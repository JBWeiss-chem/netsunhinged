"""Web application: JSON API + static three.js front end.

One WSGI application (``app``) serves both modes:

* local:       ``python -m netsunhinged serve``  (threaded stdlib server)
* production:  ``gunicorn -w 1 --threads 8 netsunhinged.server:app`` behind an HTTPS proxy.
  Exactly one worker process: sessions and analysis jobs live in memory.

Limits are set with environment variables (see ``Config``).
"""

from __future__ import annotations

import base64
import hashlib
import json
import logging
import mimetypes
import os
import re
import secrets
import threading
import time
import webbrowser
from collections import OrderedDict, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import parse_qs

import numpy as np

from . import deform, rigidity
from .cif import cif_block_names, read_cif, safe_text
from .elements import COVALENT_RADII, color
from .geometry import export_frames, path_geometry
from .model import PRESET_LABELS, MechanicalModel, ModelOptions, auto_preset, preset
from .structure import Structure

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(HERE, "static")
ROOT = os.path.dirname(HERE)
# netsunhinged/examples/ is published and installed with the package; examples_local/ next to the
# package in a source checkout (if present) holds files that must not be published (e.g. CSD
# downloads) and is never committed or copied into the Space or Docker image
EXAMPLE_DIRS = [d for d in (os.path.join(HERE, "examples"), os.path.join(ROOT, "examples_local"))
                if os.path.isdir(d)]
log = logging.getLogger("netsunhinged")


def _env(name, default, cast=float):
    try:
        return cast(os.environ.get(name, default))
    except ValueError:
        return cast(default)


class Config:
    MAX_UPLOAD = int(_env("NU_MAX_UPLOAD_MB", 5) * 1024 * 1024)   # request body size
    MAX_ATOMS = int(_env("NU_MAX_ATOMS", 6000))       # atoms after symmetry expansion
    MAX_DOF = int(_env("NU_MAX_DOF", 6000))           # unknowns of the infinitesimal analysis
    WORKERS = int(_env("NU_WORKERS", 1))              # parallel deformation analyses
    MAX_QUEUE = int(_env("NU_MAX_QUEUE", 4))          # waiting analyses
    JOB_TIMEOUT = _env("NU_JOB_TIMEOUT", 600)         # seconds per analysis
    HEAVY = int(_env("NU_HEAVY_CONCURRENCY", 2))      # concurrent load/model/export requests
    TTL = _env("NU_TTL", 3600)                        # idle seconds before a session is dropped
    MAX_SESSIONS = int(_env("NU_MAX_SESSIONS", 200))  # stored structures / jobs (LRU)
    RATE = int(_env("NU_RATE_PER_MIN", 60))           # heavy requests per client and minute
    RATE_FRAMES = int(_env("NU_RATE_FRAMES_PER_MIN", 600))  # deformed-frame requests (slider)
    TRUST_PROXY = bool(int(_env("NU_TRUST_PROXY", 0)))  # use X-Forwarded-For for rate limits
    FRAME_ANCESTORS = os.environ.get("NU_FRAME_ANCESTORS", "'self'")   # who may embed the page


class ClientError(Exception):
    """Errors whose message is safe and useful to show to the user."""

    def __init__(self, message, status=400):
        super().__init__(message)
        self.status = status


class Store:
    """Thread-safe in-memory store with idle expiry and an LRU size cap."""

    def __init__(self, ttl, cap):
        self.ttl, self.cap = ttl, cap
        self._d = OrderedDict()                       # id -> [obj, last access]
        self._lock = threading.Lock()

    def put(self, obj):
        key = secrets.token_hex(16)                   # 128-bit, unguessable
        with self._lock:
            self._d[key] = [obj, time.time()]
            while len(self._d) > self.cap:
                self._d.popitem(last=False)
        return key

    def get(self, key):
        with self._lock:
            item = self._d.get(str(key))
            if item is None:
                raise ClientError("This session has expired. Please load the structure again.", 410)
            item[1] = time.time()
            self._d.move_to_end(str(key))
            return item[0]

    def purge(self):
        cut = time.time() - self.ttl
        with self._lock:
            for k in [k for k, (_, t) in self._d.items() if t < cut]:
                del self._d[k]

    def __len__(self):
        return len(self._d)


STRUCTURES = Store(Config.TTL, Config.MAX_SESSIONS)
JOBS = Store(Config.TTL, Config.MAX_SESSIONS)
_heavy = threading.BoundedSemaphore(Config.HEAVY)
_executor = None
_pending = [0]                                        # queued + running analyses
_pending_lock = threading.Lock()
_started = threading.Event()


def _start_background():
    """Lazily start the worker pool and the janitor (after a gunicorn fork)."""
    global _executor
    if _started.is_set():
        return
    with _pending_lock:
        if _started.is_set():
            return
        _executor = ThreadPoolExecutor(max_workers=Config.WORKERS)

        def janitor():
            while True:
                time.sleep(60)
                STRUCTURES.purge()
                JOBS.purge()
        threading.Thread(target=janitor, daemon=True).start()
        _started.set()


def _r(x, nd=4):
    return np.round(np.asarray(x, float), nd).tolist()


# ------------------------------------------------------------------ payloads

def structure_payload(sid, st: Structure):
    atoms = []
    for s in range(st.N):
        role = "core" if st.core[s] else ("pruned" if st.framework[s] else "guest")
        atoms.append({
            "el": st.elements[s], "label": st.labels[s], "xyz": _r(st.cart[s]),
            "color": color(st.elements[s]), "r": COVALENT_RADII.get(st.elements[s], 1.0),
            "role": role, "env": st.env[s], "metal": bool(st.metal[s]),
            "sbu": bool(st.sbu_atom[s]), "bridging": bool(st.bridging_donor[s]),
        })
    bonds = []
    for b, (i, j, n) in enumerate(st.bonds):
        bonds.append({"i": i, "j": j, "n": list(n), "len": round(st.bond_len[b], 3),
                      "core": b in st.auto_state, "auto": st.auto_state.get(b),
                      "reason": st.auto_reason.get(b)})
    labels = {}
    for s in st.core_idx:
        d = labels.setdefault(st.labels[s], {"el": st.elements[s], "env": st.env[s], "count": 0,
                                            "bridging": False, "sbu": False, "metal": False,
                                            "sbu_internal": True, "relevant": False})
        d["count"] += 1
        d["relevant"] |= bool(st.metal[s] or st.sbu_atom[s] or st.t_center[s] or st.bridging_tot[s])
        # all bond angles of this atom lie inside one inorganic SBU (metals, mu-O, ...)
        d["sbu_internal"] &= bool(st.sbu_atom[s]) and all(
            st.sbu_atom[j] and st.sbu_id[j] == st.sbu_id[s] and (st.metal[s] or st.metal[j])
            for j, _, _ in st.core_nbrs(s))
        d["bridging"] |= bool(st.bridging_donor[s])
        d["sbu"] |= bool(st.sbu_atom[s])
        d["metal"] |= bool(st.metal[s])
    bclasses = {}
    for b, auto in st.auto_state.items():
        i, j, _ = st.bonds[b]
        key = ModelOptions.bond_key(st.labels[i], st.labels[j])
        d = bclasses.setdefault(key, {"count": 0, "auto": auto, "reason": st.auto_reason[b],
                                      "elements": "-".join(sorted((st.elements[i], st.elements[j]))),
                                      "metal": bool(st.metal[i] or st.metal[j]),
                                      "relevant": bool(st.metal[i] or st.metal[j]
                                                       or st.t_center[i] or st.t_center[j])})
        d["count"] += 1
        if d["auto"] != auto:
            d["auto"] = "mixed"
    return {
        "id": sid, "name": st.name, "cell": _r(st.cell, 4), "L": _r(st.L),
        "formula": st.formula(), "framework_formula": st.formula(st.framework),
        "framework_formula_reduced": st.formula(st.framework, reduced=True)[0],
        "z": st.formula(st.framework, reduced=True)[1],
        "stem": getattr(st, "stem", st.name),
        "has_metal": bool((st.metal & st.core).any()),
        "has_tot": bool(st.bridging_tot.any()),
        "has_carboxylate": bool(st.bridging_donor.any()),
        "n_atoms": st.N, "n_core": int(st.core.sum()), "n_guest": st.n_guest,
        "n_pruned": int((st.framework & ~st.core).sum()), "n_sbu": st.n_sbu,
        "dimensionality": st.dimensionality, "warnings": st.warnings,
        "atoms": atoms, "bonds": bonds, "labels": labels, "bond_classes": bclasses,
        "topology": st.topology(),
    }


def model_payload(st: Structure, m: MechanicalModel, rig):
    def mode_out(md):
        out = {k: v for k, v in md.items() if k not in ("vec",)}
        out["disp"] = _r(md["disp"], 4)
        return out
    atom_status, bond_status = m.dof_status()
    return {
        "atom_status": atom_status,
        "bond_status": bond_status,
        "summary": m.summary(),
        "options": m.opts.to_dict(),
        "site_body": m.site_body().tolist(),
        "bodies": [{"kind": B.kind, "linear": bool(B.linear), "n": len(B.atoms)} for B in m.bodies],
        "joints": [{"kind": j["kind"], "points": _r(j["points"], 3)} for j in m.joints
                   if j["kind"] != "rigid"],
        "bond_state": {str(b): s for b, s in m.bond_state.items()},
        "flex": m.flex.tolist(),
        "rigidity": {**{k: v for k, v in rig.items() if k not in ("macroscopic_modes", "internal_modes")},
                     "macroscopic_modes": [mode_out(md) for md in rig["macroscopic_modes"]],
                     "internal_modes": [{"disp": _r(md["disp"], 3)} for md in rig["internal_modes"]]},
    }


# ---------------------------------------------------------------- handlers

def _source(body):
    """-> (CIF text, display name or None, file stem or None) of an upload or an example."""
    text = body.get("text")
    name = body.get("name")
    stem = name
    if text is None and body.get("example"):
        if not isinstance(body["example"], str) or not body["example"].lower().endswith(".cif"):
            raise ClientError("unknown example")
        fname = os.path.basename(body["example"])
        for d in EXAMPLE_DIRS:
            p = os.path.join(d, fname)
            if os.path.isfile(p):
                with open(p, encoding="utf-8", errors="replace") as fh:
                    text = fh.read()
                break
        # examples are shown with their catalogue name (the CIF data block is often 'cmpd_2')
        stem = fname[:-4]
        name = _catalog().get("structures", {}).get(fname, {}).get("name", stem)
    if not isinstance(text, str):
        raise ClientError("no CIF given")
    if not isinstance(name, str):
        name = stem = None
    return text, name, stem


def api_blocks(body):
    """The structures (data blocks with atom sites) in a CIF, without building them."""
    text, name, _ = _source(body)
    blocks = cif_block_names(text)
    if not blocks:
        raise ValueError("no data block with _atom_site_fract_x found")
    return {"name": safe_text(name) if name else None, "blocks": blocks}


def api_load(body):
    text, name, stem = _source(body)
    block = body.get("block", 0)
    if not isinstance(block, int) or isinstance(block, bool):
        raise ClientError("invalid block number")
    st = Structure(read_cif(text, max_atoms=Config.MAX_ATOMS, block=block))
    blocks = cif_block_names(text)
    if len(blocks) > 1:
        # one of several structures in a file: it is named after its data block
        name = f"{name} · {st.name}" if name and body.get("example") else st.name
        stem = f"{stem}_{st.name}" if stem else st.name
    if name:
        st.name = safe_text(name)
    st.stem = re.sub(r"[^\w.-]+", "_", stem or st.name).strip("_")[:60] or "structure"
    return structure_payload(STRUCTURES.put(st), st)


def api_remove_guests(body):
    old = STRUCTURES.get(body.get("id"))
    st = old.without_guests()
    st.stem = getattr(old, "stem", old.name)
    return structure_payload(STRUCTURES.put(st), st)


def _options(st, body):
    if body.get("preset"):
        if body["preset"] not in ("auto", "paper", "zeolite", "rigid", "bondbend", "donor_pivots"):
            raise ClientError("unknown preset")
        return preset(st, body["preset"])
    o = body.get("options", {})
    if not isinstance(o, dict):
        raise ClientError("bad options")
    return ModelOptions.from_dict(o)


def _infinitesimal(m):
    """Infinitesimal analysis, skipped for models too large for the dense eigensolver."""
    if 6 * len(m.bodies) + 6 > Config.MAX_DOF:
        return {"skipped": True, "n_dof": 6 * len(m.bodies) + 6, "n_constraint_rows": 0, "rank": 0,
                "n_mechanisms": 0, "n_macroscopic": 0, "n_internal": 0, "smallest_nonzero_sv": [],
                "macroscopic_modes": [], "internal_modes": []}
    return rigidity.analyze(m)


def api_model(body):
    st = STRUCTURES.get(body.get("id"))
    m = MechanicalModel(st, _options(st, body))
    rig = _infinitesimal(m)
    out = model_payload(st, m, rig)
    name = body.get("preset")
    if name:
        reason = ""
        if name == "auto":
            name, reason = auto_preset(st)
        out["preset"] = {"name": name, "label": PRESET_LABELS.get(name, name), "reason": reason}
    return out


def _catalog():
    """Catalogues of all example folders, merged (categories in folder order)."""
    merged = {"categories": [], "structures": {}}
    for d in EXAMPLE_DIRS:
        path = os.path.join(d, "catalog.json")
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as fh:
                c = json.load(fh)
            merged["categories"] += [x for x in c.get("categories", []) if x not in merged["categories"]]
            merged["structures"].update(c.get("structures", {}))
    return merged


def api_examples(_):
    """Example CIFs grouped by the categories of examples/catalog.json.

    Files without a catalogue entry are listed under 'Other' with their file name."""
    files = []
    for d in EXAMPLE_DIRS:
        if os.path.isdir(d):
            files += sorted(f for f in os.listdir(d) if f.lower().endswith(".cif"))
    catalog = _catalog()
    entries = catalog.get("structures", {})
    groups = {c: [] for c in catalog.get("categories", [])}
    # catalogue order first (as curated), then any other files alphabetically
    order = {f: k for k, f in enumerate(entries)}
    files.sort(key=lambda f: (order.get(f, len(order)), f))
    for f in files:
        e = entries.get(f)
        if e:
            label = f"{e['name']} · {e['note']} · {e['id']}" if e.get("note") else f"{e['name']} · {e['id']}"
            groups.setdefault(e["category"], []).append({"file": f, "label": label})
        else:
            groups.setdefault("Other", []).append({"file": f, "label": f[:-4]})
    return {"groups": [{"category": c, "items": items} for c, items in groups.items() if items]}


class JobTimeout(Exception):
    pass


def _settings(raw):
    """Validated analysis settings (bounded work per job)."""
    raw = raw if isinstance(raw, dict) else {}

    def num(key, default, lo, hi):
        try:
            v = float(raw.get(key, default))
        except (TypeError, ValueError):
            v = default
        return min(hi, max(lo, v)) if np.isfinite(v) else default
    st = {"dmin": num("dmin", 0.9, 0.5, 0.999), "dmax": num("dmax", 1.1, 1.001, 1.5),
          "step": num("step", 0.02, 0.005, 0.1), "thr": num("thr", 0.01, 1e-12, 1.0),
          "thr2": num("thr2", 0.05, 1e-12, 1.0),
          "metric": "absolute" if raw.get("metric") == "absolute" else "relative",
          "follow_modes": bool(raw.get("follow_modes", True)),
          "extended": bool(raw.get("extended", False)),
          "paper_regimes": bool(raw.get("paper_regimes", True))}
    if (st["dmax"] - st["dmin"]) / st["step"] > 100:
        st["step"] = (st["dmax"] - st["dmin"]) / 100
    return st


def _run_job(job, st, opts, settings):
    job["status"], job["message"] = "running", "infinitesimal analysis"
    deadline = time.time() + Config.JOB_TIMEOUT
    try:
        m = MechanicalModel(st, opts)
        modes = _infinitesimal(m)["macroscopic_modes"] if settings["follow_modes"] else []
        regimes = list(deform.REGIMES) if settings["paper_regimes"] else []
        if settings["extended"]:
            regimes += list(deform.EXTENDED)
        job["message"] = "deformation scans"

        def prog(f):
            job["progress"] = f
            if time.time() > deadline:
                raise JobTimeout()

        prob, res, verdict, dt = deform.run_all(
            m, regimes=regimes, dmin=settings["dmin"], dmax=settings["dmax"], step=settings["step"],
            metric=settings["metric"], thr=settings["thr"], thr2=settings["thr2"],
            modes=modes, progress=prog)
        job["prob"] = prob
        job["raw"] = res
        job["message"] = "path geometry"
        geom = {k: path_geometry(prob, r["points"], r["pattern"]) for k, r in res.items()}
        job["result"] = {
            "verdict": verdict, "seconds": round(dt, 2),
            "regimes": {k: {"description": r["description"], "pattern": r["pattern"],
                            "classification": r["classification"], "detail": r["detail"],
                            "d": [p["d"] for p in r["points"]],
                            "E": [p["E"] for p in r["points"]],
                            "F": [p["F"] for p in r["points"]],
                            "cell": [_r(p["cell"], 4) for p in r["points"]],
                            "geometry": geom[k]}
                        for k, r in res.items()},
        }
        job["status"] = "done"
    except JobTimeout:
        job["status"] = "error"
        job["message"] = (f"The analysis exceeded the time limit of {Config.JOB_TIMEOUT / 60:.0f} min. "
                          "Use a coarser step or fewer regimes.")
    except Exception:
        ref = secrets.token_hex(4)
        log.exception("analysis failed (ref %s)", ref)
        job["status"] = "error"
        job["message"] = f"The analysis failed (internal error, reference {ref})."
    finally:
        with _pending_lock:
            _pending[0] -= 1


def api_analyze(body):
    st = STRUCTURES.get(body.get("id"))
    opts = _options(st, body)
    settings = _settings(body.get("settings"))
    if not (settings["paper_regimes"] or settings["extended"] or settings["follow_modes"]):
        raise ClientError("Select at least one kind of deformation regime.")
    _start_background()
    with _pending_lock:
        if _pending[0] >= Config.WORKERS + Config.MAX_QUEUE:
            raise ClientError("The server is busy with other analyses. Please try again in a few minutes.", 503)
        _pending[0] += 1
    job = {"status": "queued", "progress": 0.0, "message": "waiting for a free worker"}
    jid = JOBS.put(job)
    _executor.submit(_run_job, job, st, opts, settings)
    return {"job": jid}


def api_job(q):
    job = JOBS.get(q.get("job"))
    out = {"status": job["status"], "progress": job["progress"], "message": job["message"]}
    if job["status"] == "done":
        out["result"] = job["result"]
    return out


def _finished(job):
    if job.get("status") != "done":
        raise ClientError("The analysis has not finished.")
    return job


def api_frame(q):
    job = _finished(JOBS.get(q.get("job")))
    reg = job["raw"].get(q.get("regime"))
    if reg is None:
        raise ClientError("unknown regime")
    try:
        d = float(q.get("d", 1.0))
    except ValueError:
        raise ClientError("bad order parameter") from None
    pt = min(reg["points"], key=lambda p: abs(p["d"] - d))
    relax = [i for i, ch in enumerate(reg["pattern"]) if ch == "R"]
    prob = job["prob"]
    X, L = prob.positions(pt["x"], np.array(pt["cell"]), relax)
    return {"d": pt["d"], "E": pt["E"], "F": pt["F"], "xyz": _r(X, 3), "L": _r(L, 4),
            "cell": _r(pt["cell"], 4)}


def api_export(body):
    """Trajectory export. Either a scanned regime (job + regime) or explicit frames."""
    st = STRUCTURES.get(body.get("id"))
    symbols = list(st.elements)
    fmt = body.get("fmt", "extxyz")
    if fmt not in ("extxyz", "pdb", "traj"):
        raise ClientError("unknown export format")
    if body.get("job"):
        job = _finished(JOBS.get(body["job"]))
        reg = job["raw"].get(body.get("regime"))
        if reg is None:
            raise ClientError("unknown regime")
        relax = [i for i, ch in enumerate(reg["pattern"]) if ch == "R"]
        frames = []
        for p in reg["points"]:
            X, L = job["prob"].positions(p["x"], np.array(p["cell"]), relax)
            frames.append((X, L, f'd={p["d"]:.4f} E={p["E"]:.3e} F={p["F"]:.3e}'))
    else:
        raw = body.get("frames")
        if not isinstance(raw, list) or not 0 < len(raw) <= 200:
            raise ClientError("between 1 and 200 frames can be exported")
        frames = []
        for f in raw:
            if not isinstance(f, dict):
                raise ClientError("bad frame")
            X, L = np.asarray(f.get("xyz"), float), np.asarray(f.get("L"), float)
            if X.shape != (len(symbols), 3) or L.shape != (3, 3) or not (np.isfinite(X).all() and np.isfinite(L).all()):
                raise ClientError("bad frame")
            frames.append((X, L, safe_text(f.get("note", ""), 80)))
    try:
        data = export_frames(fmt, symbols, frames, st.name)
    except RuntimeError as exc:                      # e.g. ASE not installed
        raise ClientError(str(exc)) from None
    if isinstance(data, bytes):
        return {"base64": base64.b64encode(data).decode()}
    return {"text": data}


POST_ROUTES = {"/api/load": api_load, "/api/blocks": api_blocks, "/api/model": api_model,
               "/api/analyze": api_analyze,
               "/api/remove_guests": api_remove_guests, "/api/export": api_export}
GET_ROUTES = {"/api/examples": api_examples, "/api/job": api_job, "/api/frame": api_frame}
HEAVY_ROUTES = {"/api/load", "/api/model", "/api/remove_guests", "/api/export"}
# per-client request budgets per minute. /api/frame only places the bodies of an already finished
# analysis (milliseconds), but the frame slider calls it often, so it has its own, larger budget
RATE_LIMITS = {**{p: ("heavy", Config.RATE) for p in HEAVY_ROUTES | {"/api/analyze"}},
               "/api/frame": ("frames", Config.RATE_FRAMES),
               # listing the structures of a CIF only tokenises it; opening many files at once
               # sends one request per file
               "/api/blocks": ("frames", Config.RATE_FRAMES)}

# only these files are served (no directory listing, no path tricks)
STATIC_FILES = {"/": "index.html", "/index.html": "index.html", "/app.js": "app.js",
                "/style.css": "style.css",
                # three.js is served locally: no visitor IP addresses go to a third-party CDN
                "/vendor/three/three.module.min.js": "vendor/three/three.module.min.js",
                "/vendor/three/addons/controls/OrbitControls.js":
                    "vendor/three/addons/controls/OrbitControls.js"}


def _csp():
    """Content-Security-Policy; the one inline script (the import map) is allowed by hash."""
    with open(os.path.join(STATIC, "index.html"), encoding="utf-8") as fh:
        html = fh.read()
    hashes = [base64.b64encode(hashlib.sha256(m.encode()).digest()).decode()
              for m in re.findall(r"<script type=\"importmap\">(.*?)</script>", html, re.S)]
    script = " ".join(["'self'"] + [f"'sha256-{h}'" for h in hashes])
    return ("default-src 'self'; "
            f"script-src {script}; "
            "style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; "
            "object-src 'none'; base-uri 'none'; form-action 'none'; "
            f"frame-ancestors {Config.FRAME_ANCESTORS}")


CSP = _csp()
SECURITY_HEADERS = [
    ("X-Content-Type-Options", "nosniff"),
    ("Referrer-Policy", "no-referrer"),
    ("Content-Security-Policy", CSP),
    ("Permissions-Policy", "camera=(), microphone=(), geolocation=()"),
    ("Cross-Origin-Opener-Policy", "same-origin"),
]

_hits = defaultdict(deque)
_hits_lock = threading.Lock()


def _rate_ok(client, bucket="heavy", limit=None):
    now = time.time()
    with _hits_lock:
        q = _hits[(bucket, client)]
        while q and q[0] < now - 60:
            q.popleft()
        if len(q) >= (Config.RATE if limit is None else limit):
            return False
        q.append(now)
        if len(_hits) > 10000:                    # forget idle clients
            for k in [k for k, v in _hits.items() if not v or v[-1] < now - 60]:
                del _hits[k]
        return True


def _client(environ):
    if Config.TRUST_PROXY and environ.get("HTTP_X_FORWARDED_FOR"):
        return environ["HTTP_X_FORWARDED_FOR"].split(",")[0].strip()
    return environ.get("REMOTE_ADDR", "?")


def _handle(method, path, query, body, client):
    """-> (status, content type, bytes). Only ClientError messages reach the user."""
    def js(obj, status=200):
        return status, "application/json", json.dumps(obj).encode()

    if path.startswith("/api/"):
        routes = POST_ROUTES if method == "POST" else GET_ROUTES if method == "GET" else {}
        fn = routes.get(path)
        if fn is None:
            return js({"error": "not found"}, 404)
        if path in RATE_LIMITS and not _rate_ok(client, *RATE_LIMITS[path]):
            return js({"error": "Too many requests. Please wait a minute."}, 429)
        heavy = path in HEAVY_ROUTES
        if heavy and not _heavy.acquire(timeout=30):
            return js({"error": "The server is busy. Please try again shortly."}, 503)
        try:
            if method == "POST":
                try:
                    data = json.loads(body or b"{}")
                except ValueError:
                    raise ClientError("malformed request") from None
                if not isinstance(data, dict):
                    raise ClientError("malformed request")
                return js(fn(data))
            return js(fn({k: v[0] for k, v in parse_qs(query).items()}))
        except ClientError as exc:
            return js({"error": str(exc)}, exc.status)
        except ValueError as exc:                    # bad CIF etc.: our own, sanitised messages
            return js({"error": safe_text(str(exc), 300)}, 400)
        except Exception:
            ref = secrets.token_hex(4)
            log.exception("request %s failed (ref %s)", path, ref)
            return js({"error": f"Internal error (reference {ref})."}, 500)
        finally:
            if heavy:
                _heavy.release()
    if method in ("GET", "HEAD") and path in STATIC_FILES:
        fname = STATIC_FILES[path]
        with open(os.path.join(STATIC, fname), "rb") as fh:
            data = fh.read()
        ctype = mimetypes.guess_type(fname)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype.endswith("javascript"):
            ctype += "; charset=utf-8"
        return 200, ctype, data
    return 404, "text/plain; charset=utf-8", b"not found"


_REASONS = {200: "OK", 400: "Bad Request", 404: "Not Found", 410: "Gone", 413: "Payload Too Large",
            429: "Too Many Requests", 500: "Internal Server Error", 503: "Service Unavailable"}


def app(environ, start_response):
    """WSGI entry point (gunicorn: ``netsunhinged.server:app``)."""
    _start_background()
    method = environ.get("REQUEST_METHOD", "GET")
    path = environ.get("PATH_INFO", "/") or "/"
    try:
        length = int(environ.get("CONTENT_LENGTH") or 0)
    except ValueError:
        length = 0
    if length > Config.MAX_UPLOAD:
        status, ctype = 413, "application/json"
        data = json.dumps({"error": f"The file is too large (limit {Config.MAX_UPLOAD // 2**20} MB)."}).encode()
    else:
        body = environ["wsgi.input"].read(length) if length > 0 else b""
        status, ctype, data = _handle(method, path, environ.get("QUERY_STRING", ""), body,
                                      _client(environ))
    headers = [("Content-Type", ctype), ("Content-Length", str(len(data)))] + SECURITY_HEADERS
    if path.startswith("/api/"):
        headers.append(("Cache-Control", "no-store"))
    start_response(f"{status} {_REASONS.get(status, '')}", headers)
    return [b"" if method == "HEAD" else data]


def serve(host="127.0.0.1", port=8765, open_browser=True):
    """Local mode: the WSGI app on a threaded stdlib server."""
    from socketserver import ThreadingMixIn
    from wsgiref.simple_server import WSGIRequestHandler, WSGIServer, make_server

    class ThreadingWSGIServer(ThreadingMixIn, WSGIServer):
        daemon_threads = True

    class QuietHandler(WSGIRequestHandler):
        def log_message(self, *args):
            pass

    logging.basicConfig(level=logging.INFO)
    httpd = make_server(host, port, app, server_class=ThreadingWSGIServer, handler_class=QuietHandler)
    url = f"http://{host}:{port}/"
    print(f"NetsUnhinged running at {url}  (Ctrl+C to stop)")
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
