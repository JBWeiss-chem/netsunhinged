"""Geometric observables along a deformation path (for the per-regime graphs)."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .model import ModelOptions
from .rigidity import cell_params


def _angle(a, b, c):
    u, v = a - b, c - b
    cosv = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return np.degrees(np.arccos(np.clip(cosv, -1.0, 1.0)))


def _dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def _instances(m):
    """Flexible bond angles (at pivot atoms) and torsions (about effective hinge bonds)."""
    st = m.st
    angles = []
    for s in st.core_idx:
        if not m.flex[s]:
            continue
        nb = st.core_nbrs(s)
        for p in range(len(nb)):
            for q in range(p + 1, len(nb)):
                (j, nj, _), (k, nk, _) = nb[p], nb[q]
                ends = sorted((st.labels[j], st.labels[k]))
                angles.append((f"{ends[0]}–{st.labels[s]}–{ends[1]}", s, (j, nj), (k, nk)))
    torsions = []
    for b, state in m.bond_state.items():
        i, j, n = st.bonds[b]
        if state != "hinge" or (m._sbu_fixes_angles(i) and m._sbu_fixes_angles(j)):
            continue
        a = next(((k, nk) for k, nk, bb in st.core_nbrs(i) if bb != b), None)
        d = next(((k, nk) for k, nk, bb in st.core_nbrs(j) if bb != b), None)
        if a is None or d is None:
            continue
        key = ModelOptions.bond_key(st.labels[i], st.labels[j]).replace("|", "–")
        torsions.append((key, i, (j, n), a, d))
    return angles, torsions


def path_geometry(prob, points, pattern, top=8):
    """Cell parameters, volume, flexible angles and torsions along a scanned path.

    Angles/torsions are grouped by label class; for each class the instance with the
    largest change along the path is reported (change relative to d = 1)."""
    m = prob.m
    relax = [i for i, ch in enumerate(pattern) if ch == "R"]
    angles, torsions = _instances(m)
    d = [p["d"] for p in points]
    cells, vols = [], []
    A = np.zeros((len(points), len(angles)))
    T = np.zeros((len(points), len(torsions)))
    for t, p in enumerate(points):
        X, L = prob.positions(p["x"], np.array(p["cell"]), relax)

        def pos(s, n=(0, 0, 0)):
            return X[s] + np.array(n) @ L

        cells.append(cell_params(L))
        vols.append(abs(np.linalg.det(L)))
        for k, (_, s, (j, nj), (kk, nk)) in enumerate(angles):
            A[t, k] = _angle(pos(j, nj), pos(s), pos(kk, nk))
        for k, (_, i, (j, n), (a, na), (dd, nd)) in enumerate(torsions):
            T[t, k] = _dihedral(pos(a, na), pos(i), pos(j, n), pos(dd, np.add(n, nd)))
    i1 = int(np.argmin(np.abs(np.array(d) - 1.0)))
    A -= A[i1]
    T -= T[i1]
    T = (T + 180.0) % 360.0 - 180.0

    def classes(values, inst):
        groups = defaultdict(list)
        for k, item in enumerate(inst):
            groups[item[0]].append(k)
        out = []
        for key, ks in groups.items():
            sub = values[:, ks]
            best = int(np.argmax(np.abs(sub).max(axis=0)))
            out.append({"name": key, "n": len(ks), "y": np.round(sub[:, best], 3).tolist(),
                        "max": float(np.abs(sub).max())})
        out = [c for c in out if c["max"] > 0.01]
        out.sort(key=lambda c: -c["max"])
        return out[:top]

    cells = np.array(cells)
    c0 = cells[i1]
    return {
        "d": d,
        "lengths": {k: np.round(100 * (cells[:, i] / c0[i] - 1), 3).tolist()
                    for i, k in enumerate("abc")},
        "angles_cell": {k: np.round(cells[:, i + 3] - c0[i + 3], 3).tolist()
                        for i, k in enumerate(["alpha", "beta", "gamma"])},
        "volume": np.round(100 * (np.array(vols) / vols[i1] - 1), 3).tolist(),
        "flex_angles": classes(A, angles),
        "torsions": classes(T, torsions),
        "n_flex_angles": len(angles), "n_torsions": len(torsions),
    }


# ------------------------------------------------------------------ export

def export_frames(fmt, symbols, frames, name="netsunhinged"):
    """Write frames [(xyz Nx3, L 3x3, comment)] to extxyz / pdb (text) or ASE traj (bytes)."""
    if fmt == "extxyz":
        out = []
        for X, L, note in frames:
            lat = " ".join(f"{v:.6f}" for v in np.asarray(L).ravel())
            out.append(str(len(symbols)))
            out.append(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T" {note}')
            out += [f"{s:2s} {x:14.6f} {y:14.6f} {z:14.6f}" for s, (x, y, z) in zip(symbols, X)]
        return "\n".join(out) + "\n"
    if fmt == "pdb":
        out = [f"REMARK   netsunhinged trajectory {name}"]
        for f, (X, L, note) in enumerate(frames, 1):
            a, b, c, al, be, ga = cell_params(np.asarray(L))
            out.append(f"MODEL     {f:4d}")
            out.append(f"REMARK   {note}")
            out.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{al:7.2f}{be:7.2f}{ga:7.2f} P 1           1")
            for k, (s, (x, y, z)) in enumerate(zip(symbols, X), 1):
                out.append(f"ATOM  {k % 100000:5d} {s:>2s}   MOF     1    {x:8.3f}{y:8.3f}{z:8.3f}"
                           f"  1.00  0.00          {s:>2s}")
            out.append("ENDMDL")
        return "\n".join(out) + "\nEND\n"
    if fmt == "traj":
        try:
            from ase import Atoms
            from ase.io import write
        except ImportError:
            raise RuntimeError("ASE is not installed (pip install ase). Export as extended XYZ "
                               "instead - ASE reads it with ase.io.read(file, index=':').") from None
        import os
        import tempfile
        images = [Atoms(symbols=symbols, positions=X, cell=L, pbc=True) for X, L, _ in frames]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.traj")
            write(path, images)
            with open(path, "rb") as fh:
                return fh.read()
    raise ValueError(f"unknown format {fmt}")
