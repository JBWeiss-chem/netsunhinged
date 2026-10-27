"""Regression tests against the behaviour reported by Sarkisov et al. (JACS 2014).

Run with ``python -m pytest tests`` (or ``python tests/test_netsunhinged.py``).
"""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "netsunhinged", "examples")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

from netsunhinged import MechanicalModel, Structure, deform, preset, read_cif, rigidity  # noqa: E402
from netsunhinged.cif import parse_symop  # noqa: E402
import reference_structures as make_examples  # noqa: E402


def _structure(builder):
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "x.cif")
        make_examples.write_cif(path, "x", *builder())
        with open(path, encoding="utf-8") as fh:
            return Structure(read_cif(fh.read()))


def test_symop_parser():
    R, t = parse_symop("-x+1/2, y-z, 1/2+z")
    assert np.allclose(R, [[-1, 0, 0], [0, 1, -1], [0, 0, 1]])
    assert np.allclose(t, [0.5, 0, 0.5])


def test_symops_from_hm():
    from netsunhinged.cif import symops_from_hm
    assert len(symops_from_hm("Fm-3m")) == 48 * 4
    assert len(symops_from_hm("P4/mmm")) == 16
    for bad in ("P6/mmm", "P-3m1", "P6", "Fd-3m", "P21/c"):
        try:
            symops_from_hm(bad)
        except ValueError:
            continue
        raise AssertionError(f"{bad} should be rejected")


def test_cif_expansion_and_bonds():
    path = os.path.join(ROOT, "test_cifs", "exp_21789_autored.cif")
    if not os.path.isfile(path):                  # test_cifs/ holds unpublished data
        print("skip test_cif_expansion_and_bonds (test_cifs/ not present)")
        return
    with open(path, encoding="utf-8") as fh:
        st = Structure(read_cif(fh.read()))
    assert st.N == 44                     # 11 sites x 4 symmetry operations
    assert st.dimensionality == 3
    assert sum(1 for e in st.elements if e == "Cu") == 4
    # every Cu is three-coordinate (two imidazolate N, one nitrile N)
    for s in range(st.N):
        if st.elements[s] == "Cu":
            assert len(st.nbrs[s]) == 3


def test_irmof1_rigid():
    st = _structure(make_examples.irmof1)
    m = MechanicalModel(st, preset(st, "paper"))
    r = rigidity.analyze(m)
    assert r["n_macroscopic"] == 0
    assert r["n_internal"] == 24          # axial rotation of the 24 BDC linkers
    _, res, verdict, _ = deform.run_all(m, regimes=["ABC", "A", "AB-C"])
    assert verdict == "rigid"


def test_dmof_gamma_shear():
    st = _structure(make_examples.dmof)
    m = MechanicalModel(st, preset(st, "paper"))
    r = rigidity.analyze(m)
    assert r["n_macroscopic"] == 1
    assert r["macroscopic_modes"][0]["kind"] == "shear"
    _, res, _, _ = deform.run_all(m, regimes=["Gamma", "Alpha", "A"])
    assert res["Gamma"]["classification"] == "flexible"
    assert res["A"]["classification"] == "rigid"
    # bond-bending model (rigid carboxylate binding) removes the hinge
    m2 = MechanicalModel(st, preset(st, "bondbend"))
    assert rigidity.analyze(m2)["n_macroscopic"] == 0


def test_mil53_breathing():
    st = _structure(make_examples.mil53)
    m = MechanicalModel(st, preset(st, "paper"))
    _, res, verdict, _ = deform.run_all(m, regimes=["BC", "Alpha"])
    assert res["BC"]["classification"] == "flexible"      # Fig. S11, red line
    assert res["Alpha"]["classification"] == "rigid"      # Fig. S11, black line
    assert verdict == "flexible"
    # breathing: b shrinks while c grows along the mechanism path
    pts = res["BC"]["points"]
    assert pts[0]["cell"][1] < st.cell[1] and pts[0]["cell"][2] > st.cell[2]


def _mil88b_cif():
    with open(os.path.join(EXAMPLES, "MIL-88B-2OH_COD-4107051.cif"), encoding="utf-8") as fh:
        return fh.read()


def test_mil88b_modes():
    # one cell-changing mode (swelling); the ring spins stay internal although the slightly bent
    # linkers of the experimental structure couple four of them weakly to the cell
    st = Structure(read_cif(_mil88b_cif()))
    r = rigidity.analyze(MechanicalModel(st, preset(st, "auto")))
    assert r["n_macroscopic"] == 1
    assert r["macroscopic_modes"][0]["strain_share"] > 0.5
    assert r["n_internal"] == 5 and r["n_internal_weak_strain"] == 4
    # with the carboxylate C on the ring axis: six free ring spins, no strain coupling
    # (exactly on the line from the ring centre (1/2, 1/2, 0) through C2; rounding would leave a
    # residual bend that puts the last ring spin right at the null-space tolerance)
    x = 0.5 + 0.0216 * 0.1437 / 0.0697
    straight = _mil88b_cif().replace("C1 C 0.5428(13) 0.4572(13) 0.1437(13)", f"C1 C {x!r} {1 - x!r} 0.1437")
    st2 = Structure(read_cif(straight))
    r2 = rigidity.analyze(MechanicalModel(st2, preset(st2, "auto")))
    assert r2["n_macroscopic"] == 1 and r2["n_internal"] == 6 and r2["n_internal_weak_strain"] == 0


def test_mil88b_swelling():
    st = Structure(read_cif(_mil88b_cif()))
    m = MechanicalModel(st, preset(st, "auto"))
    _, res, verdict, _ = deform.run_all(m, regimes=["AB-C", "ABC"])
    assert res["AB-C"]["classification"] == "flexible"
    assert res["ABC"]["classification"] == "rigid"
    # swelling: c shrinks while a and b grow (and vice versa)
    pts = res["AB-C"]["points"]
    assert pts[-1]["cell"][0] > st.cell[0] and pts[-1]["cell"][2] < st.cell[2]
    assert pts[0]["cell"][2] > st.cell[2]


def _zif7():
    with open(os.path.join(EXAMPLES, "ZIF-7_COD-4121899.cif"), encoding="utf-8") as fh:
        return Structure(read_cif(fh.read()))


def test_no_spurious_modes_for_rigid_framework():
    # ZIF-7 with rigid Zn-N-C and N-Zn-N angles has no mechanism at all; the three rigid
    # translations must not be reported as cell-changing modes
    st = _zif7()
    r = rigidity.analyze(MechanicalModel(st, preset(st, "paper")))
    assert r["n_mechanisms"] == 0 and r["n_macroscopic"] == 0


def test_zif7_flexible_with_zn_pivots():
    from netsunhinged import ModelOptions
    st = _zif7()
    # a flexible metal must win over the rigid-SBU rule (UI keeps rigid_sbu from the preset)
    for rigid_sbu in (False, True):
        m = MechanicalModel(st, ModelOptions(flex_angle_labels={"Zn1"}, rigid_sbu=rigid_sbu))
        assert rigidity.analyze(m)["n_macroscopic"] > 0
        _, res, _, _ = deform.run_all(m, regimes=["Alpha"])
        assert res["Alpha"]["classification"] == "flexible"


def test_remove_guests():
    st = _zif7()
    assert st.n_guest == 54        # disordered water O; half-occupied copies merged
    st2 = st.without_guests()
    assert st2.n_guest == 0 and st2.N == st.N - 54
    assert not any("closer than" in w for w in st2.warnings)
    assert st2.labels == [st.labels[s] for s in range(st.N) if st.framework[s]]


def test_path_geometry_and_export():
    from netsunhinged.geometry import export_frames, path_geometry
    st = _structure(make_examples.mil53)
    m = MechanicalModel(st, preset(st, "paper"))
    prob, res, _, _ = deform.run_all(m, regimes=["BC"])
    g = path_geometry(prob, res["BC"]["points"], res["BC"]["pattern"])
    assert g["lengths"]["b"][0] < -9.9 and g["lengths"]["c"][0] > 10       # wine-rack opening
    assert g["torsions"] and g["torsions"][0]["name"].startswith("Cr")     # hinges at Cr-O
    p = res["BC"]["points"][0]
    X, L = prob.positions(p["x"], np.array(p["cell"]), [2])
    txt = export_frames("extxyz", st.elements, [(X, L, "d=0.9")] * 2)
    assert txt.count("Lattice=") == 2
    assert export_frames("pdb", st.elements, [(X, L, "")]).count("CRYST1") == 1


def test_zeolite_preset_cristobalite_collapse():
    # ideal beta-cristobalite: rigid SiO4 tetrahedra joined at flexible Si-O-Si angles can
    # collapse (tetrahedral tilting) but not expand beyond linear Si-O-Si
    st = _structure(make_examples.cristobalite)
    assert st.is_zeolite_like()
    o = preset(st, "auto")
    assert o.flex_angle_labels == {"O1"} and not o.rigid_sbu
    _, res, _, _ = deform.run_all(MechanicalModel(st, o), regimes=["ABC", "A"])
    assert res["ABC"]["classification"] == "semi-flexible"
    assert "d < 1" in res["ABC"]["detail"]
    assert res["A"]["classification"] == "rigid"
    _, res2, verdict, _ = deform.run_all(MechanicalModel(st, preset(st, "bondbend")), regimes=["ABC"])
    assert verdict == "rigid"


def test_zeolite_preset_sodalite():
    st = _structure(make_examples.sodalite)
    m = MechanicalModel(st, preset(st, "auto"))
    assert rigidity.analyze(m)["n_macroscopic"] > 0
    _, res, _, _ = deform.run_all(m, regimes=["ABC"])
    pts = res["ABC"]["points"]
    assert all(p["F"] < 0.01 for p in pts if p["d"] < 1)          # cooperative collapse
    # MOFs are not mistaken for zeolites
    assert not _structure(make_examples.irmof1).is_zeolite_like()


def test_cif_without_symops_and_disorder_merge():
    import re
    with open(os.path.join(EXAMPLES, "MOF-5_COD-7240544.cif"), encoding="utf-8") as fh:
        text = fh.read()
    ref = Structure(read_cif(text))
    assert ref.formula(ref.framework) == "C192 H96 O104 Zn32"
    # the same structure with its symmetry-operation loop removed: operations are generated
    # from the space-group symbol (Fm-3m)
    stripped = re.sub(r"loop_\s+(_space_group_symop_\S+\s+)+(('[^']*'|\S+)\s+)*?(?=loop_|_)", "",
                      text, count=1)
    assert "_space_group_symop_operation_xyz" not in stripped and "_symmetry_equiv_pos_as_xyz" not in stripped
    st = Structure(read_cif(stripped))
    assert st.formula(st.framework) == "C192 H96 O104 Zn32"
    assert any("generated 192 operations" in w for w in st.warnings)
    # DMOF-1: half-occupied carboxylate O split across a mirror -> merged to one O per site
    with open(os.path.join(EXAMPLES, "DMOF-1_COD-7038645.cif"), encoding="utf-8") as fh:
        st = Structure(read_cif(fh.read()))
    assert sum(1 for e in st.elements if e == "O") == 8
    from netsunhinged.cif import symops_from_hm
    assert len(symops_from_hm("Fm-3m")) == 192 and len(symops_from_hm("P 4/m m m")) == 16


def test_mixed_sites_and_disorder_assemblies():
    def ex(name):
        with open(os.path.join(EXAMPLES, name), encoding="utf-8") as fh:
            return Structure(read_cif(fh.read()))
    # zeolite Rho: Si/Al share each T site (occupancies 0.79/0.21) -> one T atom per site
    st = ex("Zeolite-Rho_COD-1530545.cif")
    assert st.formula(st.framework) == "O96 Si48" and int(st.t_center.sum()) == 48
    # rho-Zn(Im)2: one imidazolate is disordered with groups 3/4 (assembly B); the first
    # group of every assembly must be kept, otherwise the framework falls apart
    st = ex("rho-ZnIm2_COD-4132678.cif")
    assert st.dimensionality == 3 and st.formula(st.framework, reduced=True) == ("C6 H6 N4 Zn", 48)


def test_auto_preset_from_motifs():
    from netsunhinged.model import auto_preset

    def ex(name):
        with open(os.path.join(EXAMPLES, name), encoding="utf-8") as fh:
            return Structure(read_cif(fh.read()))
    assert auto_preset(ex("MOF-5_COD-7240544.cif"))[0] == "paper"
    assert auto_preset(ex("Zeolite-Rho_COD-1530545.cif"))[0] == "zeolite"
    name, reason = auto_preset(ex("ZIF-8_COD-7111971.cif"))      # no carboxylate, no T-O-T
    assert name == "rigid" and "ZIF" in reason


def test_extraframework_cations_in_zeolite():
    # zeolite A with K+ / Rb+ in the pores: the cations must not become part of the framework
    with open(os.path.join(EXAMPLES, "Zeolite-A_COD-1532589.cif"), encoding="utf-8") as fh:
        st = Structure(read_cif(fh.read()))
    assert st.formula(st.framework) == "Al96 O384 Si96" and st.is_zeolite_like()
    assert int(st.bridging_tot.sum()) == 384 and st.n_guest > 0
    assert all(st.elements[s] in ("K", "Rb") for s in range(st.N) if not st.framework[s])
    assert any("extra-framework" in w for w in st.warnings)


def test_symop_parser_rejects_malicious_input():
    import time
    for bad in ["9**9**9,y,z", "__import__('os').system('ls'),y,z", "x*y,y,z", "x/0,y,z",
                "(" * 10 + "x" + ")" * 10 + ",y,z", "1e9999,y,z"]:
        t = time.time()
        try:
            parse_symop(bad)
            raise AssertionError(f"accepted {bad!r}")
        except ValueError:
            pass
        assert time.time() - t < 0.1
    R, t = parse_symop("(x+1/2)*2, 1/2-y, z")
    assert np.allclose(R[0], [2, 0, 0]) and np.allclose(t, [1, 0.5, 0])


def _call(method, path, body=b"", query=""):
    import io
    from netsunhinged import server
    out = {}

    def start(status, headers):
        out["status"], out["headers"] = int(status.split()[0]), dict(headers)
    env = {"REQUEST_METHOD": method, "PATH_INFO": path, "QUERY_STRING": query,
           "CONTENT_LENGTH": str(len(body)), "wsgi.input": io.BytesIO(body), "REMOTE_ADDR": "1.2.3.4"}
    data = b"".join(server.app(env, start))
    return out["status"], out["headers"], data


def test_web_app_security():
    import json
    status, headers, _ = _call("GET", "/")
    assert status == 200 and "Content-Security-Policy" in headers
    assert headers["X-Content-Type-Options"] == "nosniff"
    for path in ["/../netsunhinged/server.py", "/static/../server.py", "/netsunhinged/server.py",
                 "/%2e%2e/README.md", "/examples/"]:
        assert _call("GET", path)[0] == 404, path
    # malformed and oversized requests, unknown and expired sessions
    assert _call("POST", "/api/load", b"{not json")[0] == 400
    from netsunhinged.server import Config
    big = b"x" * (Config.MAX_UPLOAD + 1)
    assert _call("POST", "/api/load", big)[0] == 413
    status, _, data = _call("POST", "/api/model", json.dumps({"id": "deadbeef"}).encode())
    assert status == 410 and "expired" in json.loads(data)["error"]
    # a CIF with a script in its labels: the label reaches the client sanitised
    cif = ("data_x\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n_cell_angle_alpha 90\n"
           "_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
           "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
           "'<img src=x onerror=alert(1)>' C 0 0 0\nC2 C 0.15 0 0\n")
    status, _, data = _call("POST", "/api/load", json.dumps({"text": cif, "name": "<b>x</b>"}).encode())
    assert status == 200
    payload = json.loads(data)
    assert "<" not in payload["name"] and all("<" not in a["label"] for a in payload["atoms"])


def test_cif_with_several_structures():
    import json
    from netsunhinged.cif import cif_block_names

    def text(name):
        with open(os.path.join(EXAMPLES, name), encoding="utf-8") as fh:
            return fh.read()
    # two structures in one file, behind a block without atoms (e.g. publication data)
    multi = "data_global\n_journal_name_full 'x'\n" + text("DMOF-1_COD-7038645.cif") + "\n" + text(
        "MOF-5_COD-7240544.cif")
    names = cif_block_names(multi)
    assert len(names) == 2
    assert Structure(read_cif(multi, block=1)).formula() == Structure(
        read_cif(text("MOF-5_COD-7240544.cif"))).formula()
    try:
        read_cif(multi, block=2)
        raise AssertionError("block 3 does not exist")
    except ValueError:
        pass
    # the web API lists the structures and builds the chosen one, named after its data block
    status, _, data = _call("POST", "/api/blocks", json.dumps({"text": multi, "name": "both"}).encode())
    assert status == 200 and json.loads(data)["blocks"] == names
    status, _, data = _call("POST", "/api/load", json.dumps({"text": multi, "name": "both",
                                                            "block": 1}).encode())
    payload = json.loads(data)
    assert status == 200 and payload["name"] == names[1] and payload["stem"].startswith("both_")
    assert _call("POST", "/api/load", json.dumps({"text": multi, "block": 5}).encode())[0] == 400
    assert _call("POST", "/api/load", json.dumps({"text": multi, "block": "1"}).encode())[0] == 400
    # a single-structure file keeps the name it was opened with
    status, _, data = _call("POST", "/api/load", json.dumps({"example": "ZIF-8_COD-7111971.cif"}).encode())
    assert status == 200 and "·" not in json.loads(data)["name"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
