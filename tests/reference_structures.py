"""Idealised reference frameworks used by the regression tests (P1 CIFs).

* IRMOF-1 (MOF-5): Zn4O(BDC)3, pcu, Fm-3m conventional cell.  Rigid (Sarkisov et al.).
* DMOF-type: Cu2(BDC)2(pillar), paddlewheel sql layers pillared by a linear N-C-C-N
  linker.  All carboxylate hinges are parallel (along c): flexible (wine-rack).
* MIL-53: trans corner-sharing Cr-OH chains + BDC: flexible (BC regime).
* beta-cristobalite (ideal, linear Si-O-Si) and siliceous sodalite: zeolite-type nets.

Geometries are built from standard bond lengths, not refined crystal structures; they
are test data, not literature structures (those live in examples/).
Usage: python tests/reference_structures.py OUTPUT_DIR
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def write_cif(path, name, cell, atoms):
    a, b, c = cell
    L = np.diag([a, b, c])
    inv = np.linalg.inv(L)
    lines = [f"data_{name}", f"_cell_length_a {a:.4f}", f"_cell_length_b {b:.4f}",
             f"_cell_length_c {c:.4f}", "_cell_angle_alpha 90", "_cell_angle_beta 90",
             "_cell_angle_gamma 90", "_symmetry_space_group_name_H-M 'P 1'",
             "loop_", "_symmetry_equiv_pos_as_xyz", "x,y,z",
             "loop_", "_atom_site_label", "_atom_site_type_symbol",
             "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    for label, el, x in atoms:
        f = (np.asarray(x) @ inv) % 1.0
        lines.append(f"{label} {el} {f[0]:.6f} {f[1]:.6f} {f[2]:.6f}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def ring(center, ax, w, prefix):
    """Benzene ring (para positions along ax) with 4 H; returns ipso atoms + others."""
    atoms = []
    for s in (-1, 1):
        atoms.append((f"{prefix}C2", "C", center + s * 1.39 * ax))
        for t in (-1, 1):
            p = center + s * 0.695 * ax + t * 1.204 * w
            atoms.append((f"{prefix}C3", "C", p))
            radial = (p - center) / np.linalg.norm(p - center)
            atoms.append((f"{prefix}H3", "H", p + 1.08 * radial))
    return atoms


def irmof1():
    a0 = 12.916
    A = 2 * a0
    delta = 1.94 / np.sqrt(3)
    even = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    atoms = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                c = a0 * np.array([i, j, k], float)
                vs = even if (i + j + k) % 2 == 0 else -even
                atoms.append(("O1", "O", c))
                for v in vs:
                    atoms.append(("Zn1", "Zn", c + delta * v))
                for ax_i in range(3):
                    for sgn in (1, -1):
                        e = np.zeros(3)
                        e[ax_i] = sgn
                        atoms.append(("C1", "C", c + 3.59 * e))
                        side = [v for v in vs if v[ax_i] == sgn]
                        for v in side:
                            vp = v.copy()
                            vp[ax_i] = 0
                            atoms.append(("O2", "O", c + 3.004 * e + 0.797 * vp))
                        if sgn == 1:
                            vp = side[0].copy()
                            vp[ax_i] = 0
                            w = vp / np.linalg.norm(vp)
                            atoms += ring(c + 0.5 * a0 * e, e, w, "")
    return (A, A, A), atoms


def dmof():
    zc, zo = 1.32, 1.12
    a = 2 * (2.536 + 1.49 + 1.39)
    c = 2 * zc + 2 * 2.0 + 2 * 1.15 + 1.38
    atoms = [("Cu1", "Cu", [0, 0, zc]), ("Cu1", "Cu", [0, 0, -zc])]
    ez = np.array([0, 0, 1.0])
    for e in (np.array([1.0, 0, 0]), np.array([0, 1.0, 0])):
        for sgn in (1, -1):
            atoms.append(("C1", "C", sgn * 2.536 * e))
            for t in (1, -1):
                atoms.append(("O1", "O", sgn * 1.95 * e + t * zo * ez))
        atoms += ring(0.5 * a * e, e, ez, "")
    z = zc + 2.0
    atoms += [("N1", "N", [0, 0, z]), ("C4", "C", [0, 0, z + 1.15]),
              ("C4", "C", [0, 0, z + 1.15 + 1.38]), ("N1", "N", [0, 0, z + 2 * 1.15 + 1.38])]
    return (a, a, c), atoms


def mil53(theta_deg=35.0):
    """Cr(OH)(BDC): trans corner-sharing CrO4(OH)2 chains along a, wine-rack in bc."""
    th = np.radians(theta_deg)
    half = 1.7                      # half Cr...Cr distance along the chain
    a = 4 * half
    reach = 2.45 + 1.49 + 1.39      # rod axis -> ring centre
    # neighbouring rods are 2*reach apart along +-u; rod 2 sits at (b/2, c/2)
    b, c = 4 * reach * np.cos(th), 4 * reach * np.sin(th)
    ex = np.array([1.0, 0, 0])
    atoms = []
    for rod in (np.zeros(3), np.array([0, b / 2, c / 2])):
        for k in range(2):
            atoms.append(("Cr1", "Cr", rod + 2 * half * k * ex))
        for k, sgn in ((0, 1), (1, -1)):        # carboxylate pairs at x = 1.7 and 5.1
            x = rod + (half + 2 * half * k) * ex
            u = np.array([0, np.cos(th), sgn * np.sin(th)])
            v = np.array([0, -sgn * np.sin(th), np.cos(th)])
            atoms.append(("O1", "O", x + np.sqrt(1.95 ** 2 - half ** 2) * v))   # mu-OH
            for s in (1, -1):
                uu = s * u
                atoms.append(("C1", "C", x + 2.45 * uu))
                for t in (1, -1):
                    atoms.append(("O2", "O", x + 1.12 * t * ex + 1.86 * uu))
                if s == 1:
                    atoms += ring(x + reach * uu, uu, ex, "")
    return (a, b, c), atoms


def cristobalite():
    """Ideal beta-cristobalite SiO2 (Fd-3m, linear Si-O-Si, Si-O = 1.61 A)."""
    a = 8 * 1.61 / np.sqrt(3)
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    bonds = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / 4.0
    atoms = []
    for t in fcc:
        atoms.append(("Si1", "Si", t * a))
        atoms.append(("Si1", "Si", (t + 0.25) * a))
        for b in bonds:
            atoms.append(("O1", "O", (t + b / 2) * a))
    return (a, a, a), atoms


def sodalite(a=8.83, d_sio=1.61):
    """Siliceous sodalite (SOD net). T atoms on the ideal net, O relaxed to regular
    SiO4 tetrahedra (Si-O = 1.61 A, O-O = 2.63 A) by least squares."""
    from scipy.optimize import least_squares
    base = [(.25, 0, .5), (.75, 0, .5), (.5, .25, 0), (.5, .75, 0), (0, .5, .25), (0, .5, .75)]
    T = np.array(base + [tuple((np.array(p) + .5) % 1) for p in base])
    L = a * np.eye(3)
    # T-T edges (3.12 A) and O start positions at the edge midpoints
    edges = []
    for i in range(len(T)):
        for j in range(len(T)):
            for n in np.array(np.meshgrid(*[[-1, 0, 1]] * 3)).reshape(3, -1).T:
                d = np.linalg.norm((T[j] + n - T[i]) @ L)
                if 0.1 < d < 3.4 and (i, j, tuple(n)) not in [(e[1], e[0], tuple(-np.array(e[2]))) for e in edges]:
                    edges.append((i, j, tuple(n)))
    x0 = np.array([((T[i] + T[j] + np.array(n)) / 2) @ L for i, j, n in edges])
    rng = np.random.default_rng(1)
    x0 += 0.15 * rng.standard_normal(x0.shape)          # break the planar start
    tet = {i: [] for i in range(len(T))}              # O around each T: (edge, sign of image)
    for k, (i, j, n) in enumerate(edges):
        tet[i].append((k, np.zeros(3)))
        tet[j].append((k, -np.array(n, float)))

    def res(x):
        X = x.reshape(-1, 3)
        r = []
        for t, lst in tet.items():
            P = [X[k] + s @ L for k, s in lst]
            r += [np.linalg.norm(p - T[t] @ L) - d_sio for p in P]
            r += [np.linalg.norm(P[p] - P[q]) - d_sio * np.sqrt(8 / 3)
                  for p in range(4) for q in range(p + 1, 4)]
        return r
    X = least_squares(res, x0.ravel()).x.reshape(-1, 3)
    atoms = [("Si1", "Si", t @ L) for t in T] + [("O1", "O", x) for x in X]
    return (a, a, a), atoms


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    write_cif(os.path.join(out, "IRMOF-1_ideal.cif"), "IRMOF1_ideal", *irmof1())
    write_cif(os.path.join(out, "DMOF_ideal.cif"), "DMOF_ideal", *dmof())
    write_cif(os.path.join(out, "MIL-53_ideal.cif"), "MIL53_ideal", *mil53())
    write_cif(os.path.join(out, "cristobalite_ideal.cif"), "beta_cristobalite_ideal", *cristobalite())
    write_cif(os.path.join(out, "SOD_silica_ideal.cif"), "SOD_silica_ideal", *sodalite())
    print("written to", out)
