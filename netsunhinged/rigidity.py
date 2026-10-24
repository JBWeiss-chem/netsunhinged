"""Infinitesimal rigidity of the periodic body-joint framework.

Unknowns: for every body a translation u (3) and rotation w (3), plus the symmetric cell
strain eps (6).  A joint constraint (atom shared by bodies 1 and 2) linearises to

    (p2 - p1) eps + u1 - u2 - [p1]x w1 + [p2]x w2 = 0

with p the atom position relative to the body centroid (body centroids follow the cell
affinely).  The null space of this compatibility matrix contains the mechanisms: deformations
that keep every body rigid and every joint closed to first order.  They cost nothing within the
model only; in the material, the released angles and torsions, non-bonded interactions and
guests still resist them.  Modes with a non-zero strain component change the unit
cell: they are the candidate large-amplitude ("breathing", "shear") flexibility modes
of the framework; strain-free modes are local mechanisms such as linker rotations.
"""

from __future__ import annotations

import numpy as np

from .model import MechanicalModel

STRAIN_NAMES = ["xx", "yy", "zz", "yz", "xz", "xy"]


def cell_params(L: np.ndarray):
    a, b, c = (np.linalg.norm(v) for v in L)
    al = np.degrees(np.arccos(np.dot(L[1], L[2]) / (b * c)))
    be = np.degrees(np.arccos(np.dot(L[0], L[2]) / (a * c)))
    ga = np.degrees(np.arccos(np.dot(L[0], L[1]) / (a * b)))
    return np.array([a, b, c, al, be, ga])


def strain_matrix(e6):
    xx, yy, zz, yz, xz, xy = e6
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def skew(v):
    """[v]x for an array of vectors (..., 3) -> (..., 3, 3)."""
    z = np.zeros(v.shape[:-1])
    return np.stack([
        np.stack([z, -v[..., 2], v[..., 1]], -1),
        np.stack([v[..., 2], z, -v[..., 0]], -1),
        np.stack([-v[..., 1], v[..., 0], z], -1)], -2)


_STRAIN_BASIS = []
for _i, _j in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]:
    _E = np.zeros((3, 3))
    _E[_i, _j] = _E[_j, _i] = 1.0
    _STRAIN_BASIS.append(_E)


def compatibility_matrix(m: MechanicalModel):
    """Sparse compatibility (rigidity) matrix, 3 rows per joint constraint."""
    from scipy.sparse import coo_matrix
    nb = len(m.bodies)
    ncol = 6 * nb + 6
    I, J, V = [], [], []
    r = 0

    def put(block, row0, col0):
        for a in range(block.shape[0]):
            for b in range(block.shape[1]):
                if block[a, b] != 0.0:
                    I.append(row0 + a), J.append(col0 + b), V.append(block[a, b])

    eye = np.eye(3)
    for b0, k0, b1, k1, dn, s in m.constraints:
        p0 = m.bodies[b0].pos[k0] - m.bodies[b0].centroid
        p1 = m.bodies[b1].pos[k1] - m.bodies[b1].centroid
        if b0 == b1:
            put(-skew(p0) + skew(p1), r, 6 * b0 + 3)
        else:
            put(eye, r, 6 * b0)
            put(-eye, r, 6 * b1)
            put(-skew(p0), r, 6 * b0 + 3)
            put(skew(p1), r, 6 * b1 + 3)
        d = p1 - p0
        put(np.stack([d @ E for E in _STRAIN_BASIS], 1), r, 6 * nb)
        r += 3
    for bi, B in enumerate(m.bodies):
        if B.linear:
            put(B.axis[None, :] * 3.0, r, 6 * bi + 3)
            r += 1
    return coo_matrix((V, (I, J)), shape=(r, ncol)).tocsr()


def mode_displacements(m: MechanicalModel, vec: np.ndarray):
    """Map a generalized velocity (bodies + strain) to per-site displacements and strain."""
    st = m.st
    nb = len(m.bodies)
    eps = strain_matrix(vec[6 * nb:6 * nb + 6])
    disp = np.zeros((st.N, 3))
    for s, lst in m.occurrences.items():
        b, k, n = lst[0]
        B = m.bodies[b]
        u, w = vec[6 * b:6 * b + 3], vec[6 * b + 3:6 * b + 6]
        p = B.pos[k] - B.centroid
        dx = B.centroid @ eps + u + np.cross(w, p)
        disp[s] = dx - (np.array(n) @ st.L) @ eps
    for s, (k, t) in st.anchor.items():
        lst = m.occurrences.get(k)
        if not lst:
            continue
        b, kk, n = lst[0]
        B = m.bodies[b]
        u, w = vec[6 * b:6 * b + 3], vec[6 * b + 3:6 * b + 6]
        img = np.array(n) - np.array(t)
        p = (st.frac[s] + img) @ st.L - B.centroid
        dx = B.centroid @ eps + u + np.cross(w, p)
        disp[s] = dx - (img @ st.L) @ eps
    # guests follow the cell affinely
    for s in np.where(~st.framework)[0]:
        disp[s] = st.cart[s] @ eps
    return disp, eps


def describe_strain(L, eps):
    """Relative cell-parameter changes for a small strain eps (3x3)."""
    h = 1e-6
    p0 = cell_params(L)
    p1 = cell_params(L + h * L @ eps)
    dp = (p1 - p0) / h
    rel = np.concatenate([dp[:3] / p0[:3], np.radians(dp[3:])])   # lengths relative, angles rad
    vol = float(np.trace(eps))
    return rel, vol


def classify_strain(rel, vol):
    names = ["a", "b", "c", "alpha", "beta", "gamma"]
    mag = np.abs(rel)
    tot = mag.max() if mag.max() > 0 else 1.0
    parts = [f"{names[i]} {'+' if rel[i] > 0 else '-'}{mag[i] / tot:.2f}"
             for i in np.argsort(-mag) if mag[i] / tot > 0.15]
    lengths = mag[:3].max() / tot
    angles = mag[3:].max() / tot
    if angles > 0.5 and lengths < 0.3:
        kind = "shear"
    elif abs(vol) / tot > 0.3 and np.all(np.sign(rel[:3][mag[:3] / tot > 0.3]) == np.sign(vol)) \
            and (mag[:3] / tot > 0.3).sum() == 3:
        kind = "isotropic-like expansion"
    elif lengths > 0.3 and np.ptp(np.sign(rel[:3][mag[:3] / tot > 0.3])) > 0:
        kind = "breathing (anisotropic, opposite length changes)"
    elif lengths > 0.3:
        kind = "uniaxial / Poisson-type"
    else:
        kind = "mixed"
    return kind, ", ".join(parts)


def analyze(m: MechanicalModel, tol: float = 1e-6, max_internal_modes: int = 12):
    A = compatibility_matrix(m)
    nb = len(m.bodies)
    ncol = 6 * nb + 6
    # column scaling: strain columns carry lengths of the order of the cell size
    scale = np.ones(ncol)
    Ls = np.linalg.norm(m.st.L, axis=1).mean()
    scale[6 * nb:] = 1.0 / Ls
    from scipy.sparse import diags
    As = A @ diags(scale)
    lam, V = np.linalg.eigh((As.T @ As).toarray())
    sv = np.sqrt(np.clip(lam, 0, None))
    smax = sv.max() if sv.size and sv.max() > 0 else 1.0
    null_mask = sv < tol * smax
    nullity = int(null_mask.sum())
    N = V[:, null_mask] * scale[:, None]           # back to physical variables
    # remove the 3 trivial translations
    T = np.zeros((ncol, 3))
    for bi in range(nb):
        T[6 * bi:6 * bi + 3] = np.eye(3)
    T, _ = np.linalg.qr(T)
    if N.shape[1]:
        N, _ = np.linalg.qr(N)            # orthonormal basis of the physical null space
        N = N - T @ (T.T @ N)
        # what survives the projection has singular value ~1, pure translations ~0; the
        # threshold must be absolute (a relative one keeps the translations when the
        # framework has no genuine mechanism)
        U, s2, _ = np.linalg.svd(N, full_matrices=False)
        N = U[:, s2 > 1e-6]
    n_mech = N.shape[1]

    # strain content of the mechanism space
    S = N[6 * nb:]
    macro = []
    if n_mech:
        # strain (x cell size) has the same units as the atom displacements
        Us, ss, Vs = np.linalg.svd(S * Ls, full_matrices=False)
        k_macro = int((ss > 1e-4).sum())
        for i in range(k_macro):
            vec = N @ Vs[i]
            e6 = vec[6 * nb:]
            norm = np.abs(e6).max()
            vec = vec / norm
            disp, eps = mode_displacements(m, vec)
            rel, vol = describe_strain(m.st.L, eps)
            kind, text = classify_strain(rel, vol)
            macro.append({"strain6": vec[6 * nb:].tolist(), "eps": eps.tolist(),
                          "cell_rel": rel.tolist(), "volume_rel": vol,
                          "kind": kind, "description": text,
                          "disp": disp.tolist(), "vec": vec})
        # internal (strain-free) mechanisms: orthogonal complement in the mechanism space
        P = np.eye(n_mech) - Vs[:k_macro].T @ Vs[:k_macro]
        Uq, sq, _ = np.linalg.svd(P)
        Vint = Uq[:, sq > 0.5]
        n_internal = Vint.shape[1]
    else:
        k_macro, n_internal, Vint = 0, 0, np.zeros((0, 0))

    internal = []
    for i in range(min(n_internal, max_internal_modes)):
        vec = N @ Vint[:, i]
        disp, _ = mode_displacements(m, vec)
        mx = np.linalg.norm(disp, axis=1).max()
        if mx < 1e-9:
            continue
        internal.append({"disp": (disp / mx).tolist()})
    for mode in macro:
        d = np.array(mode["disp"])
        mode["max_disp"] = float(np.linalg.norm(d, axis=1).max())

    nz = sv[~null_mask]
    return {
        "n_dof": ncol,
        "n_constraint_rows": int(A.shape[0]),
        "rank": int(ncol - nullity),
        "n_mechanisms": int(n_mech),
        "n_macroscopic": int(k_macro),
        "n_internal": int(n_internal),
        "smallest_nonzero_sv": [float(x / smax) for x in np.sort(nz)[:5]],
        "macroscopic_modes": macro,
        "internal_modes": internal,
    }
