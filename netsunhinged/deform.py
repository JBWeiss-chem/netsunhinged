"""Finite deformation test (Sarkisov, Martin, Haranczyk & Smit, JACS 2014, 136, 2228).

The unit cell is perturbed by an order parameter d (a* = d a, gamma* = d gamma, ...).
Every rigid body keeps its shape; its centroid (fractional coordinates) and orientation
are relaxed, together with the cell parameters that the regime allows to relax, by
minimising the joint mismatch energy

    E = (1/N_c) sum_c |x_c(body 1) - x_c(body 2)|^2      [A^2 per joint constraint]

which plays the role of the harmonic bond energy of the paper (Eq. 4 of the SI).  E is a
geometric joint mismatch, not a physical energy of the material.  A framework is classified as
flexible with respect to a regime if the mismatch stays ~0 over the whole range of d, i.e. if the
rigid bodies can follow the cell change without opening any joint.

The classification uses the relative mismatch F = E / E_affine, which is scale-free.  Absolute E
grows with (d - 1)^2 and with the distance between the joints, so a fixed threshold on E means
something different for every structure and every range of d (and the paper's thresholds refer to
its own energy function, not to A^2).  E is reported alongside for plotting.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import spsolve
from scipy.spatial.transform import Rotation

from .cif import cell_matrix
from .model import MechanicalModel
from .rigidity import cell_params, skew

# Table S1 of the SI: P = perturbed, R = relaxed, fixed otherwise.
PARAMS = ["a", "b", "c", "alpha", "beta", "gamma"]
REGIMES = {
    "ABC":   ("PPP---", "Linear compression/expansion"),
    "A":     ("P-----", "Young (a)"),
    "B":     ("-P----", "Young (b)"),
    "C":     ("--P---", "Young (c)"),
    "Alpha": ("---P--", "Shear (alpha)"),
    "Beta":  ("----P-", "Shear (beta)"),
    "Gamma": ("-----P", "Shear (gamma)"),
    "AB":    ("PR----", "Poisson: a perturbed, b relaxed"),
    "AC":    ("P-R---", "Poisson: a perturbed, c relaxed"),
    "BC":    ("-PR---", "Poisson: b perturbed, c relaxed"),
    "AB-C":  ("PPR---", "a,b perturbed, c relaxed"),
    "AC-B":  ("PRP---", "a,c perturbed, b relaxed"),
    "BC-A":  ("RPP---", "b,c perturbed, a relaxed"),
}
EXTENDED = {
    "A*":     ("PRRRRR", "a perturbed, all others relaxed"),
    "B*":     ("RPRRRR", "b perturbed, all others relaxed"),
    "C*":     ("RRPRRR", "c perturbed, all others relaxed"),
    "Alpha*": ("RRRPRR", "alpha perturbed, all others relaxed"),
    "Beta*":  ("RRRRPR", "beta perturbed, all others relaxed"),
    "Gamma*": ("RRRRRP", "gamma perturbed, all others relaxed"),
}


def _right_jacobian(theta):
    """Right Jacobian of SO(3) for rotation vectors (n x 3) -> (n x 3 x 3)."""
    t = np.linalg.norm(theta, axis=1)
    K = skew(theta)
    K2 = K @ K
    small = t < 1e-6
    ts = np.where(small, 1.0, t)
    c1 = np.where(small, 0.5, (1 - np.cos(ts)) / ts ** 2)
    c2 = np.where(small, 1.0 / 6.0, (ts - np.sin(ts)) / ts ** 3)
    return np.eye(3)[None] - c1[:, None, None] * K + c2[:, None, None] * K2


class DeformationProblem:
    def __init__(self, m: MechanicalModel):
        self.m = m
        st = m.st
        self.L0 = st.L
        self.p0 = np.array(st.cell, float)
        nb = len(m.bodies)
        self.nb = nb
        self.s0 = np.array([B.centroid for B in m.bodies]) @ np.linalg.inv(self.L0) if nb else np.zeros((0, 3))
        C = m.constraints
        self.nc = len(C)
        self.cb0 = np.array([c[0] for c in C], int)
        self.cb1 = np.array([c[2] for c in C], int)
        self.cp0 = np.array([m.bodies[c[0]].pos[c[1]] - m.bodies[c[0]].centroid for c in C]).reshape(-1, 3)
        self.cp1 = np.array([m.bodies[c[2]].pos[c[3]] - m.bodies[c[2]].centroid for c in C]).reshape(-1, 3)
        self.cdn = np.array([c[4] for c in C], float).reshape(-1, 3)
        self.lin = np.array([bi for bi, B in enumerate(m.bodies) if B.linear], int)
        self.lin_axis = np.array([m.bodies[bi].axis for bi in self.lin]).reshape(-1, 3)
        self.w_pin = 1.0
        self.w_gauge = 1e-3

    # x = [s (3nb), theta (3nb), relaxed cell params (scaled factors)]
    def unpack(self, x, relax_idx, cellp):
        nb = self.nb
        s = x[:3 * nb].reshape(nb, 3)
        th = x[3 * nb:6 * nb].reshape(nb, 3)
        p = cellp.copy()
        if len(relax_idx):
            p[relax_idx] = self.p0[relax_idx] * x[6 * nb:]
        return s, th, p

    def residuals(self, x, relax_idx, cellp):
        s, th, p = self.unpack(x, relax_idx, cellp)
        L = cell_matrix(*p)
        R = Rotation.from_rotvec(th).as_matrix() if self.nb else np.zeros((0, 3, 3))
        x0 = s[self.cb0] @ L + np.einsum("nij,nj->ni", R[self.cb0], self.cp0)
        x1 = s[self.cb1] @ L + np.einsum("nij,nj->ni", R[self.cb1], self.cp1)
        r = (x0 - x1 - self.cdn @ L).ravel()
        pins = self.w_pin * np.einsum("ni,ni->n", th[self.lin], self.lin_axis) if len(self.lin) else np.zeros(0)
        gauge = self.w_gauge * (s - self.s0).mean(0) * np.linalg.norm(self.L0, axis=1)
        return np.concatenate([r, pins, gauge])

    def jacobian(self, x, relax_idx, cellp):
        nb, nc = self.nb, self.nc
        s, th, p = self.unpack(x, relax_idx, cellp)
        L = cell_matrix(*p)
        R = Rotation.from_rotvec(th).as_matrix()
        Jr = _right_jacobian(th)
        rows, cols, vals = [], [], []
        ridx = np.arange(nc)[:, None] * 3 + np.arange(3)[None, :]      # nc x 3

        def add_block(blocks, col0):
            # blocks: nc x 3 x 3 (residual comp i, variable comp j), columns col0 (nc) + j
            rr = np.repeat(ridx[:, :, None], 3, axis=2)
            cc = np.repeat((col0[:, None] + np.arange(3)[None, :])[:, None, :], 3, axis=1)
            rows.append(rr.ravel()), cols.append(cc.ravel()), vals.append(blocks.ravel())

        LT = np.broadcast_to(L.T, (nc, 3, 3))
        add_block(LT, 3 * self.cb0)
        add_block(-LT, 3 * self.cb1)
        D0 = -np.einsum("nij,njk,nkl->nil", R[self.cb0], skew(self.cp0), Jr[self.cb0])
        D1 = -np.einsum("nij,njk,nkl->nil", R[self.cb1], skew(self.cp1), Jr[self.cb1])
        add_block(D0, 3 * nb + 3 * self.cb0)
        add_block(-D1, 3 * nb + 3 * self.cb1)
        nrow = 3 * nc
        # relaxed cell parameters (numerical derivative of the cell matrix)
        for q, pi in enumerate(relax_idx):
            h = 1e-6
            pp, pm = p.copy(), p.copy()
            pp[pi] += h * self.p0[pi]
            pm[pi] -= h * self.p0[pi]
            dL = (cell_matrix(*pp) - cell_matrix(*pm)) / (2 * h)
            d = (s[self.cb0] @ dL - s[self.cb1] @ dL - self.cdn @ dL).ravel()
            rows.append(np.arange(nrow)), cols.append(np.full(nrow, 6 * nb + q)), vals.append(d)
        # spin pins of linear bodies
        if len(self.lin):
            r0 = nrow + np.arange(len(self.lin))
            for k in range(3):
                rows.append(r0), cols.append(3 * nb + 3 * self.lin + k)
                vals.append(self.w_pin * self.lin_axis[:, k])
            nrow += len(self.lin)
        # gauge
        Ln = np.linalg.norm(self.L0, axis=1)
        for k in range(3):
            rows.append(np.full(nb, nrow + k)), cols.append(3 * np.arange(nb) + k)
            vals.append(np.full(nb, self.w_gauge * Ln[k] / nb))
        nrow += 3
        ncol = 6 * nb + len(relax_idx)
        return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(nrow, ncol)).tocsr()

    def energy(self, x, relax_idx, cellp):
        r = self.residuals(x, relax_idx, cellp)[:3 * self.nc]
        return float((r ** 2).sum() / max(self.nc, 1))

    def solve(self, x0, relax_idx, cellp, max_iter=300):
        """Sparse Levenberg-Marquardt.

        The problem is rank deficient (internal mechanisms, gauge), which the damping
        handles; the normal equations are solved with a sparse direct factorisation."""
        if self.nb == 0:
            return x0, 0.0
        x = x0.copy()
        r = self.residuals(x, relax_idx, cellp)
        f = r @ r
        lam = 1e-3
        for _ in range(max_iter):
            J = self.jacobian(x, relax_idx, cellp)
            A = (J.T @ J).tocsc()
            g = J.T @ r
            if np.abs(g).max() < 1e-13:
                break
            dA = A.diagonal()
            scale = np.maximum(dA, 1e-6 * max(dA.max(), 1e-12))
            # the damped normal matrix is symmetric positive definite: up to a few thousand
            # unknowns a dense Cholesky factorisation is several times faster than sparse LU
            Ad = A.toarray() if A.shape[0] <= 3000 else None
            improved = False
            for _ in range(30):
                try:
                    if Ad is not None:
                        step = cho_solve(cho_factor(Ad + np.diag(lam * scale)), -g)
                    else:
                        step = spsolve(A + diags(lam * scale), -g)
                except Exception:
                    lam *= 10
                    continue
                xn = x + step
                rn = self.residuals(xn, relax_idx, cellp)
                fn = rn @ rn
                if np.isfinite(fn) and fn < f:
                    improved = True
                    break
                lam *= 10
            if not improved:
                break
            rel = (f - fn) / max(f, 1e-30)
            x, r, f = xn, rn, fn
            lam = max(lam / 10, 1e-12)
            if rel < 1e-10 or f < 1e-24:
                break
        return x, self.energy(x, relax_idx, cellp)

    def x_init(self, n_relax):
        return np.concatenate([self.s0.ravel(), np.zeros(3 * self.nb), np.ones(n_relax)])

    # ---------------------------------------------------- run one regime
    def run_regime(self, pattern, dvals, progress=None, strain=None):
        """Scan the order parameter.  Returns list of dicts (d, E, cell, x).

        With ``strain`` (3x3) the cell is L0 (I + (d-1) strain) instead of a cell-parameter
        pattern (used to follow an infinitesimal mode to finite amplitude)."""
        if strain is None:
            pert = [i for i, ch in enumerate(pattern) if ch == "P"]
            relax = np.array([i for i, ch in enumerate(pattern) if ch == "R"], int)
        else:
            pert, relax = [], np.zeros(0, int)
        dvals = np.asarray(dvals)
        out = {}
        x1 = self.x_init(len(relax))
        out[1.0] = {"d": 1.0, "E": 0.0, "E_affine": 0.0, "F": 0.0, "cell": self.p0.tolist(), "x": x1}
        # a tiny, reproducible kick of the body orientations: high-symmetry structures (e.g.
        # linear T-O-T in ideal cristobalite) otherwise sit on a saddle point from which the
        # minimiser cannot find symmetry-breaking mechanisms such as tetrahedral tilts
        kick = np.zeros_like(x1)
        kick[3 * self.nb:6 * self.nb] = 1e-3 * np.random.default_rng(12345).standard_normal(3 * self.nb)
        for branch in (dvals[dvals > 1 + 1e-9], dvals[dvals < 1 - 1e-9][::-1]):
            x = x1 + kick
            for d in branch:
                if strain is None:
                    cellp = self.p0.copy()
                    cellp[pert] = self.p0[pert] * d
                else:
                    Ld = self.L0 + (d - 1.0) * self.L0 @ strain
                    cellp = cell_params(Ld)
                # energy without any relaxation (bodies follow the cell affinely)
                Ea = self.energy(x1, relax, cellp)
                x, E = self.solve(x, relax, cellp)
                _, _, p = self.unpack(x, relax, cellp)
                F = E / Ea if Ea > 1e-10 else 0.0
                out[float(d)] = {"d": float(d), "E": E, "E_affine": Ea, "F": F,
                                 "cell": p.tolist(), "x": x}
                if progress:
                    progress()
        return [out[k] for k in sorted(out)]

    # ------------------------------------------------- coordinates
    def positions(self, x, cellp, relax_idx=()):
        """Cartesian positions of all sites for solution x (home-cell images)."""
        m, st = self.m, self.m.st
        s, th, p = self.unpack(x, np.asarray(relax_idx, int), np.asarray(cellp, float))
        L = cell_matrix(*p)
        R = Rotation.from_rotvec(th).as_matrix()
        X = np.zeros((st.N, 3))
        for site, lst in m.occurrences.items():
            b, k, n = lst[0]
            B = m.bodies[b]
            X[site] = s[b] @ L + R[b] @ (B.pos[k] - B.centroid) - np.array(n) @ L
        for site, (k, t) in st.anchor.items():
            lst = m.occurrences.get(k)
            if not lst:
                continue
            b, _, n = lst[0]
            B = m.bodies[b]
            img = np.array(n) - np.array(t)
            pq = (st.frac[site] + img) @ st.L - B.centroid
            X[site] = s[b] @ L + R[b] @ pq - img @ L
        for site in np.where(~st.framework)[0]:
            X[site] = st.frac[site] @ L
        return X, L


def classify(E, d, thr=0.01, thr2=0.05):
    """Paper criteria: flexible / semi-flexible (limited range or low) / rigid.

    ``E`` is the relative mismatch F = E / E_affine along the scan (thresholds 0.01 / 0.05).
    The undeformed point d = 1 (F = 0 by definition) is where the interval search starts."""
    E = np.asarray(E)
    d = np.asarray(d)
    if np.all(E < thr):
        return "flexible", "mismatch below threshold over the whole range"
    below_lo = d < 1 - 1e-9
    below_hi = d > 1 + 1e-9
    if below_lo.any() and np.all(E[below_lo] < thr):
        return "semi-flexible", "mismatch below threshold for d < 1 only (one-sided)"
    if below_hi.any() and np.all(E[below_hi] < thr):
        return "semi-flexible", "mismatch below threshold for d > 1 only (one-sided)"
    i1 = int(np.argmin(np.abs(d - 1)))
    lo = hi = i1
    while lo - 1 >= 0 and E[lo - 1] < thr:
        lo -= 1
    while hi + 1 < len(E) and E[hi + 1] < thr:
        hi += 1
    if d[hi] - d[lo] >= 0.039:
        return "semi-flexible", f"mismatch below threshold only for {d[lo]:.2f} <= d <= {d[hi]:.2f}"
    if np.all(E < thr2):
        return "semi-flexible", "low but non-zero mismatch over the whole range"
    return "rigid", ""


def classify_regime(points, thr=0.01, thr2=0.05):
    """(classification, detail) from F (the web app repeats this in app.js to reclassify)."""
    return classify([p["F"] for p in points], [p["d"] for p in points], thr, thr2)


def mode_strain(mode, L):
    """Strain of an infinitesimal mode scaled so the largest relative cell change is 1."""
    from .rigidity import describe_strain
    eps = np.array(mode["eps"])
    rel, _ = describe_strain(L, eps)
    return eps / np.abs(rel).max()


def run_all(m: MechanicalModel, regimes=None, dmin=0.90, dmax=1.10, step=0.02,
            thr=0.01, thr2=0.05, modes=(), progress=None):
    """Run the deformation regimes (and optionally follow infinitesimal ``modes``).

    Returns (problem, results, verdict, seconds); results maps regime name to
    {pattern, description, points, classification, detail}."""
    prob = DeformationProblem(m)
    all_reg = {**REGIMES, **EXTENDED}
    regimes = list(REGIMES) if regimes is None else list(regimes)
    n = int(round((dmax - dmin) / step))
    dvals = np.round(np.linspace(dmin, dmax, n + 1), 6)
    total = (len(regimes) + len(modes)) * (len(dvals) - 1)
    done = [0]

    def tick():
        done[0] += 1
        if progress:
            progress(done[0] / max(total, 1))

    results = {}
    t0 = time.time()
    for name in regimes:
        pattern, desc = all_reg[name]
        pts = prob.run_regime(pattern, dvals, tick)
        cls, detail = classify_regime(pts, thr, thr2)
        results[name] = {"pattern": pattern, "description": desc, "points": pts,
                         "classification": cls, "detail": detail}
    for k, mode in enumerate(modes):
        eps = mode_strain(mode, m.st.L)
        pts = prob.run_regime(None, dvals, tick, strain=eps)
        cls, detail = classify_regime(pts, thr, thr2)
        results[f"Mode {k + 1}"] = {"pattern": "------", "strain": eps.tolist(),
                                    "description": f"follow infinitesimal mode {k + 1}: "
                                                   f"{mode['description']}",
                                    "points": pts, "classification": cls, "detail": detail}
    verdict = overall(results)
    return prob, results, verdict, time.time() - t0


def overall(results):
    cls = [r["classification"] for r in results.values()]
    if "flexible" in cls:
        return "flexible"
    if "semi-flexible" in cls:
        return "semi-flexible"
    return "rigid"
