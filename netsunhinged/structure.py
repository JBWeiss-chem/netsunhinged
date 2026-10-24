"""Periodic structure: bonding, framework detection, pruning, rings, chemistry rules."""

from __future__ import annotations

from collections import deque

import numpy as np
from scipy.spatial import cKDTree

from .cif import CifData, cell_matrix
from .elements import COVALENT_RADII, is_metal

ZERO = (0, 0, 0)


def add(n, m):
    return (n[0] + m[0], n[1] + m[1], n[2] + m[2])


def neg(n):
    return (-n[0], -n[1], -n[2])


def sub(n, m):
    return (n[0] - m[0], n[1] - m[1], n[2] - m[2])


def periodic_pairs(frac: np.ndarray, L: np.ndarray, cutoff: float):
    """All pairs (i, j, n) with |x_j + n L - x_i| < cutoff, each unordered pair once.

    Returns arrays i, j, n (k x 3 ints) and distances.
    """
    N = len(frac)
    if N == 0:
        return np.zeros(0, int), np.zeros(0, int), np.zeros((0, 3), int), np.zeros(0)
    V = abs(np.linalg.det(L))
    widths = np.array([V / np.linalg.norm(np.cross(L[(k + 1) % 3], L[(k + 2) % 3])) for k in range(3)])
    nmax = np.ceil(cutoff / widths).astype(int) + 1
    rng = [np.arange(-m, m + 1) for m in nmax]
    images = np.array(np.meshgrid(*rng, indexing="ij")).reshape(3, -1).T
    cart = frac @ L
    sup = (cart[None, :, :] + (images @ L)[:, None, :]).reshape(-1, 3)
    tree = cKDTree(sup)
    res = tree.query_ball_point(cart, cutoff)
    I, J, NN, D = [], [], [], []
    for i, hits in enumerate(res):
        for h in hits:
            img, j = divmod(h, N)
            n = images[img]
            if j < i:
                continue
            if j == i:
                nt = tuple(n)
                if nt == ZERO or nt < ZERO:  # keep one orientation of self-image pairs
                    continue
            d = np.linalg.norm(sup[h] - cart[i])
            I.append(i), J.append(j), NN.append(n), D.append(d)
    return (np.array(I, int), np.array(J, int),
            np.array(NN, int).reshape(-1, 3), np.array(D))


def _span_rank(vectors):
    if not vectors:
        return 0
    return int(np.linalg.matrix_rank(np.array(vectors, float), tol=1e-6))


def plane_deviation(points: np.ndarray) -> float:
    """Max distance of points from their best-fit plane."""
    if len(points) <= 3:
        return 0.0
    c = points - points.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return float(np.abs(c @ vt[-1]).max())


def line_deviation(points: np.ndarray) -> float:
    """RMS distance of points from their best-fit line."""
    if len(points) <= 2:
        return 0.0
    c = points - points.mean(0)
    s = np.linalg.svd(c, compute_uv=False)
    return float(np.sqrt((s[1:] ** 2).sum() / len(points)))


class Structure:
    """A periodic crystal structure with a bond graph and chemistry annotations."""

    def __init__(self, cif: CifData, bond_tol: float = 0.45, metal_metal_tol: float = 0.20,
                 max_ring: int = 8):
        self.name = cif.name
        self.cell = cif.cell
        self.L = cell_matrix(*cif.cell)
        self.bond_tol, self.metal_metal_tol, self.max_ring = bond_tol, metal_metal_tol, max_ring
        self.labels = [a.label for a in cif.atoms]
        self.elements = [a.element for a in cif.atoms]
        self.occupancy = [a.occupancy for a in cif.atoms]
        self.frac = np.array([a.frac for a in cif.atoms]).reshape(-1, 3)
        self.cart = self.frac @ self.L
        self.N = len(self.elements)
        self.warnings = list(cif.warnings)
        self.metal = np.array([is_metal(e) for e in self.elements], bool)
        self._find_bonds(bond_tol, metal_metal_tol)
        self._detach_extraframework_cations()
        self._find_components()
        self._prune()
        self._find_rings(max_ring)
        self._environments()
        self._find_sbus()
        self._auto_bond_states()

    # ---------------------------------------------------------- bonding
    def _find_bonds(self, tol, mm_tol):
        radii = np.array([COVALENT_RADII.get(e, 1.5) for e in self.elements])
        cutoff = 2 * radii.max() + tol if self.N else 0.0
        I, J, NN, D = periodic_pairs(self.frac, self.L, cutoff)
        self.bonds = []          # list of (i, j, n) : i bonded to j at image n
        self.bond_len = []
        overlaps = 0
        for i, j, n, d in zip(I, J, NN, D):
            ei, ej = self.elements[i], self.elements[j]
            if ei in ("H", "D") and ej in ("H", "D"):
                continue
            t = mm_tol if (self.metal[i] and self.metal[j]) else tol
            if d < 0.6:
                overlaps += 1
                continue
            if d < radii[i] + radii[j] + t:
                self.bonds.append((int(i), int(j), tuple(int(x) for x in n)))
                self.bond_len.append(float(d))
        if overlaps:
            self.warnings.append(f"{overlaps} atom pairs closer than 0.6 A were ignored "
                                 "(unresolved disorder?).")
        self.nbrs = [[] for _ in range(self.N)]   # site -> list of (j, n, bond_id)
        for b, (i, j, n) in enumerate(self.bonds):
            self.nbrs[i].append((j, n, b))
            self.nbrs[j].append((i, neg(n), b))

    EXTRAFRAMEWORK_CATIONS = {"Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba"}

    def _detach_extraframework_cations(self):
        """In tetrahedral (zeolite-type) frameworks, alkali / alkaline-earth cations sit in the
        pores and coordinate framework O; they are not part of the framework. Their bonds are
        removed so that they become guests. MOFs (no T-O-T bridges) are not affected."""
        cat = {s for s in range(self.N) if self.elements[s] in self.EXTRAFRAMEWORK_CATIONS}
        if not cat:
            return
        nb = [[j for j, _, _ in self.nbrs[s] if j not in cat] for s in range(self.N)]
        tcen = {s for s in range(self.N) if s not in cat and len(nb[s]) == 4
                and all(self.elements[j] == "O" for j in nb[s])}
        n_tot = sum(1 for s in range(self.N) if self.elements[s] == "O" and len(nb[s]) == 2
                    and all(j in tcen for j in nb[s]))
        if n_tot == 0:
            return
        keep = [k for k, (i, j, _) in enumerate(self.bonds) if i not in cat and j not in cat]
        self.bonds = [self.bonds[k] for k in keep]
        self.bond_len = [self.bond_len[k] for k in keep]
        self.nbrs = [[] for _ in range(self.N)]
        for b, (i, j, n) in enumerate(self.bonds):
            self.nbrs[i].append((j, n, b))
            self.nbrs[j].append((i, neg(n), b))
        els = sorted({self.elements[s] for s in cat})
        self.warnings.append(f"{len(cat)} extra-framework cations ({', '.join(els)}) in a tetrahedral "
                             "framework are treated as guests.")

    # -------------------------------------------------------- components
    def _find_components(self):
        self.component = -np.ones(self.N, int)
        self.unwrap = [ZERO] * self.N
        self.comp_dim = []
        cid = 0
        for s in range(self.N):
            if self.component[s] >= 0:
                continue
            self.component[s] = cid
            self.unwrap[s] = ZERO
            q = deque([s])
            per = []
            while q:
                i = q.popleft()
                for j, n, _ in self.nbrs[i]:
                    exp = add(self.unwrap[i], n)
                    if self.component[j] < 0:
                        self.component[j] = cid
                        self.unwrap[j] = exp
                        q.append(j)
                    elif exp != self.unwrap[j]:
                        per.append(sub(exp, self.unwrap[j]))
            self.comp_dim.append(_span_rank(per))
            cid += 1
        dims = np.array(self.comp_dim)
        self.framework = np.array([dims[c] > 0 for c in self.component], bool) if self.N else np.zeros(0, bool)
        if not self.framework.any():
            self.warnings.append("No periodic (extended) bonded network found - is this a "
                                 "molecular crystal? All atoms are treated as framework.")
            self.framework[:] = True
        fdims = {self.comp_dim[c] for c in set(self.component[self.framework])}
        self.dimensionality = max(fdims) if fdims else 0
        if self.dimensionality < 3:
            self.warnings.append(f"The bonded framework is {self.dimensionality}-periodic "
                                 "(chains/layers). The model has no constraints in directions without "
                                 "bonds and reports them as flexible, although interlayer or interchain "
                                 "interactions may hold them in place.")
        self.n_guest = int((~self.framework).sum())

    # ----------------------------------------------------------- pruning
    def _prune(self):
        """Remove dangling atoms (H, terminal groups): they add no constraints to the model.

        anchor[s] = (k, t): pruned atom s is carried by core atom k; s (image 0) sits next
        to k at image t, i.e. relative to k(image 0) s is at image -t.
        """
        alive = self.framework.copy()
        deg = np.zeros(self.N, int)
        for i, j, n in self.bonds:
            if alive[i] and alive[j]:
                deg[i] += 1
                deg[j] += 1
        parent = {}
        q = deque([s for s in range(self.N) if alive[s] and deg[s] <= 1])
        while q:
            s = q.popleft()
            if not alive[s] or deg[s] > 1:
                continue
            alive[s] = False
            for j, n, _ in self.nbrs[s]:
                if alive[j]:
                    parent[s] = (j, n)
                    deg[j] -= 1
                    if deg[j] <= 1:
                        q.append(j)
        self.core = alive
        self.anchor = {}
        for s in parent:
            k, t = s, ZERO
            while k in parent:
                j, n = parent[k]
                k, t = j, add(t, n)
            if self.core[k]:
                self.anchor[s] = (k, t)
        if self.framework.any() and not self.core.any():
            self.warnings.append("Pruning removed the whole framework (no closed loops).")
        self.core_idx = np.where(self.core)[0]

    def core_nbrs(self, s):
        return [(j, n, b) for j, n, b in self.nbrs[s] if self.core[j]]

    # -------------------------------------------------------------- rings
    def _find_rings(self, max_ring):
        """Smallest non-periodic ring through every core bond (size <= max_ring)."""
        rings = {}
        for b, (i, j, n) in enumerate(self.bonds):
            if not (self.core[i] and self.core[j]):
                continue
            # BFS from (j, n) back to (i, 0) without using bond b directly
            start, goal = (j, n), (i, ZERO)
            prev = {start: None}
            q = deque([(start, 1)])
            found = None
            while q and found is None:
                (a, na), depth = q.popleft()
                if depth >= max_ring:
                    continue
                for c, m, bb in self.nbrs[a]:
                    if not self.core[c]:
                        continue
                    node = (c, add(na, m))
                    if node == goal and not (a == j and na == n and bb == b):
                        prev[node] = (a, na)
                        found = node
                        break
                    if node not in prev and node != goal:
                        prev[node] = (a, na)
                        q.append((node, depth + 1))
            if found is None:
                continue
            path = []
            node = found
            while node is not None:
                path.append(node)
                node = prev[node]
            # canonical key: set of (site, image) shifted so the minimal site is at its min image
            key = self._canon(path)
            rings.setdefault(key, path)
        self.rings = []
        for path in rings.values():
            pts = np.array([(self.frac[s] + np.array(n)) @ self.L for s, n in path])
            organic = not any(self.metal[s] for s, _ in path)
            self.rings.append({
                "atoms": path,
                "planar": plane_deviation(pts) < 0.15,
                "organic": organic,
            })

    @staticmethod
    def _canon(nodes):
        nodes = list(nodes)
        smin = min(s for s, _ in nodes)
        best = None
        for s, n in nodes:
            if s != smin:
                continue
            key = tuple(sorted((a, sub(m, n)) for a, m in nodes))
            if best is None or key < best:
                best = key
        return best

    # ------------------------------------------------------- chemistry
    def _environments(self):
        self.env = []
        for s in range(self.N):
            nb = sorted(self.elements[j] for j, _, _ in self.nbrs[s])
            self.env.append(f"{self.elements[s]}({','.join(nb)})" if nb else self.elements[s])

    def _find_sbus(self):
        """Inorganic building units: metals plus their donor atoms, connected components."""
        sbu_atom = np.zeros(self.N, bool)
        for s in range(self.N):
            if not self.core[s]:
                continue
            if self.metal[s] or any(self.metal[j] for j, _, _ in self.core_nbrs(s)):
                sbu_atom[s] = True
        self.sbu_atom = sbu_atom
        self.sbu_id = -np.ones(self.N, int)
        sid = 0
        for s in np.where(sbu_atom)[0]:
            if self.sbu_id[s] >= 0:
                continue
            self.sbu_id[s] = sid
            q = deque([s])
            while q:
                i = q.popleft()
                for j, _, _ in self.core_nbrs(i):
                    if sbu_atom[j] and self.sbu_id[j] < 0 and (self.metal[i] or self.metal[j]):
                        self.sbu_id[j] = sid
                        q.append(j)
            sid += 1
        self.n_sbu = sid

        # donors of bridging groups (carboxylate-like): D bonded to metal, with a single
        # non-metal neighbour Z that carries >= 2 metal-bound donors.
        self.bridging_donor = np.zeros(self.N, bool)
        for s in range(self.N):
            if not self.core[s] or self.metal[s]:
                continue
            nb = self.core_nbrs(s)
            nm = [j for j, _, _ in nb if not self.metal[j]]
            if len(nm) != 1 or not any(self.metal[j] for j, _, _ in nb):
                continue
            z = nm[0]
            donors = [j for j, _, _ in self.core_nbrs(z)
                      if not self.metal[j] and any(self.metal[k] for k, _, _ in self.core_nbrs(j))]
            if len(donors) >= 2:
                self.bridging_donor[s] = True

        # zeolite-type frameworks: tetrahedral centres T (4 neighbours, all O) and the O
        # atoms bridging two of them (T-O-T)
        self.t_center = np.zeros(self.N, bool)
        for s in self.core_idx:
            nb = self.core_nbrs(s)
            self.t_center[s] = len(nb) == 4 and all(self.elements[j] == "O" for j, _, _ in nb)
        self.bridging_tot = np.zeros(self.N, bool)
        for s in self.core_idx:
            nb = self.core_nbrs(s)
            self.bridging_tot[s] = (self.elements[s] == "O" and len(nb) == 2
                                    and all(self.t_center[j] for j, _, _ in nb))

    def is_zeolite_like(self):
        """Corner-sharing tetrahedral oxide framework without organic carbon."""
        core_c = any(self.elements[s] == "C" for s in self.core_idx)
        return bool(self.bridging_tot.any()) and not core_c

    def _star_points(self, s, img=ZERO):
        pts = [(self.frac[s] + np.array(img)) @ self.L]
        for j, n, _ in self.core_nbrs(s):
            pts.append((self.frac[j] + np.array(add(img, n))) @ self.L)
        return np.array(pts)

    def _auto_bond_states(self):
        """Default torsional state of each core bond: 'hinge' (rotatable) or 'locked'."""
        in_planar_ring = set()
        for r in self.rings:
            if r["planar"] and r["organic"]:
                atoms = r["atoms"]
                for k in range(len(atoms)):
                    a, na = atoms[k]
                    c, nc = atoms[(k + 1) % len(atoms)]
                    for jj, n, b in self.nbrs[a]:
                        if jj == c and add(na, n) == nc:
                            in_planar_ring.add(b)
        self.auto_state = {}
        self.auto_reason = {}
        radii = COVALENT_RADII
        for b, (i, j, n) in enumerate(self.bonds):
            if not (self.core[i] and self.core[j]):
                continue
            if self.metal[i] or self.metal[j]:
                self.auto_state[b], self.auto_reason[b] = "hinge", "metal-ligand bond"
                continue
            if b in in_planar_ring:
                self.auto_state[b], self.auto_reason[b] = "locked", "planar ring"
                continue
            single = radii[self.elements[i]] + radii[self.elements[j]]
            short = self.bond_len[b] < single - 0.10
            pi_ = self._star_points(i)
            pj = self._star_points(j, n)
            linear = line_deviation(pi_) < 0.1 or line_deviation(pj) < 0.1
            planar = plane_deviation(np.vstack([pi_, pj])) < 0.2
            if short and (planar or linear):
                self.auto_state[b], self.auto_reason[b] = "locked", "multiple/conjugated bond"
            else:
                self.auto_state[b], self.auto_reason[b] = "hinge", "single bond"

    # ----------------------------------------------------------- topology
    def topology(self):
        """Reduced net: inorganic nodes and organic linkers with periodic edges.

        Units are lists of member atoms (site, image) in one consistent unwrapped frame, so a
        unit's position is the mean of its members for any (also deformed) geometry.
        Edges (u, v, t): unit v translated by lattice vector t is bonded to unit u.
        Infinite (periodic) SBUs such as rods are split into one node per metal."""
        unit_of = {}            # site -> (unit index, image of the site in the unit frame)
        units = []

        def grow(seed, allowed, via):
            """BFS over allowed sites; returns members or None if the component is periodic."""
            img = {seed: ZERO}
            q, periodic = deque([seed]), False
            while q:
                i = q.popleft()
                for j, n, _ in self.core_nbrs(i):
                    if not allowed(j) or not via(i, j):
                        continue
                    e = add(img[i], n)
                    if j not in img:
                        img[j] = e
                        q.append(j)
                    elif img[j] != e:
                        periodic = True
            return img, periodic

        # node atoms: metals and atoms bonded only to metals (mu-O, mu-OH ...); donor atoms
        # that also bond to the organic part belong to the linker (standard deconstruction)
        # tetrahedral centres (Si, Al, P ... in zeolites) are nodes as well
        center = self.metal | self.t_center
        node_atom = np.zeros(self.N, bool)
        for s in self.core_idx:
            nb = self.core_nbrs(s)
            node_atom[s] = center[s] or (len(nb) > 0 and all(center[j] for j, _, _ in nb))
        # inorganic nodes
        for s in self.core_idx:
            if not node_atom[s] or s in unit_of:
                continue
            img, periodic = grow(s, lambda j: node_atom[j], lambda i, j: True)
            if not periodic:
                k = len(units)
                units.append({"kind": "node", "members": [[int(a), list(n)] for a, n in img.items()]})
                for a, n in img.items():
                    unit_of[a] = (k, n)
                continue
            # rod / layer: one node per metal; donors join their first metal neighbour
            for m in img:
                if center[m]:
                    unit_of[m] = (len(units), ZERO)
                    units.append({"kind": "node", "members": [[int(m), [0, 0, 0]]]})
            for d in img:
                if center[d]:
                    continue
                for j, n, _ in self.core_nbrs(d):
                    if center[j] and j in unit_of:
                        unit_of[d] = (unit_of[j][0], neg(n))   # d relative to metal j
                        break
        # organic linkers
        for s in self.core_idx:
            if s in unit_of:
                continue
            img, _ = grow(s, lambda j: not node_atom[j], lambda i, j: True)
            k = len(units)
            units.append({"kind": "linker", "members": [[int(a), list(n)] for a, n in img.items()]})
            for a, n in img.items():
                unit_of[a] = (k, n)
        # edges between units
        edges = {}
        for i, j, n in self.bonds:
            if i not in unit_of or j not in unit_of:
                continue
            (u, ni), (v, nj) = unit_of[i], unit_of[j]
            t = sub(add(ni, n), nj)
            if u == v and t == ZERO:
                continue
            key = (u, v, t) if (u, v, t) <= (v, u, neg(t)) else (v, u, neg(t))
            edges[key] = True
        # each building block is drawn in the colour of a representative atom: one of its most
        # common non-H element (ties: metals first, then heavier elements); the front end
        # colours it through the current colour setting
        from collections import Counter
        from .elements import COVALENT_RADII
        for u in units:
            els = Counter(self.elements[a] for a, _ in u["members"] if self.elements[a] not in ("H", "D"))
            top = max(els.items(), key=lambda kv: (kv[1], is_metal(kv[0]), COVALENT_RADII.get(kv[0], 0)))[0]
            u["rep"] = next(int(a) for a, _ in u["members"] if self.elements[a] == top)
            u["label"] = "".join(f"{e}{c if c > 1 else ''}" for e, c in sorted(els.items()))
        return {"units": units, "edges": [[u, v, list(t)] for u, v, t in edges]}

    # --------------------------------------------------------------- misc
    def without_guests(self):
        """New structure containing only the bonded framework (guests/solvent removed)."""
        from .cif import CifAtom
        keep = np.where(self.framework)[0]
        atoms = [CifAtom(self.labels[s], self.elements[s], self.frac[s], self.occupancy[s])
                 for s in keep]
        warnings = []
        if any(a.occupancy < 0.99 for a in atoms):
            warnings.append("Structure contains partially occupied sites; they are kept as full atoms.")
        warnings.append(f"{self.N - len(keep)} guest atoms were removed.")
        cif = CifData(name=self.name, cell=self.cell, atoms=atoms, symops=["x,y,z"],
                      warnings=warnings)
        return Structure(cif, self.bond_tol, self.metal_metal_tol, self.max_ring)

    def formula(self, mask=None, reduced=False):
        """Hill-ordered formula of the cell (or of the atoms in mask).

        With reduced=True returns (formula divided by the common factor Z, Z)."""
        from collections import Counter
        from functools import reduce
        from math import gcd
        idx = range(self.N) if mask is None else np.where(mask)[0]
        c = Counter(self.elements[i] for i in idx)
        z = reduce(gcd, c.values(), 0) if reduced and c else 1
        order = sorted(c, key=lambda e: (e != "C", e != "H", e))
        text = " ".join(f"{e}{c[e] // z if c[e] // z > 1 else ''}" for e in order)
        return (text, z) if reduced else text
