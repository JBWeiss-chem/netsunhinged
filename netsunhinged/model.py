"""Mechanical (body-joint) model of a framework.

The framework is described by *rigid units*: sets of atom copies (site, image) whose
mutual distances are fixed.  The default (bond-bending) model fixes all bond lengths
and bond angles:

* every core atom X contributes its *star* {X} + neighbours (fixed angles at X);
* an atom whose bond angles are declared **flexible** contributes only its bonds
  (it becomes a pivot / ball joint);
* every bond is a **hinge** (free torsion) unless it is **locked**, in which case the
  two stars are merged (fixed torsion);
* optionally, every inorganic SBU is made completely rigid (Sarkisov et al.: all sites of
  an SBU within a cutoff are connected by rigid bars).

Units sharing >= 3 non-collinear atoms are fused into rigid *bodies*.  Bodies sharing two
atoms are joined by a hinge, bodies sharing one atom by a ball joint.  All shared atoms
become joint constraints of the periodic body-joint framework analysed in
:mod:`netsunhinged.rigidity` and :mod:`netsunhinged.deform`.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field

import numpy as np

from .structure import ZERO, Structure, add, line_deviation, neg, sub


@dataclass
class ModelOptions:
    flex_angle_labels: set = field(default_factory=set)    # CIF labels with flexible angles
    flex_angle_sites: set = field(default_factory=set)     # individual sites
    bond_overrides: dict = field(default_factory=dict)     # "lab1|lab2" -> 'hinge' | 'locked'
    rigid_sbu: bool = False
    sbu_cutoff: float = 5.0

    @staticmethod
    def bond_key(l1, l2):
        return "|".join(sorted((l1, l2)))

    @classmethod
    def from_dict(cls, d):
        """Options from JSON (app export or API request); malformed input raises ValueError."""
        try:
            labels = set(str(x) for x in d.get("flex_angle_labels", []))
            sites = set(int(x) for x in d.get("flex_angle_sites", []))
            overrides = {str(k): str(v) for k, v in dict(d.get("bond_overrides", {})).items()}
            cutoff = float(d.get("sbu_cutoff", 5.0))
        except (TypeError, ValueError):
            raise ValueError("malformed model options") from None
        if any(v not in ("hinge", "locked") for v in overrides.values()):
            raise ValueError("bond overrides must be 'hinge' or 'locked'")
        if not 0.0 <= cutoff <= 20.0:
            raise ValueError("the SBU cutoff must lie between 0 and 20 Å")
        return cls(flex_angle_labels=labels, flex_angle_sites=sites, bond_overrides=overrides,
                   rigid_sbu=bool(d.get("rigid_sbu", False)), sbu_cutoff=cutoff)

    def to_dict(self):
        return {
            "flex_angle_labels": sorted(self.flex_angle_labels),
            "flex_angle_sites": sorted(self.flex_angle_sites),
            "bond_overrides": self.bond_overrides,
            "rigid_sbu": self.rigid_sbu,
            "sbu_cutoff": self.sbu_cutoff,
        }


PRESET_LABELS = {
    "paper": "Carboxylate hinges (Sarkisov–Smit)",
    "zeolite": "Rigid tetrahedra (T–O–T pivots)",
    "rigid": "All angles rigid",
}


def auto_preset(st: Structure):
    """Starting model chosen from structural motifs (not from a material class).

    Returns (preset name, human-readable reason)."""
    n_carb = int(st.bridging_donor.sum())
    n_tot = int(st.bridging_tot.sum())
    if st.is_zeolite_like():
        return "zeolite", (f"{n_tot} bridging T–O–T oxygens of corner-sharing tetrahedra found: "
                           "TO4 tetrahedra rigid, T–O–T angles flexible.")
    if n_carb:
        return "paper", (f"{n_carb} carboxylate-type oxygens found: rigid SBUs, flexible M–O–C "
                         "angles (the paper's O···O hinge).")
    return "rigid", ("No carboxylate or T–O–T motif found, so nothing is made flexible "
                     "automatically. For azolate frameworks (ZIFs) the soft DOFs are usually "
                     "the angles at the metal (e.g. Zn) or at the donor N: tick them in the table.")


def preset(st: Structure, name: str) -> ModelOptions:
    """Predefined DOF models.

    * ``paper``     Sarkisov, Martin, Haranczyk, Smit (JACS 2014): rigid SBUs, rigid
                    linkers, hinges where carboxylate-like groups bind the SBU (flexible
                    M-O-C angle -> rotation about the O...O axis) and free axial rotation of
                    single bonds such as C(carboxylate)-C(ring).
    * ``zeolite``   rigid TO4 tetrahedra, flexible T-O-T angles at bridging O (rigid-unit model).
    * ``rigid``     all bond angles fixed (rigid SBUs if the framework has metals).
    * ``bondbend``  all bond lengths and angles fixed, torsions about single bonds free.
    * ``donor_pivots``  like ``paper`` but also flexible angles at every donor atom.
    * ``auto``      chosen from the structural motifs, see :func:`auto_preset`.
    """
    if name == "bondbend":
        return ModelOptions()
    if name == "rigid":
        return ModelOptions(rigid_sbu=bool((st.metal & st.core).any()))
    if name == "paper":
        labels = {st.labels[s] for s in np.where(st.bridging_donor)[0]}
        return ModelOptions(flex_angle_labels=labels, rigid_sbu=True)
    if name == "donor_pivots":
        donors = {st.labels[s] for s in np.where(st.sbu_atom & ~st.metal)[0]}
        return ModelOptions(flex_angle_labels=donors, rigid_sbu=True)
    if name == "zeolite":
        # rigid TO4 tetrahedra (angles at T fixed) joined at flexible T-O-T angles: the
        # rigid-unit model of zeolites / silicates
        bridges = {st.labels[s] for s in np.where(st.bridging_tot)[0]}
        return ModelOptions(flex_angle_labels=bridges, rigid_sbu=False)
    if name == "auto":
        return preset(st, auto_preset(st)[0])
    raise ValueError(f"unknown preset {name}")


@dataclass
class Body:
    atoms: list            # list of (site, image)
    pos: np.ndarray        # unwrapped cartesian positions (k x 3)
    centroid: np.ndarray
    linear: bool
    axis: np.ndarray       # principal axis (used if linear)
    kind: str              # 'sbu', 'linker', 'mixed'


def _canon(unit):
    """Canonical form of a unit modulo lattice translation."""
    smin = min(s for s, _ in unit)
    best = None
    for s, n in unit:
        if s != smin:
            continue
        key = tuple(sorted((a, sub(m, n)) for a, m in unit))
        if best is None or key < best:
            best = key
    return best


class MechanicalModel:
    def __init__(self, st: Structure, opts: ModelOptions):
        self.st = st
        self.opts = opts
        self.flex = np.zeros(st.N, bool)
        for s in st.core_idx:
            if st.labels[s] in opts.flex_angle_labels or s in opts.flex_angle_sites:
                self.flex[s] = True
        self.bond_state = {}
        for b, auto in st.auto_state.items():
            i, j, _ = st.bonds[b]
            self.bond_state[b] = opts.bond_overrides.get(
                ModelOptions.bond_key(st.labels[i], st.labels[j]), auto)
        self._build_units()
        self._fuse()
        self._build_constraints()
        self._find_joints()

    # ----------------------------------------------------------- units
    def _star(self, s, img=ZERO):
        st = self.st
        u = [(s, img)]
        if not self.flex[s]:
            u += [(j, add(img, n)) for j, n, _ in st.core_nbrs(s)]
        return u

    def _xyz(self, s, n):
        return (self.st.frac[s] + np.array(n)) @ self.st.L

    def _build_units(self):
        st = self.st
        self.n_sbu_released = 0
        self.released_sbu = set()
        units = []
        for s in st.core_idx:
            if self.flex[s]:
                for j, n, _ in st.core_nbrs(s):
                    units.append([(s, ZERO), (j, n)])
            else:
                units.append(self._star(s))
        for b, state in self.bond_state.items():
            if state != "locked":
                continue
            i, j, n = st.bonds[b]
            units.append(self._star(i) + self._star(j, n))
        if self.opts.rigid_sbu:
            units += self._sbu_units()
        # dedupe modulo translation
        seen = {}
        for u in units:
            u = list(dict.fromkeys(u))
            if len(u) < 2:
                continue
            seen.setdefault(_canon(u), u)
        self.units = list(seen.values())

    def _sbu_units(self):
        """For every SBU atom X: all atoms of the same SBU within cutoff (bonded path)."""
        st = self.st
        cut = self.opts.sbu_cutoff
        # an explicitly flexible metal frees its coordination sphere: such SBUs are not
        # rigidified (otherwise the rigid-SBU rule would silently undo the user's choice)
        released = {st.sbu_id[s] for s in st.core_idx if st.metal[s] and self.flex[s]}
        self.released_sbu = released
        self.n_sbu_released = len(released)
        units = []
        for s in st.core_idx:
            if not st.sbu_atom[s] or st.sbu_id[s] in released:
                continue
            x0 = self._xyz(s, ZERO)
            seen = {(s, ZERO)}
            q = deque([(s, ZERO, 0)])
            while q:
                a, na, depth = q.popleft()
                if depth >= 8:
                    continue
                for j, n, _ in st.core_nbrs(a):
                    if not st.sbu_atom[j] or st.sbu_id[j] != st.sbu_id[s]:
                        continue
                    if not (st.metal[a] or st.metal[j]):
                        continue
                    node = (j, add(na, n))
                    if node in seen:
                        continue
                    if np.linalg.norm(self._xyz(*node) - x0) > cut + 2.5:
                        continue
                    seen.add(node)
                    q.append((node[0], node[1], depth + 1))
            unit = [nd for nd in seen if np.linalg.norm(self._xyz(*nd) - x0) <= cut]
            if len(unit) >= 2:
                units.append(unit)
        return units

    # ----------------------------------------------------------- fusion
    def _noncollinear(self, nodes):
        if len(nodes) < 3:
            return False
        pts = np.array([self._xyz(s, n) for s, n in nodes])
        return line_deviation(pts) > 0.12

    def _fuse(self):
        """Union rigid units sharing >= 3 non-collinear atoms into bodies."""
        bodies = [list(u) for u in self.units]
        for _ in range(50):
            index = defaultdict(list)   # site -> [(body, image)]
            for bi, B in enumerate(bodies):
                for s, n in B:
                    index[s].append((bi, n))
            parent = list(range(len(bodies)))
            offset = [ZERO] * len(bodies)   # body shifted by offset lives in root frame

            def find(x):
                path = []
                while parent[x] != x:
                    path.append(x)
                    x = parent[x]
                root = x
                # path compression with accumulated offsets
                for p in reversed(path):
                    q = parent[p]
                    if q != root:
                        offset[p] = add(offset[p], offset[q])
                    parent[p] = root
                return root

            def off_to_root(x):
                find(x)
                return offset[x] if parent[x] != x else ZERO

            merged = False
            for ai, A in enumerate(bodies):
                shared = defaultdict(list)
                for s, n in A:
                    for bi, m in index[s]:
                        if bi == ai:
                            continue
                        shared[(bi, sub(n, m))].append((s, n))
                for (bi, t), nodes in shared.items():
                    if len(nodes) < 3 or not self._noncollinear(nodes):
                        continue
                    ra, rb = find(ai), find(bi)
                    oa, ob = off_to_root(ai), off_to_root(bi)
                    if ra == rb:
                        continue   # already fused (possibly with a translate: periodic body)
                    # B shifted by t coincides with A; in root frames:
                    # rb-frame shifted by (oa + t - ob) lies in ra-frame
                    parent[rb] = ra
                    offset[rb] = sub(add(oa, t), ob)
                    merged = True
            if not merged:
                break
            groups = defaultdict(list)
            for bi in range(len(bodies)):
                groups[find(bi)].append(bi)
            new = []
            for root, members in groups.items():
                nodes = {}
                for bi in members:
                    o = off_to_root(bi)
                    for s, n in bodies[bi]:
                        nodes[(s, add(n, o))] = True
                new.append(list(nodes))
            # drop bodies contained in others (modulo translation)
            bodies = self._drop_subsets(new)
        self.bodies = []
        for B in self._drop_subsets(bodies):
            pos = np.array([self._xyz(s, n) for s, n in B])
            c = pos.mean(0)
            lin = line_deviation(pos) < 0.1
            axis = np.zeros(3)
            if len(pos) >= 2:
                _, _, vt = np.linalg.svd(pos - c, full_matrices=False)
                axis = vt[0]
            sbu = [self.st.sbu_atom[s] for s, _ in B]
            kind = "sbu" if all(sbu) else ("mixed" if any(sbu) and any(self.st.metal[s] for s, _ in B)
                                          else "linker")
            self.bodies.append(Body(B, pos, c, lin, axis, kind))

    @staticmethod
    def _drop_subsets(bodies):
        bodies = sorted(bodies, key=len, reverse=True)
        index = defaultdict(list)    # site -> [(kept body, image)]
        keep, keep_sets = [], []
        for B in bodies:
            s0, n0 = B[0]
            contained = any(all((s, add(n, sub(m, n0))) in keep_sets[ki] for s, n in B)
                            for ki, m in index[s0])
            if contained:
                continue
            keep.append(B)
            keep_sets.append(set(B))
            for s, n in B:
                index[s].append((len(keep) - 1, n))
        return keep

    # ------------------------------------------------------- constraints
    def _build_constraints(self):
        """Every atom copy shared by several bodies gives a 3D coincidence constraint."""
        occ = defaultdict(list)    # site -> [(body, local index, image)]
        for bi, B in enumerate(self.bodies):
            for k, (s, n) in enumerate(B.atoms):
                occ[s].append((bi, k, n))
        self.occurrences = occ
        cons = []
        for s, lst in occ.items():
            b0, k0, n0 = lst[0]
            for b1, k1, n1 in lst[1:]:
                cons.append((b0, k0, b1, k1, sub(n0, n1), s))
        self.constraints = cons
        nb = len(self.bodies)
        self.n_linear = sum(1 for B in self.bodies if B.linear)
        # simple Maxwell-type count (6 dof per body, 3 per constraint, cell 6)
        self.maxwell = 6 * nb - self.n_linear - 3 * len(cons) + 6 - 3

    def _find_joints(self):
        """Classify body-body connections for display: hinge axes and ball joints."""
        st = self.st
        pair_shared = defaultdict(list)
        for s, lst in self.occurrences.items():
            for a in range(len(lst)):
                for b in range(a + 1, len(lst)):
                    (b0, _, n0), (b1, _, n1) = lst[a], lst[b]
                    if b0 == b1:
                        continue
                    key = (b0, b1, sub(n0, n1)) if b0 < b1 else (b1, b0, sub(n1, n0))
                    node = (s, n0) if b0 < b1 else (s, n1)
                    pair_shared[key].append(node)
        self.joints = []
        for (b0, b1, t), nodes in pair_shared.items():
            nodes = list(dict.fromkeys(nodes))
            if len(nodes) >= 3 and self._noncollinear(nodes):
                kind = "rigid"
            elif len(nodes) >= 2:
                kind = "hinge"
            else:
                kind = "ball"
            # shift so that the first shared atom is its home-cell copy
            home = np.array(nodes[0][1]) @ st.L
            pts = [self._xyz(s, n) - home for s, n in nodes]
            self.joints.append({"bodies": [b0, b1], "kind": kind,
                                "sites": [s for s, _ in nodes],
                                "points": [p.tolist() for p in pts]})

    # ------------------------------------------------------------ export
    def _sbu_fixes_angles(self, s):
        """True if the rigid-SBU rule fixes every bond angle at site s."""
        st = self.st
        if not (self.opts.rigid_sbu and st.sbu_atom[s]) or st.sbu_id[s] in self.released_sbu:
            return False
        return all(st.sbu_atom[j] and st.sbu_id[j] == st.sbu_id[s] and (st.metal[s] or st.metal[j])
                   for j, _, _ in st.core_nbrs(s))

    def dof_status(self):
        """Effective state of every atom label and bond class, for the DOF tables.

        Atoms: 'flexible' | 'rigid' | 'sbu' (angles fixed by the rigid SBU).
        Bonds: 'hinge' | 'locked' | 'sbu' (torsion fixed because both ends are SBU-fixed).
        """
        st = self.st
        atoms = {}
        for s in st.core_idx:
            lab = st.labels[s]
            if self.flex[s]:
                state = "flexible"
            elif self._sbu_fixes_angles(s):
                state = "sbu"
            else:
                state = "rigid"
            atoms.setdefault(lab, state)
        bonds = {}
        for b, state in self.bond_state.items():
            i, j, _ = st.bonds[b]
            key = ModelOptions.bond_key(st.labels[i], st.labels[j])
            if self._sbu_fixes_angles(i) and self._sbu_fixes_angles(j):
                state = "sbu"
            bonds.setdefault(key, state)
        return atoms, bonds

    def site_body(self):
        """Body index carrying each site (first occurrence), -1 for non-core atoms."""
        out = -np.ones(self.st.N, int)
        for s, lst in self.occurrences.items():
            out[s] = lst[0][0]
        for s, (k, _) in self.st.anchor.items():
            out[s] = out[k]
        return out

    def summary(self):
        kinds = defaultdict(int)
        for j in self.joints:
            kinds[j["kind"]] += 1
        return {
            "n_units": len(self.units),
            "n_bodies": len(self.bodies),
            "n_linear_bodies": self.n_linear,
            "n_constraints": len(self.constraints),
            "n_hinges": kinds["hinge"],
            "n_ball_joints": kinds["ball"],
            "n_rigid_joints": kinds["rigid"],
            "maxwell_count": self.maxwell,
            "n_flex_angle_atoms": int(self.flex.sum()),
            "n_sbu_released": self.n_sbu_released,
            "n_locked_bonds": sum(1 for v in self.bond_state.values() if v == "locked"),
            "n_hinge_bonds": sum(1 for v in self.bond_state.values() if v == "hinge"),
        }
