"""Minimal, dependency-free CIF reader.

Reads cell parameters, symmetry operations and the atom_site loop of the first
data block that contains atoms, and expands the asymmetric unit to the full
unit cell (P1).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np


@dataclass
class CifAtom:
    label: str          # label of the asymmetric-unit atom it was generated from
    element: str
    frac: np.ndarray    # fractional coordinates wrapped to [0, 1)
    occupancy: float = 1.0


@dataclass
class CifData:
    name: str
    cell: tuple          # a, b, c, alpha, beta, gamma (Angstrom, degrees)
    atoms: list
    symops: list
    warnings: list = field(default_factory=list)


# ---------------------------------------------------------------- tokenizer

def _tokenize(text: str):
    """Yield CIF tokens; text fields (;...;) and quoted strings become one token."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("\n".join(buf)).strip()
            continue
        j, n = 0, len(line)
        while j < n:
            ch = line[j]
            if ch.isspace():
                j += 1
            elif ch == "#":
                break
            elif ch in "'\"":
                k = j + 1
                # a closing quote must be followed by whitespace or EOL
                while k < n and not (line[k] == ch and (k + 1 == n or line[k + 1].isspace())):
                    k += 1
                yield line[j + 1:k]
                j = k + 1
            else:
                k = j
                while k < n and not line[k].isspace():
                    k += 1
                yield line[j:k]
                j = k
        i += 1


def _parse_blocks(text: str):
    """Return list of (block_name, items, loops). items: tag->value; loops: list of (tags, rows)."""
    blocks = []
    tokens = list(_tokenize(text))
    i = 0
    cur = None
    while i < len(tokens):
        t = tokens[i]
        tl = t.lower()
        if tl.startswith("data_"):
            cur = (t[5:], {}, [])
            blocks.append(cur)
            i += 1
        elif cur is None:
            i += 1
        elif tl == "loop_":
            i += 1
            tags = []
            while i < len(tokens) and tokens[i].startswith("_"):
                tags.append(tokens[i].lower())
                i += 1
            values = []
            while i < len(tokens):
                v = tokens[i]
                vl = v.lower()
                if v.startswith("_") or vl == "loop_" or vl.startswith("data_"):
                    break
                values.append(v)
                i += 1
            if tags:
                nrow = len(values) // len(tags)
                rows = [values[r * len(tags):(r + 1) * len(tags)] for r in range(nrow)]
                cur[2].append((tags, rows))
        elif t.startswith("_"):
            if i + 1 < len(tokens):
                cur[1][tl] = tokens[i + 1]
            i += 2
        else:
            i += 1
    return blocks


def _num(s: str) -> float:
    """Parse a CIF number, dropping the standard uncertainty '(..)'."""
    s = s.strip()
    if s in (".", "?", ""):
        raise ValueError("missing number")
    return float(re.sub(r"\(.*?\)", "", s))


# ------------------------------------------------------------ symmetry ops

_UNSAFE_TEXT = re.compile(r"[<>&\"`\x00-\x1f]")


def safe_text(text: str, limit: int = 80) -> str:
    """Strip characters that could break out of HTML/attributes; CIF text is untrusted."""
    return _UNSAFE_TEXT.sub("_", str(text))[:limit]


_SYMOP_TOKEN = re.compile(r"\s*(?:(\d+\.\d*|\.\d+|\d+)|([xyz])|([-+*/()]))", re.I)


def _parse_affine(expr: str):
    """Parse one component of a symmetry operation ('-x+1/2', 'x-y', '2x', '1/2+z' ...)
    into (coefficients of x, y, z; constant) with a small recursive-descent parser.

    No eval, no powers, no products of variables, bounded length and nesting: a crafted
    CIF can neither execute code nor make the parser hang."""
    if len(expr) > 40:
        raise ValueError("expression too long")
    toks, pos = [], 0
    while pos < len(expr):
        m = _SYMOP_TOKEN.match(expr, pos)
        if not m or m.end() == pos:
            raise ValueError("unexpected character")
        num, var, op = m.groups()
        if num is not None:
            toks.append(("n", float(num)))
        elif var is not None:
            toks.append(("v", "xyz".index(var.lower())))
        else:
            toks.append(("o", op))
        pos = m.end()
        if expr[pos:].strip() == "":
            break
    k = [0]

    def peek():
        return toks[k[0]] if k[0] < len(toks) else (None, None)

    def take():
        tok = peek()
        k[0] += 1
        return tok

    def factor(depth):
        kind, val = take()
        if kind == "n":
            out = (np.zeros(3), val)
            if peek()[0] == "v":                      # implicit product such as '2x'
                _, idx = take()
                vec = np.zeros(3)
                vec[idx] = val
                out = (vec, 0.0)
            return out
        if kind == "v":
            vec = np.zeros(3)
            vec[val] = 1.0
            return vec, 0.0
        if kind == "o" and val in "+-":
            vec, c = factor(depth)
            return (-vec, -c) if val == "-" else (vec, c)
        if kind == "o" and val == "(" and depth < 3:
            out = expr_(depth + 1)
            if take() != ("o", ")"):
                raise ValueError("unbalanced parentheses")
            return out
        raise ValueError("syntax error")

    def term(depth):
        vec, c = factor(depth)
        while peek() in (("o", "*"), ("o", "/")):
            _, op = take()
            v2, c2 = factor(depth)
            if op == "*":
                if vec.any() and v2.any():
                    raise ValueError("product of variables")
                # (a.X + c)(b.X + c2) with a or b zero
                vec, c = vec * c2 + v2 * c, c * c2
            else:
                if v2.any() or c2 == 0:
                    raise ValueError("division by a variable or zero")
                vec, c = vec / c2, c / c2
        return vec, c

    def expr_(depth):
        vec, c = term(depth)
        while peek() in (("o", "+"), ("o", "-")):
            _, op = take()
            v2, c2 = term(depth)
            vec, c = (vec + v2, c + c2) if op == "+" else (vec - v2, c - c2)
        return vec, c

    vec, c = expr_(0)
    if k[0] != len(toks):
        raise ValueError("trailing characters")
    if np.abs(vec).max() > 3 or abs(c) > 3:
        raise ValueError("coefficient out of range")
    return vec, c


def parse_symop(op: str):
    """'-x+1/2, y, z+1/2' -> (R 3x3, t 3)."""
    parts = [p.strip() for p in op.split(",")]
    if len(parts) != 3:
        raise ValueError(f"bad symmetry operation '{op[:60]}'")
    R = np.zeros((3, 3))
    t = np.zeros(3)
    for row, expr in enumerate(parts):
        try:
            R[row], t[row] = _parse_affine(expr)
        except ValueError as exc:
            raise ValueError(f"bad symmetry operation '{op[:60]}' ({exc})") from None
    return R, t


_CENTERING = {
    "P": [(0, 0, 0)],
    "I": [(0, 0, 0), (.5, .5, .5)],
    "F": [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)],
    "C": [(0, 0, 0), (.5, .5, 0)],
    "A": [(0, 0, 0), (0, .5, .5)],
    "B": [(0, 0, 0), (.5, 0, .5)],
}


def _rot(axis, n):
    """n-fold rotation about an axis of an orthogonal lattice (integer matrix)."""
    a = np.array(axis, float)
    a /= np.linalg.norm(a)
    th = 2 * np.pi / n
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.round(np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K).astype(int)


def symops_from_hm(symbol: str):
    """Symmetry operations from a Hermann-Mauguin symbol for *symmorphic* space groups of
    orthogonal lattices (cubic, tetragonal, orthorhombic), e.g. 'Fm-3m', 'P 4/m m m'.

    Groups with screw axes or glide planes need origin information that the symbol alone
    does not give unambiguously; for those a ValueError is raised."""
    s = symbol.replace(" ", "").replace("_", "")
    if not s or s[0] not in _CENTERING:
        raise ValueError(f"unsupported space-group symbol '{safe_text(symbol, 20)}'")
    lat, rest = s[0], s[1:]
    tokens = re.findall(r"-?[1-6](?:/[a-z])?|[a-z]", rest)
    if "".join(tokens) != rest or any(re.search(r"[abcdegn]", t) for t in tokens) \
            or re.search(r"[1-6][1-6]", rest):
        raise ValueError(f"space group '{safe_text(symbol, 20)}' has screw axes or glide planes")
    if len(tokens) == 3 and tokens[1].lstrip("-").startswith("3"):
        dirs = [(1, 0, 0), (1, 1, 1), (1, 1, 0)]                 # cubic
    elif tokens and tokens[0].lstrip("-")[0] == "4":
        dirs = [(0, 0, 1), (1, 0, 0), (1, 1, 0)]                 # tetragonal
    elif len(tokens) == 3:
        dirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]                 # orthorhombic
    else:
        raise ValueError(f"unsupported space-group symbol '{safe_text(symbol, 20)}'")
    # 3- and 6-fold axes exist only along the cubic body diagonal of an orthogonal lattice;
    # trigonal/hexagonal symbols (P6/mmm, P-3m1 ...) would give wrong integer matrices
    for k, tok in enumerate(tokens):
        if re.search(r"[36]", tok) and not (dirs[1] == (1, 1, 1) and k == 1):
            raise ValueError(f"space group '{safe_text(symbol, 20)}' is trigonal or hexagonal")
    gens = []
    for tok, ax in zip(tokens, dirs):
        for part in tok.split("/"):
            if part == "m":
                gens.append(-_rot(ax, 2))                          # mirror normal to ax
            elif part.startswith("-"):
                gens.append(-_rot(ax, int(part[1:])))             # rotoinversion
            elif part != "1":
                gens.append(_rot(ax, int(part)))
    group ={tuple(np.eye(3, dtype=int).ravel())}
    frontier = [np.eye(3, dtype=int)]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                p = g @ h
                key = tuple(p.ravel())
                if key not in group:
                    group.add(key)
                    new.append(p)
        frontier = new
        if len(group) > 48:                                       # largest point group: m-3m
            raise ValueError(f"unsupported space-group symbol '{safe_text(symbol, 20)}'")
    ops = []
    for key in sorted(group):
        R = np.array(key, float).reshape(3, 3)
        for t in _CENTERING[lat]:
            ops.append((R, np.array(t, float)))
    return ops


# ---------------------------------------------------------------- elements

_ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)")


def _element_from(type_symbol: str | None, label: str) -> str:
    for src in (type_symbol, label):
        if not src or src in (".", "?"):
            continue
        s = src.strip()
        s = s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()
        m = _ELEMENT_RE.match(s)
        if m:
            sym = m.group(1)
            from .elements import COVALENT_RADII
            if sym in COVALENT_RADII:
                return sym
            if sym[0] in COVALENT_RADII:
                return sym[0]
    raise ValueError(f"cannot determine element for atom '{safe_text(label, 32)}'")


# ------------------------------------------------------------------ reader

def cell_matrix(a, b, c, alpha, beta, gamma) -> np.ndarray:
    """Lattice vectors as rows (a along x, b in the xy plane)."""
    al, be, ga = np.radians([alpha, beta, gamma])
    ax = np.array([a, 0.0, 0.0])
    bx = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 1e-12))
    return np.array([ax, bx, [cx, cy, cz]])


def read_cif(text: str, merge_tol: float = 0.25, disorder_merge: float = 0.7,
             max_atoms: int | None = None) -> CifData:
    blocks = _parse_blocks(text)
    block = None
    for b in blocks:
        if any("_atom_site_fract_x" in tags for tags, _ in b[2]):
            block = b
            break
    if block is None:
        raise ValueError("no data block with _atom_site_fract_x found")
    name, items, loops = block
    name = safe_text(name)
    warnings = []

    try:
        cell = tuple(_num(items[k]) for k in (
            "_cell_length_a", "_cell_length_b", "_cell_length_c",
            "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))
    except KeyError as exc:
        raise ValueError(f"missing cell parameter {exc}") from None

    # symmetry operations
    symops = None
    for tags, rows in loops:
        for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz",
                    "_space_group_symop.operation_xyz"):
            if key in tags:
                col = tags.index(key)
                symops = [r[col] for r in rows]
                break
        if symops:
            break
    if not symops:
        for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
            if key in items:
                symops = [items[key]]
    # placeholders such as '?' or '.' carry no information
    symops = [s for s in (symops or []) if s.strip() not in ("?", ".", "")]
    ops = [parse_symop(s) for s in symops]
    if not ops:
        hm = next((items[k] for k in ("_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt")
                   if items.get(k, "?") not in ("?", ".")), None)
        if hm is None or hm.replace(" ", "") == "P1":
            ops = [parse_symop("x,y,z")]
            warnings.append("No symmetry operations found - assuming P1.")
        else:
            try:
                if any(abs(ang - 90.0) > 1e-3 for ang in cell[3:]):
                    raise ValueError("the cell is not orthogonal")
                ops = symops_from_hm(hm)
            except ValueError as exc:
                raise ValueError(f"The CIF lists no symmetry operations and they cannot be derived "
                                 f"safely from the space-group symbol ({exc}). Please re-export the "
                                 f"CIF with symmetry operations (e.g. from Mercury or VESTA).") from None
            warnings.append(f"No symmetry operations in the CIF - generated {len(ops)} operations "
                            f"from the space-group symbol {safe_text(hm, 20)}.")

    # atom sites
    tags, rows = next((tg, r) for tg, r in loops if "_atom_site_fract_x" in tg)

    def col(name):
        return tags.index(name) if name in tags else None

    c_label, c_type = col("_atom_site_label"), col("_atom_site_type_symbol")
    c_x, c_y, c_z = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    c_occ, c_dis = col("_atom_site_occupancy"), col("_atom_site_disorder_group")
    c_calc = col("_atom_site_calc_flag")
    c_asm = col("_atom_site_disorder_assembly")

    # disorder: within each disorder assembly keep only the lowest-numbered group (SHELX PART);
    # assemblies often use their own numbers (e.g. A: 1/2, B: 3/4)
    def _group(r):
        if c_dis is None or r[c_dis] in (".", "?"):
            return None
        try:
            return abs(int(r[c_dis]))          # negative parts: special-position disorder
        except ValueError:
            return None

    def _assembly(r):
        return r[c_asm] if c_asm is not None and r[c_asm] not in (".", "?") else ""
    first_group = {}
    for r in rows:
        g = _group(r)
        if g:
            first_group[_assembly(r)] = min(g, first_group.get(_assembly(r), g))

    asym = []
    skipped_disorder = 0
    for r in rows:
        label = r[c_label] if c_label is not None else (r[c_type] if c_type is not None else "X")
        label = safe_text(label, 32)
        if c_calc is not None and r[c_calc].lower() == "dum":
            continue
        try:
            frac = np.array([_num(r[c_x]), _num(r[c_y]), _num(r[c_z])])
        except ValueError:
            continue
        occ = 1.0
        if c_occ is not None:
            try:
                occ = _num(r[c_occ])
            except ValueError:
                occ = 1.0
        g = _group(r)
        if g and g != first_group.get(_assembly(r), g):
            skipped_disorder += 1
            continue
        el = _element_from(r[c_type] if c_type is not None else None, label)
        asym.append((label, el, frac, occ))
    if skipped_disorder:
        warnings.append(f"Skipped {skipped_disorder} atoms of minor disorder components "
                        "(the first component of each disorder assembly is kept).")
    if any(o < 0.99 for *_, o in asym):
        warnings.append("Structure contains partially occupied sites; they are kept as full atoms.")

    # expand by symmetry, merge duplicates. Partially occupied symmetry copies of the same
    # site that lie closer than a bond (disorder across a mirror or axis) are one disordered
    # atom: they are merged at their average position. A grid of buckets keeps the search
    # local (linear in the number of atoms).
    L = cell_matrix(*cell)
    reach = max(merge_tol, disorder_merge)
    widths = abs(np.linalg.det(L)) / np.array([np.linalg.norm(np.cross(L[(k + 1) % 3], L[(k + 2) % 3]))
                                               for k in range(3)])
    ngrid = np.maximum(1, np.floor(widths / reach)).astype(int)
    grid: dict = {}
    atoms: list[CifAtom] = []
    n_merged = n_mixed = 0
    for label, el, frac, occ in asym:
        for R, t in ops:
            f = (R @ frac + t) % 1.0
            f[f > 1 - 1e-6] = 0.0
            cellidx = np.floor(f * ngrid).astype(int) % ngrid
            dup = False
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        key = tuple((cellidx + (dx, dy, dz)) % ngrid)
                        for a in grid.get(key, ()):
                            df = f - a.frac
                            df -= np.round(df)
                            d = np.linalg.norm(df @ L)
                            if a.element != el:
                                # substitutional disorder (e.g. Si/Al on one zeolite T site):
                                # the site is represented by its majority element
                                if d < merge_tol and min(occ, a.occupancy) < 0.99:
                                    if occ > a.occupancy:
                                        a.element, a.label, a.occupancy = el, label, occ
                                    n_mixed += 1
                                    dup = True
                                    break
                                continue
                            if d < merge_tol:
                                dup = True
                            elif occ < 0.99 and a.label == label and d < disorder_merge:
                                a.frac = (a.frac + df / 2) % 1.0
                                n_merged += 1
                                dup = True
                            if dup:
                                break
                        if dup:
                            break
                    if dup:
                        break
                if dup:
                    break
            if not dup:
                a = CifAtom(label, el, f, occ)
                atoms.append(a)
                grid.setdefault(tuple(cellidx), []).append(a)
                if max_atoms and len(atoms) > max_atoms:
                    raise ValueError(f"The structure has more than {max_atoms} atoms per unit cell, "
                                     "which exceeds the limit of this server.")
    if n_mixed:
        warnings.append(f"{n_mixed} mixed-occupancy sites (e.g. Si/Al) are represented by their majority element.")
    if n_merged:
        warnings.append(f"Merged {n_merged} partially occupied symmetry copies of disordered sites "
                        f"(closer than {disorder_merge} Å) at their average positions.")
    return CifData(name=name, cell=cell, atoms=atoms, symops=symops, warnings=warnings)
