"""Carbon skeletons: the graph of the pi-centres, with 2D coordinates for drawing.

A skeleton can be read from SMILES or a molfile (needs RDKit), from an XYZ file
(carbon positions; bonds by distance), from an edge list, or built with the
functions at the end of this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

S3 = math.sqrt(3.0)


@dataclass
class Skeleton:
    """Graph of the pi-centres. Atoms are 0..n-1, bonds are pairs (u, v) with u < v."""

    n: int
    bonds: list
    pos: np.ndarray | None = None
    name: str = ""
    source_atoms: list | None = field(default=None, repr=False)  # index in the input molecule

    def __post_init__(self):
        seen = set()
        clean = []
        for u, v in self.bonds:
            u, v = int(u), int(v)
            if u == v or not (0 <= u < self.n and 0 <= v < self.n):
                raise ValueError(f"invalid bond ({u}, {v})")
            e = (min(u, v), max(u, v))
            if e not in seen:
                seen.add(e)
                clean.append(e)
        self.bonds = clean
        if self.pos is not None:
            self.pos = np.asarray(self.pos, float).reshape(self.n, 2)

    # ------------------------------------------------------------------ basic graph data
    @property
    def m(self):
        return len(self.bonds)

    def neighbours(self):
        adj = [[] for _ in range(self.n)]
        for u, v in self.bonds:
            adj[u].append(v)
            adj[v].append(u)
        return adj

    def components(self):
        adj = self.neighbours()
        comp = [-1] * self.n
        out = []
        for s in range(self.n):
            if comp[s] >= 0:
                continue
            comp[s] = len(out)
            stack, members = [s], [s]
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if comp[v] < 0:
                        comp[v] = len(out)
                        stack.append(v)
                        members.append(v)
            out.append(sorted(members))
        return out

    def subskeleton(self, atoms):
        atoms = sorted(atoms)
        idx = {a: k for k, a in enumerate(atoms)}
        bonds = [(idx[u], idx[v]) for u, v in self.bonds if u in idx and v in idx]
        pos = None if self.pos is None else self.pos[atoms]
        src = atoms if self.source_atoms is None else [self.source_atoms[a] for a in atoms]
        return Skeleton(len(atoms), bonds, pos, self.name, src)

    def bipartition(self):
        """Colour classes (starred, unstarred) of an alternant skeleton."""
        adj = self.neighbours()
        col = [-1] * self.n
        for s in range(self.n):
            if col[s] >= 0:
                continue
            col[s] = 0
            stack = [s]
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if col[v] < 0:
                        col[v] = 1 - col[u]
                        stack.append(v)
                    elif col[v] == col[u]:
                        raise ValueError(
                            "the skeleton has an odd ring (non-alternant); the closure "
                            "threshold theory needs an alternant (bipartite) skeleton")
        rows = [a for a in range(self.n) if col[a] == 0]
        cols = [a for a in range(self.n) if col[a] == 1]
        return rows, cols

    def biadjacency(self):
        """T with H = -[[0, T], [T^T, 0]], and for each bond its (row, column) entry."""
        rows, cols = self.bipartition()
        if len(rows) != len(cols):
            raise ValueError(
                f"unequal colour classes ({len(rows)} and {len(cols)} atoms): "
                "the skeleton has no Kekule structure and a zero eigenvalue at any coupling")
        ri = {a: i for i, a in enumerate(rows)}
        ci = {a: j for j, a in enumerate(cols)}
        T = np.zeros((len(rows), len(cols)))
        entry = []
        for u, v in self.bonds:
            i, j = (ri[u], ci[v]) if u in ri else (ri[v], ci[u])
            T[i, j] = 1.0
            entry.append((i, j))
        return T, rows, cols, entry

    def layout(self):
        """Coordinates for drawing; computed once if the input had none."""
        if self.pos is None:
            self.pos = _layout(self)
        return self.pos

    def __repr__(self):
        nm = f" {self.name!r}" if self.name else ""
        return f"<Skeleton{nm}: {self.n} atoms, {self.m} bonds>"


# ====================================================================== readers

def _need_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "reading SMILES or molfiles needs RDKit: pip install 'gapclosure[rdkit]'") from exc


def _pi_centres_rdkit(mol):
    """Carbon atoms that take part in a pi bond: aromatic, or with a double bond to carbon."""
    from rdkit import Chem
    centres = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C":
            if a.GetIsAromatic() or any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()):
                raise ValueError(
                    f"atom {a.GetIdx()} ({a.GetSymbol()}) is part of the pi system; only "
                    "hydrocarbons are covered (all Coulomb integrals equal)")
            continue
        if any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
            raise ValueError("triple bonds are not covered")
        if a.GetIsAromatic() or any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()):
            centres.append(a.GetIdx())
    return centres


def _from_rdkit(mol, name):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0 or mol.GetConformer().Is3D():
        AllChem.Compute2DCoords(mol)
    centres = _pi_centres_rdkit(mol)
    if not centres:
        raise ValueError("no pi-centres found")
    idx = {a: k for k, a in enumerate(centres)}
    bonds = []
    for b in mol.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if u in idx and v in idx:
            bonds.append((idx[u], idx[v]))
    conf = mol.GetConformer()
    pos = np.array([[conf.GetAtomPosition(a).x, conf.GetAtomPosition(a).y] for a in centres])
    pos = _normalise(pos, bonds)
    return Skeleton(len(centres), bonds, pos, name, centres)


def from_smiles(smiles: str, name: str = "") -> Skeleton:
    """Pi skeleton of a conjugated hydrocarbon given as SMILES. Needs RDKit."""
    _need_rdkit()
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
    return _from_rdkit(mol, name or smiles)


def from_molfile(path: str, name: str = "") -> Skeleton:
    """Pi skeleton from a MDL molfile / SDF (first record). Needs RDKit."""
    _need_rdkit()
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(path), removeHs=True)
    if mol is None:
        raise ValueError(f"RDKit could not read {path}")
    return _from_rdkit(mol, name or str(path))


def from_xyz(path_or_text: str, name: str = "", cutoff: float = 1.6) -> Skeleton:
    """Pi skeleton from an XYZ geometry: carbon atoms with three or fewer carbon
    neighbours closer than `cutoff` (angstrom), and carbons bonded by distance.
    Saturated carbons (four neighbours including H) are dropped. Coordinates are
    projected on the best plane for drawing."""
    text = path_or_text
    if "\n" not in path_or_text:
        with open(path_or_text) as f:
            text = f.read()
        name = name or str(path_or_text)
    lines = [ln.split() for ln in text.strip().splitlines()]
    try:
        nat = int(lines[0][0])
        body = lines[2:2 + nat]
    except (ValueError, IndexError):
        body = [ln for ln in lines if len(ln) >= 4]
    sym = [ln[0].capitalize() for ln in body]
    X = np.array([[float(t) for t in ln[1:4]] for ln in body])
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    heavy_or_h = D < np.where(np.array(sym)[:, None] == "H", 1.25, cutoff)
    np.fill_diagonal(heavy_or_h, False)
    centres = []
    for a, s in enumerate(sym):
        if s != "C":
            continue
        nb = np.nonzero(heavy_or_h[a] | heavy_or_h[:, a])[0]
        if len(nb) <= 3:
            centres.append(a)
    idx = {a: k for k, a in enumerate(centres)}
    bonds = [(idx[a], idx[b]) for a in centres for b in centres if a < b and D[a, b] < cutoff]
    P = X[centres] - X[centres].mean(0)
    _, _, Vt = np.linalg.svd(P, full_matrices=False)
    pos = _normalise(P @ Vt[:2].T, bonds)
    return Skeleton(len(centres), bonds, pos, name or "xyz", centres)


def from_edges(bonds, n: int | None = None, pos=None, name: str = "") -> Skeleton:
    """Skeleton from a list of bonds (pairs of atom indices starting at 0)."""
    bonds = [tuple(e) for e in bonds]
    if n is None:
        n = 1 + max(max(e) for e in bonds)
    return Skeleton(n, bonds, pos, name)


def from_points(points, name: str = "", tol: float = 1e-3) -> Skeleton:
    """Skeleton whose atoms are the given points and whose bonds join points at distance one."""
    P = np.asarray(points, float)
    D = np.linalg.norm(P[:, None] - P[None], axis=-1)
    bonds = [(u, v) for u in range(len(P)) for v in range(u + 1, len(P)) if abs(D[u, v] - 1) < tol]
    return Skeleton(len(P), bonds, P, name)


def _normalise(pos, bonds):
    """Scale coordinates so that the median bond has length one."""
    if bonds:
        L = np.median([np.linalg.norm(pos[u] - pos[v]) for u, v in bonds])
        if L > 0:
            pos = pos / L
    return pos - pos.mean(0)


def _layout(sk):
    """Coordinates for a skeleton read without any: RDKit depiction if available,
    else networkx Kamada-Kawai, else a spectral layout."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from .kekule import one_kekule
        rw = Chem.RWMol()
        for _ in range(sk.n):
            rw.AddAtom(Chem.Atom(6))
        M = one_kekule(sk)
        for u, v in sk.bonds:
            bt = Chem.BondType.DOUBLE if M is not None and (u, v) in M else Chem.BondType.SINGLE
            rw.AddBond(u, v, bt)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        AllChem.Compute2DCoords(mol)
        c = mol.GetConformer()
        P = np.array([[c.GetAtomPosition(a).x, c.GetAtomPosition(a).y] for a in range(sk.n)])
        return _normalise(P, sk.bonds)
    except Exception:
        pass
    try:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(range(sk.n))
        G.add_edges_from(sk.bonds)
        p = nx.kamada_kawai_layout(G)
        return _normalise(np.array([p[a] for a in range(sk.n)]), sk.bonds)
    except Exception:
        pass
    A = np.zeros((sk.n, sk.n))
    for u, v in sk.bonds:
        A[u, v] = A[v, u] = 1
    L = np.diag(A.sum(1)) - A
    w, V = np.linalg.eigh(L)
    return _normalise(V[:, 1:3], sk.bonds)


# ====================================================================== builders

def _hexagon_points(cx, cy):
    """Pointy-top hexagon of unit side: its left and right bonds are vertical."""
    return [(cx + math.cos(math.radians(30 + 60 * k)), cy + math.sin(math.radians(30 + 60 * k)))
            for k in range(6)]


def _dedupe(points):
    out = []
    for p in points:
        if not any(abs(p[0] - q[0]) < 1e-6 and abs(p[1] - q[1]) < 1e-6 for q in out):
            out.append(p)
    return out


def acene(h: int) -> Skeleton:
    """Linear acene with h >= 1 rings (h = 1 benzene, 2 naphthalene, ...)."""
    if h < 1:
        raise ValueError("an acene has at least one ring")
    return acene_chain([h], name=["benzene", "naphthalene", "anthracene", "tetracene",
                                  "pentacene"][h - 1] if h <= 5 else f"[{h}]acene")


def acene_chain(units, name: str = "") -> Skeleton:
    """Acenes joined end to end through four-membered rings.

    units[k] is the number of rings of the k-th acene; 0 stands for a single bond
    (an ethylene unit). acene_chain([1, 1]) is biphenylene, acene_chain([0, 1])
    benzocyclobutadiene, acene_chain([hA, hB]) the [hA, hB] skeletons of the
    acene-cyclobutadiene-acene class, acene_chain([1]*N) linear [N]phenylene."""
    pts = []
    x0 = 0.0
    for h in units:
        if h < 0:
            raise ValueError("ring counts must be >= 0")
        if h == 0:
            pts += [(x0, 0.5), (x0, -0.5)]
            width = 0.0
        else:
            for k in range(h):
                pts += _hexagon_points(x0 + S3 / 2 + k * S3, 0.0)
            width = h * S3
        x0 += width + 1.0
    nm = name or ("biphenylene" if list(units) == [1, 1] else
                  "benzocyclobutadiene" if list(units) in ([0, 1], [1, 0]) else
                  "acene chain " + "-".join(map(str, units)))
    sk = from_points(_dedupe(pts), nm)
    return Skeleton(sk.n, sk.bonds, sk.pos - sk.pos.mean(0), nm)


def annulene(L: int) -> Skeleton:
    """Ring of L atoms."""
    R = 1 / (2 * math.sin(math.pi / L))
    pos = [(R * math.cos(2 * math.pi * k / L), R * math.sin(2 * math.pi * k / L)) for k in range(L)]
    return Skeleton(L, [(k, (k + 1) % L) for k in range(L)], pos, f"[{L}]annulene")


def benzenoid(cells, name: str = "") -> Skeleton:
    """Benzenoid from hexagon cells in axial coordinates (i, j)."""
    pts = []
    for i, j in cells:
        pts += _hexagon_points(S3 * (i + j / 2), 1.5 * j)
    return from_points(_dedupe(pts), name or f"benzenoid {list(cells)}")


# ---------------------------------------------------------------- phenylenes (trees of hexagons)

class _Builder:
    """Planar construction by fusing regular polygons onto existing bonds."""

    def __init__(self):
        ang = math.pi / 2 - np.arange(6) * math.pi / 3
        self.pos = [(math.cos(a), math.sin(a)) for a in ang]
        self.bonds = [(i, (i + 1) % 6) for i in range(6)]
        self.rings = [list(range(6))]

    def fuse(self, u, v, k):
        """Fuse a regular k-gon onto bond (u, v), away from the ring that contains it.
        Returns the new ring as a path [u, new atoms..., v]."""
        P = np.array(self.pos, float)
        host = next(r for r in self.rings if u in r and v in r)
        c_host = P[host].mean(0)
        pu, pv = P[u], P[v]
        mid, d = (pu + pv) / 2, pv - pu
        nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
        if np.dot(nrm, mid - c_host) < 0:
            nrm = -nrm
        L = np.linalg.norm(d)
        c = mid + L / (2 * math.tan(math.pi / k)) * nrm
        R = L / (2 * math.sin(math.pi / k))
        a_u = math.atan2(*(pu - c)[::-1])
        a_v = math.atan2(*(pv - c)[::-1])
        step = 2 * math.pi / k
        wrap = lambda x: (x + math.pi) % (2 * math.pi) - math.pi  # noqa: E731
        sgn = 1 if abs(wrap(a_u + step * (k - 1) - a_v)) < 1e-6 else -1
        start = len(self.pos)
        new = list(range(start, start + k - 2))
        for i in range(k - 2):
            a = a_u + sgn * step * (i + 1)
            self.pos.append((c[0] + R * math.cos(a), c[1] + R * math.sin(a)))
        path = [u] + new + [v]
        self.bonds += list(zip(path[:-1], path[1:]))
        self.rings.append(path)
        return path

    def skeleton(self, name):
        P = np.array(self.pos)
        D = np.linalg.norm(P[:, None] - P[None], axis=-1) + 10 * np.eye(len(P))
        if D.min() < 0.5:
            import warnings
            warnings.warn("rings overlap in the plane; the drawing will overlap, "
                          "the computation is unaffected", stacklevel=3)
        return Skeleton(len(P), self.bonds, P - P.mean(0), name)


def phenylene_tree(links, name: str = "") -> Skeleton:
    """Phenylene whose hexagons form a tree.

    Hexagon 0 is placed first. Each link (i, a, j) attaches a new hexagon j to the
    existing hexagon i through a four-membered ring fused on bond a of hexagon i;
    the bonds of a hexagon are numbered 0..5 around the ring, and bond 0 of a new
    hexagon is the one it shares with its four-membered ring. Links must be given
    so that i already exists, e.g. [(0, 3, 1), (1, 2, 2)]."""
    b = _Builder()
    ring_atoms = {0: list(range(6))}            # hexagon -> atoms c0..c5 (bond p = c_p c_{p+1})
    for link in links:
        i, a, j = link[:3]
        if i not in ring_atoms or j in ring_atoms:
            raise ValueError(f"link {link}: hexagon {i} must exist and {j} must be new")
        c = ring_atoms[i]
        sq = b.fuse(c[a % 6], c[(a + 1) % 6], 4)      # [c_a, s1, s2, c_{a+1}]
        hx = b.fuse(sq[2], sq[1], 6)                   # [s2, p1, p2, p3, p4, s1]
        # bond 0 of the new hexagon is (s1, s2), then around the ring
        ring_atoms[j] = [hx[5], hx[0], hx[1], hx[2], hx[3], hx[4]]
    return b.skeleton(name or f"phenylene tree {list(links)}")


def phenylene(word: str = "", helical: bool = False) -> Skeleton:
    """[N]phenylene chain from its L/A word for the N-2 interior hexagons.

    'L' continues linearly (para), 'A' turns (angular). Angular turns alternate in
    side (zigzag) unless helical=True. The threshold depends only on the word;
    the side only changes the drawing. phenylene('') is biphenylene."""
    word = word.upper()
    if set(word) - set("LA"):
        raise ValueError("use L (linear) and A (angular)")
    links, port, side = [], 3, 1
    links.append((0, 3, 1))
    for k, letter in enumerate(word, start=1):
        if letter == "L":
            port = 3
        else:
            # the ring labels reverse at each step, so a fixed port gives a zigzag
            port = 2 if side > 0 else 4
            if helical:
                side = -side
        links.append((k, port, k + 1))
    n = len(word) + 2
    shape = "" if "AA" not in word else (" helical" if helical else " zigzag")
    nm = {"": "biphenylene"}.get(word, f"[{n}]phenylene {word}{shape}")
    return phenylene_tree(links, nm)


def benzocyclobutadiene() -> Skeleton:
    return acene_chain([0, 1])


def biphenylene() -> Skeleton:
    return phenylene("")


def triangular_phenylene() -> Skeleton:
    """Triangular [4]phenylene: a central hexagon with three outer ones."""
    return phenylene_tree([(0, 0, 1), (0, 2, 2), (0, 4, 3)], "triangular [4]phenylene")


NAMED = {
    "benzene": lambda: acene(1),
    "naphthalene": lambda: acene(2),
    "anthracene": lambda: acene(3),
    "tetracene": lambda: acene(4),
    "pentacene": lambda: acene(5),
    "benzocyclobutadiene": benzocyclobutadiene,
    "biphenylene": biphenylene,
    "linear[3]phenylene": lambda: phenylene("L"),
    "angular[3]phenylene": lambda: phenylene("A"),
    "linear[4]phenylene": lambda: phenylene("LL"),
    "zigzag[4]phenylene": lambda: phenylene("AA"),
    "helical[4]phenylene": lambda: phenylene("AA", helical=True),
    "bent[4]phenylene": lambda: phenylene("LA"),
    "triangular[4]phenylene": triangular_phenylene,
}


def named(name: str) -> Skeleton:
    key = name.lower().replace(" ", "").replace("-", "")
    for k, f in NAMED.items():
        if k.replace("-", "") == key:
            sk = f()
            sk.name = k
            return sk
    raise KeyError(f"unknown name {name!r}; known: {', '.join(NAMED)}")
