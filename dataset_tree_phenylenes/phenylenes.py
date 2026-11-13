"""Tree phenylenes: L/A classes and the distinct molecules in each class."""
import itertools
import networkx as nx
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")


def _canon(T, lab, v, parent=None):
    return "(" + lab.get(v, "") + "".join(sorted(_canon(T, lab, c, v) for c in T[v] if c != parent)) + ")"


def canon_tree(T, lab):
    """Canonical string of a tree whose degree-2 vertices carry L or A."""
    return min(_canon(T, lab, c) for c in nx.center(T))


def classes(N):
    """All L/A classes with N hexagons: {canonical string: (tree, labels)}."""
    out = {}
    trees = [nx.empty_graph(1)] if N == 1 else nx.nonisomorphic_trees(N)
    for T in trees:
        if max((d for _, d in T.degree), default=0) > 3:
            continue
        deg2 = [v for v in T if T.degree(v) == 2]
        for labs in itertools.product("LA", repeat=len(deg2)):
            lab = dict(zip(deg2, labs))
            out.setdefault(canon_tree(T, lab), (nx.Graph(T), dict(lab)))
    return out


def realisations(T, lab):
    """Every way to place the class in the plane: lists of links (i, port, j),
    hexagons numbered in BFS order, the new hexagon's bond 0 shared with its square."""
    root = nx.center(T)[0]
    order, parent = [root], {root: None}
    for v in order:
        for c in sorted(T[v]):
            if c not in parent:
                parent[c] = v
                order.append(c)
    idx = {v: k for k, v in enumerate(order)}
    options = []                     # per hexagon: list of [(child, port), ...]
    for v in order:
        kids = [c for c in sorted(T[v]) if c != parent[v]]
        d = T.degree(v)
        if parent[v] is None:
            if d == 0:
                opts = [[]]
            elif d == 1:
                opts = [[(kids[0], 0)]]
            elif d == 2:
                opts = [list(zip(kids, (0, 3)))] if lab[v] == "L" else \
                       [list(zip(kids, (0, 2))), list(zip(kids, (2, 0)))]
            else:
                opts = [list(zip(kids, p)) for p in itertools.permutations((0, 2, 4))]
        else:
            if d == 1:
                opts = [[]]
            elif d == 2:
                opts = [[(kids[0], 3)]] if lab[v] == "L" else [[(kids[0], 2)], [(kids[0], 4)]]
            else:
                opts = [list(zip(kids, (2, 4))), list(zip(kids, (4, 2)))]
        options.append(opts)
    for choice in itertools.product(*options):
        yield [(idx[v], p, idx[c]) for v, opt in zip(order, choice) for c, p in opt]


def graph(links, n_hex):
    """Bonds of the molecule; hexagon k has atoms 6k..6k+5, bond p = (p, p+1)."""
    bonds = [(6 * k + p, 6 * k + (p + 1) % 6) for k in range(n_hex) for p in range(6)]
    for i, a, j in links:
        bonds.append((6 * i + a, 6 * j))
        bonds.append((6 * i + (a + 1) % 6, 6 * j + 1))
    return bonds


def graph_key(bonds, n):
    """Canonical SMILES of the carbon graph with all bonds single: equal iff isomorphic."""
    m = Chem.RWMol()
    for _ in range(n):
        m.AddAtom(Chem.Atom(6))
    for u, v in bonds:
        m.AddBond(u, v, Chem.BondType.SINGLE)
    return Chem.MolToSmiles(m.GetMol())


def molecules(T, lab):
    """Distinct molecules of a class: {graph key: links}."""
    n = T.number_of_nodes()
    out = {}
    for links in realisations(T, lab):
        key = graph_key(graph(links, n), 6 * n)
        out.setdefault(key, links)
    return out
