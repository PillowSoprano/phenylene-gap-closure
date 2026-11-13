"""Closure thresholds of all phenylenes whose hexagons form a tree, N = 1..NMAX.

One row per L/A class (the threshold depends only on the class), plus one row per
molecule. Run:  python build.py 12
"""
import csv
import json
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
import networkx as nx
from rdkit import Chem, RDLogger

import gapclosure as gc
from gapclosure.kekule import one_kekule
from gapclosure.tree import PhenyleneTree
import phenylenes as E

RDLogger.DisableLog("rdApp.*")
BITS = 40


# ---------------------------------------------------------------- Kekule count on the tree
def kekule_count(tree):
    """Number of Kekule structures: the transfer of PhenyleneTree with all signs positive."""
    def pair(branches, port):
        factors = [1] * 6
        even_scale = odd_scale = 1
        for edge, det, deriv in branches:
            factors[edge] = det + deriv          # bond kept, or replaced by the two connecting bonds
            if edge % 2 == 0:
                odd_scale *= det
            else:
                even_scale *= det
        even = even_scale * factors[0] * factors[2] * factors[4]
        odd = odd_scale * factors[1] * factors[3] * factors[5]
        if port is None:
            return even + odd, 0
        return even + odd, (odd if port % 2 else even)

    states = {}
    for i in reversed(tree.order):
        br = [(tree.adj[i][j], *states[j]) for j in tree.children[i]]
        states[i] = pair(br, tree.adj[i].get(tree.parent[i]))
    return states[tree.root][0]


# ---------------------------------------------------------------- per molecule
def smiles(sk, mapped=False):
    """Aromatic SMILES of the carbon skeleton; with atom maps = atom index + 1 if mapped."""
    M = one_kekule(sk)
    m = Chem.RWMol()
    for a in range(sk.n):
        at = Chem.Atom(6)
        if mapped:
            at.SetAtomMapNum(a + 1)
        m.AddAtom(at)
    for u, v in sk.bonds:
        m.AddBond(u, v, Chem.BondType.DOUBLE if (u, v) in M else Chem.BondType.SINGLE)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def clash(sk):
    P = sk.pos
    D = np.linalg.norm(P[:, None] - P[None], axis=-1) + 10 * np.eye(sk.n)
    return D.min() < 0.5


def word_of(T, lab):
    """L/A word of a chain (L/A notation), the smaller of the two reading directions."""
    if max((d for _, d in T.degree), default=0) > 2:
        return None
    if T.number_of_nodes() <= 2:
        return ""
    ends = [v for v in T if T.degree(v) == 1]
    path = nx.shortest_path(T, *ends)
    w = "".join(lab[v] for v in path[1:-1])
    return min(w, w[::-1])


def one_class(args):
    key, T, lab = args
    import warnings
    warnings.simplefilter("ignore")
    N = T.number_of_nodes()
    mols = E.molecules(T, lab)
    entries, rep = [], None
    for gkey, links in mols.items():
        sk = gc.phenylene_tree(links)
        T0, *_ = sk.biadjacency()
        gap = float(2 * np.linalg.svd(T0, compute_uv=False)[-1])
        entries.append(dict(links=links, smiles=smiles(sk), gap=gap, overlaps=bool(clash(sk))))
        if rep is None or (rep[1].get("overlaps") and not entries[-1]["overlaps"]):
            rep = (sk, entries[-1])
    sk, rep_entry = rep
    t = time.time()
    if N == 1:
        r = gc.closure_threshold(sk)
        lo = hi = Fraction(1)
        tree = PhenyleneTree(1, [])
    else:
        r = gc.closure_threshold(sk, method="tree", bits=BITS)
        tree, _ = gc.recognize_phenylene_tree(sk)
        lo, hi = r.bracket
    secs = time.time() - t
    K = kekule_count(tree)
    ASC = int(r.ASC)
    deg = dict(T.degree)
    row = dict(
        N=N, la_tree=key, word=word_of(T, lab),
        n_L=sum(1 for v in lab.values() if v == "L"), n_A=sum(1 for v in lab.values() if v == "A"),
        n_branch=sum(1 for d in deg.values() if d == 3),
        K=K, ASC=ASC, minority=(K - ASC) // 2,
        rho=r.rho, rho_lower=str(lo), rho_upper=str(hi), q=r.q,
        n_molecules=len(entries),
        gap_min=min(e["gap"] for e in entries), gap_max=max(e["gap"] for e in entries),
        representative_smiles=rep_entry["smiles"], seconds=secs,
    )
    detail = dict(
        la_tree=key, N=N, rho=r.rho, rho_lower=str(lo), rho_upper=str(hi), method=r.method,
        representative=dict(
            links=[list(x) for x in rep_entry["links"]],
            smiles_mapped=smiles(sk, mapped=True),
            coords=[[round(float(x), 4), round(float(y), 4)] for x, y in sk.pos],
            closing_pattern=[[u, v, int(s)] for (u, v), s in r.pattern.items()],
            closing_orbital=[round(float(w), 6) for w in r.orbital] if r.orbital is not None else None,
        ),
    )
    mol_rows = [dict(la_tree=key, N=N, smiles=e["smiles"], huckel_gap=e["gap"],
                     overlaps_in_plane=e["overlaps"], links=json.dumps(e["links"])) for e in entries]
    return row, detail, mol_rows


ROW_FIELDS = None
MOL_FIELDS = ["id", "la_tree", "N", "smiles", "huckel_gap", "overlaps_in_plane", "links"]


def _fmt(v, digits=None):
    return repr(v) if isinstance(v, float) else ("" if v is None else v)


def main(nmax):
    tag = f"N1-{nmax}"
    files = dict(rows=f"tree_phenylenes_{tag}.csv", mols=f"tree_phenylene_molecules_{tag}.csv",
                 pats=f"tree_phenylene_patterns_{tag}.jsonl")
    for f in files.values():
        open(f, "w").close()
    first = True
    with Pool(2) as pool:
        for N in range(1, nmax + 1):
            t = time.time()
            C = E.classes(N)
            jobs = [(k, T, lab) for k, (T, lab) in sorted(C.items())]
            out = pool.map(one_class, jobs, chunksize=4)
            out.sort(key=lambda o: (-o[0]["rho"], o[0]["la_tree"]))
            with open(files["rows"], "a", newline="") as fr, open(files["mols"], "a", newline="") as fm, \
                    open(files["pats"], "a") as fp:
                wr = csv.DictWriter(fr, fieldnames=["id"] + list(out[0][0]))
                wm = csv.DictWriter(fm, fieldnames=MOL_FIELDS)
                if first:
                    wr.writeheader()
                    wm.writeheader()
                    first = False
                for k, (row, det, mr) in enumerate(out, 1):
                    cid = f"PT{N:02d}-{k:05d}"
                    wr.writerow({k2: _fmt(v) for k2, v in dict(id=cid, **row).items()})
                    fp.write(json.dumps(dict(id=cid, **det)) + "\n")
                    for m in mr:
                        m["id"] = cid
                        wm.writerow({k2: _fmt(m[k2]) for k2 in MOL_FIELDS})
            print(f"N={N:2d}: {len(C):5d} classes, {sum(len(o[2]) for o in out):6d} molecules, "
                  f"{time.time() - t:6.1f} s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
