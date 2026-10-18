import numpy as np, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from gapclosure import phenylene, bipartite, BondModel

def build(choices, cc=1.40, ch=1.08):
    g = phenylene(choices)
    P = np.array(g['pos']) * cc
    deg = np.zeros(g['n'], int); nb = [[] for _ in range(g['n'])]
    for u, v in g['edges']: deg[u] += 1; deg[v] += 1; nb[u].append(v); nb[v].append(u)
    H = []
    for i in range(g['n']):
        if deg[i] == 2:
            d = 2*P[i] - P[nb[i][0]] - P[nb[i][1]]; d /= np.linalg.norm(d); H.append(P[i] + ch*d)
    xyz = np.vstack([np.c_[P, np.zeros(len(P))], np.c_[np.array(H), np.zeros(len(H))]])
    sym = ['C']*g['n'] + ['H']*len(H)
    return g, sym, xyz

def atom_str(sym, xyz):
    return '; '.join(f'{s} {x:.8f} {y:.8f} {z:.8f}' for s, (x, y, z) in zip(sym, xyz))
