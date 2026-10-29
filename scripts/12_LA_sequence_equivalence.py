import sys, itertools, json, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, networkx as nx
from networkx.algorithms import isomorphism as iso
from gapclosure import *
def inc(g):
    T,L,R,bix=bipartite(g); terms=kekule_terms(T); G=nx.Graph()
    for k in range(len(bix)): G.add_node(('e',k),kind='edge')
    for m,(s,p) in enumerate(terms):
        G.add_node(('m',m),kind='plus' if s>0 else 'minus')
        for k,(i,j) in enumerate(bix):
            if p[i]==j: G.add_edge(('m',m),('e',k))
    return G
def la(ch): return ''.join('L' if c=='p' else 'A' for c in ch)
for N in (4,5,6):
    t=time.time(); groups={}
    for ch in itertools.product('pLR',repeat=N-2):
        groups.setdefault(la(ch),[]).append(ch)
    ok=True; reps={}
    for seq,chs in sorted(groups.items()):
        Gs=[inc(phenylene(list(c))) for c in chs]
        for c,G in zip(chs[1:],Gs[1:]):
            same=iso.GraphMatcher(Gs[0],G,node_match=lambda a,b:a['kind']==b['kind']).is_isomorphic()
            ok&=same
            if not same: print('  FAIL',N,seq,chs[0],c)
        reps[seq]=chs[0]
    # different sequences (up to reversal) not equivalent?
    print(f'[{N}]: {len(groups)} L/A sequences, {3**(N-2)} graphs; all same-sequence pairs equivalent: {ok} ({time.time()-t:.0f}s)')
