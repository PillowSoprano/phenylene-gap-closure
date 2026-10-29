import sys; from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import *
rng=np.random.default_rng(0)
def detw(nodes, edges, w):
    nodes=list(nodes); idx={v:i for i,v in enumerate(nodes)}
    # bipartition by BFS colouring of whole graph passed in col
    rows=[v for v in nodes if col[v]==0]; cols=[v for v in nodes if col[v]==1]
    if len(rows)!=len(cols): return 0.0
    ri={v:i for i,v in enumerate(rows)}; ci={v:i for i,v in enumerate(cols)}
    M=np.zeros((len(rows),len(cols)))
    for (u,v) in edges:
        if u in ri and v in ci: M[ri[u],ci[v]]=w[(u,v)]
        elif v in ri and u in ci: M[ri[v],ci[u]]=w[(u,v)]
    return np.linalg.det(M) if len(rows) else 1.0
for ch in (['L'],['p'],['p','L'],['L','R'],['L','L','p']):
    g=phenylene(ch); E=[tuple(e) for e in g['edges']]
    adj={i:[] for i in range(g['n'])}
    for u,v in E: adj[u].append(v); adj[v].append(u)
    col={0:0}; q=[0]
    for u in q:
        for v in adj[u]:
            if v not in col: col[v]=1-col[u]; q.append(v)
    hexes=[r for r in g['rings'] if len(r)==6]
    for k in range(1,len(hexes)):
        A=set().union(*hexes[:k]); B=set(range(g['n']))-A
        conn=[e for e in E if (e[0] in A)!=(e[1] in A)]
        assert len(conn)==2
        (a,c),(b,d)=[(u,v) if u in A else (v,u) for u,v in conn]
        ratios=[]
        for trial in range(5):
            w={e:rng.uniform(0.5,1.5) for e in E}
            EA=[e for e in E if e[0] in A and e[1] in A]; EB=[e for e in E if e[0] in B and e[1] in B]
            ab=[e for e in EA if set(e)=={a,b}][0]; cd=[e for e in EB if set(e)=={c,d}][0]
            def dA(): w2=dict(w); w2[ab]=0; return (detw(A,EA,w)-detw(A,EA,w2))/w[ab]
            def dB(): w2=dict(w); w2[cd]=0; return (detw(B,EB,w)-detw(B,EB,w2))/w[cd]
            wac=w[[e for e in conn if set(e)=={a,c}][0]]; wbd=w[[e for e in conn if set(e)=={b,d}][0]]
            rhs=detw(A,EA,w)*detw(B,EB,w)-wac*wbd*dA()*dB()
            ratios.append(detw(range(g['n']),E,w)/rhs)
        print(ch,'cut after hexagon',k,'det_G / formula =',np.round(ratios,12))
