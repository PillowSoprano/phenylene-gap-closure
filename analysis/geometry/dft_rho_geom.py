import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import json, numpy as np, networkx as nx, time, sys
exec(open(_os.path.join(_HERE,'run.py')).read().split("cases=")[0])
def tR(R): return 1+(1.397-R)/0.2756
d=json.load(open(_os.path.join(_HERE,'..','dft_isomers.json')))
for m in sys.argv[1:]:
    for st in ('S0',):
        G=nx.Graph(); Rb={}
        for a,b,R in d[m][st]['CC_bonds']: G.add_edge(a,b); Rb[frozenset((a,b))]=R
        T,*_=Tmat(G,{e:tR(R) for e,R in Rb.items()}); T1,*_=Tmat(G,{e:1.0 for e in Rb})
        t=time.time(); print(m,st,'rho_geom %.4f'%rho(T),'rho_eq %.4f'%rho(T1),'%.0fs'%(time.time()-t),flush=True)
