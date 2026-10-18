import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import json, time, numpy as np
from gapclosure import *
mols={'benzene':hexagon(),'naphthalene':naphthalene(),'benzocyclobutadiene':benzocyclobutadiene(),
 'biphenylene':phenylene([]),'linear[3]':phenylene(['p']),'angular[3]':phenylene(['L']),
 'linear[4]':phenylene(['p','p']),'zigzag[4]':phenylene(['L','L']),'helical[4]':phenylene(['L','R']),
 'bent[4]':phenylene(['p','L']),'triangular[4]':triangular4()}
out={}
for k,g in mols.items():
    check_geometry(g); T,*_=bipartite(g); t=time.time()
    r,yz=rho_rohn(T); K,ASC,neg=structure_counts(T)
    sv=np.linalg.svd(T,compute_uv=False)
    loc=rho_local(T,restarts=60)
    out[k]=dict(n_atoms=g['n'],K=K,ASC=int(ASC),neg=neg,gap=2*sv[-1],rho=r,rho_local=loc,
               normbound=sv[-1]/sv[0],secs=time.time()-t)
    print(k,out[k],flush=True)
json.dump(out,open(ROOT/'results'/'thresholds.json','w'),indent=1)
