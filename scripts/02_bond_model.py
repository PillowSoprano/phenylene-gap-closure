import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, json, time
from gapclosure import *
mols={'benzocyclobutadiene':benzocyclobutadiene(),'biphenylene':phenylene([]),
 'linear[3]':phenylene(['p']),'angular[3]':phenylene(['L']),
 'linear[4]':phenylene(['p','p']),'zigzag[4]':phenylene(['L','L']),'helical[4]':phenylene(['L','R']),
 'bent[4]':phenylene(['p','L']),'triangular[4]':triangular4()}
out={}
for name,g in mols.items():
    T,*_=bipartite(g); r,(y,z)=rho_rohn(T)
    out[name]={'rho':r,'y':y.tolist(),'z':z.tolist()}
    for K in (21.0,43.0,60.0):
        bm=BondModel(g,K=K); rng=np.random.default_rng(7)
        dw=bm.t0*r*(y[bm.ii]*z[bm.jj])/bm.a
        dr,Er=bm.relax()
        starts=[dw,dr+dw,dr+0.7*dw,dr+1.3*dw,dr+1.6*dw]+[dr+rng.normal(0,0.06,bm.m) for _ in range(20)]
        t=time.time(); c=bm.closure(starts,dr)
        out[name][f'K{int(K)}']=dict(E_relax=Er,dE_close=c.fun-Er,gap_relax=bm.gap(dr)/bm.t0,
            delta_relax=dr.tolist(),delta_close=c.x.tolist(),
            maxrel_close=float(np.abs(bm.a*c.x/bm.t0).max()),
            dE_witness=bm.energy(dw)[0]-Er)
        print(name,K,{k:(round(v,4) if isinstance(v,float) else '...') for k,v in out[name][f'K{int(K)}'].items()},f'{time.time()-t:.0f}s',flush=True)
    json.dump(out,open(ROOT/'results'/'bond_model.json','w'),indent=1)
