import sys, json, itertools, time; from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import *
def bounds(g):
    T,L,R,bix=bipartite(g); n=len(T); Ti=np.linalg.inv(T); sv=np.linalg.svd(T,compute_uv=False)
    L1=sv[-1]/sv[0]; L2=1/max(abs(np.linalg.eigvals(np.abs(Ti)@T)))
    terms=kekule_terms(T); sets=[(s,frozenset((i,p[i]) for i in range(n))) for s,p in terms]
    U=1.0
    for s0,M0 in sets:
        if s0>0: continue
        coef=np.zeros(n+1)
        for s,M in sets: coef[len(M-M0)]+=s
        r=[x.real for x in np.roots(coef[::-1]) if abs(x.imag)<1e-9 and 0<x.real<1]
        q=max(r); U=min(U,(1-q)/(1+q))
    K=len(sets); ASC=abs(sum(s for s,_ in sets))
    return dict(K=K,ASC=ASC,gap=2*sv[-1],L_gap=L1,L_BS=L2,U_Kek=U)
res={}
for N in (2,3,4,5,6):
    t=time.time(); cls={}
    for ch in itertools.product('pLR',repeat=N-2):
        seq=''.join('L' if c=='p' else 'A' for c in ch); key=min(seq,seq[::-1])
        cls.setdefault(key,[]).append(ch)
    for key,chs in sorted(cls.items()):
        B=[bounds(phenylene(list(c))) for c in chs]
        for f in ('K','ASC','L_BS','U_Kek'):
            vals=[b[f] for b in B]; assert max(vals)-min(vals)<1e-9 if f in ('K','ASC','U_Kek') else True, (key,f,vals)
        b=B[0]; gaps=[x['gap'] for x in B]; lg=[x['L_gap'] for x in B]; lbs=[x['L_BS'] for x in B]
        res[f'{N}:{key or "-"}']=dict(N=N,seq=key,isomers=len(chs),K=b['K'],ASC=b['ASC'],gap_min=min(gaps),gap_max=max(gaps),
             L_gap_min=min(lg),L_gap_max=max(lg),L_BS_min=min(lbs),L_BS_max=max(lbs),U_Kek=b['U_Kek'])
        r=res[f'{N}:{key or "-"}']
        print(f"[{N}] {key or '-':5s} n={len(chs):2d} K={b['K']} ASC={b['ASC']} gap {min(gaps):.3f}-{max(gaps):.3f} lower={max(max(lg),max(lbs)):.4f} (gap {max(lg):.4f}, BS {max(lbs):.4f}) upper={b['U_Kek']:.4f}",flush=True)
    print(f'  ({time.time()-t:.0f}s)')
json.dump(json.loads(json.dumps(res,default=float)),open(ROOT/'results'/'bounds_by_LA_sequence.json','w'),indent=1)
