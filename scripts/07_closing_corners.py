# All singular corners of the box at eps = rho, by meet-in-the-middle: det(corner) = A[c1] . B[c2]
import sys, json, itertools, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from gapclosure import phenylene, triangular4, benzocyclobutadiene, bipartite, kekule_terms
mols={'benzocyclobutadiene':benzocyclobutadiene(),'biphenylene':phenylene([]),'linear[3]':phenylene(['p']),'angular[3]':phenylene(['L']),
 'linear[4]':phenylene(['p','p']),'zigzag[4]':phenylene(['L','L']),'helical[4]':phenylene(['L','R']),
 'bent[4]':phenylene(['p','L']),'triangular[4]':triangular4()}
rhos=json.load(open(ROOT/'results'/'thresholds.json'))
res={}
for name,g in mols.items():
    t=time.time(); rho=rhos[name]['rho']
    T,L,R,bix=bipartite(g); terms=kekule_terms(T); m=len(bix); eix={e:k for k,e in enumerate(bix)}
    TE=[(s,[eix[(i,p[i])] for i in range(len(T))]) for s,p in terms]
    h=m//2; m2=m-h
    def side(lo,hi):
        k=hi-lo; idx=np.arange(1<<k); S=((idx[:,None]>>np.arange(k))&1)*2-1
        cols=[]
        for s,ix in TE:
            sub=[e-lo for e in ix if lo<=e<hi]
            cols.append(np.prod(1+rho*S[:,sub],axis=1) if sub else np.ones(1<<k))
        return np.array(cols).T
    A=side(0,h)*np.array([s for s,_ in TE]); B=side(h,m)
    sing=[]; second=np.inf; mn=np.inf; ch=1<<10
    for a0 in range(0,A.shape[0],ch):
        D=A[a0:a0+ch]@B.T
        mn=min(mn,D.min())
        z=np.abs(D)<1e-9
        second=min(second,np.abs(D[~z]).min())
        for i,j in zip(*np.nonzero(z)): sing.append((a0+i) | (j<<h))
    S=[np.array([((c>>e)&1)*2-1 for e in range(m)]) for c in sing]
    adj=[(a,b) for a,b in itertools.combinations(range(len(sing)),2) if np.sum(S[a]!=S[b])==1]
    def is_yz(Sv):
        for ring in g['rings']:
            pr=1
            for u,v in zip(ring,ring[1:]+ring[:1]):
                k=[kk for kk,e in enumerate(g['edges']) if set(e)=={u,v}][0]; pr*=Sv[k]
            if pr!=1: return False
        return True
    free=sorted(set(e for a,b in adj for e in np.nonzero(S[a]!=S[b])[0])) if adj else []
    res[name]=dict(m=m,n_sing=len(sing),adjacent_pairs=len(adj),yz=[bool(is_yz(s)) for s in S],
                   min_det=float(mn),smallest_nonzero=float(second),free_bonds=[g['edges'][int(k)] for k in free])
    print(name,res[name],f'{time.time()-t:.0f}s',flush=True)
json.dump(res,open(ROOT/'results'/'closing_corners.json','w'),indent=1)
