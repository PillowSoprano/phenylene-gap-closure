import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import json, numpy as np
d=json.load(open(_os.path.join(_HERE,'..','dft_isomers.json')))
def ops(X):
    c=X.mean(0); Y=X-c; u,s,vt=np.linalg.svd(Y); ax=vt  # principal axes (rows)
    out={}
    for name,diag in (('inv',[-1,-1,-1]),('C2_a',[1,-1,-1]),('C2_b',[-1,1,-1]),('C2_c',[-1,-1,1]),('m_a',[-1,1,1]),('m_b',[1,-1,1]),('m_c',[1,1,-1])):
        M=ax.T@np.diag(diag)@ax; out[name]=(Y@M.T)+c
    return out
for m in d:
    S=d[m]['S0']; nC=sum(1 for s,_ in S['geometry'] if s=='C')
    X=np.array([c for s,c in S['geometry'][:nC]])
    best=[]
    for name,Z in ops(X).items():
        D=np.linalg.norm(Z[:,None]-X[None],axis=2); perm=D.argmin(1); err=D.min(1).max()
        if err<0.05 and len(set(perm))==nC and not np.all(perm==np.arange(nC)): best.append((name,err,perm))
    for st in [k for k in d[m] if k.startswith('T1')]+['S0']:
        B={frozenset((a,b)):R for a,b,R in d[m][st]['CC_bonds']}
        for name,err,perm in best:
            dev=max(abs(R-B[frozenset((perm[a],perm[b]))]) for e,R in B.items() for a,b in [tuple(e)])
            print(f'{m:12s}{st:8s}{name:5s} S0-match {err:.3f} A  max |dR| partner bonds {dev:.4f} A')
