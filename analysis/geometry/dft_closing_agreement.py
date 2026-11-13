import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import json, numpy as np, networkx as nx, itertools, sys
exec(open(_os.path.join(_HERE,'run.py')).read().split("cases=")[0])
def tR(R): return 1+(1.397-R)/0.2756
d=json.load(open(_os.path.join(_HERE,'..','dft_isomers.json')))
m=sys.argv[1]
G=nx.Graph(); R0={}; R1={}
for a,b,R in d[m]['S0']['CC_bonds']: G.add_edge(a,b); R0[frozenset((a,b))]=R
for a,b,R in d[m]['T1_pert']['CC_bonds']: R1[frozenset((a,b))]=R
T,si,ui,col=Tmat(G,{e:tR(R) for e,R in R0.items()})
n=T.shape[0]; Ti=np.linalg.inv(T); signs=np.array(list(itertools.product([-1,1],repeat=n))); res=[]
for z in signs[signs[:,0]==1]:
    M=T[None]*signs[:,:,None]*z[None,None,:]; ev=np.linalg.eigvals(Ti[None]@M)
    re=np.where(np.abs(ev.imag)<1e-9,ev.real,0); k=np.abs(re).max(1)
    for idx in np.where(k>k.max()-1e-9)[0]:
        lam=re[idx][np.abs(re[idx]).argmax()]; res.append((abs(lam),lam,signs[idx].copy(),z.copy()))
best=max(r[0] for r in res); pats=[r for r in res if r[0]>best-1e-7]
print(m,'rho %.4f'%(1/best),'patterns',len(pats))
dt={e:-(R1[e]-R0[e])/0.2756 for e in R0}; tot=sum(abs(v) for v in dt.values())
seen=set()
for _,lam,y,z in pats:
    s={}
    for e in R0:
        i,j=tuple(e)
        if col[i]==1: i,j=j,i
        s[e]=-np.sign(lam)*y[si[i]]*z[ui[j]]
    key=tuple(sorted((tuple(sorted(e)),s[e]) for e in s))
    if key in seen: continue
    seen.add(key)
    agree=[e for e in R0 if np.sign(dt[e])==s[e]]
    print('  pattern: agree %d/%d bonds, weight %.0f%%'%(len(agree),len(R0),100*sum(abs(dt[e]) for e in agree)/tot))
