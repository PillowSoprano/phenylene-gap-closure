import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import numpy as np, networkx as nx
exec(open(_os.path.join(_HERE,'run2.py')).read().split("# one-minority closed form")[0])
from scipy.optimize import brentq
R0,x,y,Rref=1.523,0.189,0.2756,1.397; K_=2/(x*y)
def tR(R): return 1+(Rref-R)/y
def relax(G,occ):
    E=[frozenset(e) for e in G.edges()]; nodes=sorted(G); idx={v:i for i,v in enumerate(nodes)}; n=len(nodes); nh=n//2
    R={e:Rref for e in E}
    for it in range(5000):
        H=np.zeros((n,n))
        for e in E:
            i,j=tuple(e); H[idx[i],idx[j]]=H[idx[j],idx[i]]=-tR(R[e])
        w,v=np.linalg.eigh(H); occv=np.zeros(n); occv[:nh]=2
        if occ=='T1': occv[nh-1]=1; occv[nh]=1
        p={e:float(np.sum(occv*v[idx[tuple(e)[0]]]*v[idx[tuple(e)[1]]])) for e in E}
        Rn={e:R0-x*p[e] for e in E}
        if max(abs(Rn[e]-R[e]) for e in E)<1e-12: break
        R={e:0.6*R[e]+0.4*Rn[e] for e in E}
    return R,p
G=phenylene2([]); w1={frozenset(e):1.0 for e in G.edges()}
Kk=kekule(G,w1); Mm=[M for s,M in Kk if s<0][0]; Mp=[M for s,M in Kk if s>0]
# Pauling bond orders
Kall=len(Kk); pa={e:sum(e in M for s,M in Kk)/Kall for e in w1}
for lab,occ in [('S0','S0'),('T1','T1')]:
    R,p=relax(G,occ); t={e:tR(R[e]) for e in R}
    cs=[(len(M-Mm),np.prod([t[e] for e in M-Mm])/np.prod([t[e] for e in Mm-M])) for M in Mp]
    q=brentq(lambda q: sum(c*q**k for k,c in cs)-1,1e-9,1-1e-12) if sum(c for k,c in cs)>1 else None
    rho=(1-q)/(1+q) if q else 0
    # first-order: log c_a ~ (x/y) [sum_{Ma\M-} p - sum_{M-\Ma} p] relative to pref
    lin=[(x/y)*(sum(p[e] for e in M-Mm)-sum(p[e] for e in Mm-M)) for M in Mp]
    print(lab,'c_a',[round(c,3) for k,c in cs],' ln c_a',[round(np.log(c),3) for k,c in cs],' first-order',[round(l,3) for l in lin],' rho(T0=relaxed)=%.4f'%rho)
bridge=[e for e in Mm if all(e not in M for M in Mp)]
R,p=relax(G,'S0'); print('bridging bonds: Huckel bond order %.3f, Pauling %.2f, R=%.3f A'%(p[bridge[0]],pa[bridge[0]],R[bridge[0]]))
print('M- bonds Pauling orders',sorted(round(pa[e],2) for e in Mm),' mean over all bonds',round(np.mean(list(pa.values())),3))
for M in Mp: print('  M_a Pauling mean of its double bonds %.3f'%np.mean([pa[e] for e in M]))
print('  M- Pauling mean %.3f'%np.mean([pa[e] for e in Mm]))
