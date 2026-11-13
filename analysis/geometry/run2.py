import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import numpy as np, itertools, networkx as nx, time, sys
sys.path.insert(0,_HERE)
exec(open(_os.path.join(_HERE,'run.py')).read().split("cases=")[0])

def phenylene2(modes):
    G=nx.Graph(); nxt=0
    def add_hex(start):
        h=list(range(start,start+6))
        for i in range(6): G.add_edge(h[i],h[(i+1)%6])
        return h
    h=add_hex(0); nxt=6; attach=(h[3],h[4])
    for mode in list(modes)+['END']:
        h2=add_hex(nxt); nxt+=6
        a,b=attach; G.add_edge(b,h2[0]); G.add_edge(a,h2[1])
        attach={'L':(h2[3],h2[4]),'A':(h2[2],h2[3]),'B':(h2[4],h2[5])}.get(mode)
    return G

# Kekule structures with signs via det expansion
def kekule(G,w):
    T,si,ui,col=Tmat(G,w)
    S={i:v for v,i in si.items()}; U={j:v for v,j in ui.items()}
    n=T.shape[0]; out=[]
    def rec(r,used,perm):
        if r==n:
            p=list(perm); s=1
            # sign of permutation
            seen=[False]*n
            for i in range(n):
                if not seen[i]:
                    j=i;L=0
                    while not seen[j]: seen[j]=True;j=p[j];L+=1
                    if L%2==0: s=-s
            out.append((s,frozenset(frozenset((S[i],U[p[i]])) for i in range(n))))
            return
        for c in range(n):
            if c not in used and T[r,c]!=0:
                rec(r+1,used|{c},perm+[c])
    rec(0,frozenset(),[])
    tot=sum(s for s,_ in out)
    if tot<0: out=[(-s,M) for s,M in out]
    return out

# one-minority closed form for biphenylene / BCB at relaxed couplings
from scipy.optimize import brentq
for name,G in [('BCB',benzocyclobutadiene()),('biphenylene',phenylene2([]))]:
    wr,d=relax(G); w1={frozenset(e):1.0 for e in G.edges()}
    K=kekule(G,w1); Mm=[M for s,M in K if s<0]; assert len(Mm)==1; Mm=Mm[0]
    cs=[]
    for s,M in K:
        if s<0: continue
        num=np.prod([wr[e] for e in M-Mm]); den=np.prod([wr[e] for e in Mm-M])
        cs.append((len(M-Mm),num/den))
    f=lambda q: sum(c*q**k for k,c in cs)-1
    q=brentq(f,1e-9,1-1e-12); r=(1-q)/(1+q)
    Tr,*_=Tmat(G,wr)
    print(name,'c_a (k,c):',[(k,round(c,4)) for k,c in cs],' formula rho=%.5f  enumeration rho=%.5f'%(r,rho(Tr)))
    # which bonds belong only to M-
    only=[tuple(sorted(e)) for e in Mm if all(e not in M for s,M in K if s>0)]
    print('   bonds only in M-:',only,' lengths change dA:',[round(d[frozenset(e)],3) for e in only])

# zigzag vs helical [4]
t=time.time()
for modes in (['A','A'],['A','B']):
    G=phenylene2(modes); w1={frozenset(e):1.0 for e in G.edges()}
    T1,*_=Tmat(G,w1); g1=2*np.linalg.svd(T1,compute_uv=False).min()
    wr,d=relax(G); Tr,*_=Tmat(G,wr); gr=2*np.linalg.svd(Tr,compute_uv=False).min()
    print(modes,'gap_equal=%.3f'%g1,'gap_relaxed=%.3f eV'%gr, flush=True)
    print('   rho_relaxed=%.5f'%rho(Tr),' elapsed %.0fs'%(time.time()-t),flush=True)
