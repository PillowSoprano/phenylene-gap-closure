import numpy as np, itertools
# biphenylene: ring1 0..5, ring2 6..11, four-ring 0-1-7-6
E=[(i,(i+1)%6) for i in range(6)]+[(6+i,6+(i+1)%6) for i in range(6)]+[(0,6),(1,7)]
n=12
col=[0]*n
for v in range(6): col[v]=v%2
col[6]=1;col[7]=0
for v in range(6,12): col[v]=(col[6]+(v-6))%2
for a,b in E: assert col[a]!=col[b]
A=[v for v in range(n) if col[v]==0];B=[v for v in range(n) if col[v]==1]
ai={v:i for i,v in enumerate(A)};bi={v:i for i,v in enumerate(B)}
def Tmat(t):
    T=np.zeros((6,6))
    for (a,b),w in zip(E,t):
        if col[a]==1:a,b=b,a
        T[ai[a],bi[b]]=w
    return T
# Kekule structures with signs
adj={v:[] for v in range(n)}
for k,(a,b) in enumerate(E): adj[a].append((b,k));adj[b].append((a,k))
Ms=[]
for perm in itertools.permutations(range(6)):
    ks=[]
    ok=True
    for i,j in enumerate(perm):
        a,b=A[i],B[j]
        k=[k for (w,k) in adj[a] if w==b]
        if not k: ok=False;break
        ks.append(k[0])
    if ok:
        s=np.linalg.det(np.eye(6)[list(perm)])
        Ms.append((round(s),frozenset(ks)))
signs=[s for s,_ in Ms];maj=max(set(signs),key=signs.count)
mino=[M for s,M in Ms if s!=maj]
print('K',len(Ms),'ASC',abs(sum(signs)),'minority',len(mino))
M0=mino[0]
print('minority structure double bonds:',sorted(E[k] for k in M0))
# rho at minority corner
def det_corner(q): return np.linalg.det(Tmat([1.0 if k in M0 else q for k in range(len(E))]))
s0=np.sign(np.linalg.det(Tmat(np.ones(len(E)))))
lo,hi=0.0,1.0
for _ in range(80):
    m=(lo+hi)/2
    if np.sign(det_corner(m))==s0: hi=m
    else: lo=m
q=lo;print('rho(biphenylene)=%.10f q=%.10f'%((1-q)/(1+q),q))
# HSSH bond-order/bond-length iteration
def H(t):
    Hm=np.zeros((n,n))
    for (a,b),w in zip(E,t): Hm[a,b]=Hm[b,a]=-w
    return Hm
def bondorders(t,occ):
    e,C=np.linalg.eigh(H(t))
    p=np.zeros(len(E))
    for k,o in occ.items():
        for m,(a,b) in enumerate(E): p[m]+=o*C[a,k]*C[b,k]
    return p,e,C
R0,x,y,Rref=1.523,0.189,0.2756,1.397
def t_of_R(R): return 1.0+(Rref-R)/y
def relax(occ):
    R=np.full(len(E),Rref)
    for it in range(500):
        p,e,C=bondorders(t_of_R(R),occ)
        Rn=R0-x*p
        if np.max(abs(Rn-R))<1e-10: break
        R=0.5*R+0.5*Rn
    return R,p
S0={k:2 for k in range(6)}; T1={**{k:2 for k in range(5)},5:1,6:1}
R_S,_=relax(S0);R_T,_=relax(T1)
# L2 direction at S0 geometry: p_LUMO - p_HOMO
_,e,C=bondorders(t_of_R(R_S),S0)
dp=np.array([C[a,6]*C[b,6]-C[a,5]*C[b,5] for a,b in E])
print('gap S0 geometry %.4f'%(e[6]-e[5]))
print('%-8s %6s %8s %8s %9s'%('bond','inM0','R_S0','dR_T1','dp(L-H)'))
for k,(a,b) in enumerate(E):
    print('%-8s %6s %8.4f %+8.4f %+9.4f'%(f'{a}-{b}','yes' if k in M0 else '',R_S[k],R_T[k]-R_S[k],dp[k]))
dR=R_T-R_S
corner=np.array([1.0 if k in M0 else -1.0 for k in range(len(E))])
print('sign agreement corner vs -dR_T1:',sum(np.sign(-dR)==corner),'/',len(E))
print('sign agreement corner vs dp:',sum(np.sign(dp)==corner),'/',len(E))
print('corr(-dR,corner)=%.3f corr(dp,corner)=%.3f'%(np.corrcoef(-dR,corner)[0,1],np.corrcoef(dp,corner)[0,1]))
# T1 gap closes? check HOMO-LUMO gap at T1 geometry
_,eT,_=bondorders(t_of_R(R_T),S0);print('gap at T1 geometry %.4f'%(eT[6]-eT[5]))
# relative coupling change at T1 vs rho
tS,tT=t_of_R(R_S),t_of_R(R_T);rel=(tT-tS)/tS;print('max |rel coupling change| S0->T1 = %.4f'%abs(rel).max())
