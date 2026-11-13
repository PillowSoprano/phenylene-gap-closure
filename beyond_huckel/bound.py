"""Rigorous upper bound on the Hueckel-SSH adiabatic T1 energy from the closing corner:
E_T1^ad <= E_S0(R_c) - E_S0(R_S) <= (K/2)|R_c - R_S|^2, R_c realizing the (rescaled) closing corner."""
import json,numpy as np
R0,x,y,Rref=1.523,0.189,0.2756,1.397; K=2/(x*y)
def t_of_R(R): return 1.0+(Rref-R)/y
def R_of_t(t): return Rref-y*(t-1.0)
D=json.load(open('patterns_N2-9.json')); H={r['id']:r for r in json.load(open('t1_scan.json'))}
rows=[]
for d in D:
    cp=d['cp']; E=[(a,b) for a,b,_ in cp]; s=np.array([c for *_,c in cp],float); n=6*d['N']; nh=n//2
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    def ener(R,occ):
        t=t_of_R(R); Hm=np.zeros((n,n)); Hm[ia,ib]=-t; Hm[ib,ia]=-t
        e,C=np.linalg.eigh(Hm); w=np.zeros(n); w[:nh]=2
        if occ: w[nh-1]=1; w[nh]=1
        return float(w@e)+K/2*np.sum((R-R0)**2),C,w,e
    R=np.full(len(E),Rref)
    for it in range(3000):
        _,C,w,e=ener(R,0); p=(C[ia]*C[ib])@w; Rn=R0-x*p
        if abs(Rn-R).max()<1e-12: break
        R=0.5*R+0.5*Rn
    RS=R; ES,_,_,eS=ener(RS,0); tS=t_of_R(RS)
    tc=1+d['rho']*s                     # closing corner (any positive multiple is also singular)
    sc=(tc@tS)/(tc@tc); Rc=R_of_t(sc*tc)
    Ec,_,_,ec=ener(Rc,0)
    gap_c=ec[nh]-ec[nh-1]
    U_exact=Ec-ES                        # E_S0 at the singular corner = E_T1 there
    U_quad=K/2*np.sum((Rc-RS)**2)
    rows.append((d['id'],d['N'],d['rho'],H[d['id']]['ET1_adiabatic'],U_exact,U_quad,gap_c))
A=np.array([r[1:] for r in rows])
print('max gap at corner (should be 0):',A[:,5].max())
print('bound holds (ET1 <= U_exact <= U_quad):',np.all(A[:,2]<=A[:,3]+1e-9),np.all(A[:,3]<=A[:,4]+1e-9))
def rank(a): return np.argsort(np.argsort(a))
for N in range(4,10):
    m=A[:,0]==N
    print('N=%d n=%3d ET1/U_exact mean %.3f [%.3f,%.3f]  spearman(U_exact,ET1)=%.3f spearman(rho,ET1)=%.3f'%(N,m.sum(),np.mean(A[m,2]/A[m,3]),np.min(A[m,2]/A[m,3]),np.max(A[m,2]/A[m,3]),
          np.corrcoef(rank(A[m,3]),rank(A[m,2]))[0,1],np.corrcoef(rank(A[m,1]),rank(A[m,2]))[0,1]))
