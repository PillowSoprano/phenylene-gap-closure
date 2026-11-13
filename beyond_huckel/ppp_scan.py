"""HSSH geometry + PPP (Ohno) electrons. S0: RHF. T1: UHF (S_z=1), several starts, lowest kept."""
import json, sys, numpy as np, collections, time
from multiprocessing import Pool
import os
NTRIAL=int(os.environ.get('NTRIAL','3'))
t0=float(os.environ.get('T0','2.4')); U=float(os.environ.get('U','11.13'))
POT=os.environ.get('POT','ohno'); S0MODE=os.environ.get('S0MODE','RHF'); T1MODE=os.environ.get('T1MODE','UHF'); R0,x,y,Rref=1.523,0.189,0.2756,1.397
K=2*t0/(x*y)                      # eV/A^2, reproduces R = R0 - x p
def t_of_R(R): return t0*(1.0+(Rref-R)/y)
def setup(d,coords):
    cp=d['cp']; E=[(a,b) for a,b,_ in cp]; s=np.array([c for *_,c in cp],float)
    X=np.array(coords)*1.40; n=len(X)
    r=np.linalg.norm(X[:,None]-X[None],axis=2)
    if POT=='ohno': V=U/np.sqrt(1+(U*r/14.397)**2)       # Ohno, diagonal = U
    else: V=14.397/(r+14.397/U)                          # Mataga-Nishimoto, diagonal = U
    return E,s,V,n
def scf(E,V,n,R,na,nb,D0=None,maxit=400,tol=1e-8):
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    h=np.zeros((n,n)); tt=t_of_R(R); h[ia,ib]=-tt; h[ib,ia]=-tt
    if D0 is None:
        e,C=np.linalg.eigh(h); Da=C[:,:na]@C[:,:na].T; Db=C[:,:nb]@C[:,:nb].T
    else: Da,Db=D0
    Voff=V.copy(); np.fill_diagonal(Voff,0)
    for it in range(maxit):
        nt=np.diag(Da)+np.diag(Db)
        J=np.diag(Voff@(nt-1))
        Fa=h+J+np.diag(U*np.diag(Db))-Voff*Da
        Fb=h+J+np.diag(U*np.diag(Da))-Voff*Db
        ea,Ca=np.linalg.eigh(Fa); eb,Cb=np.linalg.eigh(Fb)
        Dan=Ca[:,:na]@Ca[:,:na].T; Dbn=Cb[:,:nb]@Cb[:,:nb].T
        err=max(abs(Dan-Da).max(),abs(Dbn-Db).max())
        Da=0.5*Da+0.5*Dan; Db=0.5*Db+0.5*Dbn
        if err<tol: break
    Eel=0.5*(np.sum(Da*(h+Fa))+np.sum(Db*(h+Fb)))
    p=(Da+Db)[ia,ib]
    gap=None
    if na==nb: gap=ea[na]-ea[na-1]
    return Eel,p,(Da,Db),gap
def relax(E,V,n,na,nb,D0=None):
    R=np.full(len(E),Rref); D=D0; tol=1e-4
    for g in range(300):
        Eel,p,D,gap=scf(E,V,n,R,na,nb,D,tol=tol)
        Rn=R0-x*p
        dmax=np.max(abs(Rn-R))
        if dmax<2e-6 and tol<=1e-8: break
        tol=max(1e-8,min(1e-4,0.01*dmax))
        R=0.3*R+0.7*Rn
    Etot=Eel+K/2*np.sum((R-R0)**2)
    return R,Etot,D,gap
def scf_rohf(E,V,n,R,nc,C0=None,maxit=500,tol=1e-8,mix=0.5):
    """High-spin triplet ROHF: nc doubly occupied orbitals, two singly occupied (alpha).
    Guest-Saunders effective Fock, density damping on the orbital update."""
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    h=np.zeros((n,n)); tt=t_of_R(R); h[ia,ib]=-tt; h[ib,ia]=-tt
    if C0 is None: _,C=np.linalg.eigh(h)
    else: C=C0
    Voff=V.copy(); np.fill_diagonal(Voff,0)
    Dc=C[:,:nc]@C[:,:nc].T; Do=C[:,nc:nc+2]@C[:,nc:nc+2].T
    for it in range(maxit):
        Da=Dc+Do; Db=Dc
        nt=np.diag(Da)+np.diag(Db); J=np.diag(Voff@(nt-1))
        Fa=h+J+np.diag(U*np.diag(Db))-Voff*Da
        Fb=h+J+np.diag(U*np.diag(Da))-Voff*Db
        Fm=C.T@((Fa+Fb)/2)@C; Fam=C.T@Fa@C; Fbm=C.T@Fb@C
        c=slice(0,nc); o=slice(nc,nc+2); v=slice(nc+2,n)
        Fe=Fm.copy()
        Fe[c,o]=Fbm[c,o]; Fe[o,c]=Fbm[o,c]
        Fe[o,v]=Fam[o,v]; Fe[v,o]=Fam[v,o]
        e,U_=np.linalg.eigh(Fe); Cn=C@U_
        Dcn=Cn[:,:nc]@Cn[:,:nc].T; Don=Cn[:,nc:nc+2]@Cn[:,nc:nc+2].T
        err=max(abs(Dcn-Dc).max(),abs(Don-Do).max())
        if not (np.all(np.isfinite(Dcn)) and np.all(np.isfinite(Don))): raise FloatingPointError('rohf diverged')
        C=Cn; Dc=(1-mix)*Dc+mix*Dcn; Do=(1-mix)*Do+mix*Don
        # re-orthonormalize the damped densities through their eigenvectors
        w,Q=np.linalg.eigh(Dc); Cc=Q[:,-nc:]
        P=np.eye(n)-Cc@Cc.T; w2,Q2=np.linalg.eigh(P@Do@P); Co=Q2[:,-2:]
        Rest=np.eye(n)-Cc@Cc.T-Co@Co.T; w3,Q3=np.linalg.eigh(Rest); Cv=Q3[:,-(n-nc-2):]
        C=np.hstack([Cc,Co,Cv]); Dc=Cc@Cc.T; Do=Co@Co.T
        if err<tol: break
    Da=Dc+Do; Db=Dc
    nt=np.diag(Da)+np.diag(Db); J=np.diag(Voff@(nt-1))
    Fa=h+J+np.diag(U*np.diag(Db))-Voff*Da; Fb=h+J+np.diag(U*np.diag(Da))-Voff*Db
    Eel=0.5*(np.sum(Da*(h+Fa))+np.sum(Db*(h+Fb)))
    return Eel,(Da+Db)[ia,ib],C,(Da,Db)
def relax_rohf(E,V,n,nc,mix=0.5):
    R=np.full(len(E),Rref); C=None; tol=1e-4
    for g in range(300):
        Eel,p,C,D=scf_rohf(E,V,n,R,nc,C,tol=tol,mix=mix)
        Rn=R0-x*p; dmax=np.max(abs(Rn-R))
        if dmax<2e-6 and tol<=1e-8: break
        tol=max(1e-8,min(1e-4,0.01*dmax)); R=0.3*R+0.7*Rn
    return R,Eel+K/2*np.sum((R-R0)**2),D,None
def job(args):
    d,coords=args
    E,s,V,n=setup(d,coords); nh=n//2
    RS,ES,DS,gapS=relax(E,V,n,nh,nh)
    bs_lower=0.0
    if S0MODE=='BS':
        # broken-symmetry singlet: alpha biased to one sublattice, beta to the other
        col=np.zeros(n); col[0]=1; adj={i:[] for i in range(n)}
        for a,b in E: adj[a].append(b); adj[b].append(a)
        seen={0}; st=[0]
        while st:
            u=st.pop()
            for v in adj[u]:
                if v not in seen: col[v]=1-col[u]; seen.add(v); st.append(v)
        sgn=2*col-1
        D0=(DS[0]+np.diag(0.3*sgn),DS[1]-np.diag(0.3*sgn))
        RB,EB,DB,_=relax(E,V,n,nh,nh,D0)
        if EB<ES-1e-6: bs_lower=ES-EB; RS,ES,DS=RB,EB,DB
    best=None; rng=np.random.default_rng(1)
    if T1MODE=='ROHF':
        for mix in (0.5,0.2,0.05):
            try:
                with np.errstate(all='raise'):
                    RT,ET,DT,_=relax_rohf(E,V,n,nh-1,mix); best=(RT,ET,DT); break
            except (np.linalg.LinAlgError,FloatingPointError):
                continue
        if best is None:
            return dict(id=d['id'],N=d['N'],la=d['la_tree'],rho=d['rho'],ET1=None,failed=True)
    for trial in (range(NTRIAL) if T1MODE!='ROHF' else []):
        if trial==0: D0=None
        else:
            # localized start: bias alpha excess on a random half of the atoms
            w=rng.normal(size=n); w=np.exp(1.5*w)
            Da=DS[0]+0.0; Db=DS[1]+0.0
            Da=Da+np.diag(w/w.sum())*1.0; Db=Db-np.diag(w/w.sum())*1.0
            D0=(Da,Db)
        RT,ET,DT,_=relax(E,V,n,nh+1,nh-1,D0)
        if best is None or ET<best[1]-1e-7: best=(RT,ET,DT)
    RT,ET,DT=best
    dR=RT-RS
    spin=np.diag(DT[0])-np.diag(DT[1])
    return dict(id=d['id'],N=d['N'],la=d['la_tree'],rho=d['rho'],ET1=ET-ES,gapS0_ppp=float(gapS),bs_lowering=float(bs_lower),
                agree=float(np.mean(np.sign(-dR)==s)),agree_w=float(np.sum(np.abs(dR)*(np.sign(-dR)==s))/np.sum(np.abs(dR))),
                corr=float(np.corrcoef(-dR,s)[0,1]),spin=spin.tolist(),dR=dR.tolist())
if __name__=='__main__':
    D=json.load(open('patterns_N2-9.json')); C=json.load(open('coords_N2-9.json'))
    Nmax=int(sys.argv[1]) if len(sys.argv)>1 else 9
    Nmin=int(sys.argv[2]) if len(sys.argv)>2 else 2
    D=[d for d in D if Nmin<=d['N']<=Nmax]
    tag=os.environ.get('TAG','base')
    fn='ppp_%s_N%d-%d.json'%(tag,Nmin,Nmax)
    t=time.time(); out=[]
    with Pool(int(os.environ.get('NP','2'))) as p:
        for k,r in enumerate(p.imap_unordered(job,[(d,C[d['id']]) for d in D],chunksize=2)):
            out.append(r)
            if (k+1)%10==0 or k+1==len(D):
                json.dump(out,open(fn,'w')); print('%s %d/%d %.0fs'%(tag,k+1,len(D),time.time()-t),flush=True)
    print('done',fn)
