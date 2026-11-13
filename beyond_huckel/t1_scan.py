"""HSSH: compare the L-infinity closing pattern (dataset) with T1 relaxation and the L2 gap gradient."""
import json, numpy as np, collections
R0,x,y,Rref=1.523,0.189,0.2756,1.397
K=2/(x*y)
def t_of_R(R): return 1.0+(Rref-R)/y
def run(d):
    cp=d['cp']; E=[(a,b) for a,b,_ in cp]; s=np.array([c for *_,c in cp],float)
    n=6*d['N']; nh=n//2
    def H(t):
        M=np.zeros((n,n))
        for (a,b),w in zip(E,t): M[a,b]=M[b,a]=-w
        return M
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    def solve(R,occ):
        e,C=np.linalg.eigh(H(t_of_R(R)))
        w=np.zeros(n); w[:nh]=2
        if occ=='T1': w[nh-1]=1; w[nh]=1
        p=(C[ia,:]*C[ib,:])@w
        return p,e,C,w
    def relax(occ):
        R=np.full(len(E),Rref)
        for it in range(2000):
            p,e,C,w=solve(R,occ); Rn=R0-x*p
            if np.max(abs(Rn-R))<1e-11: break
            R=0.6*R+0.4*Rn
        p,e,C,w=solve(R,occ)
        Etot=float(w@e)+K/2*np.sum((R-R0)**2)
        return R,p,e,C,Etot
    RS,pS,eS,CS,ES=relax('S0'); RT,pT,eT,CT,ET=relax('T1')
    dR=RT-RS
    dp=CS[ia,nh]*CS[ib,nh]-CS[ia,nh-1]*CS[ib,nh-1]
    agreeT=np.mean(np.sign(-dR)==s); agreeL2=np.mean(np.sign(dp)==s)
    wT=np.sum(np.abs(dR)*(np.sign(-dR)==s))/np.sum(np.abs(dR))
    corrT=np.corrcoef(-dR,s)[0,1]; corrL2=np.corrcoef(dp,s)[0,1]; corrTL2=np.corrcoef(-dR,dp)[0,1]
    gapS=eS[nh]-eS[nh-1]
    _,eTS,_,_=solve(RT,'S0'); gapT=eTS[nh]-eTS[nh-1]
    return dict(id=d['id'],N=d['N'],la=d['la_tree'],rho=d['rho'],agreeT=agreeT,agreeT_w=wT,agreeL2=agreeL2,
                corrT=corrT,corrL2=corrL2,corrTL2=corrTL2,gapS0=gapS,gap_at_T1geom=gapT,ET1_adiabatic=ET-ES,
                bad_bonds=[[a,b,int(c),float(-dR[k])] for k,((a,b),c) in enumerate(zip(E,s)) if np.sign(-dR[k])!=c])
D=json.load(open('patterns_N2-9.json'))
out=[run(d) for d in D]
json.dump(out,open('t1_scan.json','w'),indent=1)
by=collections.defaultdict(list)
for r in out: by[r['N']].append(r)
print(' N  n  agreeT(min/mean)  agreeT_w(mean)  agreeL2(mean)  corrT(min/mean)  corrTL2(mean)')
for N in sorted(by):
    xs=by[N]
    print('%2d %3d   %.2f/%.3f        %.3f          %.3f        %.2f/%.3f       %.3f'%(N,len(xs),min(r['agreeT'] for r in xs),np.mean([r['agreeT'] for r in xs]),
          np.mean([r['agreeT_w'] for r in xs]),np.mean([r['agreeL2'] for r in xs]),min(r['corrT'] for r in xs),np.mean([r['corrT'] for r in xs]),np.mean([r['corrTL2'] for r in xs])))
rho=np.array([r['rho'] for r in out]);E1=np.array([r['ET1_adiabatic'] for r in out]);g=np.array([r['gapS0'] for r in out])
print('spearman-ish corr: corr(rho,ET1)=%.3f corr(gapS0,ET1)=%.3f corr(rho,gapS0)=%.3f'%(np.corrcoef(rho,E1)[0,1],np.corrcoef(g,E1)[0,1],np.corrcoef(rho,g)[0,1]))
for N in sorted(by):
    xs=by[N]
    if len(xs)>3:
        r_=np.array([r['rho'] for r in xs]);e_=np.array([r['ET1_adiabatic'] for r in xs]);g_=np.array([r['gapS0'] for r in xs])
        print(' N=%d corr(rho,ET1)=%.3f corr(gap,ET1)=%.3f'%(N,np.corrcoef(r_,e_)[0,1],np.corrcoef(g_,e_)[0,1]))
