"""Checks for competitors.py: finite-difference test of gap gradient and S0 Hessian; bootstrap CI of Spearman differences."""
import json, numpy as np
exec(open('competitors.py').read().split("D=json.load")[0])
D=json.load(open(BH+'patterns_N2-9.json'))
d=[q for q in D if q['N']==3][0]
cp=d['cp']; E=[(a,b) for a,b,_ in cp]; m=len(E); n=6*d['N']; nh=n//2
ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
def ES0(R):
    M=np.zeros((n,n)); t=t_of_R(R); M[ia,ib]=-t; M[ib,ia]=-t
    e=np.linalg.eigvalsh(M); return 2*e[:nh].sum()+K/2*np.sum((R-R0)**2), e[nh]-e[nh-1]
R=np.full(m,Rref)
for it in range(5000):
    M=np.zeros((n,n)); t=t_of_R(R); M[ia,ib]=-t; M[ib,ia]=-t
    e,C=np.linalg.eigh(M); p=2*np.sum(C[ia,:nh]*C[ib,:nh],axis=1); Rn=R0-x*p
    if np.max(abs(Rn-R))<1e-13: break
    R=0.6*R+0.4*Rn
h=1e-4; Hn=np.zeros((m,m)); gn=np.zeros(m)
for i in range(m):
    ei=np.zeros(m); ei[i]=h
    gn[i]=(ES0(R+ei)[1]-ES0(R-ei)[1])/(2*h)
    for j in range(m):
        ej=np.zeros(m); ej[j]=h
        Hn[i,j]=(ES0(R+ei+ej)[0]-ES0(R+ei-ej)[0]-ES0(R-ei+ej)[0]+ES0(R-ei-ej)[0])/(4*h*h)
# analytic, same formulas as competitors.py
pH=C[ia,nh-1]*C[ib,nh-1]; pL=C[ia,nh]*C[ib,nh]; gam=-2*(pH-pL)/y
Co=C[:,:nh]; Cv=C[:,nh:]; eo=e[:nh]; ev=e[nh:]
V=-(Co[ia][:,:,None]*Cv[ib][:,None,:]+Co[ib][:,:,None]*Cv[ia][:,None,:])
Htt=4*np.einsum('eia,fia->ef',V,V/(eo[:,None]-ev[None,:])[None]); F=Htt/y**2+K*np.eye(m)
print('molecule',d['id'],'max |gamma_an-gamma_fd| = %.2e (scale %.2f)'%(abs(gam-gn).max(),abs(gn).max()))
print('max |F_an-F_fd| = %.2e (scale %.2f)'%(abs(F-Hn).max(),abs(Hn).max()))
# bootstrap
rows={r['id']:r for r in json.load(open('competitors.json'))}
BHm={'HSSH':{r['id']:r['ET1_adiabatic'] for r in json.load(open(BH+'t1_scan.json'))},
     'ROHF Ohno':{r['id']:r['ET1'] for r in json.load(open(BH+'ppp_rohf_ohno_N2-9.json'))}}
rng=np.random.default_rng(0)
for mod,ME in BHm.items():
    ids=[k for k,r in rows.items() if r['N']==9]
    y_=np.array([ME[k] for k in ids]); P={p:np.array([rows[k][p] for k in ids]) for p in ['rho','g0','lin_K']}
    for comp in ['g0','lin_K']:
        diffs=[]
        for b in range(2000):
            s=rng.integers(0,len(ids),len(ids))
            diffs.append(spearmanr(P['rho'][s],y_[s]).correlation-spearmanr(P[comp][s],y_[s]).correlation)
        lo,hi=np.percentile(diffs,[2.5,97.5])
        print('%s N=9: rho minus %s = %.3f, 95%% bootstrap CI [%.3f, %.3f]'%(mod,comp,np.mean(diffs),lo,hi))
