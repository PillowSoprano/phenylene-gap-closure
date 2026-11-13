"""Compare rho with gap-based predictors of the adiabatic T1 energy (Table 1, Section 6).
Hueckel-SSH model with Kwapisz-Stolarczyk parameters, as in beyond_huckel/t1_scan.py.
Predictors computed at the relaxed S0 geometry:
  g0        vertical Hueckel gap
  g0^2/P    local closure energy with an isotropic Hessian F ~ fI  (P = sum (pH-pL)^2)
  g0-sP     first-order relaxation with elastic Hessian only (F = K I)
  g0-s/2    first-order relaxation with full S0 Hessian: E ~ g0 - 1/2 gamma^T F^-1 gamma
  Uloc      local closure energy g0^2 / (2 gamma^T F^-1 gamma)
  1D        one-dimensional model, Eq. (9), with Uloc:  g0 (1 - g0/(4 Uloc))
"""
import json, numpy as np, collections, sys

def _rank(v):
    v=np.asarray(v,float); o=np.argsort(v,kind='mergesort'); r=np.empty(len(v)); r[o]=np.arange(len(v))
    for u in np.unique(v):
        i=np.where(v==u)[0]
        if len(i)>1: r[i]=r[i].mean()
    return r
class _S:
    def __init__(s,c): s.correlation=c
def spearmanr(a,b): return _S(np.corrcoef(_rank(a),_rank(b))[0,1])
BH='../beyond_huckel/'
R0,x,y,Rref=1.523,0.189,0.2756,1.397
K=2/(x*y)
def t_of_R(R): return 1.0+(Rref-R)/y

def analyse(d):
    cp=d['cp']; E=[(a,b) for a,b,_ in cp]; m=len(E)
    n=6*d['N']; nh=n//2
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    def eig(R):
        M=np.zeros((n,n)); t=t_of_R(R); M[ia,ib]=-t; M[ib,ia]=-t
        return np.linalg.eigh(M)
    R=np.full(m,Rref)
    for it in range(5000):
        e,C=eig(R); p=2*np.sum(C[ia,:nh]*C[ib,:nh],axis=1); Rn=R0-x*p
        if np.max(abs(Rn-R))<1e-12: break
        R=0.6*R+0.4*Rn
    e,C=eig(R)
    g0=e[nh]-e[nh-1]
    pH=C[ia,nh-1]*C[ib,nh-1]; pL=C[ia,nh]*C[ib,nh]
    dg_dt=2*(pH-pL)            # d(eL-eH)/dt_e, since H_ab=-t
    gam=-dg_dt/y               # d gap / d R_e
    P=np.sum((pH-pL)**2)
    # S0 Hessian wrt t (electronic): 2 * sum_i occ sum_a virt 2 V^e_ia V^f_ia/(e_i-e_a)
    Co=C[:,:nh]; Cv=C[:,nh:]; eo=e[:nh]; ev=e[nh:]
    V=-(Co[ia][:,:,None]*Cv[ib][:,None,:]+Co[ib][:,:,None]*Cv[ia][:,None,:])   # m x occ x virt
    den=eo[:,None]-ev[None,:]
    Htt=4*np.einsum('eia,fia->ef',V,V/den[None])
    F=Htt/y**2+K*np.eye(m)
    s=gam@np.linalg.solve(F,gam)
    sK=gam@gam/K
    return dict(g0=g0,g2P=g0**2/P,lin_K=g0-0.5*sK,lin_F=g0-0.5*s,Uloc=g0**2/(2*s),oneD=g0*(1-g0/(4*(g0**2/(2*s)))),minF=np.linalg.eigvalsh(F).min())

D=json.load(open(BH+'patterns_N2-9.json'))
T=json.load(open(BH+'t1_scan.json')); T={r['id']:r for r in T}
models={'HSSH':{k:v['ET1_adiabatic'] for k,v in T.items()}}
for tag,f in [('PPP-UHF Ohno',None),('ROHF Ohno','ppp_rohf_ohno_N2-9.json'),('ROHF MN','ppp_rohf_mn_N2-9.json'),('ROHF U=8','ppp_rohf_u8_N2-9.json')]:
    if f: models[tag]={r['id']:r['ET1'] for r in json.load(open(BH+f))}
rows=[]
for d in D:
    a=analyse(d); a.update(id=d['id'],N=d['N'],rho=d['rho']); rows.append(a)
json.dump(rows,open('competitors.json','w'),indent=1)
preds=['rho','g0','g2P','lin_K','lin_F','Uloc','oneD']
print('min eigenvalue of S0 Hessian over set: %.3f'%min(r['minF'] for r in rows))
for N_set,lab in [([8],'N=8'),([9],'N=9'),(list(range(4,10)),'4<=N<=9')]:
    sub=[r for r in rows if r['N'] in N_set]
    print('\n%s (%d classes)'%(lab,len(sub)))
    print('%-14s'%'model'+''.join('%9s'%p for p in preds))
    for mname,ME in models.items():
        y_=np.array([ME[r['id']] for r in sub])
        out=[]
        for p in preds:
            v=np.array([r[p] for r in sub]); out.append(spearmanr(v,y_).correlation)
        print('%-14s'%mname+''.join('%9.3f'%c for c in out))
