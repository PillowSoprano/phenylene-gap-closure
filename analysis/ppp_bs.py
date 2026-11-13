"""PPP-ROHF triplet with symmetry-broken starting orbitals (exact geometry). Keeps the lowest of the plain start and NSTART perturbed starts."""
import os, sys, json, numpy as np, warnings
os.environ.setdefault('T1MODE','ROHF'); os.environ.setdefault('POT','ohno'); warnings.filterwarnings('ignore')
sys.path.insert(0,'../beyond_huckel')
import ppp_scan as P
NSTART=int(os.environ.get('NSTART','2'))
def relax_rohf_from(E,V,n,nc,C0,mix=0.5):
    R=np.full(len(E),P.Rref); C=C0; tol=1e-4
    for g in range(300):
        Eel,p,C,D=P.scf_rohf(E,V,n,R,nc,C,tol=tol,mix=mix)
        Rn=P.R0-P.x*p; dmax=np.max(abs(Rn-R))
        if dmax<2e-6 and tol<=1e-8: break
        tol=max(1e-8,min(1e-4,0.01*dmax)); R=0.3*R+0.7*Rn
    return R,Eel+P.K/2*np.sum((R-P.R0)**2)
def t1_s0(d,coords,seed=0):
    E,s,V,n=P.setup(d,coords); nh=n//2
    RS,ES,DS,gapS=P.relax(E,V,n,nh,nh)
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    h=np.zeros((n,n)); h[ia,ib]=-P.t0; h[ib,ia]=-P.t0
    rng=np.random.default_rng(seed); ets=[]
    for k in range(NSTART+1):
        pert=np.zeros((n,n)) if k==0 else rng.normal(scale=0.05,size=(n,n))
        _,C0=np.linalg.eigh(h+(pert+pert.T)/2)
        for mix in (0.5,0.2,0.05):
            try:
                with np.errstate(all='raise'):
                    RT,ET=relax_rohf_from(E,V,n,nh-1,C0,mix); ets.append(ET); break
            except (np.linalg.LinAlgError,FloatingPointError): continue
    return min(ets)-ES, (ets[0]-ES) if ets else None, ES
if __name__=='__main__':
    from multiprocessing import Pool
    D=[d for d in json.load(open('../beyond_huckel/patterns_N2-9.json')) if int(sys.argv[1])<=d['N']<=int(sys.argv[2])]
    if len(sys.argv)>4: D=D[int(sys.argv[3]):int(sys.argv[4])]
    C=json.load(open('../beyond_huckel/coords_N2-9.json'))
    def job(d):
        best,plain,ES=t1_s0(d,C[d['id']]); return dict(id=d['id'],N=d['N'],rho=d['rho'],ET1_bs=best,ET1_plain=plain)
    done=set()
    try: done={json.loads(l)['id'] for l in open('ppp_bs_classes.jsonl')}
    except FileNotFoundError: pass
    D=[d for d in D if d['id'] not in done]
    with Pool(4) as p, open('ppp_bs_classes.jsonl','a') as f:
        for r in p.imap_unordered(job,D,chunksize=1):
            f.write(json.dumps(r)+'\n'); f.flush()
    print('finished',len(D))
