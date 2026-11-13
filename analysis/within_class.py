"""Molecules of the same L/A class share rho exactly (Section 7). Do they share the T1 energy? (Table S7)
Hueckel-SSH (Kwapisz-Stolarczyk), every molecule with N<=9 from the tree-phenylene dataset."""
import sys, csv, json, numpy as np, collections, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'../gapclosure/src')
from gapclosure.skeleton import phenylene_tree
exec(open('competitors.py').read().split("def analyse")[0])   # constants, spearmanr
def hssh(sk):
    E=sk.bonds; n=sk.n; nh=n//2; m=len(E)
    ia=np.array([a for a,_ in E]); ib=np.array([b for _,b in E])
    def solve(R,occ):
        M=np.zeros((n,n)); t=t_of_R(R); M[ia,ib]=-t; M[ib,ia]=-t
        e,C=np.linalg.eigh(M); w=np.zeros(n); w[:nh]=2
        if occ=='T1': w[nh-1]=1; w[nh]=1
        return (C[ia,:]*C[ib,:])@w, e, w
    def relax(occ):
        R=np.full(m,Rref)
        for it in range(5000):
            p,e,w=solve(R,occ); Rn=R0-x*p
            if np.max(abs(Rn-R))<1e-11: break
            R=0.6*R+0.4*Rn
        p,e,w=solve(R,occ)
        return float(w@e)+K/2*np.sum((R-R0)**2), e
    ES,eS=relax('S0'); ET,_=relax('T1')
    return ET-ES, eS[nh]-eS[nh-1]
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 9
cls={r['id']:r for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylenes_N1-12.csv'))}
out=[]
for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylene_molecules_N1-12.csv')):
    N=int(r['N'])
    if N<2 or N>NMAX: continue
    sk=phenylene_tree(json.loads(r['links']))
    et,g0=hssh(sk)
    out.append(dict(id=r['id'],N=N,smiles=r['smiles'],rho=float(cls[r['id']]['rho']),gap_eq=float(r['huckel_gap']),
                    overlaps=r['overlaps_in_plane']=='True',ET1=et,g0=g0))
json.dump(out,open('within_class_N%d.json'%NMAX if False else 'within_class_N'+str(NMAX)+'.json','w'),indent=0)
print('molecules',len(out))
for N in range(4,NMAX+1):
    xs=[o for o in out if o['N']==N]
    by=collections.defaultdict(list)
    for o in xs: by[o['id']].append(o)
    multi=[v for v in by.values() if len(v)>1]
    et=np.array([o['ET1'] for o in xs])
    # share of T1 variance between classes (eta^2)
    grand=et.mean(); ssb=sum(len(v)*(np.mean([o['ET1'] for o in v])-grand)**2 for v in by.values()); sst=((et-grand)**2).sum()
    rng_T=[np.ptp([o['ET1'] for o in v]) for v in multi]; rng_g=[np.ptp([o['g0'] for o in v]) for v in multi]
    rho=np.array([o['rho'] for o in xs]); g=np.array([o['g0'] for o in xs])
    print('N=%d molecules %4d classes %3d multi %3d | eta2(class)=%.3f | within-class T1 range median %.4f max %.4f | gap range median %.4f | total T1 SD %.4f | Spearman over molecules: rho %.3f gap %.3f'%(
        N,len(xs),len(by),len(multi),ssb/sst,np.median(rng_T) if rng_T else 0,max(rng_T) if rng_T else 0,np.median(rng_g) if rng_g else 0,et.std(),spearmanr(rho,et).correlation,spearmanr(g,et).correlation))
z=[o for o in out if o['id']=='PT04-00002']
for o in z: print('AA [4]phenylene', o['smiles'], 'gap_eq %.3f  g0 %.4f  ET1 %.4f |beta|'%(o['gap_eq'],o['g0'],o['ET1']))
