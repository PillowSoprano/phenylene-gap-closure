"""Check that rho depends only on the L/A class: tree procedure on every molecule (actual graph), N given."""
import sys, csv, json, time, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'../gapclosure/src')
from gapclosure.skeleton import phenylene_tree
from gapclosure.threshold import closure_threshold
from fractions import Fraction
Nlo,Nhi=int(sys.argv[1]),int(sys.argv[2])
cls={r['id']:r for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylenes_N1-12.csv'))}
t=time.time(); bad=0; n=0; multi=set()
from collections import Counter
cnt=Counter(r['id'] for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylene_molecules_N1-12.csv')))
for r in csv.DictReader(open('../dataset_tree_phenylenes/tree_phenylene_molecules_N1-12.csv')):
    N=int(r['N'])
    if not (Nlo<=N<=Nhi) or cnt[r['id']]<2: continue
    res=closure_threshold(phenylene_tree(json.loads(r['links'])),method='tree',bits=40)
    lo,hi=Fraction(cls[r['id']]['rho_lower']),Fraction(cls[r['id']]['rho_upper'])
    ok = not (Fraction(res.upper)<lo or Fraction(res.lower)>hi) if hasattr(res,'upper') else abs(res.rho-float(cls[r['id']]['rho']))<1e-9
    n+=1; bad+= (not ok); multi.add(r['id'])
    if not ok: print('MISMATCH',r['id'],r['smiles'],res.rho,cls[r['id']]['rho'])
print('N %d-%d: molecules in multi-molecule classes %d (classes %d), mismatches %d, %.0fs'%(Nlo,Nhi,n,len(multi),bad,time.time()-t))
