"""Broken-symmetry S0 check for classes whose ROHF triplet dropped by >0.2 eV (and PT09-00092)."""
import os, sys, json, numpy as np, warnings
os.environ['S0MODE']='BS'; os.environ['T1MODE']='ROHF'; os.environ['POT']='ohno'; warnings.filterwarnings('ignore')
sys.path.insert(0,'../beyond_huckel'); import ppp_scan as P
from multiprocessing import Pool
rows={json.loads(l)['id']:json.loads(l) for l in open('ppp_bs_classes.jsonl')}
ids=[k for k,r in rows.items() if r['ET1_plain']-r['ET1_bs']>0.2]
D={d['id']:d for d in json.load(open('../beyond_huckel/patterns_N2-9.json'))}; C=json.load(open('../beyond_huckel/coords_N2-9.json'))
def job(k):
    r=P.job((D[k],C[k])); return k,r['bs_lowering']
if __name__=='__main__':
    with Pool(4) as p: res=p.map(job,ids)
    out=[]
    for k,bl in res:
        r=rows[k]; out.append(dict(id=k,ET1_plain=r['ET1_plain'],ET1_bs=r['ET1_bs'],S0_bs_lowering=bl,ET1_bs_vs_bsS0=r['ET1_bs']+bl))
        print('%s  T1 symmetric %.3f  T1 lowest %.3f  S0 BS lowering %.3f  -> T1 lowest vs BS S0 %.3f eV'%(k,r['ET1_plain'],r['ET1_bs'],bl,r['ET1_bs']+bl))
    json.dump(out,open('s0_bs_check.json','w'),indent=1)
