"""Compare predictors of the PPP-UHF adiabatic T1 energy. Usage: python3 an_ppp2.py file1.json [file2.json ...]"""
import json,sys,numpy as np
from scipy.stats import spearmanr
H={r['id']:r for r in json.load(open('t1_scan.json'))}
def sp(a,b): return float(spearmanr(a,b).correlation)   # ties (classes of equal rho) get average ranks
for fn in sys.argv[1:]:
    P=[r for r in json.load(open(fn)) if r.get('ET1') is not None]; print('==',fn,len(P),'classes (failed dropped)')
    for N in sorted({r['N'] for r in P}):
        xs=[r for r in P if r['N']==N]
        if len(xs)<4: continue
        E=np.array([r['ET1'] for r in xs]); rho=np.array([r['rho'] for r in xs]); gH=np.array([H[r['id']]['gapS0'] for r in xs]); gP=np.array([r['gapS0_ppp'] for r in xs])
        bs=np.mean([r.get('bs_lowering',0)>0 for r in xs])
        print(' N=%d n=%d  spearman rho %.3f | Huckel gap %.3f | PPP gap %.3f | weighted agree %.3f | BS-lower frac %.2f'%(N,len(xs),sp(rho,E),sp(gH,E),sp(gP,E),np.mean([r['agree_w'] for r in xs]),bs))
