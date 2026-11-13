"""Numbers for Sections 6 and 7 (PPP-ROHF Ohno lowest triplet; PT09-00092 excluded as HF-unstable)."""
import json, csv, numpy as np, collections
exec(open('competitors.py').read().split("def analyse")[0])
rows={json.loads(l)['id']:json.loads(l) for l in open('ppp_bs_classes.jsonl')}
comp={r['id']:r for r in json.load(open('competitors.json'))}
T={r['id']:r['ET1_adiabatic'] for r in json.load(open(BH+'t1_scan.json'))}
pg={r['id']:r['gapS0_ppp'] for r in json.load(open(BH+'ppp_rohf_ohno_N2-9.json'))}
EXCL={'PT09-00092'}
preds=['rho','g0','lin_K','lin_F','g2P']
for N in (8,9):
    ids=[k for k in rows if rows[k]['N']==N and k not in EXCL]
    yH=[T[k] for k in ids]; yP=[rows[k]['ET1_bs'] for k in ids]
    print('N=%d (%d)'%(N,len(ids)))
    print('  HSSH ',' '.join('%s %.3f'%(p,spearmanr([comp[k][p] for k in ids],yH).correlation) for p in preds))
    print('  PPP  ',' '.join('%s %.3f'%(p,spearmanr([comp[k][p] for k in ids],yP).correlation) for p in preds),'PPPgap %.3f'%spearmanr([pg[k] for k in ids],yP).correlation)
rng=np.random.default_rng(0)
ids=[k for k in rows if rows[k]['N']==9 and k not in EXCL]
for lab,Y in [('HSSH',np.array([T[k] for k in ids])),('PPP',np.array([rows[k]['ET1_bs'] for k in ids]))]:
    P={p:np.array([comp[k][p] for k in ids]) for p in preds}
    for c in ['g0','lin_K']:
        d=[]
        for b in range(2000):
            s=rng.integers(0,len(ids),len(ids)); d.append(spearmanr(P['rho'][s],Y[s]).correlation-spearmanr(P[c][s],Y[s]).correlation)
        print('N=9 %s rho-%s %.3f CI [%.3f, %.3f]'%(lab,c,np.mean(d),*np.percentile(d,[2.5,97.5])))
print('class level by N (rho vs Hueckel gap):')
for N in range(3,10):
    ids=[k for k in rows if rows[k]['N']==N and k not in EXCL]
    if len(ids)<3: continue
    yH=[T[k] for k in ids]; yP=[rows[k]['ET1_bs'] for k in ids]
    print('  N=%d n=%3d HSSH rho %.3f gap %.3f | PPP rho %.3f gap %.3f'%(N,len(ids),spearmanr([rows[k]['rho'] for k in ids],yH).correlation,spearmanr([comp[k]['g0'] for k in ids],yH).correlation,spearmanr([rows[k]['rho'] for k in ids],yP).correlation,spearmanr([comp[k]['g0'] for k in ids],yP).correlation))
# equal-rho groups (different classes) with lowest PPP triplet
for N in (8,9):
    g=collections.defaultdict(list)
    for k in rows:
        if rows[k]['N']==N and k not in EXCL: g[round(rows[k]['rho'],9)].append(k)
    for v in g.values():
        if len(v)>1:
            print('N=%d equal-rho group of %d classes: rho %.4f, PPP T1 range %.3f eV, HSSH T1 range %.4f |beta|, gap range %.3f'%(N,len(v),rows[v[0]]['rho'],np.ptp([rows[k]['ET1_bs'] for k in v]),np.ptp([T[k] for k in v]),np.ptp([comp[k]['g0'] for k in v])))
    allv=[rows[k]['ET1_bs'] for k in rows if rows[k]['N']==N and k not in EXCL]; print('   N=%d PPP T1 SD %.3f range %.3f'%(N,np.std(allv),np.ptp(allv)))
n4={k:rows[k]['ET1_bs'] for k in rows if rows[k]['N']==4}; print('N=4 PPP lowest T1:',{k:round(v,3) for k,v in n4.items()})
