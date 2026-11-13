import json, numpy as np, collections
exec(open('competitors.py').read().split("def analyse")[0])
rows={}
for l in open('ppp_bs_classes.jsonl'):
    r=json.loads(l); rows[r['id']]=r
ref={r['id']:r['ET1'] for r in json.load(open(BH+'ppp_rohf_ohno_N2-9.json'))}
comp={r['id']:r for r in json.load(open('competitors.json'))}
ids=list(rows)
plain=np.array([rows[k]['ET1_plain'] for k in ids]); bs=np.array([rows[k]['ET1_bs'] for k in ids]); rf=np.array([ref[k] for k in ids])
print('classes',len(ids),' symmetric start vs ppp_rohf_ohno_N2-9.json max |diff| %.2e'%abs(plain-rf).max())
drop=plain-bs
for N in range(2,10):
    m=np.array([rows[k]['N']==N for k in ids])
    print('N=%d: %3d classes, lowered >0.005 eV: %3d, max lowering %.3f eV'%(N,m.sum(),(drop[m]>0.005).sum(),drop[m].max()))
preds=['rho','g0','lin_K','lin_F','g2P']
for N_set,lab in [([8],'N=8'),([9],'N=9'),(list(range(4,10)),'4<=N<=9')]:
    sub=[k for k in ids if rows[k]['N'] in N_set]
    for tag,E in [('symmetric start',{k:rows[k]['ET1_plain'] for k in sub}),('lowest (broken symmetry)',{k:rows[k]['ET1_bs'] for k in sub})]:
        y_=np.array([E[k] for k in sub])
        print('%-9s %-26s'%(lab,tag)+'  '.join('%s %.3f'%(p,spearmanr([comp[k][p] for k in sub],y_).correlation) for p in preds))
# bootstrap rho - lin_K for N=9 broken symmetry
sub=[k for k in ids if rows[k]['N']==9]; y_=np.array([rows[k]['ET1_bs'] for k in sub]); P={p:np.array([comp[k][p] for k in sub]) for p in ['rho','g0','lin_K']}
rng=np.random.default_rng(0)
for c in ['g0','lin_K']:
    d=[]
    for b in range(2000):
        s=rng.integers(0,len(sub),len(sub)); d.append(spearmanr(P['rho'][s],y_[s]).correlation-spearmanr(P[c][s],y_[s]).correlation)
    print('N=9 broken symmetry: rho minus %s = %.3f, 95%% CI [%.3f, %.3f]'%(c,np.mean(d),*np.percentile(d,[2.5,97.5])))
# PPP gap column
pg={r['id']:r['gapS0_ppp'] for r in json.load(open(BH+'ppp_rohf_ohno_N2-9.json'))}
for N in (8,9):
    sub=[k for k in ids if rows[k]['N']==N]
    print('N=%d PPP gap vs lowest T1: %.3f'%(N,spearmanr([pg[k] for k in sub],[rows[k]['ET1_bs'] for k in sub]).correlation))
