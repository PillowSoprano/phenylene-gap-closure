import json, numpy as np, collections
exec(open('competitors.py').read().split("def analyse")[0])
out=[json.loads(l) for l in open('within_class_mols.jsonl')]
print('N  mols cls multi | eta2  | within T1 range med/max | within g0 range med | slope dT1/dg0 within | Spearman(mols): rho  g0')
for N in range(4,10):
    xs=[o for o in out if o['N']==N]; by=collections.defaultdict(list)
    for o in xs: by[o['id']].append(o)
    multi=[v for v in by.values() if len(v)>1]
    et=np.array([o['ET1'] for o in xs]); grand=et.mean()
    eta=sum(len(v)*(np.mean([o['ET1'] for o in v])-grand)**2 for v in by.values())/((et-grand)**2).sum()
    dT=np.concatenate([[o['ET1']-np.mean([u['ET1'] for u in v]) for o in v] for v in multi])
    dg=np.concatenate([[o['g0']-np.mean([u['g0'] for u in v]) for o in v] for v in multi])
    slope=(dT@dg)/(dg@dg); r=np.corrcoef(dT,dg)[0,1]
    print('%d %5d %4d %4d | %.3f | %.4f / %.4f | %.4f | %.3f (r=%.2f) | %.3f %.3f'%(N,len(xs),len(by),len(multi),eta,
        np.median([np.ptp([o['ET1'] for o in v]) for v in multi]),max(np.ptp([o['ET1'] for o in v]) for v in multi),
        np.median([np.ptp([o['g0'] for o in v]) for v in multi]),slope,r,
        spearmanr([o['rho'] for o in xs],et).correlation,spearmanr([o['g0'] for o in xs],et).correlation))
# overlaps excluded at N=9
xs=[o for o in out if o['N']==9 and not o['overlaps']]; et=np.array([o['ET1'] for o in xs])
print('N=9 planar only (%d mols): Spearman rho %.3f g0 %.3f'%(len(xs),spearmanr([o['rho'] for o in xs],et).correlation,spearmanr([o['g0'] for o in xs],et).correlation))
