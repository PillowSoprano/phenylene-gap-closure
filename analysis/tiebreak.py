"""Why the gap wins at N=4,5 and rho wins from N=6 on: split the ranking into a between-class and a within-class part."""
import json, numpy as np, collections
exec(open('competitors.py').read().split("def analyse")[0])
out=[json.loads(l) for l in open('within_class_mols.jsonl')]
print('N mols cls | mols/cls | share of T1 variance within classes | Spearman over molecules: rho, gap, rho-then-gap | class means: rho, mean gap | pairs of molecules: same class / ordered wrongly by gap')
for N in range(4,10):
    xs=[o for o in out if o['N']==N]; by=collections.defaultdict(list)
    for o in xs: by[o['id']].append(o)
    et=np.array([o['ET1'] for o in xs]); rho=np.array([o['rho'] for o in xs]); g=np.array([o['g0'] for o in xs])
    within=sum(np.sum((np.array([o['ET1'] for o in v])-np.mean([o['ET1'] for o in v]))**2) for v in by.values())/np.sum((et-et.mean())**2)
    lex=rho+1e-6*g/g.max()            # rho first, gap breaks ties
    cm=[(v[0]['rho'],np.mean([o['g0'] for o in v]),np.mean([o['ET1'] for o in v])) for v in by.values()]
    cr=spearmanr([c[0] for c in cm],[c[2] for c in cm]).correlation; cg=spearmanr([c[1] for c in cm],[c[2] for c in cm]).correlation
    # discordant pairs for gap, split by same/different class
    n=len(xs); same=diff=0; tot_same=tot_diff=0
    for i in range(n):
        for j in range(i+1,n):
            s=xs[i]['id']==xs[j]['id']
            disc=(g[i]-g[j])*(et[i]-et[j])<0
            if s: tot_same+=1; same+=disc
            else: tot_diff+=1; diff+=disc
    print('%d %4d %3d | %.2f | %.4f | %.3f %.3f %.3f | %.3f %.3f | gap discordant: same-class %d/%d, cross-class %d/%d (%.1f%%)'%(
        N,n,len(by),n/len(by),within,spearmanr(rho,et).correlation,spearmanr(g,et).correlation,spearmanr(lex,et).correlation,cr,cg,same,tot_same,diff,tot_diff,100*diff/max(tot_diff,1)))
