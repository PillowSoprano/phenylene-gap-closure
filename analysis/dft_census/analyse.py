"""B3LYP//xTB triplet energies of all classes with six to eight rings: Table S5, the B3LYP//xTB row of
Table 1, and the numbers of Section 6 and SI Section S7 (bootstrap intervals, pairs, <S^2>).
Input: results_N6-8.json, written by phenylene_dft.py (run in phenylene_dft_colab.ipynb).
Run ../competitors.py first (it writes ../competitors.json)."""
import itertools, json, os, sys
import numpy as np
from scipy.stats import spearmanr

here = os.path.dirname(os.path.abspath(__file__)); AN = os.path.join(here, '..'); BH = os.path.join(here, '..', '..', 'beyond_huckel')
sys.path.insert(0, AN)
from longest_linear_run import longest_L

R = json.load(open(os.path.join(here, 'results_N6-8.json')))
cal, M = R['calibration'], R['molecules']
D = {d['id']: d for d in json.load(open(os.path.join(BH, 'patterns_N2-9.json')))}
C = {r['id']: r for r in json.load(open(os.path.join(AN, 'competitors.json')))}
P = {json.loads(l)['id']: json.loads(l) for l in open(os.path.join(AN, 'ppp_bs_classes.jsonl'))}
PG = {r['id']: r['gapS0_ppp'] for r in json.load(open(os.path.join(BH, 'ppp_rohf_ohno_N2-9.json')))}
sp = lambda x, y: spearmanr(x, y)[0]

print('calibration (GPU against stored CPU energies):',
      f"max |dE| = {max(max(abs(c['dE_S0_vs_stored']), abs(c['dE_T1_vs_stored'])) for c in cal.values()):.1e} Eh")
ids = sorted(M)
print(f"{len(ids)} molecules, all converged: {all(M[i]['b3lyp_S0']['converged'] and M[i]['b3lyp_T1']['converged'] for i in ids)}, "
      f"<S2> {min(M[i]['b3lyp_T1']['S2'] for i in ids):.3f}-{max(M[i]['b3lyp_T1']['S2'] for i in ids):.3f}, "
      f"triplet from displaced start {sum(M[i]['xtb_T1_used'] == 'xtb_T1_pert' for i in ids)}, "
      f"xTB minus B3LYP//xTB {np.mean([M[i]['xtb_ET1_eV'] - M[i]['ET1_eV'] for i in ids]):.2f} eV, "
      f"E_T1 {min(M[i]['ET1_eV'] for i in ids):.2f}-{max(M[i]['ET1_eV'] for i in ids):.2f} eV")

pred = {'rho': lambda i: D[i]['rho'], 'Hueckel gap': lambda i: C[i]['g0'], 'elastic': lambda i: C[i]['lin_K'],
        'full Hessian': lambda i: C[i]['lin_F'], 'g0^2/P': lambda i: C[i]['g2P'], 'longest L run': lambda i: -longest_L(D[i]['la_tree']),
        'PPP gap': lambda i: PG[i], 'PPP E_T1': lambda i: P[i]['ET1_bs']}
rng = np.random.default_rng(0)
print('\nSpearman correlation with the B3LYP//xTB triplet energy (Table S5)')
print('N  n    ' + '  '.join(f'{k:>13s}' for k in pred))
for N in (6, 7, 8):
    s = [i for i in ids if M[i]['N'] == N]; E = np.array([M[i]['ET1_eV'] for i in s])
    X = {k: np.array([f(i) for i in s], float) for k, f in pred.items()}
    print(f"{N}  {len(s):3d}  " + '  '.join(f'{sp(X[k], E):13.3f}' for k in pred))
    for other in ('Hueckel gap', 'longest L run'):
        bs = []
        for _ in range(2000):
            b = rng.choice(len(s), len(s)); bs.append(sp(X['rho'][b], E[b]) - sp(X[other][b], E[b]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        print(f"      rho minus {other}: {sp(X['rho'], E) - sp(X[other], E):+.3f}, bootstrap 95% [{lo:+.3f}, {hi:+.3f}]")

E = {i: M[i]['ET1_eV'] for i in ids}; Lm = {i: longest_L(D[i]['la_tree']) for i in ids}
n = w = 0
for a, b in itertools.combinations(ids, 2):
    if M[a]['N'] != M[b]['N'] or abs(E[a] - E[b]) <= 0.05: continue
    sr, sg = np.sign(D[a]['rho'] - D[b]['rho']), np.sign(C[a]['g0'] - C[b]['g0'])
    if sr == 0 or sg == 0 or sr == sg: continue
    n += 1; w += sr == np.sign(E[a] - E[b])
print(f"\npairs with the same N that rho and the gap order oppositely, |dE| > 0.05 eV: {n}; B3LYP//xTB follows rho in {100 * w / n:.0f}%")
n = r = g = 0
for a, b in itertools.combinations(ids, 2):
    if M[a]['N'] != M[b]['N'] or Lm[a] != Lm[b] or abs(E[a] - E[b]) <= 0.05: continue
    s = np.sign(E[a] - E[b]); n += 1
    r += np.sign(D[a]['rho'] - D[b]['rho']) == s; g += np.sign(C[a]['g0'] - C[b]['g0']) == s
print(f"pairs with the same N and longest L run, |dE| > 0.05 eV: {n}; rho correct {100 * r / n:.0f}%, gap correct {100 * g / n:.0f}%")
