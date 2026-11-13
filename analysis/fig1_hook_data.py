"""Fig. 1: lowest PPP-ROHF triplet energies of all classes with nine rings against rho and against the Hueckel gap
at equal couplings, with the two marked molecules A (PT09-00009) and B (PT09-00036).
Prints the data as pgfplots coordinate lists."""
import json, os, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
D = {d['id']: d for d in json.load(open(os.path.join(here, '..', 'beyond_huckel', 'patterns_N2-9.json')))}
P = {json.loads(l)['id']: json.loads(l) for l in open(os.path.join(here, 'ppp_bs_classes.jsonl'))}
EXCLUDED = {'PT09-00092'}          # triplet below the restricted ground state (Hartree-Fock instability)

def gap(d):
    n = 6 * d['N']; H = np.zeros((n, n))
    for a, b, _ in d['cp']: H[a, b] = H[b, a] = -1.0
    e = np.linalg.eigvalsh(H); return e[n // 2] - e[n // 2 - 1]

rows = [(i, D[i]['rho'], gap(D[i]), P[i]['ET1_bs']) for i in P if P[i]['N'] == 9 and i not in EXCLUDED]
print(f'{len(rows)} classes with nine rings')
for tag, i in (('A', 'PT09-00009'), ('B', 'PT09-00036')):
    r = next(x for x in rows if x[0] == i)
    print(f'{tag} = {i}: rho {r[1]:.3f}, gap {r[2]:.4f}, E_T1 {r[3]:.3f} eV')
print('E_T1 vs rho:', ' '.join(f'({r[1]:.4f},{r[3]:.3f})' for r in rows))
print('E_T1 vs gap:', ' '.join(f'({r[2]:.4f},{r[3]:.3f})' for r in rows))
