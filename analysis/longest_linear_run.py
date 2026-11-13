"""Longest run of linear rings as a predictor of the triplet energy (Table 1, Section 6).

For every class the predictor is the number of rings in the longest chain of adjacent rings marked L in the
L/A tree. Prints its Spearman correlation with the Hueckel-SSH and the lowest PPP-ROHF triplet energies for
eight and nine rings, and, for pairs of nine-ring classes with the same longest run whose PPP triplet energies
differ by more than 0.05 eV, how often rho and the Hueckel gap at the relaxed ground state order them correctly.
Run competitors.py first (it writes competitors.json)."""
import itertools, json, os
import numpy as np
from scipy.stats import spearmanr

here = os.path.dirname(os.path.abspath(__file__)); BH = os.path.join(here, '..', 'beyond_huckel')


def la_tree_graph(s):
    """Labels and adjacency of the L/A tree written as nested parentheses, e.g. '(A()(L()))'."""
    lab, adj, stack = [], [], []
    for c in s:
        if c == '(':
            k = len(lab); lab.append(''); adj.append([])
            if stack: adj[k].append(stack[-1]); adj[stack[-1]].append(k)
            stack.append(k)
        elif c == ')':
            stack.pop()
        else:
            lab[stack[-1]] += c
    return lab, adj


def longest_L(s):
    lab, adj = la_tree_graph(s)
    def far(x):                     # farthest L ring from x within its run of L rings
        d = {x: 0}; q = [x]
        while q:
            u = q.pop()
            for w in adj[u]:
                if lab[w] == 'L' and w not in d: d[w] = d[u] + 1; q.append(w)
        y = max(d, key=d.get); return y, d[y]
    return max([far(far(v)[0])[1] + 1 for v, l in enumerate(lab) if l == 'L'], default=0)


if __name__ == '__main__':
    D = {d['id']: d for d in json.load(open(os.path.join(BH, 'patterns_N2-9.json')))}
    P = {json.loads(l)['id']: json.loads(l) for l in open(os.path.join(here, 'ppp_bs_classes.jsonl'))}
    H = {r['id']: r for r in json.load(open(os.path.join(BH, 't1_scan.json')))}
    C = {r['id']: r for r in json.load(open(os.path.join(here, 'competitors.json')))}
    for N in (8, 9):
        ids = [i for i, d in D.items() if d['N'] == N and i in P and i != 'PT09-00092']
        L = [-longest_L(D[i]['la_tree']) for i in ids]
        print(f"N={N} ({len(ids)} classes)  longest L run: Hueckel-SSH {spearmanr(L, [H[i]['ET1_adiabatic'] for i in ids])[0]:.3f}, "
              f"PPP {spearmanr(L, [P[i]['ET1_bs'] for i in ids])[0]:.3f}")
    ids = [i for i, d in D.items() if d['N'] == 9 and i in P and i != 'PT09-00092']
    E = {i: P[i]['ET1_bs'] for i in ids}; Lm = {i: longest_L(D[i]['la_tree']) for i in ids}
    n = r = g = 0
    for a, b in itertools.combinations(ids, 2):
        if Lm[a] != Lm[b] or abs(E[a] - E[b]) <= 0.05: continue
        s = np.sign(E[a] - E[b]); n += 1
        r += np.sign(D[a]['rho'] - D[b]['rho']) == s; g += np.sign(C[a]['g0'] - C[b]['g0']) == s
    print(f"N=9, same longest L run, PPP energies differ by > 0.05 eV: {n} pairs, rho correct {100 * r / n:.0f}%, gap correct {100 * g / n:.0f}%")
