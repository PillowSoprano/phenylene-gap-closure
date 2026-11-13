"""Hueckel-SSH lowest triplet of every class (N <= 9): symmetric start (as in t1_scan.py) and four perturbed starts.
Writes hssh_bs_classes.json with ET1_plain and ET1_lowest in |beta|."""
import json, os, sys, numpy as np
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__))
R0, X, Y, RREF = 1.523, 0.189, 0.2756, 1.397
K = 2 / (X * Y)

def run(d):
    E = [(a, b) for a, b, _ in d['cp']]; n = 6 * d['N']; nh = n // 2
    ia = np.array([a for a, _ in E]); ib = np.array([b for _, b in E])
    def solve(R, occ):
        H = np.zeros((n, n)); t = 1 + (RREF - R) / Y; H[ia, ib] = -t; H[ib, ia] = -t
        e, C = np.linalg.eigh(H); w = np.zeros(n); w[:nh] = 2
        if occ == 'T1': w[nh - 1] = 1; w[nh] = 1
        return (C[ia, :] * C[ib, :]) @ w, float(w @ e)
    def relax(occ, R):
        for _ in range(6000):
            p, _ = solve(R, occ); Rn = R0 - X * p
            if np.max(abs(Rn - R)) < 1e-11: break
            R = 0.6 * R + 0.4 * Rn
        p, eel = solve(R, occ)
        return eel + K / 2 * np.sum((R - R0) ** 2), R
    ES, RS = relax('S0', np.full(len(E), RREF))
    ET0, _ = relax('T1', np.full(len(E), RREF))
    rng = np.random.default_rng(1); ets = [ET0]
    for k in range(4):
        ets.append(relax('T1', RS + rng.normal(0, 0.02, len(E)))[0])
    return dict(id=d['id'], N=d['N'], rho=d['rho'], ET1_plain=ET0 - ES, ET1_lowest=min(ets) - ES)

if __name__ == '__main__':
    D = json.load(open(os.path.join(here, '..', 'beyond_huckel', 'patterns_N2-9.json')))
    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, len(D))
    out_path = os.path.join(here, 'hssh_bs_classes.json')
    res = json.load(open(out_path)) if os.path.exists(out_path) else []
    done = {r['id'] for r in res}
    with Pool(os.cpu_count()) as p:
        res += p.map(run, [d for d in D[lo:hi] if d['id'] not in done])
    json.dump(res, open(out_path, 'w'), indent=0)
    print(len(res))
