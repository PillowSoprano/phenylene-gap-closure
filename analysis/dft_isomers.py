"""UB3LYP/6-31G* S0 and T1 of phenylenes (Sections 3, 7 and 8, Table S9), at the same level as beyond_huckel/dft_t1.py.

Jobs, in order of priority:
  1. zigzag [4]phenylene (links [[0,0,1],[0,2,2],[1,4,3]]): S0, T1 from the S0 geometry, T1 from a perturbed geometry.
     beyond_huckel/dft_t1.py covers the other member of class AA, helical [4]phenylene (gap 0.703).
  2. helical [4]phenylene: S0 (to save the geometry), T1 from a perturbed geometry.
  3. the other five synthesized phenylenes: S0 (geometry saved), T1 from a perturbed geometry.
The perturbed start displaces every atom by up to 0.02 A (also out of plane) so that the optimizer can leave the
symmetric structure. For each molecule the lowest T1 is the one to use.

Usage (from analysis/):
    python dft_isomers.py                 # all jobs, results appended to dft_isomers.json
    python dft_isomers.py zigzag4 helical4  # selected molecules
Already finished steps are skipped when the script is restarted.
Set DFT_TEST=1 for a quick check of the code path (STO-3G, few optimizer steps).
"""
import os, sys, json, time
import numpy as np
from pyscf import gto, dft, lib
from pyscf.geomopt.geometric_solver import optimize

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'gapclosure', 'src'))
from gapclosure.skeleton import phenylene_tree

TEST = os.environ.get('DFT_TEST') == '1'
BASIS = 'sto-3g' if TEST else '6-31g*'
MAXSTEPS = 4 if TEST else 100
lib.num_threads(os.cpu_count())
_ONLY = os.environ.get('DFT_STEPS')          # e.g. DFT_STEPS=T1_pert: run only these steps, S0 read from the base file
_TAG = ('_' + '_'.join(sys.argv[1:])) if sys.argv[1:] else ''   # one file per process when run in parallel
_BASE = os.path.join(here, ('dft_isomers_test' if TEST else 'dft_isomers') + _TAG + '.json')
OUT = _BASE[:-5] + '_' + _ONLY.replace(',', '_') + '.json' if _ONLY else _BASE

MOLECULES = {   # name: (class id, links)
    'zigzag4':     ('PT04-00002', [[0, 0, 1], [0, 2, 2], [1, 4, 3]]),
    'helical4':    ('PT04-00002', [[0, 0, 1], [0, 2, 2], [1, 2, 3]]),
    'biphenylene': ('PT02-00001', [[0, 0, 1]]),
    'linear3':     ('PT03-00002', [[0, 0, 1], [0, 3, 2]]),
    'angular3':    ('PT03-00001', [[0, 0, 1], [0, 2, 2]]),
    'linear4':     ('PT04-00004', [[0, 0, 1], [0, 3, 2], [1, 3, 3]]),
    'triangular4': ('PT04-00001', [[0, 0, 1], [0, 2, 2], [0, 4, 3]]),
}
STEPS = {'zigzag4': ['S0', 'T1_sym', 'T1_pert'], 'helical4': ['S0', 'T1_pert']}
DEFAULT_STEPS = ['S0', 'T1_pert']


def build(links):
    sk = phenylene_tree(links)
    X = np.asarray(sk.pos, float)
    X = X / np.median([np.linalg.norm(X[a] - X[b]) for a, b in sk.bonds]) * 1.40
    nb = {i: [] for i in range(sk.n)}
    for a, b in sk.bonds:
        nb[a].append(b); nb[b].append(a)
    atoms = [('C', (x, y, 0.0)) for x, y in X]
    for i in range(sk.n):
        if len(nb[i]) == 2:
            v = 2 * X[i] - X[nb[i][0]] - X[nb[i][1]]; v /= np.linalg.norm(v)
            h = X[i] + 1.08 * v
            atoms.append(('H', (h[0], h[1], 0.0)))
    return sk, atoms


def perturb(atoms, seed, amp=0.02):
    rng = np.random.default_rng(seed)
    return [(s, tuple(np.array(c) + rng.uniform(-amp, amp, 3))) for s, c in atoms]


def run(atoms, spin):
    mol = gto.M(atom=atoms, basis=BASIS, spin=spin, verbose=0)
    mf = dft.UKS(mol).density_fit(); mf.xc = 'b3lyp'; mf.conv_tol = 1e-9
    mol_eq = optimize(mf, maxsteps=MAXSTEPS)
    mf2 = dft.UKS(mol_eq).density_fit(); mf2.xc = 'b3lyp'; mf2.conv_tol = 1e-10
    e = mf2.kernel()
    s2 = mf2.spin_square()[0]
    geo = [(mol_eq.atom_symbol(i), [float(v) for v in mol_eq.atom_coord(i, unit='Angstrom')]) for i in range(mol_eq.natm)]
    return float(e), float(s2), geo


def cc_lengths(sk, geo):
    X = np.array([c for _, c in geo[:sk.n]])
    return [[int(a), int(b), float(np.linalg.norm(X[a] - X[b]))] for a, b in sk.bonds]


def main(names):
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for name in names:
        cid, links = MOLECULES[name]
        sk, atoms = build(links)
        r = res.setdefault(name, dict(id=cid, links=links))
        if _ONLY and 'S0' not in r:
            r['S0'] = json.load(open(_BASE))[name]['S0']
        for step in (_ONLY.split(',') if _ONLY else STEPS.get(name, DEFAULT_STEPS)):
            if step in r:
                continue
            t = time.time()
            if step == 'S0':
                e, s2, geo = run(atoms, 0)
            else:
                if 'S0' not in r:
                    raise RuntimeError(f'{name}: S0 must come first')
                start = [(s, tuple(c)) for s, c in r['S0']['geometry']]
                if step == 'T1_pert':
                    start = perturb(start, seed=1)
                e, s2, geo = run(start, 2)
            r[step] = dict(E_hartree=e, S2=s2, geometry=geo, CC_bonds=cc_lengths(sk, geo), seconds=time.time() - t)
            if step != 'S0':
                r[step]['ET1_eV'] = (e - r['S0']['E_hartree']) * 27.211386
            json.dump(res, open(OUT, 'w'), indent=1)
            msg = f"{name} {step}: E = {e:.8f} Eh, <S2> = {s2:.3f}, {time.time() - t:.0f} s"
            if 'ET1_eV' in r[step]:
                msg += f", ET1 = {r[step]['ET1_eV']:.3f} eV"
            print(msg, flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or list(MOLECULES))
