"""Triplet energies of the two hook molecules (Fig. 1) with GFN2-xTB geometries and UB3LYP/6-31G* single points.

Protocol, for each molecule:
  1. S0: GFN2-xTB optimization from the planar skeleton with small random out-of-plane noise.
  2. T1: GFN2-xTB optimization (--uhf 2) from the S0 geometry, once as is and once with random
     displacements of up to 0.02 A; the lower of the two is kept.
  3. UB3LYP/6-31G* (density fitting) single points: S0 at the S0 geometry, T1 at the T1 geometry.
     ET1 = E(T1 @ T1 geometry) - E(S0 @ S0 geometry)  (adiabatic, xTB geometries).
Calibration molecules (biphenylene, angular [3], helical [4]phenylene) have full UB3LYP/6-31G*
optimized ET1 from dft_isomers.json, so the error of the protocol can be read off directly.

Needs: xtb on PATH (conda install -c conda-forge xtb), pyscf, numpy, the gapclosure package.
Run the nine-ring molecules one after the other (each B3LYP single point needs ~11 GB of temporary disk):
         caffeinate -i sh -c 'python dft_hook_ab.py A > dft_hook_A.log 2>&1; python dft_hook_ab.py B > dft_hook_B.log 2>&1' &
Restart-safe: finished steps are kept in dft_hook_ab.json.
"""
import json, os, subprocess, sys, time
import numpy as np
from pyscf import gto, dft, lib

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'gapclosure', 'src'))
from gapclosure.skeleton import phenylene_tree

lib.num_threads(os.cpu_count())
_TAG = ('_' + '_'.join(sys.argv[1:])) if sys.argv[1:] else ''   # one file per process when run in parallel
OUT = os.path.join(here, 'dft_hook_ab' + _TAG + '.json')
WORK = os.path.join(here, 'dft_hook_ab_work')
HARTREE_EV = 27.211386

MOLECULES = {   # name: (class id, links); calibration molecules first
    'biphenylene': ('PT02-00001', [[0, 0, 1]]),
    'angular3':    ('PT03-00001', [[0, 0, 1], [0, 2, 2]]),
    'helical4':    ('PT04-00002', [[0, 0, 1], [0, 2, 2], [1, 2, 3]]),
    'A':           ('PT09-00009', [[0, 0, 1], [0, 2, 2], [0, 4, 3], [1, 2, 4], [1, 4, 5], [2, 2, 6], [4, 2, 7], [5, 2, 8]]),
    'B':           ('PT09-00036', [[0, 0, 1], [0, 2, 2], [1, 2, 3], [1, 4, 4], [2, 2, 5], [2, 4, 6], [3, 3, 7], [4, 3, 8]]),
}


def build(links, seed=0):
    sk = phenylene_tree(links)
    X = np.asarray(sk.pos, float)
    X = X / np.median([np.linalg.norm(X[a] - X[b]) for a, b in sk.bonds]) * 1.40
    nb = {i: [] for i in range(sk.n)}
    for a, b in sk.bonds:
        nb[a].append(b); nb[b].append(a)
    rng = np.random.default_rng(seed)
    atoms = [('C', (x, y, rng.normal(0, 0.05))) for x, y in X]
    for i in range(sk.n):
        if len(nb[i]) == 2:
            v = 2 * X[i] - X[nb[i][0]] - X[nb[i][1]]; v /= np.linalg.norm(v)
            h = X[i] + 1.08 * v
            atoms.append(('H', (h[0], h[1], rng.normal(0, 0.05))))
    return atoms


def write_xyz(path, atoms, comment=''):
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for s, (x, y, z) in atoms:
            f.write(f"{s} {x:.8f} {y:.8f} {z:.8f}\n")


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    atoms = [(l.split()[0], tuple(float(v) for v in l.split()[1:4])) for l in lines[2:2 + n]]
    energy = None
    for tok in lines[1].split():
        try:
            energy = float(tok); break
        except ValueError:
            pass
    return atoms, energy


def xtb_opt(atoms, tag, uhf):
    d = os.path.join(WORK, tag); os.makedirs(d, exist_ok=True)
    write_xyz(os.path.join(d, 'start.xyz'), atoms)
    cmd = ['xtb', 'start.xyz', '--gfn', '2', '--opt', 'tight', '--uhf', str(uhf)]
    with open(os.path.join(d, 'xtb.out'), 'w') as log:
        subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count())))
    geo, e = read_xyz(os.path.join(d, 'xtbopt.xyz'))
    return geo, e


def displace(atoms, seed, amp=0.02):
    rng = np.random.default_rng(seed)
    return [(s, tuple(np.array(c) + rng.uniform(-amp, amp, 3))) for s, c in atoms]


def b3lyp_sp(atoms, spin):
    # Density fitting keeps a three-index tensor in a temporary file (about 11 GB for the
    # nine-ring molecules). Run the big molecules one at a time and check the free space first.
    import gc, shutil, tempfile
    mol = gto.M(atom=atoms, basis='6-31g*', spin=spin, verbose=0)
    need = mol.nao * (mol.nao + 1) / 2 * 4300 * 8 * 1.15
    free = shutil.disk_usage(lib.param.TMPDIR or tempfile.gettempdir()).free
    if free < need:
        raise SystemExit(f"only {free / 1e9:.1f} GB free, density fitting needs about {need / 1e9:.1f} GB")
    mf = (dft.RKS(mol) if spin == 0 else dft.UKS(mol)).density_fit()
    mf.xc = 'b3lyp'; mf.conv_tol = 1e-9; mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:
        mf = mf.newton(); e = mf.kernel()
    s2 = float(mf.spin_square()[0]) if spin else 0.0
    ok = bool(mf.converged)
    del mf
    gc.collect()          # releases the temporary density-fitting file
    return float(e), s2, ok


def save(res):
    json.dump(res, open(OUT, 'w'), indent=1)


def main(names):
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for name in names:
        cid, links = MOLECULES[name]
        r = res.setdefault(name, dict(id=cid))
        t0 = time.time()
        if 'xtb_S0' not in r:
            geo, e = xtb_opt(build(links), f'{name}_S0', 0)
            r['xtb_S0'] = dict(E=e, geometry=geo); save(res)
            print(f"{name}: xTB S0 done, {time.time() - t0:.0f} s", flush=True)
        s0 = [(s, tuple(c)) for s, c in r['xtb_S0']['geometry']]
        for k, start in (('xtb_T1_sym', s0), ('xtb_T1_pert', displace(s0, 1))):
            if k not in r:
                geo, e = xtb_opt(start, f'{name}_{k}', 2)
                r[k] = dict(E=e, geometry=geo); save(res)
                print(f"{name}: {k} done, E = {e:.6f} Eh, {time.time() - t0:.0f} s", flush=True)
        t1key = min(('xtb_T1_sym', 'xtb_T1_pert'), key=lambda k: r[k]['E'])
        r['xtb_T1_used'] = t1key
        r['xtb_ET1_eV'] = (r[t1key]['E'] - r['xtb_S0']['E']) * HARTREE_EV
        if 'b3lyp_S0' not in r:
            e, s2, ok = b3lyp_sp(s0, 0)
            r['b3lyp_S0'] = dict(E=e, converged=ok); save(res)
            print(f"{name}: B3LYP S0 single point {e:.8f} Eh, {time.time() - t0:.0f} s", flush=True)
        if 'b3lyp_T1' not in r:
            e, s2, ok = b3lyp_sp([(s, tuple(c)) for s, c in r[t1key]['geometry']], 2)
            r['b3lyp_T1'] = dict(E=e, S2=s2, converged=ok); save(res)
        r['ET1_eV'] = (r['b3lyp_T1']['E'] - r['b3lyp_S0']['E']) * HARTREE_EV
        save(res)
        print(f"{name}: ET1 = {r['ET1_eV']:.3f} eV (B3LYP//xTB), xTB alone {r['xtb_ET1_eV']:.3f} eV, "
              f"<S2> = {r['b3lyp_T1']['S2']:.3f}, T1 from {t1key}, {time.time() - t0:.0f} s", flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or list(MOLECULES))
