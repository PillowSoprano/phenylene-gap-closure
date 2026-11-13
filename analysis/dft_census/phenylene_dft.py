"""Triplet energies of tree-like phenylenes from UB3LYP/6-31G* single points at GFN2-xTB geometries.

Same protocol as analysis/dft_hook_ab.py of the repository, for whole sets of L/A classes:
  1. S0: GFN2-xTB optimization from the planar census geometry with small random out-of-plane noise.
  2. T1: GFN2-xTB optimization (--uhf 2) from the S0 geometry, once as is and once with random
     displacements of up to 0.02 A; the lower of the two is kept.
  3. B3LYP/6-31G* single points with density fitting: S0 (restricted) at the S0 geometry,
     T1 (unrestricted) at the T1 geometry. E_T1 = E(T1) - E(S0), adiabatic at xTB geometries.
The molecule of each class is the one used in the census (beyond_huckel/coords_N2-9.json).
Runs on a GPU with gpu4pyscf, or on the CPU with PySCF if no GPU is available.
Every finished step is written to the results file at once, so an interrupted run resumes where it stopped.
"""
import json, os, shutil, subprocess, time, urllib.request
import numpy as np

REPO = 'https://raw.githubusercontent.com/PillowSoprano/phenylene-gap-closure/main/'
HARTREE_EV = 27.211386
CALIBRATION = {   # stored CPU results at the same xTB geometries (analysis/dft_hook_ab_*.json)
    'biphenylene': 'analysis/dft_hook_ab_biphenylene_angular3_helical4.json',
    'angular3': 'analysis/dft_hook_ab_biphenylene_angular3_helical4.json',
    'helical4': 'analysis/dft_hook_ab_biphenylene_angular3_helical4.json',
    'A': 'analysis/dft_hook_ab_A.json',
    'B': 'analysis/dft_hook_ab_B.json',
}

try:
    import cupy  # noqa: F401
    from gpu4pyscf.dft import rks as _rks, uks as _uks
    GPU = True
except Exception:
    from pyscf import dft as _dft
    GPU = False
from pyscf import gto, lib


def fetch(path, cache='repo_data'):
    os.makedirs(cache, exist_ok=True)
    local = os.path.join(cache, path.replace('/', '__'))
    if not os.path.exists(local):
        urllib.request.urlretrieve(REPO + path, local)
    return json.load(open(local))


def load_classes(Ns):
    pats = [d for d in fetch('beyond_huckel/patterns_N2-9.json') if d['N'] in Ns]
    coords = fetch('beyond_huckel/coords_N2-9.json')
    out = []
    for d in sorted(pats, key=lambda d: (d['N'], d['id'])):
        if d['id'] == 'PT09-00092':          # excluded in the paper (Hartree-Fock instability)
            continue
        out.append(dict(id=d['id'], N=d['N'], la_tree=d['la_tree'], rho=d['rho'],
                        bonds=[(a, b) for a, b, _ in d['cp']], pos=coords[d['id']]))
    return out


def build(rec, seed=0):
    X = np.asarray(rec['pos'], float) * 1.40          # census coordinates have unit bond length
    nb = {i: [] for i in range(len(X))}
    for a, b in rec['bonds']:
        nb[a].append(b); nb[b].append(a)
    rng = np.random.default_rng(seed)
    atoms = [('C', (x, y, rng.normal(0, 0.05))) for x, y in X]
    for i in range(len(X)):
        if len(nb[i]) == 2:
            v = 2 * X[i] - X[nb[i][0]] - X[nb[i][1]]; v /= np.linalg.norm(v)
            h = X[i] + 1.08 * v
            atoms.append(('H', (h[0], h[1], rng.normal(0, 0.05))))
    return atoms


def write_xyz(path, atoms):
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n\n")
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


def xtb_opt(atoms, workdir, uhf, xtb='xtb'):
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    write_xyz(os.path.join(workdir, 'start.xyz'), atoms)
    with open(os.path.join(workdir, 'xtb.out'), 'w') as log:
        subprocess.run([xtb, 'start.xyz', '--gfn', '2', '--opt', 'tight', '--uhf', str(uhf)], cwd=workdir,
                       stdout=log, stderr=subprocess.STDOUT, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count()), OMP_STACKSIZE='1G'))
    return read_xyz(os.path.join(workdir, 'xtbopt.xyz'))


def displace(atoms, seed, amp=0.02):
    rng = np.random.default_rng(seed)
    return [(s, tuple(np.array(c) + rng.uniform(-amp, amp, 3))) for s, c in atoms]


def b3lyp_sp(atoms, spin):
    """B3LYP/6-31G* single point with density fitting; returns energy (Eh), <S^2>, converged."""
    mol = gto.M(atom=atoms, basis='6-31g*', spin=spin, verbose=0)
    if GPU:
        mf = (_rks.RKS(mol, xc='b3lyp') if spin == 0 else _uks.UKS(mol, xc='b3lyp')).density_fit()
    else:
        mf = (_dft.RKS(mol) if spin == 0 else _dft.UKS(mol)).density_fit(); mf.xc = 'b3lyp'
    mf.grids.level = 3; mf.conv_tol = 1e-9; mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:                       # rare: retry with a level shift, then second order on the CPU
        mf.level_shift = 0.2; mf.max_cycle = 300; e = mf.kernel(dm0=mf.make_rdm1())
        if not mf.converged:
            cpu = mf.to_cpu() if GPU else mf
            cpu.level_shift = 0; cpu = cpu.newton(); e = cpu.kernel(); mf = cpu
    s2 = 0.0
    if spin:
        m = mf.to_cpu() if GPU and hasattr(mf, 'to_cpu') else mf
        s2 = float(m.spin_square()[0])
    ok = bool(mf.converged)
    del mf
    if GPU:
        import cupy
        cupy.get_default_memory_pool().free_all_blocks()
    return float(e), s2, ok


def _save(res, path):
    tmp = path + '.tmp'
    json.dump(res, open(tmp, 'w'), indent=1)
    os.replace(tmp, path)


def calibrate(results_path, names=('biphenylene', 'angular3', 'helical4', 'A', 'B')):
    """Single points at the stored xTB geometries; compare with the stored CPU energies."""
    res = json.load(open(results_path)) if os.path.exists(results_path) else {}
    cal = res.setdefault('calibration', {})
    for name in names:
        if name in cal:
            continue
        ref = fetch(CALIBRATION[name])[name]
        t0 = time.time()
        s0 = [(s, tuple(c)) for s, c in ref['xtb_S0']['geometry']]
        t1 = [(s, tuple(c)) for s, c in ref[ref['xtb_T1_used']]['geometry']]
        e0, _, ok0 = b3lyp_sp(s0, 0)
        e1, s2, ok1 = b3lyp_sp(t1, 2)
        et1 = (e1 - e0) * HARTREE_EV
        cal[name] = dict(E_S0=e0, E_T1=e1, S2=s2, converged=ok0 and ok1, ET1_eV=et1,
                         dE_S0_vs_stored=e0 - ref['b3lyp_S0']['E'], dE_T1_vs_stored=e1 - ref['b3lyp_T1']['E'],
                         ET1_stored_eV=ref['ET1_eV'], seconds=time.time() - t0)
        _save(res, results_path)
        print(f"{name:12s} E_T1 {et1:.4f} eV (stored {ref['ET1_eV']:.4f}),  dE(S0) {e0 - ref['b3lyp_S0']['E']:+.1e} Eh, "
              f"dE(T1) {e1 - ref['b3lyp_T1']['E']:+.1e} Eh,  {time.time() - t0:.0f} s", flush=True)
    return cal


def run(classes, results_path, workroot='xtb_work', xtb='xtb'):
    res = json.load(open(results_path)) if os.path.exists(results_path) else {}
    mols = res.setdefault('molecules', {})
    for k, rec in enumerate(classes):
        cid = rec['id']
        r = mols.setdefault(cid, dict(id=cid, N=rec['N'], la_tree=rec['la_tree'], rho=rec['rho']))
        if 'ET1_eV' in r:
            continue
        t0 = time.time()
        if 'xtb_S0' not in r:
            geo, e = xtb_opt(build(rec), os.path.join(workroot, cid, 'S0'), 0, xtb)
            r['xtb_S0'] = dict(E=e, geometry=geo); _save(res, results_path)
        s0 = [(s, tuple(c)) for s, c in r['xtb_S0']['geometry']]
        for key, start in (('xtb_T1_sym', s0), ('xtb_T1_pert', displace(s0, 1))):
            if key not in r:
                geo, e = xtb_opt(start, os.path.join(workroot, cid, key), 2, xtb)
                r[key] = dict(E=e, geometry=geo); _save(res, results_path)
        t1key = min(('xtb_T1_sym', 'xtb_T1_pert'), key=lambda q: r[q]['E'])
        r['xtb_T1_used'] = t1key
        r['xtb_ET1_eV'] = (r[t1key]['E'] - r['xtb_S0']['E']) * HARTREE_EV
        t1 = time.time()
        if 'b3lyp_S0' not in r:
            e, _, ok = b3lyp_sp(s0, 0)
            r['b3lyp_S0'] = dict(E=e, converged=ok); _save(res, results_path)
        if 'b3lyp_T1' not in r:
            e, s2, ok = b3lyp_sp([(s, tuple(c)) for s, c in r[t1key]['geometry']], 2)
            r['b3lyp_T1'] = dict(E=e, S2=s2, converged=ok); _save(res, results_path)
        r['ET1_eV'] = (r['b3lyp_T1']['E'] - r['b3lyp_S0']['E']) * HARTREE_EV
        r['seconds'] = dict(xtb=t1 - t0, dft=time.time() - t1)
        _save(res, results_path)
        print(f"[{k + 1}/{len(classes)}] {cid} N={rec['N']}  E_T1 {r['ET1_eV']:.3f} eV  (xTB {r['xtb_ET1_eV']:.3f})  "
              f"<S2> {r['b3lyp_T1']['S2']:.3f}  xTB {t1 - t0:.0f} s, DFT {time.time() - t1:.0f} s", flush=True)
    return res
