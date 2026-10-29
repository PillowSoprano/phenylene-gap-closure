import sys, json, time, os, numpy as np
from pyscf import gto, dft, scf, mcscf, lib
from chem import *
from gapclosure import bipartite
lib.num_threads(int(os.environ.get("NT","2")))
HART = 27.211386245988
name, ch = sys.argv[1], json.loads(sys.argv[2])
g, sym, _ = build(ch); nC = g['n']
X0 = np.load(f'{name}_eq.npy')
D = json.load(open('directions.json'))[name]
E_LEVELS = [0.5, 1.0, 1.5, 2.0, 3.0]
H_of = {}                                           # each H follows its carbon
for k in range(nC, len(sym)):
    H_of[k] = int(np.argmin(np.linalg.norm(X0[:nC]-X0[k], axis=1)))

def bond_lengths(X): return np.array([np.linalg.norm(X[u]-X[v]) for u, v in g['edges']])
r_eq = bond_lengths(X0)

def displace(target_dr, iters=8):
    X = X0.copy(); target = r_eq + target_dr
    for _ in range(iters):
        C = X[:nC, :2]; r = bond_lengths(X)
        B = np.zeros((len(g['edges']), 2*nC))
        for k, (u, v) in enumerate(g['edges']):
            e = (C[u]-C[v])/r[k]; B[k, 2*u:2*u+2] = e; B[k, 2*v:2*v+2] = -e
        dx = np.linalg.lstsq(B, target-r, rcond=None)[0].reshape(nC, 2)
        X[:nC, :2] += dx
        for hk, ck in H_of.items(): X[hk, :2] += dx[ck]
    ach = bond_lengths(X) - r_eq
    cos = float(ach@target_dr/np.linalg.norm(ach)/np.linalg.norm(target_dr))
    return X, float(np.abs(ach-target_dr).max()), cos

points = [('ref', None, 0.0, np.zeros(len(r_eq)))]
# random controls first: the efficiency comparison depends only on them (and P@1.5)
for i, r in enumerate(D['R']):
    for sgn in (+1, -1): points.append((f'R{i}{"+" if sgn>0 else "-"}@1.5', 'R', 1.5, sgn*np.sqrt(1.5)*np.array(r)))
for lab, key in [('P', 'P'), ('BLA', 'BLA')]:
    for Em in E_LEVELS: points.append((f'{lab}@{Em}', lab, Em, np.sqrt(Em)*np.array(D[key])))

outf = f'{name}_points.json'
res = json.load(open(outf)) if os.path.exists(outf) else {}
dm_ref = None
def follow(mf, max_iter=4):
    """Follow internal instabilities until the SCF solution is stable (or max_iter)."""
    stable = False
    for _ in range(max_iter):
        mo1, _, stable, _ = mf.stability(return_status=True)
        if stable: break
        mf.kernel(dm0=mf.make_rdm1(mo1, mf.mo_occ))
    return bool(stable)
def mkmol(X): return gto.M(atom=atom_str(sym, X), basis='def2-svp', verbose=0)
for lab, kind, Em, dr in points:
    if lab in res and (res[lab].get('stab') or kind not in (None, 'P')):
        if lab == 'ref' and dm_ref is None:
            m0 = dft.RKS(mkmol(X0)).density_fit(); m0.xc = 'b3lyp'; m0.kernel(); dm_ref = m0.make_rdm1()
        continue
    t = time.time(); X, err, cos = displace(dr); mol = mkmol(X)
    mf = dft.RKS(mol).density_fit(); mf.xc = 'b3lyp'; mf.conv_tol = 1e-9
    e = mf.kernel(dm0=dm_ref) if dm_ref is not None else mf.kernel()
    if lab == 'ref': dm_ref = mf.make_rdm1()
    no = mol.nelectron//2; gap = (mf.mo_energy[no]-mf.mo_energy[no-1])*HART
    rec = dict(kind=kind, Em=Em, E=e, gap=gap, maxerr=err, cos=cos, maxdr=float(np.abs(dr).max()))
    if kind in (None, 'P'):
        # broken-symmetry singlet from rotated HOMO/LUMO
        C = mf.mo_coeff.copy(); h, l = C[:, no-1].copy(), C[:, no].copy()
        Ca, Cb = C.copy(), C.copy()
        Ca[:, no-1] = (h+l)/np.sqrt(2); Ca[:, no] = (h-l)/np.sqrt(2)
        Cb[:, no-1] = (h-l)/np.sqrt(2); Cb[:, no] = (h+l)/np.sqrt(2)
        occ = np.zeros_like(mf.mo_occ); occ[:no] = 1
        molT = mol.copy(); molT.spin = 2; molT.build()
        tm = dft.UKS(molT).density_fit(); tm.xc = 'b3lyp'; tm.conv_tol = 1e-9
        et = tm.kernel()
        # diradical character y0 from broken-symmetry UHF natural orbitals (Yamaguchi / Nakano PUHF formula)
        hf = scf.RHF(mol).density_fit(); hf.conv_tol = 1e-9; hf.kernel()
        Ch = hf.mo_coeff.copy(); hh, ll = Ch[:, no-1].copy(), Ch[:, no].copy(); Cha, Chb = Ch.copy(), Ch.copy()
        Cha[:, no-1] = (hh+ll)/np.sqrt(2); Cha[:, no] = (hh-ll)/np.sqrt(2)
        Chb[:, no-1] = (hh-ll)/np.sqrt(2); Chb[:, no] = (hh+ll)/np.sqrt(2)
        # two broken-symmetry starting points: HOMO/LUMO mixing, and an antiferromagnetic guess that puts
        # opposite spin density on the two carbon sublattices; keep the lower UHF solution
        D = hf.make_rdm1()
        sl = mol.aoslice_by_atom()
        mA = np.zeros(mol.nao); mB = np.zeros(mol.nao)
        Lset = set(bipartite(g)[1])
        for at in range(nC):
            (mA if at in Lset else mB)[sl[at, 2]:sl[at, 3]] = 1
        delta = 0.5*(mA[:, None]*D*mA[None, :] - mB[:, None]*D*mB[None, :])
        guesses = [uh_tmp.make_rdm1((Cha, Chb), (occ, occ)) for uh_tmp in [scf.UHF(mol)]] + [(D/2 + delta/2, D/2 - delta/2)]
        uhs = []
        for dm0 in guesses:
            u_ = scf.UHF(mol).density_fit(); u_.conv_tol = 1e-9; u_.max_cycle = 200; u_.kernel(dm0=dm0); uhs.append(u_)
        uh = min(uhs, key=lambda u_: u_.e_tot); euh = uh.e_tot; uhf_stable = None
        uhf_candidates = [float(u_.e_tot) for u_ in uhs]
        Da, Db = uh.make_rdm1(); Sm = mol.intor('int1e_ovlp')
        w, V = np.linalg.eigh(Sm); Sh = V @ np.diag(np.sqrt(w)) @ V.T
        nocc_nat = np.sort(np.linalg.eigvalsh(Sh @ (Da+Db) @ Sh))[::-1]
        nH, nL = nocc_nat[no-1], nocc_nat[no]; Tq = (nH-nL)/2; y0 = 1 - 2*Tq/(1+Tq**2)
        # broken-symmetry UKS: from the HOMO/LUMO-mixed guess and from the stable UHF density; keep the lower
        cands = []
        for dm0 in (um_guess := dft.UKS(mol).make_rdm1((Ca, Cb), (occ, occ)), uh.make_rdm1()):
            um = dft.UKS(mol).density_fit(); um.xc = 'b3lyp'; um.conv_tol = 1e-9; um.kernel(dm0=dm0)
            cands.append((um.e_tot, um.spin_square()[0], bool(um.converged)))
        ebs, s2, bs_conv = min(cands)
        rec.update(E_BS=ebs, S2_BS=float(s2), E_T=et, y0=float(y0), n_LUNO=float(nL), S2_UHF=float(uh.spin_square()[0]),
                   E_UHF=float(euh), E_RHF=float(hf.e_tot), bs_converged=bs_conv, uhf_converged=bool(uh.converged), uhf_stable=uhf_stable, E_UHF_candidates=uhf_candidates, E_BS_candidates=[float(x[0]) for x in cands], stab=True)
    res[lab] = rec
    json.dump(res, open(outf, 'w'), indent=1)
    print(f'{lab:12s} Em={Em:4.1f} dE={(e-res["ref"]["E"])*HART:7.3f} eV gap={gap:6.3f} cos={cos:.4f} maxerr={err:.1e}'
          + (f' y0={rec["y0"]:.4f} S2={rec["S2_BS"]:.3f} dBS={(ebs-e)*HART:+.3f}' if 'y0' in rec else '') + f'  ({time.time()-t:.0f}s)', flush=True)
print('done', name)
