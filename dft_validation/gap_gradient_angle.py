"""Angle between the predicted pattern x* (Hueckel bond-length model) and the B3LYP gap gradient."""
import json, numpy as np, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from gapclosure import phenylene, BondModel
HART = 27.211386245988
out = {}
for name, ch in [('biphenylene', []), ('angular3', ['L'])]:
    try: G = json.load(open(f'{name}_gapgrad.json'))
    except FileNotFoundError: continue
    if 'grad_eV_per_A' not in G: continue
    gB = np.array(G['grad_eV_per_A'])                       # d(gap)/d(bond length), eV/A, B3LYP
    g = phenylene(ch); bm = BondModel(g, K=43.0); d, _ = bm.relax()
    U, s, Vt = np.linalg.svd(bm.Tw(d)); gM = -2*bm.a*U[:, -1][bm.ii]*Vt[-1][bm.jj]   # model gap gradient
    hh = 1e-5; I = np.eye(bm.m)
    F = np.column_stack([(bm.energy(d+hh*I[k])[1]-bm.energy(d-hh*I[k])[1])/(2*hh) for k in range(bm.m)]); F = (F+F.T)/2
    Fi = np.linalg.inv(F)
    xs = -Fi @ gM                                            # predicted pattern (direction)
    xB = -Fi @ gB                                            # best direction for the B3LYP gap at fixed model energy
    cosF = float(xs @ F @ xB / np.sqrt(xs @ F @ xs) / np.sqrt(xB @ F @ xB))
    cosE = float(gM @ gB / np.linalg.norm(gM) / np.linalg.norm(gB))
    D = json.load(open('directions.json'))[name]; P = np.array(D['P'])
    assert abs(np.dot(P, xs)/np.linalg.norm(P)/np.linalg.norm(xs) - 1) < 1e-8
    pts = json.load(open(f'{name}_points.json'))
    lin = -float(gB @ (np.sqrt(0.5)*P))                      # predicted gap decrease at model energy 0.5 eV
    meas = pts['ref']['gap'] - pts['P@0.5']['gap']
    R = [np.array(r) for r in D['R']]
    cr = [abs(float(r @ F @ xB)/np.sqrt(r @ F @ r)/np.sqrt(xB @ F @ xB)) for r in R]
    out[name] = dict(cos_F=cosF, cos_euclid=cosE, grad_B=gB.tolist(), grad_model=gM.tolist(),
                     lin_pred_dgap_P05=lin, meas_dgap_P05=meas, cos_F_random=cr, h=G['h_A'])
    print(f'{name}: cos_F(x*, best B3LYP direction) = {cosF:.3f}; Euclidean cos(gap gradients) = {cosE:.3f}')
    print(f'   random patterns |cos_F|: median {np.median(cr):.3f}, max {max(cr):.3f}')
    print(f'   linear check at P@0.5: predicted gap decrease {lin:.3f} eV, computed {meas:.3f} eV')
json.dump(out, open('gapgrad_analysis.json', 'w'), indent=1)
