"""Evaluate the global polar minorant and a posteriori T1/closure ratio bound."""
from pathlib import Path
import json
import numpy as np
from ssh_core import Model,linear_phenylene
HERE=Path(__file__).resolve().parent
def bound(model,b,tc):
    eb,res=model.energy_gradient(b);enorm=np.linalg.norm(res)
    G=model.global_growth_matrix(b);mu=float(np.linalg.eigvalsh(G)[0]);g=float(2*model.svd(b)[1][-1])
    out=dict(mu=mu,reference_gap=g,stationarity_residual=float(enorm))
    if mu<=0:return dict(**out,status='this minorant is not positive definite')
    # Valid at an arbitrary reference b, avoiding an exact-stationarity assumption.
    if enorm/mu>=g/2:lower=0.
    elif (2+enorm)/mu<=g/2:lower=g-(2+enorm)**2/(2*mu)
    else:lower=mu*g*g/8-enorm*g/2
    upper=model.energy_gradient(tc)[0]-eb+enorm**2/(2*mu)
    out.update(status='positive global minorant',T1_lower=float(lower),closure_upper=float(upper),
               ratio_lower=float(lower/upper),ground_energy_bracket=[float(eb-enorm**2/(2*mu)),float(eb)])
    return out
if __name__=='__main__':
    rows=[]
    for r in json.loads((HERE/'linear_scaling.json').read_text()):
        model=linear_phenylene(r['N']);rows.append(dict(family='linear',N=r['N'],**bound(model,np.array(r['tS']),np.array(r['t_energy']))))
    for r in json.loads((HERE/'objective_comparison.json').read_text()):
        model=Model(r['edges']);rows.append(dict(family='selected_tree',id=r['id'],**bound(model,np.array(r['tS']),np.array(r['t_energy']))))
    # Check the identity-derived global inequality well beyond the local basin.
    rng=np.random.default_rng(1115);tests=0;smallest_margin=float('inf')
    for r in json.loads((HERE/'objective_comparison.json').read_text()):
        model=Model(r['edges']);b=np.array(r['tS']);fb,grad=model.energy_gradient(b);G=model.global_growth_matrix(b)
        for scale in [.01,.1,1.,3.]:
            for _ in range(10):
                delta=rng.normal(size=len(b))*scale
                margin=model.energy_gradient(b+delta)[0]-fb-grad@delta-.5*delta@G@delta
                assert margin>=-1e-9;smallest_margin=min(smallest_margin,margin);tests+=1
    out=dict(status='formula evaluations and numerical identity checks passed',arithmetic='floating point; proof in PROOFS.md',
             rows=rows,global_inequality_checks=tests,smallest_check_margin=float(smallest_margin))
    (HERE/'global_lower_bound.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows:print(r.get('id',r.get('N')),r['mu'],r.get('ratio_lower'),flush=True)
