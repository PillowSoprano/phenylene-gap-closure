from pathlib import Path
import json,time
import numpy as np
from ssh_core import linear_phenylene,KAPPA
HERE=Path(__file__).resolve().parent
out=[];start=time.time();rng=np.random.default_rng(1109)
for n in [2,3,4,6,8,12,16,24,32,48,64]:
 model=linear_phenylene(n);ts,rs=model.relax();es=rs['energy']
 ts2,rs2=model.relax(start=ts+rng.normal(0,.02,len(ts)))
 if rs2['energy']<es-1e-8:ts,rs,es=ts2,rs2,rs2['energy']
 tt,rt=model.relax(True,start=ts)
 tt2,rt2=model.relax(True,start=ts+rng.normal(0,.02,len(ts)))
 if rt2['energy']<rt['energy']:tt,rt=tt2,rt2
 candidates=[model.close(ts,st)for st in [ts,tt,ts+rng.normal(0,.03,len(ts))]]
 td,rd=min(candidates,key=lambda z:z[1]['distance'])
 energy_candidates=[model.close(ts,st,energy=True)for st in [td,tt]]
 te,re=min(energy_candidates,key=lambda z:z[1]['energy'])
 s=model.svd(ts)[1];C=model.electronic_curvature(ts)
 row=dict(N=n,ET1=rt['energy']-es,U_distance=rd['energy']-es,U_energy=re['energy']-es,
          ratio_distance=(rt['energy']-es)/(rd['energy']-es),ratio_energy=(rt['energy']-es)/(re['energy']-es),
          gapS=float(2*s[-1]),sigma_second=float(s[-2]),Hessian_min=float(KAPPA-np.linalg.eigvalsh(C)[-1]),
          tS=ts.tolist(),tT=tt.tolist(),t_distance=td.tolist(),t_energy=te.tolist(),S0=rs,T1=rt,
          distance_point=rd,energy_point=re)
 out.append(row);(HERE/'linear_scaling.json').write_text(json.dumps(out,indent=2)+'\n')
 print(n,row['gapS'],row['ratio_distance'],row['ratio_energy'],row['Hessian_min'],rd['sigma_min'],re['sigma_min'],time.time()-start,flush=True)

