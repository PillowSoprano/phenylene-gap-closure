"""Positive central closure via the unique-minority matching polynomial, not sigma_min."""
from pathlib import Path
import json,time
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize
from ssh_core import acene_pair,forcing_minority
HERE=Path(__file__).resolve().parent

def matchings(edges,indices):
    vertices={v for k in indices for v in edges[k]};adj={v:[]for v in vertices}
    for k in indices:
        a,b=edges[k];adj[a].append((b,k));adj[b].append((a,k))
    def rec(free):
        if not free:return [[]]
        a=min(free,key=lambda a:sum(b in free for b,k in adj[a]));out=[]
        for b,k in adj[a]:
            if b in free:out.extend([[k]+x for x in rec(free-{a,b})])
        return out
    return rec(vertices)

class AceneClosure:
    def __init__(self,n):
        self.n=n;self.model=acene_pair(n,n);e=len(self.model.edges);unit=5*n+1
        self.M=np.zeros(e);self.M[forcing_minority(self.model)]=1
        self.parts=[]
        for ids in [range(unit),range(unit,2*unit)]:
            mm=matchings(self.model.edges,ids);assert len(mm)==n+1
            I=np.zeros((len(mm),e))
            for j,ms in enumerate(mm):I[j,ms]=1
            self.parts.append(I)

    def phi(self,logt):
        value=-self.M@logt;grad=-self.M.copy()
        for I in self.parts:
            z=I@logt;l=logsumexp(z);value+=l;grad+=np.exp(z-l)@I
        return float(value),grad

    def close(self,ts,start,energy=False):
        def obj(x):
            t=np.exp(x)
            if energy:f,g=self.model.energy_gradient(t,True)
            else:f=.5*np.sum((t-ts)**2);g=t-ts
            return f,t*g
        if energy:bounds=[(np.log(1e-8),np.log(20))]*len(ts)
        else:bounds=[(np.log(t),np.log(20))if m else(np.log(1e-8),np.log(t))for t,m in zip(ts,self.M)]
        constraint=dict(type='eq',fun=lambda x:self.phi(x)[0],jac=lambda x:self.phi(x)[1])
        r=minimize(obj,np.log(np.maximum(start,1e-8)),jac=True,method='SLSQP',constraints=[constraint],bounds=bounds,
                   options=dict(ftol=1e-11,maxiter=1500))
        t=np.exp(r.x);phi,gphi=self.phi(r.x);_,g=obj(r.x);lam=-g@gphi/(gphi@gphi)
        s,n=self.model.signed_sigma_gradient(t);g_t=g/t;phi_t=gphi/t;lamt=-g_t@phi_t/(phi_t@phi_t)
        sigma_multiplier=-g_t@n/(n@n)
        return t,dict(success=bool(r.success),message=str(r.message),nit=int(r.nit),log_partition_residual=phi,
             sigma_min=float(abs(s)),sigma_second=float(self.model.svd(t)[1][-2]),distance=float(np.linalg.norm(t-ts)),
             energy=float(self.model.energy_gradient(t)[0]),kkt_log_residual=float(np.linalg.norm(g+lam*gphi)),
             kkt_residual=float(np.linalg.norm(g_t+lamt*phi_t)),sigma_multiplier=float(sigma_multiplier),
             min_coupling=float(t.min()),max_coupling=float(t.max()))

if __name__=='__main__':
    rng=np.random.default_rng(1114);out=[];beg=time.time()
    for n in [1,2,4,8,16,32,64]:
        cl=AceneClosure(n);model=cl.model;ts,rs=model.relax();es=rs['energy']
        tt,rt=model.relax(True,start=ts+rng.normal(0,.02,len(ts)))
        # The analytic minority ray always supplies a positive central closing point.
        from scipy.optimize import brentq
        sign=2*cl.M-1
        rho=brentq(lambda e:cl.phi(np.log(1+e*sign))[0],0,.4)
        tc=1+rho*sign;tc*=tc@ts/(tc@tc)
        td,rd=cl.close(ts,tc)
        trials=[cl.close(ts,td,energy=True)]
        for _ in range(3):trials.append(cl.close(ts,td*np.exp(rng.normal(0,.08,len(ts))),energy=True))
        te,re=min(trials,key=lambda z:z[1]['energy'])
        row=dict(n=n,atoms=model.n,rho=rho,edges=model.edges,minority_matching=np.flatnonzero(cl.M).tolist(),
            tS=ts.tolist(),tT=tt.tolist(),t_distance=td.tolist(),t_energy=te.tolist(),S0=rs,T1=rt,
            distance_point=rd,energy_point=re,energy_trials=[r for t,r in trials],
            gapS=float(2*model.svd(ts)[1][-1]),ET1=rt['energy']-es,U_distance=rd['energy']-es,U_energy=re['energy']-es,
            ratio_distance=(rt['energy']-es)/(rd['energy']-es),ratio_energy=(rt['energy']-es)/(re['energy']-es),
            distance_sign_agreement=int(np.sum(np.sign(td-ts)==sign)),energy_sign_agreement=int(np.sum(np.sign(te-ts)==sign)))
        out.append(row);(HERE/'acene_matching_closure.json').write_text(json.dumps(out,indent=2)+'\n')
        print(n,'gap',row['gapS'],'E',row['ET1'],'U',row['U_energy'],'ratio',row['ratio_energy'],'lambda',re['sigma_multiplier'],'KKT',re['kkt_residual'],'distance',rd['distance'],'seconds',time.time()-beg,flush=True)
