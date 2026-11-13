"""Real bipartite Hückel--SSH in coupling coordinates; all stored bonds are independent."""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize

R0, X, Y, RREF = 1.523, .189, .2756, 1.397
KAPPA = 2*Y/X
REST = 1+(RREF-R0)/Y

@dataclass
class Model:
    edges: list
    kappa: float = KAPPA
    rest: object = REST

    def __post_init__(self):
        self.n = max(max(e) for e in self.edges)+1
        adj=[[] for _ in range(self.n)]
        for a,b in self.edges: adj[a].append(b);adj[b].append(a)
        colors={}
        for seed in range(self.n):
            if seed in colors:continue
            colors[seed]=0; stack=[seed]
            while stack:
                a=stack.pop()
                for b in adj[a]:
                    if b not in colors:colors[b]=1-colors[a];stack.append(b)
                    assert colors[b]!=colors[a]
        aa=[a for a in range(self.n)if colors[a]==0];bb=[b for b in range(self.n)if colors[b]==1]
        assert len(aa)==len(bb);self.m=len(aa)
        ai={a:i for i,a in enumerate(aa)};bi={b:i for i,b in enumerate(bb)}
        self.rr=np.array([ai[a]if colors[a]==0 else ai[b]for a,b in self.edges])
        self.cc=np.array([bi[b]if colors[a]==0 else bi[a]for a,b in self.edges])
        self.rest=np.broadcast_to(self.rest,(len(self.edges),)).copy()
        self.refsign=np.linalg.slogdet(self.matrix(np.ones(len(self.edges))))[0]
        if self.refsign==0:self.refsign=1

    def matrix(self,t):
        out=np.zeros((self.m,self.m));out[self.rr,self.cc]=t;return out

    def svd(self,t):return np.linalg.svd(self.matrix(t),full_matrices=False)

    def energy_gradient(self,t,triplet=False):
        u,s,vt=self.svd(t);stop=self.m-int(triplet)
        p=np.sum(u[self.rr,:stop]*vt[:stop,self.cc].T,axis=1)
        return float(-2*s[:stop].sum()+self.kappa/2*np.sum((t-self.rest)**2)),self.kappa*(t-self.rest)-2*p

    def signed_sigma_gradient(self,t):
        u,s,vt=self.svd(t)
        orientation=self.refsign*np.linalg.slogdet(u)[0]*np.linalg.slogdet(vt)[0]
        return float(orientation*s[-1]),orientation*u[self.rr,-1]*vt[-1,self.cc]

    def relax(self,triplet=False,start=None):
        if start is None:start=np.ones(len(self.edges))
        # Unrestricted coupling coordinates, matching the given harmonic model.
        result=minimize(lambda t:self.energy_gradient(t,triplet),start,jac=True,method='BFGS',
                        options=dict(gtol=2e-10,maxiter=3000))
        return result.x,dict(energy=float(result.fun),gradient_norm=float(np.linalg.norm(result.jac)),
                            success=bool(result.success),message=str(result.message),nit=int(result.nit))

    def close(self,ts,start,energy=False):
        if energy:
            fun=lambda t:self.energy_gradient(t,True) # E1=E0 on the singular surface; smooth for simple nullity.
        else:fun=lambda t:(.5*np.sum((t-ts)**2),t-ts)
        cons=dict(type='eq',fun=lambda t:self.signed_sigma_gradient(t)[0],
                  jac=lambda t:self.signed_sigma_gradient(t)[1])
        res=minimize(fun,start,jac=True,method='SLSQP',constraints=[cons],
                     bounds=[(1e-9,None)]*len(ts),options=dict(ftol=1e-12,maxiter=1000))
        t=res.x;f,normal=self.signed_sigma_gradient(t);_,grad=fun(t)
        lam=-grad@normal/(normal@normal)
        tangent=grad+lam*normal
        u,s,vt=self.svd(t)
        return t,dict(success=bool(res.success),message=str(res.message),nit=int(res.nit),
                      sigma_min=float(s[-1]),sigma_second=float(s[-2]),
                      kkt_residual=float(np.linalg.norm(tangent)),distance=float(np.linalg.norm(t-ts)),
                      energy=float(self.energy_gradient(t)[0]),min_coupling=float(t.min()),multiplier=float(lam))

    def electronic_curvature(self,t):
        u,s,vt=self.svd(t);rows=[]
        for i in range(self.m):
            for j in range(i+1,self.m):
                b=u[self.rr,i]*vt[j,self.cc]-u[self.rr,j]*vt[i,self.cc]
                rows.append(np.sqrt(2/(s[i]+s[j]))*b)
        B=np.array(rows)
        return B.T@B

    def global_growth_matrix(self,t,theta=.5):
        """Global anchored quadratic minorant from the polar variational formula."""
        u,s,vt=self.svd(t)
        B=(u[self.rr,:,None]*vt.T[self.cc,None,:]).reshape(len(t),-1).T
        denominators=(theta*s[:,None]+(1-theta)*s[None,:]).ravel()
        return self.kappa*np.eye(len(t))-2*(B.T/denominators)@B

def linear_phenylene(n):
    edges=[(6*j+i,6*j+(i+1)%6)for j in range(n)for i in range(6)]
    for j in range(n-1):edges.extend([(6*j+3,6*(j+1)),(6*j+4,6*(j+1)+1)])
    return Model(edges)

def acene_pair(left,right):
    def unit(n):
        parent=list(range(6*n))
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for i in range(n-1):
            parent[root(6*i+3)]=root(6*(i+1)+1)
            parent[root(6*i+4)]=root(6*(i+1))
        labels={v:i for i,v in enumerate(sorted({root(i)for i in range(6*n)}))}
        edges=sorted({tuple(sorted((labels[root(6*j+i)],labels[root(6*j+(i+1)%6)])))for j in range(n)for i in range(6)})
        return edges,(labels[root(0)],labels[root(1)]),len(labels)
    e1,p1,n1=unit(left);e2,p2,n2=unit(right)
    edges=e1+[(a+n1,b+n1)for a,b in e2]+[(p1[j],p2[j]+n1)for j in range(2)]
    return Model(edges)

def penalty_minimum(model,start,eta=1):
    """Minimize F1+eta*g on the two determinant-sign branches; local numerical solutions."""
    out=[]
    for sign in [1,-1]:
        def objective(t):
            f,g=model.energy_gradient(t,True);s,n=model.signed_sigma_gradient(t)
            return f+2*eta*sign*s,g+2*eta*sign*n
        cons=dict(type='ineq',fun=lambda t:sign*model.signed_sigma_gradient(t)[0],
                  jac=lambda t:sign*model.signed_sigma_gradient(t)[1])
        r=minimize(objective,start,jac=True,method='SLSQP',constraints=[cons],
                   bounds=[(1e-9,None)]*len(start),options=dict(ftol=1e-11,maxiter=1000))
        out.append(dict(sign=sign,energy=float(r.fun),gap=float(2*model.svd(r.x)[1][-1]),
                        success=bool(r.success),t=r.x.tolist(),nit=int(r.nit)))
    return out

def forcing_minority(model):
    """Unique matching containing the two bridge bonds of acene_pair()."""
    forced=[len(model.edges)-2,len(model.edges)-1]
    used={v for e in forced for v in model.edges[e]};free=set(range(model.n))-used
    adj=[[]for _ in range(model.n)]
    for k,(a,b)in enumerate(model.edges):adj[a].append((b,k));adj[b].append((a,k))
    def solve(vertices):
        if not vertices:return [[]]
        a=min(vertices,key=lambda a:sum(b in vertices for b,k in adj[a]));out=[]
        for b,k in adj[a]:
            if b in vertices:
                out.extend([[k]+tail for tail in solve(vertices-{a,b})])
        return out
    matches=solve(free);assert len(matches)==1
    return forced+matches[0]

