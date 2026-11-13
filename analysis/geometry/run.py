import numpy as np, itertools, networkx as nx

def benzocyclobutadiene():
    G=nx.cycle_graph(6); # ring 0..5
    G.add_edges_from([(0,6),(6,7),(7,1)])
    return G

def phenylene(modes):
    # modes: list for internal hexagons, 'L' or 'A'
    G=nx.Graph(); hexes=[]
    def add_hex(start):
        h=list(range(start,start+6))
        for i in range(6): G.add_edge(h[i],h[(i+1)%6])
        return h
    h=add_hex(0); hexes.append(h); nxt=6
    attach=(h[3],h[4])  # first hexagon: any bond away from nothing
    seq=list(modes)+['END']
    for k,mode in enumerate(seq):
        h2=add_hex(nxt); nxt+=6
        a,b=attach
        # four ring a-b-h2[0]-h2[1]-a
        G.add_edge(b,h2[0]); G.add_edge(a,h2[1])
        if mode=='L': attach=(h2[3],h2[4])
        elif mode=='A': attach=(h2[2],h2[3])
    return G

def Tmat(G,w):
    col=nx.bipartite.color(G)
    S=[v for v in G if col[v]==0]; U=[v for v in G if col[v]==1]
    assert len(S)==len(U)
    si={v:i for i,v in enumerate(S)}; ui={v:i for i,v in enumerate(U)}
    T=np.zeros((len(S),len(U)))
    for (i,j) in G.edges():
        if col[i]==1: i,j=j,i
        T[si[i],ui[j]]=w[frozenset((i,j))]
    return T,si,ui,col

def rho(T):
    n=T.shape[0]; Ti=np.linalg.inv(T); best=0
    signs=np.array(list(itertools.product([-1,1],repeat=n)))
    ys=signs; zs=signs[signs[:,0]==1]
    for z in zs:
        # T∘ y z^T : entry T_ij y_i z_j
        M=T[None,:,:]*ys[:,:,None]*z[None,None,:]
        A=Ti[None]@M
        ev=np.linalg.eigvals(A)
        re=np.where(np.abs(ev.imag)<1e-9,np.abs(ev.real),0)
        best=max(best,re.max())
    return 1/best

def relax(G,t0=2.5,al=4.1,fs=43.,pref=2/3,it=500):
    E=[frozenset(e) for e in G.edges()]
    d={e:0.0 for e in E}
    nodes=sorted(G); idx={v:i for i,v in enumerate(nodes)}
    for _ in range(it):
        H=np.zeros((len(nodes),)*2)
        for e in E:
            i,j=tuple(e); t=t0-al*d[e]; H[idx[i],idx[j]]=H[idx[j],idx[i]]=-t
        w,v=np.linalg.eigh(H); occ=v[:,:len(nodes)//2]
        new={}
        for e in E:
            i,j=tuple(e); p=2*np.dot(occ[idx[i]],occ[idx[j]])
            new[e]=2*al*(pref-p)/fs
        diff=max(abs(new[e]-d[e]) for e in E)
        for e in E: d[e]=0.5*d[e]+0.5*new[e]
        if diff<1e-12: break
    return {e:t0-al*d[e] for e in E}, d

cases={'benzocyclobutadiene':benzocyclobutadiene(),'biphenylene':phenylene([]),
       'linear[3]':phenylene(['L']),'angular[3]':phenylene(['A'])}
for name,G in cases.items():
    w1={frozenset(e):1.0 for e in G.edges()}
    T1,*_=Tmat(G,w1)
    wr,d=relax(G)
    Tr,*_=Tmat(G,wr)
    gap1=2*np.linalg.svd(T1,compute_uv=False).min()
    gapr=2*np.linalg.svd(Tr,compute_uv=False).min()
    print(f"{name:20s} atoms={G.number_of_nodes()} bonds={G.number_of_edges()} rho_equal={rho(T1):.5f} rho_relaxed={rho(Tr):.5f} gap_relaxed={gapr:.3f}eV  t range {min(wr.values()):.3f}-{max(wr.values()):.3f} eV  dmin/dmax {min(d.values()):+.3f}/{max(d.values()):+.3f} A")
