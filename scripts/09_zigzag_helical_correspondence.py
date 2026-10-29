import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np, networkx as nx, sympy as sp
from networkx.algorithms import isomorphism as iso
from gapclosure import phenylene, triangular4, benzocyclobutadiene, bipartite, kekule_terms
gz, gh = phenylene(['L','L']), phenylene(['L','R'])
for nm,g in (('zigzag LL',gz),('helical LR',gh)):
    P=np.array(g['pos']); C=[P[r].mean(0) for r in g['rings']]
    hexes=[c for c,r in zip(C,g['rings']) if len(r)==6]
    print(nm,'terminal hexagon centre distance',np.linalg.norm(hexes[0]-hexes[-1]).round(3), 'nrings',len(g['rings']))
def inc(g):
    T,L,R,bix=bipartite(g); terms=kekule_terms(T); G=nx.Graph()
    for k in range(len(bix)): G.add_node(('e',k),kind='edge')
    for m,(s,p) in enumerate(terms):
        G.add_node(('m',m),kind='plus' if s>0 else 'minus')
        for k,(i,j) in enumerate(bix):
            if p[i]==j: G.add_edge(('m',m),('e',k))
    return G,T,bix,terms
Gz,Tz,bz,tz=inc(gz); Gh,Th,bh,th=inc(gh)
gm=iso.GraphMatcher(Gz,Gh,node_match=lambda a,b:a['kind']==b['kind'])
print('isomorphic',gm.is_isomorphic())
mp={a[1]:b[1] for a,b in gm.mapping.items() if a[0]=='e'}
# count all correspondences
n_all=sum(1 for _ in iso.GraphMatcher(Gz,Gh,node_match=lambda a,b:a['kind']==b['kind']).isomorphisms_iter())
print('number of correspondences',n_all)
# verify polynomial identity symbolically
x=sp.symbols('x0:30')
def poly(T,bix,terms,sub):
    return sp.expand(sum(s*sp.Mul(*[sub[[k for k,(i,j) in enumerate(bix) if i==ii and j==p[ii]][0]] for ii in range(len(T))]) for s,p in terms))
pz=poly(Tz,bz,tz,x); ph=poly(Th,bh,th,[x[[k for k in mp if mp[k]==kk][0]] for kk in range(30)])
print('polys equal under relabelling:',sp.expand(pz-ph)==0, ' n terms',len(pz.args))
# is the correspondence induced by a graph isomorphism of the skeletons? (it cannot be: different graphs)
Gz0=nx.Graph(gz['edges']); Gh0=nx.Graph(gh['edges']); print('skeleton graphs isomorphic:',nx.is_isomorphic(Gz0,Gh0))
# does the bond map send adjacent bonds to adjacent bonds?
def adjb(g,a,b): return len(set(g['edges'][a])&set(g['edges'][b]))==1
L=nx.line_graph(Gz0)
bad=sum(1 for a in range(30) for b in range(a+1,30) if adjb(gz,a,b)!=adjb(gh,mp[a],mp[b]))
print('pairs where adjacency not preserved',bad)
json.dump(dict(zigzag_edges=gz['edges'],helical_edges=gh['edges'],map={int(k):int(v) for k,v in mp.items()},
  zigzag_pos=gz['pos'],helical_pos=gh['pos'],zig_rings=gz['rings'],hel_rings=gh['rings']),open(ROOT/'results'/'zigzag_helical_correspondence.json','w'))
for k in range(30): print(k, gz['edges'][k], '->', mp[k], gh['edges'][mp[k]])
