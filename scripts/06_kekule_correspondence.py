import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import networkx as nx
from networkx.algorithms import isomorphism as iso
from gapclosure import phenylene, bipartite, kekule_terms
def inc_graph(g):
    T, L, R, bix = bipartite(g); terms = kekule_terms(T)
    G = nx.Graph()
    for k, e in enumerate(bix): G.add_node(('e', k), kind='edge')
    for m, (s, p) in enumerate(terms):
        G.add_node(('m', m), kind='plus' if s > 0 else 'minus')
        for k, (i, j) in enumerate(bix):
            if p[i] == j: G.add_edge(('m', m), ('e', k))
    return G
Gz, Gh = inc_graph(phenylene(['L','L'])), inc_graph(phenylene(['L','R']))
gm = iso.GraphMatcher(Gz, Gh, node_match=lambda a, b: a['kind'] == b['kind'])
print('signed Kekule structures of zigzag[4] and helical[4] identical up to relabelling of bonds:', gm.is_isomorphic())
Gl, Gb = inc_graph(phenylene(['p','p'])), inc_graph(phenylene(['p','L']))
print('control (linear[4] vs bent[4]):', iso.GraphMatcher(Gl, Gb, node_match=lambda a, b: a['kind'] == b['kind']).is_isomorphic())
