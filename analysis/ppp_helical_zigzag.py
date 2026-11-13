"""Section 7: lowest PPP-ROHF (Ohno) triplets of helical and zigzag [4]phenylene, both on the regular
planar layout of the gapclosure package and with the multi-start protocol of ppp_bs.py."""
import os, sys, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..', 'gapclosure', 'src'))
os.chdir(here)
import ppp_bs as B
from gapclosure.skeleton import phenylene_tree
for name, links in (('helical', [[0, 0, 1], [0, 2, 2], [1, 2, 3]]), ('zigzag', [[0, 0, 1], [0, 2, 2], [1, 4, 3]])):
    sk = phenylene_tree(links); pos = np.asarray(sk.pos, float)
    pos = pos / np.median([np.linalg.norm(pos[a] - pos[b]) for a, b in sk.bonds])
    d = dict(id=name, N=4, rho=0, cp=[[a, b, 1] for a, b in sk.bonds])
    lowest, plain, _ = B.t1_s0(d, pos.tolist())
    print(f'{name} [4]phenylene: lowest triplet {lowest:.3f} eV (plain start {plain:.3f} eV)')
