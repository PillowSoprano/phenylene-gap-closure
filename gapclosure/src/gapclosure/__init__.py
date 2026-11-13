"""gapclosure: Hueckel gap-closure thresholds of conjugated hydrocarbons.

    >>> import gapclosure as gc
    >>> r = gc.closure_threshold("biphenylene")
    >>> round(r.rho, 6)
    0.188829
    >>> r.equation
    'q^2 + 2q^4 + q^6 = 1'
"""
from .skeleton import (Skeleton, acene, acene_chain, annulene, benzenoid, benzocyclobutadiene,
                       biphenylene, from_edges, from_molfile, from_points, from_smiles, from_xyz,
                       named, phenylene, phenylene_tree, triangular_phenylene)
from .kekule import exact_det, signed_structures
from .threshold import ClosureResult, closure_threshold, rohn
from .tree import PhenyleneTree, recognize as recognize_phenylene_tree

__version__ = "0.1.0"

__all__ = [
    "Skeleton", "ClosureResult", "closure_threshold", "PhenyleneTree", "recognize_phenylene_tree",
    "acene", "acene_chain", "annulene", "benzenoid", "benzocyclobutadiene", "biphenylene",
    "phenylene", "phenylene_tree", "triangular_phenylene", "named",
    "from_smiles", "from_molfile", "from_xyz", "from_edges", "from_points",
    "signed_structures", "exact_det", "rohn",
]
