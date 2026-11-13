"""The closure threshold rho, a closing pattern, and the closing orbital.

rho is the smallest eps such that some choice of bond couplings in
[1 - eps, 1 + eps] (relative to equal couplings) closes the Hueckel HOMO-LUMO
gap, that is, makes det T = 0. A closing pattern is such a choice at eps = rho:
every bond is raised (+rho) or lowered (-rho). The closing orbital is the zero
mode of the Hueckel matrix at the closing pattern.

closure_threshold() picks the method:

1. no minority Kekule structure: the gap cannot close for eps < 1 (rho = 1);
2. exactly one minority structure M0: the one-minority theorem, sum over the other structures
   M of q^{|M0 \\ M|} = 1 with q = (1 - rho)/(1 + rho); the pattern raises the
   bonds of M0 and lowers all others;
3. a phenylene whose hexagons form a tree: the tree theorem, exact rational bisection
   in O(N) operations per step, any size;
4. otherwise, if small enough: Rohn's vertex criterion (exact, exponential);
5. otherwise: a certified bracket [lower, upper] with a pattern attaining upper.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product

import numpy as np

from .kekule import TooManyStructures, exact_det, signed_structures, theorem_a_root
from .skeleton import Skeleton
from . import tree as _tree


@dataclass
class ClosureResult:
    skeleton: Skeleton
    rho: float
    lower: float
    upper: float
    method: str
    pattern: dict                  # bond (u, v) -> +1 (raised) or -1 (lowered); 0 = not changed
    orbital: np.ndarray | None     # weight of the closing orbital on each atom (sums to one)
    K: int | None = None           # number of Kekule structures
    ASC: int | None = None         # algebraic structure count |det T|
    minority: int | None = None    # number of minority structures
    equation: str | None = None    # closure equation in q (the one-minority theorem)
    minority_structure: list | None = None
    gap: float | None = None       # Hueckel gap at equal couplings, in units of |beta|
    notes: list = field(default_factory=list)
    bracket: tuple | None = None   # exact (lower, upper) as Fractions, when the tree theorem was used

    @property
    def exact(self):
        return self.upper - self.lower <= 1e-9 * max(1.0, self.upper)

    @property
    def q(self):
        return (1 - self.rho) / (1 + self.rho)

    @property
    def changes(self):
        """Relative coupling change of every bond at the closing pattern."""
        return {b: s * self.rho for b, s in self.pattern.items()}

    def couplings(self):
        """Bond couplings (relative to the reference) at the closing pattern."""
        return {b: 1 + s * self.rho for b, s in self.pattern.items()}

    def to_dict(self):
        sk = self.skeleton
        d = dict(name=sk.name, atoms=sk.n, bonds=sk.m, rho=self.rho, lower=self.lower,
                 upper=self.upper, exact=self.exact, q=self.q, method=self.method,
                 K=self.K, ASC=self.ASC, minority_structures=self.minority,
                 equation=self.equation, huckel_gap=self.gap,
                 closing_pattern=[dict(atoms=[int(u), int(v)], sign=int(s),
                                       relative_change=float(s * self.rho))
                                  for (u, v), s in self.pattern.items()],
                 closing_orbital=None if self.orbital is None else [float(w) for w in self.orbital],
                 bracket=None if self.bracket is None else [str(x) for x in self.bracket],
                 notes=self.notes)
        if sk.source_atoms is not None:
            d["input_atom_index"] = [int(a) for a in sk.source_atoms]
        return d

    def to_json(self, **kw):
        return json.dumps(self.to_dict(), **kw)

    def plot(self, ax=None, **kw):
        from .plot import plot_closing
        return plot_closing(self, ax=ax, **kw)

    def summary(self):
        sk = self.skeleton
        lines = [f"{sk.name or 'skeleton'}: {sk.n} atoms, {sk.m} bonds"]
        if self.K is not None:
            lines.append(f"  Kekule structures K = {self.K}, ASC = {self.ASC}, "
                         f"minority structures = {self.minority}")
        if self.exact:
            lines.append(f"  rho = {self.rho:.10f}   (q = {self.q:.10f})")
        else:
            lines.append(f"  {self.lower:.10f} <= rho <= {self.upper:.10f}")
        lines.append(f"  method: {self.method}")
        if self.equation:
            lines.append(f"  closure equation: {self.equation}")
        if self.gap is not None:
            lines.append(f"  Hueckel gap at equal couplings: {self.gap:.6f} |beta|")
        if self.pattern and self.rho < 1:
            up = sum(1 for s in self.pattern.values() if s > 0)
            down = sum(1 for s in self.pattern.values() if s < 0)
            lines.append(f"  closing pattern: {up} bonds raised, {down} lowered by {100 * self.rho:.2f}%")
        lines += ["  note: " + n for n in self.notes]
        return "\n".join(lines)

    def __str__(self):
        return self.summary()


# ====================================================================== helpers

def _weighted(T, entry, sigma, eps):
    W = T.copy()
    for (i, j), s in zip(entry, sigma):
        W[i, j] = 1 + s * eps
    return W


def _det_sign(T, entry, sigma, eps):
    s, _ = np.linalg.slogdet(_weighted(T, entry, sigma, eps))
    return s


def _corner_root(T, entry, sigma, lo, hi, base_sign):
    """Root of det along the corner sigma between lo (same sign as at eps = 0) and hi."""
    for _ in range(80):
        mid = (lo + hi) / 2
        if _det_sign(T, entry, sigma, mid) == base_sign:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _orbital(sk, T, rows, cols, entry, sigma, eps):
    W = _weighted(T, entry, sigma, eps)
    U, s, Vt = np.linalg.svd(W)
    w = np.zeros(sk.n)
    w[rows] = U[:, -1] ** 2 / 2
    w[cols] = Vt[-1] ** 2 / 2
    return w


def _equation(exponents):
    c = Counter(exponents)
    terms = []
    for k in sorted(c):
        coef = "" if c[k] == 1 else str(c[k])
        terms.append(f"{coef}q" if k == 1 else f"{coef}q^{k}")
    return " + ".join(terms) + " = 1"


def rohn(T, chunk=1 << 15):
    """Exact regularity radius by Rohn's vertex criterion: 2^(2n-1) eigenproblems.
    Returns (rho, y, z); the closing pattern is sigma_ij = -y_i z_j."""
    n = len(T)
    Ti = np.linalg.inv(T)
    best, arg = 1.0, None
    ys = np.array(list(product([-1, 1], repeat=n)), float)
    zs = np.array(list(product([-1, 1], repeat=n - 1)), float).reshape(-1, n - 1)
    zs = np.hstack([np.ones((len(zs), 1)), zs])
    A = np.einsum("ij,yj,jk->yik", Ti, ys, T)
    zper = max(1, chunk // len(ys))
    for zb in range(0, len(zs), zper):
        Z = zs[zb:zb + zper]
        M = (A[None] * Z[:, None, None, :]).reshape(-1, n, n)
        lam = np.linalg.eigvals(M)
        lam = np.where(np.abs(lam.imag) < 1e-9, lam.real, -np.inf)
        lmax = lam.max(axis=1)
        k = int(np.argmax(lmax))
        if lmax[k] > max(1 / best, 1 + 1e-12):
            best = 1 / lmax[k]
            iz, iy = divmod(k, len(ys))
            arg = (ys[iy].copy(), Z[iz].copy())
    return best, arg


def _local_upper(T, restarts=30, seed=0):
    """Best (y, z) found by sign-flip hill climbing; 1/lambda is a valid upper bound."""
    n = len(T)
    Ti = np.linalg.inv(T)
    rng = np.random.default_rng(seed)

    def lam(y, z):
        ev = np.linalg.eigvals(Ti @ (T * np.outer(y, z)))
        ev = ev.real[np.abs(ev.imag) < 1e-9]
        return ev.max() if len(ev) else -np.inf

    best, arg = -np.inf, None
    for _ in range(restarts):
        y = rng.choice([-1.0, 1.0], n)
        z = rng.choice([-1.0, 1.0], n)
        v = lam(y, z)
        improved = True
        while improved:
            improved = False
            for vec in (y, z):
                for k in range(n):
                    vec[k] *= -1
                    w = lam(y, z)
                    if w > v + 1e-12:
                        v, improved = w, True
                    else:
                        vec[k] *= -1
        if v > best:
            best, arg = v, (y.copy(), z.copy())
    return (1 / best, arg) if best > 1 else (1.0, None)


# ====================================================================== main entry

def resolve(mol) -> Skeleton:
    """Accept a Skeleton, a known name (e.g. 'biphenylene') or a SMILES string."""
    if isinstance(mol, Skeleton):
        return mol
    if isinstance(mol, str):
        from .skeleton import from_smiles, named
        try:
            return named(mol)
        except KeyError:
            return from_smiles(mol)
    raise TypeError("expected a Skeleton, a name or a SMILES string")


def closure_threshold(mol, *, bits: int = 40, exact_limit: int = 10,
                      kekule_limit: int = 200_000, method: str = "auto") -> ClosureResult:
    """Closure threshold, closing pattern and closing orbital of a conjugated hydrocarbon.

    mol          Skeleton, a known name, or SMILES (SMILES needs RDKit)
    bits         bisection steps for phenylene trees (bracket width 2^-bits)
    exact_limit  largest n = atoms/2 for Rohn's exponential method
    kekule_limit stop enumerating Kekule structures beyond this many
    method       'auto', 'theorem-a', 'tree', 'rohn' or 'bounds'
    """
    sk = resolve(mol)
    comps = sk.components()
    if len(comps) > 1:
        parts = [(c, _single(sk.subskeleton(c), bits, exact_limit, kekule_limit, method)) for c in comps]
        c, r = min(parts, key=lambda t: t[1].rho)
        pattern = {b: 0 for b in sk.bonds}
        for (u, v), s in r.pattern.items():
            pattern[(min(c[u], c[v]), max(c[u], c[v]))] = s
        orb = np.zeros(sk.n)
        if r.orbital is not None:
            orb[c] = r.orbital
        notes = r.notes + [f"{len(comps)} separate pi systems; the threshold is that of the "
                           f"one with atoms {c[:6]}{'...' if len(c) > 6 else ''}"]
        return ClosureResult(sk, r.rho, r.lower, r.upper, r.method, pattern, orb,
                             gap=min(p[1].gap for p in parts), notes=notes)
    return _single(sk, bits, exact_limit, kekule_limit, method)


def _single(sk, bits, exact_limit, kekule_limit, method):
    T, rows, cols, entry = sk.biadjacency()
    n = len(T)
    det = exact_det(T)
    gap = float(2 * np.linalg.svd(T, compute_uv=False)[-1]) if n else 0.0
    if det == 0:
        return ClosureResult(sk, 0.0, 0.0, 0.0, "det T = 0: the gap is closed at equal couplings",
                             {b: 0 for b in sk.bonds}, None, ASC=0, gap=0.0)
    base = 1 if det > 0 else -1
    res = dict(ASC=abs(det), gap=gap)

    # ---- Kekule structures and sign classes
    structs = None
    if method in ("auto", "theorem-a"):
        try:
            structs = signed_structures(sk, limit=kekule_limit)
        except TooManyStructures:
            if method == "theorem-a":
                raise
    if structs is not None:
        minor = [M for s, M in structs if s < 0]
        res.update(K=len(structs), minority=len(minor))
        if not minor:
            return ClosureResult(sk, 1.0, 1.0, 1.0,
                                 "no minority Kekule structure: the gap cannot close for eps < 1",
                                 {b: 0 for b in sk.bonds}, None, **res)
        if len(minor) == 1:
            M0 = minor[0]
            exps = [len(M0 - M) for s, M in structs if s > 0]
            q = theorem_a_root(exps)
            rho = (1 - q) / (1 + q)
            sigma = {b: (1 if b in M0 else -1) for b in sk.bonds}
            orb = _orbital(sk, T, rows, cols, entry, [sigma[b] for b in sk.bonds], rho)
            return ClosureResult(sk, rho, rho, rho, "single minority structure (exact)", sigma,
                                 orb, equation=_equation(exps), minority_structure=sorted(M0), **res)
        if method == "theorem-a":
            raise ValueError(f"{len(minor)} minority structures; the one-minority theorem needs exactly one")

    # ---- phenylene tree
    if method in ("auto", "tree"):
        rec = _tree.recognize(sk)
        if rec is None and method == "tree":
            raise ValueError("not a phenylene whose hexagons form a tree")
        if rec is not None:
            tr, bitmap = rec
            lo, hi = tr.bracket(bits)
            if hi == 1 and tr.regular(Fraction(1) - Fraction(1, 2 ** bits)):
                return ClosureResult(sk, 1.0, float(lo), 1.0, "phenylene tree (exact)",
                                     {b: 0 for b in sk.bonds}, None, **res)
            mask = tr.corner(hi)
            sign = {}
            for bit, bonds in bitmap.items():
                for b in bonds:
                    sign[b] = 1 if mask >> bit & 1 else -1
            sigma = [sign[b] for b in sk.bonds]
            rho = _corner_root(T, entry, sigma, float(lo), float(hi), base)
            orb = _orbital(sk, T, rows, cols, entry, sigma, rho)
            r = ClosureResult(sk, rho, float(lo), float(hi), "phenylene tree (exact)",
                              dict(zip(sk.bonds, sigma)), orb, bracket=(lo, hi), **res)
            r.notes.append(f"exact bracket {lo} < rho <= {hi} (width 2^-{bits}); rho is the root of "
                           "the closing corner inside it")
            return r

    # ---- Rohn's vertex criterion
    if method == "rohn" or (method == "auto" and n <= exact_limit):
        rho, (y, z) = rohn(T)
        sigma = [-int(y[i] * z[j]) for i, j in entry]
        orb = _orbital(sk, T, rows, cols, entry, sigma, rho)
        return ClosureResult(sk, rho, rho, rho, "vertex enumeration (Rohn)",
                             dict(zip(sk.bonds, sigma)), orb, **res)

    # ---- certified bracket
    Ti = np.linalg.inv(T)
    lower = 1 / max(abs(np.linalg.eigvals(np.abs(Ti) @ np.abs(T))))
    upper, arg = _local_upper(T)
    sigma = None if arg is None else [-int(arg[0][i] * arg[1][j]) for i, j in entry]
    if structs is not None:
        for s, M0 in structs:
            if s > 0:
                continue
            sg = [1 if b in M0 else -1 for b in sk.bonds]
            if _det_sign(T, entry, sg, upper) != base:
                r_c = _corner_root(T, entry, sg, 0.0, upper, base)
                if r_c < upper:
                    upper, sigma = r_c, sg
    if sigma is None:
        return ClosureResult(sk, upper, lower, upper, "bounds", {b: 0 for b in sk.bonds}, None, **res,
                             notes=["no closing corner found below eps = 1"])
    orb = _orbital(sk, T, rows, cols, entry, sigma, upper)
    return ClosureResult(sk, upper, lower, upper, "bounds", dict(zip(sk.bonds, sigma)), orb, **res,
                         notes=["lower bound from rho(|T^-1||T|), upper bound from the best corner "
                                "found; the pattern attains the upper bound",
                                f"exact methods: the one-minority theorem, the tree theorem, or Rohn up to {exact_limit} "
                                "starred atoms (raise exact_limit at exponential cost)"])
