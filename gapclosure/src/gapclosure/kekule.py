"""Signed Kekule structures.

det T is the signed sum of the Kekule structures; a structure's sign is the sign
of its permutation. The sign that wins in det T is the majority sign; the other
structures are minority structures.
"""
from __future__ import annotations

from fractions import Fraction


class TooManyStructures(RuntimeError):
    pass


def exact_det(T) -> int:
    """Exact determinant of a 0/1 (or integer) matrix by fraction-free elimination."""
    A = [[int(round(x)) for x in row] for row in T]
    n = len(A)
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    A[k], A[r] = A[r], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1] if n else 1


def _perm_sign(p):
    s, seen = 1, [False] * len(p)
    for a in range(len(p)):
        if not seen[a]:
            L, b = 0, a
            while not seen[b]:
                seen[b] = True
                b = p[b]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def signed_structures(sk, limit: int = 200_000):
    """All Kekule structures as (sign, frozenset of bonds), sign +1 for the majority.

    Raises TooManyStructures beyond `limit`."""
    T, rows, cols, entry = sk.biadjacency()
    n = len(rows)
    bond_of = {e: b for b, e in zip(sk.bonds, entry)}
    opts = [[j for j in range(n) if T[i, j]] for i in range(n)]
    out = []
    perm = [-1] * n
    used = [False] * n

    def rec(done):
        if done == n:
            if len(out) >= limit:
                raise TooManyStructures(f"more than {limit} Kekule structures")
            out.append((_perm_sign(perm), frozenset(bond_of[(i, perm[i])] for i in range(n))))
            return
        # row with fewest free columns
        best, bopts = None, None
        for i in range(n):
            if perm[i] < 0:
                o = [j for j in opts[i] if not used[j]]
                if best is None or len(o) < len(bopts):
                    best, bopts = i, o
                    if len(o) <= 1:
                        break
        if not bopts:
            return
        for j in bopts:
            perm[best] = j
            used[j] = True
            rec(done + 1)
            used[j] = False
        perm[best] = -1

    rec(0)
    total = sum(s for s, _ in out)
    flip = -1 if total < 0 else 1
    return [(s * flip, M) for s, M in out]


def one_kekule(sk):
    """One Kekule structure as a set of bonds, or None."""
    try:
        T, rows, cols, entry = sk.biadjacency()
    except ValueError:
        return None
    n = len(rows)
    adj = [[j for j in range(n) if T[i, j]] for i in range(n)]
    match_c = [-1] * n

    def aug(i, seen):
        for j in adj[i]:
            if not seen[j]:
                seen[j] = True
                if match_c[j] < 0 or aug(match_c[j], seen):
                    match_c[j] = i
                    return True
        return False

    for i in range(n):
        if not aug(i, [False] * n):
            return None
    bond_of = {e: b for b, e in zip(sk.bonds, entry)}
    return {bond_of[(match_c[j], j)] for j in range(n)}


def theorem_a_root(exponents) -> Fraction | float:
    """Root q in (0, 1) of sum_k q^k = 1 (one-minority theorem), by bisection to double precision."""
    exps = list(exponents)
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if sum(mid ** k for k in exps) < 1:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-17:
            break
    return (lo + hi) / 2
