"""Tree theorem: exact regularity test for phenylenes whose hexagons form a tree.

The decision procedure makes one
pass from the leaves to the root with exact integer arithmetic, keeping two
extreme states per hexagon. Ports 0..5 are successive bonds of a hexagon; a link
(i, a, j, b) joins port a of hexagon i and port b of hexagon j through a
four-membered ring.
"""
from fractions import Fraction
from itertools import product
from math import gcd


def hex_pair(weights, branches, port):
    """(D, dD/dt_port) at one hexagon. Branch tuples are (port, connector product, D, U)."""
    factors = list(weights)
    even_scale = odd_scale = 1
    for edge, connector, det, derivative in branches:
        factors[edge] = weights[edge] * det - connector * derivative
        if edge % 2 == 0:
            odd_scale *= det
        else:
            even_scale *= det
    even = even_scale * factors[0] * factors[2] * factors[4]
    odd = odd_scale * factors[1] * factors[3] * factors[5]
    if port is None:
        return even + odd, 0
    # The parent port has no child, so the division is exact, even for signed terms.
    derivative = (odd if port % 2 else even) // weights[port]
    return even + odd, derivative


class PhenyleneTree:
    def __init__(self, size, links, root=0):
        if size < 1 or not 0 <= root < size:
            raise ValueError("Invalid tree size or root.")
        self.size = size
        self.root = root
        adj = [{} for _ in range(size)]
        for i, a, j, b in links:
            if not (0 <= i < size and 0 <= j < size and i != j and j not in adj[i]):
                raise ValueError("Invalid or repeated link.")
            for node, neighbor, port in [(i, j, a), (j, i, b)]:
                if port not in range(6) or any((port - p) % 6 in (0, 1, 5) for p in adj[node].values()):
                    raise ValueError("Ports on a ring must have disjoint endpoints.")
                adj[node][neighbor] = port
        if sum(map(len, adj)) != 2 * (size - 1):
            raise ValueError("Links must form a tree.")
        self.adj = adj
        self.order = [root]
        self.parent = {root: None}
        for i in self.order:
            for j in adj[i]:
                if j not in self.parent:
                    self.order.append(j)
                    self.parent[j] = i
        if len(self.order) != size:
            raise ValueError("Links are disconnected.")
        self.children = {i: [j for j in adj[i] if j != self.parent[i]] for i in self.order}
        self.connector_bit = {j: 6*size+k for k,j in enumerate(self.order[1:])}
        # Uniform messages on both orientations of each cut. The complement at
        # the reference point gives the same positive test functional at all radii.
        unit = {}
        def unit_pair(i, skip):
            branches = [(p, 1, *unit[j, i]) for j, p in adj[i].items() if j != skip]
            return hex_pair((1,) * 6, branches, adj[i].get(skip))
        for i in reversed(self.order[1:]):
            j = self.parent[i]
            unit[i, j] = unit_pair(i, j)
        for i in self.order:
            for j in self.children[i]:
                unit[i, j] = unit_pair(i, j)
        self.unit = unit
        self.asc = unit_pair(root, None)[0]
        assert self.asc > 0

    def regular(self, eps):
        """True iff every determinant is positive throughout the closed eps-box."""
        return self._decide(eps, False)[0]

    def corner(self, eps):
        """A high-bond mask of a nonpositive corner, or None if the box is regular.

        Bits 6*i+p label hexagon bonds. Each remaining bit (self.connector_bit)
        controls the two connecting bonds to one child simultaneously.
        """
        if not 0 < eps < 1:
            raise ValueError("Corner extraction uses 0 < eps < 1.")
        return self._decide(eps, True)[1]

    def weighted_det(self, ring_weights, connector_products):
        """Evaluate the signed sum for positive integer bond weights."""
        states = {}
        for i in reversed(self.order):
            branches = [(self.adj[i][j], connector_products[j], *states[j])
                        for j in self.children[i]]
            states[i] = hex_pair(ring_weights[i], branches, self.adj[i].get(self.parent[i]))
        return states[self.root][0]

    def _complete_corner(self, high_mask, fixed, low, center, high):
        # A nonpositive partial witness used uniform weights outside the subtree.
        # Each remaining bond or connector product is a separate affine variable.
        def det_at(mask, occupied):
            def value(bit):
                if occupied >> bit & 1:
                    return high if mask >> bit & 1 else low
                return center
            rings = [tuple(value(6*i+p) for p in range(6)) for i in range(self.size)]
            connectors = {j: value(bit)**2 for j,bit in self.connector_bit.items()}
            return self.weighted_det(rings,connectors)
        for bit in range(7*self.size-1):
            if not fixed >> bit & 1:
                fixed |= 1 << bit
                if det_at(high_mask, fixed) > 0:
                    high_mask |= 1 << bit
        assert det_at(high_mask,fixed) <= 0
        return high_mask

    def _decide(self, eps, recover):
        eps = Fraction(eps)
        if not 0 <= eps <= 1:
            raise ValueError("Relative radius must be in [0,1].")
        if eps == 0:
            return True, None
        if eps == 1:
            return False, None
        d, b = eps.denominator, eps.numerator
        low, high = d - b, d + b
        weights = [(w, sum((x == high) << k for k,x in enumerate(w)) if recover else 0)
                   for w in product((low, high), repeat=6)]
        products = (low * low, high * high)
        extremes = {}
        fixed = {}
        for i in reversed(self.order):
            parent = self.parent[i]
            port = self.adj[i].get(parent)
            if parent is not None:
                outside_det, outside_derivative = self.unit[parent, i]
            choices = []
            if recover:
                fixed[i] = 63 << (6*i)
            for j in self.children[i]:
                bit = self.connector_bit[j]
                if recover:
                    fixed[i] |= (1 << bit) | fixed[j]
                choices.append([((self.adj[i][j], w, det, deriv),
                                 mask | (1 << bit) if recover and w == high*high else mask)
                                for w in products for det, deriv, mask in extremes[j]])
            min_vec = max_vec = None
            min_test = max_test = None
            for selected in product(*choices):
                branches = [pair for pair,mask in selected]
                branch_mask = sum(mask for pair,mask in selected) if recover else 0
                for t, local_mask in weights:
                    det, deriv = hex_pair(t, branches, port)
                    mask = branch_mask | (local_mask << (6*i)) if recover else 0
                    if parent is None:
                        if det <= 0:
                            return False, mask if recover else None
                        continue
                    test = outside_det * det - d * outside_derivative * deriv
                    if test <= 0:
                        corner = self._complete_corner(mask, fixed[i], low, d, high) if recover else None
                        return False, corner
                    # Normalizing by this positive functional puts all pairs on
                    # one affine line. Its extrema are ordered by deriv / test.
                    if min_vec is None or deriv * min_test < min_vec[1] * test:
                        min_vec, min_test = (det, deriv, mask), test
                    if max_vec is None or deriv * max_test > max_vec[1] * test:
                        max_vec, max_test = (det, deriv, mask), test
            if parent is not None:
                # Only positive scales are discarded, so subsequent signs agree.
                rays = []
                for det, deriv, mask in (min_vec, max_vec):
                    common = gcd(det, deriv)
                    pair = det // common, deriv // common, mask
                    if not rays or pair[:2] != rays[0][:2]:
                        rays.append(pair)
                extremes[i] = rays
        return True, None

    def bracket(self, bits=35):
        """Return exact dyadic endpoints (lower, upper), lower < rho <= upper."""
        if bits < 1:
            raise ValueError("Use at least one bisection.")
        lower, upper = Fraction(0), Fraction(1)
        for _ in range(bits):
            mid = (lower + upper) / 2
            if self.regular(mid):
                lower = mid
            else:
                upper = mid
        return lower, upper


# ====================================================================== recognition

def _six_cycles(adj):
    """All 6-cycles of a graph of maximum degree 3, as tuples of atoms in ring order."""
    found = set()
    out = []
    for s in range(len(adj)):
        stack = [(s, [s])]
        while stack:
            u, path = stack.pop()
            if len(path) == 6:
                if s in adj[u]:
                    key = frozenset(path)
                    if key not in found:
                        found.add(key)
                        out.append(tuple(path))
                continue
            for v in adj[u]:
                if v > s and v not in path:
                    stack.append((v, path + [v]))
    return out


def recognize(sk):
    """If `sk` is a phenylene whose hexagons form a tree, return (tree, bitmap),
    where bitmap[bit] lists the bonds of `sk` that the corner bit controls.
    Otherwise return None."""
    adj = sk.neighbours()
    if sk.n < 6 or sk.n % 6 or max(map(len, adj)) > 3:
        return None
    cycles = _six_cycles(adj)
    hexbonds = {(min(c[k], c[(k + 1) % 6]), max(c[k], c[(k + 1) % 6])) for c in cycles for k in range(6)}
    connectors = [b for b in sk.bonds if b not in hexbonds]
    # removing the connecting bonds must leave disjoint hexagons covering every atom
    hadj = [[] for _ in range(sk.n)]
    for u, v in hexbonds:
        hadj[u].append(v)
        hadj[v].append(u)
    if any(len(a) != 2 for a in hadj):
        return None
    ring_of = [-1] * sk.n
    rings = []
    for s in range(sk.n):
        if ring_of[s] >= 0:
            continue
        order = [s]
        prev, cur = None, s
        while True:
            nxt = hadj[cur][0] if hadj[cur][0] != prev else hadj[cur][1]
            if nxt == s:
                break
            order.append(nxt)
            prev, cur = cur, nxt
            if len(order) > 6:
                return None
        if len(order) != 6:
            return None
        for a in order:
            ring_of[a] = len(rings)
        rings.append(order)
    port_of = []
    for c in rings:
        port_of.append({frozenset((c[p], c[(p + 1) % 6])): p for p in range(6)})
    pairs = {}
    for u, v in connectors:
        i, j = ring_of[u], ring_of[v]
        if i == j:
            return None
        key = (min(i, j), max(i, j))
        pairs.setdefault(key, []).append((u, v) if i < j else (v, u))
    links = []
    for (i, j), cs in pairs.items():
        if len(cs) != 2:
            return None
        (x1, y1), (x2, y2) = cs
        a = port_of[i].get(frozenset((x1, x2)))
        b = port_of[j].get(frozenset((y1, y2)))
        if a is None or b is None:
            return None
        links.append((i, a, j, b, cs))
    if len(links) != len(rings) - 1:
        return None
    try:
        tree = PhenyleneTree(len(rings), [lk[:4] for lk in links])
    except (ValueError, AssertionError):
        return None
    bitmap = {}
    for i, c in enumerate(rings):
        for p in range(6):
            u, v = c[p], c[(p + 1) % 6]
            bitmap[6 * i + p] = [(min(u, v), max(u, v))]
    for i, a, j, b, cs in links:
        child = j if tree.parent.get(j) == i else i
        bitmap[tree.connector_bit[child]] = [(min(u, v), max(u, v)) for u, v in cs]
    return tree, bitmap
