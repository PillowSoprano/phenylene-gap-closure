# gapclosure

How far must the bond couplings of a conjugated hydrocarbon move before its Hückel HOMO–LUMO gap closes?

The answer is the **closure threshold ρ**. It is the smallest relative change ε such that some couplings within [1−ε, 1+ε] of the reference make det T = 0, where H = −[[0, T], [Tᵀ, 0]]. `gapclosure` takes a molecule and returns three things:

- **ρ**, exact wherever a theorem allows it;
- a **closing pattern**: one choice of couplings at ε = ρ that closes the gap, with every bond either raised by ρ or lowered by ρ;
- the **closing orbital**: the zero mode at that pattern.

It also draws them in one picture:

![closing pattern of [5]phenylene LAL](docs/closing_LAL.png)

Red bonds are raised and blue bonds are lowered. The circles show the weight of the closing orbital on each atom.

## Install

```bash
pip install gapclosure              # numpy + matplotlib
pip install "gapclosure[rdkit]"     # also read SMILES and molfiles
```

## Use

```bash
gapclosure biphenylene --plot biphenylene.png
gapclosure "c1ccc2c(c1)-c1cc3c(cc1-2)-c1ccccc1-3" --json     # SMILES, needs RDKit
gapclosure --phenylene LALA                                  # [6]phenylene from its L/A word
gapclosure --acenes 2,3                                      # naphthalene–cyclobutadiene–anthracene
gapclosure --xyz geometry.xyz                                # carbons and bonds read from an optimized geometry
```

```python
import gapclosure as gc

r = gc.closure_threshold("biphenylene")
r.rho            # 0.18882895203249...
r.equation       # 'q^2 + 2q^4 + q^6 = 1'
r.pattern        # {(u, v): +1 or -1}  raised / lowered at eps = rho
r.orbital        # weight of the closing orbital on each atom
r.plot()         # the red-blue figure (matplotlib Axes)
r.to_json()

sk = gc.phenylene_tree([(0, 0, 1), (0, 2, 2), (0, 4, 3)])   # any tree of hexagons
gc.closure_threshold(sk)
```

### Input

| Input | How |
|---|---|
| a name | `"biphenylene"`, `"angular[3]phenylene"`, `"triangular[4]phenylene"`, … (`gapclosure --list`) |
| SMILES, molfile | `gc.from_smiles`, `gc.from_molfile` (RDKit). Sp² carbons form the π system; substituents such as methyl groups are ignored. |
| XYZ geometry | `gc.from_xyz`: carbons with at most three neighbours, bonded by distance |
| edge list | `gc.from_edges([(0, 1), (1, 2), ...])` |
| builders | `gc.acene(h)`, `gc.acene_chain([hA, hB, ...])`, `gc.phenylene("LA")`, `gc.phenylene_tree(links)`, `gc.benzenoid(cells)`, `gc.annulene(L)` |

The theory covers alternant (bipartite) hydrocarbons with equal Coulomb integrals. A molecule with heteroatoms or odd rings in the π system raises an error.

## How ρ is computed

`closure_threshold` uses the first method that applies.

| Case | Method | Cost | Result |
|---|---|---|---|
| no Kekulé structure of the minority sign | the gap cannot close for ε < 1 | Kekulé enumeration | ρ = 1 |
| exactly one minority structure M₀ | One-minority theorem: Σ q^{\|M₀∖M\|} = 1 over the other structures, q = (1−ρ)/(1+ρ) | Kekulé enumeration | exact, with the equation |
| phenylene whose hexagons form a tree | Tree theorem: one exact integer pass from the leaves to the root per bisection step | O(N) operations on O(N log d)-bit integers per step | exact rational bracket of width 2⁻⁴⁰, and the root of the closing pattern inside it |
| other skeletons with at most 10 starred atoms | Rohn's vertex criterion | 2^(2n−1) eigenproblems | exact |
| anything else | lower bound 1/ρ(\|T⁻¹\|\|T\|), upper bound from the best pattern found | polynomial | certified bracket |

Phenylene trees are recognized from the graph alone, so a SMILES string or a relabelled edge list goes through the tree theorem as well. The method can be forced with `method="theorem-a" | "tree" | "rohn" | "bounds"`.

Speed on one core: [4]phenylenes take 0.03 s, a 30-ring phenylene about 1 s, and a 100-ring phenylene about 10 s. For comparison, Rohn's method needs minutes already at 24 atoms.

## Tests

```bash
pip install "gapclosure[test,rdkit]" && pytest
```

The tests reproduce the published thresholds of all phenylenes up to N = 4 to 10⁻¹⁰, match the tree theorem against Rohn's method, and reproduce the acene–cyclobutadiene–acene table.

## References

- X. Yu, The distance to gap closure sets the triplet energy of phenylenes (code and data: doi:10.5281/zenodo.23448791). Its Supporting Information contains the one-minority theorem, the tree theorem, the closing patterns and the [N]phenylene values.
- J. Rohn, Systems of linear interval equations, Linear Algebra Appl. 126 (1989) 39–78.

MIT licence.
