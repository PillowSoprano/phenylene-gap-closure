# Hückel gap-closure thresholds of all tree-like phenylenes, N ≤ 12

This dataset covers every phenylene with N ≤ 12 hexagons whose hexagons form a tree. Each pair of neighbouring hexagons is joined by a four-membered ring.

For each phenylene it gives the **closure threshold ρ**: the smallest relative change of the bond couplings that closes the Hückel HOMO–LUMO gap. It also gives the bracket that certifies ρ, a closing pattern, the closing orbital, and the counts of Kekulé structures.

## What a row is

The threshold depends only on the **L/A class** of a phenylene. This is the tree of hexagons in which every hexagon with two neighbours is marked:
- **L** (linear) when its two four-membered rings sit on opposite bonds;
- **A** (angular) when they sit on bonds two apart.

A hexagon with three neighbours always carries its four-membered rings on alternate bonds.

One class usually contains several molecules. For example, zigzag and helical [4]phenylene are not isomorphic as graphs, but both belong to class `AA` and have the same ρ. They differ in their Hückel gap.

The dataset therefore has three files:

| File | One row per | Rows |
|---|---|---|
| `tree_phenylenes_N1-12.csv` | L/A class | 9 232 |
| `tree_phenylene_molecules_N1-12.csv` | molecule, distinct up to graph isomorphism | 147 673 |
| `tree_phenylene_patterns_N1-12.jsonl` | L/A class: closing pattern and closing orbital on a representative molecule | 9 232 |

The number of classes for each N is 1, 1, 2, 4, 8, 18, 44, 108, 282, 761, 2096, 5907. The number of molecules is 1, 1, 2, 5, 12, 37, 123, 446, 1689, 6693, 27034, 111630. This is the number of catafusenes, OEIS A002216 (the Harary–Read numbers), as it should be: tree-like phenylenes correspond one to one to catacondensed hexagonal systems, helicenes included.

## Columns

### `tree_phenylenes_N1-12.csv`

| Column | Meaning |
|---|---|
| `id` | `PTnn-kkkkk`: N, then rank by decreasing ρ within N |
| `N` | number of hexagons |
| `la_tree` | canonical string of the L/A tree. Each hexagon is written `(label children...)`, the label is L, A, or empty (one or three neighbours), and the string is rooted at the centre of the tree |
| `word` | for chains, the L/A word of the N−2 inner hexagons in L/A notation (the smaller of its two reading directions); empty for branched trees and for N ≤ 2 |
| `n_L`, `n_A`, `n_branch` | hexagons marked L, marked A, and with three neighbours |
| `K` | number of Kekulé structures |
| `ASC` | algebraic structure count, \|det T\| |
| `minority` | number of minority Kekulé structures, (K − ASC)/2 |
| `rho` | closure threshold. This is the root of the closing pattern's determinant inside the certified bracket, so it is an upper bound on ρ and lies within 2⁻⁴⁰ of it |
| `rho_lower`, `rho_upper` | exact rationals with rho_lower < ρ ≤ rho_upper (tree theorem, 40 bisection steps) |
| `q` | (1 − ρ)/(1 + ρ) |
| `n_molecules` | molecules in the class |
| `gap_min`, `gap_max` | smallest and largest Hückel gap (units of \|β\|) among the molecules of the class, at equal couplings |
| `representative_smiles` | one molecule of the class, chosen so that its rings do not overlap in the plane whenever such a molecule exists |
| `seconds` | time taken by the tree theorem for this class (one core) |

For N = 1 (benzene) there is no minority structure, so the gap cannot close for ε < 1. The file records ρ = 1.

### `tree_phenylene_molecules_N1-12.csv`

| Column | Meaning |
|---|---|
| `id`, `la_tree`, `N` | the class of the molecule |
| `smiles` | canonical SMILES (RDKit) |
| `huckel_gap` | Hückel gap at equal couplings, units of \|β\| |
| `overlaps_in_plane` | true if the drawing with regular polygons has two atoms closer than half a bond length. These are helicene-like molecules; they exist, but not as flat drawings |
| `links` | the construction: a list of (i, p, j), where hexagon j is attached through a four-membered ring on bond p of hexagon i. Pass it to `gapclosure.phenylene_tree(links)` |

### `tree_phenylene_patterns_N1-12.jsonl`

Each line holds one class: `id`, `la_tree`, `N`, `rho`, `rho_lower`, `rho_upper`, and `representative`. The representative has:
- `links`;
- `smiles_mapped`: SMILES whose atom map number is the atom index plus one;
- `coords`: 2D coordinates, with bond length one;
- `closing_pattern`: a list of [u, v, s]. The bond u–v is raised by ρ when s = +1 and lowered by ρ when s = −1;
- `closing_orbital`: the weight of the zero mode on each atom at the closing pattern; the weights sum to one.

## How it was computed

1. **Classes.** Every tree with N vertices and maximum degree 3 (networkx `nonisomorphic_trees`), with every L/A labelling of its degree-2 vertices, deduplicated by a canonical string.
2. **Molecules.** Every planar placement of each class: the side of each angular turn and the order of the branches at each three-neighbour hexagon. These are deduplicated by RDKit canonical SMILES of the carbon graph.
3. **ρ.** The decision procedure of the tree theorem runs with exact integer arithmetic, and bisection on dyadic ε gives the bracket. The closing pattern is the corner returned at `rho_upper`, and `rho` is the root of the determinant along that corner. Software: `gapclosure` 0.1.0.
4. **K.** The same pass from the leaves to the root, with all signs positive.

## Checks

- The ρ values of all nine skeletons with N ≤ 4 agree with exact vertex enumeration (Rohn's criterion) to 10⁻¹⁰, or better.
- For all 22 chain classes with N ≤ 6, K, ASC and the range of Hückel gaps agree with the earlier phenylene paper. ρ also lies between that paper's lower bound and its upper bound.
- The Kekulé counts from the tree pass agree with direct enumeration on 47 molecules with N ≤ 6.
- The molecule counts reproduce OEIS A002216 for N ≤ 12.
- On 300 randomly chosen classes with N ≥ 3:
  - the closing pattern makes T singular at `rho`: the ratio σ_min/σ_max is at most 3·10⁻¹⁶;
  - none of 300 random corners per class at `rho_lower` has det T ≤ 0;
  - another molecule of the same class, run separately through the tree theorem, gives the same ρ in all 291 classes that have one.
- Only biphenylene has a single minority structure, as the one-minority theorem requires.

The build takes about 22 minutes on two cores (`python build.py 12`, log in `build_N12.log`).

## Licence and citation

Data: CC BY 4.0. Code: MIT.

Cite the repository (doi:10.5281/zenodo.23448791) and the paper; see `../CITATION.cff`.
