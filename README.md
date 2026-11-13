# The distance to gap closure sets the triplet energy of phenylenes: code and data

Code and data for the paper *The distance to gap closure sets the triplet energy of phenylenes* (X. Yu). Every number, table and figure in the paper and its Supporting Information (SI) comes from a script in this repository.

> Author: Xiyao Yu (School of Chemistry, Chemical Engineering and Biotechnology, Nanyang Technological University, Singapore; ORCID [0009-0006-2673-5269](https://orcid.org/0009-0006-2673-5269)) · Archive: [doi:10.5281/zenodo.23448791](https://doi.org/10.5281/zenodo.23448791) · Code: MIT (`LICENSE`) · Data and figures: CC BY 4.0 (`LICENSE-DATA`)

The closure threshold ρ of a carbon skeleton is the smallest relative change of the bond couplings that closes the Hückel HOMO–LUMO gap. The paper shows that ρ ranks the triplet energies of phenylenes better than the gap and the usual gap-based predictors, and explains why: the triplet relaxes toward the nearest geometry with a closed gap, and its energy is a nearly fixed fraction of the cost of getting there.

## Layout

| Folder | Contents |
|---|---|
| `gapclosure/` | Python package: ρ, closing patterns and closing orbitals of any conjugated hydrocarbon; exact algorithm for phenylenes whose rings form a tree; drawings |
| `dataset_tree_phenylenes/` | ρ, closing pattern and Kekulé counts of every tree-like phenylene with N ≤ 12 (9 232 L/A classes, 147 673 molecules); see its README |
| `kekule_thresholds/` | Signed Kekulé structures, exact thresholds and closing patterns of the small phenylenes (SI Sections S1–S4, Appendices A and B), and the B3LYP steepest-descent check (Fig. S3) |
| `beyond_huckel/` | Hückel–SSH and PPP triplet energies of all 467 classes with N ≤ 9; nearest closed-gap geometry; UB3LYP triplets of small phenylenes |
| `analysis/` | Table 1, lowest PPP triplets, classes and molecules within classes, symmetry breaking, DFT of the isomers and of molecules A and B |
| `analysis/geometry/` | Thresholds and closing-pattern agreement at given geometries (Sections 3, 4 and 8) |
| `closure_energy/` | Lower bound on the triplet energy, acene pairs, cheapest closed-gap geometry (Section 5, SI Section S6) |
| `figures/` | Figure files of the paper. Figs. 1, 2, 4b–c and 5 are TikZ/pgfplots files with the plotted numbers written in them |

## Where each result comes from

Run every script from its own folder (`cd analysis; python census_final.py`); paths between folders are relative.

### Main text

| Item | Script | Output or data |
|---|---|---|
| Fig. 1: lowest PPP triplet of the 281 nine-ring classes against ρ and against the gap; molecules A and B | `analysis/fig1_hook_data.py` | prints the coordinates in `figures/hook_rho.tex`, `figures/hook_gap.tex` |
| Fig. 1, Section 6: B3LYP//xTB triplets of A and B, 1.62 and 1.44 eV | `analysis/dft_hook_ab.py` | `analysis/dft_hook_ab_A.json`, `dft_hook_ab_B.json` |
| Fig. 2: mechanism sketch | drawn from the one-dimensional model of Section 5 | `figures/scheme.tex` |
| Fig. 3: graph and signed Kekulé structures of benzocyclobutadiene | `kekule_thresholds/figures/plot_figs.py` | `figures/fig_bcb_graph.pdf`, `fig_bcb_kekule.pdf` |
| Fig. 4a: Kekulé structures of biphenylene | `kekule_thresholds/figures/plot_figs.py` | `figures/fig_kekule.pdf` |
| Fig. 4b–c, Section 3: biphenylene, ρ = 0.310 at the relaxed ground state and 0.111 at the relaxed triplet; weights c_a | `analysis/geometry/mech.py`; bond lengths and triplet relaxation in `beyond_huckel/biphenylene.py` | printed; `beyond_huckel/biphenylene_output.txt` |
| Section 3: ρ = 0.310 as bond-length changes (≤ 0.09 Å), bond lengths 1.388–1.488 Å | `analysis/biphenylene_bond_scale.py` | printed |
| Sections 3 and 4: bridging bond 1.508 Å in UB3LYP; ρ = 0.350 at the UB3LYP ground state | `analysis/dft_isomers.py`, `analysis/geometry/dft_rho_geom.py biphenylene` | `analysis/dft_isomers.json` |
| Section 5, Fig. 5: E_T1/U between 0.73 and 0.87 (mean 0.82); r between 0.45 and 0.62; Eq. (9) within 11 %, correlation 0.985; Spearman 0.983 (ρ) and 0.980 (U) | `beyond_huckel/l2nearest.py`, then `analysis/section5_ratio.py` | `beyond_huckel/l2nearest.npy` |
| Section 5: the cheapest closed-gap geometry lies 1–2 % below U | `closure_energy/cheapest_closure.py` | reads `beyond_huckel/l2nearest_points.json` |
| Section 5: two acenes joined through a four-membered ring, E_T1/U from 0.87 to 0.003, closure cost about 1.8 \|β\| | `closure_energy/acene_matching_closure.py` | `closure_energy/acene_matching_closure.json` |
| Table 1, Section 6: Spearman correlations of ρ and five other predictors; bootstrap intervals | `analysis/competitors.py`, then `analysis/census_final.py`; checks in `analysis/check_competitors.py` | `analysis/census_final.txt` |
| Section 6: lowest PPP triplet of every class | `analysis/ppp_bs.py`, `analysis/ppp_bs_analyse.py`, `analysis/s0_bs_check.py` | `analysis/ppp_bs_classes.jsonl`, `ppp_bs_summary.txt`, `s0_bs_check.json` |
| Section 6: the gap orders 11 % of cross-class pairs wrongly at N = 9 | `analysis/tiebreak.py` | printed |
| Section 6: where A and B close their gap (90 % of the closing orbital on a linear [3] unit, 89 % of the triplet change, 99 % sign agreement) | `analysis/hook_ab_closure.py ../gapclosure/src ../beyond_huckel/patterns_N2-9.json` | `analysis/hook_ab_closure.json` |
| Section 7: ρ depends only on the L/A class, checked on every molecule | `analysis/class_invariance.py Nmin Nmax` | printed |
| Section 7: zigzag and helical [4]phenylene, ρ = 0.2821 at relaxed geometries; PPP 2.11 and 2.05 eV; UB3LYP 1.94 and 1.81 eV | `analysis/relaxed_rho_AA.py`, `analysis/ppp_helical_zigzag.py`, `analysis/dft_isomers.py` | `analysis/relaxed_rho_AA.out`, `analysis/dft_isomers.json` |
| Section 8: symmetric and broken Hückel–SSH triplet of angular [3]phenylene (0.022 \|β\|, 86 % of bonds, 97 % of the change) | `analysis/angular3_symmetry_hssh.py` | printed |
| Section 8: symmetry of the UB3LYP triplets; agreement of their bond-length changes with the closing patterns | `analysis/geometry/dft_symmetry_check.py`, `analysis/geometry/dft_closing_agreement.py <molecule>` | printed |

### Supporting Information

| Item | Script | Output or data |
|---|---|---|
| Table S1: c_a of biphenylene at four sets of reference couplings | `analysis/geometry/mech.py`, `analysis/geometry/dft_rho_geom.py biphenylene` | printed |
| S2: closed forms of the one-minority theorem | `kekule_thresholds/scripts/15_closed_forms.py` | `kekule_thresholds/results/15_closed_forms.log` |
| S2: two-sided bounds | `kekule_thresholds/scripts/10_bounds.py` | `kekule_thresholds/results/bounds.json` |
| S3: L/A sequence determines the signed Kekulé polynomial; gluing formula (Fig. S1) | `kekule_thresholds/scripts/12_LA_sequence_equivalence.py`, `13_gluing_formula_check.py`, `figures/plot_figs.py` | `kekule_thresholds/results/12_*.log`, `13_*.log` |
| Table S2: exact thresholds, their certification, Kekulé counts, single-bond changes | `kekule_thresholds/scripts/01_thresholds.py`, `08_certify_thresholds.py`, `03_tables.py` | `kekule_thresholds/results/thresholds.json`, `table1.json`, `08_*.log` |
| S4.1: frontier-sign bound ρ_F | `kekule_thresholds/scripts/18_frontier_sign_bound.py` | `kekule_thresholds/results/frontier_sign_bound.json` |
| Table S3: bounds by L/A sequence for N ≤ 6 | `kekule_thresholds/scripts/11_bounds_by_LA_sequence.py` | `kekule_thresholds/results/bounds_by_LA_sequence.json` |
| S4.2: linear phenylene polymer gapless, zigzag and helical chains gapped | `kekule_thresholds/scripts/14_polymer_bands.py` | `kekule_thresholds/results/14_polymer_bands.log` |
| S4.3, Fig. S2: closing patterns and closing orbitals; terminal elimination | `kekule_thresholds/scripts/07_closing_corners.py` (cross-check `05_closing_patterns.py`), `17_terminal_elimination.py`; drawn by `figures/export_figdata.py` and `plot_figs.py` | `kekule_thresholds/results/closing_corners.json`, `terminal_elimination.json` |
| S5: exact thresholds of all tree-like phenylenes | `dataset_tree_phenylenes/build.py`, package `gapclosure` | `dataset_tree_phenylenes/` |
| S6: upper bound from the closing pattern (the triplet reaches 0.19 to 0.60 of it) | `beyond_huckel/bound.py` | `beyond_huckel/bound_output.txt` |
| S6: lower bound (0.006 against 0.85 for linear [64]phenylene) | `closure_energy/global_lower_bound.py`, `closure_energy/linear_scaling.py` | `closure_energy/global_lower_bound.json`, `linear_scaling.json` |
| S7: Hückel–SSH triplets of all classes | `beyond_huckel/t1_scan.py` | `beyond_huckel/t1_scan.json` |
| S7: perturbed starts lower the Hückel–SSH triplet in 12 of 467 classes, by at most 0.022 \|β\| | `analysis/hssh_bs_scan.py` | `analysis/hssh_bs_classes.json` |
| Table S5: Spearman correlation by N; Table S8: classes of equal ρ | `analysis/census_final.py` | `analysis/census_final.txt` |
| Table S6: PPP triplets in the symmetric convention | `beyond_huckel/ppp_scan.py`; `beyond_huckel/an_ppp2.py ppp_scan_N8.json ppp_base_N9-9.json ppp_rohf_ohno_N2-9.json ppp_rohf_mn_N2-9.json ppp_rohf_u8_N2-9.json` | `beyond_huckel/ppp_*.json` |
| Table S7: molecules within classes | `analysis/within_class.py` (or `within_class_chunk.py N start stop`), `analysis/within_class_analyse.py`, `analysis/tiebreak.py` | `analysis/within_class_mols.jsonl`, `within_class_summary.txt` |
| Table S9: UB3LYP triplets of seven synthesized phenylenes | `beyond_huckel/dft_t1.py`, `analysis/dft_isomers.py` | `beyond_huckel/dft_t1_*.json`, `analysis/dft_isomers.json` |
| Table S4: B3LYP//xTB protocol, molecules A and B | `analysis/dft_hook_ab.py` | `analysis/dft_hook_ab_*.json` |
| S7, Fig. S3: steepest descent of the B3LYP gap, cos_F = 0.989 | `kekule_thresholds/dft_validation/` (see below) | `kekule_thresholds/dft_validation/*_points.json`, `gapgrad_analysis.json` |
| S7: nearest closed-gap point | `beyond_huckel/l2nearest.py` | `beyond_huckel/l2nearest.npy` |
| Appendix A, Fig. S4: bond correspondence of zigzag and helical [4]phenylene | `kekule_thresholds/scripts/06_kekule_correspondence.py`, `09_zigzag_helical_correspondence.py` | `kekule_thresholds/results/zigzag_helical_correspondence.json` |
| Appendix B: linear phenylenes up to N = 160 | `kekule_thresholds/scripts/16_linear_chains.py` | `kekule_thresholds/results/linear_chains.json` |

`kekule_thresholds/scripts/02_bond_model.py` relaxes the bond-length model used by `17_terminal_elimination.py` and by the B3LYP check.

## Reproducing

```bash
pip install -r requirements.txt
pip install -e gapclosure
gunzip -k dataset_tree_phenylenes/*.gz        # class_invariance.py and within_class*.py read the molecule table
```

Thresholds and closing patterns of the small phenylenes (run times on two cores):

```bash
cd kekule_thresholds
python scripts/01_thresholds.py          # ~9 min
python scripts/02_bond_model.py          # ~15 min
python scripts/03_tables.py              # ~1 min
for s in 05 06 07 08 09 10 11 12 13 14 15 16 17 18; do python scripts/${s}_*.py; done
python figures/export_figdata.py         # ~7 min; writes figures/fig_data/*.csv
python figures/plot_figs.py              # <1 min; needs LaTeX with cm-super
```

Triplet energies and the census:

```bash
cd beyond_huckel
python t1_scan.py                                         # Hückel–SSH triplets, N <= 9
TAG=base NP=4 python ppp_scan.py 8 2                      # PPP-UHF, Ohno, N <= 8 (stored as ppp_scan_N8.json)
TAG=base NP=4 python ppp_scan.py 9 9                      # N = 9
T1MODE=ROHF TAG=rohf_ohno NP=4 python ppp_scan.py 9 2     # PPP-ROHF, Ohno; TAG=rohf_mn POT=mn and TAG=rohf_u8 U=8 give the other two
NSTART=30 python l2nearest.py                             # nearest closed-gap point, a few minutes
IDS="PT02-00001 PT04-00001 PT08-00001 PT08-00049" NSTART=30 python l2nearest.py   # points for cheapest_closure.py
cd ../analysis
python competitors.py                # writes competitors.json, read by the scripts below
python ppp_bs.py 2 9                 # lowest PPP triplets, several hours
python census_final.py; python ppp_bs_analyse.py; python tiebreak.py; python within_class_analyse.py
python section5_ratio.py; python fig1_hook_data.py
cd geometry; python mech.py; python dft_closing_agreement.py biphenylene
```

`within_class.py` covers 2 315 molecules and takes several hours; `within_class_chunk.py` runs it in slices and appends to `within_class_mols.jsonl`.

The DFT calculations need PySCF; `analysis/dft_isomers.py`, `analysis/dft_hook_ab.py` and `beyond_huckel/dft_t1.py` also need geomeTRIC, and `dft_hook_ab.py` needs the `xtb` program (GFN2-xTB). Optimizations of the four-ring phenylenes take hours on one core, and the density-fitting tensor of a nine-ring molecule is about 11 GB. All DFT results are included as JSON, so the analysis scripts run without them.

The B3LYP steepest-descent check (B3LYP/def2-SVP, density fitting, geometries by pyberny) runs from `kekule_thresholds/dft_validation/`:

```bash
python opt.py biphenylene '[]'           # geometry optimization, ~5 min
python opt.py angular3 '["L"]'           # ~25 min
python directions.py                     # model directions (writes directions.json)
NT=2 python run_points.py biphenylene '[]'   # ~1.5 h
NT=2 python run_points.py angular3 '["L"]'   # ~5 h
python evaluate.py
NT=2 python gap_gradient.py biphenylene '[]'      # 28 single points, ~15 min
NT=2 python gap_gradient.py angular3 '["L"]'     # 44 single points, ~60 min
python gap_gradient_angle.py             # cosine between x* and the B3LYP gap gradient
```

`run_points.py` saves after every geometry and skips points already present, so it can be stopped and restarted. The files in `dft_validation/logs/` are the console output of these runs.

The TikZ figure files in `figures/` need TikZ and pgfplots.

## Citation

See `CITATION.cff`.
