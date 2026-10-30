# Signed Kekulé structures and Hückel gap-closure thresholds in phenylenes

Code and data for the paper of the same title.

> Author: Xiyao Yu (School of Chemistry, Chemical Engineering and Biotechnology, Nanyang Technological University, Singapore; ORCID [0009-0006-2673-5269](https://orcid.org/0009-0006-2673-5269); XYU029@e.ntu.edu.sg) · Archive: [doi:10.5281/zenodo.23448791](https://doi.org/10.5281/zenodo.23448791) · Code: MIT (`LICENSE`) · Data and figures: CC BY 4.0 (`LICENSE-DATA`)

The closure threshold ρ of a carbon skeleton is the smallest relative change of the bond couplings that closes the Hückel HOMO–LUMO gap; it is the regularity radius of an interval matrix and depends only on the signed Kekulé structures, which for phenylenes depend only on the L/A sequence. The repository computes ρ exactly for the [N]phenylenes with N ≤ 4, certifies each value, lists every coupling change that attains it, brackets ρ for N ≤ 6, compares the threshold with closure energies in a bond-length model, and tests the bond-length pattern predicted by that model with B3LYP.

## Layout

| Path | Contents |
|---|---|
| `src/gapclosure.py` | Molecular graphs, Kekulé structures and algebraic structure count, exact closure threshold (Poljak–Rohn formula and full vertex enumeration), Hückel + σ-spring bond-length model |
| `scripts/` | Numbered scripts that produce every number in the paper (see below) |
| `results/` | Outputs of `scripts/` (JSON tables and log files) |
| `dft_validation/` | B3LYP test: scripts, optimized geometries, all computed points, run logs |
| `figures/` | The figure files of the paper, the script that draws them, and `fig_data/` (every number plotted in Figs. 2, 5 and 7, as CSV) |

## Where each result comes from

Section, theorem, table and figure numbers refer to the paper.

| Paper item | Produced by | Output |
|---|---|---|
| Closed forms q² = 1/φ (benzocyclobutadiene) and q³ + q = 1 (biphenylene) of Theorem 2, compared with the computed ρ (Sections 3.1, 5.1) | `scripts/15_closed_forms.py` | `results/15_closed_forms.log` |
| Exact thresholds ρ (Table 1) | `scripts/01_thresholds.py` | `results/thresholds.json` |
| Certification of ρ at ρ(1 ± 10⁻⁸) by corner enumeration (Section 5.1) | `scripts/08_certify_thresholds.py` | `results/08_certify_thresholds.log` |
| All closing patterns (Table 1, Lemma 4, Section 5.3), including linear [3]phenylene | `scripts/07_closing_corners.py` | `results/07_closing_corners.log`, `results/closing_corners.json` |
| Closing patterns of sign-vector type (cross-check of the above) | `scripts/05_closing_patterns.py` | `results/05_closing_patterns.log` |
| Linear [3]phenylene contains the spectrum of benzocyclobutadiene; equal gap and threshold (Section 5.3) | `scripts/15_closed_forms.py` | `results/15_closed_forms.log` |
| Zigzag and helical [4]phenylene: bond correspondence (Theorem 4, Appendix B, Fig. 6) | `scripts/06_kekule_correspondence.py`, `scripts/09_zigzag_helical_correspondence.py` | `results/06_kekule_correspondence.log`, `results/09_zigzag_helical_correspondence.log`, `results/zigzag_helical_correspondence.json` |
| Two-sided bounds of Proposition 3 for the skeletons of Table 1 | `scripts/10_bounds.py` | `results/10_bounds.log`, `results/bounds.json` |
| Bounds by L/A sequence for N ≤ 6 (Table 2) | `scripts/11_bounds_by_LA_sequence.py` | `results/11_bounds_by_LA_sequence.log`, `results/bounds_by_LA_sequence.json` |
| Theorem 4 checked for all phenylenes with N ≤ 6 (Section 4) | `scripts/12_LA_sequence_equivalence.py` | `results/12_LA_sequence_equivalence.log` |
| Gluing formula (Lemma 2), numerical check | `scripts/13_gluing_formula_check.py` | `results/13_gluing_formula_check.log` |
| Linear phenylene polymer gapless; long zigzag and helical chains gapped (Section 5.2) | `scripts/14_polymer_bands.py` | `results/14_polymer_bands.log` |
| Closure energies in the bond-length model, f_σ = 21, 43, 60 eV Å⁻² (Section 6, Appendix A) | `scripts/02_bond_model.py` | `results/bond_model.json` |
| Columns of Tables 1 and 3; Hessian eigenvalue range (Section 6) | `scripts/03_tables.py` | `results/table1.json` (Table 1), `results/table2.json` (Table 3) |
| Spearman correlations and p values (Section 6, Table 3); 1.8 Å, 0.12 Å and 113 eV (Section 7); isotropic ceiling 5.3 and 6.7 (Appendix C) | `scripts/04_statistics.py` | `results/04_statistics.log` |
| B3LYP efficiencies (Appendix C, Fig. 7) | `dft_validation/` (see below), evaluated by `dft_validation/evaluate.py` | `dft_validation/*_points.json` |
| B3LYP gap gradient and its angle to the predicted pattern, cos_F = 0.989 (Section 6, Appendix C) | `dft_validation/gap_gradient.py`, `dft_validation/gap_gradient_angle.py` | `dft_validation/*_gapgrad.json`, `dft_validation/gapgrad_analysis.json` |
| Figures 1–7 | `figures/export_figdata.py` → `figures/fig_data/*.csv` → `figures/plot_figs.py` | see the table below |

| Figure | File(s) in `figures/` |
|---|---|
| 1 | `fig_bcb_graph.pdf`, `fig_bcb_kekule.pdf` |
| 2 | `fig1_hook.pdf` |
| 3 | `fig_kekule.pdf` |
| 4 | `fig_la.pdf` |
| 5 | `fig2_closing.pdf` |
| 6 | `fig_corresp.pdf` |
| 7 | `fig_dft.pdf` |

The file names are those used in the manuscript source.

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01_thresholds.py          # ~9 min on 2 cores
python scripts/02_bond_model.py          # ~15 min
python scripts/03_tables.py              # ~1 min
python scripts/04_statistics.py          # <1 min
python scripts/05_closing_patterns.py    # ~7 min
python scripts/06_kekule_correspondence.py
python scripts/07_closing_corners.py     # ~1.5 min
python scripts/08_certify_thresholds.py  # ~3 min
python scripts/09_zigzag_helical_correspondence.py  # ~4.5 min
python scripts/10_bounds.py
python scripts/11_bounds_by_LA_sequence.py
python scripts/12_LA_sequence_equivalence.py
python scripts/13_gluing_formula_check.py
python scripts/14_polymer_bands.py
python scripts/15_closed_forms.py        # <1 s
python figures/export_figdata.py         # ~7 min; writes figures/fig_data/*.csv
python figures/plot_figs.py              # <1 min; needs LaTeX with cm-super; writes figures/*.pdf
```

The B3LYP test (B3LYP/def2-SVP, density fitting, PySCF 2.14) is run from inside `dft_validation/`:

```bash
cd dft_validation
python opt.py biphenylene '[]'           # geometry optimization, ~5 min
python opt.py angular3 '["L"]'           # ~25 min
python directions.py                     # model directions (writes directions.json)
NT=2 python run_points.py biphenylene '[]'   # ~1.5 h
NT=2 python run_points.py angular3 '["L"]'   # ~5 h
python evaluate.py                       # efficiencies of predicted, bond-equalizing and random patterns
NT=2 python gap_gradient.py biphenylene '[]'      # 28 single points, ~15 min
NT=2 python gap_gradient.py angular3 '["L"]'     # 44 single points, ~60 min
python gap_gradient_angle.py             # cosine between x* and the B3LYP gap gradient
```

`run_points.py` saves after every geometry and skips points already present, so it can be stopped and restarted. The files in `dft_validation/logs/` are the console output of these runs.

## Citation

See `CITATION.cff`. A journal reference will be added once the paper is published.
