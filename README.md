# Signed Kekulé structures and Hückel gap-closure thresholds in phenylenes

Code, data and manuscript source for the paper of the same title.

> Author: Xiyao Yu (School of Chemistry, Chemical Engineering and Biotechnology, Nanyang Technological University, Singapore; ORCID [0009-0006-2673-5269](https://orcid.org/0009-0006-2673-5269); XYU029@e.ntu.edu.sg) · Status: manuscript in preparation · Code: MIT (`LICENSE`) · Data, figures and text: CC BY 4.0 (`LICENSE-DATA`)

The closure threshold ρ of a carbon skeleton is the smallest relative change of the bond couplings that closes the Hückel HOMO–LUMO gap; it is the regularity radius of an interval matrix and depends only on the signed Kekulé structures, which for phenylenes depend only on the L/A sequence. The repository computes ρ exactly for the [N]phenylenes with N ≤ 4, certifies each value, lists every coupling change that attains it, brackets ρ for N ≤ 6, compares the threshold with closure energies in a bond-length model, and tests the bond-length pattern predicted by that model with B3LYP.

## Layout

| Path | Contents |
|---|---|
| `src/gapclosure.py` | Molecular graphs, Kekulé structures and algebraic structure count, exact closure threshold (Poljak–Rohn formula and full vertex enumeration), Hückel + σ-spring bond-length model |
| `scripts/` | Numbered scripts that produce every number in the paper (see below) |
| `results/` | Outputs of `scripts/` (JSON tables and log files) |
| `dft_validation/` | B3LYP test: scripts, optimized geometries, all computed points |
| `notebooks/` | Exploration notebooks v0.2 (thresholds, phenylene series, bond-length model) and v0.3 (local closure estimate); both executed |
| `paper/` | `paper.tex`, `refs.bib`, compiled `paper.pdf`, figures, `fig_data/` (every plotted number as CSV) and the figure scripts |

## Where each result comes from

| Paper item | Produced by | Output |
|---|---|---|
| Closed forms q² = 1/φ, q³ + q = 1 | `notebooks/bond_gap_research_seed_v0.2.ipynb` (§2, §5); compared with the computed ρ in `scripts/08_certify_thresholds.py` | notebook output |
| Exact thresholds ρ (Table 1) | `scripts/01_thresholds.py` | `results/thresholds.json` |
| Certification of ρ at ρ(1 ± 10⁻⁸) by corner enumeration | `scripts/08_certify_thresholds.py` | `results/08_certify_thresholds.log` |
| All closing patterns (Table 1, Section 4.3), including linear [3]phenylene | `scripts/07_closing_corners.py` | `results/07_closing_corners.log`, `results/closing_corners.json` |
| Closing patterns of sign-vector type (cross-check) | `scripts/05_closing_patterns.py` | `results/05_closing_patterns.log` |
| Helical and zigzag [4]phenylene: bond correspondence of Proposition 1 (Appendix B) | `scripts/06_kekule_correspondence.py`, `scripts/09_zigzag_helical_correspondence.py` | `results/06_kekule_correspondence.log`, `results/09_zigzag_helical_correspondence.log`, `results/zigzag_helical_correspondence.json` |
| Two-sided bounds of Proposition 2 for the skeletons of Table 1 | `scripts/10_bounds.py` | `results/10_bounds.log`, `results/bounds.json` |
| Bounds by L/A sequence for N ≤ 6 (Table 2) | `scripts/11_bounds_by_LA_sequence.py` | `results/11_bounds_by_LA_sequence.log`, `results/bounds_by_LA_sequence.json` |
| Theorem 2 checked for all phenylenes with N ≤ 6 | `scripts/12_LA_sequence_equivalence.py` | `results/12_LA_sequence_equivalence.log` |
| Gluing formula (Lemma 2), numerical check | `scripts/13_gluing_formula_check.py` | `results/13_gluing_formula_check.log` |
| Linear phenylene polymer gapless; long zigzag and helical chains gapped (Section 5.2) | `scripts/14_polymer_bands.py` | `results/14_polymer_bands.log` |
| Closure energies in the bond-length model, f_σ = 21, 43, 60 eV Å⁻² | `scripts/02_bond_model.py` | `results/bond_model.json` |
| Tables 1 and 2, Hessian eigenvalue range | `scripts/03_tables.py` | `results/table1.json`, `results/table2.json` |
| Spearman correlations and p values; 1.8 Å, 0.12 Å, 113 eV; isotropic baseline 5.3/6.7; model ratios 6.7/8.7 | `scripts/04_review_checks.py` | `results/04_review_checks.log` |
| Linear [3]phenylene contains the spectrum of benzocyclobutadiene | `notebooks/bond_gap_research_seed_v0.2.ipynb` (§6) | notebook output |
| B3LYP efficiencies (Section 5.2) | `dft_validation/` (see below), evaluated by `dft_validation/evaluate.py` | `dft_validation/*_points.json` |
| B3LYP gap gradient and its angle to the predicted pattern (cos_F = 0.989) | `dft_validation/gap_gradient.py`, `dft_validation/gap_gradient_angle.py` | `dft_validation/*_gapgrad.json`, `dft_validation/gapgrad_analysis.json` |
| Figures 1 and 3 | `paper/export_figdata.py` → `paper/fig_data/*.csv` → `paper/plot_figs.py` | `paper/fig1_hook.pdf`, `paper/fig2_closing.pdf` |
| Figures 2, 4 and 5 | `paper/plot_new_figs.py` (after `scripts/09_…`) | `paper/fig_kekule.pdf`, `paper/fig3_dft_2panel.pdf`, `paper/fig_corresp.pdf` |

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01_thresholds.py          # ~9 min on 2 cores
python scripts/02_bond_model.py          # ~15 min
python scripts/03_tables.py              # ~1 min
python scripts/04_review_checks.py       # <1 min
python scripts/05_closing_patterns.py    # ~7 min
python scripts/06_kekule_correspondence.py
python scripts/07_closing_corners.py     # ~1.5 min
python scripts/08_certify_thresholds.py  # ~3 min
python scripts/09_zigzag_helical_correspondence.py
python scripts/10_bounds.py
python scripts/11_bounds_by_LA_sequence.py
python scripts/12_LA_sequence_equivalence.py
python scripts/13_gluing_formula_check.py
python scripts/14_polymer_bands.py
python paper/export_figdata.py           # ~5 min; then
python paper/plot_figs.py                # needs LaTeX with cm-super
python paper/plot_new_figs.py
cd paper && latexmk -pdf paper.tex
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

`run_points.py` saves after every geometry and skips points already present, so it can be stopped and restarted.

## Citation

See `CITATION.cff`. A journal reference will be added once the paper is published.
