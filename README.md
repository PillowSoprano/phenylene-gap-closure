# Many small bond changes close the HOMO–LUMO gap of phenylenes long before one large change does

Code, data and manuscript source for the paper of the same title.

> Author: Xiyao Yu (School of Chemistry, Chemical Engineering and Biotechnology, Nanyang Technological University, Singapore; ORCID [0009-0006-2673-5269](https://orcid.org/0009-0006-2673-5269); XYU029@e.ntu.edu.sg) · Status: manuscript in preparation · Code: MIT (`LICENSE`) · Data, figures and text: CC BY 4.0 (`LICENSE-DATA`)

Closing the Hückel gap of biphenylene through a single bond needs that bond to shorten by about 1.8 Å; letting all fourteen bonds change, some up and some down, needs at most 0.12 Å per bond. The repository computes how large such cooperative changes must be (the closure threshold ρ, exact for [N]phenylenes with N ≤ 4), where they must sit, what they cost in a bond-length model, and tests the predicted bond pattern in a preregistered B3LYP calculation.

## Layout

| Path | Contents |
|---|---|
| `src/gapclosure.py` | Molecular graphs, Kekulé structures and algebraic structure count, exact closure threshold (Rohn criterion and full vertex enumeration), Hückel + σ-spring bond-length model |
| `scripts/` | Numbered scripts that produce every number in the paper (see below) |
| `results/` | Outputs of `scripts/` (JSON tables and log files) |
| `dft_validation/` | Preregistered B3LYP test: protocol (`PREREGISTRATION.md` with SHA-256 change log), scripts, optimized geometries, all computed points |
| `dft_validation/superseded/` | Results replaced during the work, kept because the change log refers to them |
| `notebooks/` | Exploration notebooks v0.2 (thresholds, phenylene series, bond-length model) and v0.3.2 (local closure estimate, B3LYP evaluation); both executed |
| `paper/` | `paper.tex`, `refs.bib`, compiled `paper.pdf`, figures, `fig_data/` (every plotted number as CSV) and the two figure scripts |

## Where each result comes from

| Paper item | Produced by | Output |
|---|---|---|
| Closed forms q² = 1/φ, q³ + q = 1; Rohn criterion vs. vertex enumeration | `notebooks/bond_gap_research_seed_v0.2.ipynb` (§2, §5) | notebook output |
| Exact thresholds ρ (Table 1) | `scripts/01_thresholds.py` | `results/thresholds.json` |
| Closure energies in the bond-length model, K = 21, 43, 60 eV Å⁻² | `scripts/02_bond_model.py` | `results/bond_model.json` |
| Tables 1 and 2, Hessian eigenvalue range | `scripts/03_tables.py` | `results/table1.json`, `results/table2.json` |
| Spearman correlations and p values; 1.8 Å, 0.12 Å, 113 eV; model gap crossing at 0.66 eV; isotropic ceiling 5.3/6.7; model ratios 6.7/8.7 | `scripts/04_review_checks.py` | `results/04_review_checks.log` |
| Number of closing patterns (bent, triangular, zigzag, helical [4]) | `scripts/05_closing_patterns.py` | `results/05_closing_patterns.log` |
| Angular [3]phenylene: two mirror closing patterns | `notebooks/bond_gap_v0.3_baseline_audit.ipynb` (§1c) | notebook output |
| Helical and zigzag [4]phenylene share their signed Kekulé structures | `scripts/06_kekule_correspondence.py` | `results/06_kekule_correspondence.log` |
| Linear [3]phenylene contains the spectrum of benzocyclobutadiene | `notebooks/bond_gap_research_seed_v0.2.ipynb` (§6) | notebook output |
| B3LYP test (Section 6, criteria 1 and 2) | `dft_validation/` (see below), evaluated by `dft_validation/evaluate.py` and notebook v0.3.2 §6 | `dft_validation/*_points.json` |
| Figures 1–3 | `paper/export_figdata.py` → `paper/fig_data/*.csv` → `paper/plot_figs.py` | `paper/fig*.pdf` |

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01_thresholds.py          # ~9 min on 2 cores
python scripts/02_bond_model.py          # ~15 min
python scripts/03_tables.py              # ~1 min
python scripts/04_review_checks.py       # <1 min
python scripts/05_closing_patterns.py    # ~7 min
python scripts/06_kekule_correspondence.py
python paper/export_figdata.py           # ~5 min; then
python paper/plot_figs.py                # needs LaTeX with cm-super
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
python evaluate.py                       # applies the preregistered criteria
```

`run_points.py` saves after every geometry and skips points already present, so it can be stopped and restarted. The protocol and pass criteria were written to `PREREGISTRATION.md` before any displaced geometry was computed; every later change is appended there with the SHA-256 of the file after the change (`PREREGISTRATION.sha256`). Appendix B of the paper lists the four changes.

## Citation

See `CITATION.cff`. A journal reference will be added once the paper is published.
