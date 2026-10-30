# Plotted data (Figs. 2, 5 and 7)

Every number plotted in Figs. 2, 5 and 7 of the paper, as CSV. The files are written by `../export_figdata.py` and read by `../plot_figs.py`, which does no further calculation for these figures. The B3LYP values (`dft_*.csv`) are taken from `dft_validation/*_points.json`.

Units: molecular coordinates are schematic, with unit bond length (not optimized geometries). Energies are in eV; the Hückel gap in `hook_curves.csv` is in units of |β|.

| File | Columns | Contents |
|---|---|---|
| `hook_atoms.csv` | atom, x, y | Biphenylene, atom coordinates (Fig. 2) |
| `hook_bonds.csv` | atom1, atom2, cooperative_sign, single_bond_route | Each bond; cooperative_sign = +1 (strengthened) or −1 (weakened) at the threshold; single_bond_route = 1 marks the bond increased in the single-bond route |
| `hook_curves.csv` | route, largest_relative_change_percent, huckel_gap_abs_beta | Hückel gap along the single-bond and the cooperative route |
| `hook_summary.csv` | quantity, value | ρ of biphenylene; relative increase needed in the single-bond route |
| `closing_<skeleton>_atoms.csv` | atom, x, y, closing_orbital_weight | Weight of the closing orbital on each atom (sums to 1) (Fig. 5) |
| `closing_<skeleton>_bonds.csv` | atom1, atom2, coupling_change | Relative change of each coupling in the closing pattern (±ρ) |
| `closing_summary.csv` | skeleton, rho, sigma_min_at_pattern | Exact threshold; smallest singular value at the closing pattern (zero to machine precision) |
| `dft_ladders.csv` | molecule, pattern, model_energy_eV, b3lyp_deformation_energy_eV, b3lyp_gap_eV, gap_decrease_eV, y0, S2_broken_symmetry | Predicted pattern (P) and bond-equalizing control (BLA) at increasing model energy; the reference point appears once in each ladder; y0 and ⟨S²⟩ are given on the P ladder only (Fig. 7, Appendix C) |
| `dft_random.csv` | molecule, label, b3lyp_deformation_energy_eV, b3lyp_gap_eV, gap_decrease_eV | The 16 random patterns (8 directions × 2 signs), all at a model energy of 1.5 eV |
