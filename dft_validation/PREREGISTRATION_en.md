# Preregistration: chemical validation of the predicted bond-length direction (operationalization of Section 5.2 of v0.3.1)

> **English translation.** The authoritative version is the Chinese original, `PREREGISTRATION.md`, whose successive SHA-256 hashes are listed in `PREREGISTRATION.sha256`. This translation was added after the study was completed and is not part of that hash chain. It follows the original sentence by sentence; notes in [square brackets] are the translator's and are not in the original.

Date of registration: 2026-10-18. This file was written before any DFT calculation on a displaced geometry; at this point only the equilibrium geometry optimizations of the two molecules had been run. Any later change to this file is recorded, with its reason, in the "Change log" at the end.

## Prediction to be tested

In the Hückel + σ-spring model (t0 = 2.5 eV, α = 4.1 eV/Å, K = 43 eV/Å² [K is the σ force constant, written k_σ in the paper]):

- biphenylene has a larger gap than angular [3]phenylene, but a lower closure cost;
- the mechanism is g0²/P: the frontier orbitals of biphenylene are more concentrated on bonds (P ≈ 0.45 against 0.27).

What is tested: whether the "cooperative bond-length direction that lowers the gap at least energy cost" given by this model is still more effective than other directions of equal cost in the real molecules as described by B3LYP/def2-SVP; and whether biphenylene is indeed pushed toward open-shell character more easily than angular [3]phenylene.

## Methods (fixed)

- **Electronic structure:** RKS B3LYP/def2-SVP, density fitting (PySCF default auxiliary basis). Geometry optimization with pyberny.
- **Geometry:** each molecule optimized to its equilibrium structure, used as the reference.
- **Directions (defined in the space of C–C bond lengths of the model; bonds numbered as in the graph construction):**
  - **P (predicted):** at the relaxed geometry of the model, x* ∝ −F⁻¹a, where a is the gradient of the model gap with respect to bond lengths and F is the Hessian of the model energy.
  - **C_BLA (totally symmetric control):** −δ_relax, i.e. the direction that undoes the equilibrium bond-length alternation of the model (making bond lengths more uniform).
  - **C_rand (random control):** 8 Gaussian random directions (numpy seed 2026), each computed with both signs.
- **Amplitude:** each direction normalized by the model elastic energy ½xᵀFx. P and C_BLA at E_m = 0.5, 1.0, 1.5, 2.0, 3.0 eV; C_rand only at E_m = 1.5 eV.
- **Mapping to Cartesian displacements:** only carbon atoms are moved (in plane), solved by minimum-norm least squares with the Wilson B matrix of C–C bond stretches; each hydrogen is translated with its carbon (keeping the C–H bond length). Reported: the linearization residual, and the cosine between the measured C–C bond-length changes after displacement and the target direction.
- **Observables:**
  - All geometries: RKS total energy (relative to the equilibrium geometry, denoted ΔE_DFT) and KS HOMO–LUMO gap (denoted Δgap = gap_ref − gap).
  - Every geometry on the P and C_BLA ladders, and the reference: natural-orbital occupation numbers from CASSCF(2,2)/def2-SVP (starting from RHF orbitals), with y0 taken as the occupation of the weakly occupied natural orbital; ⟨S²⟩ and E_BS − E_RKS of the broken-symmetry UB3LYP singlet (HOMO/LUMO mixing initial guess); UB3LYP triplet energy, giving ΔE_ST = E_T − E_S (the lower of RKS and BS).
- All geometries are single-point calculations at **fixed displacement, without further relaxation**.

## Pass criteria (locked)

**Criterion 1 (directional efficiency):** at E_m = 1.5 eV, define the efficiency η = Δgap / ΔE_DFT. For each random direction take the larger η of its two signs ("best-sign random control"), giving 8 values; their median is η_rand.
**Pass condition:** η_P ≥ 1.5 × η_rand in both molecules.
η_P / η_BLA is recorded as supplementary information and does not enter the decision.

**Criterion 2 (open-shell tendency):** along each molecule's predicted direction P, express y0 as a function of ΔE_DFT (linear interpolation between ladder points, including the reference point ΔE = 0).
- 2a: on a grid with spacing 0.25 eV (including 0) over the common ΔE_DFT range of the two curves, require y0(biphenylene) > y0(angular [3]) at every grid point.
- 2b: the lowest ΔE_DFT at which ⟨S²⟩ of BS-UB3LYP first exceeds 0.1. If it occurs in both molecules, biphenylene must be lower; if it occurs in only one, that one must be biphenylene; if it occurs in neither, 2b is recorded as "inconclusive" and does not overturn 2a.
**Pass condition:** 2a passes and 2b does not contradict it.
The increase of y0 relative to the reference point is also reported as supplementary information.

## Conclusion levels (locked)

- Criteria 1 and 2 both pass: the chemical story can be written (mechanism + reversal).
- Criterion 1 passes, criterion 2 fails: write only the mechanism, not the reversal.
- Criterion 1 fails: the conclusion stops at the model level; archive.

## Known limitations (declared in advance)

- The single-reference B3LYP description is unreliable for strongly open-shell systems; CASSCF(2,2) contains only static correlation; both are used only to compare directions and orderings.
- Displacements are fixed without relaxation, so ΔE_DFT contains costs such as angle strain that the model lacks; this is one of the things being tested.
- Only one pair of molecules; the conclusions are not extrapolated to other skeletons.

## Change log

(None yet) [this placeholder line was left in the original when the entries below were appended]

- 2026-10-18 (before any result on a displaced geometry): a CASSCF(2,2) single point at the biphenylene equilibrium geometry takes about 8 minutes. To limit the total time, along each ladder the CASSCF orbitals of the previous point are used as the initial guess (warm start); the first point uses the orbitals of the reference geometry. This is only a computational detail; the definitions of the observables and the criteria are unchanged. SHA-256 of the original file: 7eefd93c9cef475122b38cbecf198c8c8af748190897ef6bfb27e3fbeb5b2d4d.
- 2026-10-18 (still no result on any displaced geometry; only the CASSCF(2,2) y0 = 0.058 of the biphenylene reference geometry has been obtained): even with the warm start a single ladder point still takes more than 20 minutes, and angular [3] in full would take more than 10 hours. Instead, y0 is computed from broken-symmetry UHF/def2-SVP natural orbitals with the Yamaguchi formula (y = 1 − 2T/(1+T²), T = (n_HONO − n_LUNO)/2, i.e. the PUHF-type diradical index commonly used by Nakano et al.), replacing the CASSCF(2,2) y0. The form and thresholds of criterion 2 are unchanged; only the source of y0 changes. The CASSCF result for the biphenylene reference geometry is kept in biphenylene_points_casref.json as supplementary information. [The file is now in `superseded/`.]
- 2026-10-18 (data seen so far: RKS/BS/UHF results for the biphenylene reference geometry and the P ladder at 0.5–2.0 eV, and the results for the angular [3] reference geometry): UHF for the angular [3] reference geometry converged to an almost RHF-like solution (⟨S²⟩ = 0.013, y0 ≈ 0), whereas UHF for biphenylene shows strong spin-symmetry breaking (⟨S²⟩ = 1.52). This shows that the HOMO/LUMO mixing guess does not necessarily find the lowest UHF solution, so this difference in y0 risks being an artifact of convergence. Correction: both UHF and broken-symmetry UKS undergo internal stability analysis and are reconverged along unstable directions until stable (at most 4 rounds). All ladder points of both molecules are recomputed with the new procedure, and the old results are discarded. This is a correction to how the observables are implemented; the criteria are unchanged. But it was made after part of the data had been seen, and is recorded as such.
- 2026-10-18 (no new result since the previous correction): stability analysis of B3LYP UKS is too expensive (a single point had not finished after 20 minutes) and is not feasible. Instead: UHF still undergoes stability analysis; broken-symmetry UKS is converged separately from the HOMO/LUMO mixing guess and from the density of the stable UHF solution, keeping the lower energy. The criteria are unchanged.
- 2026-10-18 (under the new procedure only the biphenylene reference geometry has been obtained, with the same result as before the correction: y0 = 0.0857): to limit the total time, the C_BLA ladder is computed with RKS only (energy and gap). Open-shell quantities on the BLA ladder do not enter any criterion and were only supplementary information.
- 2026-10-18 (obtained under the new procedure: the biphenylene reference and P@0.5, with the same results as before; UHF stability analysis on the angular [3] reference geometry had not finished after 1 hour on a single thread): UHF stability analysis is not feasible for angular [3]. Instead, UHF is converged separately from two symmetry-broken initial guesses, keeping the lower energy: (i) HOMO/LUMO mixing; (ii) an antiferromagnetic guess, i.e. opposite spin densities placed on the two carbon sublattices (the usual way to find symmetry-broken UHF solutions in alternant hydrocarbons). BS-UKS is still started from the HOMO/LUMO mixing guess and from the selected UHF solution, keeping the lower. All P-ladder points of both molecules are recomputed uniformly with this procedure. The criteria are unchanged.
