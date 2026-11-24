# Tier 2: density matrix of a small particle under a few-cycle pulse

Code: `src/enzkin/densmat.py`. Model: Un, Sarkar & Sivan, refs. [16] (formalism) and [62] (results). Validation: `scripts/validate_densmat.py` → `results/tier2/densmat_validation.json`. Units as everywhere in `enzkin`: eV, fs, nm; fields in V/nm; "P" means P/ε₀ in V/nm.

## 1. Equations of [16] in oscillator form

[16] Eqs. (13)–(15) evolve the populations f(E, t) and the dipole-weighted coherences X(E, E′) between conduction-band states of different energy, with the real field (no rotating-wave or cycle average):

```
∂f(E)/∂t |_photon = −(iV/2ħ) E_loc ∫ρ(E′) [X(E,E′) − X(E′,E)] dE′                                  (13)
∂X(E,E′)/∂t + η̄ [X(E,E′) − X(E′,E)] + i (E−E′)/ħ · X(E,E′) = −(i/ħ) |D|² E_loc [f(E) − f(E′)]      (14)
P = (V/4) ∫∫ ρ(E)ρ(E′) [X(E,E′) + X(E′,E)] dE dE′                                                    (15)
```

η̄ = [η(E) + η(E′)]/2, η = e–e + e–ph + e–impurity out-scattering rates of the current f. Because X(E′,E) = X(E,E′)*, write X = a + ib. Then ȧ = Δ b and ḃ = −2η̄ b − Δ a − |D|²E_loc(f − f′)/ħ with Δ = (E − E′)/ħ, so each pair of energies is a damped oscillator. With Q_p = (e/ε₀) ρ_lo ρ_hi V a_p dE² for the pair p = (E_lo < E_hi):

```
Q̈_p + 2 η̄_p Q̇_p + Δ_p² Q_p = S_p (f_lo − f_hi) E_loc,        P/ε₀ = Σ_p Q_p
S_p = (e²/ε₀) (Δ_p/ħ) W(E_lo,E_hi) ρ_lo ρ_hi dE²,   W = V|D|²/e²
d(ρ f dE)/dt |_hi = − d(ρ f dE)/dt |_lo = E_loc Q̇_p / [(e/ε₀)(E_hi − E_lo)]
```

The last line is (13): the work E_loc Q̇_p done on a pair moves electrons from the lower to the upper level, or back when it is negative (stimulated emission). Number is conserved pair by pair and the electron energy changes by exactly the work E_loc·Ṗ; the code holds both to round-off (`tests/test_densmat.py`). The sphere's field is [16] Eq. (17): E_loc = [3ε_h E_inc − P/ε₀]/(ε∞ + 2ε_h).

Collisions on the populations are the WP1 terms (`ElectronElectronB1`, `ElectronPhonon`), refreshed together with η every 0.1 fs; the lattice stays at 300 K during the tens of femtoseconds simulated. As in [16], the field-driven coupling between different coherences is dropped.

## 2. Departures from the printed formulas

**Dipole strength, [16] Eq. (16): unusable as printed, reconstructed.** Evaluated as printed, Eq. (16) is negative over most of the (E, E′) plane and misses the oscillator-strength sum rule by a factor of order −10⁵ at E_F for L = 4 nm (`scripts/validate_densmat.py`, entry `sum_rule_at_EF`). The figure of the same quantity ([16] Fig. 1b) is sound, so the formula was rederived:

- States of a box of width L along the field: plane waves across, standing waves along z. A dipole transition keeps the transverse momentum, so k_z′² = k′² − k² + k_z². The squared box element is 64 k_z²k_z′²/[L²(k_z² − k_z′²)⁴] for odd level differences and zero for even ones.
- Averaging over the direction of k and summing over final states in dE′ gives

  W(k, k′) = (64π / L k k′) ∫ dk_z k_z² k_z′ / [(k_z − k_z′)⁴ (k_z + k_z′)⁴],  k_z from √max(k² − k′², 0) to k.

  Away from the diagonal this reproduces [16] Fig. 1b: log₁₀|D|² = −58.64, −58.67, −59.98 C²m² at E = 0.5, 1.5, 2 E_F (E′ = E_F, L = 4 nm) against about −58.6, −58.3, −60.0 read from the figure.
- At k_z′ → k_z the integrand diverges as (k_z − k_z′)⁻⁴. This is the continuum limit of a sum over discrete level differences s = ±1, ±3, …, which has no term at s = 0. It is replaced by [(k_z − k_z′)² + (c_reg/L)²]⁻². With c_reg = 1 the second moment equals the discrete sum (Σ_odd 1/s² = π²/8), which keeps the oscillator strength; the diagonal singularity of W becomes logarithmic, as [16] states for its own formula, and the scale 1/L is the one that appears in the logarithms of the printed Eq. (16).

**Sum rule.** For any dispersion, Σ_k′ (E_k′ − E_k)|z_kk′|² = (ħ²/2) ∂²E/∂p_z², whose angular average is ħ² w_opt(E)/2m₀ with the optical weight of [15] Eq. 9. The reconstructed W gives 0.84 of the bulk ω_p² for the cold 4 nm particle (0.90 at E_F), 0.91 at 8 nm and 0.95 at 16 nm (on grids fine enough to resolve the kernel, whose width in energy shrinks as 1/L). So the particle's effective plasma frequency is 8 % below the bulk value; this is a property of the quasi-continuum box model, not a numerical error.

**Weak-field permittivity, [16] Eq. (18).** As printed it lacks one factor 1/ħ (the prefactor must be V/2ε₀ħ²) and has −2iη̄ω in the denominator. With the e^{−iωt} convention used in `enzkin` a passive medium needs

  ε(ω) = ε∞ + Σ_p S_p (f_lo − f_hi) / (Δ_p² − ω² − 2i η̄_p ω).

**ε∞** is not printed in [16, 62]. It is fitted to Re ε(175 THz) = −1.6 of [16] Fig. 2b: ε∞ = 3.32 with the [14, 15] parameter set.

## 3. Validation (no parameter other than ε∞ adjusted)

Parameter set of [14, 15] (`ito_sarkar2023`), L = 4 nm, 6 eV band, 192 points, pulse of [16, 62]: 5 fs, 175 THz.

| Quantity | [16] / [62] | enzkin.densmat |
| --- | --- | --- |
| ε(175 THz), 300 K | −1.6 + 1.4i | −1.60 (fitted) + 1.27i |
| Zero of Re ε | 215 THz | 214 THz |
| Energy given to the electrons, E₀ = 1.2 V/nm | 0.63 eV/nm³ | 0.637 |
| Instantaneous absorbed power, E₀ = 1.2 V/nm | −0.14 … +0.33 eV/nm³/fs | −0.147 … +0.343 |
| Local-field spectral peak, E₀ = 0, 0.4, 0.6, 0.8, 1.2, 1.8, 2.4 V/nm | 179, 178, 177, 176, 173, 170, 168 THz | 177.3, 176.2, 174.9, 173.3, 170.3, 167.4, 166.0 |
| Local-field spectral width, same fields | 66, 70, 74, 79, 90, 105, 121 THz | 62, 67, 73, 78, 91, 107, 120 |
| ε(160 THz) before and after the pump, E₀ = 1.2 V/nm ([62] Fig. 4) | −2.6 + 1.85i → about −1.0 + 2.25i | −2.51 + 1.70i → −0.92 + 2.16i |
| 1/η_ee at E_F, 5 fs after the peak, E₀ = 1.8 V/nm ([62] Sec. 4.2) | about 12 fs | 14 fs |
| Effective temperature at that point | "∼4500 K" | 6050 K |

The last row is the only disagreement. The absorbed energy at 1.2 V/nm agrees to 1 %, so the energies at 1.8 V/nm should agree as well, and 1.27 eV/nm³ corresponds to 6050 K by energy matching. The e–e time of 12–14 fs is what the collision integral gives near 6000–6500 K (20 fs at 4500 K, `scripts/check_un2026.py`), so the "∼4500 K" of [62] does not fit its own rate; with 6000 K the earlier puzzle of working-log §4.2 disappears.

## 4. Closures built on the same oscillators

| Level | Populations | Polarisation | State per cell |
| --- | --- | --- | --- |
| L0 | kinetic, with collisions | all pair oscillators | N + N(N−1) numbers |
| L1 ("thermal model" of [62]) | Fermi–Dirac at T_e, dU/dt = E_loc·Ṗ − P_eph(T_e) | all pair oscillators, η(E; T_e) tabulated | 1 + N(N−1) |
| L2 (few-pole law) | T_e as in L1 | K fixed overlapping poles, strengths S_k(T_e) ≥ 0 fitted to χ(ω; T_e), widths γ_k + [η̄(T_e) − η̄(T₀)] | 1 + 2(K+1) |
| L3 (Drude, plan Tier 2a) | T_e as in L1 | one pole at zero frequency, bulk ω_p²(T_e), damping fitted to the cold particle at the carrier | 3 |

All are driven by the instantaneous work E_loc·Ṗ. In L2 the temperature dependence of the width matters: a fit with fixed widths reproduces χ(ω; T_e) equally well at every frozen temperature but errs in time, because the electron temperature changes within the memory time of the coherences (`docs/gates/GT2_fewcycle_closure_report.md`).

## 5. Limits of the model

- Absorption is by the finite size of the particle (Landau-type) only. Phonon-, impurity- and Umklapp-assisted absorption, which is what absorbs in a film, is absent ([16] says so). The model must not be applied to a film cell as it stands.
- A box stands in for a sphere, in the quasi-continuum limit. At 4 nm only four to five subbands are occupied at E_F, so the true level structure is discrete; the regularisation constant c_reg moves ε″ by about 5 % per 10 %.
- No field-driven coupling between coherences, hence no intra-cycle (Kane-mass) nonlinearity, although the field displaces the Fermi sphere by about k_F at 1.2 V/nm (`scripts/field_regime.py`).
- Intraband only; lattice fixed at 300 K.
