# Notes on Un, Sarkar & Sivan, Nanophotonics 15, e70272 (2026) — ref. [62]

What in the paper's printed formulas had to be corrected to use or check them, and what was checked against enzkin. Script: `scripts/check_un2026.py` → `results/wp1/check_un2026.json`. Digest: `references/LITERATURE_SUMMARY.md`, entry [62].

## Corrections to printed formulas (Supporting Information)

| Eq. | As printed | Used here | Why |
| --- | --- | --- | --- |
| S9 (impurity rate) | [(1 − q_TF²/2k²)² − 1]⁻¹ | [(1 + q_TF²/2k²)² − 1]⁻¹ | The angular integral of (q² + q_TF²)⁻² over the elastic sphere, q² = 2k²(1 − cos θ), gives 2/[q_TF²(4k² + q_TF²)] = (1/2k⁴)[(1 + q_TF²/2k²)² − 1]⁻¹. The printed form is negative for q_TF < 2k (e.g. at E_F, q_TF = 3.8 nm⁻¹, 2k_F = 7.1 nm⁻¹). Checked numerically in `tests/test_collisions.py`. |
| S22 (Fermi–Dirac product) | cosh(x/2)/[cosh(x₁/2)cosh(x₂/2)cosh(x₃/2)] | the same divided by 4 | With f = e^{−x/2}/(2cosh(x/2)) the two products sum to cosh(x/2)/[4 cosh cosh cosh]. Without the 1/4 the coefficients come out four times those of Fig. S1. |
| S23 | tanh((𝓔₁ − 𝓔₂)/k_BT̃_e) | tanh((𝓔₁ − 𝓔₂)/2k_BT̃_e) | Second-order expansion of cosh(x/2)/cosh((x+u)/2) in x gives 1 − ½tanh(u/2)x + ¼tanh²(u/2)x²; (S26) has the /2 already. |
| S26 (A_n) | ∬ d𝓔₁d𝓔₂ … without a temperature factor | multiplied by 1/(k_BT̃_e)² | A_n must be dimensionless; with the factor A₀ → π² and A₂ → 1 at low T as the paper states and Fig. S1 shows. |
| S25 (K) | 3/(32ħ)(m*/ħ²q_TF)³(e²/ε₀ε_b)²(1+2Cμ)³/k_μ [Φ] | the same times (2π)⁻³ | The T → 0 limit of (S19), whose prefactor is (3/4)(2π/ħ)((2π)²/2)(2π)⁻⁶(1/2q_TF³) = 3/[16(2π)³ħ]·(1/q_TF³)…; enzkin's `fermi_liquid_K` is this limit. As printed K would be 248 times too large. |
| S2 vs S7 | (S2) has (3/4)(2π/ħ)(1/2)Σ…, (S7) has (3/4)(2π/ħ)Σ… | not used | The rate should be the functional derivative of the collision term, so the two should carry the same factor. They differ by 2 as printed. (S19), the rate actually plotted, carries c_ee = 3/16 in enzkin's convention. Open question for the authors. |
| S3–S4 vs S6, S9, S18, S25 | screening with ε∞ in S3–S4, with ε_b elsewhere | ε_b = 4 | Neither constant is in Table S1. |
| Table S1 | ρ = 7120 kg/m² | kg/m³ | unit typo |

## Checks against enzkin

1. **e–e kernel and prefactor.** (S19) is the B1 integral of `collisions/ee.py` with c_ee = 3/16, term by term.
2. **Generalised Fermi-liquid coefficients** (independent implementation, 600-point integration grid, 4 eV band): A₀(E_F) = 9.46, 8.51, 6.06, 4.65, 3.73 and A₂(E_F) = 0.95, 0.83, 0.53, 0.40, 0.32 at 300, 1000, 3000, 5000, 7000 K, against Fig. S1 (by eye) 9.4, 8.6, 6.3, 4.9, 4.0 and 0.95, 0.85, 0.55, 0.42, 0.35. The remaining 5–7 % at high T may come from the unprinted band top or screening constant. The formula reproduces the B1 rate of a Fermi–Dirac distribution within 3 % at E_F and 1.5E_F, and within 30 % at 0.5E_F at 300 K.
3. **e–e magnitude.** B1 gives 1/η(E_F) = 20 fs at 4500 K and 12 fs at about 6500 K; the paper quotes 12 fs at about 4500 K. Its T̃_e(t) is not printed, so this 1.5–2× difference cannot be traced.
4. **e–ph rate.** With D = 8.6 eV, enzkin gives 0.0028, 0.015, 0.041, 0.069 fs⁻¹ at 0.1, 1, 2, 3 E_F (300 K); Fig. 5c shows about 0.001, 0.01, 0.03, 0.06. With D = 17.2 eV enzkin's rate is six times higher. The low-energy difference comes from the quasi-elastic (sub-cell) phonons that enzkin includes (`wp1_kinetic_model.md` §5).
5. **Damping budget and consequences** (bulk film, 300 K, ε∞ = 3.45):

| | e–e (E_F) | e–ph (E_F) | impurity (E_F) | Drude Γ | ε(230 THz) | G at T_e = 1000 K | Re ε recovery, 100 fs pulse, 10⁸ J/m³ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [14, 15] set (D 17.2 eV, γ_imp 0.0352 fs⁻¹ fitted to [15]) | 0.0004 | 0.0588 | 0.0352 | 0.175 fs⁻¹ | 0.172 + 0.399i | 2.6×10¹⁶ W m⁻³ K⁻¹ | 0.52 ps |
| [62] Table S1 (D 8.6 eV, Brooks–Herring 0.05 nm⁻³) | 0.0004 | 0.0147 | 0.0167 | 0.074 fs⁻¹ | 0.142 + 0.169i | 6.6×10¹⁵ W m⁻³ K⁻¹ | 1.9 ps |

Rates in fs⁻¹. Measured recovery in ITO films is 60–400 fs ([1, 6, 21, 60]). The linear permittivity fixes only the total damping; the e–ph share sets the cooling. The two sets come from the same group and differ by four in that share, so the project treats D (or G) and the impurity model as inferred, not fixed (plan revision 2.1, decision D15).

## What is not adopted

- The density-matrix photon–electron term and its finite-size (box) dipoles: they describe absorption in a nanoparticle, not in a film.
- The paper's thermal model is not implemented as such; the plan's Tier 2a uses its idea (drive the electrons by the instantaneous J·E) inside the WP3 film solver.
