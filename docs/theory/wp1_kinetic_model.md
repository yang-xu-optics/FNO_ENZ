# WP1 kinetic model: equations as implemented

This note is the specification of the `enzkin` package (`src/enzkin`). It states the equations exactly as coded, the discrete forms that make the conservation laws hold to round-off, and every place where the implementation departs from a printed formula. Source digests and equation numbers are in [`references/LITERATURE_SUMMARY.md`](../../references/LITERATURE_SUMMARY.md); reference numbers `[n]` are those of the project plan.

## 1. Units and grid

Working units are eV, fs, nm, K (`constants.py`): ħ = 0.6582 eV fs, e²/ε₀ = 18.095 eV nm, mₑ = 5.6856 eV fs² nm⁻², densities in nm⁻³ (1 nm⁻³ = 10²⁷ m⁻³), fields in V/nm.

The energy grid is uniform and cell-centred, Eᵢ = (i + ½)ΔE, i = 0…N−1, ΔE = E_max/N with E_max = 4 eV (`grid.py`). Two reasons:

- cell centres avoid k = 0, where the screened-Coulomb prefactor 1/k diverges;
- on a uniform grid E₁ + E₂ = E₃ + E₄ is the index identity i₁ + i₂ = i₃ + i₄, so two-body energy conservation is exact, not interpolated.

Moments are midpoint sums: n = Σρᵢfᵢ ΔE, u = ΣEᵢρᵢfᵢ ΔE.

## 2. Band (`band.py`)

Kane dispersion ħ²k²/2m₀ = γ(E) = E(1 + CE), with m₀ = 0.3964 mₑ, C = 0.4191 eV⁻¹ [14]:

| Quantity | Expression |
| --- | --- |
| wavevector | k = √(γ/h₂ₘ), h₂ₘ = ħ²/2m₀ |
| density of states (both spins) | ρ = (1 + 2CE) k / (2π² h₂ₘ) |
| group velocity | v = 2h₂ₘ k / (ħ(1 + 2CE)) |
| optical weight | w_opt = (1 + 8Cγ/3)/(1 + 2CE)³ |

For n = 1.5 nm⁻³: E_F = 0.8805 eV, k_F = 3.541 nm⁻¹, T_F = 10 218 K, Sommerfeld coefficient 12.67 J m⁻³ K⁻² ([14] quotes 0.88 eV, 3.54 nm⁻¹, 12.7).

## 3. State and equation (`model.py`)

State y = [f(E₀)…f(E_{N−1}), T_ph, U_abs, U_pair].

```
df/dt        = (df/dt)_exc + Q_ee[f] + Q_eph[f, T_ph] + S_inter
C_ph dT_ph/dt = P_eph                       (power leaving the electrons through phonons)
dU_abs/dt     = p_abs + p_field,inter       (energy taken from the optical field)
dU_pair/dt    = interband pair power        (energy stored in electron-hole pairs)
```

The last two are bookkeeping variables. Carrying them in the state makes the audit

    U_e + C_ph T_ph + U_pair − U_abs = const

a linear invariant of the ODE system. Runge–Kutta methods preserve linear invariants exactly, so the audit and (without interband terms) the electron number hold to ~10⁻¹³ relative for any step size.

## 4. Electron–electron collisions (`collisions/ee.py`)

Both levels are Pauli-blocked two-body integrals, Q(E₁) = Σ W(1,2;3,4)[f₃f₄(1−f₁)(1−f₂) − f₁f₂(1−f₃)(1−f₄)] over E₂, E₃ with E₄ = E₁ + E₂ − E₃. After multiplication by ρ(E₁) both kernels are symmetric under (1,2)↔(3,4) and under (1↔2, 3↔4). Consequences, all checked in `tests/test_collisions.py`: number and energy rates vanish to round-off; every Fermi–Dirac function is an exact fixed point; entropy production is non-negative.

### B1: momentum-resolved, statically screened Coulomb

Golden rule with |V(q)|² = [e²/(ε₀ε_b(q² + κ²))]²; for isotropic f the angular integrals are analytic and k dk = m₀(1 + 2CE) dE/ħ²:

```
Q(E1) = P/k1 · Σ_{E2,E3} ΔE² (1+2CE2)(1+2CE3)(1+2CE4) Φ [ ... ]
Φ     = A(q_hi) − A(q_lo)  if positive, else 0;   A(q) = qκ/(q²+κ²) + atan(q/κ)
q_lo  = max(|k1−k3|, |k2−k4|),   q_hi = min(k1+k3, k2+k4)
P     = c_ee /(ħ(2π)³) · (e²/(ε₀ε_b))² · (m₀/(ħ²κ))³
κ²    = (e²/(ε₀ε_b)) ∫ρ(−∂f/∂E)dE          (recomputed from f at every call; ε_b = 4)
```

This is [14] Eq. 15–17 with three corrections:

1. **Screening.** [14] Eq. 15 prints an extra factor 4π³ in q_TF²; a footnote of [16] says it should be removed. Without it κ = 3.83 nm⁻¹ at 300 K.
2. **Dimensions.** [14] Eq. 16–17 are printed with 1/ħ² where 1/ħ is needed; [14] Eq. 19 (the Fermi-liquid coefficient K) is consistent with the form above.
3. **Prefactor.** A spin-summed direct-term derivation gives c_ee = 1. The prefactor printed in [14] corresponds to c_ee = 3/16, which is the default so results compare with [14, 15]; it gives K = 0.0651 eV⁻² fs⁻¹. c_ee is a material parameter (the "scattering prefactor" of model M1).

κ² is evaluated on cell edges, Σ ρ(edge)(fᵢ − fᵢ₊₁), which is exact summation by parts.

Evaluation: the antiderivative A is tabulated at k_a + k_b and |k_a − k_b| along diagonals (d = i₂ − i₄) so the innermost loop is a contiguous min/max/multiply reduction (numba, parallel over E₁). Cost is O(N³): 1.8 ms at N = 128, 13 ms at 256, 80 ms at 512 on 4 CPU cores.

A second implementation of the same discrete integral, `collisions/ee_jax.py`, is written in JAX: for each shift d the kernel is an N×N matrix built from two vectors and the sums are two matrix–vector products. It is batched over cells, differentiable, and meant for GPUs; it agrees with the numba kernel to 10⁻¹¹ and keeps the same invariants (`tests/test_ee_jax.py`). On a CPU it is about half as fast as numba.

### B0: random-k, constant matrix element [40]

```
Q(E1) = K0 · Σ_{E2,E3} ΔE² ρ2 ρ3 ρ4 [ ... ]
```

With a = ρf and b = ρ(1−f): gain = (1−f)·corr(b, a∗a), loss = f·corr(a, b∗b), evaluated with FFT convolutions in O(N log N) (0.25 ms, almost independent of N up to 512). The FFT result equals the direct triple sum to 10⁻¹⁴. K₀ is set so that B0 and B1 share the Fermi-liquid coefficient, K₀ρ_F³/2 = K(κ), and follows κ[f]; a B0–B1 comparison therefore tests the kernel's shape, not its strength.

The two kernels agree near the Fermi level and part ways with excitation energy: at μ + 0.5 eV the B0 rate is twice the B1 rate, because B0 grows with the cube of the Kane density of states while the Coulomb kernel suppresses large momentum transfer.

### RTA baseline

`ElectronElectronRTA`: −(f − f_FD)/τ(E) with the Fermi-liquid τ(E) and f_FD matched in density and energy. As [40] notes and [39] quantifies, an energy-dependent τ on a single distribution violates conservation; `project_conserving` removes the violating component along h(1−h) and E·h(1−h). The same projection is the hard constraint intended for the WP2 network.

## 5. Electron–phonon collisions (`collisions/eph.py`)

Acoustic deformation-potential scattering with a Debye spectrum [14 Eq. 5–11]: D = 17.2 eV, ρ_m = 7120 kg m⁻³, v_s = 6400 m s⁻¹, E_D = 60.3 meV, lattice in equilibrium at T_ph. The same authors later use D = 8.6 eV with the same formula, together with Brooks–Herring impurity scattering ([62] Table S1; preset `ito_un2026()`); G scales as D², so the two sets differ by four in the electron cooling rate. Which value describes a given film is a WP4 inference question (see `docs/theory/un2026_notes.md`).

Multiplying [14] Eq. 9 by ρ(E) shows that the kernel is symmetric in the two electron states,

    d(ρf)/dt |_E = κ_ep (1+2CE) ∫dε ε² (1+2CE′)[…],   κ_ep = m₀²D²/(4π³ħ³ρ_m(ħv_s)⁴),

so the term is implemented as a master equation between pairs (Eᵢ, Eⱼ = Eᵢ + sΔE):

```
T_ij = W_ij [ (n_s + 1) f_j (1 − f_i) − n_s f_i (1 − f_j) ]          n_s = n_B(sΔE, T_ph)
W_ij = κ_ep (1+2CE_i)(1+2CE_j) ∫ ε² hat_s(ε) dε   over   ħv_s|k_j − k_i| ≤ ε ≤ min(E_D, ħv_s(k_i + k_j))
```

hat_s is the linear hat function centred on sΔE, so phonon energies between grid nodes are represented and ∫ε³ (the energy-loss moment) is integrated exactly.

**Departure from [14].** The printed cut-off ε ≤ min(E_D, 2ħv_s k(E)) depends on one state only; with it the term does not conserve electron number and overestimates the coupling (G ≈ 8×10¹⁶ W m⁻³ K⁻¹). The symmetric triangle condition above differs at O(ε/E) and gives: exact number conservation; Fermi–Dirac at T_ph as an exact fixed point (detailed balance); the lattice power Σ T_ij ε_s known exactly; and G(T→0) = D²m₀²k_B k_F⁴(1+2CE_F)²/(π³ρ_mħ³) = 2.51×10¹⁶ W m⁻³ K⁻¹, which is [14] Eq. 27. The numerical G is 2.44, 3.33, 6.67 ×10¹⁶ at T_e = 500, 2000, 5000 K (N = 512) and changes by 0.7 % between N = 256 and 1024. Neglecting momentum conservation raises G 14-fold ([14] says "~30"; its own formula gives 16.7).

**Collision rate for the damping.** 1/τ_eph(E) of [14] Eq. 11 counts every scattering event, including quasi-elastic ones whose phonon energy is below one grid cell. Those do not move population between cells, so the master equation omits them, but they relax momentum and dominate the rate (the integrand ε²n_B(ε) ≈ ε k_BT is weighted to small ε). `ElectronPhonon.rate` therefore adds the intra-cell (s = 0) term and integrates the Bose factor over each hat function instead of sampling it at the node:

```
1/τ_i = (1/ρ_i) Σ_{s≥0} [ B_s(i,i+s) + W_s(i,i+s) f_{i+s} + B_s(i−s,i) + W_s(i−s,i)(1 − f_{i−s}) ],   B_s = κ_ep m m′ ∫ε² n_B(ε,T_ph) hat_s dε
```

The rate is 0.0585 fs⁻¹ at E_F and 300 K, proportional to T_ph, and its optical average is the same to three digits from N = 64 to 1024. (The first implementation sampled n_B at the nodes and had no s = 0 term; it gave a rate three times too small at N = 128 and was the reason the damping was not grid-converged.)

## 6. Photo-excitation (`excitation.py`)

[14] Eq. 4–6: electrons move up by one photon energy with a constant matrix element, and the overall rate B(t) is fixed by ∫Eρ(∂f/∂t)_exc dE = p_abs. A transition from cell j lands at E_j + ħω = E_j + (m + w)ΔE and is split between cells j+m and j+m+1 with weights (1−w, w), so number is conserved for any photon energy. Two variants: `absorption_only` (as printed in [14, 15]) and `net` (absorption minus stimulated emission, Pauli factors combining to f(E−ħω) − f(E)).

## 7. Drives (`drive.py`, `tmm.py`)

| Drive | p_abs(t) | Use |
| --- | --- | --- |
| `PowerDrive` | prescribed | tests; kinetics-only comparisons |
| `FieldDrive` | (ω/2)ε₀ε″[f]·\|E_loc(t)\|², in-film intensity defined as \|E_loc\|²/2Z₀ | WP1 dataset |
| `ThinFilmDrive` | same, with ⟨\|E_loc\|²⟩ from a transfer-matrix solution of the layer stack at the instantaneous ε[f] | adiabatic self-consistent scheme of [15]; replaced by FDTD in WP3 |

## 8. Optics (`optics.py`)

```
ω_p²[f] = e²/(ε₀m₀) ∫ρ f w_opt dE                         [15] Eq. 10 (by parts: e²/(3ε₀)∫ρv²(−∂f/∂E)dE)
ε(ω)    = ε_∞ − e²/(ε₀m₀) ∫ρ f w_opt /(ω + iη(E))² dE       [15] Eq. 9
η(E)    = 1/τ_ee(E) + 1/τ_eph(E) + γ_imp [+ η_imp^BH(E) if n_imp > 0]
```

The optional Brooks–Herring term (`collisions/eimp.py`) is the out-scattering rate of elastic screened-Coulomb scattering on n_imp ionised donors of charge Z, screened by ε_b with the κ of the current f:
η_imp(E) = (2π n_imp/ħ)(Ze²/ε₀ε_b)²(2π)⁻³ (π/k²)(dk/dE) [(1 + κ²/2k²)² − 1]⁻¹. This is [62] Eq. S9 with the sign inside the bracket corrected (printed 1 − κ²/2k²); the test `test_brooks_herring_rate_matches_elastic_sphere_integral` checks it against the angular integral. Being elastic, it changes η(E) but not f(E).

To first order in η/ω this is a Drude term with Γ = 2⟨η⟩ (`drude_equivalent`). ε_∞ and γ_imp are not printed in [14, 15]; matching [15]'s plotted ε(230 THz, 300 K) = 0.172 + 0.399i gives ε_∞ = 3.451 and γ_imp = 0.0352 fs⁻¹ (`enzkin.calibrate`), with the zero crossing at 224 THz ([15]: ≈225 THz).

## 9. Interband channels (`interband.py`)

| Model | Term | Published form | Kinetic lift used here |
| --- | --- | --- | --- |
| M2 | `AugerLikeRelaxation` | ṅ_c = η[n₀(T) − n_c], U̇_c = −E_g ṅ_c, η = 1/120 fs⁻¹, E_g = 3.9 eV [22] | number and energy rates injected along the Fermi–Dirac tangent directions h(1−h), E·h(1−h); n₀(T) from charge neutrality with a parabolic valence band |
| M3 | `TwoPhotonAbsorption` | β I² N₁(1 − N₂/N₂max), gap 2.67 eV [23] | pairs generated at β I²/(2ħω)·(1 − f(E*)), E* fixed by a vertical transition; Pauli blocking and its lifting by intraband absorption arise from f itself |

The valence-band mass is printed in neither paper and is a declared free parameter (default 1 mₑ). Two-photon absorption is closed for 2ħω < E_g: at E_g = 2.67 eV the channel opens only below 929 nm.

## 10. Thermal baseline, model M0 (`ttm.py`)

U̇_e = p_abs − P_eph(T_e, T_ph), C_ph Ṫ_ph = P_eph, with f a Fermi–Dirac function at all times. U_e(T_e), μ(T_e) and P_eph come from the same grid and the same electron–phonon operator as the kinetic model, so a kinetic and a thermal run differ only by non-thermal electrons.

## 11. Time integration (`integrate.py`)

Fourth-order Runge–Kutta, with each step limited to h ≤ 1/λ, where λ is the larger of the fastest collision rate and the rate that would change any f(E) by 2 % in one step. Near equilibrium the rates are below 0.1 fs⁻¹; a strongly excited distribution under B0 scatters at several per fs, and a fixed 0.5 fs step is unstable there. Outputs land on a uniform sampling grid.

## 12. Known limitations

- **Single band of finite width.** The grid ends at E_max (4 eV as stated in [14]). A band of that width holds at most about 2.7 eV per electron, i.e. ≈5×10⁸ J/m³ above the cold state for n = 1.5 nm⁻³. As the absorbed energy approaches that, the heat capacity collapses, the matched temperature runs away and the excitation term eventually has nothing to excite. In the [15] validation this shows as a 16–32 % overshoot of T_e,max at 50–75 GW/cm², removed (to 5 %) by E_max = 6 eV. Keep ΔE fixed when widening the band (N = 192 at 6 eV corresponds to N = 128 at 4 eV). Beyond that, real ITO has higher bands and photoemission, which this model does not.
  The convergence sweep (`docs/gates/G1b_report.md`) turned this into a rule. At a fixed 31 meV step the permittivity transient changes by less than 1 % when the band is widened, as long as the occupation of the top cell never exceeds 4×10⁻⁵ during the run; above 2×10⁻⁴ it changes by more than 1 %. In temperature terms a 4 eV band is safe to about 3500 K and a 6 eV band to about 6700 K. At the hottest conditions tested (≈8500 K) the transient still changes by 3.1 % from 6 to 8 eV and by 1.3 % from 8 to 10 eV, and the error of a narrow band is always upward. Since the Kane dispersion is an extrapolation at those energies, results above roughly 7000 K depend at the percent level on states the model does not describe.
- **Energy step.** With ΔE = 31 meV the cold Fermi edge (k_B T = 26 meV at 300 K) is one cell wide. The resulting quadrature error is a common scale factor on the transient of T_e, ω_p² and ε: about 1 % for a temperature rise of 130–190 K and 0.4 % or less above 1800 K; halving ΔE divides it by about 3. The cold-state permittivity also shifts by a constant between grids (10⁻³ in Re ε at 128 points), so ε∞ and γ_imp must be calibrated on the grid in use.
- Static Thomas–Fermi screening; no exchange; no dynamic screening.
- e–e collisions enter the optical damping at full weight as in [15]; [16] argues only Umklapp e–e processes should.
- Isotropic f only: the current-carrying (l = 1) part of the distribution is not evolved, so the choice between the two current equations of WP3 is not yet decided by the kinetics (see LITERATURE_SUMMARY Part I §6).
- Cycle-averaged absorption; valid for pumps ≳ 20 fs.
- Acoustic phonons only; no Umklapp.
- Interband terms are moment-level lifts of lumped models, not microscopic impact-ionisation or two-photon kernels.
