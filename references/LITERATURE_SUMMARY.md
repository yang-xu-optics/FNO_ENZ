# Literature summary: theory cornerstone for the ENZ kinetic-closure project

Built 2026-11-17 from the 61 references of `docs/PROJECT_PLAN.md`; [62] (Un, Sarkar & Sivan, Nanophotonics 2026) added 2026-11-18. Reference numbers `[n]` are the plan's.

**How to use this file.** Part I is the synthesis: the physics the project stands on, the equations and parameter values to implement, where printed formulas are wrong, and where the literature disagrees with the project plan. Part II holds one digest per paper with equations, numbers and their locations. When Part I and a digest differ, the digest is closer to the source; when a digest and the paper differ, the paper wins.

**How it was produced, and how far to trust it.** 53 papers were read in full text: 52 from PDF (with arXiv TeX source where it exists; [62] with its Supporting Information, supplied by the project lead) and [20] from Europe PMC full-text XML; the 97-page [43] was read selectively, as its digest states. [35] and [57] were available as abstracts only. Seven paywalled papers ([1], [3], [4], [7], [32], [33], [34]) could not be obtained; their digests are reconstructions from the abstract (where one could be retrieved) and from how the other papers in the library cite them, and are flagged "NOT READ". Supplementary materials were mostly unavailable, and for papers with an arXiv and a journal version the arXiv version was read. Each digest states exactly what was read; the appendix lists it for all 62. The reading was done by ten parallel reading passes with a common template; statements marked "(my reading)", "(my calculation)" or "(background)" inside digests are the reader's inferences, not the paper's. Items in Part I marked **[verified in code]** were re-derived and checked numerically while building the WP1 solver (`src/enzkin`, `tests/`).

---

# Part I. Synthesis

## 1. The physical picture

1. **ITO near its ENZ wavelength is a degenerate electron gas in a non-parabolic band.** Typical films: n ≈ 1–1.5×10²¹ cm⁻³, band-edge mass 0.26–0.40 mₑ, Kane non-parabolicity C ≈ 0.42 eV⁻¹, E_F ≈ 0.9–1 eV (T_F ≈ 10⁴ K), ε∞ ≈ 3.5–4, zero crossing at 1200–1550 nm [1, 5, 6, 9, 14, 23].
2. **Sub-gap light heats the conduction electrons; hot electrons are heavier on average; the plasma frequency drops and the ENZ point redshifts** [5–9, 14, 15]. Because Re ε ≈ 0 initially, the *relative* change of ε, and hence of the index, is of order one [1, 2]. The mass that matters is the transport mass, ω_p² ∝ ∫ρ f w_opt dE ([15] Eq. 9–10), identical to the averaged-mass formulas of [5, 6, 8] after integration by parts **[verified in code]**.
3. **The electron heat capacity is small** (γ_e ≈ 12.7 J m⁻³ K⁻², five times less than gold), so electrons reach several thousand kelvin at tens of GW/cm² and a sizeable fraction of T_F at higher intensity; the chemical potential falls and can turn negative [14, 15] **[verified in code]**.
4. **Electron–electron scattering is fast** (Fermi-liquid coefficient K ≈ 0.065 eV⁻² fs⁻¹, thermalisation in tens of fs, about ten times faster than in gold) [14]. A two-temperature description is therefore adequate for pulses of ≳100 fs, and non-thermal electrons matter for the *rise* under ≲30 fs pumps [15]. Under a 5 fs pump the rate rises to ≈1/12 fs at ≈4500 K [62], and a generalised Fermi-liquid formula with a finite band bottom reproduces the full integral (A₀ < π², a negative linear term) [62] **[verified in code against B1]**. The plan's phrase "attributes the real-permittivity change to non-thermal electrons [14, 15]" overstates this: at 220 fs the thermal and kinetic permittivities coincide.
5. **Electron–phonon energy transfer is slow and restricted by momentum conservation**: k_F ≈ q_D/4, so only phonons with q ≤ 2k participate; G ≈ 2.5×10¹⁶ W m⁻³ K⁻¹ at low T_e, rising with T_e [14] **[verified in code]**. Recovery of the response is electron cooling on a few hundred fs [1, 5, 6, 9, 15, 21]. **The e–ph strength is not settled**: the same group's later parameter set [62] halves the deformation potential (D = 8.6 instead of 17.2 eV, G four times smaller) and moves the missing damping into Brooks–Herring impurity scattering; with it a 100 fs transient recovers in 1.9 ps instead of 0.5 ps, slower than any measured ITO film (60–400 fs) **[verified in code]**.
6. **"Instantaneous" rise means ε′ follows the absorbed energy**, i.e. the running integral of the pump, not the electron temperature of a thermalised gas [5, 15]. In [15] the 10–90 % rise is ≈210–230 fs for a 220 fs pump and 37 fs for a 30 fs pump. No microscopic model in the library predicts a rise *faster* than the pump. Under a single-cycle pump [62] finds the same cumulative behaviour, decorated by a 2ω ripple from coherences, and reproduces it with a thermal model driven by the instantaneous absorbed power Ṗ·E.
7. **At ≳1 TW/cm² and sub-10 fs pumps the response changes sign or oscillates** [20, 22, 23, 24]. Three explanations compete (Section 5). They differ on whether the conduction-band density changes.
8. **What is measured is not ε(t).** Every experiment reads the response through a film: transit-limited frequency shifts [17], interface phase terms that reverse the shift's sign [18], a saturated Berreman mirror whose *reflectivity* rise was fitted [19], thickness-dependent red shifts [20], and non-collinear geometric smearing of 100 fs–1 ps [60]. This is hypothesis H3, and the literature supports it piecemeal but nobody has tested it with a distribution-resolved model.
9. **The ultimate response time is acknowledged as open** by the reviews themselves ("still an open research question" [26]; "unexpectedly fast… sparked ongoing theoretical investigations" [29]).

## 2. Governing equations adopted in WP1 (with corrections to the printed forms)

State per spatial cell: f(E,t), T_ph(t). Units and discretisation are documented in `docs/theory/wp1_kinetic_model.md`.

**Band** ([14] Eq. 1–2; [15] Eq. 7–8):

$$\frac{\hbar^2k^2}{2m_0}=E(1+CE),\qquad \rho(E)=\frac{1+2CE}{2\pi^2}\Big(\frac{2m_0}{\hbar^2}\Big)^{3/2}\sqrt{E(1+CE)},\qquad m_0=0.3964\,m_e,\ C=0.4191\ \mathrm{eV^{-1}}.$$

**Boltzmann equation** ([14] Eq. 3):

$$\partial_t f=\Big(\partial_t f\Big)_{exc}+Q_{ee}[f]+Q_{e\text{-}ph}[f,T_{ph}]+S_{inter}[f,\dots],\qquad C_{ph}\,\dot T_{ph}=-\int E\rho\,Q_{e\text{-}ph}\,dE .$$

**Excitation** ([14] Eq. 4–6; [15] Eq. 2–3), upward transitions by one photon, normalised to the absorbed power density:

$$\Big(\partial_t f\Big)_{exc}=B(t)\big[\rho(E-\hbar\omega)f(E-\hbar\omega)(1-f(E))-\rho(E+\hbar\omega)f(E)(1-f(E+\hbar\omega))\big],\quad \int E\rho\,(\partial_t f)_{exc}dE=p_{abs}=\tfrac{\omega}{2}\varepsilon_0\varepsilon''|E_{loc}|^2 .$$

The matrix element D_J is never specified in [14, 15]; the plan's wording "phonon- and impurity-assisted absorption [14]" is an interpretation. [16] calls this normalisation "ad hoc" and replaces it with an explicit dipole element for a nanoparticle. Stimulated emission is omitted in [14, 15]; [15] notes this fails when electrons pile up at the band top.

**Electron–electron, screened Coulomb (level B1)** ([14] Eq. 15–17):

$$Q_{ee}(E_1)=\frac{P}{k_1}\!\int\! dE_2\,dE_3\,(1+2CE_2)(1+2CE_3)(1+2CE_4)\,\Phi\,\big[f_3f_4(1-f_1)(1-f_2)-f_1f_2(1-f_3)(1-f_4)\big],\quad E_4=E_1+E_2-E_3,$$

$$\Phi=\Big[\frac{q\kappa}{q^2+\kappa^2}+\arctan\frac{q}{\kappa}\Big]_{q_{lo}}^{q_{hi}},\quad q_{lo}=\max(|k_1-k_3|,|k_2-k_4|),\ q_{hi}=\min(k_1+k_3,k_2+k_4),$$

$$P=c_{ee}\,\frac{1}{\hbar(2\pi)^3}\Big(\frac{e^2}{\varepsilon_0\varepsilon_b}\Big)^2\Big(\frac{m_0}{\hbar^2\kappa}\Big)^3,\qquad \kappa^2=\frac{e^2}{\varepsilon_0\varepsilon_b}\int\rho\,(-\partial_E f)\,dE,\qquad \varepsilon_b=4 .$$

Corrections **[verified in code]**: (i) [14] Eq. 15 carries a spurious factor 4π³ in q_TF² — [16] says so in a footnote; without it κ = 3.83 nm⁻¹. (ii) [14] Eq. 16–17 are printed with one power of ħ too many; [14] Eq. 19 for K is dimensionally right. (iii) The printed prefactor corresponds to c_ee = 3/16; a spin-summed direct-term golden-rule derivation gives c_ee = 1. With c_ee = 3/16, K = 0.065 eV⁻² fs⁻¹. Prefactors across the literature differ by factors of 3 ([40]: textbook Fermi-liquid prefactor had to be cut 2.8× to match the full integral; [37]: D_e = 0.016–0.021 eV⁻¹), so **c_ee is a fit parameter**, as the plan's M1 row assumes. [62] prints the energy-space rate (S19) with exactly the B1 kernel and c_ee = 3/16 **[verified in code]**; its population equation (S2) carries an extra ½ relative to its rate (S7) as printed, which sharpens rather than settles the question.

**Electron–electron, random-k (level B0)** ([40] Eq. 14–17):

$$Q_{ee}(E_1)=K_0\!\int\! dE_2\,dE_3\,\rho_2\rho_3\rho_4\,\big[f_3f_4(1-f_1)(1-f_2)-f_1f_2(1-f_3)(1-f_4)\big].$$

With a = ρf, b = ρ(1−f) this is one convolution and one correlation: gain = (1−f)·corr(b, a∗a), loss = f·corr(a, b∗b), O(N log N) by FFT on a uniform grid. [40] itself is silent on numerics; the reduction is ours and holds to 10⁻¹⁴ against direct summation **[verified in code]**.

**Electron–phonon** ([14] Eq. 9–11), acoustic deformation potential with D = 17.2 eV, ρ_m = 7120 kg m⁻³, v_s = 6400 m s⁻¹, E_D = 60.3 meV:

$$Q_{e\text{-}ph}(E)=\frac{D^2}{4\pi\rho_m(\hbar v_s)^4}\frac{m_0}{\hbar k(E)}\int_0^{\epsilon_{max}}\!\!\epsilon^2d\epsilon\Big\{(n_B+1)\big[(1-f(E))f(E+\epsilon)m(E+\epsilon)-(1-f(E-\epsilon))f(E)m(E-\epsilon)\big]+n_B\big[\dots\big]\Big\},$$

with m(E) = 1+2CE. Correction **[verified in code, and independently by the reading pass]**: with the printed cut-off ε_max = min(E_D, 2ħv_s k(E)) the term does not conserve electron number and gives G ≈ 8×10¹⁶; with the pair-symmetric triangle condition ħv_s|k−k′| ≤ ε ≤ min(E_D, ħv_s(k+k′)) it conserves number exactly and reproduces the published 2.5×10¹⁶ W m⁻³ K⁻¹ ([14] Eq. 27). Other flags in [14]: Eq. 26 is garbled; the "~30×" reduction of G by momentum conservation is 14–17× by the paper's own formula; "Umklapp treated as in [14]" in the plan has no counterpart — [14] never mentions Umklapp (only a footnote in [16] does).

**Electron–impurity (damping only)** ([62] Eq. S6, S9): Brooks–Herring scattering on ionised donors, elastic, so it changes η(E) but not f(E). Out-scattering rate η_imp(E) = (2π n_imp/ħ)(Ze²/ε₀ε_b)²(2π)⁻³(π/k²)(dk/dE)[(1 + κ²/2k²)² − 1]⁻¹; the printed bracket has a minus sign that makes the rate negative for κ < 2k **[verified in code: closed form against the angular integral]**. [62] uses n_imp = 0.05 nm⁻³, Z = 1 together with D = 8.6 eV; [15] used an unprinted constant rate. Implemented as an option (`enzkin/collisions/eimp.py`, `ito_un2026()`).

**Optics from f** ([15] Eq. 9–10):

$$\varepsilon(\omega)=\varepsilon_\infty-\frac{e^2}{\varepsilon_0m_0}\int\frac{\rho f\,w_{opt}}{(\omega+i\eta(E))^2}dE,\quad w_{opt}=\frac{1+\tfrac83CE(1+CE)}{(1+2CE)^3},\quad \eta=\tau_{ee}^{-1}+\tau_{e\text{-}ph}^{-1}+\tau_{imp}^{-1}.$$

To first order in η/ω this is a Drude term with damping Γ = 2⟨η⟩. Neither ε∞ nor τ_imp is printed in [15]; fitting its own 300 K curve (ε = 0.172 + 0.399i at 230 THz) with the printed band parameters gives ε∞ = 3.45 and τ_imp⁻¹ = 0.035 fs⁻¹ **[verified in code]**; reading ε∞ = 4 off the same curve instead implies ħω_p = 1.87 eV, 8 % above what the printed n, m₀, C give (1.73 eV). This gap is unexplained and needs the authors. [16] argues that only Umklapp e–e collisions should damp the polarisation, whereas [15] adds τ_ee⁻¹ at full weight.

**Two-temperature limit (model M0)** ([14] Eq. 20–21, [15] App. B): C_e(T_e)Ṫ_e = p_abs − G(T_e−T_ph), C_ph Ṫ_ph = G(T_e−T_ph), with C_e and μ(T_e) from the full Fermi integrals and C_ph = 2.54×10⁶ J m⁻³ K⁻¹. Terminology trap: in [38] "extended TTM" means adding non-thermal electrons; the thermal-only model is the plain TTM, so the plan's M0 label "extended two-temperature model [38]" should read "two-temperature model with full Fermi integrals".

**What reproducing [15] showed** (`scripts/validate_un2023.py`, `docs/WORKLOG.md` §1.5) **[verified in code]**:

- The modelled stack is air → glass → 40 nm ITO → 100 nm Au with the pump at 60° *in air* (35.3° in the glass). With it the self-consistent absorptivity matches [15] within 1–6 % from 2.5 to 75 GW/cm².
- Peak electron temperature matches within 2–7 % up to 22 GW/cm². At 50–75 GW/cm² it matches (to 5 %) only if the conduction band is taken 6 eV wide; with the 4 eV stated in [14] the band saturates and T_e overshoots by 16–32 %. The band width is an unprinted, consequential parameter.
- [15]'s ε″(t) is matched to 1–6 % with the Thomas–Fermi wavevector frozen at its 300 K value, and is 5–32 % too high if κ follows the hot distribution as [14] Eq. 15 literally prescribes.
- Δε′ comes out 20–28 % smaller than in [15] at every intensity, in line with the unexplained 8 % plasma-frequency gap noted above.

## 3. ITO parameter sets in the library

| Source | Film | n (10²¹ cm⁻³) | band-edge m*/mₑ | C (eV⁻¹) | ε∞ | Drude ω_p, γ | zero crossing |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [14, 15] | 40 nm on glass under 100 nm Au (model) | 1.5 | 0.3964 | 0.4191 | not printed (3.45 fitted here) | ħω_p = 1.73 eV from these values | 225 THz (1332 nm) |
| [1] via [9] | 310 nm on glass | — | 0.4 (from [7]) | 0.4191 | 3.8055 | 2π×473 THz, 0.0468 ω_p | 1237–1240 nm |
| [5] | generic | 1.1 | — | — | ≈4 | — | — |
| [6] | generic | 1.13 | 0.26 (edge), 0.49 (average) | 0.41 | 3.66 | — | 1325 nm |
| [10] | 23 nm | — | 0.26 | 0.418 | — | 1.9 eV, 0.17 eV | 1240 nm |
| [13] | — | — | — | — | 3.99 | 3.16×10¹⁵, 0.185×10¹⁵ rad/s | 1200 nm |
| [17] | 2×310 nm clamped | — | — | — | ≈3.8 (reader's fit) | ≈2.95×10¹⁵, ≈1.4×10¹⁴ s⁻¹ (reader's fit) | ≈1240 nm |
| [22] | 310 nm on glass | not stated | "modified Kane" | not stated | not stated | ≈3.03×10¹⁵, 0.11–0.14×10¹⁵ rad/s (from figure) | ≈1300 nm |
| [23] | 270 nm on glass | 1.0 | 0.35 | — | 3.9 | γ = 1.24×10¹⁴ rad/s | 1240 nm |
| [62] Table S1 | 4 nm sphere (model); D fitted to the bulk permittivity of [7]'s film | 1.5 | 0.3964 | 0.4191 | not printed | D = 8.6 eV, Brooks–Herring n_imp = 0.05 nm⁻³: Γ ≈ 0.074 fs⁻¹ in bulk with ε∞ = 3.45 (my calculation, against 0.175 fs⁻¹ for [15]) | — (sphere resonance ≈180 THz) |

The band-edge mass is the least settled parameter (0.26 vs 0.35–0.40); a comment in the [16] source even says 0.4 "is too high (0.4 instead of 0.26)". WP1 samples it rather than fixing it.

## 4. Time scales

| Process | Value | Source |
| --- | --- | --- |
| Optical cycle at 1.2–1.3 µm | 4.0–4.3 fs | — |
| Momentum relaxation 1/Γ | 5–10 fs | [9, 15, 23] |
| e–e thermalisation | 10–30 fs near E_F + ħω/2; [9] estimates 4.2 fs; [60] assumes ≈30 fs | [14, 15, 9, 60] |
| e–ph collision time at E_F | ≈16 fs (number-conserving form) | [14], code |
| Electron cooling / response recovery | 360 fs [1]; 170–400 fs [6, 60]; 62→170 fs rising with intensity [21]; 10–30 fs (transmission) vs hundreds of fs (reflection) at 2 TW/cm², 6 fs pump [20] | |
| Auger-like lag | 120 fs (assigned, not fitted) | [22] |
| Inferred reflectivity rise under 225 fs pump | 1–10 fs "equally good", no error bar | [19] |
| e–e thermalisation under a 5 fs pump (4 nm sphere, density matrix) | ≈12 fs at "1 TW/cm²" (T̃_e ≈ 4500 K), rate ∝ \|E₀\|²; ~1 fs projected at 10 TW/cm². enzkin B1 gives 20 fs at 4500 K, 12 fs at ≈6500 K | [62] |
| e–ph rate at E_F, 300 K | 0.059 fs⁻¹ with D = 17.2 eV [14]; ≈0.01 fs⁻¹ in [62] (0.015 in enzkin with its D = 8.6 eV) | [14, 62], code |

Open conflict: [15] predicts *faster* T_e decay at higher intensity; [21] measures *slower* relaxation with intensity; [20] reports recovery ten times faster at ten times the intensity. [15] reconciles its own prediction with slower *reflection* recovery through the ENZ shift exceeding the resonance width — an electromagnetic effect, again H3.

## 5. Mechanisms for the high-intensity anomalies

| | [22] Auger-like (M2) | [23] two-photon absorption with state filling (M3) | [24] extreme temperature (M0) |
| --- | --- | --- | --- |
| Published model | ṅ_c = η[n₀(T) − n_c], U̇_c = U̇₀ − E_g ṅ_c; η = 1/120 fs⁻¹ assigned; E_g = 3.9 eV; energy input prescribed, no lattice equation | three lumped levels; rate β I² N₁(1 − N₂/N₂max) plus one-photon σ I N₂; fundamental gap 2.67 eV | theory only: fixed density, ω_p(T_e) falls while Γ = Γ₀[1 + (T_e/T₀)^m] rises |
| Driver | absorbed energy (hot-electron population), first-order lag | instantaneous I², gated by depletion of the band bottom | T_e alone |
| Changes n_c? | yes | yes (≈+65 % in the fit) | no (transmittance flips; Re n does not without a density change) |
| Fixed fluence, shorter pulse | fluence-set below ≈120 fs (inference) | stronger; no reversal at 100 fs in the authors' simulation | weakly dependent |
| Data | 310 nm ITO, 1300 nm, 136 fs, 4.6–153 mJ/cm²; anomalies above 55 mJ/cm² | 270 nm ITO, 800 nm, 6.4 fs, 2.9–11 TW/cm²; "4.9 TW/cm²" is a fit crossover, full reversal at ≈8–10 TW/cm² | 10 nm AZO cavity, 10 fs, 1425 nm (simulated) |

Three consequences for the project:

- **Two-photon absorption is energetically closed at the plan's wavelengths.** With E_g = 2.67 eV it needs λ < 929 nm; two 1300 nm photons carry 1.91 eV. The M2-versus-M3 test can only be run on 800 nm-pumped data ([20, 23]), i.e. in the few-cycle regime the plan defers to Tier 2. [22] and [23] also assume different gaps (3.9 vs 2.67 eV). **[TPA threshold verified in code]**
- **A fixed-density thermal model can already flip the sign of the transmittance** [24]. Hypothesis H2 ("sign reversal requires a channel that changes carrier density") must be stated for ω_p or Re ε, not for ΔT/T. M0 also needs the damping law with (m, T₀) free; the plan lists only "coupling constants".
- **Lag is a weak discriminator.** [23]'s reversal also peaks ≈150 fs after a 6.4 fs pump; [22]'s η was set to the observed delay. The energetic argument and the pulse-duration scaling at fixed fluence are the usable tests.

[25] (CdO, 4 µm) attributes sub-cycle time refraction to two-photon interband absorption, Δn ∝ −∫I²dt; [10]'s "invertible nonlinearity" is the sample-dependent direction of the ENZ shift under one-photon interband pumping, not an intensity-driven reversal.

## 6. The electromagnetic side

**Constitutive equation.** The papers give two forms, and only two:

$$\text{(a)}\ \ \dot J+\gamma J=\varepsilon_0\,\omega_p^2(t)\,E\qquad\qquad \text{(b)}\ \ \dot J+\Big(\gamma-\frac{d\ln\omega_p^2}{dt}\Big)J=\varepsilon_0\,\omega_p^2(t)\,E$$

(a) is used by the FDTD two-temperature baseline [9], by the hydrodynamic models [11, 12] (which add −(ṅ/n)Ṗ for a density change), and is the oscillator-density model of [30] (D, B, P, Ṗ continuous across a temporal step). (b) is the main form of [31] (J/ω_p² continuous, J jumps); at fixed density it equals the plan's WP3 law with +(ṁ*/m*)J. Both give the same new frequency after a step but different time-reflected and transmitted amplitudes, and only (a) has a zero-frequency mode. **None of [26, 27, 30, 31] identifies (a) with a density change and (b) with a mass change — that assignment is the project's.** It is physically reasonable (new carriers arrive without drift; heating at conserved momentum lowers the velocity), but which form heated ITO obeys depends on whether the heating collisions conserve momentum: e–e scattering does, assisted absorption does not. It has to be derived from the l = 1 moment of the kinetic equation in WP3; when benchmarking against [9] the ṁ* term must be off.

**What the film does to the measurement (H3).** [17]: the shift measures Δn/Δt during transit. [18]: interface phase terms reverse the sign of the shift at high angle; its quasi-static transfer-matrix model with ω_p(t) is exactly the adiabatic baseline for the H3 test. [19]: α = 1/2 fs⁻¹ is the rise of |r|² of a saturated ITO/Au mirror from a one-parameter, non-dispersive r(t). [20]: thickness-dependent red shift. [62] finds its single-nanosphere theory at odds with [20]'s 40 fs transmission decay and size of change, and names field and temperature gradients and transport in the 300 nm film as one possible cause (asserted), next to interband Kerr and multiphoton terms. [60]: in a non-collinear geometry the pump front sweeps the spot in (w/v) sin θ ≈ 0.35–1 ps for a 600 µm pump — an effect a 1D solver cannot represent. [12]: the 35–50 nm nonlocal blueshift was simulated with m₀* = 0.033 mₑ, n ~ 10²⁰ cm⁻³; with realistic ITO values the pressure coefficient is ≈30× smaller, so it is not the expected correction for WP3.

**Energy cost.** Δε ~ 1 in a TCO needs >10 J/cm³ absorbed ([27] Eq. 5; ≈0.6 eV per electron by [5]'s formula); the 2ω time-crystal term is ≈100× weaker than the slow grating term [27]; a photonic time crystal needs Δn = 0.1–1 within one cycle [28].

## 7. Experiments as datasets

| Ref. | Sample | Pump | Observable | Model used by authors | Raw data |
| --- | --- | --- | --- | --- | --- |
| [1] | 310 nm ITO/glass | 150 fs, ≤250 GW/cm², 1240 nm | R, T vs intensity at 30°; n₂, β by z-scan vs angle/wavelength | two-temperature | not read (paywalled) |
| [17] | 620 nm (2×310) | 1235 nm, 54–485 GW/cm² | probe frequency shift up to 14.9 THz | split-step NLSE, fitted phase | on request |
| [18] | ITO films (arXiv v1) | — | shift vs angle, TM/TE | quasi-static TMM with ω_p(t) | none stated |
| [19] | 40 nm ITO/glass + 100 nm Au | 225 fs field-FWHM (159 fs intensity), ≈124 GW/cm², probe 60°, pumps 54°/66° | time-diffraction fringes | fitted r(t) | on request |
| [20] | ITO, 800 nm ≈6 fs pump, ≈2 TW/cm² | — | red/blue shifts, recovery | FDTD (SI) | none stated |
| [21] | ITO and AZO films | ≈60 fs, 20–200 GW/cm², 0–30° | TE-probe transmission decay τ_r | hydrodynamic + TTM form | in paper/SI |
| [22] | 310 nm ITO/glass | 1300 nm, 136 fs, 4.6–153 mJ/cm² | complex r(t), t(t) by FROG retrieval | (n_c, U_c) rate model | on request |
| [23] | 270 nm ITO/glass | 800 nm, 6.4 fs, ≤11 TW/cm² | ΔR/R, ΔT/T at 1225 nm | three-level rate model + TMM | on request (MATLAB code too) |
| [24] | (simulation) | — | — | TTM, Python code on Figshare 10.6084/m9.figshare.33283908 | public |
| [62] | (simulation) 4 nm ITO sphere | 175 THz, 5 fs, E₀ = 0.6–2.4 V/nm | ε(160 THz, τ_d), f(𝓔,t), e–e and e–ph rates | density matrix + thermal model | on request |
| Herbst et al. 2026 (not in library, cited by [62]) | ITO nanocrystals | two-cycle, ≈10.7 fs, near their ENZ (≈2 µm) | transmitted field by fieldoscopy (90 as resolution) | — | not checked |
| [2] | 900 nm AZO | 785 nm, 100 fs | n₂ | χ⁽³⁾ | open DOI 10.17861/d82f29dc-5c47-4c5a-b893-c1ea93ab5224 |

Corrections to the plan's anomaly table: [21] reports relaxation-time ratios of 2.7 (ITO) and 2.9 (AZO), i.e. +171 % and +185 % (">200 %" holds only relative to the extrapolated zero-intensity value); the "40 %" ENZ-mode figure is AZO (ITO ≈ +23 %). [19]'s 225 fs is a field-envelope FWHM (159 fs in intensity). [15] models the earlier saturable-mirror experiment (pump 60°, probe 65°), not the double-slit experiment [19]. [24] is a theory paper, not an observation. The plan's film range "20–310 nm" excludes [17]'s 620 nm stack.

## 8. Machine learning: what exists and what it implies for WP2

| Ref. | Learns | Constraints | Resolution transfer | In-solver stability | Speed-up |
| --- | --- | --- | --- | --- | --- |
| [46] RelaxNet | Q[f] as ν(E − f), input M − f, fixed-grid MLP | conservation soft (drifts of 10⁻⁴–3 % visible in its figures); positivity and Maxwellian fixed point architectural | none | shock tests with on-the-fly data | not quantified |
| [47] | Q[f] for *linear* radiative transfer, DeepONet | **mass conserved exactly by construction** (10⁻⁷–10⁻⁸); entropy left open | none | — | — |
| [48] FourierSpecNet | Q[f], exactly bilinear low-rank spectral form, 3072 parameters, RK3 outside | none enforced | train 64², test 16²–128², error unchanged | 200–500 steps | ×0.8 (N=16), ×2.7 (64), ×70 (512) vs fast spectral |
| [49] | **flow map**, standard FNO | soft penalty (energy error 10⁻³–10⁻¹ vs 10⁻⁶ spectral) | asserted, untested | 100-step BGK rollouts | ×2 (2V), ×1231 (3V endpoint) |
| [50] | generalized collision kernel from MD | conservation, H-theorem, positivity by construction | — | — | — |
| [51] | ∂ₓq closure, FNO, trained through 100-step rollouts | none | none | 5000 steps; **offline-trained version fails**; fails just outside trained amplitudes | not reported |
| [52] | ∂ₓq from T_e(x), FNO | none stated | train 6× coarser in x, 10–20× in t: no error growth | 3000 implicit steps | ×40 vs SNB reduced model (not vs PIC) |
| [53, 54] | density flow map (TDDFT) | soft particle number | — | drifts; 10 fs stable with laser history | ×8 |

Implications:

- **The novelty claim survives, reworded.** No paper in the group mentions Pauli blocking, Fermi–Dirac statistics, Uehling–Uhlenbeck or electron–phonon coupling. But a hard conservation constraint exists ([47], one invariant, linear operator) and an exact equilibrium fixed point exists ([46]). What is new: two invariants (number and energy) with density-of-states weights, a nonlinear Pauli-blocked operator (cubic in f, which [48]'s bilinear ansatz cannot represent), a Fermi–Dirac fixed point *with* exact conservation, an entropy constraint, and a resolution-transferable architecture. The plan's "[46, 47, 49] use soft or architectural constraints" should become "[46] enforces conservation by a loss penalty; [47] enforces one invariant exactly for a linear operator".
- **[49] is not a collision-operator learner** (the plan cites it as one), and the convolution structure of spectral Boltzmann solvers [48] is over Fourier indices of a v − v* kernel, not a consequence of an energy delta function. The energy-convolution structure used for B0 here is a separate fact.
- **"Train at 256, evaluate at 512 with <2× error growth" has no support in [42, 43]**; the only cross-resolution numbers are in [45] Table 1, where it holds only if the training grid resolves the solution. [43] Sec. 4.4 states that FFT layers cover Neumann boundaries via the cosine transform — the citable basis for the even-extension layer.
- **Train through the solver.** [51]: offline training fails on deployment; [56]: 0.021 single-call error becomes 0.073 after one composition. Unrolled training in the plan is necessary, not optional.
- **Keep Maxwell conventional.** Learned field solvers reach 12–39 % error [55] and accumulate 5–9 % in 160 fs [56]; [57] reports NMSE "as low as 0.24".
- **Uniform vs clustered grid.** The FFT form of B0 and the exact index identity i₁+i₂ = i₃+i₄ need a uniform energy grid; the plan's WP2 grid is clustered near E_F. Either the network works on the uniform grid, or the mapping layer must carry the conservation weights.
- **Cost claims.** [39] and [41] give no usable speed-up for a kinetic closure: [39] has no timings; [41]'s DMD (10–16×) extrapolates one autonomous trajectory and cannot serve as a closure.

**Baselines.** A relaxation-time model with τ(E) on a single distribution violates number and energy conservation ([40] notes only constant τ conserves; [39] quantifies ≈10⁻³ electrons per atom and ≈2 kK); report it with a conserving projection or use the two-subsystem form of [39]. [38]'s eTTM conserves only at the per-mille level and was validated against Boltzmann on a single snapshot.

## 9. Inference and metrology

[61]: use simulator gradients when available; otherwise likelihood or likelihood-ratio surrogates. WP4's deterministic, differentiable forward model with an explicit noise model has a tractable likelihood, so Hamiltonian Monte Carlo is ordinary likelihood-based inference; neural posterior estimation is justified mainly by amortisation across identifiability trials. [61] does not treat Bayes factors or identifiability.

[58]: a network trained on simulated FROG traces was worse than PCGPA/ptychography on measured data until refined through a differentiable forward model. [60] is a beam-deflection intensity cross-correlation in 240 nm AZO, not a FROG gate; phase "is implicitly lost"; fitted τ_p = 95 fs, τ_R = 400 fs, instantaneous rise assumed. [59] (ITO, ENZ ≈1510 nm, non-collinear four-wave-mixing gate) shows the group's laboratory capability. **None of [58–60] treats double-blind retrieval or a complex gate; the closest is [22]'s retrieval of r(t).** The plan's "related ENZ pulse metrology [58–60]" should be narrowed accordingly.

## 10. Where the literature contradicts the project plan

This list refers to revision 1 of the plan (2026-11-16). Revision 2 (2026-11-17) incorporates every item; it is kept here as the record of why the plan changed. Items 12–13 come from [62] and are incorporated in revision 2.1 (2026-11-18).

1. Two-photon absorption cannot operate at 1.1–1.6 µm for E_g = 2.67 eV → M3 is testable only on 800 nm data (Tier 2 regime).
2. H2 must be phrased for Re ε or ω_p; transmittance sign reversal does not need a density change [24].
3. M0 is a two-temperature model (not "extended") and needs a temperature-dependent damping law to produce any reversal.
4. The density-versus-mass reading of the two current equations is the project's extrapolation, not a result of [30, 31].
5. [14] contains no Umklapp treatment; its absorption term has an unspecified matrix element; the full e–e integral is in [14], not [15].
6. [49] learns a flow map; [46, 47] are not Fourier operators; [47] already has a hard constraint.
7. The 256→512 criterion is unsupported by the cited FNO papers.
8. [24] is theory; [21]'s ">200 %" and "40 %" need the qualifiers above; [19]'s pump is 159 fs in intensity; [15] models a different experiment from [19].
9. [58–60] are not double-blind FROG precedents.
10. [36]'s tiny non-thermal power fraction is a CW result and must not be used for femtosecond pumping.
11. Citation fixes: [16] has three authors (Un, Sarkar, Sivan); [22] was retitled and published (Optica 13, 1943); [10], [13], [29], [39], [40], [46], [48], [50]–[52], [54]–[56], [59], [60] are now published (see `references/README.md`); [47]'s authors are Lee, Schotthöfer, Xiao, Krumscheid, Frank; [50]'s second author is J. W. Burby.
12. Revision 2 says sub-8 fs data need the density matrix of [16] and that it is too expensive in a film. [62] shows that for a 5 fs pump a thermal model driven by the instantaneous absorbed power reproduces the density-matrix permittivity, and [16]'s finite-size absorption does not transfer to films anyway: the first few-cycle step is a field-resolved J·E drive of the WP3 cells (Tier 2a).
13. Revision 2 fixes the e–ph coupling at the [14] value; the authors' own later set [62] has a four times weaker coupling and an energy-dependent impurity damping. D (or G) and the impurity model must be inferred.

## 11. What must come from the authors or from papers not in the library

- [15]: ε∞, the e–impurity rate, the gold permittivity and incidence medium, the conduction-band width actually used, whether the screening wavevector was recomputed from f, grid and time step; and why its plotted ω_p is 8 % above the value implied by its printed parameters. Its sample is that of Tirole et al., Phys. Rev. Applied 2022 (not in the library).
- [22]: Supplementary Information (n₀(T) formula, Kane parameters, valence density of states, U_max vs fluence); raw r(t), t(t).
- [23]: the MATLAB rate-equation code (units of β and σ are inconsistent as printed); which pulse duration the model used.
- [24]: Supplementary Sections 1–13; the Figshare code is public.
- [1]: the paper and its supplement (two-temperature parameters, R/T data vs angle).
- [62]: T̃_e(t) for the runs of Figs. 4–5; ε∞ and the dielectric constant used for screening (ε∞ in S3–S4, ε_b elsewhere); the band top and grid; whether the ½ of (S2) or the rate (S7) without it was integrated; which intensity convention ("1 TW/cm²" for E₀ = 1.8 V/nm); how D = 8.6 eV and n_imp were fitted and whether they apply to films such as [15]'s.
- Valence-band mass and density of states for ITO (needed by M2 and M3; printed in neither [22] nor [23]).

---

# Part II. Per-paper digests


## ENZ nonlinearity: foundations and models
### [1] Alam et al., "Large optical nonlinearity of indium tin oxide in its epsilon-near-zero region" — Science 352, 795 (2016)

**Source read:** abstract only (verbatim OpenAlex copy in `references/text/01_alam2016.abstract.txt`). Main text and Supplementary Materials NOT READ (paywalled, no PDF). Everything beyond the abstract is reconstructed from how papers in this library describe it and from background knowledge; verify before use. Each item below carries a provenance tag: **[abs]** abstract; **[lib: X]** stated about this paper by library paper X (Baxter = [9], Zhou = [17], Li = [10], Un = [15], Galiffi = [26], Ball = [60]); **[bg]** my recollection of the paper, not checked against any source here.
**Role in this project:** The calibration dataset: step 2 of the WP3 validation ladder (the only step at which parameters may be tuned) and the WP4 dataset that fixes the intraband response of M0/M1.

**What it does.** It measures the intensity-dependent refractive index of an ITO film pumped at its zero-permittivity wavelength with femtosecond pulses. The film acquires "an ultrafast and large intensity-dependent refractive index in the region of the spectrum where the real part of its permittivity vanishes": the real index changes by 0.72 ± 0.025, which is 170% of the linear index, and the change is reversible with a recovery time of about 360 fs [abs]. The response is quantified as an effective Kerr coefficient from z-scan measurements as a function of pump wavelength and angle of incidence, as intensity-dependent reflectance/transmittance, and as a degenerate pump–probe transient [lib: Baxter].

**Model / method.** No equation could be copied (paper not read). What the library attributes to it:

- Linear response: a Drude fit of the film [lib: Baxter, written in Baxter's Eq. 1 form]

$$\varepsilon(\omega)=\varepsilon_\infty-\frac{\omega_p^2}{\omega^2+i\gamma\omega}$$

- Nonlinear response reported as Kerr-like, $n=n_0+n_2I$, with $\mathrm{Re}(n_2)$ and $\beta=(4\pi/\lambda)\,\mathrm{Im}(n_2)$ "extracted via z-scan measurements" [lib: Baxter]. These are effective, pulse- and geometry-specific coefficients.
- Dynamics described with a two-temperature model (TTM): pump energy heats the conduction electrons, the hot distribution in the non-parabolic band has a red-shifted plasma frequency, and the index recovers through electron–phonon coupling [lib: Zhou, citing this paper and one other for the TTM]. The approach "uses fitting parameters, such as electron-phonon coupling strength G determined from time-resolved index change measurements" [lib: Li, said of this paper and two others together] and involves a temperature-dependent electron heat capacity [lib: Li]. The electron-dynamics model is separate from the optical-propagation model [lib: Baxter].
- [bg] The Supplementary Materials contain the TTM (with a non-thermal energy reservoir feeding electrons and lattice, the same structure as Baxter's Eq. 4) and a Kane-band $\omega_p(T_e)$. The equations and all TTM parameter values (coupling constant, heat capacities, relaxation times, electron density, band-edge mass, non-parabolicity, Fermi energy) are **not available** from any source in this library.

**Numbers.**

| Quantity | Value | Provenance |
| --- | --- | --- |
| Film | ITO, 310 nm, on glass | lib: Baxter Sec. 3 |
| Drude fit (room temperature) | $\omega_p=2\pi\times473$ THz, $\gamma=0.0468\,\omega_p$, $\varepsilon_\infty=3.8055$ | lib: Baxter Sec. 3, "reported in Ref. [14]" |
| Zero-permittivity wavelength | 1240 nm (Baxter quotes $\lambda_{ZC}=1237$ nm from $2\pi\sqrt{\varepsilon_\infty}c/\omega_p$) | lib: Baxter Secs. 2, 3 |
| Derived from the Drude fit (my arithmetic) | $\hbar\omega_p=1.96$ eV; $\gamma=1.39\times10^{14}$ s$^{-1}$ ($1/\gamma=7.2$ fs); $\mathrm{Re}\,\varepsilon=0$ at 1242 nm when damping is kept; at 1240 nm $\varepsilon=0.01+0.35i$, $n=0.42+0.41i$ | my arithmetic; $n_0=0.42$ agrees with 0.72 being 170% |
| Pulse duration | 150 fs | lib: Baxter Sec. 4 |
| z-scan result | maximum $\mathrm{Re}(n_2)$ and minimum $\beta$ at pump wavelength 1240 nm and 60° from the surface normal: $n_2\approx0.11$ cm²/GW, $\beta\approx-7500$ cm/GW (Baxter prints the unit of $\beta$ once as cm²/GW; typo) | lib: Baxter Sec. 3 |
| Intensity scan | 1240 nm pump at 30°, 0 to 250 GW/cm² (incident, from the air side as modelled by Baxter): $T$ from ≈0.12 to ≈0.4, $R$ from ≈0.2 to ≈0.04; absorptance also reported | lib: Baxter Sec. 3 |
| Index change | $\Delta n=0.72\pm0.025$ = 170% of linear index | abs |
| Related figure of merit | $\Delta n\approx0.72$ with $\alpha=27{,}362$ cm$^{-1}$ | lib: Ball (cites this paper and one other) |
| Transient (degenerate pump–probe, transmission) | rise time < 200 fs for the 150 fs pulse; fall time ≈ 360 fs; both defined 10%–90% of the maximum change | lib: Baxter Sec. 4; recovery also abs, Galiffi |
| TTM parameters, $N$, $m^\ast$, $C$, $E_F$, substrate thickness, beam size, repetition rate | not available | — |

Unverified recollections [bg] (do not copy without opening the paper): commercially sourced film on 1.1 mm float glass; TM-polarized pump; z-scan over roughly 950–1500 nm and 0–60°; enhancement of $n_2$ and $\beta$ at 1240 nm relative to 970 nm by factors 43 and 53 at normal incidence and 1837 and 2377 at 60°; $\beta\approx-7.3\times10^3$ cm/GW at 60°; $\Delta n=0.72$ obtained by transfer-matrix inversion of the 30° $R$, $T$ data near 150 GW/cm², above which the change saturates.

**Key results.**
- Order-unity, reversible index change in a sub-wavelength film at ENZ: $\Delta n=0.72\pm0.025$, recovery ≈ 360 fs (measured) [abs].
- The effective $n_2$ peaks at the zero-permittivity wavelength and grows strongly with angle of incidence up to 60° for the polarization used (measured) [lib: Baxter].
- Transmittance roughly triples and reflectance falls about five-fold between 0 and 250 GW/cm² at 30° (measured) [lib: Baxter]; the film goes from metallic to dielectric (Baxter's interpretation).
- Rise is pulse-limited (< 200 fs with 150 fs pulses); recovery ≈ 360 fs (measured) [lib: Baxter].

**Assumptions and limitations.**
- $n_2$ and $\beta$ are effective coefficients for one film, one pulse duration and one geometry. Baxter calls the Kerr description "a rough approximation" (higher orders, non-instantaneous); Un states that such fitted values "strongly depend on the pulse duration".
- z-scan averages over the transverse beam profile, the pulse and the film depth (my reading).
- 150 fs pulses cannot resolve a rise faster than about the pulse; nothing here constrains few-fs dynamics (my reading).
- The TTM assumes thermalized electrons; parameters were fitted to the transient [lib: Li].

**Use for this project.**
- Fit targets, in order of directness: $R(I)$, $T(I)$, $A(I)$ at 1240 nm, 30°; the pump–probe transmission transient; the z-scan maps $n_2(\lambda,\theta)$, $\beta(\lambda,\theta)$. Compare raw $R$, $T$ against incident intensity; do not fit the derived $n_2$.
- Raw data: no data-availability statement is known to me (paper not read). Plan to request from the authors or digitise.
- Adopt the Drude fit above as the linear starting point for this sample; step 1 of the ladder should reproduce the low-intensity $T\approx0.12$, $R\approx0.2$ at 30°.
- Baseline to beat (Baxter's FDTD+TTM, with a fitted damping slope and fitted electron–phonon prefactor): $n_2\approx0.16$ cm²/GW and $\beta\approx-13000$ cm/GW against 0.11 and −7500; $T$ 0.13→0.4, $R$ 0.19→0.02; rise 170 fs, fall 385 fs. Un notes that fit needed damping rising linearly with $T_e$ and an electron–phonon coupling reduced by a factor of 8.
- Possible mismatch with the plan: the WP4 table lists "reflection and transmission against intensity and angle". As described by Baxter, the intensity scan of $R$, $T$ is at one angle (30°); the angle and wavelength dependence is in the z-scan coefficients. Check the paper before designing the likelihood.
- What limits the speed, per the same group later [lib: Zhou]: rise "limited by the thermalization time of the conduction band electrons owing to electron-electron scattering" and by the pump envelope; recovery by electron–phonon coupling. With 150 fs pulses this dataset cannot separate the two rise contributions, so it calibrates magnitude and recovery, not rise (H1, H3).
- Open: the paper's own material and TTM parameters must be read from the Supplementary Materials before any claim of "same parameters as Alam".

---
### [2] Caspani et al., "Enhanced nonlinear refractive index in epsilon-near-zero materials" — Phys. Rev. Lett. 116, 233901 (2016)

**Source read:** arXiv:1603.03581v2 full text + TeX source (accepted manuscript); figures checked in the PDF. "Read from figure" marks approximate values.
**Role in this project:** Background for order-unity index change in a second oxide (AZO), the origin of the Kerr-type $n_2$ description the project replaces, and an openly available secondary dataset.

**What it does.** A 900 nm Al-doped ZnO film is pumped at 785 nm and probed from 1150 to 1550 nm. The pumped permittivity is retrieved from $\Delta R$, $\Delta T$ by an inverse transfer-matrix method and converted to $\chi^{(3)}$ and $n_2$. $\mathrm{Re}\,n_2$ peaks at the ENZ wavelength and nonlinear absorption changes sign there; the paper attributes this to the linear index alone, with a nearly constant $\chi^{(3)}$. No microscopic mechanism or time constant is given.

**Model / method.** Non-degenerate Kerr index (Eqs. 1–3), with $n_{r,i}$ the linear index at the probe, $n_r^{\rm pump}$ the real index at the pump, $D=n_r^{\rm pump}(n_r^2+n_i^2)$, $\beta_2=4\pi n_{2i}/\lambda$, $\delta n=n_2I$:

$$n_2=\frac{3}{2\varepsilon_0 c}\frac{\chi^{(3)}}{n_r^{\rm pump}(n_r+in_i)},\qquad n_{2r}=\frac{3}{2\varepsilon_0 c}\frac{n_r\chi^{(3)}_r+n_i\chi^{(3)}_i}{D},\qquad n_{2i}=\frac{3}{2\varepsilon_0 c}\frac{n_r\chi^{(3)}_i-n_i\chi^{(3)}_r}{D}$$

$\chi^{(3)}$ from the slope of the pumped permittivity $\varepsilon_{nl}$ against pump intensity $I_p$, where that dependence is linear (Eq. 4):

$$\chi^{(3)}(\omega_p,\omega_{\rm probe})=\frac{n_r^{\rm pump}\varepsilon_0c}{3}\frac{\partial\varepsilon_{nl}(\omega_{\rm probe},I_p)}{\partial I_p}$$

Eqs. 5, 6 are Eqs. 2, 3 with fitted constants $A$, $B$ for $\chi^{(3)}_r$, $\chi^{(3)}_i$.

**Numbers.**
- Sample: oxygen-deprived AZO, 900 nm, pulsed laser deposition. Substrate, electron density, mass, damping: not stated.
- ENZ wavelength "∼1300 nm" (dashed line at ≈1.31 µm in Figs. 1c, 5). Unpumped $\varepsilon\approx0.14+0.17i$ at 1258 nm (Fig. 4a intercepts, read from figure). $n_r=0.09$ at 1390 nm.
- Pump: 785 nm, 100 fs, 100 Hz, normal incidence, waist 125 µm; $I_p$ = 435, 870, 1300 GW/cm² (not stated whether incident or in-film).
- Probe: 100 fs, 1150–1550 nm, < 10° incidence, crossed polarization, waist 45 µm, bandwidth 15 nm, below 1 GW/cm². Delay "optimised to maximise the nonlinear effect"; no delay scan shown.
- $A=4\times10^{-20}$, $B=1\times10^{-20}$ m²/V² (after Eq. 6).
- $n_{2r}$: peak ≈ 3.5×10⁻¹³ cm²/W near 1.31 µm, six-fold above its value at 1152 nm (≈0.6×10⁻¹³; Fig. 5a, read from figure). $\beta_2$: ≈0 near 1.19 µm, minimum ≈ −3.2×10⁻⁸ cm/W at 1.34–1.39 µm (Fig. 5b, read from figure).
- $\delta n_r=0.4$ at 1390 nm at 1300 GW/cm². Probe carrier shift up to 17.5 ± 1.6 nm (Fig. 5d).
- No damage up to 2 TW/cm² (785 nm, 100 Hz).

**Key results.**
- Six-fold enhancement of $n_{2r}$ at ENZ (measured).
- Nonlinear absorption turns negative (saturable) at the ENZ wavelength (measured; predicted from the linear index with $\chi^{(3)}_r\sim\chi^{(3)}_i$ constant).
- Index change comparable to the linear index, so the power-series expansion of the polarization "breaks down" (asserted).

**Assumptions and limitations.**
- Phenomenological, non-dispersive $\chi^{(3)}$.
- The retrieval treats the 900 nm film as homogeneous although pump absorption varies with depth (my reading).
- One delay and 100 fs pulses: values are specific to this pulse duration (my reading).

**Use for this project.**
- Data-availability statement: "All data relevant to this work may be obtained at DOI: 10.17861/d82f29dc-5c47-4c5a-b893-c1ea93ab5224." It is the only data repository stated in the papers of this group that could be read.
- Usable as an out-of-sample test if AZO parameters are added; Secondo et al. [6] fit these data with a thermal effective-mass model and $\tau_{el}=170$ fs.
- Pitfalls: [6] quotes this sample's ENZ point as 1330 nm; $n_2$ is not a material constant, and Eq. 4 holds only while $\varepsilon_{nl}(I_p)$ is linear.

---
### [3] Reshef et al., "Nonlinear optical effects in epsilon-near-zero media" — Nat. Rev. Mater. 4, 535 (2019)

**Source read:** NOT READ — paywalled, no PDF and no abstract could be retrieved. Reconstructed from citing papers in this library and background knowledge; verify before use.
**Role in this project:** Context only: the general review of ENZ nonlinear optics by the group that produced [1], cited in the plan as background.

**What it does.** (background, not checked against the paper) A review of why a vanishing permittivity enhances nonlinear optical effects, which homogeneous and structured media provide it, and what had been demonstrated by 2019: intensity-dependent index, harmonic generation, wave mixing, frequency translation, and nanoantennas coupled to ENZ films.

**Model / method.** No equation could be copied. Two relations usually associated with this literature (background, not from this paper as read):
- $\Delta n=\Delta\varepsilon/2\sqrt{\varepsilon}$, which grows as $\varepsilon\to0$. Li et al. [10] state the same relation in their introduction.
- When $\Delta n$ is comparable to $n$, the index must be computed non-perturbatively from the field-dependent permittivity, not as $n_0+n_2I$ (same group, Opt. Lett. 42, 3225 (2017)).

**Numbers.** None attested in this library.

**Key results.** As used by library papers:
- Cited for the statement that in a highly nonlinear low-index medium $\Delta n/(n_1+\Delta n)$ can approach unity (Zhou et al. [17]).
- Cited as having noted that doped ZnO slightly outperforms ITO (Secondo et al. [6], Sec. 5).
- Otherwise cited as a general review ([5], [9], [10]).

**Use for this project.**
- Introductory citation only. Take no parameter or equation from this digest.
- Open the paper if comparative $n_2$ values or its statements on mechanisms and time constants are needed; none are known here.

---
### [4] Kinsey et al., "Near-zero-index materials for photonics" — Nat. Rev. Mater. 4, 742 (2019)

**Source read:** NOT READ — paywalled, no PDF and no abstract could be retrieved. Reconstructed from citing papers in this library and background knowledge; verify before use.
**Role in this project:** Context only: the materials-oriented review of near-zero-index media, cited in the plan as background.

**What it does.** (background, not checked against the paper) A review of near-zero-index materials (epsilon-near-zero, mu-near-zero, and both): homogeneous materials such as transparent conducting oxides and nitrides, structured media, their linear effects (field enhancement, slow light, wavelength stretching), their nonlinear and time-varying uses, and optical loss as the practical limit.

**Model / method.** No equation could be copied.

**Numbers.** None attested in this library.

**Key results.** As used by library papers:
- Cited for the statement that the ENZ condition may be met "in bulk materials near resonances or through free carriers, as well as in nanostructured materials as an effective property by mixing both metals and dielectrics" (Secondo et al. [6], Sec. 1).
- Cited for materials whose refractive index can be drastically modified by an optical pump (Li et al. [10]).
- Otherwise cited as a general review ([5], [16], [19], [26]).
- (background, unverified) The review separates two ways of pumping a conducting oxide: intraband (hot electrons, heavier average mass, red shift of the plasma frequency) and interband (added carriers, blue shift), with different recovery channels.

**Use for this project.**
- Introductory citation only. Take no parameter or equation from this digest.
- If the intraband/interband sign distinction is quoted in support of H2, cite it from a paper that was read (for example [10]) or open this review first.

---
### [5] Khurgin et al., "Fast and slow nonlinearities in epsilon-near-zero materials" — Laser Photonics Rev. 15, 2000291 (2021)

**Source read:** arXiv:2007.05613v1, full text (33-page Word-generated PDF; no TeX source, so equations were transcribed from the PDF pages and checked numerically where noted). The arXiv title is "Fast and slow nonlinearities in ENZ materials". The published version was not available and may differ. Equation numbers are those of arXiv v1; its text cites "Figure 5" for the plot labelled Figure 6.
**Role in this project:** Analytic form of the "established picture" (plan refs [5–8]): closed-form Kane-band limits that WP1's $\omega_p^2[f]$ and models M0/M1 must reproduce at low excitation, and the two time constants (momentum and energy relaxation) behind the project's central question.

**What it does.** It asks how the "slow" hot-electron nonlinearity of transparent conducting oxides (index switching) and the "fast" instantaneous one (harmonic generation) are related. A two-level density-matrix model gives the bound-electron (interband) reference; a single non-parabolic band with a momentum relaxation rate gives the free-electron (intraband) case. Both intraband terms come from non-parabolicity, and their ratio is $2\tau_{el}/\tau_s$, the analogue of $2T_1/T_2$. Estimates for ITO: $\chi^{(3)}_{fast}\approx2.5\times10^{-19}$ m²/V² and $\chi^{(3)}_{slow}$ up to $10^{-17}$ m²/V².

**Model / method.** Fields are $E(t)=E_1e^{-i\omega_1t}+E_2e^{-i\omega_2t}+c.c.$; $\Delta\omega$ is the beat frequency, $\bar\omega$ the mean frequency; the process is degenerate four-wave mixing, $\omega_3=\omega_1-\omega_2+\omega_1$.

*Interband reference* (Eq. 8; $N$ density of two-level systems, $\mu_{21}=ez_{21}$, $T_1$ population lifetime, $\gamma=1/T_2$ dephasing rate, $\Phi$ an irrelevant phase), and the fast-part estimate (Eq. 9):

$$\chi^{(3)}\sim\left(1+\frac{2T_1/T_2-1}{1-i\Delta\omega T_1}\right)\frac{N\mu_{21}^4/\varepsilon_0\hbar^3}{\left[(\omega_{21}-\bar\omega)^2+\gamma^2\right]^{3/2}}e^{j\Phi}=\chi^{(3)}_{fast}+\chi^{(3)}_{slow},\qquad \chi^{(3)}_{fast}\sim\chi^{(1)}\frac{e^2z_{21}^2}{\left[(\omega_{21}-\bar\omega)^2+\gamma^2\right]\hbar^2}$$

*Intraband: carrier motion* (Eqs. 10, 11; $\gamma$ is now the momentum scattering rate):

$$\frac{dk_z}{dt}=-\frac{e}{\hbar}E_z-\gamma(k_z-k_{z0}),\qquad \delta k(t)=\frac{e}{\hbar(i\omega_1-\gamma)}E_1e^{-i\omega_1t}+\frac{e}{\hbar(i\omega_2-\gamma)}E_2e^{-i\omega_2t}+c.c.$$

*Transport mass* (Eq. 14), the angle average of $\hbar^{-1}\partial v_z/\partial k_z$ for an isotropic band, derived in [6]; it differs from $m_c^{-1}=\hbar^{-1}dv/dk$ and $m_{opt}^{-1}=\hbar^{-1}v/k$ and stays finite and positive when $dv/dk\le0$:

$$m_t^{-1}=\hbar^{-1}\left[\frac23\frac{v(k)}{k}+\frac13\frac{dv(k)}{dk}\right]=\frac{1}{3\hbar k^2}\frac{d}{dk}\left(k^2v(k)\right)$$

*Kane two-band dispersion* (Eqs. 15–18; $m^\ast$ band-edge mass, $E_{gap}$ gap). As printed $E(0)=E_{gap}/2$, although the text says energy is measured from the band bottom:

$$E(k)=\frac{E_{gap}}{2}\sqrt{1+\frac{2\hbar^2k^2}{m^\ast E_{gap}}}=\frac{E_{gap}}{2}\sqrt{1+(k/k_0)^2},\qquad k_0=\sqrt{m^\ast E_{gap}/2\hbar^2}$$

$$v_z(k)=v_{sat}\frac{k/k_0}{\sqrt{1+(k/k_0)^2}},\qquad v_{sat}=\frac12\frac{E_{gap}}{\hbar k_0}=\frac{\hbar k_0}{m^\ast},\qquad m_t^{-1}=m^{\ast-1}\frac{2(k/k_0)^2+3}{3\left[1+(k/k_0)^2\right]^{3/2}}$$

*Derivatives of the inverse mass* (Eqs. 19, 20):

$$\alpha=\frac{dm_t^{-1}}{dk}=\frac{m^{\ast-1}}{k_0}\frac{k}{k_0}\frac{2(k/k_0)^2+5}{3\left[1+(k/k_0)^2\right]^{5/2}},\qquad \beta=\frac{d^2m_t^{-1}}{dk^2}=\frac{m^{\ast-1}}{k_0^2}\frac{4(k/k_0)^4+14(k/k_0)^2-5}{3\left[1+(k/k_0)^2\right]^{7/2}}$$

($\alpha$ as printed; differentiating Eq. 18 gives this expression with a minus sign, since the mass grows with $k$. The paper then uses $\alpha_F=|\alpha(k_F)|$.)

*Current and Fermi-sea averages* (Eqs. 22, 23; $N$ carrier density, $k_F$ Fermi wavevector; the $\alpha$ term cancels because $\alpha(k)$ is odd):

$$J(t)\approx-Ne\hbar\left[\langle m_t^{-1}\rangle\delta k(t)+\frac16\langle\beta\rangle\delta k^3(t)\right],\qquad \langle m_t^{-1}\rangle=\frac{3}{k_F^2}\int_0^{k_F}m_t^{-1}(k)k^2dk=\frac{m^{\ast-1}}{\left[1+(k_F/k_0)^2\right]^{1/2}}$$

$$\langle\beta\rangle=\frac{3}{k_F^2}\int_0^{k_F}\beta k^2dk=\frac{m^{\ast-1}}{k_0^2}\left\{\frac{4\sinh^{-1}(k_F/k_0)}{(k_F/k_0)^3}-\frac{4(k_F/k_0)^4+11(k_F/k_0)^2+4}{(k_F/k_0)^2\left[1+(k_F/k_0)^2\right]^{1/2}}\right\}$$

(Both as printed. The prefactor must be $3/k_F^3$ for the dimensions and for the stated result. The $\langle\beta\rangle$ bracket as printed gives −10.9 at $k_F/k_0=1.6$; numerical integration of Eq. 20 gives +0.266, which is reproduced if the last exponent is 5/2 instead of 1/2, and then $f_\beta=0.50$ as the paper states. My check.)

*Linear response* (Eq. 25): $\chi^{(1)}_{fc}(\omega)=-\dfrac{Ne^2\langle m_t^{-1}\rangle}{\varepsilon_0\omega(\omega+i\gamma)}$, $\varepsilon(\omega)=\varepsilon_\infty+\chi^{(1)}_{fc}(\omega)$.

*Fast intraband nonlinearity* (Eqs. 28–30), with $f_\beta(k_F)=k_0^2\langle\beta\rangle\langle m_t\rangle$:

$$\chi^{(3)}_{fast}=\frac{Ne^4\langle\beta\rangle}{2\varepsilon_0\hbar^2(\bar\omega^2+\gamma^2)^{3/2}}\approx\chi^{(1)}_{fc}(\omega)\frac{e^2f_\beta}{2\hbar^2k_0^2(\bar\omega^2+\gamma^2)}\approx\varepsilon_\infty\frac{e^2}{2\hbar^2\omega^2k_0^2}f_\beta\ \text{(near ENZ)}$$

(As printed. The first form is one power of frequency short of the second; the second is dimensionally consistent and reproduces the quoted number.)

*Slow intraband nonlinearity.* Absorbed power and hot-carrier energy density (Eqs. 31, 32; $\tau_{el}$ electron–lattice energy relaxation time):

$$\delta P(t)=J(t)E(t)=\frac{2Ne^2\gamma\langle m_t^{-1}\rangle}{\bar\omega^2+\gamma^2}\left[E_1^2+E_2^2+2E_1E_2\cos[(\omega_1-\omega_2)t]\right],\qquad \frac{dU_{hot}(t)}{dt}=\delta P(t)-\frac{U_{hot}}{\tau_{el}}$$

A fraction $f_{hot}$ of carriers is promoted across the Fermi level; each gains $\delta E=U/f_{hot}N$ and shifts by $\delta k=\delta E/\hbar v_F$ (Eqs. 34–36):

$$\delta m_t^{-1}(t)=\alpha_F\frac{U_{hot}(t)}{f_{hot}N\hbar v_F},\qquad \delta\varepsilon(t)=\delta\chi^{(1)}(t)\approx f_{hot}\frac{\delta m_t^{-1}}{\langle m_t^{-1}\rangle}\chi^{(1)}=\alpha_F\frac{U_{hot}(t)}{\langle m_t^{-1}\rangle N\hbar v_F}\chi^{(1)}$$

Resulting susceptibility (Eqs. 38–42), with $\tau_s=\gamma^{-1}$, $\alpha_F/v_F=f_\alpha(k_F)/\hbar k_0^2$ and $f_\alpha(k_F)=\dfrac{2(k_F/k_0)^2+5}{3[1+(k_F/k_0)^2]^2}$:

$$\chi^{(3)}_{slow}=\frac{\gamma\tau_{el}}{1-i\Delta\omega\tau_{el}}\frac{2\alpha_F}{\hbar v_F}\frac{e^2}{\bar\omega^2+\gamma^2}\chi^{(1)}=\frac{2\tau_{el}/\tau_s}{1-i\Delta\omega\tau_{el}}\chi^{(1)}\frac{e^2f_\alpha}{\hbar^2k_0^2(\bar\omega^2+\gamma^2)}\approx\frac{2\tau_{el}/\tau_s}{1-i\Delta\omega\tau_{el}}\chi^{(3)}_{fast}$$

*Link between the two mechanisms* (Eq. 44, from the oscillator sum rule; $P_{cv}=m_0z_{cv}E_{gap}/\hbar$): $\dfrac{m_0E_{gap}}{2\hbar^2k_0^2}=1+\dfrac{2P_{cv}^2}{m_0E_{gap}}$, hence $k_0^{-1}\sim z_{cv}$.

*Figure of merit* (Eq. 45; $L_a\sim c\chi^{(1)}/\gamma$ absorption length, $n_g$ group index): $FOM=\chi^{(3)}_{slow}L_aB\approx\dfrac{n_ge^2cf_\alpha}{2\pi\hbar^2k_0^2\bar\omega^2}$, $B=1/2\pi\tau_{el}$.

**Numbers.**
- ITO example (Sec. 4): $N=1.1\times10^{21}$ cm⁻³, $k_F=(3\pi^2N)^{1/3}=3.2$ nm⁻¹, $k_0\approx2$ nm⁻¹ ("obtained by fitting", source not stated), $k_F/k_0\approx1.6$, $f_\beta\approx0.5$, $f_\alpha\approx f_\beta/2$ (Fig. 6; Eq. 40 gives 0.27), $\chi^{(1)}_{fc}\sim-\varepsilon_\infty\approx-4$, photon energy 1 eV.
- Not stated: band-edge mass, $E_{gap}$, Drude damping value, ENZ wavelength, any film or pulse parameter.
- $\chi^{(3)}_{fast}\approx2.5\times10^{-19}$ m²/V² (Eq. 30). $\chi^{(3)}_{slow}$ "about two order of magnitude higher", up to $10^{-17}$ m²/V² (Sec. 5).
- Bound electrons: $z_{21}\sim0.1$ nm, detuning ∼1 eV, $\chi^{(1)}\sim10$, giving $\sim10^{-19}$ m²/V² (Sec. 2; printed as "$10^{19}$").
- $\tau_s$: "on the scale of 10 fs or less" for ITO, AZO. $\tau_{el}$: "100's of femtoseconds", $100\ \text{fs}<\tau_{el}<1$ ps. Electron–electron collisions: "fast (100fs or less)", yet "it may take longer than 100 fs to establish thermal equilibrium". Interband $T_1$: nanoseconds.
- Validity: $T_e\ll E_F/k_B\approx15{,}000$ K; Fermi energy "approaching 1 eV".
- No fluence, energy-cost, $n_2$ or $\Delta n$ estimate appears in this version.

**Key results.**
- Fast term: third-order expansion of the band velocity ($\langle\beta\rangle$), no characteristic time; also gives third-harmonic generation (derived).
- Slow term: absorbed energy moves carriers to heavier states ($\alpha_F$); it follows $U_{hot}$, so it decays with $\tau_{el}$ and has 3 dB bandwidth $1/2\pi\tau_{el}$ (derived).
- $\chi^{(3)}_{slow}/\chi^{(3)}_{fast}\approx2\tau_{el}/\tau_s$ at low beat frequency (derived).
- $\delta\varepsilon$ is independent of $f_{hot}$; the authors conclude that computing the electron temperature or the exact distribution "is unnecessary" while $T_e\ll E_F/k_B$ (derived within the linearization).
- Fast intraband and fast interband magnitudes are similar because $k_0^{-1}\sim z_{cv}$ (derived, order of magnitude).
- Intraband loss is unavoidable; the product of nonlinearity, absorption length and bandwidth is about the same for any mechanism at a given wavelength (asserted).

**Assumptions and limitations.**
- Perturbative: third order in field, linear in absorbed energy. No saturation or order-unity $\delta\varepsilon$.
- One isotropic Kane band, zero-temperature Fermi sphere, constant $\gamma$ unaffected by heating, single-exponential energy loss, no lattice heating.
- All promoted carriers are treated as sitting at $k_F$ with slope $\alpha_F$. This needs $m_t^{-1}(E)$ to be linear over $E_F\pm\hbar\omega$, while $\hbar\omega\approx E_F\approx1$ eV (my reading).
- Bulk, local, uniform internal field; no film or propagation.
- Signs are not tracked ($\alpha_F$ is a modulus); treat Eqs. 35–41 as magnitudes. Physically $\langle m_t^{-1}\rangle$ falls and $\mathrm{Re}\,\varepsilon$ rises (my reading).
- "Interband" here means bound-electron nonlinearity of conventional materials, not carrier generation in the oxide.

**Use for this project.**
- $\varepsilon_0\omega_p^2=Ne^2\langle m_t^{-1}\rangle$ with Eq. 14 is the plan's $\omega_p^2[f]$ integrated by parts, for any isotropic $f$ (my derivation). The $T=0$ value $\langle m_t^{-1}\rangle=m^{\ast-1}[1+(k_F/k_0)^2]^{-1/2}$ is a unit test for WP1.
- Mapping to the plan's band, with energy from the band bottom: $\hbar^2k^2/2m^\ast=E(1+CE)$, $C=1/E_{gap}$, $k_0^2=m^\ast/2\hbar^2C$ (my derivation). The paper gives neither $C$ nor $m^\ast$. Taking $C\approx0.42$ eV⁻¹ (the value Baxter et al. [9] use for ITO; not from this paper), $k_0=2$ nm⁻¹ implies $m^\ast\approx0.26\,m_0$, which matches the ITO band fit in [6].
- Low-excitation target: from Eqs. 36, 39, $|\delta\varepsilon|\approx|\chi^{(1)}|f_\alpha\langle m_t\rangle U_{hot}/N\hbar^2k_0^2$ (my rearrangement). To first order the solver's $\delta\varepsilon$ should depend on absorbed energy density only.
- Energy scale (my arithmetic; not in the paper and beyond its linear range): with $k_F/k_0=1.6$, $|\chi^{(1)}|=4$ and $\hbar^2k_0^2/m^\ast=E_{gap}/2\approx1.2$ eV taken from the ITO fit in [6], $\delta\varepsilon\approx1.6\,(U_{hot}/N)/\text{eV}$. So $\delta\varepsilon\sim1$ needs about 0.6 eV per electron, about 100 J/cm³ at $N=1.1\times10^{21}$ cm⁻³.
- Speed: this model has no thermalization delay; the kinetic rise is the pump integral and recovery is $\tau_{el}$. Zhou et al. [17] instead say the rise is limited by electron–electron thermalization. H1 and H3 test exactly this difference.
- The fast term is absent from a cycle-averaged Boltzmann model. By Eq. 42 it overtakes the slow term for modulation faster than about $2/\tau_s$ (my inference, by analogy with the paper's interband statement); relevant to Tier 2.
- Do not copy Eq. 23 or the first form of Eq. 28 as printed.

---
### [6] Secondo et al., "Absorptive loss and band non-parabolicity as a physical origin of large nonlinearity in epsilon-near-zero materials" — Opt. Mater. Express 10, 1545 (2020)

**Source read:** arXiv:2004.00134v1, full text with Appendices A–D (Word-generated PDF, no TeX source; equations transcribed from the PDF pages).
**Role in this project:** The simplest thermal-electron model (starting point of M0): Fermi–Dirac electrons in a non-parabolic band, with $\mu$ and $T_e$ fixed by number and energy conservation. It also supplies an ITO band parameterization.

**What it does.** It gives a fit-parameter-free model of the intensity-dependent index of Drude ENZ oxides pumped below the gap: absorbed energy moves electrons to heavier states, the average effective mass rises, the plasma frequency red-shifts. The numerical Fermi–Dirac version reproduces the AZO index-change spectra of Caspani et al. [2].

**Model / method.**

Momentum-dependent mass (Eq. 5, polar-angle average of $dv_z/dk_z$) and generalized Drude formula (Eq. 6; $\Gamma$ scattering rate, $N$ carrier density):

$$m^\ast(k)=\hbar\left[\frac23\frac{v(k)}{k}+\frac13\frac{dv(k)}{dk}\right]^{-1},\qquad \varepsilon(\omega)=\varepsilon_\infty-\frac{Nq^2}{\varepsilon_0(\omega^2+i\omega\Gamma)m^\ast_{avg}}$$

Average mass (Eq. 8; $\rho$ density of states, $f$ Fermi–Dirac) and the resulting change at the screened plasma frequency (Eq. 11):

$$m^\ast_{avg}=\frac1N\int\frac{f(E,\mu,T_e)}{m^\ast(E)}\rho(E)\,dE,\qquad \delta\varepsilon\big|_{\omega_p}\approx\varepsilon_\infty\frac{\delta m^\ast_{avg}(I)}{m^\ast_{avg}}$$

(Eq. 8 as printed; the left side must be $1/m^\ast_{avg}$ to agree with Eq. 7. Read that way, Eqs. 29–30 with Table 2 reproduce the $m^\ast_{avg}$ column of Table 1 and the 1325 nm ENZ wavelength. My check.)

Analytic estimate with no electron temperature (Eqs. 14, 15, 19), assuming $m^\ast_{avg}\approx m^\ast(k_F)$, saturated velocity at $k_F$, and $\delta k_i=\delta E_i/\hbar v_F$ per absorbed energy $\delta E_i$:

$$\delta\varepsilon(\omega_p)=\varepsilon_\infty\frac{\delta U_{abs}}{N\hbar v_Fk_F},\qquad \delta U_{abs}=I_i(1-e^{-\alpha d})\tau/d,\qquad n_2(\omega_p)=\varepsilon_\infty\frac{\alpha\tau}{2n'(\omega_p)NE_F}$$

$I_i=(1-R)I$ is the intensity just inside the film, $d$ thickness, $\alpha=\alpha_0n_g$ the slow-light-enhanced absorption. $\tau$ is the accumulation time: $\tau_{el}$ for long pulses, $\tau_p$ when $\tau_p<\tau_{el}$.

Full model (Eqs. 21–23), solved for $\mu(I)$, $T_e(I)$, then inserted in Eq. 6; absorptance $A$ from a transfer-matrix calculation:

$$\int f\rho\,dE=N,\qquad \int f(E,\mu(I),T_e(I))E\rho\,dE-\int f(E,\mu_0,T_0)E\rho\,dE=\delta U_{abs}=AI\tau/d$$

Band (Appendix B, Eqs. 29, 30): hyperbola $E=\sqrt{a^2+a^2k^2/b^2}-a$,

$$m^\ast(E)=\left[\frac23\frac{a^2}{\hbar^2b^2(E+a)}+\frac13\frac{a^4}{\hbar^2b^2(E+a)^3}\right]^{-1},\qquad DOS(E)=\frac{b^3}{\pi^2a^3}\sqrt{(E+a)^2-a^2}\,(E+a)$$

For $\tau_p\sim\tau_{el}$ (Eq. 31): $\delta U_{abs}=I_i\eta\alpha_0n_g\tau_p$, $\eta=\max[(\gamma(t)/2)(1+\mathrm{erf}(t/\tau_p))]$, $\gamma(t)=e^{-t/\tau_{el}}$ for $t\ge0$, 1 otherwise. $\Gamma=q/(\mu_em^\ast_{avg}(I=0))$ from the optical mobility, held constant.

**Numbers.**
- Table 2 ($a$ in eV, $b$ in Å⁻¹): ITO 1.22, 0.204 (fit range 0–2 eV); doped ZnO 1.35, 0.211. CdO, GaN and Ga₂O₃ are also listed.
- Implied (my arithmetic): band-edge mass $\hbar^2b^2/a$ = 0.26 $m_0$ (ITO), 0.25 $m_0$ (ZnO); Kane $C=1/2a$ = 0.41 eV⁻¹ (ITO), 0.37 eV⁻¹ (ZnO); ITO $E_F\approx1.06$ eV at the Table 1 density.
- Table 1 (ENZ at 1325 nm). ITO: $\varepsilon_\infty=3.66$, $N=11.3\times10^{20}$ cm⁻³, $m^\ast_{avg}=0.49\,m_0$, $dm^\ast(k_F)/dE=0.36$ (printed unit eV⁻¹; consistent with $m_0$/eV), $\mu_e$ 20–35 cm²/Vs. ZnO: 3.3, 9.19, 0.44 $m_0$, 0.36, 10–40 cm²/Vs.
- AZO estimate (Sec. 3): $\varepsilon_\infty\approx3.5$, $\alpha\sim10^5$ m⁻¹, $\tau=\tau_p=100$ fs, $n'\approx0.3$, $N\approx10^{21}$ cm⁻³, $E_F\approx4.5$ eV, giving $n_2\sim8\times10^{-17}$ m²/W; measured $3\times10^{-17}$ m²/W [2].
- AZO fit (Fig. 3): pump 785 nm, 100 fs FWHM, 0.4 and 0.9 TW/cm²; $\tau_{el}=170$ fs (measured AZO recovery, Clerici et al. 2017); scattering times "10s of (fs)". $m^\ast_{avg}$ rises about 20% at 1 TW/cm² (≈0.435 to ≈0.53 $m_0$, $T_e$ up to ≈5000 K; read from Fig. 3b).

**Key results.**
- Measured $\delta n(\lambda)$ reproduced near ENZ without free parameters (simulated against measured).
- $m^\ast_{avg}$ is linear in intensity, $\delta n$ is not; saturation is built in, $\delta n_{max}=\sqrt{\varepsilon_\infty}-n_0$ (derived).
- Speed: energy accumulates over $\tau$ and leaves to the lattice with $\tau_{el}$; "sub-picosecond rise and fall time" (asserted).

**Assumptions and limitations.**
- Instantaneous thermalization, one $T_e$ for the film, quasi-static peak response; constant $\Gamma$ and density.
- Eq. 19's estimate inserts $E_F\approx4.5$ eV where the derivation has $\hbar v_Fk_F$, about 1.6 eV for this band; that would raise the estimate about threefold (my reading).
- Typos: implicit hyperbola printed with $(E-a)^2$; stray factor $E$ in the inline formula before Eq. 23.

**Use for this project.**
- Adopt Eqs. 21–23 with Eqs. 29–30 as the minimal M0 and as the Fermi–Dirac limit of WP1's $\omega_p^2[f]$; Eq. 5 is the transport mass of [5].
- The ITO fit implies band-edge mass 0.26 $m_0$, below the plan's 0.3–0.4 range; 0.49 $m_0$ is the Fermi-sea average. Baxter et al. [9] use 0.4 $m_e$ as band-edge mass. Decide which mass WP1 samples.
- Constant $\Gamma$ conflicts with [7] and [9], which let damping rise with $T_e$.

---
### [7] Wang et al., "Extended Drude model for intraband-transition-induced optical nonlinearity" — Phys. Rev. Applied 11, 064062 (2019)

**Source read:** abstract only (verbatim OpenAlex copy in `references/text/07_wang2019.abstract.txt`). Main text NOT READ (paywalled, no PDF). Statements tagged **[lib: X]** are how library paper X describes this paper (Baxter = [9], Sarkar = [14], Un = [15], Gurung = [21], Choi = [24]); verify before use.
**Role in this project:** Part of the "established picture" (plan refs [5–8]): a thermal-level model in which both the plasma frequency and the Drude damping depend on electron temperature, relevant to M0 and to the damping $\gamma[f]$ in the WP3 constitutive law.

**What it does.** It presents an "extended Drude model" intended "to fully predict the spectrum" of the nonlinearity caused by intraband transitions. It attributes "the origin of the nonlinearity to be due to both band nonparabolicity and temperature-dependent mobility". The model "matches experiment well" for ITO [abs].

**Model / method.** Equations not available. What the library attributes to it:
- The Drude description of the oxide and the statement that the temperature dependence of the plasma frequency follows from the energy dependence of the effective mass in the conduction band [lib: Baxter Sec. 2].
- The non-parabolic density of states, quoted as Baxter's Eq. 30 with this paper as source ($m_0^\ast$ band-edge mass, $C$ non-parabolicity):

$$D(E)=\frac{1}{2\pi^2}\left(\frac{2m_0^\ast}{\hbar^2}\right)^{3/2}(E+CE^2)^{1/2}(1+2CE)$$

- It is classed among the "thermal models" of the permittivity as a function of electron temperature [lib: Sarkar].
- With a companion paper (Wang et al., Nanophotonics 9, 4287 (2020)): the authors incorporated "multiple scattering mechanisms" and "found that the temperature-dependent damping rate was primarily influenced by ionized impurities and acoustic phonon scattering, increasing at higher electron temperatures" [lib: Gurung, citing both papers together].

**Numbers.**
- $m_0^\ast=0.4\,m_e$ for ITO [lib: Baxter SI Sec. 1, attributed to this paper].
- Interband background permittivity $\varepsilon_b=4$ [lib: Sarkar, attributed to this paper and one other].
- "an almost-linear relationship was found between an effective $\gamma$ and the pulse intensity by fitting the measured spectra to the Drude model" [lib: Baxter Sec. 3].
- Electron density, $C$, $\varepsilon_\infty$, damping, ENZ wavelength, film thickness, pulse parameters, time constants: not available.

**Key results.**
- Non-parabolicity and temperature-dependent mobility both contribute to the nonlinearity (asserted in the abstract; modelled and compared with experiment per the abstract).
- The imaginary part of the ITO permittivity grows under illumination (measured) [lib: Un].
- The response was described as saturable [lib: Un, citing this paper with two others].
- Experimental fitting in this paper and its companion points to a damping that rises more steeply than linearly at high electron temperature [lib: Choi, citing both together].

**Assumptions and limitations.**
- Thermal-electron description, by the library's classification. Whether it treats dynamics, time constants or speed limits is unknown.

**Use for this project.**
- Evidence that damping must change with excitation. This conflicts with the constant damping of [5] and [6] and supports a distribution-dependent $\gamma[f]$ in WP3 and a $T_e$-dependent damping in M0.
- Sarkar et al. [14] used the room-temperature ITO permittivity measured here to fix their deformation potential, so this paper's sample underlies the parameter set of [14, 15] that WP1 starts from (my reading).
- Open the paper for the mobility model $\mu(T_e)$ and its ITO parameters before using it as a baseline. Do not adopt $m_0^\ast=0.4\,m_e$ as a band-edge mass unchecked: the ITO band fit in [6] implies 0.26 $m_0$ at the band edge and about 0.49 $m_0$ averaged over the Fermi sea.

---
### [8] Secondo et al., "Deterministic modeling of hybrid nonlinear effects in epsilon-near-zero thin films" — Appl. Phys. Lett. 120, 031103 (2022)

**Source read:** published PDF, full (cover page + 8 article pages; the file contains no supplementary material).
**Role in this project:** Thermal-level, permittivity-only reference for combined intraband and interband pumping (H2, M2–M4); its averaged effective mass is an equivalent form of the project's $\omega_p^2[f]$.

**What it does.** Builds a parameter-free ("deterministic") model of the transient Drude permittivity of a degenerately doped oxide under a sub-gap pump (heats carriers), an above-gap pump (adds carriers), or both, using a first-principles band structure in place of a non-parabolicity parameter. It is tested on pump–probe reflection from Ga:ZnO, not ITO, then used to predict AZO/GZO, GaN, CdO and ITO films.

**Model / method.** Absorbed energy density from a sub-gap pump, with its decay (Eqs. 1, 2):

$$\Delta U_{IR}(t)=A(\lambda_{IR},t)\frac{I_{0IR}\,t_{IR}}{d}\,\eta_{IR}(t),\qquad \eta(t)=\frac12\Big(\mathrm{erf}\Big(\frac{t}{\sigma}\Big)+1\Big)\begin{cases}1,& t<0\\ e^{-t/\tau_{IR}},& t\ge0\end{cases}$$

$A$ is the absorptivity from a transfer-matrix calculation, $I_{0IR}$ the incident intensity, $t_{IR}$ the pulse FWHM, $d$ the film thickness, $\sigma$ "the time of the Gaussian IR pump pulse", $\tau_{IR}$ the electron–phonon relaxation time.

Carriers from an above-gap pump and the energy they add (Eqs. 5, 6); $\eta_{UV}$ is Eq. 2 with $\tau_{UV}$:

$$\Delta N=A(\lambda_{UV},t)\frac{I_{0UV}\,t_{UV}}{d}\frac{1}{\hbar\omega}\,\eta_{UV}(t),\qquad \Delta U_{UV}\approx\Delta N\Big(\mu+\frac{\Delta N}{\rho(\mu)}\Big)$$

Number and energy conservation fix $\mu$ and $T_e$ of a Fermi–Dirac distribution $f$ at each time (Eqs. 7, 8):

$$N=N_0+\Delta N=\int_0^\infty f(\mu,T_e,E)\rho(E)\,dE,\qquad U=U_0+\Delta U_{IR}+\Delta U_{UV}=\int_0^\infty f(\mu,T_e,E)\rho(E)E\,dE$$

Average optical mass from the band $E(k)$ (Eqs. 3, 4), inserted with $N$ in the Drude permittivity:

$$m^*_{avg}(T_e)=\Big(\frac1N\int_0^\infty\frac{f(\mu,T_e,E)\rho(E)}{m^*(E)}dE\Big)^{-1},\qquad m^*(E)=\Big(\frac{2}{3\hbar^2k}\frac{dE}{dk}+\frac{1}{3\hbar^2}\frac{d^2E}{dk^2}\Big)^{-1}$$

$E(k)$ comes from HSE hybrid DFT tuned to the experimental gap; the band minimum is assumed at Γ, spherical and non-degenerate. There is no current equation and no field propagation.

**Numbers.**
- Sample: 225 nm Ga:ZnO on 24 nm ZnO on c-plane sapphire; $\varepsilon(\lambda_{ENZ}=1700\ \mathrm{nm})=0+i0.41$.
- Pumps: 1890 nm (intraband) and 325 nm (interband), 60 fs; incident peak irradiance 13, 33, 83 GW/cm² and 1, 2, 6, 14 GW/cm². Probe: white light 1150–1600 nm, reflection at 30°.
- Time constants: thermalisation "on the order of 10s fs"; $\tau_{IR}\sim200$ fs "generally". Values used for GZO: not stated.
- ITO rows of Table I (hypothetical films crossing zero at 1.8, 2.0, 2.2 µm): carrier density 4.44, 3.46, 2.77 × 10²⁰ cm⁻³; effective mass 0.39, 0.38, 0.37. In₂O₃ band from HSE with $E_g=2.68$ eV (Table II).
- Predicted for 100 nm ITO crossing at 1.8 µm (read from Fig. 3): $\Delta\varepsilon\approx+0.46$ at 600 GW/cm² pumped at 787 nm; $\Delta\varepsilon\approx-1.75$ at 50 GW/cm² pumped 0.1 eV above the band edge.

**Key results.**
- Peak $\Delta R/R$ against probe wavelength and irradiance is reproduced "without fitting parameters": about −0.2 at 83 GW/cm² (infrared pump) and +0.45 at 14 GW/cm² (ultraviolet pump), read from Fig. 2c,d (measured and simulated).
- $\Delta R/R$ changes sign at short probe wavelengths (near 1200–1250 nm) for both pumps; attributed to impedance matching with air, not to a change of mechanism (measured; explanation asserted).
- Interband pumping is more efficient but slower than intraband (asserted, with Fig. 2a,b). Figures of merit (Eqs. 9, 10) rank GaN best for intraband and CdO for interband use (simulated).

**Assumptions and limitations.**
- Thermalised electrons, single-exponential decay, no lattice temperature equation. The authors call the model "kinetic", but it tracks only $\mu$ and $T_e$ (my reading).
- Absorption uniform through the film, $d<1/\alpha$; holes, Moss–Burstein shift and gap renormalisation neglected; expected to fail at high ultraviolet irradiance through band filling (authors).
- A change of Drude damping is not discussed (my reading).

**Use for this project.**
- Eqs. 3–4 are the integrated-by-parts form of the plan's $\omega_p^2[f]=\frac{e^2}{3\varepsilon_0}\int g v^2(-\partial f/\partial E)dE$: $\omega_p^2=e^2N/(\varepsilon_0m^*_{avg})$ for any isotropic $E(k)$ and any $f$ (my derivation). Use as a unit test, and to replace the Kane band by a DFT band.
- Eqs. 5–8 are a ready thermal baseline for M2–M4.
- The sign change with probe wavelength is electromagnetic; relevant to H3 and WP5.
- The plan lists [8] among two-temperature models. It is an energy-balance model validated on GZO.
- Data "available from the corresponding author upon reasonable request."

---
### [9] Baxter et al., "Understanding the nonlinear optical response of epsilon near zero materials in the time-domain" — arXiv:2110.14806v1 (2021); published as Adv. Photonics Res. 4, 2200280 (2023)

**Source read:** arXiv v1 full text with its Supplemental Information (SI) + TeX source; figures checked in the PDF. The published version was not available and may differ. Equation and figure numbers follow the arXiv v1 PDF (SI equations are 11–36).
**Role in this project:** The FDTD + two-temperature baseline for WP3 and the thermal-electron reference for M0 in WP4; its $\omega_p(T_e)$ is the Fermi–Dirac limit of the project's $\omega_p^2[f]$.

**What it does.** A two-temperature model (TTM) is solved in every cell of a 3D FDTD grid and coupled to Maxwell's equations through a Drude current whose plasma frequency depends on the local electron temperature via the non-parabolic band. Two parameters are tuned to the 310 nm ITO film of Alam, De Leon and Boyd (Science 2016). The model then reproduces that paper's z-scan $n_2$ map, intensity-dependent reflectance and transmittance, and pump–probe rise and fall times, and is applied to a gold-antenna metasurface on ITO.

**Model / method.**

Drude permittivity and Kane band (Eqs. 1, 2), with $\omega_p=\sqrt{Ne^2/\varepsilon_0 m^\ast}$ and zero crossing $\lambda_{\rm ZC}=2\pi\sqrt{\varepsilon_\infty}\,c/\omega_p$:

$$\varepsilon(\omega)=\varepsilon_\infty-\frac{\omega_p^2}{\omega^2+i\gamma\omega},\qquad \frac{\hbar^2k^2}{2m_0^\ast}=E+CE^2$$

Plasma frequency (Eq. 3), derived in SI Sec. 1 from the Boltzmann equation in the relaxation-time approximation for an isotropic band, with $\partial f/\partial\mathbf k\to\partial f_{FD}/\partial\mathbf k$ in the force term:

$$\omega_p^2(T_e)=\frac{e^2}{3m_0^\ast\varepsilon_0\pi^2}\left(\frac{2m_0^\ast}{\hbar^2}\right)^{3/2}\int_{0}^{\infty}dE\,(E+CE^2)^{3/2}(1+2CE)^{-1}\left(-\frac{\partial f_{FD}(E,T_e)}{\partial E}\right)$$

$f_{FD}=(\exp((E-\mu)/k_BT_e)+1)^{-1}$; $E$ is measured from the band edge. The SI copies (Eqs. 21, 23) print $(E+CE)^{3/2}$ and $+\partial f_{FD}/\partial E$ (as printed; typos, Eq. 3 is the consistent form).

Extended TTM with athermal energy density $U$ in J/m³ (Eq. 4a–c). Written down but not used:

$$\frac{\partial U}{\partial t}=-\frac{U}{\tau_{ee}(T_e)}-\frac{U}{\tau_{ep}(T_e,T_l)}+P(\mathbf r,t)$$
$$C_e(T_e)\frac{\partial T_e}{\partial t}=-g_{ep}(T_e,T_l)(T_e-T_l)+\frac{U}{\tau_{ee}(T_e)},\qquad C_l(T_l)\frac{\partial T_l}{\partial t}=g_{ep}(T_e,T_l)(T_e-T_l)+\frac{U}{\tau_{ep}(T_e,T_l)}$$

TTM actually solved (Eq. 5a,b), local in each cell, with no heat-diffusion term:

$$C_e(T_e)\frac{\partial T_e(t,\mathbf r)}{\partial t}=-g_{ep}(T_e,T_l)\,(T_e-T_l)+P(\mathbf r,t),\qquad C_l(T_l)\frac{\partial T_l(t,\mathbf r)}{\partial t}=g_{ep}(T_e,T_l)\,(T_e-T_l)$$

Maxwell and Drude current (Eqs. 6–8):

$$\mu_0\frac{\partial\mathbf H}{\partial t}=-\nabla\times\mathbf E,\qquad \varepsilon_\infty\varepsilon_0\frac{\partial\mathbf E}{\partial t}=\nabla\times\mathbf H-\mathbf J,\qquad \frac{\partial\mathbf J}{\partial t}=-\gamma\mathbf J+\varepsilon_0\omega_p^2(T_e)\,\mathbf E$$

Coupling: the absorbed power $P=\mathbf J\cdot\mathbf E$ ("given by Poynting's theorem"; no cycle average is stated) drives Eq. 5a, and $\omega_p(T_e)$ from Eq. 3 enters Eq. 8. Bound electrons are a constant $\varepsilon_\infty$.

Damping, introduced in Sec. 3 as a fit (Eq. 10). It depends on $T_e$ only, not on $T_l$:

$$\gamma(T_e)=\gamma_0\Big(1+\frac{T_e}{T_0}\Big)$$

TTM coefficients (SI Sec. 3, Eqs. 30, 33, 34). $\mu(T_e)$ is found numerically from $n(T_e)=\int_0^\infty f_{FD}D\,dE=\int_0^{E_F}D\,dE$ (Eqs. 35, 36):

$$D(E)=\frac{1}{2\pi^2}\Big(\frac{2m_0^*}{\hbar^2}\Big)^{3/2}(E+CE^2)^{1/2}(1+2CE),\qquad C_e(T_e)=\int_0^\infty\frac{\partial f_{FD}(E,\mu,T_e)}{\partial T_e}D(E)\,E\,dE$$
$$g_{ep}(T_e)=\frac{A_{ep}\pi\hbar k_B}{D(E_F)}\int_{-\infty}^\infty D^2(E)\,\frac{\partial f_{FD}(E,\mu,T_e)}{\partial T_e}\,dE$$

(Eq. 34 as printed. With $A_{ep}$ in eV² and $\partial f/\partial T_e$ it does not have units of W m⁻³ K⁻¹. The cited source, Lin, Zhigilei and Celli 2008, has $-\partial f/\partial E$ and prefactor $\pi\hbar k_B\lambda\langle\omega^2\rangle/g(\varepsilon_F)$ — background, not from this paper. Check before implementing.)

Why athermal electrons are dropped (SI Sec. 2): the Fermi-liquid lifetime $\tau=\frac{128}{\pi^2\sqrt3\,\omega_p}\big(\frac{E_F}{E-E_F}\big)^2\approx\frac{7.5}{\omega_p}\big(\frac{E_F}{E-E_F}\big)^2$ (Eq. 25) is applied to an initial step perturbation $\delta f(E,0)=A\{f_{FD}(E-E_{pump})[1-f_{FD}(E)]-f_{FD}(E)[1-f_{FD}(E+E_{pump})]\}$ (Eq. 28).

Numerical scheme. In-house parallel 3D FDTD with Bloch-periodic in-plane boundaries for oblique incidence. The TTM is integrated by Runge–Kutta in the same code; "at each step in the FDTD simulation" $\mathbf E$, $\mathbf H$, $\mathbf J$, $T_e$, $T_l$ are updated together, so slow thermal and fast field variables share one time step. $C_e$, $g_{ep}$ and $\omega_p^2$ are precomputed and fitted by fifth-order polynomials in $T_e$. Pump–probe uses two runs: the pump run stores $T_e(\mathbf r,t)$ and $T_l(\mathbf r,t)$, and a weak probe is then propagated through that stored time-varying medium.

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| Sample | 310 nm ITO on glass (Alam 2016 film); glass index not stated | Sec. 3 |
| Room-temperature Drude | $\omega_p=2\pi\times473$ THz, $\gamma=0.0468\,\omega_p$, $\varepsilon_\infty=3.8055$, $\lambda_{\rm ZC}=1237$ nm | Sec. 3 |
| Band | $m_0^*=0.4\,m_e$, $C=0.4191$ eV⁻¹, $E_F\approx1$ eV; electron density not stated | SI Sec. 1, 3 |
| Lattice heat capacity | $C_l=2.54\times10^{6}$ J m⁻³ K⁻¹, constant | SI Sec. 3 |
| Electron–phonon prefactor | $A_{ep}=5.25\times10^{-4}$ eV², fitted to the fall time | SI Sec. 3 |
| Damping slope | $T_0=15000$ K, fitted to the R/T saturation | Eq. 10 |
| $\omega_p(T_e)$ | about $2.97\to1.9\times10^{15}$ rad/s for 300 → 20 000 K (read from plot) | Fig. 1 |
| Athermal decay | exponential fit $t_0=4.2$ fs for ITO ($E_F=E_{pump}=1.0$ eV, 300 K, called an upper bound); 30 fs for gold | Figs. 9, 10 |
| Pulses | p-polarized; 150 fs FWHM (stated for Sec. 4 and SI Sec. 4); incident peak intensity | Sec. 3, 4 |
| FDTD cell size, time step | not stated | — |

**Key results.**
- $n_2$ map (simulated against measured, Fig. 3): effective $n_2$ peaks at 1250 nm and 60° (steps 50 nm, 10°), with $n_2\approx0.16$ cm²/GW and $\beta\approx-13000$ cm/GW; Alam: 0.11 cm²/GW and −7500 (printed "cm²/GW"; unit typo). $n_2=\overline{\Delta n}/\overline I$ (Eq. 9) comes from the space–time mean $T_e$, not from a simulated z-scan.
- R and T against intensity at 1240 nm, 30° (Fig. 4): T rises 0.13 → 0.4 and R falls 0.19 → 0.02 over 0 → 250 GW/cm²; measured 0.12 → 0.4 and 0.2 → 0.04. This needed Eq. 10; with constant $\gamma$ saturation sets in only above 250 GW/cm².
- Pump–probe transmission at 1240 nm, 30° (Fig. 5): 10–90 % rise 170 fs and fall 385 fs; measured < 200 fs and ≈ 360 fs for a 150 fs pulse. Peak $\Delta T\approx0.18$; the pump intensity is not stated.
- Effective index change in that run: real part ≈ +0.5, imaginary part ≈ −0.2, a metal-to-dielectric transition within the pulse (Fig. 6a; simulated).
- Largest $\Delta n$ at 2 GW/cm²: pump at 1228 nm, 35°, $\Delta n\approx0.02$; chirp is irrelevant. The best wavelength rises to ≈ 1360 nm at 110 GW/cm² as the ENZ point redshifts (Figs. 11, 12; simulated).
- Gold-bar metasurface on the film at 2 GW/cm²: $\Delta n$ five times the bare film and $\Delta T<0$, the sign seen by Alam et al. 2018 (Fig. 7; simulated).
- "Including hydrodynamic terms has very little effect on the nonlinearity when solved in conjunction with the TTM" (asserted; tests not shown).

**Assumptions and limitations.**
- Electrons are a Fermi–Dirac distribution at all times. The 4 fs justification assumes an instantaneous perturbation and a single-exponential fit (authors).
- $\gamma(T_e)$ is a first-order fit, and $g_{ep}$ uses a formula for metals with a fitted prefactor (authors). Agreement with Alam's data is therefore partly by construction.
- No interband transitions, no bound-electron nonlinearity, one sample (authors).
- Eq. 8 is Drude with a time-varying coefficient $\varepsilon_0\omega_p^2(t)$ on $\mathbf E$; there is no term in $\dot m^*/m^*$ acting on $\mathbf J$ (my reading).
- No heat transport between cells (my reading).
- The intensity normalising $n_2$ is "time-averaged pulse intensity" in the main text and intensity "in the ITO layer" in an SI footnote (my reading: ambiguous).

**Use for this project.**
- Adopt Eq. 3 as the thermal check of the plan's $\omega_p^2[f]$. With $g=D(E)$ and $v=\hbar k/[m_0^*(1+2CE)]$ the two are identical for $f=f_{FD}$; the Fig. 1 curve is the target. (My estimate, not in the paper: $E_F=1$ eV with these band parameters gives $n\approx1.9\times10^{21}$ cm⁻³ and $\omega_p\approx2.9\times10^{15}$ rad/s at $T=0$.)
- Baseline for WP3 step 2: Figs. 4 and 5 with the table above. Cell size, time step, glass index and the pulse length of Sec. 3 are missing and must be chosen or requested. arXiv v1 has no data-availability statement; field and temperature movies are in its supplementary material.
- The plan's WP3 current equation carries $(\gamma+\dot m^*/m^*)J$; Eq. 8 does not. Switch that term off to benchmark against this paper.
- The claim that thermalisation takes ≈ 4 fs is the opposite premise to H1 and M1.
- The plan cites [9] for an *extended* TTM in M0. The paper writes that form (Eq. 4) but solves the simple one (Eq. 5), with Fermi integrals and $\mu(T_e)$ reduced to polynomial fits.
- The probe is propagated through the stored time-varying film, so electromagnetic effects on the probe exist at temperature level. This is the comparison H3 needs.
- $P=\mathbf J\cdot\mathbf E$ is instantaneous here and cycle-averaged in the plan; check that the difference is negligible for a 150 fs pump.

---
### [10] Li et al., "Invertible optical nonlinearity in epsilon-near-zero materials" — arXiv:2203.09600v1 (2022); published as Phys. Rev. Research 5, 013198 (2023)

**Source read:** arXiv v1 full text + TeX source; figures checked in the PDF. The Supplemental Material (set-up, Figs. S1–S4, Table S2 with fitted Drude and band parameters) is not in the library; the published version was not checked. Equation numbers follow the arXiv PDF.
**Role in this project:** Supports the WP5 probe-wavelength scan and bears on H2: a carrier-density increase competes with band non-parabolicity in setting the sign of the plasma-frequency shift.

**What it does.** Broadband-probe pump–probe transmission on 23 nm ITO films with ENZ wavelengths of 1220, 1400 and 1500 nm, pumped at the ENZ wavelength (intraband) or at 360 nm (interband). Near-infrared pumping always redshifts the ENZ resonance. Pumping at 360 nm redshifts it in the 1220 nm film and blueshifts it in the other two.

**Model / method.** Drude permittivity with $\omega_p=\sqrt{ne^2/\epsilon_0m^*}$ (Eq. 2) and a band with non-parabolicity $1/E_g$ (Eq. 3):

$$\frac{k^2\hbar^2}{2m^{*}}=E(k)+\frac{E^2}{E_g}$$

Electrons are a Fermi–Dirac distribution $f_0(\mu,T)$ at all times. The expressions for $\omega_p$, $n(\mu,T)$ and $U(\mu,T)$ are taken from Guo et al., Nat. Photon. 10, 267 (2016) and are not reproduced.

- Intraband pump: $n$ constant, $U(\mu,T)=U(\mu,300\,\mathrm K)+dU$, with $dU$ the absorbed energy density.
- Interband pump: the density rises by the absorbed photon density (Eqs. 5, 6). How $T$ is fixed here is not stated in the main text.

$$n(\mu,300\,\mathrm K)+dn=\frac{1}{\pi^2}\int_0^\infty dE\left(\frac{2m}{\hbar^2}\Big(E+\frac{E^2}{E_g}\Big)\right)^{1/2}\frac{m}{\hbar^2}\Big(1+\frac{2E}{E_g}\Big)f_0(\mu,T),\qquad dn=\frac{FA\cos\theta}{\pi r^2h\,E_{photon}}$$

($F$ pulse energy, $A$ absorption, $\theta$ incidence angle, $r$ beam radius, $h$ film thickness.)

- Optics: transfer-matrix method; observable $\Delta A=-\log(1+\Delta T/T_0)$ (Eq. 4). $\omega_p$ is fitted to one spectrum per delay. There is no current equation and no time-domain propagation.

Mechanism of the inversion: added carriers raise $\omega_p$ (blueshift), but they also raise the Fermi level, which in a non-parabolic band raises the average effective mass (redshift). Large non-parabolicity (small $E_g$) or a low Fermi level lets the mass term win.

**Numbers.**
- Commercial film: 23 nm ITO on glass, ENZ at 1240 nm, $\omega_p=1.9$ eV, $\gamma_p=0.17$ eV (Fig. 1a); $\epsilon_\infty$ not stated.
- Sputtered, annealed films: 23 nm, ENZ at 1220, 1400, 1500 nm (Fig. 1c); fitted parameters are in Table S2 (not available).
- Model illustration: $E_g=2.39$ eV, $m^*=0.26\,m_0$; $E_g$ scanned 2–2.8 eV (Fig. 3).
- Pump at the ENZ wavelength or at 360 nm (the band gap of these films); probe about 1000–1600 nm (figure axes). Pulse durations and angles: not stated in the main text.
- Near-infrared pump, 1220 nm film: $\Delta A=100$ mOD at 1000 µJ/cm², "corresponding to $\Delta T/T=35\%$" (as printed; Eq. 4 with a base-10 logarithm gives about 26 %).
- Decay of $dU$ (Figs. 4f, 5f). Near-infrared pump: 149, 135, 114 fs for the 1220, 1400, 1500 nm films. 360 nm pump: 142 fs (1220 nm); 119 and 504 fs (1400 nm); 130 and 499 fs (1500 nm). The ≈ 500 fs component is attributed to recombination.

**Key results.**
- Every transient spectrum is bipolar about the ENZ wavelength: for a redshift $\Delta A<0$ on the short-wavelength side and $>0$ on the long side; for a blueshift the pattern is reversed (measured, Figs. 4, 5).
- Under 360 nm pumping the 1220 nm film, which has the largest non-parabolicity, redshifts; the 1400 and 1500 nm films blueshift (measured; modelled with one fitted $\omega_p$ per spectrum).

**Assumptions and limitations.**
- Thermalised electrons only; the rise and any non-thermal stage are not modelled (my reading).
- The same number of carriers is assumed excited in every sample; the Burstein–Moss shift is neglected (authors).

**Use for this project.**
- WP5: probe on both sides of the ENZ wavelength. The sign flip across it marks a resonance shift; only reversal of the whole pattern marks a change of mechanism.
- H2 and M2–M4: a carrier-density increase does not by itself give a blueshift.
- The inversion here comes from one-photon above-gap pumping and sample choice. It is not the intensity-driven sign reversal in the plan's anomaly table.
- The non-parabolicity $1/E_g=0.418$ eV⁻¹ matches [9]; $m^*=0.26\,m_0$ does not ($0.4\,m_e$ there).
- No data-availability statement in arXiv v1.

---
### [11] Scalora et al., "Electrodynamics of conductive oxides: intensity-dependent anisotropy, reconstruction of the effective dielectric constant, and harmonic generation" — arXiv:2003.01059 (2020); published as Phys. Rev. A 101, 053828 (2020)

**Source read:** arXiv PDF, full (no TeX source; equations read from the page images). Published version not checked. Gaussian units.
**Role in this project:** The hydrodynamic–Maxwell level of theory that WP3 defers (nonlocal and hot-electron terms), and a source of constitutive forms for a changing mass versus a changing density.

**What it does.** Extends the model of [12] for a 20 nm ITO layer. Free- and bound-electron polarizations are integrated with the vector Maxwell equations in the time domain, and an effective permittivity is reconstructed as $\varepsilon_{jj}\approx1+4\pi P_j/E_j$ averaged over the layer. Nonlocal terms make the response anisotropic and blueshift the ENZ crossing; hot electrons redshift it. Simulation only, apart from an ellipsometry fit.

**Model / method.** Free-electron polarization (Eq. 2) in scaled coordinates $\varsigma=y/\lambda_0$, $\xi=z/\lambda_0$, $\tau=ct/\lambda_0$, with $\lambda_0=1\ \mu$m and $\tilde\gamma_f=\gamma_f\lambda_0/c$:

$$\ddot{\mathbf P}_f+\tilde\gamma_f\dot{\mathbf P}_f=\frac{n_{0,f}e^2\lambda_0^2}{m_0^*c^2}\mathbf E-\frac{e\lambda_0}{m_0^*c^2}(\nabla\cdot\mathbf P_f)\mathbf E+\Big(\sum_{l=1,2,3}(-\tilde\Lambda)^l(\mathbf E\cdot\mathbf E)^l\Big)\mathbf E+\frac{e\lambda_0}{m_0^*c^2}\dot{\mathbf P}_f\times\mathbf H$$
$$\qquad+\frac{3E_F}{5m_0^*c^2}\Big(\nabla(\nabla\cdot\mathbf P_f)+\frac12\nabla^2\mathbf P_f\Big)-\frac{1}{n_{0,f}e\lambda_0}\Big[(\nabla\cdot\dot{\mathbf P}_f)\dot{\mathbf P}_f+(\dot{\mathbf P}_f\cdot\nabla)\dot{\mathbf P}_f\Big]$$

In order: Drude drive; Coulomb term from continuity; hot-electron series; magnetic Lorentz force; pressure and viscosity, with $E_F=\frac{\hbar^2}{2m_0^*}(3\pi^2n_{0f})^{2/3}$; convection. The hot-electron series expands $1/m_f^*(T_e)$ for $m_f^*(T_e)\approx m_e(0.033+aK_BT_e)$, with $T_e$ proportional to absorbed energy and $\tilde\Lambda=aK_B\,n_{0,f}e^2\lambda_0^2\sigma_0\tau_0/(m_0^*c^2)$ ($\sigma_0$ a conductivity, $\tau_0$ the pulse duration). It is instantaneous in $|\mathbf E|^2$, and the left-hand side has no $\dot m^*$ term. Bound electrons obey a nonlinear Lorentz equation (Eq. 3).

Changing carrier density (Eqs. 8, 10, unscaled; mass changes neglected there):

$$\dot n=-\frac1e\nabla\cdot\dot{\mathbf P}_f+\Omega\frac{\partial}{\partial t}\mathbf E\cdot\mathbf E,\qquad \ddot{\mathbf P}_f-\frac{\dot n}{n}\dot{\mathbf P}_f+(\dot{\mathbf P}_f\cdot\nabla)\Big(\frac{\dot{\mathbf P}_f}{ne}\Big)+\gamma\dot{\mathbf P}_f=\frac{ne^2}{m}\mathbf E+\frac{e}{mc}\dot{\mathbf P}_f\times\mathbf B-\frac{e\nabla p}{m}$$

with $p=p_0(n/n_0)^{5/3}$, $p_0=n_0E_F$. The resulting Eq. 11 adds $+\frac{\Omega e^2\lambda_0^2}{m_0^*c^2}(\mathbf E\cdot\mathbf E)\mathbf E$ (a blueshift) and $+\frac{\Omega}{n_0}\dot{\mathbf P}_f\frac{\partial}{\partial\tau}(\mathbf E\cdot\mathbf E)$. Its pressure prefactor is $\frac{5E_F}{3m_0^*c^2}$, against $\frac{3E_F}{5m_0^*c^2}$ in Eq. 2 (as printed; inconsistent).

**Numbers.**
- 20 nm ITO (in vacuum in Fig. 1), ENZ near 1240 nm; ≈ 100 fs p-polarized pulses at 60°; incident peak 1 MW/cm² (linear), 10, 20, 25 GW/cm².
- $n_0\sim10^{20}$ cm⁻³, $m_0^*\sim0.033\,m_e$, $aK_B=1.75\times10^{-6}$ per K (Fig. 9), nonlocal parameter $\eta=3E_F/(5m_0^*c^2)\approx10^{-5}$. Drude–Lorentz fit values: not stated.

**Key results.** All simulated.
- With nonlocal terms the layer-averaged $\varepsilon_{yy}$ and $\varepsilon_{zz}$ have ENZ crossings that are blueshifted and "mismatched by nearly 10 nm" (Fig. 4b); at 20 GW/cm² their minima are about 100 nm apart (Fig. 7b).
- $\varepsilon_{yy}(z)$ oscillates with a period of a few nanometres across the layer (Fig. 3b); extra absorption resonances appear near 700 and 900 nm (Fig. 1a).
- Hot-electron terms up to seventh order are needed for stability at 10–20 GW/cm², although the in-film intensity is enhanced only two to three times. 25 GW/cm² is the threshold of optical bistability (Fig. 7a).
- The local and nonlocal linear $|\langle\varepsilon_{zz}\rangle|$ already differ by order unity (Fig. 5).

**Assumptions and limitations.**
- The hot-electron response is "apparently instantaneous". Proposed cure: $m_f^*(T_e)\approx m_0^*(1+aT_e+bT_e^2)$ with $T_e(t)=\Lambda\int_{-\infty}^t\dot{\mathbf P}_f\cdot\mathbf E\,dt'$ (authors). No relaxation, no lattice, no electron distribution.
- $n_0$ and $m_0^*$ are chosen to match second-harmonic amplitudes. $m_0^*$ is about ten times below the 0.26–0.40 $m_e$ of [8–10]; the ratio still gives $\omega_p\approx3.1\times10^{15}$ rad/s (my calculation).
- The density source $\Omega\,\mathbf E\cdot\mathbf E$ follows the field and vanishes after the pulse (my reading).
- (Background, not from this paper) the usual hydrodynamic coefficient is $\beta^2=\tfrac35v_F^2=6E_F/(5m^*)$, which matches neither Eq. 2 nor Eq. 11.

**Use for this project.**
- This paper puts a mass change only into the coefficient of $\mathbf E$ and a density change into an extra $-(\dot n/n)\dot{\mathbf P}_f$ term. The plan's WP3 equation has $+(\dot m^*/m^*)J$ and no $\dot n$ term. The conventions differ; derive the form from the kinetic equation.
- If nonlocal terms are added after Gate 3, keep the structure of Eq. 2 but re-derive the coefficient and use measured $n$ and $m^*$.
- With nonlocal terms the permittivity retrieved from a 20 nm layer differs from the local one, so the size of an "index change" depends on the reference.

---
### [12] Rodríguez-Suné et al., "Study of second and third harmonic generation from an indium tin oxide nanolayer: influence of nonlocal effects and hot electrons" — arXiv:1910.00313 (2019); published as APL Photonics 5, 010801 (2020)

**Source read:** arXiv PDF, full (no TeX source; equations read from the page images). Published version not checked. Gaussian units.
**Role in this project:** Source of the plan's statement that nonlocal effects blueshift the ENZ resonance by tens of nanometres in the thinnest films (WP3, deferred correction); first statement of the hydrodynamic model detailed in [11].

**What it does.** Measures second- and third-harmonic generation (SHG, THG) from a 20 nm ITO film near its ENZ wavelength against pump wavelength and angle, and compares with a time-domain hydrodynamic–Maxwell model with nonlocal terms and a temperature-dependent effective mass. SHG spectra and angular curves are reproduced with the free-electron mass and density as the only free parameters.

**Model / method.** Free-electron polarization with a temperature-dependent mass (Eq. 3), where $E_F=\frac{\hbar^2}{2m^*}(3\pi^2n_0)^{2/3}$:

$$\ddot{\mathbf P}_f+\gamma_f\dot{\mathbf P}_f=\frac{n_{0,f}e^2}{m_f^*(T_e)}\mathbf E-\frac{e}{m_f^*(T_e)}\mathbf E(\nabla\cdot\mathbf P_f)+\frac{e}{m_f^*(T_e)}\dot{\mathbf P}_f\times\mathbf H+\frac{3E_F}{5m_f^*(T_e)}\Big(\nabla(\nabla\cdot\mathbf P_f)+\nabla^2\mathbf P_f\Big)-\frac{1}{n_{0,f}e}\Big[(\nabla\cdot\dot{\mathbf P}_f)\dot{\mathbf P}_f+(\dot{\mathbf P}_f\cdot\nabla)\dot{\mathbf P}_f\Big]$$

The viscosity term is $\nabla^2\mathbf P_f$ here and $\tfrac12\nabla^2\mathbf P_f$ in [11] (as printed in each). The mass enters only as $1/m_f^*(T_e)$ on the right-hand side; there is no $\dot m^*$ term.

Mass and its expansion (Eqs. 5, 7), after assuming $\mathbf J=\sigma_0\mathbf E$ and dropping the volume integral:

$$m_f^*(T_e)\approx m_0^*+\alpha K_BT_e=m_0^*+\alpha K_B\Lambda\int\!\!\int\mathbf J\cdot\mathbf E\,d\mathbf r^3dt,\qquad \frac{n_{0,f}e^2}{m_f^*(T_e)}\mathbf E\approx\frac{n_{0,f}e^2}{m_0^*}\mathbf E-\tilde\Lambda(\mathbf E\cdot\mathbf E)\mathbf E+\tilde\Lambda^2(\mathbf E\cdot\mathbf E)^2\mathbf E$$

Bound electrons follow a nonlinear Lorentz oscillator (Eq. 4); the local permittivity is Drude plus one Lorentz term (Eq. 2). Fields are propagated by a spectral method and the material equations by a predictor–corrector scheme; grid and step sizes are not stated.

**Numbers.**
- Sample: 20 nm ITO, RF-sputtered (In₂O₃/SnO₂ 90/10 wt %) on fused silica, annealed in vacuum at 600 °C. Re ε crosses zero near 1260 nm (ellipsometry); the SHG enhancement appears near 1240 nm.
- Laser: 1140–1600 nm, 1 kHz, 50–102 fs FWHM, incident power density 1–2 GW/cm².
- Measured SHG efficiency at 60°: peak ≈ 0.9 × 10⁻⁷ reflected and ≈ 0.6 × 10⁻⁷ transmitted (read from Fig. 2a).
- Simulation: 1 GW/cm², ≈ 100 fs plane-wave pulses, $m_0^*=0.033\,m_e$, $n_{0,f}\sim10^{20}$ cm⁻³. Drude–Lorentz fit values: not stated.

**Key results.**
- The blueshift, exact wording. Abstract: "nonlocal effects, which blueshift the epsilon-near-zero resonance by tens of nanometers". Body: the chosen mass and density "combine to yield a Fermi energy that imparts a blueshift of nearly 50nm to the SHG peak spectral position. The magnitude of $\tilde\Lambda$ is chosen so that its effect is to redshift the SHG spectra by approximately 15nm, for a net blueshift of ~35nm" (simulated transmitted SHG spectra, Fig. 3c).
- Transmitted and reflected SHG maxima are about 20 nm apart, even without nonlocal terms (measured and simulated).
- THG at 1240 nm has a minimum near 30° incidence (measured and simulated); whether bound or hot electrons cause it is left open.
- Field-intensity enhancement inside the layer is "a mere factor of two" (simulated).

**Assumptions and limitations.**
- The hot-electron response is instantaneous in $|\mathbf E|^2$; a linear $m^*(T_e)$ replaces the two-temperature model, justified for "a few thousand degrees Kelvin" (authors). No distribution, no relaxation.
- The blueshift is a simulated shift of the SHG peak under chosen $m_0^*$ and $n_0$. It is not a measured shift of the linear ENZ resonance (my reading).
- With $m_0^*=0.033\,m_e$ and $n=10^{20}$ cm⁻³ their formula gives $E_F\approx2.4$ eV and $v_F\approx5\times10^{6}$ m/s. With $m^*\approx0.4\,m_e$ and $n\approx10^{21}$ cm⁻³ it gives about 0.9 eV and $0.9\times10^{6}$ m/s, so the pressure coefficient $3E_F/5m^*$ is about 30 times smaller (my calculation).

**Use for this project.**
- The plan's "blueshift of tens of nanometres in the thinnest films" is correctly attributed, but one 20 nm film was studied and the number depends on the parameters above. Recompute the shift with measured ITO parameters before using it as the expected size of the correction.
- The authors warn that ellipsometry data from a very thin layer "may intrinsically contain shifts due to nonlocal effects" and should not be transferred to another thickness. Relevant to per-sample parameters in WP4.
- No data-availability statement in the arXiv version.

---
### [13] Tamashevich et al., "Field theory description of the non-perturbative optical nonlinearity of epsilon-near-zero media" — arXiv:2308.03181v2 (2023); published as APL Photonics 9, 016105 (2024)

**Source read:** arXiv v2 full text + TeX source (appendices A–D, the path-integral derivations, skimmed). Published version not checked. Equation numbers follow the arXiv v2 PDF.
**Role in this project:** Context only: a closed-form phenomenological $\varepsilon(I)$ for ITO. It is not one of the mechanisms M0–M4.

**What it does.** Derives the intensity-dependent permittivity of a Drude medium from a Huttner–Barnett–Bechler Lagrangian by path integrals, then replaces the perturbative Kerr interaction by a Born–Infeld Lagrangian to obtain a non-perturbative $\varepsilon(I)$. A saturation term is added by hand.

**Model / method.** Non-perturbative and saturated permittivities (Eqs. 21, 22):

$$\varepsilon^{(np)}(I)=\varepsilon_{lin}+\frac a2\left(\frac{I/I_c}{\sqrt{1-I/I_c}}\right),\qquad \varepsilon^{(np)}(I)=\varepsilon_{lin}+\frac a2\left(\frac{I/I_c}{\sqrt{1-I/I_c+(I/I_0)^2}}\right)$$

with $a=4n_0n_2I_c$ from the low-intensity Kerr limit, $I_c=2cn_0b^2$ a free critical intensity ($b$ is the Born–Infeld parameter), $I_0=|\delta\varepsilon/2n_0n_2|$ and $\delta\varepsilon=\varepsilon_{sat}-\varepsilon_{lin}$.

**Numbers.**
- ITO Drude fit (Sec. V B): $\varepsilon_\infty=3.99$, $\omega_p=3.16\times10^{15}$, $\gamma=0.185\times10^{15}$, at $\omega=1.57\times10^{15}$ (printed "Hz"; angular frequency by context), $\lambda_{ENZ}=1200$ nm.
- Fig. 4: $\varepsilon=0+0.45i$, $\chi^{(3)}=7.58\times10^{-17}$ m²/V², $n_0=0.4743$, $n_2=0.9517$ cm²/GW, $I_c=300$ GW/cm². Fig. 6: $I_c=100$ GW/cm², $\varepsilon_{sat}=1.252+0.606i$.

**Key results.**
- The Born–Infeld form and the square-root expansion model of Reshef et al. (2017) agree at low intensity and deviate above about 20 GW/cm² (Fig. 4; derived).
- The saturated form follows the Reshef et al. model with $\chi^{(3)}$ to $\chi^{(7)}$ (Fig. 6). The comparison is with that fitted curve, not with raw data.

**Assumptions and limitations.**
- Scalar field, intensity dependence only: no time dependence or relaxation; $I_c$ is real, so nonlinear loss cannot be represented; no nonlocal or surface terms (authors).
- No electron distribution, temperature or carrier density, and no current equation (my reading).

**Use for this project.**
- At most a two-parameter fitting form for steady-state $\Delta n(I)$ curves. It cannot address rise time, recovery or sign reversal.

---

## Kinetic theory of low-electron-density Drude materials
### [14] Sarkar et al., "Electronic and thermal response of low-electron-density Drude materials to ultrafast optical illumination" — Phys. Rev. Applied 19, 014005 (2023)

**Source read:** arXiv:2207.14528v1 (29 Jul 2022), full text incl. appendices + TeX source; figure values read from the vector paths of the PDF. Published version not seen. arXiv title: "The electronic and thermal response of low electron density Drude materials to ultrafast optical illumination"; Sarkar and Un are equal first authors. Equation numbers are arXiv-v1 numbers (they match how [15, 16] cite this paper: Eq. 3 = BE, 15 = q_TF, 17 = e-e rate).
**Role in this project:** Equation-level basis of the WP1 solver (model M1; B1 e-e kernel; e-ph term) and of the two-temperature parameters used in M0.

**What it does.** Writes a Boltzmann equation (BE) for the isotropic conduction-electron distribution f(𝓔,t) of ITO in a Kane band, with full (non-RTA) e-e and e-ph collision integrals and a phonon temperature. Derives energy-dependent collision rates and two-temperature-model (TTM) parameters, and solves one example (220 fs, 2.5 GW/cm²). Headline: momentum conservation cuts e-ph coupling strongly because k_F < q_D/2; e-e scattering is ~10× faster than in Au; the electron heat capacity is ~5× smaller than in Au, so electrons heat more and cool faster.

**Model / method.** Band and DOS (per volume, spin included):

$$\frac{\hbar^2k^2}{2m_e^*}=\mathcal{E}(1+C\mathcal{E})\quad(1),\qquad \rho_e(\mathcal{E})=\frac{1+2C\mathcal{E}}{2\pi^2}\Big(\frac{2m_e^*}{\hbar^2}\Big)^{3/2}\sqrt{\mathcal{E}(1+C\mathcal{E})}\quad(2)$$

$$\frac{\partial f}{\partial t}=\Big(\frac{\partial f}{\partial t}\Big)_{exc}+\Big(\frac{\partial f}{\partial t}\Big)_{e\text{-}ph}+\Big(\frac{\partial f}{\partial t}\Big)_{e\text{-}e}\quad(3)$$

No field acceleration, drift or gradients; no interband transitions; impurity collisions "contribute to the permittivity [but] have no effect on the electron distribution" and are omitted.

*Excitation* (ω = pump frequency; D_J = squared matrix element, explicit form **not stated**: "be it Landau damping, surface/phonon-assisted absorption, etc."):

$$\Big(\frac{\partial f(\mathcal{E})}{\partial t}\Big)_{exc}=A(t)\Big[D_J(\mathcal{E}-\hbar\omega,\mathcal{E})\rho_e(\mathcal{E}-\hbar\omega)f(\mathcal{E}-\hbar\omega)(1-f(\mathcal{E}))-D_J(\mathcal{E},\mathcal{E}+\hbar\omega)\rho_e(\mathcal{E}+\hbar\omega)f(\mathcal{E})(1-f(\mathcal{E}+\hbar\omega))\Big]\quad(5)$$

$$P_{abs}\equiv\int\mathcal{E}\rho_e\Big(\frac{\partial f}{\partial t}\Big)_{exc}d\mathcal{E}=\frac{\omega}{2}\varepsilon_0\varepsilon''|\mathbf{E}(t)|^2\ (6),\qquad A(t)=\frac{\varepsilon_0\varepsilon''|\mathbf{E}(t)|^2}{2\hbar\int D_J(\mathcal{E}-\hbar\omega,\mathcal{E})\rho_e(\mathcal{E})\rho_e(\mathcal{E}-\hbar\omega)f(\mathcal{E}-\hbar\omega,t)(1-f(\mathcal{E},t))d\mathcal{E}}\ (7)$$

E is the local field (taken uniform), ε″ the imaginary part of the ITO permittivity. (5)-(7) conserve electron number and inject exactly P_abs (my check).

*Electron-phonon* (acoustic, deformation potential D, linear dispersion 𝓔_ph = ħv_ph q, Debye cut-off, phonons thermal at T_ph, n_B Bose function). The momentum-space golden-rule term (8) has coupling πD²q²/(ρVω_ph^q). Its angular integration gives (9): ∫d cosθ δ(𝓔_{k±q} − 𝓔_k ∓ ħω_ph) = m_e*/(ħ²kq) if |m_e*v_ph/(ħk) ∓ q/(2k)| ≤ 1, else 0; the authors simplify the condition to q ≲ 2k and obtain

$$\Big(\frac{\partial f(\mathcal{E})}{\partial t}\Big)_{e\text{-}ph}=\frac{D^2}{4\pi\rho(\hbar v_{ph})^4}\sqrt{\frac{m_e^*}{2\mathcal{E}(1+C\mathcal{E})}}\int_0^{\mathcal{E}_{ph}^{max}(\mathcal{E})}\mathcal{E}_{ph}^2\,d\mathcal{E}_{ph}\Big[(n_B+1)\big((1-f(\mathcal{E}))f(\mathcal{E}+\mathcal{E}_{ph})(1+2C(\mathcal{E}+\mathcal{E}_{ph}))-(1-f(\mathcal{E}-\mathcal{E}_{ph}))f(\mathcal{E})(1+2C(\mathcal{E}-\mathcal{E}_{ph}))\big)+n_B\big((1-f(\mathcal{E}))f(\mathcal{E}-\mathcal{E}_{ph})(1+2C(\mathcal{E}-\mathcal{E}_{ph}))-(1-f(\mathcal{E}+\mathcal{E}_{ph}))f(\mathcal{E})(1+2C(\mathcal{E}+\mathcal{E}_{ph}))\big)\Big]\quad(10)$$

(the source prints the last factor as "(1 − f(𝓔+𝓔_ph)f(𝓔)…", a missing parenthesis),

$$\mathcal{E}_{ph}^{max}(\mathcal{E})=\min(\mathcal{E}_D,2\hbar v_{ph}k(\mathcal{E}))\quad(11)$$

Rates are functional derivatives, 1/τ(𝓔) = δ(∂f/∂t)/δf(𝓔):

$$\frac{1}{\tau_{e\text{-}ph}(\mathcal{E})}=\frac{D^2}{4\pi\rho(\hbar v_{ph})^4}\sqrt{\frac{m_e^*}{2\mathcal{E}(1+C\mathcal{E})}}\int_0^{\mathcal{E}_{ph}^{max}}\mathcal{E}_{ph}^2d\mathcal{E}_{ph}\Big[(n_B+f(\mathcal{E}+\mathcal{E}_{ph}))(1+2C(\mathcal{E}+\mathcal{E}_{ph}))+(n_B+1-f(\mathcal{E}-\mathcal{E}_{ph}))(1+2C(\mathcal{E}-\mathcal{E}_{ph}))\Big]\quad(12)$$

$$\frac{1}{\tau_{e\text{-}ph}(\mathcal{E})}\approx\frac{D^2k_BT_{ph}}{4\pi\rho(\hbar v_{ph})^2}\frac{m_e^*}{\hbar k(\mathcal{E})}(1+2C\mathcal{E})\times\{4k^2(\mathcal{E})\ \text{for}\ k<q_D/2;\ \ q_D^2\ \text{for}\ k\ge q_D/2\}\quad(13)$$

Appendix A.1 (phase space): an electron of momentum k exchanges only phonons with q ≲ 2k; all phonons are available only above 𝓔_min,D = (1/2C)[√(1 + C𝓔_D²/(2m_e*v_ph²)) − 1] ≈ 2.46 eV ≈ 2.8 𝓔_F. Appendix A.2 (shape): expanding (10) in 𝓔_ph gives three terms (A1): (𝓔_ph^max)⁴/4 · ∂_𝓔[f(1−f)], (𝓔_ph^max)⁵/5 · ∂²_𝓔 f, and 96(m_e*v_ph²)²𝓔(1+C𝓔)(1+2C𝓔) f(1−f) H(𝓔)H(𝓔_min,D−𝓔), all times the prefactor of (10) × (1+2C𝓔); the third exists only because of the cut-off and "ensures electron number conservation".

*Electron-electron* (static Thomas-Fermi screening):

$$U(\mathbf{q})=\frac{e^2}{4\pi\varepsilon_0\varepsilon_b}\int d^3r\frac{e^{-q_{TF}r}}{r}e^{-i\mathbf{q}\cdot\mathbf{r}}=\frac{e^2}{4\pi\varepsilon_0\varepsilon_b}\frac{1}{q^2+q_{TF}^2}\ (14),\qquad q_{TF}^2=-\frac{4\pi^3e^2}{\varepsilon_0\varepsilon_b}\int d\mathcal{E}\frac{\partial f}{\partial\mathcal{E}}\rho_e(\mathcal{E})\ (15)$$

$$\Big(\frac{\partial f(\mathcal{E})}{\partial t}\Big)_{e\text{-}e}=-\frac{3}{16\hbar^2}\frac{1}{(2\pi)^3}\Big(\frac{e^2}{\varepsilon_0\varepsilon_b}\Big)^2\Big(\frac{m_e^*}{\hbar^2q_{TF}}\Big)^3\frac{\hbar}{\sqrt{2m_e^*(\mathcal{E}+C\mathcal{E}^2)}}\int d\mathcal{E}_1d\mathcal{E}_2d\mathcal{E}_3\,(1+2C\mathcal{E}_1)(1+2C\mathcal{E}_2)(1+2C\mathcal{E}_3)\Big[\frac{qq_{TF}}{q^2+q_{TF}^2}+\tan^{-1}\frac{q}{q_{TF}}\Big]_{q^{min}}^{q^{max}}\delta(\mathcal{E}+\mathcal{E}_1-\mathcal{E}_2-\mathcal{E}_3)\big[(1-f_3)(1-f_2)f_1f-(1-f)(1-f_1)f_2f_3\big]\quad(16)$$

with q^min = max(|k−k₂|,|k₁−k₃|), q^max = min(k+k₂,k₁+k₃). The rate 1/τ_e-e(𝓔) (17) is the same expression without the leading minus sign and with last bracket [(1−f₃)(1−f₂)f₁ + (1−f₁)f₂f₃]. For a Fermi-Dirac f:

$$\frac{1}{\tau_{e\text{-}e}(\mathcal{E})}\approx K\Big\{(\pi k_BT_e)^2+(\mathcal{E}-\mathcal{E}_F)^2+\frac{4}{3}\frac{C(\mathcal{E}-\mathcal{E}_F)}{(1+2C\mathcal{E}_F)^3}(\pi k_BT_e)^2+\cdots\Big\}\quad(18)$$

$$K=\frac{3}{32\hbar}\frac{1}{(2\pi)^3}\Big(\frac{e^2}{\varepsilon_0\varepsilon_b}\Big)^2\Big(\frac{m_e^*}{\hbar^2q_{TF}}\Big)^3\frac{1}{k_F}\Big[\frac{2k_Fq_{TF}}{4k_F^2+q_{TF}^2}+\tan^{-1}\frac{2k_F}{q_{TF}}\Big](1+2C\mathcal{E}_F)^3\quad(19)$$

*TTM reduction.* The "effective" T_e of a non-thermal f is defined by energy matching, ∫𝓔ρ_e f d𝓔 = ∫𝓔ρ_e f^T(𝓔,T_e) d𝓔 (20). Then

$$C_e\frac{dT_e}{dt}=P_{abs}-G_{e\text{-}ph}(T_e-T_{ph})\ (21),\qquad C_{ph}\frac{dT_{ph}}{dt}=G_{e\text{-}ph}(T_e-T_{ph})\ (22)$$

$$\mu(T_e)\approx\mathcal{E}_F-\frac{(\pi k_BT_e)^2}{6}\frac{1}{\rho_e(\mathcal{E}_F)}\frac{d\rho_e}{d\mathcal{E}}\Big|_{\mathcal{E}_F}\quad(24)$$

$$C_e\equiv\frac{d\mathcal{U}_e}{dT_e}\approx\frac{\pi^2k_B^2}{3}T_e[\rho_e+\mu\rho_e^{(1)}]+\frac{7\pi^4k_B^4}{90}T_e^3[3\rho_e^{(2)}+\mu\rho_e^{(3)}]+\frac{d\mu}{dT_e}\Big(\mu\rho_e+\frac{(\pi k_BT_e)^2}{6}[2\rho_e^{(1)}+\mu\rho_e^{(2)}]+\frac{7(\pi k_BT_e)^4}{360}[2\rho_e^{(3)}+\mu\rho_e^{(4)}]\Big)\quad(25)$$

(all ρ_e^{(n)} at μ). G_e-ph is defined by ∫𝓔ρ_e(∂f^T/∂t)_e-ph d𝓔 ≈ −G_e-ph(T_e − T_ph). Limits:

$$G_{e\text{-}ph}=\frac{D^2m_e^{*2}k_B\mathcal{E}_D^4}{16(\pi\hbar)^3\rho(\hbar v_{ph})^4}(1+2C\mu(T_e))^2\ \text{(no momentum conservation)}\ (27),\qquad G_{e\text{-}ph}(T_e\to0)=\frac{D^2m_e^{*2}k_Bk_F^4}{\pi^3\rho\hbar^3}(1+2C\mathcal{E}_F)^2\ (k_F<q_D/2)\ (28)$$

The general integral (26) and its closed form in polylogarithms (A2)-(A8) are long; see flags below.

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| m_e* ; C | 0.3964 m_e ; 0.4191 eV⁻¹ | after Eq. 1 |
| n_e ; 𝓔_F ; band width 𝓔_max | ~1.5×10²⁷ m⁻³ ; 0.88 eV ; 4 eV | Sec. II |
| ε_b (interband background, used in screening) | 4 | after Eq. 15 |
| D ; ρ ; v_ph | 17.2 eV (fitted to room-temperature permittivity) ; 7120 kg m⁻³ ; 6400 m/s | after Eq. 8 |
| 𝓔_D = k_BT_D | ≈ 0.06 eV | Sec. II B |
| k_F ; q_D | 3.54 nm⁻¹ ; 14.3 nm⁻¹ (Au: 11.57 ; 6.85) | Sec. II B |
| C_ph | 2.54×10⁶ J m⁻³ K⁻¹ | Sec. III B |
| γ_e (C_e ≈ γ_eT_e up to ~1300 K; (25) good to ~3500 K) | 12.7 J m⁻³ K⁻² (Au 67.6) | Sec. III B |
| G_e-ph | "~3×10¹⁶" W m⁻³ K⁻¹ in text; 2.5×10¹⁶ for T_e → 0 (Fig. 3c dotted line; Eq. 28 evaluates to the same); ≈7×10¹⁶ at 5000 K (Fig. 3c, read ±10 %); Au 2.5×10¹⁶ | Sec. III B |
| G_e-ph/C_e | 6.6 ps⁻¹ (300 K), 2.19 (1000 K), 1.49 (2000 K), 1.45 (2450 K); Au 1.23 → 0.15 | Fig. 4b |
| μ(T_e) | 0.88 eV → ≈0.50 eV at 5000 K | Fig. 3a |
| 1/τ_e-ph (caption: T_e = 500 K, T_ph = 300 K) | 0.060 fs⁻¹ at 𝓔_F; 0.27 fs⁻¹ plateau above 2.8𝓔_F; C = 0: 0.030; no momentum conservation: 0.22-0.42 | Fig. 1b |
| 1/τ_e-e (same caption) | 5.3×10⁻² (𝓔 = 0), 4.7×10⁻⁴ (𝓔_F), 2.2×10⁻² fs⁻¹ (2𝓔_F); C = 0 lower by ≈1.5 | Fig. 1c |
| K | not quoted | - |

Worked example (Sec. III): I_inc = I₀exp[−4 ln2 (t/τ_p)²], λ = 1300 nm, τ_p = 220 fs, I₀ = 2.5 GW/cm², "bulk ITO"; how I_inc maps to the local |E|² and ε″ is not stated. Results (Fig. 4a): T_e,max = 2464 K at t = +125 fs (BE), 2460 K (TTM); T_e = 2027/2262/1699/900/463 K at 0/250/500/1000/1500 fs (BE); T_ph = 320.9 K at 1650 fs. Au with the same heat source: T_e,max = 1217 K, T_ph = 312.8 K. Text: rise ~2100 K (ITO) vs ~900 K (Au), ratio ≈ √(67.6/12.7). Numerical method (grids, time step) not stated.

**Key results.**
- Momentum conservation lowers 1/τ_e-ph(𝓔_F) ~5× and G_e-ph "~30×" (derived/simulated); G_e-ph then *increases* with T_e, unlike in noble metals.
- Non-parabolicity raises 1/τ_e-ph ~2× and 1/τ_e-e by (1+2C𝓔_F)^{3/2}(1+C𝓔_F)^{−5/4} ≈ 1.5 (derived).
- TTM with T_e-dependent C_e, μ, G reproduces BE temperatures for a 220 fs pulse (simulated).
- f(𝓔) shows steps at 𝓔_F + nħω that e-e scattering smooths within the pulse (simulated).

**Assumptions and limitations.**
- Stated: isotropic f(𝓔); acoustic phonons only; thermal phonons; uniform field; no heat loss to the environment; parameters sample-dependent.
- Consistency checks (all mine; numerical ones use the parameter table, a 4 eV band and my own integration):
  - Eq. 14: the Fourier transform of the Yukawa potential is 4π/(q²+q_TF²), so the printed right-hand side has a spurious 1/4π. The factor (e²/ε₀ε_b)² in (16) is consistent with the correct form.
  - Eq. 15 (Thomas-Fermi wavevector) is dimensionally consistent (m⁻²) but carries an extra 4π³. [16] states "the factor 4π³ in that equation should be removed". Without it q_TF = 3.83 nm⁻¹ at 300 K; with it 42.6 nm⁻¹.
  - Eqs. 16-17 (e-e prefactor) as printed have units of s⁻¹/(J s): one power of ħ too many. Eq. 19 is consistent (J⁻² s⁻¹) and corresponds to the prefactor 3/(16ħ)·(2π)⁻³·(e²/ε₀ε_b)²(m_e*/ħ²q_TF)³·(1/k(𝓔)). With that prefactor and q_TF without 4π³, Eq. 17 reproduces Fig. 1c to a uniform factor 1.14 if T_e = 300 K is used; at the captioned T_e = 500 K the value at 𝓔_F would be 1.1×10⁻³, not 4.7×10⁻⁴ fs⁻¹. K then evaluates to ≈0.065 eV⁻² fs⁻¹.
  - Eq. 12 reproduces Fig. 1b within 5 %.
  - Eq. 10 with the cut-off (11) evaluated at 𝓔 for all four terms does not conserve electron number and gives G_e-ph ≈ 8×10¹⁶ W m⁻³ K⁻¹. Applying the cut-off symmetrically to each pair of states, 𝓔_ph/(ħv_ph) ≤ k(𝓔) + k(𝓔±𝓔_ph), which is what (9) implies before the simplification, conserves number and gives 2.4, 2.65, 3.3, 4.4, 6.7×10¹⁶ at 350, 1000, 2000, 3000, 5000 K, matching Figs. 3c and 4b.
  - Eq. 26 as printed has m_e* where (27), (28) and (A2)-(A5) have m_e*², and does not reduce to (27). Use (10) numerically or (A2)-(A5).
  - γ_e as printed, π²k_B²[ρ_e + μρ_e^{(1)}]/3, evaluates to 26 J m⁻³ K⁻². The quoted 12.7 equals π²k_B²ρ_e(𝓔_F)/3, which (25) gives once its dμ/dT_e term is included.
  - (𝓔_D/2ħv_phk_F)⁴ = 16.4, not "~30"; Fig. 3c also shows a ratio ≈16.

**Use for this project.**
- Adopt (1)-(3), (5)-(7), (10)-(12), (16)-(17) with the corrections above and the parameter table as the WP1 nominal set.
- Validation targets: Fig. 1b,c rate values; G_e-ph(T_e), C_e(T_e) (my exact Fermi-integral values: 3.80×10³, 1.23×10⁴, 2.26×10⁴, 3.98×10⁴ J m⁻³ K⁻¹ at 300, 1000, 2000, 5000 K); T_e(t) of Fig. 4a. That curve coincides with the 2.5 GW/cm² film case of [15] to reading accuracy, so the absorbed power must be taken from [15]'s film calculation.
- No data-availability statement; no code.
- Conflicts with the plan: Umklapp is not mentioned in this paper; the excitation term is a generic, power-normalised matrix element, not specifically phonon- or impurity-assisted absorption; the full e-e integral is given here, not in [15]. The paper's closing inference that falling G_e-ph/C_e explains slower reflectivity recovery at high intensity is revised in [15].
- Open: D_J and the impurity collision rate are not specified ([15] says its e-imp expression is "the same as in" this paper, but arXiv v1 has none; check the published version); optical-phonon scattering is excluded by assumption.

---
### [15] Un et al., "Electronic-based model of the optical nonlinearity of low-electron-density Drude materials" — Phys. Rev. Applied 19, 044043 (2023)

**Source read:** arXiv:2210.08504v3 (30 Mar 2023), full text incl. appendices + TeX source; published version not seen. Figure values come from the vector paths of the authors' figure files; reading uncertainty is about a line width (±10 K on T_e, ±0.5 K on T_ph, ±0.005 on ε′, ±0.001 on ε″, ±0.003 on absorptivity and reflection, ±3 fs).
**Role in this project:** First validation target of the WP1 solver (220 fs and 30 fs cases); the adiabatic-permittivity film model that WP3 and H3 are tested against; source of the TTM/eTTM baselines (M0) and of M1.

**What it does.** Couples the Boltzmann equation (BE) of [14] self-consistently to a time-dependent permittivity and to the field in a 40 nm ITO film, in the adiabatic approximation (model name ANTHEM). Simulates 2.5-75 GW/cm² pumps at 220 fs, two 30 fs cases, and probe reflection. Headline: T_e reaches the Fermi temperature and μ turns negative; Re ε rises by up to ~2 while absorptivity falls from 0.9 to 0.3; a TTM with the [14] parameters reproduces the BE at 220 fs, but at 30 fs the non-thermal part of f is needed to get the rise of ε′ right.

**Model / method.**

$$\frac{\partial f}{\partial t}=\Big(\frac{\partial f}{\partial t}\Big)_{exc}+\Big(\frac{\partial f}{\partial t}\Big)_{e\text{-}ph}+\Big(\frac{\partial f}{\partial t}\Big)_{e\text{-}e}+\Big(\frac{\partial f}{\partial t}\Big)_{e\text{-}imp}\quad(1)$$

The excitation term (2) is Eq. 5 of [14] with A(t) renamed B(t; ω_pump), fixed by

$$\int\mathcal{E}\rho_e(\mathcal{E})\Big(\frac{\partial f}{\partial t}\Big)_{exc}d\mathcal{E}=p_{abs}(t)\ (3),\qquad p_{abs}(t)=\frac{\omega_{pump}}{2}\varepsilon_0\varepsilon''(t;\omega_{pump})|\mathbf{E}(t)|^2\ (6)$$

The real field is written E(t)e^{−iωt} + c.c. with slowly varying envelope E(t). Adiabatic Maxwell equations for the envelope:

$$\nabla\times\mathbf{H}=-i\omega_{pump}\varepsilon_0\varepsilon(t;\omega_{pump})\mathbf{E},\qquad\nabla\times\mathbf{E}=-i\omega_{pump}\mu_0\mathbf{H}\quad(4)$$

$$\varepsilon(t;\omega_{pump})=\varepsilon_\infty+\lim_{\mathbf{q}\to0}\frac{2e^2}{\varepsilon_0q^2}\int\frac{d^3k}{(2\pi)^3}\frac{f_{\mathbf{k}+\mathbf{q}}-f_{\mathbf{k}}}{\mathcal{E}_{\mathbf{k}+\mathbf{q}}-\mathcal{E}_{\mathbf{k}}-\hbar\omega_{pump}-i\hbar(\eta_{\mathbf{k}+\mathbf{q}}+\eta_{\mathbf{k}})/2}\quad(5)$$

$$\varepsilon(t;\omega_{pump})=\varepsilon_\infty-\frac{e^2}{\varepsilon_0m_e^*}\int\frac{\rho_e(\mathcal{E})f(\mathcal{E},t)}{(\omega_{pump}+i\eta(\mathcal{E},t))^2}\frac{1+8C\mathcal{E}(1+C\mathcal{E})/3}{(1+2C\mathcal{E})^3}d\mathcal{E}\quad(9)$$

with η = τ⁻¹_e-e + τ⁻¹_e-ph + τ⁻¹_e-imp (Matthiessen), each rate the functional derivative of its collision term evaluated on the current f. For energy-independent η this is ε = ε_∞ − ω_p²/(ω + iη)² with

$$\omega_p^2[f]=\frac{e^2}{3\pi^2\varepsilon_0m_e^*}\int\Big[\frac{2m_e^*}{\hbar^2}\mathcal{E}(1+C\mathcal{E})\Big]^{3/2}(1+2C\mathcal{E})^{-1}\Big(-\frac{\partial f}{\partial\mathcal{E}}\Big)d\mathcal{E}\quad(10)$$

*Self-consistent loop (Sec. II B):* f(𝓔,t) → η(𝓔,t) → ε(t) by (9) → local pump envelope in the ITO by solving (4) with the transfer-matrix method for the layer stack → p_abs by (6), spatially averaged over the film → B(t) by (3) → excitation term → advance f. The phonon temperature follows C_ph dT_ph/dt = −∫𝓔ρ_e(∂f/∂t)_e-ph d𝓔.

*Probe:* ε(t; ω_pump, ω_probe) is (5) with ħω_probe in the denominator (12); the reflected envelope is r(t; ω_probe) × incident envelope, and R(τ) = ∫|E_probe,ref|²dt / ∫|E_probe,inc|²dt.

*Adiabatic approximation (App. A).* Starts from the non-local-in-time relation D(t) = ε₀ε_∞E(t) + ∫_{−∞}^{t}R(t,t′)E(t′)dt′ for the real fields (A1), with a memory kernel that decays at rate η (A2). Validity conditions: (i) f changes more slowly than η, printed as |(df/dt)_exc| < η/ħ in the text and < η in the footnote (only the footnote form is dimensionally right); (ii) τ_pump ≫ 2π/ω_pump; (iii) τ_pump ≫ 1/η. With η ≈ 0.1 fs⁻¹, (i) holds for local field < 2.7×10⁹ V/m, "an incident intensity of ~200 GW/cm²" in this geometry. Without the approximation, Maxwell + BE needs a time step < 0.2 fs.

*eTTM (App. B).* f = f^T(μ(T_e),T_e) + f^NT:

$$\frac{\partial f^{NT}}{\partial t}=\Big(\frac{\partial f}{\partial t}\Big)_{exc}+\Big(\frac{\partial f^{NT}}{\partial t}\Big)_{e\text{-}e}+\Big(\frac{\partial f^{NT}}{\partial t}\Big)_{e\text{-}ph}\ (B1a),\qquad\frac{d\mathcal{U}^{NT}}{dt}=p_{abs}+\Big(\frac{d\mathcal{U}^{NT}}{dt}\Big)_{e\text{-}e}+\Big(\frac{d\mathcal{U}^{NT}}{dt}\Big)_{e\text{-}ph}\ (B1b)$$

$$C_e\frac{dT_e}{dt}=\Big(\frac{d\mathcal{U}^T}{dt}\Big)_{e\text{-}ph}-\Big(\frac{d\mathcal{U}^{NT}}{dt}\Big)_{e\text{-}e}\ (B2),\qquad C_{ph}\frac{dT_{ph}}{dt}=-\Big(\frac{d\mathcal{U}^T}{dt}\Big)_{e\text{-}ph}-\Big(\frac{d\mathcal{U}^{NT}}{dt}\Big)_{e\text{-}ph}\ (B3)$$

with C_e = ∫𝓔ρ_e ∂f^T/∂T_e d𝓔 (B4), μ(T_e) from number conservation (B5), (∂f^NT/∂t)_e-e = −f^NT/τ_e-e (B6), (dU^T/dt)_e-ph = −G_e-ph(T_e − T_ph), and (dU^NT/dt)_e-ph = (dU/dt)_e-ph − (dU^T/dt)_e-ph computed from the full e-ph term (B9, B10); the last is found negligible. Setting τ_e-e → 0 gives the TTM (B12) = Eqs. 21-22 of [14]. Three temperatures are distinguished: *extracted* (energy-matched to the BE f; called "effective" in [14]), *effective* (TTM), *instantaneous* (eTTM thermal part).

*Gurzhi comparison.* Earlier ITO models used

$$\frac{1}{\tau_{e\text{-}e}(T_e)}=\frac{\omega_{pump}^2}{4\pi^2\omega_p(T_e)}\Big[1+\Big(\frac{2\pi k_BT_e}{\hbar\omega_{pump}}\Big)^2\Big]\quad(B7)$$

with ω_p from (10); it grows super-quadratically because ω_p falls with T_e. The functional-derivative rate of [14] scales as (1+2Cμ(T_e))³[(πk_BT_e)² + (𝓔−μ(T_e))²] and grows only about linearly because μ drops. Fig. 9a, at 𝓔 = μ + ħω_pump/2: this work 0.020, 0.038, 0.052, 0.061 fs⁻¹ at 2000, 4000, 6000, 8000 K; Gurzhi 0.036, 0.100, 0.212, 0.379 fs⁻¹. For e-ph, the earlier formula (B8; Γ₀ = 0.5342 fs⁻¹, T_D = 900 K) gives a similar τ_e-ph < 10 fs, but an RTA for the e-ph loss of non-thermal energy is called "improper": e-ph collisions transfer thermal energy and do not thermalise f.

**Numbers.**
- Material: m_e* = 0.3964 m_e, C = 0.4191 eV⁻¹, C_ph = 2.54×10⁶ J m⁻³ K⁻¹, T_F = 10204 K (Fig. 2 caption), T_D = 700 K (implied by the right axis of Fig. 3). Interband threshold ~3.2 eV. **ε_∞: value not stated. The e-imp collision term and rate: no expression and no value given**; the text says they are "the same as in" [14] and three other references, but arXiv v1 of [14] omits the term.
- Sample (Fig. 1; "similar to" the saturable-mirror experiment of Tirole et al., Phys. Rev. Applied 2022, not the double-slit paper [19]): 40 nm ITO "deposited on glass and covered by a 100 nm gold film"; the schematic draws SiO₂ on the incidence side with a 40 nm label. Pump at 60°, probe at 65°, both p-polarised; whether the angles are in air or glass is not stated. ENZ (ε′ = 0 at 300 K) at ~225 THz. Glass and gold are linear.
- Pulses: Gaussian field envelope exp[−2 ln2 (t/τ)²], τ_pump = τ_probe = 220 fs (intensity FWHM), ω/2π = 230 THz for both, I₀ = |E_pump,0|²/2·√(ε₀/μ₀) is the *incident* peak intensity, 2.5-75 GW/cm². Short-pulse cases: τ_pump = 30 fs with I₀ = 161 GW/cm² (same pulse energy as 22 GW/cm² at 220 fs) and 22 GW/cm²; the probe stays 220 fs. T₀ = 300 K.
- Linear values at 230 THz, 300 K: ε = 0.172 + 0.399i (Fig. 4); absorptivity 0.905 (Fig. 6c); probe reflection 0.071 (Fig. 7a).

Validation table, 220 fs pump (BE/ANTHEM; t = 0 is the pump peak; "n.p." = not plotted):

| I₀ (GW/cm²) | T_e,max (K) @ t (fs) | TTM T_e,max (K) | T_ph(990 fs) / T_ph,max (K) | ε′_max (Δε′) @ t | ε″_max (Δε″) @ t | Absorptivity | Probe R_max @ τ | R decay rate (ps⁻¹) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.5 | 2464 @ 124 | 2472 | 319.2 / 321.0 | n.p. | n.p. | 0.777 | 0.237 @ 160 | 3.46 |
| 5 | 3313 @ 117 | 3315 | 332.5 / 336.2 | 0.765 (0.593) @ 130 | 0.495 (0.096) @ 163 | 0.668 | 0.378 @ 154 | 3.00 |
| 10 | 4348 @ 109 | 4331 | 352.7 / 359.0 | n.p. | n.p. | 0.545 | 0.516 @ 146 | 2.57 |
| 22 | 5902 @ 105 | 5850 | 389.6 / 399.6 | 1.476 (1.305) @ 118 | 0.547 (0.149) @ 351 | 0.418 | 0.645 @ 151 | 2.09 |
| 50 | 8393 @ 108 | 8247 | 456.3 / 472.3 | 2.013 (1.842) @ 125 | 0.594 (0.195) @ 574 | 0.318 | 0.734 @ 152 | 1.75 |
| 75 | 10262 @ 115 | 10003 | 507.2 / 528.2 | 2.315 (2.143) @ 134 | 0.630 (0.231) @ 707 | 0.281 | n.p. | n.p. |

Sources: Figs. 2a, 3, 4a-b, 6, 7. T_ph,max is the paper's T₀ + U_abs/C_ph. At 50 and 75 GW/cm² ε″ has an early local maximum (0.517 at −37 fs; 0.515 at −75 fs) and a dip (0.513 at +54 fs; 0.491 at +77 fs) before the late maximum. Further: T_e(990 fs) = 914, 1273, 1641, 2083, 2680, 3108 K; 10-90 % rise of ε′ = 203-232 fs; μ/𝓔_F at 50 GW/cm² is negative from +24 to +218 fs, minimum −0.195; ω_p²/ω_p²(300 K) falls to 0.554 at 50 GW/cm²; dT_e/dt = −2.4 K/fs at 2000 K and −10.0 K/fs at 6000 K (T_ph = 300 K). The decay rate is 1/(τ_open − τ_peak) of the markers in Fig. 7a; the open "half-maximum" markers lie 53-61 % of the way from the pre-pump reflection up to the peak, so the baseline is not the pre-pump value (my reading). Extrapolated T_ph,max = 680, 1061, 1821, 2582 K at 125, 250, 500, 750 GW/cm², crossing the 1900 K damage line near 525 GW/cm².

30 fs pump (Figs. 8, 10, 11):

| Case | T_e,max (K) @ t | TTM T_e,max | T_ph(450 fs) (K) | ε′_max (Δε′) @ t | 10-90 % rise of ε′ | ε″_max (Δε″) @ t | Probe R_max @ τ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 161 GW/cm² | 7470 @ 26 fs | 7160 | 366.0 | 1.717 (1.545) @ 51 fs | 37 fs | 0.546 (0.147) @ 354 fs; ≈0.501 at 0-20 fs | 0.653 @ 170 fs |
| 22 GW/cm² | 3187 @ 27 fs | 3161 | 315.3 | 0.682 (0.510) @ 39 fs | 33 fs | 0.482 (0.083) @ 69 fs | 0.296 @ 130 fs |
| 22 GW/cm², 220 fs | 5902 @ 106 fs | 5850 | 368.2 | 1.476 @ 118 fs | 211 fs | 0.547 @ 354 fs | 0.645 @ 151 fs |

Peak T_ph for the 30 fs cases is not shown (curves end at 450 fs). No decay-rate markers are given; with a half-way-to-pre-pump criterion I read 1.9, 3.4 and 1.8 ps⁻¹ for the three rows (the authors' own marker gives 2.09 for the third). 161 GW/cm²: time-resolved absorption falls 0.905 → 0.261 within the pulse (Fig. 10c); f(𝓔_F + ħω/2) peaks at 0.130 (0.085 for 22 GW/cm²). The eTTM instantaneous T_e peaks at 6957 K at +67 fs versus extracted 7470 K at +26 fs; with a thermal-only permittivity the eTTM overshoots (T_e > 10⁴ K per the text).

**Key results.**
- T_e,max ≈ T_F at 75 GW/cm²; μ < 0 for I₀ ≥ 50 GW/cm² (simulated).
- Δε′ tracks T_e and hence ∫p_abs dt; ε′ peaks 118-134 fs after the peak of a 220 fs pump and rises 12-fold at 50 GW/cm². ε″ rises less and non-monotonically (τ⁻¹_e-e up, ω_p² down, then phonon heating) (simulated).
- Absorptivity falls because the ENZ resonance shifts; ITO is "not a saturable absorber" (simulated/argued).
- Probe-reflection recovery slows with intensity (3.5 → 1.75 ps⁻¹) although |dT_e/dt| grows, because the ENZ shift exceeds the resonance width; stated to match the saturable-mirror experiment of Tirole et al. quantitatively (simulated; no overlay with data shown).
- For 30 fs pumps, rise of ε′ is pulse-limited only if f^NT enters the permittivity; with a thermal-only permittivity it is τ_e-e-limited (tens of fs) (simulated).
- The [14] parameters explain the ad hoc fixes in the model of Baxter et al. [9] (G_e-ph reduced 8×, damping taken linear in T_e, no e-ph loss from non-thermal electrons) (argued).

**Assumptions and limitations.**
- Stated: adiabatic, quasi-monochromatic; uniform absorbed power across the film; large spot; no interband, multiphoton, emission or photoemission channels, hence limited to ≤ 75 GW/cm² at 220 fs; ε_∞ constant.
- (my reading) Eq. 6 is the usual cycle-averaged Joule power only if |E| is the peak amplitude of the real field; with the envelope convention above the envelope is half of that.
- (my reading) ε ∝ 1/(ω+iη)² corresponds to a Drude rate γ ≈ 2η. A fit of this form to the 300 K curve in the Fig. 1 inset gives ε_∞ = 4.00, ħω_p = 1.87 eV, η = 0.075 fs⁻¹ (rms 0.002). Eq. 10 with n = 1.5×10²⁷ m⁻³ and the stated m_e*, C gives ħω_p = 1.73 eV, which would put ENZ at 225 THz only for ε_∞ ≈ 3.45. This 8 % gap in ω_p is unexplained. Comparing η with 1/τ_e-ph(𝓔_F) = 0.060 fs⁻¹ from [14] suggests an e-imp rate of order 0.01-0.02 fs⁻¹.
- (my reading) T_ph,max − 300 K equals absorptivity × I₀τ cos60°/(d·C_ph) to 0.1 K, so the pulse energy is normalised as I₀τ, without the Gaussian factor 1.064; the top axis of Fig. 6 uses yet another scale (12.1 eV/nm³ at 75 GW/cm²). Normalisation-free targets are the absorbed energy densities U_abs = C_ph(T_ph,max − T₀) = 0.333, 0.574, 0.935, 1.58, 2.73, 3.62 eV/nm³ for 2.5-75 GW/cm².
- Text versus figures: 161 → 22 GW/cm² at 30 fs is said to halve ΔT_e and Δε′; the figure gives factors 0.40 and 0.33.
- Numerical method (energy grid, time step) and D_J are not stated.

**Use for this project.**
- Use the two tables as WP1/WP3 acceptance targets. Test the kinetics first: drive the solver with a p_abs(t) whose time integral is U_abs and compare T_e(t) and T_ph(t); the time profile of p_abs is published only for the 161 GW/cm² case (absorption versus time, Fig. 10c). The optical targets (ε(300 K), absorptivity versus I₀, ε(t)) additionally need ε_∞, the e-imp rate, D_J and the stack details, none of which is given.
- Ask the authors for these and obtain the published version: [16] cites an "Eq. (B6)" of this paper containing a q_TF-dependent bracket, which does not exist in arXiv v3. No data-availability statement.
- The plan says [14, 15] attribute the real-permittivity change to non-thermal electrons. That is too strong: at 220 fs the thermal (TTM) permittivity matches the BE; the non-thermal part matters for the rise under ≤ 30 fs pumps. "Instantaneous" here means ε′ follows the absorbed energy without extra delay, not a 1-10 fs rise under a 220 fs pump.
- H3 baseline: the probe model here is r(t) at the carrier frequency times the incident envelope, in a uniformly heated film.

---
### [16] Un et al., "Quantum-optical theory of the few-femtosecond nonlinear optical response of Drude metals with a nonparabolic conduction band" — Phys. Rev. B 113, 165411 (2026)

**Source read:** arXiv:2601.10985v2 (26 Feb 2026), full text + TeX source; figure values from the vector paths of the PDF (±1 THz, ±0.01). The arXiv version lists three authors, Ieng-Wai Un, Subhajit Sarkar and Yonatan Sivan (the project plan cites Un alone), and spells the title "few femtosecond … non-parabolic". The journal reference could not be checked against these files.
**Role in this project:** Tier-2 stretch model of WP1 (coherent level for few-cycle pumps; state two-dimensional in energy); it also corrects one formula each of [14] and [15] and shows what the cycle-averaged Boltzmann description leaves out.

**What it does.** Derives equations of motion for the conduction-electron density matrix in energy space for a Kane band: populations f(𝓔,t) plus coherences between states of different energy, driven by the real field with no rotating-wave or cycle averaging. Applies them to an ITO nanosphere in the quasi-static limit under a 5 fs, 175 THz pulse. This paper concentrates on the formulation and on coherence signatures; strength and speed of the permittivity response are in a companion Letter (its Ref. 1, Un, Sarkar, Sivan, "Origin, strength, and speed of the nonlinear optical response of transparent conducting oxides to single-cycle optical pulses"). The Letter is reference [62] of this library.

**Model / method.** Start (momentum space): reduced density matrix ρ with ρ_kk = f_k, dipole coupling −d_z E(t), Kane dispersion ħk(𝓔) = √(2m_e*𝓔(1+C𝓔)) (12). Dipole elements follow from momentum elements of a one-dimensional infinitely deep well of width L, (d_z)_kk′ = −2ieħ (p_z)_kk′/[(p²)_kk − (p²)_k′k′] (6). Relaxation operator: populations relax by the same e-e, e-ph and e-imp collision integrals as in [14, 15] (7a); coherences decay in a relaxation-time form corrected to preserve the continuity equation,

$$\mathcal{R}_{\mathbf{k}\mathbf{k}'}=-\eta_{\mathbf{k}\mathbf{k}'}\big(\rho_{\mathbf{k}\mathbf{k}'}-h_{\mathbf{k}\mathbf{k}'}\rho_{\mathbf{k}'\mathbf{k}}\big),\qquad \eta_{\mathbf{k}\mathbf{k}'}=(\eta_{\mathbf{k}}+\eta_{\mathbf{k}'})/2,\quad \eta_{\mathbf{k}}=\delta\mathcal{R}_{\mathbf{k}\mathbf{k}}/\delta f_{\mathbf{k}}\quad(8)$$

(h is a phase factor, equal to 1 here). The field-driven coupling between different coherences (last term of Eq. 5b) is dropped, i.e. "coherent optical nonlinearities are negligible". After angular averaging, with ρ_e(𝓔) = (k²/π²)(dk/d𝓔) and V the particle volume, the state is f(𝓔,t) and the dipole-weighted coherence \overline{d_{𝓔′𝓔}ρ_{𝓔𝓔′}}(𝓔,𝓔′,t) ≡ (4π)⁻²∫dΩ_k dΩ_k′ (d_z)_k′k ρ_kk′:

$$\Big(\frac{\partial f(\mathcal{E},t)}{\partial t}\Big)_{photon\text{-}e}=-\frac{iV}{2\hbar}E(t)\int\rho_e(\mathcal{E}')\big(\overline{d_{\mathcal{E}'\mathcal{E}}\rho_{\mathcal{E}\mathcal{E}'}}-\overline{d_{\mathcal{E}\mathcal{E}'}\rho_{\mathcal{E}'\mathcal{E}}}\big)d\mathcal{E}'\quad(13)$$

$$\frac{\partial(\overline{d_{\mathcal{E}'\mathcal{E}}\rho_{\mathcal{E}\mathcal{E}'}})}{\partial t}+\eta_{\mathcal{E}\mathcal{E}'}\big(\overline{d_{\mathcal{E}'\mathcal{E}}\rho_{\mathcal{E}\mathcal{E}'}}-\overline{d_{\mathcal{E}\mathcal{E}'}\rho_{\mathcal{E}'\mathcal{E}}}\big)+i\Big(\frac{\mathcal{E}-\mathcal{E}'}{\hbar}\Big)\overline{d_{\mathcal{E}'\mathcal{E}}\rho_{\mathcal{E}\mathcal{E}'}}=-\frac{i}{\hbar}|D(k(\mathcal{E}),k(\mathcal{E}'),L)|^2E(t)\big[f(\mathcal{E},t)-f(\mathcal{E}',t)\big]\quad(14)$$

$$P(t)=\frac{V}{4}\int\rho_e(\mathcal{E})\rho_e(\mathcal{E}')\big(\overline{d_{\mathcal{E}'\mathcal{E}}\rho_{\mathcal{E}\mathcal{E}'}}+\overline{d_{\mathcal{E}\mathcal{E}'}\rho_{\mathcal{E}'\mathcal{E}}}\big)d\mathcal{E}\,d\mathcal{E}'\quad(15)$$

The full population equation is ∂f/∂t = (13) + (∂f/∂t)_e-e + (∂f/∂t)_e-ph + (∂f/∂t)_e-imp. Here η_𝓔𝓔′ = (η_𝓔 + η_𝓔′)/2 and η_𝓔 = η_e-e + η_e-ph + η_e-imp (functional derivatives, Matthiessen). The squared dipole element, angle-averaged (units C² m²):

$$|D(k,k',L)|^2=\frac{2e^2(2\pi L)^2}{(kL)(k'L)(k^2-k'^2)^2LV}\Big\{\frac{((kL)^2-(k'L)^2)^2-16}{4}\ln\Big(\frac{4+|(kL)^2-(k'L)^2|}{4+(kL-k'L)^2}\Big)-\big[((kL)^2-(k'L)^2)^2+1\big]\ln\Big(\frac{1+|(kL)^2-(k'L)^2|}{1+(kL-k'L)^2}\Big)-\frac{3}{4}\big[((kL)-(k'L))^2\big]^2\ln\Big|\frac{k-k'}{k+k'}\Big|\Big\}\quad(16)$$

(as printed). It is concentrated near 𝓔 = 𝓔′ with a logarithmic divergence there. Field in a sphere of radius a ≪ λ₀ in a host ε_h (quasi-electrostatic, ∂B/∂t ≈ 0; field and polarisation uniform inside):

$$\mathbf{E}_{NP}(t)=\frac{1}{\varepsilon_\infty+2\varepsilon_h}\Big[3\varepsilon_h\mathbf{E}_{inc}(t)-\frac{\mathbf{P}_{NP}(t)}{\varepsilon_0}\Big]\quad(17)$$

with P_NP from (15). Weak-field CW limit:

$$\varepsilon(\omega,T_0)=\varepsilon_\infty+\frac{V}{2\varepsilon_0\hbar}\int d\mathcal{E}\,d\mathcal{E}'\frac{\rho_e(\mathcal{E})\rho_e(\mathcal{E}')(\mathcal{E}-\mathcal{E}')|D|^2(f(\mathcal{E},T_0)-f(\mathcal{E}',T_0))}{\omega^2-(\mathcal{E}-\mathcal{E}')^2/\hbar^2-2i\eta_{\mathcal{E}\mathcal{E}'}\omega}\quad(18)$$

Diagnostics: P_photon-e(t) = ∫₀^{𝓔_max}ρ_e 𝓔 (∂f/∂t)_photon-e d𝓔 (19), stated to equal E_NP·∂P_NP/∂t, and U(t) = ∫_{−∞}^{t}P_photon-e dt′ (20). Initial state: 300 K Fermi-Dirac populations, zero coherences.

*Coherent versus incoherent.* Coherent: photon-electron coupling through the coherences, which carry phase and oscillate at ω₀; absorption and stimulated emission both arise from (13)-(14), and energy conservation holds without the Poynting normalisation A(t) of [14, 15], which the authors call "ad hoc". Incoherent: the three collision integrals acting on f, and dephasing of coherences at rate η. Neglected: spontaneous emission (App. A: ~10⁶ times weaker), interband transitions, "electron- and phonon-assisted absorption events" (absorption is of Landau-damping, finite-size type through D), and the direct field-driven coupling between coherences (the authors' "coherent optical nonlinearities").

*Relation to the Boltzmann model of [14, 15].* The paper does not derive the reduction; it refers to Sarkar, Un, Sivan, "Comparative study of models of optically-induced electron non-equilibrium in Drude metals" (2024, no journal given). The collision integrals are identical by construction. (my reading) For a quasi-monochromatic field with envelope slow compared with 1/η, eliminating the coherences from (14) gives a golden-rule rate ∝ |E|²|D|²[f(𝓔′) − f(𝓔)] with Lorentzian energy conservation, i.e. the excitation term of [14] with D_J ↔ |D|² and the prefactor fixed by the field; (18) is the corresponding Drude/Lindhard permittivity. That limit is the adiabatic, cycle-averaged model of [15].

**Numbers.**
- Material: m_e* = 0.3964 m_e, C = 0.4191 eV⁻¹; otherwise "the same material parameters as" [14, 15]. n_e, 𝓔_F, ε_∞ not stated; the text uses "an electron density of 1/nm³" in an estimate ([14]: 1.5 nm⁻³).
- Geometry: single ITO nanosphere, quasi-static; "the 4 nm sphere" with "33 electrons", i.e. 4 nm diameter (33.5 nm³); L = 4 nm in (16) (Fig. 1). This confirms the plan's description, with the caveat that the matrix elements come from a 1D box, not a sphere. Host ε_h not stated; the plotted |3ε_h/(ε_ITO+2ε_h)| is reproduced exactly by ε_h = 1 (my check).
- Nanoparticle permittivity at 300 K (Fig. 2b, Eq. 18): −1.6 + 1.4i at 175 THz; ε′ = 0 at 215 THz; ε″ larger than bulk because of finite size. Local-field factor peaks at 2.08 near 180 THz.
- Pulse: E_inc = ẑE₀exp[−2 ln2 (t/τ_p)²]cos(ω₀t), τ_p = 5 fs, ω₀/2π = 175 THz (0.723 eV, cycle ≈ 6 fs), E₀ = 0.6, 1.2, 2.4 V/nm incident peak field (Fig. 8c also 0.4, 0.8, 1.8). In vacuum these are 0.048, 0.19 and 0.76 TW/cm² (my conversion).
- Term sizes at E₀ = 0.6 V/nm (Fig. 3): photon-e term up to ±0.4 fs⁻¹ near 𝓔_F; e-e "an order of magnitude weaker", e-ph "a further order of magnitude smaller"; e-e becomes appreciable only for t ≳ 15 fs.
- E₀ = 1.2 V/nm (Fig. 6): P_photon-e oscillates between +0.33 and −0.14 eV nm⁻³ fs⁻¹; U(t) rises in steps from about −3 fs to +8 fs and settles at 0.63 eV/nm³.
- Local field versus E₀ ≈ 0, 0.4, 0.6, 0.8, 1.2, 1.8, 2.4 V/nm (Fig. 8c): spectral peak 179, 178, 177, 176, 173, 170, 168 THz; spectral width 66, 70, 74, 79, 90, 105, 121 THz.
- Phase delay of E_NP relative to its linear reference at t = +10 fs (Fig. 10b): 0.05π, 0.20π, 0.65π for 0.6, 1.2, 2.4 V/nm. Phase of E_NP relative to E_inc: between about −0.2π and −0.4π during the pulse; quasi-static CW value −0.4π.
- Level spacing in a few-nm particle ≲ 0.02-0.05 eV, comparable with η; band width ≳ 4 eV.
- Numerical cost: no grid size, time step or run time is given. Only a scaling statement: the e-e term costs O(N_k⁶) on an N_k³ momentum grid and O(N_E²) in energy space.

**Key results.**
- During a 5 fs pulse the distribution is driven almost entirely by the coherent photon-electron term; collisions act afterwards (simulated).
- Harmonic content: (∂f/∂t)_photon-e, f and the absorbed power oscillate at ~2ω₀ (period ~3 fs). With rising field this population pulsation drives ±3ω₀, ±5ω₀ components in the coherences, hence in the polarisation, and 4ω₀, 6ω₀ in the populations (Fig. 4, shown for the photon-e term) (simulated). The Outlook calls this "repeated frequency-doubling features in the polarization". No emitted-harmonic spectra or efficiencies are given.
- P_photon-e is repeatedly negative, but U(t) shows a net stepwise increase and "resembles the permittivity dynamics" of the Letter (simulated; formal link left open).
- Stimulated absorption grows super-linearly with |E₀|² because already-excited electrons absorb again (sequential excited-state absorption); stimulated emission scales as |E₀|⁴; the net absorption grows sub-linearly (simulated; scaling plots are in the Letter).
- Sign: the local-field spectrum red-shifts and broadens, "consistent with the expected increase in the real part of the ITO permittivity". No sign reversal is reported up to 2.4 V/nm.
- The phase lag builds on the trailing edge, not at the peak; the authors read this as an absorptive rather than Kerr-type nonlinearity and propose it as an experimental test (simulated/argued).
- Response time: no rise-time number is given here. The absorbed energy accumulates within the pulse, over roughly 10 fs (my reading of Fig. 6b).

**Assumptions and limitations.**
- Stated: quasi-static sphere; 1D-box matrix elements; relaxation-time dephasing; no interband transitions; no coherent wave mixing. A footnote notes that the e-e term in the population equation contains only normal processes, that these do not relax the polarisation after angular averaging, that Umklapp processes do, and that the ratio for ITO "is presently unknown; however, it is likely to be of order unity". How this enters η_e-e numerically is not stated.
- (my reading) (14) is a damped oscillator with damping 2η, so the Drude rate of (18) is γ = 2η_𝓔𝓔′.
- Eq. 21 (decomposition of f(𝓔) − f(𝓔′) into absorption/emission in/out) is printed with four plus signs; two must be negative.
- The text says "a few eV are absorbed per nm³" and then 0.03-0.15 eV per electron for 33 electrons in 33 nm³. The two statements differ by a factor ~33; Fig. 6b gives 0.63 eV/nm³ at 1.2 V/nm, about 0.4-0.6 eV per electron (my arithmetic).
- (my reading) The O(N_E²) cost must be per energy point; the full e-e term is then O(N_E³), as the plan assumes. The coherence array adds N_E² complex unknowns stepped at sub-cycle resolution.

**Use for this project.**
- Corrections to adopt in WP1: footnote to Eq. 7: "the factor 4π³ in [Eq. 15 of [14]] should be removed"; a second footnote corrects a bracket in "Eq. (B6)" of [15] to (1 + q_TF²/2k²)² − 1; no such bracket exists in arXiv v3 of [15], whose B6 is the relaxation-time relation.
- A block of the TeX source excluded from compilation writes the e-e rate with prefactor (3/4)(2π/ħ)((2π)²/2)(2π)⁻⁶(m_e*/ħ²)³(e²/ε₀ε_b)²/(2q_TF³)·(1/k) = 3/(16ħ(2π)³)…, which supports the dimensionally corrected form of Eq. 17 of [14]. It is not part of the paper; treat as a hint. A source comment also remarks that the effective mass "so far is too high (0.4 instead of 0.26)".
- Tier-2 implementation: state (f, coherence array on an 𝓔×𝓔′ grid), equations (13)-(15) plus the WP1 collision terms; validation targets are Figs. 6 and 8c.
- Future work as listed by the authors: other "particle sizes and shapes"; other "material parameters or TCOs (effective mass, non-parabolicity, electron density, ...)"; "different geometries (primarily, of layers, which is the geometry of choice in most experiments to date)"; other pump wavelengths (800 nm); "probe wavelength and their duration and chirp"; comparison with "more computationally accessible approaches (e.g., the hydrodynamic approach)" to study "complex spatio-temporal dynamics of single cycle pulses in TCOs"; adding interband transitions through a valence-band block; a formal link between U(t) and the permittivity. Films are named explicitly; propagation only through the "spatio-temporal dynamics" phrase. There is no probe pulse in this paper.
- No data-availability statement.
- The sign reversal of Ozlu et al. [23] is not addressed: interband channels are absent.

**Companion Letter.** Published as Un, Sarkar, Sivan, Nanophotonics 15, e70272 (2026); added to the library as **[62]** and digested there (2026-11-18). It supersedes the short note on its arXiv version that stood here.
---

### [62] Un, Sarkar & Sivan, "Origin, strength, and speed of the nonlinear optical response of transparent conducting oxides to single-cycle optical pulses" — Nanophotonics 15, e70272 (2026); arXiv:2608.11756

**Source read:** published version (Wiley, CC BY 4.0, received 7 June, accepted 18 August 2026), all 11 pages incl. figures, and the Supporting Information (9 pages), both supplied by the project lead on 2026-11-18; figure values read by eye from the PDF (±10 %). Not in revision 2 of the plan; added as [62]. The arXiv v1 text (12 Aug 2026) read the day before is the same paper.
**Role in this project:** the only microscopic, collision-resolved calculation of TCO dynamics under a single-cycle pump; source of the authors' current ITO parameter set (Table S1), of a generalised Fermi-liquid e–e rate, and of the "thermal model driven by the instantaneous absorbed power" that defines the cheap few-cycle path (plan Tier 2a). Companion of [16], which holds the formalism.

**What it does.** Solves the density matrix (populations and coherences) of the conduction electrons of a 4 nm ITO sphere in air, quasi-statically coupled to its own polarization, under a 5 fs, 175 THz pump, with the rigorous e–e (screened Coulomb), e–ph (deformation potential) and e–impurity (Brooks–Herring) collision integrals of [14, 15] acting on the populations and their functional derivatives as decoherence rates. A weak 5 or 20 fs probe is simulated to extract a time-dependent "instantaneous permittivity". The paper asks how strong the response is, where it comes from, and how fast it turns on and off.

**Model / method.**

Equations of motion (SI Eq. S1; d_kk′ dipole elements of a 1D infinite well of the particle size, h_kk′ a phase factor that keeps charge continuity, η_kk′ = (η_k + η_k′)/2):

$$\partial_t f_k=-\tfrac{i}{\hbar}E_{NP}(t)\sum_{k'}\big[d_{k'k}\rho_{kk'}-d_{kk'}\rho_{k'k}\big]+\big(\partial_tf_k\big)_{e\text{-}e}+\big(\partial_tf_k\big)_{e\text{-}ph}+\big(\partial_tf_k\big)_{e\text{-}imp},$$

$$\partial_t\rho_{kk'}=-\eta_{kk'}\big(\rho_{kk'}-h_{kk'}\rho_{k'k}\big)-\tfrac{i}{\hbar}(\mathcal E_k-\mathcal E_{k'})\rho_{kk'}+\tfrac{i}{\hbar}E_{NP}(t)d_{kk'}\big[f_{k'}-f_k\big],$$

solved in energy space after angular integration ([16]). Local field (main Eq. 1): E_NP = [3ε_h E_inc − P_NP/ε₀]/(ε∞ + 2ε_h), P_NP = (2/V)Σ d_k′k ρ_kk′ (S10). Collision terms: e–e (S2) with M(q) = e²/[ε₀ε∞(q² + q_TF²)] and q_TF² = −(e²/ε₀ε∞)Σ∂f/∂𝓔 (S3–S4; a footnote repeats that the 4π³ of [15] Eq. 15 is an error); e–ph (S5) with coupling πD²q²/(ρVω_ph), linear phonons, Debye cut-off; e–impurity (S6), Brooks–Herring. Rates η = functional derivatives (S7–S9). The linear permittivity (S11) and the "ANTHℰM" instantaneous permittivity (S12) of the sphere have the Lindhard-like form Σ 2|d|²(𝓔_k − 𝓔_k′)(f_k′ − f_k)/[(𝓔_k − 𝓔_k′)² − (ħω)² − 2iħ²ωη_kk′]. Coherences are split into absorption and stimulated-emission parts (S13–S16) to separate the two.

Thermal model (main Eq. 2), the paper's simple alternative:

$$C_e(\tilde T_e)\,\frac{d\tilde T_e}{dt}=-G_{e\text{-}ph}(\tilde T_e-T_{ph})+\dot{\mathbf P}_{NP}(t)\cdot\mathbf E_{NP}(t),$$

with the instantaneous (not cycle-averaged) net absorbed power, C_e from [14] Eq. 25, G from [14] Eq. 26; the coherences and P_NP are then computed from the density-matrix equation with the thermal f(T̃_e(t)) (footnote 2). For long pulses Ṗ·E averages to ω₀ε″|E|²/2, the absorption of [15].

Generalised Fermi-liquid e–e rate (S24), derived by expanding only the Fermi-Dirac product of the thermal rate (S21) to second order in (𝓔 − μ)/k_BT̃_e while keeping the screened matrix element exact and the lower band edge finite:

$$\eta^T_{e\text{-}e}(\mathcal E,\tilde T_e)\approx K(\tilde T_e)\Big[A_0(k_B\tilde T_e)^2+A_1\,k_B\tilde T_e(\mathcal E-\mu)+A_2(\mathcal E-\mu)^2\Big],$$

with K the T → 0 coefficient of [14] Eq. 19 (S25) and A_n(𝓔, T̃_e) double integrals (S26); A₀ → π², A₁ → 0, A₂ → 1 in the textbook limit. T̃_e is defined by energy matching (S20).

As printed, five formulas need corrections (details and checks in `docs/theory/un2026_notes.md`): the Fermi-Dirac product (S22) lacks a factor 1/4; (S23) has tanh((𝓔₁ − 𝓔₂)/k_BT̃_e) where (S26) and the expansion need tanh((𝓔₁ − 𝓔₂)/2k_BT̃_e); (S26) lacks 1/(k_BT̃_e)²; (S25) lacks the 1/(2π)³ of (S19); the bracket of (S9) reads 1 − q_TF²/2k² where 1 + q_TF²/2k² is needed (the printed form is negative for q_TF < 2k). The population equation (S2) carries a factor ½ that the rate (S7) does not, so the two differ by 2 as printed. Screening uses ε∞ in (S3–S4) but ε_b in (S6), (S9), (S18), (S25). ρ is printed in kg/m².

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| System | ITO sphere, 4 nm diameter, in air (ε_h = 1), quasi-static; dipole resonance ≈180 THz (≈0.7 eV) | Sec. 3, Fig. 1 |
| Pump | E_inc = E₀exp[−2ln2(t/τ_p)²]cos ω₀t, τ_p = 5 fs (intensity FWHM), ω₀/2π = 175 THz (0.72 eV, 1.71 µm, cycle 5.7 fs) | Sec. 3 |
| Fields | E₀ = 0.6, 1.2, 1.8, 2.4 V/nm incident; 1.2 V/nm "comparable to the field levels used in" [20] and Herbst et al. 2026 | Sec. 4.1 |
| Intensity convention | "1 TW/cm²" for E₀ = 1.8 V/nm, i.e. I = E₀²/Z₀ (0.86 TW/cm²); with I = E₀²/2Z₀ it is 0.43 TW/cm² (my calculation) | Sec. 4.2 |
| Table S1 | C = 0.4191 eV⁻¹; m* = 0.3964 mₑ; 𝓔_F = 0.8793 eV; n ≈ 1.5×10²⁷ m⁻³; **D = 8.6 eV** ("obtained by fitting the ITO permittivity at the bulk limit and at room temperature" of Wang et al. 2019, plan ref. [7]); ρ = 7120 (printed kg/m²); v_ph = 6400 m/s; 𝓔_D = 0.06 eV; **n_imp = 0.05 nm⁻³, Z = 1** | SI Table S1 |
| Not printed | ε∞, ε_b, band top 𝓔_max, energy grid, time step, T̃_e(t) for any run | — |
| ε(ω) of the sphere at 300 K | ε′ ≈ −3.2, ε″ ≈ 2.4 at 150 THz; ε′ = 0 near 240 THz; ε′ ≈ +1.6, ε″ ≈ 0.3 at 300 THz (by eye) | Fig. 1b |
| Instantaneous ε at 160 THz, E₀ = 1.2 V/nm | ε′: −2.6 → −1.0; ε″: 1.85 → 2.25; rise between ≈ −5 and +8 fs with a ≈3 fs ripple; thermal model and ANTHℰM nearly coincide (by eye) | Fig. 4a,b |
| Δε′ (max) vs absorbed energy U_abs | ≈0.4, 1.5, 2.7, 3.9 at U_abs ≈ 0.02, 0.2, 0.5, 1.05×10⁹ J/m³ (sublinear); without stimulated emission linear, ≈4 at 0.27×10⁹ J/m³ (by eye) | Fig. 4c |
| Δε′ vs E₀ | ≈0.3, 0.7, 1.0, 1.5, 2.6, 3.9 for E₀ ≈ 0.4–2.4 V/nm; with excited-state absorption but no stimulated emission 2.7, 6.9, 13.3 at 1.2, 1.8, 2.4 V/nm (by eye) | Fig. 4d |
| e–e rate | 1/η_e-e ≈ 12 fs, nearly energy independent, at T̃_e ≈ 4500 K ("1 TW/cm²", E₀ = 1.8 V/nm, t = 5 fs); η_e-e(𝓔_F) ≈ 0.005 → 0.17 fs⁻¹ for E₀² ≈ 0.16 → 5.8 V²/nm², roughly ∝ E₀²; at 𝓔_F + 1 eV ≈ 0.03 → 0.18 fs⁻¹ | Sec. 4.2, Fig. 5a,b, S1a |
| A₀, A₁, A₂ at 𝓔 = 𝓔_F | A₀ ≈ 9.4, 8.6, 6.3, 4.9, 4.0; A₁ ≈ −0.1, −0.3, −0.75, −0.85, −0.85; A₂ ≈ 0.95, 0.85, 0.55, 0.42, 0.35 for T̃_e = 300, 1000, 3000, 5000, 7000 K (by eye) | Fig. S1b–d |
| e–ph rate | ≈0.001 (0.1𝓔_F), ≈0.01 (𝓔_F), ≈0.06 fs⁻¹ (3𝓔_F), the same at −5, 0, +5 fs | Fig. 5c |
| Simulation window | ≈ −10 to +20 fs | Fig. 2 |

**Key results.**

- Simulated: the permittivity change is of order one (Δε′ up to ≈4, "∼250 %" relative), comparable to longer pulses of similar fluence; it accumulates during the pulse and persists afterwards, so it is absorptive and "should not be approximated by a Kerr nonlinearity", although a Kerr form can fit the initial rise.
- Simulated: during the pulse the dynamics is coherent (photon–electron term dominates; collisions act after). The population deviation is not the long-pulse non-thermal shoulder but an oscillating population inversion near 𝓔_F, attributed to momentum conservation in the photon–electron term.
- Simulated: stimulated emission grows to nearly the size of absorption at the highest field but never exceeds it; excited-state absorption makes Δε′ superlinear in E₀², stimulated emission makes the net change sublinear (sub-cubic nonlinearity), and Δε′ is linear in the gross absorbed energy only without stimulated emission.
- Simulated: the 3 fs (1/2ω₀) ripple on ε(t) is the signature of coherences; the thermal model reproduces it through the phase lag of Ṗ behind E.
- **Simulated: the thermal model (Eq. 2) gives a permittivity "nearly identical" to ANTHℰM for E₀ = 1.2 V/nm**, "despite the strong deviation from thermal equilibrium", which the authors attribute to the thermal nature of the absorption term and to fast e–e scattering.
- Simulated: the rise time of the extracted ε″ shortens with probe duration (5 vs 20 fs); ε′ is less sensitive.
- Derived and simulated: the generalised Fermi-liquid formula reproduces the exact e–e rate; the finite band bottom lowers A₀ and A₂ below π² and 1 and creates a negative A₁ that grows with T̃_e. The e–e rate grows roughly linearly with E₀²; "single femtosecond" times are projected for 10 TW/cm² (extrapolation).
- Argued: only the Umklapp fraction of e–e scattering contributes to the permittivity rise; the permittivity relaxes on the e–ph time, which is distribution-independent and slow, so e–ph does not explain the 40 fs transmission decay of [20].
- Argued: the cumulative heating limits optical photonic time crystals in TCOs (repeated modulation keeps heating the electrons); THz operation is not affected.
- Stated: the predicted response is "much larger than the experimentally measured transmission and reflection changes" of [20] and its slow decay is "at odds with the observed 40-fs decay time of the transmission". Offered reasons (asserted): electrons driven closer to transparency exposing an instantaneous interband Kerr term (Narimanov, ACS Photonics 12, 402 (2025)); multiphoton absorption (Clerici et al., Nat. Commun. 8, 15829 (2017)); field and temperature gradients, transport and spatiotemporal effects in the 300 nm film; the hydrodynamic model of Scalora et al., Phys. Rev. A 112, 013502 (2025) reproduces [20]'s data.

**Assumptions and limitations.**

- Stated: 4 nm quasi-static sphere (uniform field and temperature); intraband only (the resonance at 0.7 eV is far from interband thresholds); nonlocality neglected ("at most, a modest quantitative" effect); the field-driven coupling between coherences dropped as in [16].
- Absorption in the model is finite-size (Landau-type) through the box dipoles; phonon-, impurity- and Umklapp-assisted absorption, which is what absorbs in a 40–300 nm film, is absent (my reading). The thermal-model agreement is shown for one field and one probe frequency.
- The field displaces the Fermi sphere by roughly k_F at E₀ = 1.2 V/nm (local field ≈2E₀), the regime where the dropped coherence–coherence coupling matters (my estimate, `scripts/field_regime.py` method).
- Up to 1.05×10⁹ J/m³ is absorbed, twice what a 4 eV band can hold above the cold state (≈5×10⁸ J/m³), so the unprinted band top must be well above 4 eV for the highest fields (my calculation).
- The time dependence of T̃_e is not printed, so the "12 fs at 4500 K" statement cannot be checked against a temperature-resolved rate (my reading).

**Use for this project.**

- **Checked against enzkin (`scripts/check_un2026.py`, `results/wp1/check_un2026.json`) [verified in code]:**
  - The exact rate (S19) has the kernel and prefactor of the B1 integral with c_ee = 3/16, so their e–e term is ours.
  - The generalised Fermi-liquid coefficients, implemented independently with the corrections above, reproduce Fig. S1b–d: A₀(𝓔_F) = 9.46, 8.51, 6.06, 4.65, 3.73 against ≈9.4, 8.6, 6.3, 4.9, 4.0; A₂ = 0.95, 0.83, 0.53, 0.40, 0.32 against ≈0.95, 0.85, 0.55, 0.42, 0.35. The formula then matches the B1 rate of a Fermi-Dirac distribution within 3 % near μ (and within 30 % at 0.5𝓔_F at 300 K, where the expansion is poor, as the authors say).
  - B1 gives 1/η_e-e(𝓔_F) = 20 fs at 4500 K and ≈12 fs only at ≈6500 K. Their "≈12 fs at ≈4500 K", and Fig. 5a near 0.5𝓔_F at t = −5 fs, are 1.5–2 times faster than our rate at the stated conditions. Not resolvable from the paper: T̃_e(t), the dielectric constant in the screening and the ½ between (S2) and (S7) are open.
  - Fig. 5c matches our e–ph rate with D = 8.6 eV within 1.2–1.5× at 𝓔_F–3𝓔_F, and not with D = 17.2 eV of [14] (6× too high). The authors did change D.
  - With Table S1 (D = 8.6 eV, Brooks–Herring impurities) and ε∞ = 3.45 the bulk film has Γ = 0.074 fs⁻¹ and ε(230 THz) = 0.14 + 0.17i, against Γ = 0.175 fs⁻¹ and 0.172 + 0.399i of [15]; G = 6.6×10¹⁵ instead of 2.6×10¹⁶ W m⁻³ K⁻¹; a 100 fs pulse depositing 10⁸ J/m³ gives a Re ε recovery (1/e) of 1.9 ps instead of 0.52 ps. Measured film recoveries are 60–400 fs [1, 6, 21, 60]. **The same group's two parameter sets differ by four in the e–ph coupling, so D (or G) must be inferred, not fixed** (plan revision 2.1, decision D15).
- Few-cycle path (plan Tier 2a): the thermal-model result says the permittivity follows the instantaneous net absorbed energy even under a single-cycle pump. In a film the analogue is to drive each kinetic (or thermal) cell by the field-resolved J·E of the Maxwell solver instead of the cycle-averaged absorption term. Validation target: Fig. 4a,b (ε′ −2.6 → −1.0, ε″ 1.85 → 2.25 at 160 THz, E₀ = 1.2 V/nm) in the quasi-static sphere limit.
- e–e rates for thermal models: the generalised Fermi-liquid formula is a cheap, checked e–e rate for M0's damping law (only its Umklapp share damps, per the authors).
- Code: `enzkin.material.ito_un2026()` (Table S1) and Brooks–Herring damping (`enzkin/collisions/eimp.py`, off unless `n_imp > 0`).
- Cited, not in the library: Herbst et al., Adv. Sci. 13, e16818 (2026) (fieldoscopy of ITO nanocrystals under two-cycle pulses, a candidate few-cycle dataset); Scalora et al. 2025 (hydrodynamic model of [20]); Narimanov 2025; Sommer et al., Nature 534, 86 (2016).
- Data availability (verbatim): "The data that support the findings of this study are available from the corresponding author upon reasonable request."

---
## Time-varying ENZ experiments
### [17] Zhou et al., "Broadband frequency translation through time refraction in an epsilon-near-zero material" — Nat. Commun. 11, 2180 (2020)

**Source read:** published PDF, full (main text, Methods, all four figures), cross-checked against the PMC full-text XML. Supplementary Notes 1-6 (setup, simulation model, higher-order effects, sample holder) were NOT available; anything that lives only there is marked "not in main text".
**Role in this project:** WP3 validation-ladder step 3 and the WP4 "probe frequency shift" dataset (a prediction target, not a calibration set).

**What it does.** A degenerate pump-probe experiment with ~120 fs pulses sends a weak probe through a 620 nm ITO stack near its zero-permittivity wavelength (1240 nm) while a strong pump changes the refractive index. The transmitted probe spectrum is recorded against pump-probe delay and pump intensity. The probe red-shifts when it meets a rising index (pump lags the probe) and blue-shifts when it meets a falling index. The headline total translation is 14.9 THz at 1235 nm. The effect persists from 1000 to 1500 nm and is strongest where Re(ε) ≈ 0.

**Model / method.**

Time refraction at a temporal index boundary, for a dispersionless medium with conserved wavevector (Introduction, unnumbered):

$$ n_1 f_1 = n_2 f_2 = (n_1+\Delta n)(f_1+\Delta f) \;\;\Rightarrow\;\; \Delta f = -\frac{\Delta n\, f_1}{n_1+\Delta n} $$

f_1 is the probe frequency before the boundary, n_1 the initial index, Δn = n_2 − n_1, Δf = f_2 − f_1. Δn > 0 gives a red shift.

The paper states that this formula does not describe the film directly: the 120 fs pulses "never reside entirely within the ITO thin film", so the shift depends on the index change rate Δn/Δt met during transit, not on Δn alone (Results).

Simulation: a nonlinear Schrödinger equation solved by the split-step Fourier method in the slowly-varying-envelope approximation. "An iterative algorithm" adjusts "the approximate shape of the time-varying nonlinear phase variations induced by the index change" until the measured spectra are reproduced (Supplementary Note 2; equations not in main text). The modelling is phenomenological: the fitted object is a time-dependent phase, an effective Δn(t), per pump intensity. No two-temperature or Drude calculation is fitted. The two-temperature picture is qualitative only: hot electrons lower the plasma frequency through the momentum-dependent effective mass; the rise is limited by electron-electron thermalization and the pump envelope; recovery is by electron-phonon coupling.

**Numbers.** Experiment card (Methods unless noted). Values "read from" a figure are my readings and approximate.

| Item | Value |
| --- | --- |
| Material | Commercial ITO (PGO GmbH), 310 nm film on 1.1 mm glass |
| Sample | Two such films pressed together in a holder with tightening screws, giving "620 nm" of ITO with a 1.1 mm glass slab on each side (Supp. Note 6). Whether a gap remains between the films: not in main text |
| Permittivity | Re(ε) = 0 at 1240 nm (text). Fig. 4f axis labels: Re(ε) = −2.3, −0.8, 0.3, 1.0, 1.5 at 190, 220, 250, 280, 310 THz. Read from Fig. 1b: Im(ε) ≈ 0.2, 0.35, 0.6 at 1000, 1240, 1500 nm |
| Drude parameters, carrier density, mobility | Not in main text ("well described by the Drude model"). My least-squares Drude fit to the five axis labels: ε_∞ ≈ 3.8, ω_p ≈ 2.95×10¹⁵ rad/s; Im(ε) ≈ 0.35 at 1240 nm then gives γ ≈ 1.4×10¹⁴ s⁻¹ (±15%). An estimate, not the authors' values |
| Source | OPA pumped by an amplified Ti:sapphire laser, "~120 fs"; how measured, and whether intensity FWHM: not stated. Repetition rate: not stated |
| Pump | Degenerate with probe, 1000-1500 nm; p-polarized; 15° incidence; 25 cm lens; spot "~100 µm"; peak intensity 54-268 GW/cm² (Fig. 3), up to about 600 GW/cm² (Fig. 4). Incident or in-film: not stated |
| Probe | Split from the same pulse by a pellicle; p-polarized; 10° incidence; 10 cm lens; spot diameter ~45 µm at 1235 nm; intensity "kept low" |
| Detection | Transmitted probe through a multimode fibre (50 µm core) to an optical spectrum analyser. Observable: probe spectrum against delay t_d (−200 to +200 fs in Fig. 3) and pump intensity; each spectrum normalized individually |
| Temporal resolution | Set by the ~120 fs pulses |

Measured values (no uncertainties are quoted):

- λ_0 = 1235 nm, total translation against intensity (Fig. 3a-c): 4.5 THz at 54 GW/cm², 8.6 THz at 161 GW/cm², 10.8 THz at 268 GW/cm².
- λ_0 = 1235 nm, fixed delay (Fig. 2): red shift 9.1 THz at t_d = −60 fs and blue shift 4.0 THz at t_d = +60 fs, at the largest intensity shown (value not stated).
- λ_0 = 1235 nm, maximum: 14.9 THz = 11.1 THz red + 3.8 THz blue. In Fig. 4c this is the last point, at about 485 GW/cm² (read from figure). The text calls it "over 6%" of the carrier and the conclusion "~6.5%"; 14.9/242.7 THz = 6.1% (my arithmetic).
- λ_0 = 1235 nm, 483 GW/cm²: an additional blue-shifted peak, up to 10.6 THz (~52 nm); total above 20 THz, Δf/f_0 ~ 9%. Attributed to unspecified "higher-order nonlinear optical effects" (Supp. Note 4). The 14.9 THz figure excludes this peak although the intensity is essentially the same (my reading of Fig. 4c).
- Fig. 4a-e, maximum red / blue shift at the highest intensity plotted (read from figure): 1000 nm, −7.6 / +4.0 THz at ≈555 GW/cm²; 1184 nm, −10.8 / +5.3 THz at ≈455; 1235 nm, −11.1 / +3.8 THz at ≈485; 1305 nm, −8.1 / +4.1 THz at ≈515; 1495 nm, −5.5 / +4.8 THz at ≈595. Text: 5.4 THz red shift at 1495 nm.
- Fig. 4f (eight wavelengths, f_0 ≈ 200-300 THz): fractional red shift at ~450 GW/cm² between about 2.2% and 4.4%; total shift at maximum intensity between about 4.7% and 6.4%, peaking near 253 THz (1184 nm).
- Index change: described only as "unity-order". No Δn(t), rise time or decay time is given; "effective values of the index change rate" are in Supp. Note 2.

**Key results.**

- Red shift for t_d < 0, blue shift for t_d > 0, both near t_d ≈ 0, where the front of the probe sees a rising index and the back a falling one (measured).
- Shifts grow with intensity and shrink with |t_d| (measured).
- Red shift exceeds blue shift because the index falls more slowly than it rises (measured; attributed to electron-phonon cooling).
- Shifts peak near the zero-permittivity wavelength (measured, Fig. 4f).
- A fitted phase profile reproduces the delay-resolved spectra (simulated, Fig. 3d-f).

**Assumptions and limitations.**

- "The predictions of our model are only dependent on the envelope-averaged dynamics of the ITO" (authors).
- The film is far thinner than the pulse; the shift mixes Δn, its rate of change and the transit time (authors).
- Stated limits on the shift: linear loss, higher-order nonlinearity, dispersion and "the interplay between the pulse width and the interaction time".
- The formula assumes a dispersionless bulk medium. Interface (Fresnel) phase terms are not discussed; Bohn et al. [18] show they rival the bulk term below about 1 µm (my reading).
- Two films clamped face to face may leave a middle interface or gap (my reading).
- Fig. 4b read-offs sum to about 16 THz at 1184 nm, above the headline 14.9 THz at 1235 nm, and Fig. 4f peaks there; the text does not comment (my reading).

**Use for this project.**

- Dataset: Fig. 3a-c (spectra against delay at 54, 161, 268 GW/cm²) and Fig. 4a-e (maximum shifts against intensity at five wavelengths). Data-availability statement, verbatim: "All data supporting this study are available from the corresponding author upon request." Code: "All relevant computer codes supporting this study are available from the corresponding author upon request."
- Validation targets: the three Fig. 3 totals; the red/blue asymmetry; the peak near ENZ; sublinear growth of the red shift with intensity (Fig. 4).
- The extra blue-shifted peak at 483 GW/cm² is unexplained; it may bear on H2 (my reading).
- H3: the paper says the observable is a transit-limited rate and that the rise depends on "the energy deposition rate" and the pump envelope. A 620 nm absorbing stack also has a depth-dependent pump field (my reading; not discussed).
- Conflict with the plan: WP3 lists films of 20-310 nm. This sample is 620 nm (2 × 310 nm) between two glass slabs, at 15°/10° incidence, p-polarized.
- Obtain the Supplementary Information before use: pulse characterization, sample holder and the effective Δn(t) are there.

---
### [18] Bohn et al., "Spatiotemporal refraction of light in an epsilon-near-zero indium tin oxide layer: frequency shifting effects arising from interfaces" — Optica 8, 1532 (2021)

**Source read:** arXiv:2107.05970v1, full text, Supplement 1, TeX source and figures. The preprint has no subtitle and three authors (Bohn, Luk, Hendry); the published version adds S. Horsley. Published version not compared.
**Role in this project:** WP3 validation step 3, and the adiabatic thin-film (Fresnel) baseline for the H3 test.

**What it does.** For sub-micrometre ITO films the probe frequency shift is the time derivative of the whole transmission phase, which contains interface (Fresnel) terms as well as bulk propagation. For a TM probe at high angle, just below the ENZ frequency, the front-interface term reverses the sign of the shift (blue, not red) while transmission still rises. Measured on 407 nm and 115 nm films; reproduced by a quasi-static transfer-matrix model.

**Model / method.** Shift Δω = −dΦ/dt; bulk phase ΔΦ(t) = Δn(t) k_0 L. Three-layer transmission (1 air, 2 ITO, 3 substrate):

$$ t = \frac{1}{1 + r_{12} r_{23}\,{\rm e}^{2{\rm i} k_2 d}}\; t_{12} t_{23}\; {\rm e}^{{\rm i} k_2 d} \qquad \text{(S1)} $$

$$ \Delta \Phi_{\rm TMM} = \Delta \Phi_{\rm den}+\Delta \Phi_{12}+\Delta \Phi_{23} + \Delta \Phi_{\rm e} \qquad \text{(S2)} $$

$$ t_{12}(\tau) \propto \sqrt{\frac{\varepsilon_{\rm ITO}(\tau)}{\varepsilon_{\rm ITO}(\tau)-1}} \approx {\rm i}\sqrt{\varepsilon_{\rm ITO}(\tau)} \quad (|\varepsilon_{\rm ITO}| \ll 1) \qquad \text{(S4)} $$

d is the thickness, k_2 = n_ITO k_0 cos Θ_2, r_ij and t_ij Fresnel coefficients; (S4) is for TM near grazing incidence. Drude ITO, f_ENZ = (1/2π)√(ω_p²/ε_∞ − γ²). Only ω_p varies: a reduction of up to 10%, shaped as a 107 fs FWHM Gaussian convolved with a 150 fs exponential decay. For comparison with experiment, ω_p(t) is extracted from the measured transmission change through the same model, so the shift has "no fitting parameters". A graded-index variant scales ω_p with the computed depth profile of pump absorption (Sec. S3).

**Numbers.** Experiment card.

- 407 nm ITO (UQG Ltd) on coverslip: ε_∞ = 3.45, ω_p = 2.5×10¹⁵ rad/s, ω_p/γ = 11.5, f_ENZ = 212 THz (211.5 in text; 1417 nm) (Table S1).
- 115 nm ITO on cover glass: ε_∞ = 3.81, ω_p = 3.03×10¹⁵ rad/s, ω_p/γ = 16.9, f_ENZ = 246 THz. Sputtered at room temperature from 90/10 In₂O₃/SnO₂, annealed in forming gas at 425-525 °C (the recipe does not name the sample; the 407 nm film was bought). Carrier density, mobility: not stated.
- Laser: Ti:sapphire 800 nm, 107 fs, 1 kHz, two OPAs; OPA pulse durations not separately reported.
- Pump: TE; incoming I_0 = 400 GW/cm² (100-400 in Sec. S4.2); 480 µm FWHM diameter in air; 5° below the probe angle; 250 THz (407 nm film) or 214 THz (115 nm film).
- Probe: polarized at 45°, TM or TE selected by analyser; 250 µm FWHM; 200 THz (407 nm) or 240 THz (115 nm); angles 0 to ~70°.
- Detection: transmitted spectra against delay (Andor Shamrock 163, DU490A-1.7); observables T/T_0 and central-frequency shift.
- 407 nm: T/T_0 up to 2-3; Δf ≈ −0.7 THz at 5°, ≈ 0 at 50°, ≈ +0.3 THz at 70° (text says 75°) (Figs. 3b, 4b). Model: −0.60 to +0.45 THz. Extracted ω_p reduction up to ~10%.
- 115 nm: −0.30 to +0.60 THz (Fig. 5a).
- Model thickness scan: red shift at 0° linear in thickness, about −0.85 THz at 700 nm; blue shift at 70° saturates near +0.6 THz above ~200 nm (Fig. 5b).
- No uncertainties quoted.

**Key results.**

- Sign reversal of the shift with angle at positive ΔT; only for TM, only below the ENZ frequency (measured).
- Interface blue shift can exceed the bulk red shift for films up to ~700 nm (simulated).
- A depth-graded ω_p gives a smaller blue shift than a uniform layer and matches experiment better (simulated, Fig. S5).

**Assumptions and limitations.**

- Quasi-static: ε(t) in steady-state Fresnel formulas; probe dispersion and bandwidth ignored (my reading).
- Uniform excitation in depth except in Sec. S3, whose profile the authors call "a rough assumption".
- Intensity on the sample falls at high angle as the beam footprint stretches (authors).

**Use for this project.**

- H3: direct evidence that angle, thickness, polarization and depth profile change the measured shift at fixed material dynamics; internal reflections matter for d ≲ 200 nm.
- Use (S1) with ε(t) as the adiabatic baseline; compare the Maxwell-kinetic result with Figs. 4b, 5a, S7.
- Data: arXiv v1 has no data-availability statement (the TeX source holds a commented-out placeholder naming the University of Exeter repository, no URL). Check the published version.
- The authors attribute their smaller shifts relative to [17] to higher Drude damping and thinner films.

---
### [19] Tirole et al., "Double-slit time diffraction at optical frequencies" — Nat. Phys. 19, 999 (2023)

**Source read:** arXiv:2206.04362v2 (PDF dated 16 June 2022), main text and Supplementary Materials in full, all figures. No TeX source. The published Nature Physics version was not compared; numbers are from the arXiv version.
**Role in this project:** WP4 dataset "time-diffraction fringes", the main constraint on rise time and the anomaly behind H1 and H3.

**What it does.** Two pump pulses, separated by a tunable delay, each switch the reflectivity of an ITO/gold mirror from 0.08 to 0.6. A longer probe pulse spans both switching events ("time slits"), and its reflected spectrum shows fringes whose period is the inverse of the slit separation. The fringes stay visible much further from the carrier than a pump-limited rise would allow. A Fourier-transform model with the slit rise as the fitted shape parameter gives a 10-90% rise of 7 fs, with 1 to 10 fs fitting equally well, under a 225 fs pump.

**Model / method.** Equations are unnumbered, in the Supplementary Materials. Single-slit aperture function, α and β positive constants:

$$ f_{ss}(t) = \frac{1}{(1+e^{-\alpha t})\times(1+e^{\beta t})} $$

Double-slit amplitude reflection coefficient, taken real:

$$ r(t) = A \times f_{ss}\!\left(t-\frac{S}{2}\right) + B \times f_{ss}\!\left(t+\frac{S}{2}\right) + C $$

S is the slit separation, A and B the slit amplitudes, C the unmodulated reflection coefficient; R(t) = |r(t)|². The asymptotic case (α → +∞, β → 0) uses Heaviside functions:

$$ r(t) = A \times H\!\left(t-\frac{S}{2}\right) + B \times H\!\left(t+\frac{S}{2}\right) + C $$

Probe spectrum:

$$ I(f) = \mathcal{F}^{2}\big( r(t)\cdot E_{probe}(t-t_0) \big) $$

E_probe is a Gaussian envelope of 794 fs FWHM on a 230.2 THz carrier; t_0 is an offset between the probe peak and the centre of the slits; the transform is squared to give intensity.

How the rise time was inferred:

1. Pump-probe with a short (225 fs) probe gives R against delay (Fig. S3A), fitted with "r(t)² ∗ E_probe(t)²" (a convolution). This fixes A, B, C, S and β. The rise cannot be resolved: "the rise time is limited by the probe pulse".
2. The 794 fs probe spectra are fitted with I(f), separation and decay held at the step-1 values; "the rise time is taken as a fitting parameter" (α). t_0 is a further model parameter.
3. Result: α = 1/2 fs⁻¹, "an intensity reflectivity rise time (10-90%) of 7 fs". The data are close to the Heaviside limit, where fringe intensity falls as 1/(f − f_0)², "therefore the fitting is equally good for a reflectivity rise time (10-90%) in the range 1-10 fs".

Stated approximations: r(t) real (no phase evolution); non-dispersive ("the same temporal response at all frequencies"); the saturation plateau is not modelled.

A second, adiabatic and dispersive model computes "a negative shift in the plasma frequency from the convolution of the gaussian pump intensity with an exponential rise and decay response function". It reproduces fringes, period and red/blue asymmetry but not the spectral extent, because its rise "is fundamentally limited by the pump duration". No parameters are given.

**Numbers.** Experiment card.

| Item | Value |
| --- | --- |
| Sample | 40 nm ITO on a glass coverslip, covered by 100 nm gold. Fig. 1G labels the stack SiO2/ITO/Au and draws the beams on the glass side (my reading of the sketch) |
| ENZ frequency | 227 THz (1320 nm) |
| Drude parameters, carrier density, mobility, growth | Not stated; deferred to ref. (16), Tirole et al., "Saturable time-varying mirror based on an ENZ material" |
| Source | Light Conversion Pharos with Orpheus OPA, 230.2 THz (1300 nm). Repetition rate: not stated |
| Pumps | Two pulses, 225 ± 5 fs FWHM of the electric-field envelope (autocorrelation, Gaussian fit); intensity FWHM is √2 shorter (159 fs, my arithmetic). p-polarized; 54° and 66°. Intensity "around 124 GW/cm²", saturating; Fig. S3C plots up to about 90 GW/cm². Incident or in-film, per pump or total: not stated. Spot size: not stated |
| Probe, characterization | Same 225 fs pulse, 60°, p-polarized |
| Probe, diffraction | Spectrally filtered in a 4-f line: 794 ± 37 fs field FWHM (sum-frequency cross-correlation with the pump in 400 µm GaP, Gaussian fit); 1 THz bandwidth; 230.2 THz; 60°; p-polarized |
| Detection | Reflected probe spectrum, Princeton Instruments NIRvana camera |
| Temporal resolution | Pump-probe: limited by the 225 fs probe. Spectral method: set by how far from the carrier the fringes are tracked, about 10 THz (my reading) |

Measured and fitted values:

- Reflectivity 0.08 → 0.6, "modulation contrast of 76%"; saturates at 0.6 (Fig. S3C).
- β = 1/400 fs⁻¹: "decay time of the intensity reflectivity (defined as the time to decrease from 90% to 10%) of 625 fs (or 1/e of 330 fs)".
- α = 1/2 fs⁻¹: 10-90% rise 7 fs; acceptable range 1-10 fs. No formal uncertainty or goodness-of-fit.
- My numerical check: the quoted 7, 625 and 330 fs are those of |f_ss(t)|² alone (7.1, 623, 331 fs). For the full R = |A f_ss + C|² rising from 0.08 to 0.6, the same α and β give about 8 fs, 0.9 ps and 0.45 ps. The 10-90% rise of |f_ss|² is about 3.6/α.
- Second slit amplitude "0.93% of the first one" (as printed; Fig. 1H shows near-equal peaks, so a ratio of 0.93 is evidently meant).
- Slit separations: 800 and 500 fs (Fig. 2A, B), 2.3 ps (Figs. 1H, S3A), scan of about ±1.2 ps (Fig. 2E, F). Below 300 fs the slits merge.
- Fringe period ∝ 1/S (Fig. 2E); error bars 50 fs in separation, 0.15 THz in period.
- Up to six oscillations; out to 10 THz (~60 nm) on the red side, 4 THz on the blue side.
- Fig. S4D (S = 800 fs, log scale): model curves for response times 3.6, 7, 17, 32, 47, 61 fs. The 3.6 and 7 fs curves follow the data; longer ones fall below it towards 221 THz, by roughly an order of magnitude for 61 fs (read from figure).
- Sensitivity (my estimate, not in the paper): a logistic edge multiplies the fringe intensity at detuning Δf by [x/sinh x]², x = 2π²Δf/α. At 9 THz this is 0.99, 0.96, 0.80 and 0.09 for rises of 3.6, 7, 17 and 61 fs, consistent with Fig. S4D and the 1-10 fs degeneracy.
- Optical cycle quoted as 4.4 fs.

**Key results.**

- Spectral fringes from a temporal double slit, period ∝ 1/S (measured).
- Fringe intensity scaled by (f − f_0)² is nearly constant, the signature of a step-like rise (measured, Fig. 2C, D).
- Rise of 1-10 fs, over ten times shorter than the pump (inferred from the fit).
- Red/blue asymmetry attributed to the time-varying phase of r, i.e. time refraction (asserted; reproduced qualitatively by the adiabatic model).

**Assumptions and limitations.**

- The fitted quantity is the rise of the reflection coefficient of a resonant ITO/gold stack at 60°, not of the ITO permittivity. The authors call the measured reflectivity "a convolution of the driving pump pulses, the probe and the material response".
- The mirror is saturated, so R is a strongly nonlinear function of the permittivity change; a saturating map can steepen the edge of R(t) relative to ε(t) (my reading).
- 1-10 fs is in effect an upper bound near 10 fs; the data cannot tell 1 fs from 7 fs.
- Fixed functional form, one α for both slits, Gaussian probe.
- "A phenomenological model of the phase evolution is known" (as printed; the context implies it is not known).

**Use for this project.**

- Fit target: fringe envelope against detuning at fixed S (Figs. 2A-D, S4D) and the interferogram (Fig. 2F), not a single rise time. Forward model: complex, dispersive r(ω, t) of glass/ITO(40 nm)/Au(100 nm), p-polarized, pumps at 54° and 66°, probe at 60°.
- Data availability, verbatim: "Data and codes are available upon request."
- H3: the cleanest test. Field enhancement at the Berreman resonance, saturation, phase evolution and probe dispersion are all outside the authors' fit, and each can change the apparent rise.
- Pitfall: 225 fs and 794 fs are field-envelope FWHM; the plan's "225 fs pump" is 159 fs in intensity.
- Discrepancy with the task description: the plan's reference [15] (Un et al.) quotes pump 60° and probe 65° for the earlier "saturable time-varying mirror" experiment, its reference [8]. Here the probe is at 60° and the pumps at 54° and 66°. The 40 nm ITO on glass with 100 nm gold is confirmed.
- Missing for modelling: ITO Drude parameters, spot sizes, repetition rate, whether 124 GW/cm² is incident. See ref. (16) or request.

---
### [20] Lustig et al., "Time-refraction optics with single cycle modulation" — Nanophotonics 12, 2221 (2023)

**Source read:** PMC full-text XML (PMC11501756), main text and figure captions in full, plus the abstract file. No PDF, so figure panels were not seen; Supplementary Information not available.
**Role in this project:** the "fast rise and fast recovery under near-single-cycle pumping" anomaly; Tier-2 (sub-8 fs) target and an H3 test case.

**What it does.** An 800 nm modulator pulse, compressed from 30 fs to about 6 fs, changes the index of a 310 nm ITO film while a 40-44 fs probe near the ENZ wavelength passes through. The transmitted probe spectrum red-shifts and then blue-shifts within tens of femtoseconds. Shorter modulators give a faster and larger red shift. The blue shift implies an index recovery of 10-30 fs, which the authors leave unexplained.

**Model / method.** Phenomenological:

$$ \frac{{\rm d}}{{\rm d}t}\Delta\varepsilon(t) \approx I_{\rm modulator}(t) - \frac{1}{T_{\rm decay}}\Delta\varepsilon(t) \qquad (1) $$

$$ n(t) = {\rm Re}\left\{\sqrt{\varepsilon_{I_{\rm modulator}=0} + \Delta\varepsilon(t)}\right\} \qquad (2) $$

I_modulator(t) is the modulator intensity "in units of [permittivity/time]", T_decay the relaxation time, ε at zero modulator intensity the unmodulated permittivity. The rise is taken as instantaneous because "absorption is an induced transition process". FDTD simulations pass a probe through a slab whose index follows Eq. (2) (details in SI). The red-shift rise is fitted by an error function convolved with the probe envelope (Fig. 4B). T_decay is not formally fitted: runs with 20 fs and 100 fs are compared with the data.

**Numbers.** Experiment card.

- Sample: ITO, 0.31 µm; permittivity by ellipsometry (Fig. 1D); ENZ point 1225 nm. Growth, substrate, Drude parameters, carrier density: not in main text.
- Modulator: 800 nm; 30 fs compressed in an Ar-filled hollow-core fibre; FROG intensity FWHM 6 fs (6.8 fs in Fig. 2; "~5-6 fs" in abstract); 6 to 30 fs used; 700 µm diameter region; "~2 TW/cm²" on the sample, incident or in-film not stated. Angle, polarization, repetition rate: not stated.
- Probe: down-converted, 1200 nm (Fig. 1 caption), "near the ENZ point at 1225 nm" (text); 40 fs (Fig. 1) or 44 fs (Figs. 3, 5), FROG; five orders of magnitude weaker; much smaller spot.
- Detection: output facet imaged onto photodiodes (transmitted intensity; Fresnel reflection from the front interface) and a spectrometer (transmitted spectrum), against delay. Zero delay set at 20% of the maximum change.
- Temporal resolution: limited by the probe; the apparent response is "a convolution of the index response and the temporal envelope of the probe pulse". Probe transit time "on the order of 1 fs".
- Index change: real part +~0.5, imaginary part −~0.15 (plane-wave transfer-matrix estimate, SI); ~0.12 per probe cycle for the shortest modulator; 3-4 cycles to reach 0.5; "raised within 5-10 fs".
- Recovery: transmitted intensity and blue shift, "20-30 fs" (Sec. 6) or "10-20 fs" (Secs. 7, 8); Fresnel reflection, "several hundreds of femtoseconds".
- FDTD: 20 fs Gaussian modulator; T_decay = 20 fs resembles the data, 100 fs gives almost no blue shift; probes of 44 and 20 fs.
- AZO: index rise up to 0.15, slow recovery, almost no blue shift (SI).
- No uncertainties quoted.

**Key results.**

- Red-shift rise time follows the modulator duration down to ~6 fs (measured).
- Red and blue shifts have similar durations and overlap in delay because the probe is longer than the index transient (measured; simulated).
- Shorter modulators give a larger mean red shift (measured). The authors' explanation is electromagnetic: the film is thinner than the rise, so the shift depends on how much the index changes during transit.
- Transmission and reflection give different recovery times, attributed to slower relaxation at the surface (asserted).

**Assumptions and limitations.**

- Spatially uniform, non-dispersive Δε(t) driven linearly by absorbed energy (my reading).
- Quoted recovery time varies between 10 and 30 fs within the paper.
- Carrier-envelope phase not measured. Pump at 800 nm and ~2 TW/cm², unlike the ENZ-pumped datasets.

**Use for this project.**

- H3: the authors attribute the pulse-width dependence of the shift to finite thickness and the red/blue overlap to probe duration. The reflection/transmission mismatch is a surface-versus-bulk observable for the Maxwell-kinetic solver.
- Conflict: recovery of 10-30 fs here against 60-170 fs, lengthening with intensity, in Gurung et al. [21].
- Data availability: the article has no data-availability statement. Fig. 3 panel values need the PDF.
- Outside the cycle-averaged Boltzmann scope (~6 fs pump); stretch target only.

---
### [21] Gurung et al., "Control of ultrafast hot electron dynamics in epsilon-near-zero conductive oxide thin films" — Sci. Adv. 11, eadu8850 (2025)

**Source read:** published PDF, full (main text, Materials and Methods, Figs. 1-4), cross-checked against the PMC XML. The Supplementary Materials were NOT available; they hold the hydrodynamic equations (section S4), the fit function (eq. S23), the material parameters (table S1) and the pulse durations (fig. S3). Obtain them before use.
**Role in this project:** WP3 validation-ladder step 4 and the WP4 dataset "relaxation time against intensity and angle"; a hydrodynamic baseline for recovery dynamics.

**What it does.** Degenerate pump-probe transmission measurements on a 306 nm ITO film and a 205 nm AZO film give a "hot electron relaxation time" that lengthens with pump intensity and with the field enhancement of the radiative ENZ mode (varied through pump angle and wavelength). A hydrodynamic model derived from the Boltzmann equation attributes this to quenching of electron-phonon damping and thermalization for fast electrons, and gives a closed-form first-order relaxation rate.

**Model / method.** The only model equation in the main text is the perturbative relaxation rate:

$$ \gamma_{\rm r} \approx \gamma_{\rm th} - \frac{\gamma_{\rm th}\,(2\gamma-\gamma_{\rm th})\,e^{2} f(\theta,\omega)}{12\,\gamma\, k_{\rm B} T_{0}\, m\, \varepsilon_{0} c\,(\omega^{2}+\gamma^{2})}\; I_{0} \qquad (1) $$

γ_r = 1/τ_r is the hot-electron relaxation rate; γ the linear dephasing rate (Drude damping); 1/γ_th = τ_th the "linear thermalization time"; T_0 the lattice temperature; m the electron effective mass; ω the carrier angular frequency; I_0 the intensity of the impinging radiation; f(θ,ω) the field-intensity enhancement. Field enhancement and slowdown measure (Materials and Methods, unnumbered):

$$ f(\theta,\omega)=\left|\frac{E}{E_{\rm i}}\right|^{2}=\frac{A(\theta,\lambda)\cdot n_{\rm i}\cdot \cos(\theta)}{\varepsilon''}\,\frac{\lambda}{2\pi t}, \qquad \zeta(\theta,\lambda)=f(\theta,\lambda)\times\lambda^{2} $$

E is the field averaged over the film, E_i the incident field, n_i the index of the incident medium, ε'' the imaginary permittivity, t presumably the film thickness (not defined), A the absorptance from a linear transfer-matrix calculation with ellipsometric indices.

The full model, described in words only: Boltzmann equation with full electron-electron and electron-phonon collision integrals; Landau weak-coupling limit (large-angle collisions neglected), giving a Fokker-Planck-Landau equation; Rosenbluth potentials evaluated analytically; method of moments, giving equations for electron density, mean velocity and temperature with intensity-dependent dephasing, heating and relaxation rates. Electron-electron collisions thermalize the electrons instantaneously. The calculation is local, with a film-averaged field and uniform temperature, integrated by fourth-order Runge-Kutta. The equations are in section S4 and refs. (27) Silvestri et al., Nanophotonics 13, 2003 (2024) and (29) Marini et al., Faraday Discuss. 214, 235 (2019).

Parameters: τ_th is not predicted. It is "estimated from the intensity-dependent pump-probe experiment ... by extrapolating the relaxation time to very low intensity" and "used as a fixed parameter": 21 fs for AZO, 42 fs for ITO. Plasma frequency, γ and m are in table S1.

Definition of relaxation time: "the time required for energetic electrons to equilibrate with the lattice". Operationally, τ_r comes from fitting the normalized probe transmittance ΔT/ΔT_max "using the same dependence as the normalized electron temperature in the TTM" (eq. S23, not seen), on the stated basis that ΔT follows ΔT_e. It is the decay constant of a transmission transient, equated with the electron-temperature decay.

**Numbers.** Experiment card. Figure read-offs are mine and approximate.

| Item | Value |
| --- | --- |
| ITO | 306 nm; RF magnetron sputtering, 90:10 In₂O₃ target, 50 W, base pressure 1×10⁻⁵ mbar, Ar 60 sccm, oxygen "adjusted to 5×10⁻³ mbar", substrate at 500 °C, 20 nm/min for 900 s. Substrate material: not stated |
| AZO | 205 nm; atomic layer deposition (diethylzinc, water, trimethylaluminum) on 500 µm fused silica |
| ENZ wavelength | ITO 1190 nm; AZO 1480 nm |
| Carrier density, mobility, Drude parameters | Not in main text (section S1, table S1) |
| Source | Coherent Ti:sapphire with OPA, 1100-1600 nm, 1 kHz; "60 fs" (Fig. 1D label). Duration varies with wavelength (fig. S3, not seen) |
| Pump | 90% arm; TM-polarized; 200 mm lens; radius 1.33 times the probe's (absolute sizes in table S3); peak intensity about 20-300 GW/cm², taken as incident ("impinging"); angle 0-50° |
| Probe | 10% arm, same wavelength; TE-polarized; 150 mm lens; drawn 10° from the pump in Fig. 1D (my reading of the schematic; not in text) |
| Detection | Probe transmission against delay; step 1 µm (3.3 fs) |
| Temporal resolution | Set by the ~60 fs pulses; rise of ΔT about 50 fs (read from Fig. 2E) |

Intensity dependence (normal incidence, pump and probe at the ENZ wavelength; Fig. 2F and text):

- ITO: 62.6 ± 8 fs at ~20 GW/cm² to 169.5 ± 3 fs at ~200 GW/cm².
- AZO: 31.4 ± 7 fs at ~20 GW/cm² to 89.6 ± 3 fs at ~200 GW/cm² (in Fig. 2F these points sit near 34 and 209 GW/cm²).
- Fig. 2F points (intensity in GW/cm², τ_r in fs): ITO (20, 66), (38, 97), (52, 115), (67, 134), (95, 136), (111, 149), (133, 151), (190, 169); AZO (34, 34), (48, 45), (71, 53), (103, 64), (124, 68), (152, 76), (209, 90), (275, 96). Error bars about ±4 to ±11 fs.
- The text says "nearly three times". Its numbers give ratios of 2.7 (ITO) and 2.9 (AZO), increases of 171% and 185%; relative to τ_th, 4.0 and 4.3 times (my arithmetic). The abstract's "more than 200%" fits only the second reading, or "to more than 200% of".
- The dependence is sublinear (authors: higher-order terms, falling dephasing rate, saturating thermalization time).

ENZ-mode dependence (100 GW/cm²; Figs. 3, 4):

- τ_r normalized to its maximum: 81 fs (AZO) and 155 fs (ITO), at 30° and the ENZ wavelength.
- Angle scan 0-50° in 10° steps, τ_r/τ_max (experiment): AZO 0.70, 0.84, 0.93, 1.00, 0.96, 0.87; ITO 0.86, 0.92, 1.01, 1.06, 1.00, 0.82.
- Text: ENZ-mode excitation changes τ_r by "more than 24 fs for AZO and 29 fs for ITO"; normal-incidence values are "approximately 20 to 40% lower" than at the peak. That is about +42% for AZO and +23% for ITO (my arithmetic); the abstract's "more than 40%" is the AZO figure.
- Wavelength scan at normal incidence, τ_r/τ_max (experiment): AZO 0.52-0.70 over 1190-1580 nm, rising; ITO 0.80-0.56 over 1170-1550 nm, falling.
- For ITO the analytic ζ varies far more with angle (about 0.62 to 1.0) and lies well below the data against wavelength (0.6 to 0.35); the authors call ζ "a first-order approximation".

**Key results.**

- τ_r grows with pump intensity in both oxides (measured; reproduced by the hydrodynamic model with τ_th fixed from the same data).
- τ_r peaks near 30-35° pump incidence, following the field enhancement (measured; error bars of about ±0.1 to ±0.17 in τ_r/τ_max are not much smaller than the trend, read from Figs. 3B, 4B).
- AZO relaxes about twice as fast as ITO at equal intensity (measured).

**Assumptions and limitations.**

- Authors: local limit, no nonlocality or spill-out (justified by thickness > 100 nm); large-angle collisions neglected; instantaneous electron thermalization; thickness enters only through f; uniform temperature and film-averaged field, "accurate only for films much thinner than the wavelength".
- f comes from linear optics and is not updated as the permittivity changes during the pump (my reading).
- ΔT is mapped one-to-one onto ΔT_e, with no Maxwell propagation of the probe (my reading). H3 questions exactly this.
- The shortest AZO times (31-45 fs) are below the 60 fs pulse duration and depend on the deconvolution in eq. S23 (my reading).
- The plotted AZO model curve extrapolates to roughly 30 fs at zero intensity, not obviously 21 fs (my reading of Fig. 2F).

**Use for this project.**

- Fit targets: Fig. 2F and Figs. 3B, 4B. Better: digitise the transients (Figs. 1B, 2A, 2C, 2E) and fit ΔT(t) directly, since τ_r depends on an unseen fit function.
- Data availability, verbatim: "All data needed to evaluate the conclusions in the paper are present in the paper and/or the Supplementary Materials." No repository.
- Eq. (1) is a baseline to beat. With assumed m = 0.35 m_e and γ = 1-2×10¹⁴ s⁻¹ (background, not from this paper), its first-order term lowers γ_r by about 0.45 f at 100 GW/cm² and 1190 nm (my estimate), so the expansion parameter is not small over the measured range.
- Conflict in the library: Un et al. [15] predict that the electron-temperature decay becomes faster with intensity; this paper notes that and measures the opposite. Lustig et al. [20] report 10-30 fs recovery at ~2 TW/cm².
- H3: at oblique TM incidence the absorbed-power profile and the TE probe's transmission both change with angle; test whether the angle trend survives a full Maxwell treatment (my suggestion).

---
### [22] Harwood et al., "Competing Auger-like and intraband excitations drive anomalous modulations in indium tin oxide" — arXiv:2605.21432v2 (2026); published as Optica 13, 1943 (2026)

**Source read:** arXiv v2 full text + TeX source, main text only. The Supplementary Information (SI) is not in the library; it holds the retrieval algorithm, the formula for $n_0(T)$, the Kane parameters, the depth-resolved transfer-matrix check, the transmission retrievals and further fits. Published Optica version not read. The plan cites the v1 title, "Intraband and interband competition drives ultrafast modulations of indium tin oxide".
**Role in this project:** Defines model M2 (Auger-type interband source term, WP1/WP4), supplies the WP4 dataset of complex reflection coefficients at high fluence, and is the precedent for the complex-gate retrieval of WP5.

**What it does.** A 310 nm ITO film is pumped and probed with 136 fs pulses at 1300 nm, at pump fluences of 4.6–153 mJ/cm². The delay-resolved probe spectrum is treated as a FROG trace whose gate is the film's time-dependent complex Fresnel coefficient; this yields amplitude and phase of $r(t)$, which a transfer-matrix inversion converts into Drude $\omega_p(t)$ and $\gamma(t)$. Above about 55 mJ/cm² the $\omega_p$ dip becomes asymmetric and $\omega_p$ overshoots equilibrium, while $\gamma$ rises several-fold. The authors reproduce this with a thermal conduction-electron model in which the conduction-band density relaxes, at a fixed rate, toward its value for equilibrium with the valence band at the current electron temperature.

**Model / method.**

*Retrieval ("GRUMPY FROG").* The probe field $E(t)$ is gated by an "effective time-varying, complex and adiabatic" reflection coefficient $r(t)=|r|e^{i\varphi(t)}$, valid if the film is non-dispersive across the probe bandwidth:

$$ I(\omega,\tau) = \left|\int_{-\infty}^{\infty} E(t)\cdot r(t-\tau)\,e^{-i\omega t}\,\mathrm{d}t\right|^2 \qquad \text{(Eq. 2)} $$

$I(\omega,\tau)$ is the measured probe spectrogram and $\tau$ the delay. An iterative FROG algorithm ("well-known algorithms"; not specified) returns both $E(t)$ and $r(t)$. The normalised $r(t)$ is scaled by the static Fresnel coefficient from ellipsometry; uncertainty bands come from bootstrapping. A transfer matrix for a *uniformly modulated* Drude film then gives effective $\omega_p(t)$, $\gamma(t)$. Transmission coefficients were also retrieved (SI); reflection is preferred because the 1 mm substrate adds cross-phase modulation.

*Carrier model.* The state is the conduction-electron density $n_c(t)$ and conduction-electron energy density $U_c(t)$. There is no lattice-temperature equation and no absorption term driven by the field.

$$ \dot{n}_c = \eta\left[n_0 - n_c\right] \qquad \text{(Eq. 3)} $$

$n_0$ is "the conduction electron density assuming perfect equilibrium with the valence band at temperature $T$"; $n_c$ is the actual density; $\eta$ is called the Auger excitation rate. The same term gives excitation ($n_0>n_c$) and recombination ($n_0<n_c$). Solution: $n_c(t)=\eta\int_{-\infty}^{t}\mathrm{d}t'\,e^{-\eta(t-t')}n_0(t')$.

$$ U_0(t) = \frac{1}{2}U_{\max}\left[e^{-\beta t} + U_{\infty}/U_{\max}\right]\left[\tanh(\alpha t) + 1\right] \qquad \text{(Eq. 4)} $$

$U_0$ is the conduction-electron energy density without interband exchange; $\alpha$ is the inverse rise time set by the pump, $\beta$ the electron–phonon decay rate, $U_\infty$ a residual energy fitted to long delays (origin "not understood").

$$ \dot{U}_c(t) = \dot{U}_0(t) - E_g\,\dot{n}_c(t) \qquad \text{(Eq. 5)} $$

Each promoted electron costs $E_g$, assuming excitation from the valence-band edge.

*Closure and optics.* Assuming a Fermi–Dirac distribution, $T_c(t)$ and $\mu_c(t)$ "follow from $n_c$ and $U_c$ numerically". The plasma frequency is $\omega_p=\sqrt{n_c e^2/m^*\varepsilon_0}$ (introduction, unnumbered), with $m^*$ from a "modified Kane model": energy-dependent below the Fermi level, approximately constant above it, the high-energy value rescaled to fit. The averaging formula is not given. The damping is the computed valence-hole density "scaled by a constant"; no equation is printed.

Approximations: conduction electrons always thermal; one rate $\eta$ with no explicit temperature or density dependence (it all sits in $n_0(T)$); prescribed energy input; uniform film; adiabatic, non-dispersive gate.

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| Film / substrate | ITO 310 nm on glass (SiO₂, 1 mm) | Results; Fig. 1a |
| Carrier density, mobility, $\varepsilon_\infty$ | not stated (main text) | — |
| ENZ wavelength | "near" 1300 nm; exact value not stated | Results |
| Static $R$, $T$ at 1300 nm | ≈46 %, ≈20 % (read from figure) | Fig. 1b |
| Static $\omega_p$, $\gamma$ | ≈3.03×10¹⁵ rad/s; ≈0.11–0.14×10¹⁵ rad/s (read from figure) | Fig. 4d–i |
| Pump | 1300 nm (230 THz), 136 fs, p-polarised, 8°, 4.6–153 mJ/cm², 20 Hz | Results |
| Fluence definition, peak intensity | not stated (presumably incident). Gaussian 136 fs FWHM would give 55 mJ/cm² ≈ 0.38 TW/cm², 153 mJ/cm² ≈ 1.06 TW/cm² (my conversion) | — |
| Probe | 1300 nm, 136 fs, s-polarised, 4°, ≈0.1 mJ/cm², analyser at output | Results |
| Observable | spectrally resolved reflected and transmitted probe vs delay; temporal resolution not stated (retrieval deconvolves the probe) | Results |
| Retrieval error | 0.3 % rms between measured and retrieved spectrogram at 133 mJ/cm² | Results |
| Onset of anomalous dynamics | >55 mJ/cm² (fast partial recovery, oscillation; phase sign reversal near −100 fs) | Figs. 1d, 3 |
| $\omega_p$ change | saturates near −20 % at low fluence; above threshold exceeds equilibrium "by over 3 %" (elsewhere "over 4 %"; Fig. 4f shows ≈+2 % at 133 mJ/cm², my reading) | Results; Fig. 4f |
| $\omega_p(t)$ at 133 mJ/cm² (my reading) | min ≈2.36 at ≈−50 fs, max ≈3.08 at ≈+280 fs, second min ≈2.84 at ≈850 fs (×10¹⁵ rad/s); intraband-only model min ≈2.07 (text: −33 %) | Fig. 4f |
| $\gamma$ peaks | "exceeding 600 %" of static (caption: "700 % above"); ≈0.39 and ≈0.66×10¹⁵ rad/s at 57 and 133 mJ/cm² (my reading) | Results; Fig. 4h,i |
| $r(t)$ at 133 mJ/cm² (my reading) | $\lvert r\rvert$: 0.69 → 0.15 near −30 fs → ≈0.5 by 300 fs; $\varphi$: 1.65 → 1.08 rad near 0 fs → 1.84 rad near 1.1 ps | Fig. 2c |
| Delay of interband contribution | peak at 120 fs (measured) | Results |
| $\eta$ | 1/120 fs⁻¹ ("assign") | below Eq. 3 |
| $\alpha$, $\beta$ | 1/50 fs⁻¹; 1/300 fs⁻¹ (fixed, low-fluence value) | below Eq. 4 |
| $E_g$ | 3.9 eV (direct gap) | below Eq. 5 |
| High-energy $m^*$ | 0.75 × Kane extrapolation (fitted at highest fluence) | Results |
| Valence density of states, $U_{\max}$, $U_\infty$, hole-scattering constant, Kane parameters | not stated | — |
| Parameter uncertainties | none given | — |
| Electron temperature | intraband-only: >5 eV (also ">6 eV (70,000 K)"); with Auger: "below 1.5 eV" in equilibrium; transient spike ≈3.4 eV vs ≈5 eV at 133 mJ/cm² (my reading) | Results; Fig. 4b,c |

**Key results.**
- Measured: above ≈55 mJ/cm², $R$ and $T$ show a fast partial recovery and oscillation, asymmetric between $R$ and $T$, which signals a change of $\gamma$ as well as $\omega_p$.
- Retrieved: $\omega_p$ changes sign relative to equilibrium and $\gamma$ rises several-fold, with the same fluence threshold.
- Argued: a system in thermal equilibrium passes each temperature twice, so $\omega_p(t)$ would show two minima of equal depth; the observed unequal minima are "model-independent evidence for an out-of-equilibrium process".
- Argued: two-photon absorption "operates only during the pump pulse and is ruled out by this delay, which exceeds the pump duration".
- Argued: the threshold onset of $\gamma$ is "inconsistent with temperature-dependent scattering models", which predict smooth growth from the lowest fluence.
- Fitted: the model follows $\omega_p(t)$ at 19, 57 and 133 mJ/cm² with band parameters fixed at the highest fluence; $\gamma(t)$ follows the hole density, less closely.
- Model output: interband excitation is a heat sink that limits electron temperature; offered as the reason for the high damage threshold (asserted).

**Assumptions and limitations.**
- Authors: the uniform-film inversion makes $\gamma$ dip below its static value at long delay, an artefact removed by a depth-dependent model (SI); the $\gamma$ model omits impurity, phonon and electron–electron scattering.
- Not a two-temperature model in the usual sense: no lattice temperature and no absorbed power; the relation between $U_{\max}$ and fluence is not given (my reading).
- The energy-deposition pulse implied by $\alpha$ is sech² with FWHM $1.76/\alpha\approx 88$ fs, shorter than the 136 fs pump (my reading).
- $\eta$ is set equal to the observed 120 fs delay. With one duration and one wavelength, the data cannot separate a fluence threshold from an intensity threshold (my reading).
- A two-photon carrier population is cumulative and persists after the pulse, so the timing argument against it is weaker than stated. The stronger argument is energetic: two 1300 nm photons carry 1.91 eV, below both 3.9 eV and the 2.67 eV gap of [23] (my reading).
- The gate assumption (adiabatic, non-dispersive) is what hypothesis H3 questions; its check is only in the SI.
- Internal inconsistencies: 131 vs 133 mJ/cm²; 3 % vs 4 %; 600 % vs 700 %; 5 vs 6 eV; the text cites Fig. 4b/4c for panels 4a/4b.

**Use for this project.**
- M2 source term: $\dot n = \eta[n_0(T)-n]$ with energy sink $E_g$ per carrier (Eqs. 3, 5). The plan's "Auger coefficient" is therefore a relaxation rate ($\eta^{-1}=120$ fs), not a $Cn^2p$ coefficient. A distribution-resolved solver needs a definition of $T$ (or an impact-ionisation rate over $f(E)$) and a valence density of states; neither is given numerically, so the SI or the authors are needed.
- Dynamic damping is required at high fluence; first model: $\gamma$ proportional to hole density.
- Fit $r(t)$ (Fig. 3a,b) directly; the published $\omega_p(t)$, $\gamma(t)$ embed a uniform-film Drude inversion with unstated $\varepsilon_\infty$.
- Data availability (verbatim): "Data underlying the results presented in this paper may be obtained from the authors upon reasonable request."
- Open question: dependence on pulse duration at fixed fluence, the plan's discriminator against M3, is not measured here.
- Conflicts: gap 3.9 eV here against 2.67 eV in [23]; the threshold argument on $\gamma$ contradicts the smooth $\Gamma(T_e)$ law of [24]; the plan's pump range of 1.1–1.6 µm excludes two-photon absorption across 2.67 eV (needs λ < 929 nm), so M3 cannot act on this dataset.

---
### [23] Ozlu et al., "Nonlinearity reversal in epsilon-near-zero indium tin oxide driven by few-cycle light pulse" — arXiv:2606.11567v1 (2026)

**Source read:** arXiv v1 full text + TeX source, including the Supplementary Information (pump FROG, transmitted probe spectra, 100 fs simulation).
**Role in this project:** Defines model M3 (two-photon interband absorption gated by Pauli unblocking, WP1/WP4) and is the Tier-2 dataset (few-cycle pump) for WP4.

**What it does.** A 270 nm ITO film is pumped at 800 nm with a compressed few-cycle pulse and probed at its ENZ wavelength with a 100 fs pulse. At low intensity $\Delta R/R_0<0$ and $\Delta T/T_0>0$, the usual intraband response. As the pump intensity rises toward 10 TW/cm², a contribution of opposite sign grows about 150 fs after the pump, and at the two highest intensities both $\Delta R/R_0$ and $\Delta T/T_0$ change sign. Its amplitude grows quadratically with intensity. The authors fit a three-level rate-equation model in which two-photon absorption (TPA) from the valence band into the bottom of the conduction band becomes allowed once one-photon intraband absorption has emptied those states.

**Model / method.**

Levels: valence band $|1\rangle$, lower conduction band $|2\rangle$ ($m_2^*=0.35\,m_0$), upper conduction band $|3\rangle$ ($m_3^*=0.50\,m_0$).

$$ \frac{dN_{1}}{dt} = -\beta I^{2} N_{1} f_{2} + \gamma_{21}\left(N_{2}P_{1} - N_{2,\mathrm{eq}} P_{1,\mathrm{eq}}\right) \qquad \text{(Eq. 1)} $$

$$ \frac{dN_{2}}{dt} = +\beta I^{2} N_{1} f_{2} - \gamma_{21}\left(N_{2}P_{1} - N_{2,\mathrm{eq}} P_{1,\mathrm{eq}}\right) - \sigma I N_{2} + \gamma_{32} N_{3} \qquad \text{(Eq. 2)} $$

$$ \frac{dN_{3}}{dt} = +\sigma I N_{2} - \gamma_{32} N_{3} \qquad \text{(Eq. 3)} $$

$$ P_1 = N_{1eq} - N_1 \qquad \text{(Eq. 4)} $$

with the Pauli-blocking factor $f_{2} = \max\left(0,\, 1 - N_{2}/N_{2,\mathrm{max}}\right)$. $N_{1,2,3}$ are level populations; $P_1$, $P_{1,\mathrm{eq}}$ the instantaneous and equilibrium valence-hole populations; $I$ the pump intensity; $\beta$ the TPA coefficient; $\sigma$ the intraband (one-photon) cross-section; $\gamma_{32}$ "phonon-mediated relaxation from the upper to the lower conduction band"; $\gamma_{21}$ net recombination in detailed-balance form. Degenerate limit: at equilibrium $N_2=N_{2,\mathrm{eq}}=N_{2,\mathrm{max}}$ and $P_{1,\mathrm{eq}}=0$, so $f_2=0$ and TPA is fully blocked until level 2 is depleted.

$$ \varepsilon(\omega,t) = \varepsilon_{\infty} - \frac{\omega^{2}_{p_2}(t)}{\omega^{2}+i\gamma_{D_2}\omega} - \frac{\omega^{2}_{p_3}(t)}{\omega^{2}+i\gamma_{D_3}\omega} \qquad \text{(Eq. 5)} $$

with $\omega_{p2}(t) = \omega_{p2,\mathrm{eq}}\sqrt{N_{2}(t)/N_{2,\mathrm{eq}}}$. $\omega_{p3}$ "is initialized to zero at equilibrium and grows as carriers are promoted into level $|3\rangle$"; its formula is not printed. The values in Fig. 4 are consistent with $\omega_{p3}^2=\omega_{p2,\mathrm{eq}}^2\,(N_3/N_{2,\mathrm{eq}})(m_2^*/m_3^*)$ (my check).

Optics and fit: transfer matrix for 270 nm ITO on 1.1 mm glass at 20°, TE, convolved with the 100 fs probe envelope. Levenberg–Marquardt fit to $\Delta R/R_0$ and $\Delta T/T_0$ at all five intensities at once, with one intensity scale factor per curve free within ±10 %.

Not in the model: electron temperature, thermalisation, lattice, temperature-dependent damping, pump propagation or field enhancement in the film, coherent pump–probe coupling.

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| Film / substrate | ITO 270 nm on 1.1 mm glass | Results |
| $\varepsilon_\infty$ | 3.9 (fixed) | Table 1 |
| Carrier density | 10²⁷ m⁻³ (fixed) | Table 1 |
| Drude damping $\gamma_{D2}=\gamma_{D3}$ | 1.24×10¹⁴ rad/s (fixed) | Table 1 |
| $m_2^*$ | 0.35 $m_0$ (fixed) | Table 1 |
| Mobility | not stated ($e/m^*\gamma \approx 40$ cm²/Vs, my calculation) | — |
| ENZ wavelength | 1240 nm (Re ε = 0); probe at 1225 nm is called the ENZ wavelength; Re ε ≈ 0.1, Im ε ≈ 0.307 at the probe | Fig. 1b; Fig. 4b,c |
| Gap | electronic gap of In₂O₃ "at most 2.67 eV"; optical gap ≈4 eV (Moss–Burstein shifted); two pump photons ≈3.1 eV | Sec. 2.2 |
| Fermi level | not stated | — |
| Pump | 800 nm, spectrum 700–900 nm, "sub-8 fs" (6.4 fs FWHM by SHG-FROG; Fig. 1a is labelled "10 fs beam"), 0° incidence, s-polarised, beam size 750 µm | Methods; SI S1 |
| Pump intensities | $I/I_{\max}$ = 0.32, 0.42, 0.50, 0.83, 1.00 with $I_{\max}\approx 10$ TW/cm²; fitted values 2.85, 3.74, 4.51, 8.83, 11 TW/cm² | Fig. 2; SI S2; Fig. 3 |
| Intensity definition | peak ("peak intensities"); incident or in-film not stated, and no in-film correction is described, so read as incident (my reading); horizontal error bars ≈±10 % | Fig. 2d |
| Pulse energy, fluence, repetition rate | not stated | — |
| Probe | 1225 nm, 100 fs, 20°, s-polarised, 50 µm | Methods |
| Observable, resolution | $\Delta R/R_0$, $\Delta T/T_0$, transmitted probe spectrum vs delay; resolution set by the 100 fs probe; zero delay = delay of maximum probe frequency shift | Results |
| Crossover intensity | 4.9 TW/cm² (4.91 in Fig. 2d) | Results; Fig. 2d |
| Peak response | $\Delta R/R_0 \approx +0.7$ at the highest intensity, ≈150 fs after the pump; low-intensity $\Delta R/R_0\approx-0.45$, $\Delta T/T_0\approx+0.4$; highest-intensity $\Delta T/T_0\approx-0.65$ (figure values my reading); recovery ≈1 ps | Results; Fig. 2a,b |
| $\beta$ | 1.13×10⁻¹² "cm⁴ s" (free) | Table 1 |
| $\sigma$ | 24.0 "cm²" (free) | Table 1 |
| $\gamma_{32}$ | 1.61×10¹² s⁻¹ (621 fs) (free) | Table 1 |
| $\gamma_{21}$ | 1.27×10¹³ s⁻¹ (79 fs) (free) | Table 1 |
| $N_{3,\max}$ | 3.00 (free) | Table 1 |
| $m_3^*$ | 0.50 $m_0$ (free) | Table 1 |
| Parameter uncertainties | none given | — |
| Model, highest intensity | $\omega_{p2}/2\pi$: 480 → ≈285 THz; $\omega_{p3}/2\pi$ → ≈457 THz; Re ε: 0.1 → +0.37 → ≈−0.97 → zero near 300 fs → +0.39 near 780 fs (my reading) | Fig. 4a–c |
| Model, 100 fs pump, 1.1 TW/cm² (same fluence) | Re ε stays positive (max ≈0.73); no reversal | SI S3 |

Units (as printed; looks inconsistent because the rates $\beta I^2$ and $\sigma I$ are in s⁻¹ only if $I$ is a photon flux, which would give $\sigma I\sim10^{33}$ s⁻¹ at 10 TW/cm²). With $I$ in W/cm² the numbers give $\sigma I = 2.4\times10^{14}$ s⁻¹ and $\beta I^2 = 1.1\times10^{14}$ s⁻¹ at 10 TW/cm², i.e. order-unity transfer in about 10 fs, which matches Fig. 4. So read $\sigma$ = 24.0 cm² J⁻¹ and $\beta$ = 1.13×10⁻¹² cm⁴ J⁻² s, with normalised populations (my reading; confirm against the authors' code).

**Key results.**
- Measured: complete sign reversal of both $\Delta R/R_0$ and $\Delta T/T_0$ at the two highest intensities (≈8.3 and 10 TW/cm²); at ≈5 TW/cm² only a dip appears (Fig. 2a,b, my reading).
- Measured: the reversal appears ≈150 fs after the pump has left, so it reflects real carriers "rather than … an instantaneous virtual-transition or coherent process".
- Fitted: the peak $\Delta R/R_0$ in the 200–400 fs window, after subtracting the lowest-intensity curve, follows a negative linear plus positive quadratic term in intensity (five points). The 4.9 TW/cm² "threshold" is where this fit changes sign. Coefficients are not tabulated; from Fig. 2d, ≈−0.067 per TW/cm² and ≈+0.014 per (TW/cm²)² (my digitisation).
- Fitted: the model reproduces sign reversal, peak position and later negative trough of $\Delta R/R_0$; it underestimates $\Delta T/T_0$ at the two highest intensities and relaxes too fast beyond ≈1.2 ps.
- Simulated: at equal fluence a 100 fs pump gives no reversal, so "peak intensity, not fluence, is the critical parameter".
- Simulated: Re ε completes a cycle (up, strongly negative, back above equilibrium) within ≈300 fs.

**Assumptions and limitations.**
- Authors: fixed damping; constant relaxation rates; no thermalisation and no coherent artefact, so the sharp early decay is missed; higher-order nonlinearities ignored.
- Three lumped levels replace the distribution; "Pauli blocking" is one occupancy factor (my reading).
- $N_{3,\max}$ is fitted but absent from Eqs. 1–4, and $N_{1,\mathrm{eq}}$ is not given; the equations are incomplete as printed (my reading).
- The pulse duration used in the model is not stated. The fluence-matched pair in the SI (1.1 TW/cm², 100 fs) equals 11 TW/cm² × 10 fs, not 6.4 fs (my reading).
- Only intensity was varied in the experiment, at one duration. The claim that intensity rather than fluence controls the effect rests on the model simulation, not on data (my reading).
- The 150 fs delay is attributed to thermalisation, which the model lacks. Fig. 3 uses a shifted delay axis (peak near 480 fs against ≈130 fs in Fig. 2) and its caption swaps the panels; whether a time offset was fitted is not stated (my reading).
- From Fig. 4a the fitted model needs $N_2+N_3\approx1.65\,N_{2,\mathrm{eq}}$ just after the pulse, a 65 % rise in conduction-band density (my calculation). Three of five intensity scale factors sit near the ±10 % bound (my comparison of Figs. 2 and 3).

**Use for this project.**
- M3 source term: generation rate $\beta I(t)^2 N_v[1-f(E_\mathrm{final})]$, with final states $2\hbar\omega-E_g\approx0.43$ eV above the conduction-band minimum for $E_g$ = 2.67 eV (my arithmetic). A distribution-resolved solver supplies the unblocking from $f(E,t)$ itself, so $\sigma$, $N_{2,\max}$ and $N_{3,\max}$ are not needed. Add recombination $\gamma_{21}NP$.
- Carry $\beta$ with the unit caveat and a wide prior; $E_g$ = 2.67 eV against 3.9 eV in [22].
- Dataset: five $\Delta R/R_0$ and $\Delta T/T_0$ traces (Figs. 2, 3). Data availability (verbatim): "Source data and simulation code are available from the corresponding author upon request." Code: "MATLAB code used for rate-equation modeling and TMM fitting is available upon request."
- Tier 2 only: 6.4 fs is about 2.4 cycles at 800 nm, outside the cycle-averaged Boltzmann description.
- TPA across 2.67 eV needs pump wavelength below 929 nm. In the plan's 1.1–1.6 µm range M3 is energetically closed, so M2 against M3 can only be tested with 800 nm data (my reading).
- Plan check: "about 5 TW/cm²" is a fit crossover; full reversal needs ≈8 TW/cm² or more. The reversal is in $R$ and $T$, not in a retrieved $\omega_p$; [24] shows that a transmittance sign change does not by itself prove a density change, although here $R$ rises while $T$ falls.
- The 150 fs lag of the reversal behind a 6.4 fs pump is close to the 120 fs lag that [22] treats as the signature of Auger excitation; timing alone does not separate M2 from M3 here (my reading).

---
### [24] Choi et al., "Pathway to optical-cycle dynamic photonics: extreme electron temperatures in transparent conducting oxides" — Nat. Commun. 17, 11129 (2026); arXiv:2512.24641v2

**Source read:** arXiv v2 full text + TeX source, main text and Methods only. Supplementary Sections 1–13 are not in the library; they hold the two-temperature equations, the temperature-dependent heat capacity and electron–phonon coupling, the fitted $T_0$, the Drude parameters, the layer thicknesses and the comparison with experiment. Published version not read.
**Role in this project:** Defines the "extreme electron temperature" explanation that model M0 (thermal electrons only) must cover in WP4.

**What it does.** This is a theoretical and numerical paper; it reports no new measurement. It asks whether thermalised electrons alone can produce oscillatory, sign-reversing optical transients on the scale of a few optical cycles. A spatially resolved two-temperature model feeds an electron-temperature-dependent Drude permittivity, for a 10 nm film with aluminium-doped zinc oxide (AZO) parameters on an inverse-designed Si/SiO₂ mirror that concentrates absorption. With the carrier density held fixed, the transmittance change reverses sign several times within tens of femtoseconds, because the falling plasma frequency and the rising damping push the transmittance in opposite directions. The real refractive-index change does not reverse sign unless the carrier density also changes; the authors obtain that with thermionic injection into a second, 2 nm low-doped layer.

**Model / method.**

$$ n(\omega, T_\mathrm{e}) = \sqrt{\varepsilon(\omega, T_\mathrm{e})} = \sqrt{\varepsilon_{\infty} - \frac{\omega_\mathrm{p}(T_\mathrm{e})^2}{\omega^2 + i \omega \Gamma(T_\mathrm{e})}} \qquad \text{(Eq. 1)} $$

$$ \omega_\mathrm{p}^2(T_\mathrm{e}) = \frac{e^2}{3 m^\ast \varepsilon_0 \pi^2} \int_{0}^{\infty} \frac{\mathrm{d} E}{1 + 2 C E} \left( \frac{2m^*}{\hbar^2 (E + C E^2)} \right)^{3/2} \left( -\frac{\partial f_\mathrm{FD}(E, \mu(T_\mathrm{e}), T_\mathrm{e})}{\partial E} \right) \qquad \text{(Eq. 2)} $$

(as printed; looks inconsistent because the bracket must be $k^2 = 2m^*(E+CE^2)/\hbar^2$, with $E+CE^2$ in the numerator, for the right-hand side to reduce to $n_e e^2/[\varepsilon_0 m^*(1+2CE_F)]$ at low temperature).

$$ \Gamma(T_\mathrm{e}) = \Gamma_0 \left[ 1 + \left(\frac{T_\mathrm{e}}{T_0}\right)^m \right] \qquad \text{(Eq. 3)} $$

$f_\mathrm{FD}$ is the Fermi–Dirac distribution, $\mu$ the chemical potential, $m^*$ the band-edge mass, $C$ the non-parabolicity, $\Gamma_0$ the room-temperature damping, $m$ a power-law exponent and $T_0$ a fitted temperature scale. The baseline is $m=2$; $T_0$ is fitted to the transmission data of Lustig et al. [20] (Supplementary Section 1).

Two chemical-potential models:

$$ \mu(T_\mathrm{e}) = \varepsilon_F \left[ 1 - \frac{\pi^2}{12} \left( \frac{T_\mathrm{e}}{T_\mathrm{F}} \right)^2 \right] \qquad \text{(Eq. 4, Sommerfeld)} $$

$$ n_{e} = \int_{0}^{\infty} f_\mathrm{FD}(E,\mu(T_\mathrm{e}),T_\mathrm{e})\, g(E)\,\mathrm{d}E \qquad \text{(Eq. 5, carrier conservation)} $$

with $g(E)$ the density of states (stray commas in the source equation removed). Eq. 4 does not conserve particles at high temperature (it overestimates the density), and the authors use it on purpose as a "hypothetical carrier injection" case, only in Figs. 1d and 2a. All transmittance results use Eq. 5.

$$ \ln\left(\frac{\mathcal{T}_\mathrm{mod}(\lambda, T_\mathrm{e})}{\mathcal{T}(\lambda)}\right) = \frac{\partial \ln \mathcal{T}(\lambda)}{\partial \varepsilon'} \Delta \varepsilon'(T_\mathrm{e}) + \frac{\partial \ln \mathcal{T}(\lambda)}{\partial \varepsilon''} \Delta \varepsilon''(T_\mathrm{e}) \qquad \text{(Eq. 6)} $$

$\alpha\equiv\partial\ln\mathcal{T}/\partial\varepsilon'$ and $\beta\equiv\partial\ln\mathcal{T}/\partial\varepsilon''$ are the refractive and absorptive sensitivities of the stack. A sign change of $\Delta\mathcal{T}$ occurs when $\alpha\Delta\varepsilon'$ and $\beta\Delta\varepsilon''$ have opposite signs and exchange dominance as $T_\mathrm{e}$ rises and falls.

Thermionic injection across the bilayer interface (modified Richardson–Dushman):

$$ \frac{dN_{\text{esc}}}{dt} = A R_1 R_2 (k_{\mathrm B} T_{\mathrm e})^2 \exp\left[ -\frac{eE_f - \mu(T_{\mathrm e}) + e\phi + \phi_{\mathrm{sc}}}{k_{\mathrm B} T_{\mathrm e}} \right] \qquad \text{(Eq. 7)} $$

with $A = 4\pi m^\ast h^{-3}$, $R_{1,2}$ the semi-axes of the emission region, $\phi$ the barrier (set to 0) and $\phi_\mathrm{sc}\approx aN_\mathrm{esc}e^2/R_1$ the space-charge barrier (Eq. 8); Eq. 9 gives the time-integrated yield.

Heating: the pump spectrum is passed through a transfer matrix to obtain depth-resolved absorption, which is the source of "the coupled two-temperature equations" with temperature-dependent chemical potential, electron heat capacity and electron–phonon coupling. Those equations and functions are in Supplementary Sections 5–6 and are not available here. The optics are adiabatic: the permittivity is evaluated at the instantaneous $T_\mathrm{e}$.

"Extreme" is defined as $T_\mathrm{e}$ approaching or exceeding the Fermi temperature, where $n_\mathrm{e}\lambda_T^3\ll1$ with $\lambda_T=\sqrt{2\pi\hbar^2/m^*k_\mathrm{B}T}$ and the distribution approaches Maxwell–Boltzmann.

**Numbers.**

| Quantity | Value | Where |
| --- | --- | --- |
| Material | Drude parameters "obtained from an aluminum-doped zinc oxide (AZO) film" (not ITO) | Sec. 2 |
| Absorber | 10 nm, $n_e$ = 7.35×10²⁰ cm⁻³ (optimised) | Sec. 2 |
| Stack | absorber on a five-layer Si/SiO₂ mirror; thicknesses not stated in main text; reference: same film on fused silica | Sec. 2 |
| $\varepsilon_\infty$, $\Gamma_0$, mobility | not stated in main text (Supplementary Fig. 3a) | — |
| ENZ / design wavelength | 1425 nm | Sec. 2 |
| Fermi temperature | 11584 K (≈1.0 eV) | Sec. 2 |
| $m^*$, $C$ | $m^*=0.24\,m_0$ quoted for the non-parabolicity sweep; $C$ not stated (sweep 0.1–0.6 eV⁻¹). $m^*=0.24\,m_0$ with $C=0.24$ eV⁻¹ reproduces 11584 K at this density (my back-calculation) | Sensitivity section; Fig. 6d |
| Absorptance | >60 % at ENZ, p-polarised; three times the single film | Sec. 2 |
| Pump | Gaussian, 10 fs, 1425 nm, p-polarised, 60°, 2 mJ/cm² (incident), peak ≈200 GW/cm² | Sec. 2; Fig. 1b |
| Peak $T_\mathrm{e}$ | "close to the Fermi temperature" in the cavity; ≈11.5×10³ K against ≈8.5×10³ K for the single film (my reading); maximum at ≈45 fs; ≈15×10³ K in the bilayer | Fig. 1c; Fig. 3; Fig. 4c |
| Damping law | $m$ = 2 baseline; $T_0$ not stated in main text; sweeps $m$ = 1.8–2.4, $T_0$ = 1.0–3.0×10³ K | Eq. 3; Fig. 6e,f |
| Electron–electron time | 5–6 fs; thermalisation taken as 10–20 fs | Results |
| Probe | p-polarised, 10°–50° | Fig. 3 |
| $\Delta\mathcal{T}/\mathcal{T}$ near ENZ (1431 nm) | three cycles of ≈20 fs (4–5 optical cycles), "approximately 20 %"; at 50°: ≈+13 %, −21 %, +14 % near 25, 45, 100 fs (my reading) | Fig. 3b |
| $\Delta\mathcal{T}/\mathcal{T}$ off ENZ (1075–1260 nm with angle) | opposite sign sequence, about ten times smaller (≈1–2 %) | Fig. 3c |
| Acceptor layer | 2 nm, 0.5×10²⁰ cm⁻³; $\Delta n/n\approx2$ %, three cycles of ≈20 fs at 1.8–2.0 µm | Fig. 4 |
| Shorter pump | 5 fs, 2 mJ/cm², acceptor 0.4×10²⁰ cm⁻³: first period 9 fs | Fig. 5c |
| 800 nm pump | needs "more than two orders of magnitude" higher fluence for the same $T_\mathrm{e}$ | Sec. 2 |

**Key results.**
- Simulated: with carrier density conserved, $\Delta\mathcal{T}$ reverses sign twice before $T_\mathrm{e}$ peaks. The cause is competition between $\Delta\omega_p$ and $\Delta\Gamma$, and the oscillation is tied to the extinction coefficient $\Delta k$, not to $\Delta n$.
- Simulated: under carrier conservation $\Delta n$ stays single-peaked even at $T_\mathrm{e}\approx T_\mathrm{F}$. A sign-reversing $\Delta n$ needs the density to change, giving competing $n_e$ and $m^*$ contributions to $\omega_p$.
- Stated: oscillations in $\Delta\mathcal{T}$ "can easily be misconstrued as evidence of oscillatory $\Delta n$".
- Asserted (evidence in Supplementary Section 1, not seen): with fitted $T_0$ the model accounts for the femtosecond transmission relaxation measured in ITO by Lustig et al. [20].
- Argued: a thermal model suffices because the electron–electron time is 5–6 fs, thermalisation completes in 10–20 fs, and earlier work [9, 14, 15] predicts weak non-thermal contributions. Non-thermal electrons "may affect the earliest response during the first tens of femtoseconds".
- Argued: multiphoton interband absorption is negligible here because the AZO gap (300–350 nm) needs four to five ENZ-wavelength photons.
- Stated: "we do not claim that this is the only possible explanation"; interband effects at higher intensity and intervalley transfer are left open.

**Assumptions and limitations.**
- Authors: thermalised electrons only; $m$ and $T_0$ are sample-specific and must be calibrated; idealised thermionic boundary; validity for shortened pumps "requires further investigation"; non-adiabatic effects may arise at these speeds.
- No measurement, AZO parameters, a 10 nm film in a cavity at 60°: the results are not a prediction for a 300 nm ITO film near normal incidence (my reading).
- The claimed thermalisation time (10–20 fs) is as long as the predicted features (≈20 fs) and the 10 fs pump, so the thermal assumption is weakest where the result lies (my reading).
- The response is a single-valued function of $T_\mathrm{e}$, so features on heating recur in reverse on cooling. [22] uses the absence of that symmetry in its $\omega_p(t)$ data as evidence against any equilibrium description (my reading).
- The reading of the plan's [14, 15] as predicting weak non-thermal effects differs from the plan, which says they attribute the real-permittivity change to non-thermal electrons.

**Use for this project.**
- M0 needs three ingredients: $\omega_p(T_\mathrm{e})$ from Eq. 2 (the thermal limit of the WP1 plasma-frequency integral, with the bracket corrected), $\mu(T_\mathrm{e})$ from Eq. 5, and $\Gamma(T_\mathrm{e})$ from Eq. 3 with $m$ and $T_0$ free. Without Eq. 3, M0 cannot produce any sign reversal. The plan lists only "coupling constants" as free parameters of M0.
- Never use Eq. 4 at high temperature; it changes the particle number.
- Hypothesis H2 must be stated for $\omega_p$ or Re ε. A fixed-density thermal model reverses the sign of the transmittance, so a sign flip in $R$ or $T$ alone does not falsify or support H2. The $\omega_p$ overshoot above equilibrium in [22] is what M0 cannot give.
- Validation targets: the $\Delta\mathcal{T}$ sign map (Fig. 2b) and the three-cycle trace (Fig. 3b), once the supplementary parameters are obtained.
- Data availability (verbatim): "The data generated in this study are provided in the article and Supplementary Information and can be reproduced using the Python code published on Figshare under DOI: https://doi.org/10.6084/m9.figshare.33283908." The code is deposited at the same DOI.
- Conflict with the plan: its background table lists this paper as an observation; it is a model prediction, and its only link to experiment is the unseen supplementary fit to [20].
- Conflict with [22]: that paper calls a threshold-like rise of damping "inconsistent with temperature-dependent scattering models" such as Eq. 3.

---
### [25] Segal et al., "Sub-cycle time-refraction at optical frequencies" — arXiv:2601.05566v1 (2026); also CLEO 2026, paper FM2E.2

**Source read:** arXiv v1 PDF, full (it has no methods section or supplementary material).
**Role in this project:** Motivation for sub-cycle index modulation, and an independent duration scan at fixed energy that invokes two-photon interband absorption (M3-type physics). The material is cadmium oxide, not ITO, so it is not a WP4 fitting dataset.

**What it does.** A 4 µm probe passes through a doped CdO film with its ENZ point at 4.3 µm while an intense 800 nm "modulator" pulse changes the index. The modulator is compressed from 41 fs to 10 fs at constant energy; 10 fs is shorter than one probe cycle (13 fs). Probe transmission falls, reflection rises and the transmitted spectrum blue-shifts, more so for shorter modulators. The authors attribute this to two-photon interband absorption that adds conduction electrons.

**Model / method.** There is no quantitative model. The only relation is the assumed index trajectory:

$$ \Delta n(t) \propto -\int_{-\infty}^{t} I^{2}_{\mathrm{modulator}}(\tau)\,\mathrm{d}\tau \qquad \text{(Fig. 1 caption, unnumbered)} $$

where $I_\mathrm{modulator}$ is the modulator intensity profile from FROG. The size of $\Delta n$ is estimated "assuming that the change to the refractive index of the CdO sample can be described by a simple shift of the ENZ wavelength". The sub-cycle character is inferred, not time-resolved: the modulator is shorter than the probe cycle and the index is assumed to follow the running integral of $I^2$.

**Numbers.**
- Sample: doped CdO (dopant not stated), 2.7 µm thick as printed, on sapphire (Fig. 1c label); ENZ at 4.3 µm (Fig. 1e); band gap ≈2.4 eV. Carrier density, mobility, $\varepsilon_\infty$: not stated.
- Probe: 4 µm (13 fs period), ≈100 fs, ≈0.01 mJ.
- Modulator: 0.8 µm (≈1.5 eV), ≈1 mJ, durations 41, 25, 17 and 10 fs at constant energy. Spot size, fluence, intensity, angles, polarisations: not stated.
- Observables: transmitted and Fresnel-reflected probe power, transmitted spectrum against delay; zero delay defined as a 10 % change in transmission. Temporal resolution: not stated.
- Index change: $\Delta n\approx-0.15$ during ≈0.77 probe cycle.
- Blue shift: to below 3 µm for the 10 fs modulator, "13% shorter than without the modulator". Ratio $\omega_\mathrm{shift}/\omega_0$ ≈ 1.118, 1.121, 1.139, 1.145 for 41, 25, 17, 10 fs (Fig. 4b, my reading). No error bars.

**Key results.**
- Measured: at fixed energy the frequency shift and the reflection rise are larger for shorter modulators; the transmission drop "does not depend on the duration of the modulator pulse" (Fig. 2a caption).
- Inferred: negative $\Delta n$, picosecond recovery and a pure blue shift indicate a multi-photon interband process, "most likely two-photon absorption", unlike the intraband response of [20].

**Assumptions and limitations.**
- Transmission falls to nearly zero for every duration, so it is saturated and cannot test the $I^2$ law; the shift ratios span 2.4 % with no stated uncertainty (my reading).
- The probe is 100 fs long, so the sub-cycle rise is not resolved (my reading).
- One sentence reads "in contrast to inter-band transitions where the relaxation can be as fast as tens of fs" (as printed; must mean intra-band).
- No data-availability statement.

**Use for this project.**
- Qualitative support for the M3 discriminator: at fixed energy, a response driven by $\int I^2\mathrm{d}t$ strengthens as the pulse shortens.
- Outside the WP1 range (CdO, 4 µm probe); the plan's background table should say so.
- Frequency shift against modulation time is a qualitative check for the WP3 Maxwell-kinetic solver, if CdO parameters can be obtained.

---

## Time-varying media: context and constraints
### [26] Galiffi et al., "Photonics of time-varying media" — Adv. Photonics 4, 014002 (2022)

**Source read:** arXiv:2111.08640v2 (18 Nov 2021) + TeX source; equation and figure numbers are those of the arXiv version. Read: Sec. 1 (introduction), Sec. 2 (temporal interfaces: boundary conditions, temporal slabs, anisotropic switching, dispersion, time interfaces in spatial structures), Sec. 3 (photonic time crystals, including topology, non-Hermiticity, time-modulated surfaces), Sec. 5 (experiments, all subsections), Sec. 6. Skipped: Sec. 4 (spatiotemporal metamaterials: nonreciprocity, synthetic motion, luminal amplification, space-time metasurfaces) apart from its opening paragraph.
**Role in this project:** Background for "why it matters" and for WP3: the standard temporal boundary conditions, time refraction/reflection and time-crystal formulas, and the review's statement that the material response time of ITO is an open question.

**What it does.** A review of wave physics in media whose parameters change in time. It derives the temporal boundary conditions and scattering coefficients for a non-dispersive step, extends them to slabs and periodic modulation (photonic time crystals), summarises the dispersive case, and surveys optical experiments, with epsilon-near-zero (ENZ) oxides as the main platform.

**Model / method.**
Temporal boundary conditions at a switch at $t_0$, from integrating Maxwell's curl equations (Eq. (3)):

$$\mathbf{B}(t_0^+)=\mathbf{B}(t_0^-),\qquad \mathbf{D}(t_0^+)=\mathbf{D}(t_0^-) \qquad (3)$$

Momentum conservation for a non-dispersive step $(\varepsilon_1,\mu_1)\to(\varepsilon_2,\mu_2)$ (Eq. (5)), and forward ($T$) and backward ($R$) amplitudes of $D$ (Eq. (6)), $Z_i=\sqrt{\mu_i/\varepsilon_i}$:

$$\omega_2\sqrt{\varepsilon_2\mu_2}=\omega_1\sqrt{\varepsilon_1\mu_1} \qquad (5)$$

$$T=\frac{Z_2+Z_1}{2Z_2},\qquad R=\frac{Z_2-Z_1}{2Z_2} \qquad (6)$$

These hold for spatial homogeneity, isotropy and instantaneous response. Double switch $1\to2\to1$ lasting $\Delta t_2$ (Eqs. (14), (15)):

$$|\mathbb{T}_2|^2=1+\frac{1}{2}\frac{(Z_1^2-Z_2^2)^2}{Z_1^2Z_2^2}\sin^2(\omega_2\Delta t_2),\qquad |\mathbb{R}_2|^2=\frac{1}{2}\frac{(Z_1^2-Z_2^2)^2}{Z_1^2Z_2^2}\sin^2(\omega_2\Delta t_2) \qquad (14),(15)$$

(as printed; composing Eq. (6) twice gives the prefactor 1/4, not 1/2 — my check. $|\mathbb T_2|^2=1+|\mathbb R_2|^2$ holds either way.) Periodic modulation at frequency $\Omega$ opens momentum gaps at $\omega=\Omega/2$, inside which waves grow in time.

Dispersive medium (Sec. 2.4): general kernel $\mathbf{D}(t)=\varepsilon_0[\varepsilon_\infty\mathbf{E}(t)+\int dt'\chi(t,t')\mathbf{E}(t-t')]$, and the Drude–Lorentz equation (Eq. (12)) with $\omega_p=\sqrt{Ne^2/(m_e\varepsilon_0)}$:

$$\frac{d^2\mathbf{P}}{dt^2}+\gamma\frac{d\mathbf{P}}{dt}+\omega_0^2\mathbf{P}=\varepsilon_0\omega_p^2\mathbf{E} \qquad (12)$$

For a time-switched $N$ (citing Bakunov et al. and [30]) the review states that the boundary conditions "turn out to be the continuity of $\mathbf{B}$, $\mathbf{D}$, $\mathbf{E}$ (and therefore $\mathbf{P}$) and $d\mathbf{P}/dt$", with four new frequencies $\pm\omega_{1,2}$.

**Numbers.**
- Speed: all-optical modulation acts on sub-ps scales, "very close to the time of an optical cycle (1-10 fs regime)"; a time crystal needs $\Omega=2\omega$ (Secs. 3, 5.1).
- Depth: for an instantaneous index change the time-refraction shift is proportional to $\Delta n/n$; ITO gives order-unity $\Delta n$ (Sec. 5.1.1, citing [1]).
- Time-refraction shifts: up to 58 THz in ITO/AZO films (Bruno et al.; Fig. 14 caption: 500 nm AZO, 105 fs pump, 770 GW/cm², degenerate at 1400 nm); 11.2 THz at 4 GW/cm² with gold antennas on ITO (Pang et al., Sec. 5.1.2; the Fig. 14 caption instead lists 0.5, 1 and 2 GW/cm², 50 fs pulses, 1304 nm). Whether intensities are incident or in-film is not stated.
- ITO: penetration depth about 100 nm near ENZ; 80 nm ITO metasurface (Liu et al.); 45 % reflectance change in Kretschmann geometry (Bohn et al.); recovery time 360 fs in unsaturated ITO (citing [1]) (Secs. 5.1.1, 5.2.1).
- Pump energy for time crystals: no figure given; only "Huge energy exchanges and stringent simultaneity of switching are typically required" for switching a whole medium (Sec. 2.5).

**Key results.**
- $D$ and $B$ are continuous at a time step; frequency changes, momentum is conserved; a backward wave appears (derived).
- With dispersion and a switched $N$, $E$, $P$ and $dP/dt$ are continuous as well (reported from [30] and Bakunov et al.).
- The ITO nonlinearity is attributed to intraband transitions that redistribute "the effective mass of conducting electrons" and shift the Drude plasma frequency (asserted; illustrated with Fig. 14a, taken from [17]).
- Photocarrier and nonlinear modulation are "generally limited by the pulse duration and not material-limited" (asserted).
- In the 80 nm ITO film a shorter pump gives a larger shift, and saturation sets in at lower energies for longer pulses; "the ultimate material response time for short and intense modulation constitutes still an open research question" (reported from Liu et al.; the conclusion is the review's).

**Assumptions and limitations.**
- Eqs. (5), (6), (14), (15) assume no dispersion; the review warns they change otherwise.
- Written before time diffraction [19] and single-cycle experiments [20]; neither is covered.
- No kinetic or microscopic model of the ITO response.
- The dispersive boundary conditions are quoted for a switched $N$, while the ITO mechanism is described as an effective-mass change; the review does not connect the two (my reading).

**Use for this project.**
- Cite for the standard boundary conditions and for the explicit "open research question" on the response time.
- It does not distinguish density change from mass change; the plan's claim cannot rest on [26].
- Eqs. (5), (6) are the non-dispersive limit, in which $E$ jumps at a step; in a Drude medium $E$ is continuous, so use them only as a limiting check of the WP3 solver.
- Conflict to note: "not material-limited" (Sec. 5.1) sits beside "open research question" (Sec. 5.1.1) in the same review.
- No data-availability statement (review).

---
### [27] Hayran et al., "ħω versus ħk: dispersion and energy constraints on time-varying photonic materials and time crystals" — Opt. Mater. Express 12, 3904 (2022)

**Source read:** arXiv:2208.00459v2 (7 Sep 2022) full text + TeX source; figures checked in the PDF. The published OME version was not compared.
**Role in this project:** Motivation ("why it matters") and WP3 context: the order-of-magnitude energy cost of fast, deep index modulation in transparent conducting oxides (TCOs), the slow/fast split of the TCO nonlinearity, and the time-varying Drude–Lorentz equation with a list of which microscopic parameters a modulation changes.

**What it does.** A perspective with two arguments. First, dispersion: a causal material responds non-instantaneously, so its time-varying susceptibility depends on two times; adiabatic "time-dependent parameter" models overestimate the effect of fast modulation, and the electric field stays continuous at a temporal boundary. Second, energy: changing a refractive index in time transfers energy ($\hbar\omega$), unlike spatial structuring, which only supplies momentum ($\hbar k$). Estimates are given for generic dielectrics, Kerr media and TCOs. For TCOs the authors conclude a photonic-time-crystal response is "possible but challenging, and limited to short times", and usually masked by a slower dynamic-grating response.

**Model / method.**
Two-time constitutive relation (Eq. (1)), causal: $\chi(t,t')=0$ for $t-t'\le0$:

$$\mathbf{P}(t)=\epsilon_0\int_{-\infty}^{\infty}\chi(t,t')\,\mathbf{E}(t')\,dt' \qquad (1)$$

Time-varying Drude–Lorentz model (Eq. (2)); $f$ oscillator strength ($f=1$ for Drude), $\omega_p=\sqrt{(N_1-N_2)e^2/\epsilon_0m^*}$ with $N_1$ ($N_2$) the density in the lower (excited) level, $\hbar\omega_0=E_{\rm gap}$ (Penn gap); low-density approximation (carriers independent):

$$\frac{\partial^2\mathbf{P}(t)}{\partial t^2}+\gamma(t)\frac{\partial\mathbf{P}(t)}{\partial t}+\omega_0^2(t)\mathbf{P}(t)=f\epsilon_0\omega_p^2(t)\mathbf{E}(t) \qquad (2)$$

Green's-function form (Eq. (3)):

$$\Big[\frac{\partial^2}{\partial t^2}+\gamma(t)\frac{\partial}{\partial t}+\omega_0^2(t)\Big]G(t,t')=\delta(t-t'),\qquad \chi(t,t')=f\omega_p^2(t)\,G(t,t') \qquad (3)$$

(as printed; looks inconsistent because solving Eq. (2) puts the source-time value $\omega_p^2(t')$ in the kernel, not $\omega_p^2(t)$. With $\omega_p^2(t)$, $P$ would jump at a step in $\omega_p$, whereas the paper's own boundary-condition argument and its FDTD, which integrates Eq. (2), keep $P$ and $E$ continuous — my reading.)

Energy estimate for TCOs, "in the most optimistic case" (Eq. (5)); $\delta U_{abs}$ absorbed energy density, $N$ carrier density, $\hbar v_Fk_F\sim E_F$:

$$\delta\varepsilon\sim\varepsilon_\infty\frac{\delta U_{abs}}{N\hbar v_Fk_F} \qquad (5)$$

**Numbers.**
- Required depth: relative permittivity change in time-crystal theory is typically ~10 % and up to 200 %; relative momentum-gap width is proportional to relative modulation strength (Sec. 1). Required speed: modulation frequency comparable to or above the signal frequency; time reversal needs modulation at $2\omega$.
- Generic solids (Sec. 3): $n$ = 1.4–3.4; Moss rule $n^4E_{\rm gap}\sim95$ eV; changing $n$ by 10–20 % costs ~$10^4$ J/cm³ (field > $10^8$ V/cm), enough for breakdown. If 0.1 % is dissipated, fewer than 100 modulation cycles heat the material by 1000 K.
- Electro-/acousto-optic: $\Delta n/n$ below a fraction of a percent; electro-optic speed at best 100s of GHz. Instantaneous Kerr effect: $\Delta n/n$ saturates below 1 % at well over 1 TW/cm².
- TCOs (Sec. 3): observed $\Delta n/n$ ~100 %; response "few 100's of fs", stated to be "two-to-three orders of magnitude too slow" for $2\omega$ modulation.
- Eq. (5) inputs: $\varepsilon_\infty\sim1$, $E_F\sim1$ eV, $N$ = few × $10^{20}$ cm⁻³; result: $\delta\varepsilon\sim1$ needs more than 10 J/cm³ absorbed.
- Fig. 3 (lossless Drude time crystal): $\omega_p^2(t)=\omega_{p0}^2+\Delta\omega_p^2(\cos\omega_{\rm mod}t-1)$, $\omega_{p0}=1.225\,\omega_p$, $\Delta\omega_p=0.707\,\omega_p$, $\omega_{\rm mod}=4\omega_p=10^{13}$ rad/s. Gap/mid-gap ≈ 8 % at $\Delta\epsilon(\omega_{p0})=0.333\,\epsilon_0$ and ≈ 15 % at $0.667\,\epsilon_0$ (read from figure); both are stated to need more than 50 TW/cm³ of supplied power density, for $\epsilon_\infty=1$, $N=10^{20}$ cm⁻³, $E_F=1$ eV, $\omega_p=2.5\times10^{12}$ rad/s. Text: "tens of TW/cm³ … for THz-scale modulations". Units are power per volume.
- Third-order response of TCOs (Sec. 3): ratio of slow (dynamic-grating) to fast ($2\omega$, time-crystal) nonlinearity ~$2\tau_{\rm el}/\tau_s$, with $\tau_{\rm el}$ (electron cooling to the lattice) a few hundred fs and $\tau_s$ (momentum scattering) a few fs; slow ≈ 100 × fast. General form $2T_1/T_2$.
- Second-order alternative: $\chi^{(2)}>100$ pm/V in III–V crystals, $\Delta n\sim\tfrac12\chi^{(2)}E_{2\omega}$, $\Delta n/n$ of a few percent with < 100 GW/cm² intracavity (a few GW/cm² incident).
- Fig. 1 (FDTD, Lorentz medium, $\omega_0(t)$ modulated between about 430 and 470 THz): $\omega_p=2\pi\cdot400$ THz, $\gamma=2\pi\cdot2$ THz, signal 390 THz, modulation 8 THz (slow) or 195 THz (fast). Fig. 2 (time switch): $\omega_p/2\pi$ 400 → 500 THz, $\omega_0=2\pi\cdot430$ THz, $\gamma=2\pi\cdot0.02$ THz, signal 400 THz with 8.8 THz bandwidth; non-dispersive comparison $\chi$ 6.42 → 10.05.

**Key results.**
- For modulation at half the signal frequency the adiabatic approximation overestimates the generated harmonics by more than 100 % (simulated, Fig. 1f).
- In a dispersive medium $E$ is continuous at an instantaneous switch; a non-dispersive model overestimates time reflection and gap width (argued from Eq. (1); simulated, Fig. 2).
- $\delta\varepsilon\sim1$ in a TCO costs > 10 J/cm³ (order-of-magnitude estimate).
- For large $\Delta n/n$ "only individual or infrequent changes" of the index may be feasible; the absorbed energy a material survives limits the number of modulation periods and the gap width (asserted).
- In TCOs a $\chi^{(3)}$ experiment is "always dominated" by the slow dynamic grating; the two responses are "externally indistinguishable" (asserted).

**Assumptions and limitations.**
- All energy figures are order-of-magnitude; $\varepsilon_\infty$ is set to ~1.
- Fig. 3 is for $\omega_{\rm mod}=10^{13}$ rad/s, far below optical frequencies. The upper axis is consistent with $\delta U_{abs}\,\omega_{\rm mod}$ (107 TW/cm³ at $\delta\varepsilon=0.667$); the conversion is not stated, and the axis ticks are unevenly spaced (my reading).
- $\omega_p=2.5\times10^{12}$ rad/s does not correspond to $N=10^{20}$ cm⁻³ (about $5.6\times10^{14}$ rad/s for a free-electron mass) (my check).
- Eq. (4) of the paper (static susceptibility) carries a prefactor $2\pi$ that Eq. (3) does not produce (as printed).
- Single time-refraction events, rise time and time diffraction are not analysed.

**Use for this project.**
- Energy scale to test the kinetic model against: > 10 J/cm³ per unit $\delta\varepsilon$, i.e. about 0.1 mJ/cm² absorbed in a 100 nm film (my arithmetic).
- Eq. (2) names every modulated parameter. The text gives two routes to changing $\omega_p$: carrier density $N_1-N_2$, and effective mass through heating in a non-parabolic band. Both enter only as $\omega_p(t)$ in one equation, so [27] does not distinguish a density change from a mass change. It distinguishes $\omega_p$ modulation from $\omega_0$ modulation.
- Conflict: for a density change [27] takes the susceptibility "scaled at each time instant" by the instantaneous density, i.e. the $A(t)$ kernel, while [30] derives the $A(\tau)$ kernel and says the two differ.
- Tension with [19]: in [27] the strong TCO response integrates absorbed energy and decays with $\tau_{\rm el}$ (few hundred fs), so its rise should follow the pump; a 1–10 fs rise inferred under a 225 fs pump [19] does not fit that picture (my reading). [27] does not say the origin of the response time is unknown.
- Mapping to the model family (my reading): the weak, fast "ballistic motion … inside the non-parabolic conduction band" response belongs to the coherent Tier 2 level; the strong, slow one is carrier heating (M0/M1).
- Data availability (quoted): "Data underlying the results presented in this paper are not publicly available at this time but may be obtained from the authors upon reasonable request."

---
### [28] Saha et al., "Photonic time crystals: a materials perspective [Invited]" — Opt. Express 31, 8267 (2023)

**Source read:** published PDF, full (7 pages, 8267–8273), figures checked.
**Role in this project:** Background ("why it matters"): states the modulation depth and speed a photonic time crystal (PTC) needs, and the materials-community view of what sets rise and relaxation times in transparent conducting oxides (TCOs).

**What it does.** A short perspective, without equations or new data, on materials whose refractive index can change by a large amount fast enough for a PTC at optical frequencies. It argues that optically pumped TCOs near epsilon-near-zero (ENZ) are the most promising and names two open problems: the fixed relaxation time and the pump power.

**Numbers.**
- Requirement (Sec. 1): index modulated "by an appreciable amount (in the range 0.1-1)" on "a time scale of a single optical cycle", to get a detectable time-reflected signal and/or wide momentum gaps.
- Pump energy (Sec. 4): "millijoules of energy pumped into optical systems in femtosecond cycles, resulting terawatts of peak power" (citing Dikopoltsev et al.). Units are energy and power, not fluence or intensity; no area or volume given.
- Interband pumping of TCOs (Sec. 2.2): index decrease of 0.1 in Al-doped ZnO (AZO) and 1.6 in undoped ZnO. Rise "instantaneous (follows the pump pulse)", "as fast as a few femtoseconds". Relaxation by trap-assisted (Shockley–Read–Hall) recombination over sub-picosecond to a few picoseconds; sub-picosecond in AZO near ENZ (Fig. 2a).
- Intraband pumping of TCOs (Sec. 2.2): recovery "a few hundred femtoseconds" (citing [1], [2] and Yang et al.).
- Fig. 2b (from Lustig et al., CLEO 2021): an 8 fs pump gives a much sharper transition than a 30 fs pump.
- Fig. 1d (ITO, from [1]): Re(n) rises from about 0.4 to about 1.1, saturating above roughly 100–150 GW/cm² incident intensity (read from figure, approximate).
- Damage (Sec. 4): unity-order index changes in oxides at powers 10× below damage threshold; TiN withstands "powers over 1 TW" (as printed; no area).
- Slower alternatives (Sec. 2.1): thermal phase change in nanoseconds (GSST, Δn = 0.5–3, 30–100 ns); electrical gating in picoseconds.

**Key results.**
- Thermal, electrical and mechanical modulation are too slow for an optical PTC (asserted).
- In TCOs the excitation speed is "limited only by the temporal pulse width of the excitation laser" (asserted).
- The bottleneck is relaxation, which "mostly depends on electron recombination and lattice cooling, which are intrinsic properties of the materials" (asserted).
- Proposed ways round it: two-colour (interband plus intraband) pumping, demonstrated in AZO by Clerici et al.; and stimulated emission from an inverted conduction-band population, of which the authors write "It is not at all clear that this idea is physical" (proposed).
- ENZ operation and antenna coupling lower the power needed (asserted).

**Assumptions and limitations.**
- Qualitative; no model, no energy-density estimate, no rise-time analysis.
- The intraband mechanism is described as a sub-bandgap pump that "energizes conduction electrons, and changing the free carrier concentration" (as printed); no effective-mass or non-parabolicity argument is given.
- Fig. 2a caption cites the paper's Ref. [18] (a GSST paper) for AZO data; the citation number looks wrong (my reading).
- Data availability (quoted): "No data were generated or analyzed in the presented research."

**Use for this project.**
- Cite for the PTC requirement (Δn of 0.1–1 within one optical cycle) and for the view that the rise is pump-limited and relaxation is the limit. The project's central question tests the first half of that view.
- Conflict: [28] words the intraband response as a carrier-concentration change; [26], [27] and [29] attribute it to effective mass. This is the density-versus-mass ambiguity that WP3 must resolve.
- Conflict: intraband recovery of "a few hundred femtoseconds" versus the 10–30 fs recovery of the transmitted probe reported in [20] (see that digest).
- No statement that the origin of the response time is unknown, apart from the doubt about the stimulated-emission idea.

---
### [29] Lobet et al., "New horizons in near-zero refractive index photonics and hyperbolic metamaterials" — ACS Photonics 10, 3805 (2023)

**Source read:** arXiv:2306.01314v2 PDF (42 pages), not the ACS version. Read: abstract, introduction, "Nonlinear properties of NZI media and their application to all-optical switching", "ENZ media for time-varying photonics", conclusions. Skipped: light and thermal emission, sensing and magneto-optics sections.
**Role in this project:** Context for the Background: typical index changes, relaxation times and irradiances of ENZ oxides, and a 2023 statement that the sub-10 fs rise of ITO was unexpected.

**What it does.** A multi-author perspective. The nonlinear section treats ENZ as a dispersion (slow-light) enhancement of "slow", real-state nonlinearities; the time-varying section lists ITO experiments and modelling needs. It contains no equations relevant here.

**Numbers.**
- Mechanisms (nonlinear section): free-carrier generation (interband) blue-shifts the index curve; free-carrier redistribution (intraband) red-shifts it. Fig. 3a models the intraband case as an effective-mass change, 0.3 to 0.5 (units not stated), in a Drude medium with $\varepsilon_\infty=4$, $\tau=6$ fs, $N=8\times10^{20}$ cm⁻³.
- ENZ oxides: $\Delta n\sim0.1$–1 with ~1 ps relaxation; $n_{2,\rm eff}=\Delta n/I\sim0.1$–$5\times10^{-3}$ cm²/GW at 1100–1700 nm; 10–1000 GW/cm² used without permanent damage.
- Table 1 (near normal incidence): AZO $3.5\times10^{-4}$ cm²/GW, ~0.8 ps, pump 785 nm, probe 1258 nm, 100 fs pulses; ITO $1.80\times10^{-3}$ cm²/GW, ~1 ps, pump 1100 nm, probe 1250 nm, 150 fs pulses.
- Time refraction: up to 14.9 THz shift in a 620 nm ITO film (citing [17]). Time diffraction: rise time estimated below 10 fs (citing [19]).

**Key results.**
- The fast rise "sparked ongoing theoretical investigations on the nature of such unprecedented response times"; [14, 15] are cited for eightfold higher electron temperatures than in standard plasmonic materials (reported).
- Boundary conditions and conservation laws at temporal inhomogeneities are "a field of intense ongoing investigation", citing [30]; [31] is cited as a route beyond adiabatic approximations (asserted).

**Assumptions and limitations.**
- Values are quoted from other papers. For time crystals only "strong modulations at single-cycle timescales" is stated; no depth, rise-time or pump-energy requirement is quantified.

**Use for this project.**
- Citable statement that the origin of the sub-10 fs rise was open in 2023.
- Do not use Table 1 values as model inputs: its footnote says $n_{2,\rm eff}$ depends on wavelengths, pulse width, thickness, irradiance and angle.
- Intraband response is assigned to effective mass and interband to carrier density, but no current equation is given.
- The plan lists [29] as arXiv only; it is published (ACS Photonics 10, 3805–3820). First author confirmed.

---
### [30] Solís et al., "Time-varying materials in the presence of dispersion: plane-wave propagation in a Lorentzian medium with temporal discontinuity" — Photon. Res. 9, 1842 (2021)

**Source read:** arXiv:2103.06142v1 (10 Mar 2021) full text + TeX source. The published Photon. Res. version was not available; equation and figure numbers below are those of arXiv v1 and may differ in print. The arXiv title reads "in Presence of Dispersion" (no "the").
**Role in this project:** WP3 — the reference for temporal boundary conditions and frequency conversion when the plasma frequency of a dispersive medium is stepped by changing the *density of oscillators*; it is the "density-change" case against which the plan's mass-change current equation is contrasted, and an analytic test for the Maxwell solver.

**What it does.** A monochromatic plane wave (frequency $\omega_-$) travels in an unbounded Lorentz medium whose plasma frequency is stepped at $t=0$ from $\omega_{p-}$ to $\omega_{p+}$ by changing the volumetric density $N$ of oscillators; the resonance frequency $\omega_0$ and damping $\gamma$ stay fixed. Because the polarization obeys a second-order equation, momentum conservation gives two new frequencies, each with a forward and a backward wave. Four amplitudes need four temporal boundary conditions: continuity of $D$ and $B$ plus continuity of $P$ and $dP/dt$. Closed-form amplitudes are derived, re-derived by Laplace transform, and checked with FDTD; a transmission-line equivalent and the lossy case follow.

**Model / method.**
Constitutive equation with a time-dependent source coefficient (Eq. (3) lossless; Eq. (37a) with loss), $\omega_p=q_e\sqrt{N/(\epsilon_0 m_e)}$, $N(t)=N_-+(N_+-N_-)U(t)$, $U$ the unit step:

$$\frac{d^2P(t)}{dt^2}+\gamma\frac{dP(t)}{dt}+\omega_0^2P(t)=\epsilon_0\omega_p^2(t)E(t) \qquad (37\text{a})$$

Equivalent kernel (Eq. (5) and text below Eq. (4)), with $A(t)=\omega_p^2(t)$ and $\chi_n(t)=\omega_0^{-1}\sin(\omega_0 t)U(t)$:

$$P(t)=\int_{-\infty}^{t}h(t,\tau)E(\tau)\,d\tau,\qquad h(t,\tau)=\epsilon_0A(\tau)\chi_n(t-\tau) \qquad (5)$$

The response to a field impulse at $\tau$ depends on $N(\tau)$ only: "the new dipoles brought into the medium after the electric-field impulse … simply have no excitation to respond to". The authors state that this differs from the model of Mirmoosa et al. (their Ref. [42]), where $h(t,\tau)=\epsilon_0A(t)\chi_n(t-\tau)$.

Temporal boundary conditions (Eqs. (8), (10), from the initial-value theorem), in addition to $D(0^+)=D(0^-)$, $B(0^+)=B(0^-)$:

$$P(0^+)=P(0^-),\qquad \frac{dP}{dt}(0^+)=\frac{dP}{dt}(0^-) \qquad (8),(10)$$

$D$ and $P$ continuous imply $E$ continuous. In the non-dispersive limit $E$ and $P$ jump.

Momentum conservation $k_-(\omega_-)=k_+(\omega_+)$, lossless, $\epsilon_\infty=1$ (Eq. (12)); for vacuum before the step ($\omega_{p-}=0$) it reduces to Eq. (15):

$$\omega_-\sqrt{1+\frac{\omega_{p-}^2}{\omega_0^2-\omega_-^2}}=\omega_+\sqrt{1+\frac{\omega_{p+}^2}{\omega_0^2-\omega_+^2}} \qquad (12)$$

$$\omega_+^4-(\omega_-^2+\omega_0^2+\omega_{p+}^2)\,\omega_+^2+\omega_0^2\omega_-^2=0 \qquad (15)$$

Roots $\pm\omega_1$, $\pm\omega_2$ with $\omega_1<\omega_0<\omega_2$ when $\omega_-<\omega_0$. Fields after the step, $E_+=e^{-ikz}\sum_{l=1}^{2}(f_le^{i\omega_lt}+b_le^{-i\omega_lt})$ (Eq. (18a)), with $\epsilon_l\equiv\epsilon_+(\omega_l)$, $\omega_l=\sqrt{\epsilon_-/\epsilon_l}\,\omega_-$ and amplitudes (Eq. (20), incident amplitude 1; $+$ forward, $-$ backward):

$$[f_1,b_1]=\frac{\epsilon_2-\epsilon_-}{2\sqrt{\epsilon_-}(\epsilon_2-\epsilon_1)}(\sqrt{\epsilon_-}\pm\sqrt{\epsilon_1}),\qquad [f_2,b_2]=\frac{\epsilon_--\epsilon_1}{2\sqrt{\epsilon_-}(\epsilon_2-\epsilon_1)}(\sqrt{\epsilon_-}\pm\sqrt{\epsilon_2}) \qquad (20)$$

Transmission-line form (Eq. (36)), $L(t)=1/(\epsilon_0\omega_p^2(t))$, $C(t)=\epsilon_0(\omega_p(t)/\omega_0)^2$:

$$L(t)\frac{d^2P(t)}{dt^2}+\frac{1}{C(t)}P(t)=E(t) \qquad (36)$$

with the explicit remark "(i) $\omega_0$ is kept constant and (ii) there is no $\frac{dL(t)}{dt}\frac{dP(t)}{dt}$ term"; the step is modelled as switching a new $L_pC_p$ branch in parallel, not as changing the inductance of the existing branch. The inductor "represents mass".

**Numbers.**
- Worked example: vacuum to $\epsilon_+(\omega_-)=4$ with $\omega_0=2\omega_-$, giving $\omega_{p+}=3\omega_-$; new frequencies $\omega_1=0.54\,\omega_-$, $\omega_2=3.70\,\omega_-$ (Fig. 1).
- Gradual transition $N(t)\propto(1+\tanh Rt)/2$; FDTD converges to the step result at $R=10^5/T$, $T=2\pi/\omega_-$ (Sec. IV.A, Fig. 2a).
- Lossy examples: $\epsilon_+(\omega_-)=4-0.1i$, i.e. $\omega_{p+}\approx3\omega_-$, $\gamma=0.1\omega_-$ (Fig. 4a); underdamped $\gamma=0.5\omega_-$ (Fig. 5); the $\omega_2$ pair becomes critically damped at $\gamma=7.01\omega_-$ (Fig. 4b); overdamped $\gamma=7.3\omega_-$ (Fig. 6).
- Limit $\omega_0/\omega_-\to\infty$ at fixed $\epsilon_+(\omega_-)$: $[f_1,b_1]\to\frac{\sqrt{\epsilon_-}}{2\epsilon_1}(\sqrt{\epsilon_-}\pm\sqrt{\epsilon_1})$ (the non-dispersive result), $f_2=b_2=\frac{\epsilon_1-\epsilon_-}{2\epsilon_1}$, $\epsilon_2\to0$ (Sec. IV.A.1).

**Key results.**
- Switching $N$ leaves $P$, $dP/dt$, $D$, $B$ and hence $E$ continuous (derived).
- One input frequency becomes two, each with forward and backward waves (derived; FDTD-confirmed).
- For $\omega_0\gg\omega_-$ at fixed $\epsilon_+(\omega_-)$ the $\omega_1$ waves reproduce the non-dispersive time-refraction/reflection amplitudes; the $\omega_2$ waves see $\epsilon\to0$, carry no power, and die out with infinitesimal loss, which is how the non-dispersive jump in $E$ is recovered (derived; App. C).
- With loss the new frequencies are complex: waves decay in time, not in space (derived; FDTD-confirmed).

**Assumptions and limitations.**
- Unbounded homogeneous medium, monochromatic plane wave, prescribed step (or tanh) in $N$; $\omega_0$, $\gamma$, mass constant; $\epsilon_\infty=1$. Eq. (20) is written for general $\epsilon_-$; the examples, the Laplace derivation and the lossy case take vacuum before the step ($\omega_{p-}=0$).
- Anomalous-dispersion interval not treated; energy balance deferred to a later paper.
- For a decrease of $\omega_p$ ($\omega_{p+}<\omega_{p-}$) the equivalent added branch is "non-Foster" (negative $L_p$, $C_p$); the authors defer its discussion to a future study (Sec. V), while stating in the Conclusions that the analysis "also holds for the inverse transition to free-space, or any other transition". The intraband response of ITO is such a decrease.
- The Drude case ($\omega_0=0$) is not treated explicitly, and an effective-mass change is never considered (my reading).
- In the paper $\mathbf J$ (Eq. (22b)) is a separate source current, dropped in the "source-free case"; the material current is $dP/dt$.
- Typos in arXiv v1: the boundary condition is printed as $\frac{dP(0^+)}{dt}=\frac{P(0^-)}{dt}$; "$\sum_l(f_l+b_l)=0$" in Sec. IV.A.1 contradicts row 1 of Eq. (19), which gives 1.

**Use for this project.**
- What it establishes for the plan's claim: the governing equation and boundary conditions for a *density* step, and that a second kernel, $A(t)$ instead of $A(\tau)$, is a different model. It does not discuss a mass change. Reading the absent $\dot L\,\dot P$ term as the mass term is the project's extrapolation: with $L\propto m^*/N$, $\frac{d}{dt}(L\dot P)=E$ would add $(\dot L/L)\dot P$, which is $(\dot m^*/m^*)\dot P$ at fixed $N$ (my reading).
- Drude limit for WP3 ($\omega_0\to0$, $J=dP/dt$): $dJ/dt+\gamma J=\epsilon_0\omega_p^2(t)E$ with $J$ continuous, i.e. new carriers start at rest (my reading). The plan's equation instead makes $J/\omega_p^2$ continuous at fixed $N$.
- Solver test, Drude limit of the density step (my derivation from Eqs. (12) and (A2) with $\omega_0=0$): roots $\omega_+=0$ and $\omega_+=\sqrt{\omega_-^2-\omega_{p-}^2+\omega_{p+}^2}$; electric-field amplitudes $[f,b]=\tfrac12(1\pm\omega_-/\omega_+)$ on the non-zero root; the zero-frequency mode carries magnetic field and current but no $E$. For $\omega_{p-}=0$ this is the classic sudden-plasma result (background, not from this paper).
- Discriminator (my derivation, not in any of [26]–[31]; lossless Drude, $\epsilon_\infty=1$): a mass step at fixed $N$ ($J/\omega_p^2$ continuous, the plan's equation) gives the same $\omega_+$ but $[f,b]=\tfrac12(1\pm\omega_+/\omega_-)$ and no zero-frequency mode. Frequency shift alone cannot separate the two models; time-reflected and transmitted amplitudes can. Both amplitude formulas were checked numerically (single-mode equations integrated through a fast tanh step); they are not from the literature.
- No data-availability statement.

---
### [31] Horsley et al., "Eigenpulses of dispersive time-varying media" — Phys. Rev. Lett. 130, 203803 (2023)

**Source read:** arXiv:2208.11778v1 (24 Aug 2022) full text + TeX source. The Supplementary Information is not in the arXiv files (only an uncaptioned SI figure), and the published PRL version was not available. The Drude-model derivation, all simulation parameters and the comparison with the adiabatic approach are in the SI and were NOT read.
**Role in this project:** WP3 / H3 — a frequency-domain operator method for reflection and transmission of a pulse from a thin Drude film with a prescribed time-varying plasma frequency; a non-adiabatic semi-analytic cross-check for the FDTD solver, and the place where the choice of two-time kernel (which the plan's current equation fixes) is made explicit.

**What it does.** The current response of a time-varying dispersive medium is written with a two-time conductivity. In the frequency domain this becomes an operator in which the time argument is replaced by $-i\partial_\omega$, with a fixed ordering. The usual Fresnel and slab formulas then hold with operators (matrices on a frequency grid) in place of the wavevector and impedance. Eigenvectors of these operators are "eigenpulses" whose spectrum is unchanged by the interaction; poles are bound states. The slab result is compared with COMSOL for a time-modulated Drude film at 45° incidence.

**Model / method.**
Two-time constitutive relation (Eq. (1)); causality requires $\sigma(t,t-t')=0$ for $t'>t$:

$$\boldsymbol{j}(t)=\int_{-\infty}^{\infty}dt'\,\sigma(t,t-t')\,\boldsymbol{E}(t') \qquad (1)$$

With $\tilde{\boldsymbol j}(\omega)=\int dt\,\boldsymbol j(t)e^{i\omega t}$ (Eq. (3)), all $\partial_\omega$ standing to the left of the frequency dependence ("anti-normal" ordering):

$$\tilde{\boldsymbol{j}}(\omega)=\hat{\sigma}(-i\partial_{\omega},\omega)\,\tilde{\boldsymbol{E}}(\omega) \qquad (3)$$

Drude model with time-varying plasma frequency (Eq. (4)); $1/\gamma$ is the collision time, $\mu=1$, and $\omega_0\equiv\omega_p(-\infty)$:

$$\hat{\sigma}(-i\partial_{\omega},\omega)=\omega_p^2(-i\partial_{\omega})\,\frac{i\epsilon_0}{\omega+i\gamma} \qquad (4)$$

Homogeneous medium, propagation along $x$ (Eq. (5)), $k_0=\omega/c$, $\eta_0=\sqrt{\mu_0/\epsilon_0}$, $\hat K_p^2=i\eta_0\hat\sigma k_0+k_0^2$:

$$H(x,\omega)=\exp(\pm i\hat{K}_p x)\,H_0(0,\omega) \qquad (5)$$

Half-space, in-plane wavevector $k_\parallel$, $k_x=[k_0^2-k_\parallel^2]^{1/2}$. TE: $\hat K_s^2=i\eta_0k_0\hat\sigma+k_0^2-k_\parallel^2$, $\hat Z_s=k_x^{-1}\hat K_s$. TM: $\hat K_p^2=i\eta_0\hat\sigma k_0+k_0^2-k_\parallel^2$, $\hat Z_p=k_x^{-1}[1+i\eta_0\hat\sigma k_0^{-1}]^{-1}\hat K_p$. For both (Eqs. (6), (7); TE defined on $E$, TM on $H$):

$$\hat r=(1-\hat Z)(1+\hat Z)^{-1},\qquad \hat t=2(1+\hat Z)^{-1} \qquad (6),(7)$$

Slab of thickness $d$, TM, $\hat A_\pm=1\pm\hat Z_p$ (Eqs. (9), (10)):

$$\hat r_{\rm slab}=\left[\hat A_+e^{i\hat K_pd}\hat A_--\hat A_-e^{-i\hat K_pd}\hat A_+\right]\left[\hat A_-e^{i\hat K_pd}\hat A_--\hat A_+e^{-i\hat K_pd}\hat A_+\right]^{-1} \qquad (9)$$

$$\hat t_{\rm slab}=4\hat Z_p\left[\hat A_+e^{-i\hat K_pd}\hat A_+-\hat A_-e^{i\hat K_pd}\hat A_-\right]^{-1} \qquad (10)$$

Alternative kernel: the authors note the formalism also works with the first argument taken as $t'$, i.e. $\sigma(t',t-t')$, "with only minor modifications"; that form needs the opposite (normal) ordering (their footnote [24]).

**Numbers.**
- Discretisation: $\partial_\omega$ is an $N\times N$ matrix on $N$ frequency points (finite difference or Fourier transform); $N$ not stated.
- Fig. 2 (half-space): $N(t)=\omega_p^2(t)/\omega_0^2$ dips from 1 to about 0.2 near $\omega_0t\approx130$; eigenpulse with $|r|=1.00$ near $\omega\approx0.52\,\omega_0$; Gaussian pulse centred near $0.5\,\omega_0$ (read from figure, approximate).
- Fig. 3 (slab): TM, 45° incidence, transmitted spectra around 200 THz for delays $-250$ to $250$ fs; line-outs at $\Delta t=50,100,120,150$ fs agree with COMSOL. Thickness, $\omega_p$, $\gamma$, modulation depth and timing: not stated in the main text (SI).
- Speed-up and accuracy figures: not stated in the main text.

**Key results.**
- Operator Fresnel and slab coefficients for a dispersive time-varying Drude medium (derived).
- Eigenpulses exist: an incident spectrum that is reflected with $|r|=1$ and unchanged shape, while a Gaussian pulse is broadened and reshaped (simulated, Fig. 2).
- The slab operator result is in "excellent agreement" with COMSOL (simulated, Fig. 3).
- TE and TM operators differ because $\hat\sigma$ and $\omega$ do not commute; at normal incidence $\hat r_s=-\hat r_p$ (derived).

**Assumptions and limitations.**
- $\omega_p(t)$ is prescribed; no self-consistent nonlinearity. Only $\omega_p$ varies; $\gamma$ constant.
- The slab solution needs $\hat K_p$ uniform through the film: no depth-dependent modulation (my reading).
- No background permittivity $\epsilon_\infty$ appears in $\hat K^2$ (my reading); ITO needs it added.
- The words "mass", "density" and "carrier" do not occur in the paper; the modulation is only "a time-varying plasma frequency", plotted as "$N(t)$".

**Use for this project.**
- Time-domain content of Eq. (4) (my derivation, not written in the paper): $j(t)=\epsilon_0\omega_p^2(t)\int_{-\infty}^{t}e^{-\gamma(t-t')}E(t')\,dt'$, equivalently
$$\frac{dj}{dt}+\Big(\gamma-\frac{1}{\omega_p^2}\frac{d\omega_p^2}{dt}\Big)j=\epsilon_0\omega_p^2(t)E .$$
With $\omega_p^2\propto N/m^*$ at fixed $N$ this is exactly the plan's WP3 equation, $\dot J+(\gamma+\dot m^*/m^*)J=\epsilon_0\omega_p^2E$; across a step $j/\omega_p^2$ is continuous and $j$ jumps.
- The footnoted alternative $\sigma(t',t-t')$ gives $\dot j+\gamma j=\epsilon_0\omega_p^2(t)E$ with $j$ continuous, which is the density-step model of [30] (my derivation).
- So [31] contains both kernels but does not name their physical meaning, does not compare their predictions, and calls the difference a minor modification. The plan's sentence that this literature "shows the distinction matters [30, 31]" is stronger than either paper; the mass/density reading is the project's.
- Benchmark: run the WP3 solver with prescribed $\omega_p(t)$ and uniform modulation against Eqs. (9)–(10), in both orderings. Parameters must come from the SI or the authors.
- For H3: this is a non-adiabatic alternative to "adiabatic permittivity + Fresnel", but it cannot represent depth-inhomogeneous excitation.
- No data-availability statement in the arXiv version.

---

## Non-equilibrium electron kinetics in metals
### [32] Rethfeld et al., "Ultrafast dynamics of nonequilibrium electrons in metals under femtosecond laser irradiation" — Phys. Rev. B 65, 214303 (2002)

**Source read:** abstract only (verbatim OpenAlex abstract in `references/text/32_rethfeld2002.abstract.txt`). Body NOT READ — paywalled, no PDF. Everything beyond the abstract is reconstructed from how [14], [36], [38], [39], [40] describe the paper and from background knowledge; each item is tagged with its origin. No equation below was copied from the paper: verify against the original before implementing.
**Role in this project:** Methodological ancestor of the WP1 solver: full Boltzmann collision integrals for laser absorption, electron–electron and electron–phonon scattering without relaxation times; the template for fidelity level B1 (momentum-resolved screened interaction with screening recomputed from $f$) and for absorption in photon-energy steps.

**What it does.** (abstract) Follows the distribution function of the free-electron gas of a metal during and after a sub-picosecond pulse of moderate intensity. Absorption by inverse bremsstrahlung, electron–electron thermalization and electron–phonon coupling are each "described by a full Boltzmann collision integral without using any relaxation-time approach"; the model "is free of phenomenological parameters". A "system of time- and energy-dependent integro-differential equations" is solved numerically for aluminium. Near the damage threshold the electron–lattice energy exchange agrees with the two-temperature model despite the non-equilibrium distribution; at low excitation the non-equilibrium distribution delays the cooling of the electron gas.

**Model / method.**

Structure (abstract; [39] SI Sec. S1 writes the same scheme): electrons $f$ and phonons $g$ evolve as

$$\frac{\partial f}{\partial t} = \left.\frac{\partial f}{\partial t}\right|_\text{el-el} + \left.\frac{\partial f}{\partial t}\right|_\text{el-ph} + \left.\frac{\partial f}{\partial t}\right|_\text{laser},\qquad \frac{\partial g}{\partial t} = \left.\frac{\partial g}{\partial t}\right|_\text{ph-el}$$

([39] SI Eq. 1, citing [32]). The equations are "time- and energy-dependent" (abstract), so the distributions are treated as isotropic (my inference; [39] SI: "assuming an isotropic momentum space"). [39] SI Eq. (23) summarizes each collision term of this model family as a sum over all partners of (matrix element)² × (Pauli collision functional) × (energy delta function), with "a screened Coulomb potential and a plane wave approach" for the matrix element.

Electron–electron term (background, unverified). Momentum-conserving two-particle scattering $\boldsymbol k,\boldsymbol k_1\to\boldsymbol k_2,\boldsymbol k_3$ in a parabolic band:

$$\left.\frac{\partial f(\boldsymbol k)}{\partial t}\right|_\text{el-el} = \frac{2\pi}{\hbar}\sum_{\boldsymbol k_1}\sum_{\boldsymbol k_3} |M_{ee}(\Delta k)|^2\,\mathcal F\;\delta\big(E(\boldsymbol k)+E(\boldsymbol k_1)-E(\boldsymbol k_2)-E(\boldsymbol k_3)\big),\qquad \boldsymbol k_2=\boldsymbol k+\boldsymbol k_1-\boldsymbol k_3$$

$$\mathcal F = f(\boldsymbol k_2)f(\boldsymbol k_3)[1-f(\boldsymbol k)][1-f(\boldsymbol k_1)] - f(\boldsymbol k)f(\boldsymbol k_1)[1-f(\boldsymbol k_2)][1-f(\boldsymbol k_3)]$$

$$|M_{ee}(\Delta k)|^2=\left(\frac{e^2}{\varepsilon_0\Omega}\,\frac{1}{\Delta k^2+\kappa^2}\right)^2,\qquad \kappa^2=\frac{e^2 m_e}{\pi^2\hbar^2\varepsilon_0}\int_0^\infty f(k)\,\mathrm dk$$

$\Delta k$ is the exchanged momentum, $\Omega$ the normalization volume, $\kappa$ a static screening wave number recomputed from the instantaneous distribution. The matrix-element form matches [40] Eq. (10). For a parabolic band the screening expression equals $(e^2/\varepsilon_0)\int D(E)\,(-\partial f/\partial E)\,\mathrm dE$ after integration by parts (my check), i.e. [40] Eq. (17) with its sign and units repaired.

Electron–phonon term (background, unverified). Golden-rule sum over phonon wave vectors with absorption and emission, Pauli factors for electrons and Bose factors $g$, $g+1$; matrix element of screened electron–ion type, recalled as $|M_{ep}(q)|^2=\frac{e^2}{2\varepsilon_0\Omega}\frac{\hbar\omega_q}{q^2+\kappa^2}$; one longitudinal acoustic branch; umklapp processes neglected ([37] says of the successor model [33] that it "ignores Umklapp processes"). A mirror equation drives $g(q)$.

Absorption term (background, unverified). Inverse bremsstrahlung: an electron exchanges momentum $\Delta\boldsymbol k$ with an ion and $\ell$ photons with the field,

$$\left.\frac{\partial f(\boldsymbol k)}{\partial t}\right|_\text{absorb}=\frac{2\pi}{\hbar}\sum_{\Delta\boldsymbol k}|M_{ei}(\Delta k)|^2\sum_{\ell}J_\ell^2\!\left(\frac{e\,\boldsymbol E_L\cdot\Delta\boldsymbol k}{m_e\,\omega_L^2}\right)\Big\{f(\boldsymbol k')[1-f(\boldsymbol k)]-f(\boldsymbol k)[1-f(\boldsymbol k')]\Big\}\,\delta\big(E(\boldsymbol k)-E(\boldsymbol k')-\ell\hbar\omega_L\big)$$

with $\boldsymbol k'=\boldsymbol k+\Delta\boldsymbol k$, $J_\ell$ the Bessel function of the first kind, $\boldsymbol E_L$ and $\omega_L$ the laser field amplitude and frequency, $M_{ei}$ a screened electron–ion matrix element. The argument is $\Delta k$ times the quiver amplitude $eE_L/m_e\omega_L^2$. Terms $\ell=\pm1$ give the one-photon steps; $|\ell|\ge2$ are simultaneous multi-photon processes. Support from the library: the abstract calls this a "plasmalike absorption term"; [38] calls it the "electron-ion-photon collision term"; [39] SI states that simultaneous two-photon absorption is included "as described by Rethfeld et al."

Numerical method: integro-differential equations in energy and time (abstract). Grid, quadrature, integrator, cost and conservation diagnostics: not known (not read).

**Numbers.**
- Material: aluminium (abstract).
- Recalled, unverified (background): 630 nm, 100 fs pulse of constant intensity, field amplitudes up to about 1.4×10⁸ V/m. Check before quoting.
- Thermalization and cooling times: not available from the abstract.

**Key results.**
- Absorption by free electrons is "well described by a plasmalike absorption term", with "good agreement of calculated absorption characteristics with values experimentally found" (abstract; simulated vs measured).
- Near damage threshold the electron–lattice energy exchange follows the two-temperature model (abstract; simulated).
- At low excitation the electron gas cools with a delay; the cooling time "depends thus on excitation parameters" (abstract; simulated).
- Thermalization time is defined as "the time after which the collective behavior of laser-excited electrons equals the thermalized limit" (abstract).
- The excited distribution has steps one photon energy wide; the paper introduced the function $\Phi$ in which a Fermi function is a straight line of slope ∝ 1/T ([40] Fig. 2 caption; [14] cites [32] for the multi-photon "shoulders").
- A classical field-acceleration drive term is adequate only when the energy gained from the field far exceeds one photon energy ([36] SI, citing [32]; [14] refers to [32] for the discussion).

**Assumptions and limitations.**
- Free-electron, parabolic single band; isotropic distributions (abstract and background).
- Static screening (background).
- No transport; one material (abstract).
- Later papers of the same group call such kinetic simulations "numerically rather expensive" ([38]) and "computationally unfeasible at scale" ([39]).

**Use for this project.**
- B1 in the plan is this formulation: $|M(\Delta k)|^2$ with $\kappa[f]$. After angular integration the kernel depends on the magnitudes of three momenta and does not reduce to convolutions, unlike the constant-matrix-element integral of [40] (my reading).
- The plan cites [14, 32] for "phonon- and impurity-assisted absorption in steps of the photon energy". [32] is ion-assisted inverse bremsstrahlung with a Bessel series. Check it against the absorption term of [14]: the Sivan-group form in [36] (Eqs. A3–A9) is a power-normalized golden-rule term without an explicit momentum-supplying partner. Decide which one WP1 implements.
- Size of the Bessel argument for ITO (my estimate; $m^*=0.35\,m_e$, 1240 nm, $\Delta k\sim k_\text{F}\approx3\times10^{9}$ m⁻¹ at $10^{27}$ m⁻³): ≈0.2 at a field equivalent to 10 GW/cm² in vacuum and ≈1.8 at 1 TW/cm². One-photon absorption suffices at the low end of the plan's range; simultaneous multi-photon intraband terms are not negligible at the top.
- The screening formula needs the Kane-band generalization $\kappa^2=(e^2/\varepsilon_0)\int D\,(-\partial f/\partial E)\,\mathrm dE$.
- "Delayed cooling at low excitation" is a qualitative validation target for the electron–phonon term.
- Obtain the paper before coding: prefactors, the phonon model and the numerical scheme are unverified here.

---
### [33] Mueller et al., "Relaxation dynamics in laser-excited metals under nonequilibrium conditions" — Phys. Rev. B 87, 035139 (2013)

**Source read:** abstract only (verbatim OpenAlex abstract in `references/text/33_mueller2013.abstract.txt`). Body NOT READ — paywalled, no PDF. Details beyond the abstract are reconstructed from [37], [38], [39], [40] and background knowledge and are tagged; verify before use. Authors: B. Y. Mueller, B. Rethfeld.
**Role in this project:** The version of the Rethfeld Boltzmann model with a realistic density of states; it is the reference solver behind the validation statements of [38], [39], [40] and the usual citation for the Fermi-liquid lifetime used in the relaxation-time baselines.

**What it does.** (abstract) Applies "complete Boltzmann collision integrals to describe the transient electron distribution due to excitation, thermalization, and relaxation", with "the density of states of real materials". It extracts "the electron thermalization time and the electron-phonon coupling under nonequilibrium conditions" for aluminium, gold and nickel, and finds that both depend strongly on material and excitation type and that "a nonthermalized electron gas couples differently to the phononic system as a thermalized one".

**Model / method.**
- Collision terms: laser excitation, electron–electron and electron–phonon, all as full collision integrals (abstract). [39] SI Eq. (23) gives the generic form used in this model family,
$$\left.\frac{\partial f(k)}{\partial t}\right|_\text{coll} = \sum_{\text{all }k'} M^2(k, k')\times\mathcal{F}(f(k), f(k'), f(k\pm k'))\times\delta(E(k), E(k\pm k')),$$
with isotropic momentum space, $M^2$ "calculated using a screened Coulomb potential and a plane wave approach" and $\mathcal F$ the Pauli collision functional; it refers to [33] for details.
- Excitation: an "electron-ion-photon collision term"; it produces a step structure under Pauli blocking, includes two-photon absorption and imprints DOS features on the excited distribution ([38] Secs. 2.2, 3.3).
- Screening parameter recomputed from the distribution at each time step ([40] cites [33] for its Eq. 17).
- Single effective band, free-electron-like dispersion combined with a DOS ([37]: "a simplified single-band free-electron-like model of the electrons"); the electron–phonon matrix element "ignores Umklapp processes" ([37]).
- [39] calls its own reference "a modified form" of this model "using a random-$k$ approximation", so [33] itself presumably does not use it (my inference).
- Fermi-liquid lifetime cited to [33] (with Pines–Nozières and Kaveh–Wiser) by [38], [39], [40]:
$$\tau_\text{E} = \tau_0 \frac{E_\text{F}^2}{\left(E-E_\text{F}\right)^2 + \left(\pi k_\text{B}T_\text{e}\right)^2}$$
([40] Eq. 19; [38] calls it "a simplified version from Fermi-liquid theory"). Whether [33] prints exactly this form: not verified.
- Numerical method, grids, conservation checks: not known.

**Numbers.**
- Materials: Al, Au, Ni (abstract).
- Al electron–phonon coupling at 2000 K as quoted by [37]: $G \approx 2.9\times10^{7}$ W m⁻³ K⁻¹ for thermalized electrons, "1.5x larger than" $1.9\times10^{7}$ W m⁻³ K⁻¹ for non-thermalized electrons with the same energy (exponents as printed in [37]; almost certainly 10¹⁷, since all other Al values there are of order 10¹⁷).
- Thermalization times, pulse parameters: not available from the abstract.

**Key results.**
- Thermalization time and electron–phonon coupling strength "strongly depend on the material properties and the excitation type" (abstract; simulated).
- A non-thermal electron gas couples to phonons differently from a thermal one; for Al the coupling is lower by a factor 1.5 at equal energy (abstract; number from [37]).
- The DOS shapes the excited distribution; [38] (Sec. 3.2, citing [33]) attributes the residual per-mille conservation error of its eTTM "probably" to neglecting this.

**Assumptions and limitations.**
- One effective band; plane-wave matrix elements; no umklapp ([37], [39]).
- [37] argues that its predictions "do not match experiment" for Al because of the simple electron–phonon matrix element.

**Use for this project.**
- The energy-space formulation with an arbitrary DOS is the natural route to a Kane band: replace the DFT DOS by $D(E)$ of $\hbar^2k^2/2m_0^*=E(1+CE)$ (my reading).
- A non-thermal distribution changes the electron–lattice energy flow by tens of percent in Al. The thermal-only model M0 will misestimate the recovery unless $Q_{e\text{-}ph}[f,T_l]$ is evaluated on the actual $f$, as WP1 plans.
- Source of the Fermi-liquid $\tau(E)$ used in [38]–[40]; [40] finds the textbook prefactor 2.8 times too large (0.284 fs against a fitted 0.1 fs) when compared with the full integral.
- Umklapp is neglected here ([37]); the plan takes the electron–phonon term, including Umklapp, from [14] instead.
- Obtain the paper if B1 is to follow this implementation; none of its equations is verified here.

---
### [34] Carpene, "Ultrafast laser irradiation of metals: Beyond the two-temperature model" — Phys. Rev. B 74, 024301 (2006)

**Source read:** abstract only (verbatim OpenAlex abstract in `references/text/34_carpene2006.abstract.txt`). Body NOT READ — paywalled, no PDF. Details beyond the abstract are reconstructed from [9], [14], [15], [38], [39], [40] and background knowledge and are tagged; verify before use.
**Role in this project:** Origin of the extended two-temperature model (eTTM) that [38] refines and [39] replaces; source of the two relaxation-time formulas used in the WP2 baselines.

**What it does.** (abstract) Extends the two-temperature model so that it incorporates "the transient electronic excitation process induced by an ultrafast laser pulse and the subsequent relaxation via electron-electron and electron-phonon scattering". The extension makes the thermal response depend on photon energy and yields "distinct heating terms for the electronic and the lattice temperatures", illustrated by "explanatory simulations".

**Model / method.**
- Three subsystems: phonons, thermal electrons, and non-thermal electrons and holes created by the laser ([38] Sec. 2.1 and Fig. 1).
- A constant density of states around the Fermi edge is assumed, which allows "large parts of the model" to be solved analytically ([38]); single-photon absorption only ([39]).
- The laser-created deviation from the Fermi function has the step form later written as [38] Eq. (3), and decays as
$$\Delta f'(t, t', E) = \Delta f_\text{L}(t', E) \exp\left(-\frac{t-t'}{\tau_\text{el*-el}}-\frac{t-t'}{\tau_\text{el*-ph}}\right)$$
([38] Eq. 5, "According to" [34] and Tsibidis). Energy integrals run only over $[E_\text{F}-E_\text{pt},\,E_\text{F}+E_\text{pt}]$ ([38] Sec. 2.3).
- Relaxation times taken from this paper by later work: the Fermi-liquid prefactor
$$\tau_{0} = \frac{128}{\sqrt{3}\pi^2\omega_\text{p}}$$
([40] Eq. 20; [38] Eq. A6; both cite [34] with Pines–Nozières), and the non-thermal electron–phonon time
$$\tau_\text{el*-ph} = \tau_\text{mfp}\frac{\hbar\omega}{k_\text{B}\theta_\text{D}}$$
([39] Eq. 15, citing [34] and Groeneveld et al. 1995), which assumes that each collision, occurring every $\tau_\text{mfp}$, removes one maximum phonon energy $k_\text{B}\theta_\text{D}$ from an electron with excess energy $\hbar\omega$.
- Energy dependence of the electron–electron time in the original: recalled as the zero-temperature Fermi-liquid form $\tau_0E_\text{F}^2/(E-E_\text{F})^2$, which diverges at the Fermi edge (background, unverified; [38] says its Eq. 6 "prevents a singularity at the Fermi edge" present in Tsibidis's version).
- [15] lists [34] among studies that estimate $\tau_{e\text{-}e}$ with the Gurzhi formula ([15] Eq. B7) and that give the non-thermal distribution the fixed profile $A(t)\{f^T(E-\hbar\omega)[1-f^T(E)]-f^T(E)[1-f^T(E+\hbar\omega)]\}$ with $A(t)$ set by the non-thermal energy.

**Numbers.** None can be given from the abstract. Values that later papers compute with the formulas above: $\tau_0=0.329$ fs and $\tau_\text{el*-ph}=299.22$ fs for Al at 1.55 eV ([38] App. A.3).

**Key results.**
- Separate, delayed heating terms for electrons and lattice replace the instantaneous source of the two-temperature model (abstract; derived).
- The thermal response depends on photon energy (abstract; simulated).

**Assumptions and limitations.**
- Constant DOS; integration window of one photon energy; one-photon absorption ([38], [39]).
- Energy and particle number are conserved only approximately in this model family ([38] Fig. 6; [39] SI S3).
- [14] calls the relaxation-time closure of such models "somewhat ambiguous" in its choice of energy-averaged decay coefficients.

**Use for this project.**
- Baseline lineage: [34] → Tsibidis 2018 → [38] → [39]. Implement [39], not [34].
- The constant-DOS, fixed-chemical-potential assumptions fail for ITO: [15] states the fixed-profile ansatz "is not valid for ITO since the chemical potential decreases significantly with $T_e$".
- [9] argues the opposite simplification for ITO: with $E_\text{F}\approx\hbar\omega\approx1$ eV the thermalization time is "negligible compared to the pump pulse duration", so athermal electrons are neglected. [15] instead finds the non-thermal distribution "essential to properly capture the nonlinear response of ITO". H1 tests this disagreement.
- $\tau_0$ and $\tau_\text{el*-ph}$ are free-electron estimates; [40] had to reduce $\tau_0$ by a factor 2.8.

---
### [35] Della Valle et al., "Real-time optical mapping of the dynamics of nonthermal electrons in thin gold films" — Phys. Rev. B 86, 155139 (2012)

**Source read:** abstract only (verbatim OpenAlex abstract in `references/text/35_dellavalle2012.abstract.txt`). Body NOT READ — no PDF arrived in this session (the manifest lists an open copy at hal.science/hal-02394961). Details beyond the abstract come from [14], [15] and background knowledge and are tagged; verify before use. Authors: G. Della Valle, M. Conforti, S. Longhi, G. Cerullo, D. Brida.
**Role in this project:** Experimental precedent for separating non-thermal from thermal electron contributions in a pump–probe signal with a rate-equation model; the place where the three-equation "extended two-temperature model" with a non-thermal energy reservoir is applied to data.

**What it does.** (abstract) Measures the nonlinear optical response of thin gold films with "high temporal resolution (∼15 fs) and broad spectral coverage" and models the data "quantitatively ... without free parameters using an extended version of the two-temperature model". The combination allows "tracking the evolution from a purely nonthermal electron distribution towards a fully thermalized one".

**Model / method.**
- From the library: [14] groups this paper with Schoenlein et al. (1987), Carpene [34] and Dubi–Sivan [36] as noble-metal models in which "it is customary to add a dynamical equation for the total non-thermal energy ... usually done within the RTA, requiring a somewhat ambiguous choice of an energy-averaged decay coefficients of the non-thermal energy to the (thermal) electron and phonon subsystems".
- Recalled structure (background, unverified; not copied from the paper). A non-thermal energy density $N(t)$ is added to the two temperatures:
$$\frac{\mathrm dN}{\mathrm dt}=-aN-bN+P_a(t),\qquad C_e\frac{\mathrm dT_e}{\mathrm dt}=-G\,(T_e-T_l)+aN,\qquad C_l\frac{\mathrm dT_l}{\mathrm dt}=G\,(T_e-T_l)+bN$$
with $P_a$ the absorbed power density, $a$ the rate at which non-thermal electrons heat the thermal electron gas, $b$ their rate of direct energy transfer to the lattice, $C_e$, $C_l$ the heat capacities and $G$ the electron–phonon coupling. Symbols and signs are my recollection.
- Recalled optical model (background, unverified): the probe response is computed from the change of interband absorption of gold caused by two parts of the occupation change — a thermal smearing of the Fermi edge set by $T_e$, and a non-thermal step-like part extending one pump-photon energy above and below the Fermi level, with amplitude set by $N(t)$. [15] describes this fixed-profile ansatz, $f^{NT}(E,t)=A(t)\{f^T(E-\hbar\omega)[1-f^T(E)]-f^T(E)[1-f^T(E+\hbar\omega)]\}$ with $A(t)$ from the non-thermal energy, citing Sun et al. (1994), [34] and a later paper of the same group (Zavelani-Rossi et al. 2015), not [35] itself.

**Numbers.**
- Time resolution ≈15 fs; sample: thin gold films (abstract).
- Film thickness, pump and probe wavelengths, fluence, rates $a$ and $b$, heat capacities, $G$, extracted thermalization time: not available from the abstract.

**Key results.**
- A parameter-free extended two-temperature model reproduces the measured broadband transient response (abstract; measured vs simulated).
- The data resolve the passage from a purely non-thermal to a thermalized electron distribution (abstract; measured).

**Assumptions and limitations.**
- Non-thermal electrons enter only through one scalar, their total energy; their spectral shape is prescribed, not evolved (background and [14], [15]).
- Gold, interband probe: the non-thermal signature is read out through d-band to Fermi-level transitions (background). In ITO the near-infrared response is intraband, so this read-out does not carry over (my reading).
- [15] states that the fixed-profile ansatz is "a decent approximation for noble metals because the chemical potential is nearly $T_e$-independent" and "not valid for ITO".

**Use for this project.**
- The single rate equation for a non-thermal energy density with relaxation-rate closures belongs to this model family (Sun et al. 1994, [34], [35]), not to [38]. If WP2 wants that three-equation baseline, cite these and take the equations from the paper, not from this digest.
- Cautionary point for H1: a scalar non-thermal reservoir was sufficient for gold; [15] argues it is insufficient for ITO because the permittivity depends on the shape of $f$ and on the shifting chemical potential.
- Raw data availability: not known. Fetch the open copy (HAL) before using any number.

---
### [36] Dubi et al., "'Hot' electrons in metallic nanostructures: non-thermal carriers or heating?" — Light Sci. Appl. 8, 89 (2019)

**Source read:** arXiv:1810.00565v4 (21 Jul 2019), full text with its Supplementary Information appendix, + TeX source. Authors: Y. Dubi, Y. Sivan. Published version not compared; equation numbers are those of the arXiv PDF.
**Role in this project:** Ancestor of the Sivan-group kinetic model [14, 15] behind M1: golden-rule absorption term, relaxation-time electron–electron term, and the split of $f$ into thermal and non-thermal parts.

**What it does.** Solves a Boltzmann equation for the conduction electrons of a silver nanoparticle under continuous-wave (CW) illumination, coupled to energy balances for electrons and lattice with heat loss to the environment, and finds the unique steady state. Electron and phonon temperatures are nearly equal, and the power going into non-thermal carriers is 10⁻⁸ to 10⁻⁵ of the power that heats the particle.

**Model / method.**

$$\frac{\partial f (\mathcal E, T_e,T_{ph})}{\partial t} = \left(\frac{\partial f}{\partial t} \right)_{ex} +\left(\frac{\partial f}{\partial t} \right)_{e-e}+\left(\frac{\partial f}{\partial t} \right)_{e-ph} \tag{1}$$

$$\frac{dU_e}{dt}=W_{ex}-W_{e-ph},\qquad \frac{dU_{ph}}{dt}=W_{e-ph}-G_{ph-env}(T_{ph}-T_{env}) \tag{2, 3}$$

$W_{ex}=\int \mathcal E\rho_e(\partial f/\partial t)_{ex}\,d\mathcal E$, $W_{e-ph}=-\int \mathcal E\rho_e(\partial f/\partial t)_{e-ph}\,d\mathcal E$, $\rho_e = \frac{3 n_e}{2 \mathcal E_F} \sqrt{\mathcal E/\mathcal E_F}$. Absorption is incoherent and normalized to the absorbed power:

$$\left(\frac{\partial f}{\partial t}\right)_{ex}(\mathcal E) = \frac{N_A \phi_A(\mathcal E)}{\rho_e(\mathcal E)},\qquad \phi_A(\mathcal{E};\omega) = \int_0^\infty d\omega\, [A(\mathcal{E};\omega) - A(\mathcal{E} + \hbar \omega;\omega)] \tag{A9, A8}$$

$$A(\mathcal{E};\omega) = \frac{n_A(\omega)}{N_A} \frac{D_J(\mathcal{E},\mathcal{E} - \hbar \omega) \rho_J(\mathcal{E},\mathcal{E} - \hbar \omega)}{\int D_J(\mathcal{E},\mathcal{E} - \hbar \omega) \rho_J(\mathcal{E},\mathcal{E} - \hbar \omega) d\mathcal{E}},\quad \rho_J(\mathcal{E}_{final},\mathcal{E}_{initial}) = \left[f(\mathcal E_{initial}) \rho_e(\mathcal{E}_{initial})\right] \left[(1 - f(\mathcal E_{final}) \rho_e(\mathcal E_{final})\right] \tag{A3, A4}$$

(Eq. A4 as printed; a closing parenthesis after $f(\mathcal E_{final})$ is missing.) $N_A=\langle p_{abs}\rangle_t/\hbar\omega$ with $\langle p_{abs}\rangle_t=\omega\epsilon''\langle\vec E\cdot\vec E\rangle_t$ (Eqs. A5, A6); $D_J$ is a squared transition matrix element. Electron–phonon scattering is the Bloch–Boltzmann–Peierls integral with Debye phonons and no Umklapp (Eq. A10). Electron–electron scattering:

$$\tau_{e-e}^{-1}(\mathcal E) = K \left[\left(\pi k_B T_e\right)^2 + \left(\mathcal E - \mathcal E_F\right)^2\right] \tag{A11}$$

$$\left(\frac{\partial f}{\partial t} \right)_{e-e} = - \frac{f(\mathcal E,T_e,T_{ph}) - f^T(\mathcal E,T_e)}{\tau_{e-e}(\mathcal E)} + \mathcal{F}_{e-e}(\mathcal E),\qquad \int_0^\infty \mathcal E \rho_e(\mathcal E) \left[(\Delta_\tau f)_{e-e} + \mathcal{F}_{e-e}(\mathcal E)\right] d\mathcal E = 0 \tag{A12, A13}$$

$f^T$ is a Fermi–Dirac function at the temperature $T_e$ the electrons would reach if illumination stopped without further exchange with phonons. $\mathcal F_{e-e}$ restores energy conservation; its functional form is not given. The non-thermal part is $\Delta f^{NT}\equiv f-f^T(\mathcal E,T_e)$.

**Numbers.**
- Ag sphere, 5 nm; photon energy 2.25 eV; $\epsilon_{Ag}=-8.5+1.8i$; host permittivity 4.25; $\mathcal E_F=5.1$ eV; band width 9 eV; $n_e=5.86\times10^{28}$ m⁻³; $T_{env}=297$ K; $G_{ph-env}=5\times10^{14}$ W m⁻³ K⁻¹ (Table I).
- $K=2\times10^{14}$ eV⁻² s⁻¹ is quoted for Au and Ag, but the value used was fitted to another work's data and is not stated; the fit in Fig. 6 gives ≈100 fs at 1 eV from $\mathcal E_F$, i.e. $K\approx10^{13}$ eV⁻² s⁻¹ (my reading).
- $G_{e-ph}\sim3\times10^{16}$ W m⁻³ K⁻¹ for Ag, quoted from [37], with "negligible temperature-dependence, up to about 3000" K.
- Local $|E|^2$ = 1.4×10⁶–1.4×10⁹ (V/m)²; incident intensity 7×10¹–7×10⁴ W/cm² (Figs. 3, 4).
- Temperature 297.0003–297.3 K; $T_e-T_{ph}$ = 10⁻⁵–10⁻² K (Fig. 3).
- Non-thermal power is "more than 8 orders of magnitude" below heating power at low field and 10⁻⁵ of it at the highest field (Fig. 4).

**Key results.**
- (a) Energy partition: at CW steady state nearly all absorbed power heats electrons and lattice (simulated).
- $T_e\approx T_{ph}$, both linear in $|E|^2$ (simulated).
- Just above (below) the Fermi energy the true non-thermal part consists of holes (electrons), an electron–phonon effect (simulated).
- Changing the form of $\tau_{e-e}(\mathcal E)$ alters the distribution slightly and the temperatures not at all (simulated).

**Assumptions and limitations.**
- Parabolic band, no interband transitions, constant heat capacities and $G_{e-ph}$, no field-acceleration term (stated).
- (b) Temperature-dependent heat capacity, coupling and collision rates beyond free-electron models are not treated.
- (c) No acceleration method; the full electron–electron integral is avoided as "highly time-consuming".
- Permittivity is not updated from $f$ (stated as needed at higher intensity).

**Use for this project.**
- The 10⁻⁸ fraction is a CW number at intensities three to four orders of magnitude below the plan's lowest pump intensities (0.1–1 GW/cm²). The authors call their results "very different from those known for electron dynamics under ultrafast illumination". Do not use it in the Gate-1 argument.
- Eqs. (A3)–(A9): absorption term to compare with [14]; it conserves particle number by construction.
- Eqs. (A12)–(A13): an energy-conserving relaxation-time variant for the WP2 baseline list; particle conservation with $\tau(\mathcal E)$ is not addressed (compare [39]).
- $\Delta f^{NT}$ is the deviation from an energy-matched Fermi–Dirac function, the same quantity WP2 feeds to the network.
- Data-availability statement: none found.

---
### [37] Brown et al., "Ab initio phonon coupling and optical response of hot electrons in plasmonic metals" — Phys. Rev. B 94, 075120 (2016)

**Source read:** arXiv:1602.00625v2 (28 Jul 2016) full text + TeX source, including the supplementary permittivity tables shipped with the source. Published version not compared; equation numbers are those of the arXiv PDF.
**Role in this project:** Shows how far thermal-electron properties (heat capacity, electron–phonon coupling, electron–electron rate) depart from free-electron models at high electron temperature; informs the thermal model M0 and all temperature-level baselines.

**What it does.** Computes from density-functional theory (PBEsol+U, Wannier interpolation), without fitted parameters, the electronic heat capacity $C_e(T_e)$, the electron–phonon coupling $G(T_e)$ and the dielectric function $\epsilon(\omega,T_e)$ of Al, Ag, Au and Cu. Free-electron and earlier semi-empirical estimates fail for the noble metals above about 2000–3000 K. Electrons are assumed thermal throughout.

**Model / method.**

$$C_e(T_e) = \int_{-\infty}^{\infty} \mathrm{d}\varepsilon\, g(\varepsilon)\, \varepsilon\, \frac{\partial f(\varepsilon,T_e)}{\partial T_e} \tag{3}$$

$g(\varepsilon)$ is the electronic DOS per unit volume, $f$ the Fermi function. The Sommerfeld limit is $C_e = \frac{\pi^2 n_e k_B^2}{2\varepsilon_F} T_e$ (the text prints "$n_e = 3\pi^2 k_F^3$", as printed; it should read $k_F^3/3\pi^2$). $G(T_e)$ is a golden-rule sum over all electron and phonon states with ab initio matrix elements, Umklapp included implicitly (Eq. 8). For phonon energies $\ll k_BT_e$ and $T_e\gg T_l$ it reduces to

$$G(T_e) \approx \frac{\pi k_B}{\hbar g(\varepsilon_F)} \int_{-\infty}^\infty \mathrm{d}\varepsilon\, h(\varepsilon)\, g^2(\varepsilon)\, \frac{-\partial f(\varepsilon,T_e)}{\partial\varepsilon} \tag{10}$$

with $h(\varepsilon)$ an energy-resolved coupling strength (Eq. 11). The earlier semi-empirical formula of Lin and Zhigilei (Eq. 9) replaces $h(\varepsilon)$ by a constant, $\lambda\langle(\hbar\omega)^2\rangle$. Electron–electron scattering rate, free-electron form with the prefactor fitted to DFT rates at 300 K:

$$\tau^{-1}_{ee}(\varepsilon,T_e) \approx \frac{D_e}{\hbar} [ (\varepsilon - \varepsilon_F)^2 + (\pi k_B T_e)^2 ] \tag{15}$$

for $|\varepsilon-\varepsilon_F| \ll \varepsilon_F$ and $T_e \ll \varepsilon_F/k_B$. The electron–phonon scattering rate is taken as independent of $T_e$. $\epsilon(\omega,T_e)$ is the sum of Drude, direct interband and phonon-assisted intraband terms (Eqs. 12–14); $T_e$ enters through the occupations and through the extra linewidth $(D_e/\hbar)(\pi k_B)^2(T_e^2-T_0^2)$.

**Numbers.**
- $D_e$ (Table I): 0.017 (Al), 0.021 (Ag), 0.016 (Au), 0.020 (Cu) eV⁻¹; equivalently $D_e/\hbar\approx$ 2.4–3.2×10¹³ eV⁻² s⁻¹ (my conversion).
- Highest electron temperature that avoids melting after equilibration: ≈5700 (Al), 8300 (Ag), 7500 (Au), 6700 K (Cu) (Sec. I).
- Sommerfeld $C_e$ "accurate ... up to ∼2000 K" for all four metals and throughout for Al. At 8000 K, in 10⁵ J m⁻³ K⁻¹ (Fig. 3, by eye): Ag ≈7, Au ≈10, Cu ≈20, against ≈5 from Sommerfeld.
- $G$ in 10¹⁷ W m⁻³ K⁻¹ at $T_l=298$ K (Fig. 7, by eye): Al ≈5–6, about twice the measured ≈2.5; Ag ≈0.3; Au ≈0.2 rising to ≈0.5 at 8000 K; Cu ≈1 rising to ≈2.7.
- $h(\varepsilon)$ is 2–3 times smaller in the d-bands than at $\varepsilon_F$ (Fig. 6).
- Quoted from [33] for Al at 2000 K: $G\approx2.9\times10^{7}$ W m⁻³ K⁻¹ thermalized, "1.5x larger than" $1.9\times10^{7}$ non-thermalized (as printed; inconsistent with Fig. 7, the exponent should be 17).

**Key results.**
- (b) Free-electron $C_e$ underestimates the noble metals at high $T_e$; $G(T_e)$ rises above ≈3000 K, but less than a constant-$h$ estimate predicts (calculated).
- For Al, equilibrium $G$ is twice the measured value; the authors attribute the gap to non-thermal electrons in the experiments (asserted).

**Assumptions and limitations.**
- (a) Thermal electrons only: no result on the split of energy between thermal and non-thermal electrons (non-equilibrium extension named as ongoing work).
- $T_e$ dependence of linewidths is added analytically, not recomputed (stated).
- (c) No reduced-order or data-driven acceleration.

**Use for this project.**
- For M0, compute $C_e(T_e)$ from Eq. (3) with the Kane-band DOS and a temperature-dependent chemical potential, not from the Sommerfeld formula. In ITO non-parabolicity and $k_BT_e$ approaching $E_F\approx1$ eV play the role the d-bands play here (my reading).
- Eq. (15) is the Fermi-liquid form also used in [36], [38]–[40]; prefactors differ by a factor of three or more between sources (see the [40] digest).
- The factor-two gap for Al is a warning for WP4: a coupling constant fitted with a thermal model absorbs non-thermal effects.
- Data: tables of $\epsilon(\omega,T_e)$ are in `references/sources/37_brown2016/SI/` (`*_ReEpsilon.dat`, `*_ImEpsilon.dat` for Al, Ag, Au, Cu; first column $\hbar\omega$ in eV, first row $T_e$ in K, 300–8000 K in 100 K steps). They can test a permittivity routine in the thermal limit.

---
### [38] Uehlein et al., "Implementation of the electronic non-equilibrium in the two-temperature model" — Nanomaterials 12, 1655 (2022)

**Source read:** arXiv:2112.01851v1 (3 Dec 2021) full text + TeX source. The published Nanomaterials version was not available and was not compared; equation numbers are those of the arXiv PDF.
**Role in this project:** Non-ML baseline for WP2 (extended two-temperature model, eTTM) and the reduced limit that the WP1 solver is meant to recover; its plain TTM, Eq. (1), is the closest thing here to model M0.

**What it does.** Restates the eTTM of Carpene [34] and Tsibidis (Appl. Phys. A 124, 311 (2018)): a phonon bath, a thermal (Fermi-distributed) electron bath and a non-thermal electron–hole population created by the laser. Two changes are made: a Fermi-liquid relaxation time that stays finite at the Fermi edge, and energy integrals taken over the whole band instead of $[E_\text{F}-E_\text{pt},E_\text{F}+E_\text{pt}]$. Applied to aluminium, the eTTM gives a lower and later electron-temperature maximum than the TTM and a step-shaped distribution close to a Boltzmann result directly after absorption.

**Model / method.**

Reference TTM ($u$ internal energy density, $c$ heat capacity, $g$ electron–phonon coupling parameter, $s(t)$ absorbed laser power density):

$$\frac{\mathrm{d} u_\text{el}}{\mathrm{d} t} = c_\text{el} \frac{\partial T_\text{el}}{\partial t} = - g \left(T_\text{el}-T_\text{ph}\right) + s(t),\qquad \frac{\mathrm{d} u_\text{ph}}{\mathrm{d} t} = c_\text{ph} \frac{\partial T_\text{ph}}{\partial t} = g \left(T_\text{el}-T_\text{ph}\right) \tag{1a,b}$$

eTTM: the source $s(t)$ is replaced by energy flows out of the non-thermal electrons (index el*):

$$\frac{\mathrm{d} u_\text{el}}{\mathrm{d} t} = - g \left(T_\text{el}-T_\text{ph}\right) + \frac{\partial u_\text{el*-el}}{\partial t} =: -\Pi_\text{el-ph}(t) + \Xi_\text{el*-el}(t) \tag{2a}$$

$$\frac{\mathrm{d} u_\text{ph}}{\mathrm{d} t} = g \left(T_\text{el}-T_\text{ph}\right) + \frac{\partial u_\text{el*-ph}}{\partial t} =: \Pi_\text{el-ph}(t) + \Xi_\text{el*-ph}(t) \tag{2b}$$

There is no differential equation for a non-thermal energy or particle density. The non-thermal population is a sum over past time steps of "excitation functions". Photons of energy $E_\text{pt}$ arriving in $[t',t'+\Delta t]$ create, relative to the current thermal Fermi function $f(E)=f(T_\text{el}(t'),E)$,

$$\Delta f_\text{L}(t', E) = \delta(t') \left\{f(E-E_\text{pt}) \left[1-f(E)\right]-f(E)\left[1-f(E+E_\text{pt})\right]\right\} \tag{3}$$

with amplitude $\delta(t')$ fixed by the absorbed energy ($D$ = electronic DOS, energy zero at the Fermi edge):

$$\int_{-\infty}^{\infty} \Delta f_\text{L}(t', E)\,D(E)\,E\,\mathrm{d}E = s(t') \Delta t \tag{4}$$

Each excitation function then decays exponentially at every energy:

$$\Delta f'(t, t', E) = \Delta f_\text{L}(t', E) \exp\left(-\frac{t-t'}{\tau_\text{el*-el}}-\frac{t-t'}{\tau_\text{el*-ph}}\right) \tag{5}$$

$$\tau_\text{el*-el}(E, T) = \tau_0 \frac{E_\text{F}^2}{(E-E_\text{F})^2+(\pi k_\text{B} T_\text{el})^2} \tag{6}$$

Eq. (6) is called "a simplified version from Fermi-liquid theory" (ref. [33]); $\tau_\text{el*-ph}$ is one constant. The full distribution is

$$f^\text{tot}(t, E) = f(t, E) + \Delta f(t, E) = f(t, E) + \sum_{t' = 0,\ t' \mathrel{+}= \Delta t}^{t} \Delta f'(t, t', E) \tag{7}$$

and the energy rate delivered by the electrons created at $t'$ is

$$\frac{\partial u'(t, t')}{\partial t} = -\frac{\partial}{\partial t}\int_{-\infty}^{\infty} \Delta f'(t, t', E)\,D(E)\,E\,\mathrm{d}E = - \int_{-\infty}^{\infty} \left(-\frac{1}{\tau_\text{el*-el}}-\frac{1}{\tau_\text{el*-ph}}\right) \Delta f_\text{L}(t', E) \exp\left(-\frac{t-t'}{\tau_\text{el*-el}}-\frac{t-t'}{\tau_\text{el*-ph}}\right) D(E)\,E\,\mathrm{d}E =: \Xi'_\text{el*-el}(t, t') + \Xi'_\text{el*-ph}(t, t') \tag{8}$$

The $1/\tau_\text{el*-el}$ part goes to the thermal electrons and the $1/\tau_\text{el*-ph}$ part to the phonons. Summing over history,

$$\Xi_\text{el*-el/ph}(t) := \frac{\partial u_\text{el*-el/ph}(t)}{\partial t} = \sum_{t' = 0,\ t' \mathrel{+}= \Delta t}^{t} \Xi'_\text{el*-el/ph}(t, t') \tag{9}$$

How the relaxation times are obtained: $\tau_0 = \frac{128}{\sqrt{3}\pi^2\omega_\text{P}}$ (Eq. A6, from the plasma frequency); $\tau_\text{el*-ph}$ "according to" Carpene [34] and Groeneveld et al. (1995) with textbook data (formula not printed here; [39] gives $\tau_\text{mfp}\hbar\omega/k_\text{B}\theta_\text{D}$). Electron heat capacity from the DOS: $c_\text{el} = \int \frac{\mathrm{d}f(T_\text{el}, E)}{\mathrm{d}T}\,D(E)\,E\,\mathrm{d}E$ (Eq. A4), DOS from DFT (Lin et al. 2008); $c_\text{ph}$ constant (Dulong–Petit). No transport: homogeneously heated thin film. Laser: $s(t)=(1-R)\alpha I(t)$ with

$$I(t) = \sqrt{\frac{4\ln(2)}{\pi}} \frac{\Phi}{\tau_\text{p}}\exp\left[-4\ln(2)\left(\frac{t}{\tau_\text{p}}\right)\right] \tag{A1}$$

(as printed; looks inconsistent because a Gaussian needs $(t/\tau_\text{p})^2$ in the exponent.)

**Numbers.** (aluminium; Appendix A unless noted)
- Pulse: fluence $\Phi=15$ J m⁻² (incident), FWHM $\tau_\text{p}=50$ fs, photon energy $E_\text{pt}=1.55$ eV, absorption coefficient $\alpha=1.3123\times10^{8}$ m⁻¹, reflectivity $R=0.8682$.
- $c_\text{ph}=2.5\times10^{6}$ J K⁻¹ m⁻³ ($\rho_V=2.7\times10^{6}$ g m⁻³, $M=26.98$ g mol⁻¹).
- $\tau_0=0.329$ fs with $\hbar\omega_\text{P}=14.98$ eV; $\tau_\text{el*-ph}=299.22$ fs.
- $g=2.6148\times10^{17}$ J K⁻¹ s⁻¹ m⁻³, constant, Lin et al. value at 1000 K; fcc lattice constant $4.05\times10^{-10}$ m.
- $E_\text{F}$ used in Eq. (6), time step $\Delta t$, energy grid: not stated.
- Peak electron temperature ≈2.2×10³ K (TTM) vs ≈1.95×10³ K (eTTM); common temperature ≈350 K after ≈1.5 ps (Fig. 2, by eye).
- Peak non-thermal energy density ≈2.1×10⁸ J m⁻³, "more than two times" the thermal-electron peak; thermalization complete within 200 fs after the pulse (Fig. 3, Sec. 3.1).
- Conservation error (Eq. 10, normalized to initial value): per-mille; energy peak ≈+1.6×10⁻³ (this model) vs ≈+1.8×10⁻³ (Tsibidis); particles ≈−2.0×10⁻³ vs ≈−2.7×10⁻³ (Fig. 6, by eye).
- Density of electrons 1–1.5 eV above $E_\text{F}$: eTTM peak ≈2×10²⁶ m⁻³, about ten times the TTM peak, decaying about 50 fs after the pulse (Fig. 8).
- δ-pulse comparison (Sec. 3.3): initial 300 K; TTM temperature after excitation 36 418 K (as printed; the plotted TTM curve in Fig. 7a falls to ≈0.05 at $E_\text{F}+3$ eV, which looks like $k_\text{B}T\approx1$ eV — my reading; treat as unverified).
- Consistency check (mine): $(1-R)\alpha\Phi=2.6\times10^{8}$ J m⁻³ with the printed values, about twice the ≈1.3×10⁸ J m⁻³ final energy gain in Fig. 3 and twice the ≈2.4×10²¹ W m⁻³ peak source in Fig. 4. Do not rely on the printed fluence to reproduce the figures.

**Key results.**
- eTTM heats the thermal electrons indirectly, so the maximum of $T_\text{el}$ is lower and later than in the TTM; phonon temperature and the later relaxation are nearly unchanged (simulated, Figs. 2–4).
- $\Xi_\text{el*-ph}\ll\Xi_\text{el*-el}$: non-thermal electrons give their energy mainly to thermal electrons (simulated, Fig. 4).
- Whole-band integration removes the kink at $E_\text{F}\pm E_\text{pt}$ and reduces the conservation error, but energy and particle number are still not conserved exactly (simulated, Figs. 5–6).
- Validation against Boltzmann is limited to one snapshot: the distribution directly after a δ-like excitation of the same energy density, compared with the electron–ion–photon collision term of [33]. Both show a step of width $E_\text{pt}$; the Boltzmann step is lower, has a second step from two-photon absorption and carries DOS features that the eTTM lacks (simulated, Fig. 7).
- No comparison of time evolution with Boltzmann; "Further comparisons ... are ongoing" (stated).

**Assumptions and limitations.**
- Step shape of Eq. (3) ignores the DOS; DOS enters only the amplitude (stated).
- One-photon absorption only (stated).
- Non-thermal electrons decay in place; no secondary-electron cascade (my reading of Eq. 5).
- Constant $g$ and $\tau_\text{el*-ph}$; no transport (stated).
- Cost grows with the history sum in Eqs. (7), (9); no timings given (my reading).
- Eq. (6) gives picosecond lifetimes at the Fermi edge (stated).

**Use for this project.**
- Discrepancy with the task description: this paper has no extra rate equation for a non-thermal energy or particle density, and its relaxation times are not derived from a Boltzmann solution. The three-equation form with a non-thermal energy density belongs to the Sun/Carpene/Della Valle line ([34], [35]).
- Terminology clash with the plan: M0 is "thermal electrons only: extended two-temperature model", citing [9, 24, 38]. In [38], "extended" means adding non-thermal electrons; its thermal-only model is the plain TTM, Eq. (1), with $c_\text{el}$ from Eq. (A4).
- "Recover the extended two-temperature limit [38]" needs a definition. A workable one: Boltzmann with $Q_{ee}$ replaced by per-energy exponential decay, compared through $f^\text{tot}$ of Eq. (7).
- For ITO, prefer the DOS-weighted excitation of [39]: with $\hbar\omega\approx E_\text{F}$ the hole step reaches the band bottom. [15] states that the fixed-profile ansatz of Eq. (3) "is not valid for ITO" because the chemical potential shifts.
- Expect per-mille conservation errors, far above the plan's 10⁻⁶ target.
- Reusable diagnostic: density in an energy window (Fig. 8).
- Data-availability statement: none in the arXiv version.

---
### [39] Uehlein et al., "Capturing non-equilibrium electron dynamics in metals accurately and efficiently" — J. Appl. Phys. 138, 063103 (2025)

**Source read:** arXiv:2503.09479v1 (12 Mar 2025) main text + Supplementary Information, PDF and TeX source. The published J. Appl. Phys. version was not available and was not compared. The SI carries a different title ("Non-equilibrium electron dynamics in metals modelled accurately and efficiently"). Authors: Uehlein, Snowden, Seibel, Held, Weber, Maurer, Rethfeld.
**Role in this project:** The "athermal electron model" (AthEM) is the second non-ML baseline of WP1/WP2; its SI also shows why an energy-dependent relaxation time applied to a single distribution is not a valid baseline.

**What it does.** Splits the electrons into a Fermi-distributed thermal bath and an explicitly propagated athermal distribution $f_\text{el*}(E,t)$ (positive for electrons, negative for holes), coupled to a phonon bath. Laser excitation enters through Fermi's golden rule with the real DOS; thermalization uses a Fermi-liquid relaxation time; whatever energy and particles the athermal part loses are credited to the thermal bath, so both are conserved by construction. For Al and Au the model reproduces Boltzmann results with a full collision integral (random-k, as in [40]) for excitation, electron temperature and electron–phonon energy flow, but overestimates energy-resolved thermalization times.

**Model / method.**

Transported quantities: $f_\text{el*}(E,t)$ on an energy grid, plus three scalars — thermal-electron temperature $T_\text{el}$, phonon temperature $T_\text{ph}$ and thermal-electron particle density $n_\text{el}$ (equivalently chemical potential $\mu$). The total distribution is $f=f_\text{el}+f_\text{el*}$.

$$\frac{\partial f_\text{el*}}{\partial t} = \left.\frac{\partial f_\text{el*}}{\partial t}\right|_\text{laser} + \left.\frac{\partial f_\text{el*}}{\partial t}\right|_\text{el*-el} + \left.\frac{\partial f_\text{el*}}{\partial t}\right|_\text{el*-ph} \tag{1}$$

$$\frac{\mathrm{d}u_\text{el}}{\mathrm{d}t} = - g_\text{ep}\,(T_\text{el}-T_\text{ph}) + \left.\frac{\partial u_\text{el}}{\partial t}\right|_\text{el-el*},\qquad \frac{\mathrm{d}u_\text{ph}}{\mathrm{d}t} = g_\text{ep}\,(T_\text{el}-T_\text{ph}) + \left.\frac{\partial u_\text{ph}}{\partial t}\right|_\text{ph-el*} \tag{2a,b}$$

$$\left.\frac{\partial u_s}{\partial t}\right|_{s-\text{el*}} = - \left.\frac{\partial u_\text{el*}}{\partial t}\right|_{\text{el*}-s} = - \int \left.\frac{\partial f_\text{el*}}{\partial t}\right|_{\text{el*}-s}\,D(E)\,E\,\mathrm{d}E,\qquad s\in\{\text{el},\text{ph}\} \tag{3}$$

$$\frac{\mathrm{d}n_\text{el}}{\mathrm{d}t} = - \frac{\mathrm{d}n_\text{el*}}{\mathrm{d}t} = - \int \left(\left.\frac{\partial f_\text{el*}}{\partial t}\right|_{\text{el*-el}} + \left.\frac{\partial f_\text{el*}}{\partial t}\right|_{\text{el*-ph}}\right)\,D(E)\,\mathrm{d}E \tag{5}$$

$u$ = internal energy density, $D(E)$ = electronic DOS, $g_\text{ep}$ = electron–phonon coupling parameter. "Equations (1), (2), (3) and (5)" define the model.

Laser term, with $E_\pm=E\pm\hbar\omega$ and $V$ the unit-cell volume:

$$\left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_\text{laser} = \left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_{h^+} + \left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_{e^-} \tag{6}$$

$$\left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_{h^+} = \frac{2\pi V}{\hbar}\,D(E_+)\,|M_{E,E_+}|^2\,f(E)[1-f(E_+)] \tag{7}$$

$$\left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_{e^-} = \frac{2\pi V}{\hbar}\,D(E_-)\,|M_{E,E_-}|^2\,f(E_-)[1-f(E)] \tag{8}$$

(as printed; both terms are positive, although Eq. (7) is described as creating a hole at $E$ and holes count negatively — the hole term needs a minus sign for Eq. (4), $\int \partial_t f_\text{el*}|_\text{laser} D\,\mathrm dE=0$, to hold.) The two matrix elements are taken independent of energy and are re-solved at every time step from two conditions: equal numbers of holes and electrons (Eq. 9) and

$$s(t) = \int \left.\frac{\partial f_\text{el*}(E,t)}{\partial t}\right|_\text{laser}\,D(E)\,E\,\mathrm{d}E,\qquad s(t)=\alpha I(t) \tag{10}$$

with constant absorption coefficient $\alpha$ (SI Eqs. 12–17 give the 2×2 linear system for $\delta_{h^+},\delta_{e^-}=\frac{2\pi V}{\hbar}|M|^2$).

Thermalization (relaxation-time closure):

$$\left.\frac{\partial f_\text{el*}(E)}{\partial t}\right|_\text{el*-el} = -\frac{f_\text{el*}(E)}{\tau_\text{el*-el}(E)} + \frac{f_\text{el}(E)-f^\text{rlx}(E)}{\tau_\text{el*-el}(E)} \tag{11}$$

$$\tau_\text{el*-el}(E) = \tau_0\frac{\mu^2}{(E-\mu)^2+(\pi k_BT_{\text{el}})^2},\qquad \tau_0 = \frac{128}{\sqrt{3}\pi^2\omega_\text{p}} \tag{12, 13}$$

$f^\text{rlx}(E;T^\text{rlx},\mu^\text{rlx})$ is the Fermi function with the same energy and particle density as the total $f$, found by two-dimensional root finding at each step (SI S5). The second term of Eq. (11) is "a very small correction" at weak excitation. Athermal–phonon coupling:

$$\left.\frac{\partial f_\text{el*}}{\partial t}\right|_\text{el*-ph} = - \frac{f_\text{el*}}{\tau_\text{el*-ph}},\qquad \tau_\text{el*-ph} = \tau_\text{mfp}\frac{\hbar\omega}{k_\text{B}\theta_\text{D}} \tag{14, 15}$$

$\tau_\text{mfp}$ is called "the mean free path of the electrons" but is a time; $\theta_\text{D}$ is the Debye temperature. Eq. (15) is "a rough estimate yielding an upper limit for the energy relaxation rate". Thermal variables are advanced via (SI Eq. 26)

$$\frac{\mathrm{d}T_\text{el}}{\mathrm{d}t} = \frac{1}{c_T p_\mu - p_T c_\mu} \left(p_\mu \frac{\mathrm{d}u_\text{el}}{\mathrm{d}t} - c_\mu \frac{\mathrm{d}n_\text{el}}{\mathrm{d}t} \right),\quad c_x = \int \frac{\partial f_\text{el}}{\partial x}D\,E\,\mathrm{d}E,\quad p_x = \int \frac{\partial f_\text{el}}{\partial x}D\,\mathrm{d}E,\quad x\in\{T_\text{el},\mu\}$$

with $\mu$ from one-dimensional root finding and constant $c_\text{ph}$. Athermal–athermal collisions are set to zero (SI S2). The reference Boltzmann model uses full collision integrals for all terms (screened Coulomb, plane waves; refs. [32], [33]), e–e in the random-k approximation, and includes simultaneous two-photon absorption "as described by Rethfeld et al." [32].

**Numbers.**

| Parameter (Table I) | Al | Au |
| --- | --- | --- |
| Fermi energy (eV) | 11.15 | 9.94 |
| lattice constant (Å) | 4.04 | 4.07 |
| $\theta_\text{D}$ (K) | 428 | 165 |
| $c_\text{ph}$ (10⁶ J K⁻¹ m⁻³) | 2.5 | 2.5 |
| $\hbar\omega_\text{p}$ (eV) | 14.98 | 9.03 |
| $\tau_\text{mfp}$ (fs) | 7.12 | 10.9 |

- Laser sets (SI Table S1): A — Al, 25 fs FWHM, 3.1 eV, absorbed energy density 8.1×10⁸ J m⁻³; B — Au, 25 fs, 3.1 eV, 1.7×10⁹ J m⁻³; C — Au, 50 fs, 1.55 eV, 1.0×10⁷ J m⁻³; D — Au, 25 fs, 1.55 eV, 3.3×10⁹ J m⁻³.
- DOS: DFT (FHI-aims, PBE, 12×12×12 k-grid), −20 to 20 eV. $g_\text{ep}(T_\text{el},T_\text{ph})$ from a separate Boltzmann calculation (their ref. Held 2025); values not given.
- Derived (my calculation from Eqs. 13, 15): $\tau_0\approx0.33$ fs (Al), 0.55 fs (Au); $\tau_\text{el*-ph}\approx0.60$ ps (Al, 3.1 eV), 2.4 ps (Au, 3.1 eV), 1.2 ps (Au, 1.55 eV).
- Energy-grid size, time step, integrator, run times, speed-up factor: not stated.
- Al thermalization time of the spectral density (window 0.2 eV) for centre energies 2–3 eV above $E_\text{F}$: Boltzmann ≈4.5→1 fs, AthEM ≈6→4.5 fs, constant-τ RTA ≈5 fs flat (Fig. 5 inset, by eye).
- Single-distribution RTA with τ(E) in Al: loses up to ≈1.2×10⁻³ electrons per atom within 20 fs; quasi-temperature falls from ≈4 kK to ≈2 kK (Fig. 4).
- eTTM [38] particle error: +1.2×10⁻³ (Al), −11.8×10⁻³ (Au) electrons per atom (SI S3).
- Au, set B: electrons above $E_\text{F}$ ≈6.7×10²⁷ m⁻³ total, ≈2.9×10²⁷ m⁻³ primary; peak quasi-temperature ≈6 kK (Figs. 6, 9, by eye).
- Au, set C, spectral density at 1.2 eV, convolved with 40 fs Gaussian: AthEM peak more than twice Boltzmann; both about 10⁶ times the TTM (Fig. 12).

**Key results.**
- Primary-carrier distributions of AthEM and Boltzmann coincide for Al and Au, including d-band features; the eTTM step does not (simulated, Fig. 3).
- AthEM conserves particles and energy with τ(E); a single-distribution RTA with the same τ(E) does not (simulated, Fig. 4; derived, SI S6).
- AthEM overestimates thermalization times at every energy; Boltzmann shows slower thermalization at the Fermi edge than AthEM (simulated, Figs. 5, 8).
- Quasi-temperatures of AthEM, Boltzmann and TTM agree; the thermal bath lags by "a few femtoseconds" (simulated, Fig. 9).
- Electron–phonon energy flow is dominated by the thermal channel and matches Boltzmann (simulated, Figs. 10–11).
- Normalized photoemission transient of gold (Kühne et al. 2022) is matched by AthEM and Boltzmann, not by TTM (measured vs simulated).
- "Significantly reduced computational effort" is asserted without numbers.

**Assumptions and limitations.**
- Single band, energy-independent optical matrix element, constant $\alpha$ (stated).
- One-photon steps only; no simultaneous two-photon absorption (stated).
- No athermal–athermal collisions; no in-scattering of secondaries beyond the $f^\text{rlx}$ term (stated).
- Fermi-liquid τ(E) is a single-particle lifetime; phonons thermal (stated).
- No transport (stated as future work).
- Agreement is judged by eye; no error norms (my reading).

**Use for this project.**
- Baseline to implement: Eqs. (1)–(15) plus SI Eq. (26); state size N + 3; cost per step O(N) plus a root find (my reading).
- The DOS-weighted excitation, Eqs. (7)–(8), carries over to a Kane band; with $\hbar\omega\approx E_\text{F}$ in ITO, $D(E_-)=0$ below the band bottom truncates the hole step automatically (my reading).
- Because $n_\text{el}$ and $\mu$ are tracked, AthEM can follow the large chemical-potential shift of ITO, which the fixed-profile eTTM cannot ([15]).
- $\tau_0$ of Eq. (13) is a free-electron estimate; [40] needed 0.1 fs instead of 0.284 fs, and [15] reports a non-parabolicity factor for ITO. Fit it.
- Neglecting athermal–athermal collisions is doubtful at TW/cm² (my reading); test against B0.
- The plan's baseline "RTA with energy-dependent lifetime [40]" should be run in this two-subsystem form, or its conservation error reported.
- Speed-up numbers requested for this paper do not exist; measure them at Gate 1.
- Data availability: "available from the corresponding authors upon reasonable request." No code link.

---
### [40] Roden et al., "Thermalization of optically excited Fermi systems: electron-electron collisions in solid metals" — Phys. Plasmas 33, 053906 (2026)

**Source read:** arXiv:2601.11371v1 (16 Jan 2026) full text + TeX source. The published Phys. Plasmas version was not available and was not compared; equation numbers below are those of the arXiv PDF.
**Role in this project:** Defines fidelity level B0 of WP1 (random-k electron–electron collision integral with one constant matrix element) and supplies the constant-τ and Fermi-liquid relaxation-time baselines of WP2.

**What it does.** Derives the electron–electron Boltzmann collision integral for an isotropic, spin-degenerate, single-band metal in the random-**k** approximation, with a screened-Coulomb matrix element averaged to a single number. It then follows the thermalization of one step-shaped laser-excited distribution of a degenerate free-electron gas (no electron–phonon term, no laser during relaxation) and compares the full collision integral (FCI) with three relaxation-time approaches (RTA): constant τ, Fermi-liquid τ(E) with calculated prefactor, and Fermi-liquid τ(E) with fitted prefactor. Fermi-liquid τ(E) with fitted prefactor reproduces the FCI distribution best, but only the FCI gives the non-monotonic energy-resolved occupation time near the Fermi energy.

**Model / method.**

Momentum-resolved starting point (Fermi's golden rule; factor 1/2 avoids double counting of $\boldsymbol{p'},\boldsymbol{k'}$):

$$\left.\frac{\partial f_{\boldsymbol{p}}}{\partial t}\right|_{\text{el-el}} = \frac{2\pi}{\hbar}\,\frac{1}{2}\sum_{\boldsymbol{p'}\boldsymbol{k}\boldsymbol{k'}} \delta\left(E_{\boldsymbol{p}} + \epsilon_{\boldsymbol{k}} - E_{\boldsymbol{p'}} - \epsilon_{\boldsymbol{k'}}\right)\left[ \mathcal{F}_{\text{in}} \left|M\left(\boldsymbol{p'},\boldsymbol{k'},\boldsymbol{p},\boldsymbol{k}\right)\right|^2 - \mathcal{F}_{\text{out}} \left|M\left(\boldsymbol{p},\boldsymbol{k},\boldsymbol{p'},\boldsymbol{k'}\right)\right|^2\right] \tag{3}$$

$$\mathcal{F}_{\text{in}} = f_{\boldsymbol{p'}}f_{\boldsymbol{k'}}(1-f_{\boldsymbol{p}})(1-f_{\boldsymbol{k}}),\qquad \mathcal{F}_{\text{out}} = f_{\boldsymbol{p}}f_{\boldsymbol{k}}(1-f_{\boldsymbol{p'}})(1-f_{\boldsymbol{k'}}) \tag{4a,b}$$

The effective $|M|^2$ is the sum of four processes (same-spin direct and exchange, opposite-spin without and with spin flip). The exchange interference term is neglected (Eq. 7), so the four are treated as distinguishable and equal, each given by a Yukawa potential $V_{\text{ee}} = \frac{e^2}{4\pi\varepsilon_0 r}\mathrm{e}^{-\kappa r}$ (Eq. 9):

$$\left|M_{\text{ee}}(\Delta k)\right|^2 = \left( \frac{e^2}{\varepsilon_0\Omega_0} \frac{1}{\Delta k^2 + \kappa^2} \right)^2,\qquad \Delta k^2 = (\boldsymbol{k} - \boldsymbol{k'})^2 = (\boldsymbol{p} - \boldsymbol{p'})^2 \tag{10, 11}$$

$\Omega_0$ is the unit-cell volume, $\kappa$ the screening parameter. Random-**k** approximation: angular dependence is dropped (assumed smeared out by elastic scattering), momentum conservation is "implicitly fulfilled", and every momentum sum becomes an energy integral,

$$\sum_{\boldsymbol{k}} f_{\boldsymbol{k}} \rightarrow \int \mathrm{d}\varepsilon\, D_{\varepsilon}f_{\varepsilon},\qquad D_{\varepsilon} = \frac{\Omega_0}{2\pi^2} k^2 \frac{\mathrm{d}k}{\mathrm{d}\varepsilon} \tag{12, 13}$$

with $D_\varepsilon$ the density of states (per unit cell, by Eq. 13) and $f_\varepsilon$ the spherical average of $f_{\boldsymbol k}$. The authors call the transformation strictly valid only where Eq. (13) holds, "for example" free electrons. The delta function removes the $\varepsilon'$ integral ($\varepsilon'=\varepsilon+\Delta E$, $\Delta E = E-E'$). Final collision integral:

$$\left.\frac{\partial f_{E}}{\partial t}\right|_{\text{el-el}}= \frac{2\pi}{\hbar} \, g \int \mathrm{d}E' \int \mathrm{d}\varepsilon\; D_{E'}D_{\varepsilon}D_{\varepsilon + \Delta E}\,\mathcal{F}_E\,\left|M_{\text{ee}}\right|^2 \tag{14}$$

$$\mathcal{F}_E = f_{E'}f_{\varepsilon+\Delta E}(1-f_{E})(1-f_{\varepsilon}) - f_{E}f_{\varepsilon}(1-f_{E'})(1-f_{\varepsilon+\Delta E}) \tag{15}$$

$g=2$ is the spin-degeneracy factor, "entering through the two last terms" of Eq. (8), the opposite-spin processes (numerically $\tfrac12\times4$ equal terms, my reading). The matrix element is taken independent of energy and replaced by its average over all $\boldsymbol k,\boldsymbol k'$ in the sphere of radius $k_r = \sqrt[3]{6\pi^2/\Omega_0}$:

$$\left\langle M_{\text{ee}}^2\right\rangle = \frac{\Omega_0^2}{(2 \pi)^6} \int \mathrm{d}^3 k\int \mathrm{d}^3 k'\, \left|M_{\text{ee}}(\Delta k)\right|^2 = \frac{e^4}{16 \pi^4 \varepsilon_0^2} \left[ \frac{k_r^3}{3\kappa}\left( \frac{\pi}{2} + \arctan{\frac{4k_r^2 - \kappa^2}{4k_r\kappa}}\right) + \left(\frac{k_r^2}{2} + \frac{\kappa^2}{12} \right) \ln\left(\frac{\kappa^2}{4k_r^2 + \kappa^2}\right) + \frac{k_r^2}{3} \right] \tag{16}$$

Screening is static and recomputed from the distribution at every time step:

$$\kappa^2(t) = \frac{e^2}{\varepsilon_0} \int \mathrm{d}E\, D_E \frac{\mathrm{d}f_E(t)}{\mathrm{d}E} \tag{17}$$

(as printed; looks inconsistent because $\mathrm{d}f/\mathrm{d}E<0$ makes the right side negative, and with $D_E$ per unit cell a factor $1/\Omega_0$ is needed for units of m$^{-2}$ — use $-\mathrm{d}f/\mathrm{d}E$ and a DOS per volume.)

Relaxation-time comparison models:

$$\left.\frac{\partial f_{E}(t)}{\partial t}\right|_{\text{el-el}} = - \frac{f_E(t)-f_E^{\text{Fermi}}}{\tau} \tag{18}$$

$$\tau_\text{E} = \tau_0 \frac{E_\text{F}^2}{\left(E-E_\text{F}\right)^2 + \left(\pi k_\text{B}T_\text{e}\right)^2},\qquad \tau_{0,\,\text{calc}} = \frac{128}{\sqrt{3}\pi^2\omega_\text{p}} \tag{19, 20}$$

$f_E^{\text{Fermi}}$ is the Fermi distribution with the same energy and particle density as $f_E$; it defines the quasi-temperature $T_\text{e}$. Thermalization is measured by $\Delta_{\text{MAD}}(t)=\int \mathrm{d}E\,\left|f_E(t)-f_E^\text{Fermi}\right|$ (Eq. 21), fitted with $\exp(-t/\tau_\text{const})$ (Eq. 22).

Numerical evaluation: not described. The paper gives no energy grid, quadrature, time integrator, time step, cost or scaling, and does not mention convolutions or FFT.

**Numbers.**
- System: free-electron gas, $D(E)\propto\sqrt E$; $E_\text{F}=11.18$ eV, $m=1.65\,m_\text{e}$, $\hbar\omega_\text{p}=17.34$ eV (Sec. III B, from B. Müller's PhD thesis). Material not named. Electron density, $\Omega_0$, $k_r$, $\kappa$ and $\langle M_\text{ee}^2\rangle$ values: not stated. (My calculation: these $E_\text{F}$ and $m$ imply $n\approx3.6\times10^{29}$ m$^{-3}$, which reproduces 17.34 eV, and a $T=0$ Thomas–Fermi $\kappa\approx3.0\times10^{10}$ m$^{-1}$.)
- Initial state: produced beforehand by a full Boltzmann excitation term as in [32], [33]; wavelength 400 nm (3.1 eV); absorbed energy equivalent to a quasi-temperature of about 1800 K; reference cold distribution 300 K (Sec. III A, Fig. 2). Pulse duration and fluence: not stated.
- MAD decay time: FCI 15.2 fs; RTA constant τ: $\tau_\text{const}=15.2$ fs (fitted to FCI); RTA with $\tau_{0,\text{calc}}=0.284$ fs: 39.8 fs; RTA with $\tau_{0,\text{fit}}=0.1$ fs: 15.2 fs by construction (Fig. 3, Sec. III B).
- FCI distribution "closely resembles" equilibrium within about 30 fs; cold Fermi edge still visible at 5 and 10 fs; replicas at ±3.1 eV gone by 5 fs (Fig. 4a).
- Energy-resolved occupation time, FCI (Fig. 5, read by eye): local maximum ≈22 fs at $E_\text{F}$; discontinuity and global maxima ≈45 fs at about ±0.7 eV (text: "towards 1 eV"); a few fs beyond ±2 eV; small kinks at ±3.1 eV. Eq. (19) at $E_\text{F}$, 1800 K gives ≈53 fs ($\tau_{0,\text{fit}}$) and ≈150 fs ($\tau_{0,\text{calc}}$) (my calculation).

**Key results.**
- Full e–e collision integral reduces to a double energy integral over three DOS factors and one Pauli functional, Eq. (14) (derived).
- Constant-τ RTA matches the integrated MAD decay but keeps the step structure at all times and overestimates the hot tail (simulated, Fig. 4b).
- Fermi-liquid τ(E) with the textbook prefactor thermalizes 2.6 times too slowly (39.8 vs 15.2 fs); the prefactor must be reduced from 0.284 to 0.1 fs (simulated).
- With fitted prefactor, τ(E) relaxes faster than FCI in the first 15 fs and slower later; high-energy tails agree, the Fermi-edge region does not (simulated, Figs. 3–5).
- Occupation time is symmetric about $E_\text{F}$ for the free-electron DOS; its peaks sit where the non-thermal distribution crosses the hot Fermi distribution and in- and out-scattering balance (simulated).
- "Formally, only a constant relaxation time ensures the conservation of energy and particle density" (asserted, citing [39]).
- At intermediate times the FCI distribution looks hotter than its final temperature (simulated).

**Assumptions and limitations.**
- Random-**k**: no momentum conservation, no angular information; exchange interference neglected; one energy-independent matrix element (stated).
- Static screening only; one band; results only for a free-electron DOS (stated).
- One excitation condition (400 nm, ≈1800 K), strongly degenerate; no e–ph coupling or laser during relaxation (stated).
- No conservation diagnostics and no numerical details are reported (my reading).
- The Φ-function is printed as $\Phi=\ln(f^{-1}-1)$ but the Fig. 2 axis corresponds to $-\ln(f^{-1}-1)$, i.e. logit $f$ (my reading).
- $\tau_{0,\text{fit}}$ is fitted to this single case; its transferability is untested (my reading).

**Use for this project.**
- Adopt Eqs. (14)–(17) as B0. Eq. (13) holds for any isotropic single band, so a Kane-band $D(E)$ is admissible (my reading; the authors only name free electrons).
- FFT hypothesis: the paper is silent — it neither supports nor contradicts it. The structure of Eq. (14) does support it (my derivation): with $a=D f$, $b=D(1-f)$ and $C(u)=\int\mathrm d\varepsilon\,b(\varepsilon)\,a(\varepsilon+u)$,
$$\left.\frac{\partial f_E}{\partial t}\right|_{\text{el-el}}=\frac{2\pi}{\hbar}\,g\,\langle M_\text{ee}^2\rangle\Big[(1-f_E)\,(a*C)(E)-f_E\,(b*\tilde C)(E)\Big],\quad \tilde C(u)=C(-u),$$
one cross-correlation and two convolutions, O(N log N) on a uniform grid with zero padding, versus O(N³) by direct summation. I checked this numerically (N = 96, $\sqrt E$ DOS, step-shaped $f$): agreement with the direct triple sum to 10⁻¹⁴, particle and energy moments zero to round-off, Fermi–Dirac an exact fixed point. A matrix element depending only on $\Delta E$ would keep the convolution form; one depending on $\Delta k$ (B1) does not.
- The convolution needs a uniform energy grid, which conflicts with the clustered grid planned for WP2.
- $\kappa$ and $\langle M_\text{ee}^2\rangle$ depend on $f$ only through one scalar per step, so screening feedback is free in B0.
- Validation target: 15.2 fs MAD decay and Fig. 4a/5, but $\Omega_0$ is missing, so reproduction needs the authors' parameters or the thesis.
- Baseline pitfall: the plan lists "RTA with energy-dependent lifetime [40]" as a baseline; applied to a single distribution it violates particle and energy conservation ([39] quantifies this). Report it with that caveat or use the two-subsystem form of [39].
- Prefactors differ widely between sources: $1/(\tau_0E_\text{F}^2)$ is 2.8×10¹³ (calc) and 8.0×10¹³ eV⁻² s⁻¹ (fit) here, against $D_e/\hbar\approx$ 2.4–3.2×10¹³ eV⁻² s⁻¹ in [37] (my calculation). Treat the ITO prefactor as a fit parameter.
- Open: accuracy of B0 against a momentum-resolved integral is not tested here; that is the Gate-1 measurement.

---
### [41] Maliyov et al., "Dynamic mode decomposition of nonequilibrium electron-phonon dynamics: accelerating the first-principles real-time Boltzmann equation" — npj Comput. Mater. 10, 123 (2024)

**Source read:** arXiv:2311.07520v1 (13 Nov 2023) full text + TeX source. Published npj version not compared; equation numbers are those of the arXiv PDF.
**Role in this project:** The one data-driven, linear-algebra acceleration of a Boltzmann solver in this group: a non-neural comparison point for WP2 and input to the Gate-1 question of how much speed-up cheap methods already give.

**What it does.** Applies dynamic mode decomposition (DMD) to first-principles real-time Boltzmann simulations with electron–phonon collisions. The equation is stepped explicitly over a short sampling window; a low-rank linear propagator fitted to those snapshots extrapolates the populations to steady state or equilibrium. About 10 % of the trajectory suffices in the cases shown, for a cost reduction of 10–17 times.

**Model / method.**

$$\frac{\partial f_{\mathbf{k}}(t)}{\partial t} = \frac{\mathbf{F}}{\hbar}\cdot \nabla _\mathbf{k}f_{\mathbf{k}}(t)+\mathcal{I}[f_{\mathbf{k}}(t)] \tag{1}$$

$f_{\mathbf k}$ are band- and momentum-resolved populations, $\mathbf F$ an external field, $\mathcal I$ the electron–phonon collision integral with ab initio matrix elements (Eq. 14). Electron–electron collisions and an optical drive are absent; the materials are GaAs, GaN and graphene. DMD assumes a time-independent linear map between successive snapshots:

$$f_{\mathbf{k}}(t + \Delta t) = \mathbf{A} f_{\mathbf{k}}(t) \tag{2}$$

Snapshots $t_1\dots t_{M-1}$ form the columns of $\mathbf X_1$ and $t_2\dots t_M$ those of $\mathbf X_2$. A truncated singular value decomposition $\mathbf{X}_1 \approx \tilde{\mathbf{U}} \tilde{\mathbf{\Sigma}} \tilde{\mathbf{V}}^\dagger$ keeps $r$ singular values, and

$$\tilde{\mathbf{A}} = \tilde{\mathbf{U}}^\dagger \mathbf{A} \tilde{\mathbf{U}} = \tilde{\mathbf{U}}^\dagger \mathbf{X}_2 \tilde{\mathbf{V}} \tilde{\mathbf{\Sigma}}^{-1},\qquad \tilde{\mathbf{A}} \mathbf{W} = \mathbf{W} \mathbf{\Lambda},\qquad \mathbf{\Phi} = \mathbf{X}_2 \tilde{\mathbf{V}} \tilde{\mathbf{\Sigma}}^{-1} \mathbf{W} \tag{7, 8, 9}$$

$$\omega^{\mathrm{DMD}}_l = -i \frac{\ln{\lambda_l}}{\Delta t},\qquad \mathbf{b} = \mathbf{\Phi}^{+} f_{\mathbf{k}}(0) \tag{10, 12}$$

$$f_{\mathbf{k}}(t>t_M) \approx \sum_{l=1}^r \, b_l\,\phi^l_{\mathbf{k}} \, e^{ i \omega_l^{\mathrm{DMD}}\, t} \tag{4}$$

$\lambda_l$ are the eigenvalues of $\tilde{\mathbf A}$, $\phi^l$ the columns of $\mathbf\Phi$, $^+$ the pseudo-inverse. For relaxation to a known equilibrium, DMD is applied to $\delta f_{\mathbf{k}}(t) = f_{\mathbf{k}}(t) - f^{\mathrm{FD}}_{\mathbf{k}}$. For velocity–field curves the field is raised in steps and DMD is re-learned at each value.

**Numbers.**
- State size $N\approx10^5$–$10^6$ k-points (grids "typically greater than 100×100×100"); rank $r\approx10$; windows of order 100 snapshots.
- GaAs, 5 kV cm⁻¹, lattice at 300 K: a 1.7 ps window (170 of 1250 snapshots, total 12.5 ps) reproduces the steady-state drift velocity (Fig. 2e).
- Full velocity–field curves: cost lower "by a factor of 10.5 for GaAs and ∼16.5 for graphene" (Fig. 3).
- GaN, electrons ≈1 eV above the band edge, 400 fs of dynamics: DMD on $f$ needs a 200 fs window; DMD on $\delta f$ needs 50 fs (≈12 %) and returns the mean relaxation rate "5.23 eV fs⁻¹" within 0.8 % (as printed; looks inconsistent because ≈1 eV lost in ≈0.4 ps suggests eV ps⁻¹).
- Time step: "femtosecond", value not stated; snapshot spacing 10 fs (my calculation from 170 snapshots in 1.7 ps).

**Key results.**
- (c) Singular values decay fast: about ten modes describe 10⁵ populations (computed, Fig. 2b).
- Subtracting the known equilibrium shortens the required window four-fold (computed, Fig. 4).

**Assumptions and limitations.**
- DMD extrapolates one trajectory; it does not approximate the collision operator and is re-learned for every field value (stated for field steps).
- A linear, time-invariant propagator suits autonomous dynamics only: constant field or free relaxation (my reading).
- (a), (b): nothing on thermal versus non-thermal energy partition or on temperature-dependent material parameters.

**Use for this project.**
- As a Maxwell closure this form of DMD is unusable: the pump makes each cell non-autonomous, and the closure needs $Q[f]$ at arbitrary $f$, not the continuation of a known trajectory (my reading).
- The benchmark it sets is 10–17 times per trajectory; the WP2 criterion is 100 times per evaluation.
- The low rank suggests a projection-based reduced model (proper orthogonal decomposition of $f-f_\text{FD}$ with a Galerkin-projected B0) as a further non-ML baseline (my reading).
- Subtracting the matched Fermi–Dirac state helped here too, which supports the WP2 input choice.
- Code and data: the solver is the open-source Perturbo package; no DMD code or data link in the arXiv version.

---

## Neural operators and learned kinetic closures
### [42] Li et al., "Fourier neural operator for parametric partial differential equations" — ICLR (2021)

**Source read:** arXiv:2010.08895v3 (ICLR version), full text + TeX source.
**Role in this project:** Defines the Fourier layer that WP2 adapts to the energy axis; origin of the resolution-transfer claim behind the "train at 256, evaluate at 512" criterion.

**What it does.** Proposes a neural operator whose kernel integral is a convolution parameterised in Fourier space and evaluated by FFT. On 1-D Burgers, 2-D Darcy and 2-D Navier–Stokes it reports the lowest errors among the learned solvers compared, errors independent of grid resolution, and much faster inference than a pseudo-spectral solver.

**Model / method.** The input function $a(x)$ is lifted pointwise, $v_0(x)=P(a(x))$, updated $T$ times, and projected pointwise, $u(x)=Q(v_T(x))$. Layer update (Eq. 2) and Fourier integral operator (Eq. 4):

$$v_{t+1}(x) := \sigma\Big( W v_t(x) + \bigl(\mathcal{K}(a;\phi)v_t\bigr)(x) \Big), \qquad \bigl(\mathcal{K}(\phi)v_t\bigr)(x)= \mathcal{F}^{-1}\Bigl(R_\phi \cdot (\mathcal{F} v_t) \Bigr)(x), \qquad \forall x \in D$$

$v_t: D\to\mathbb R^{d_v}$ are hidden fields, $W\in\mathbb R^{d_v\times d_v}$ a pointwise linear map, $\sigma$ a pointwise nonlinearity, $\mathcal F$ the Fourier transform. $R_\phi$ is a complex $k_{\max}\times d_v\times d_v$ tensor acting on modes $|k_j|\le k_{\max,j}$; higher modes are discarded. Assumptions: a translation-invariant, periodic kernel $\kappa(x-y)$; a uniform grid for the FFT.

**Numbers.**
- Network (Sec. 5): 4 Fourier layers, ReLU, batch normalisation; $k_{\max,j}=16$, $d_v=64$ in 1-D; $k_{\max,j}=12$, $d_v=32$ in 2-D. 1000 training, 200 test samples; Adam, 500 epochs, learning rate 10⁻³ halved every 100.
- Burgers, viscosity 0.1, relative L2 error, trained and tested at the same resolution $s$ (Table 3): 0.0149 ($s=256$), 0.0158 (512), 0.0160 (1024), 0.0146 (2048), 0.0142 (4096), 0.0139 (8192).
- Darcy (Table 4): 0.0098–0.0109 for $s$ = 85–421.
- Navier–Stokes at 64×64 (Table 1): FNO-3D 0.0086 (viscosity 10⁻³), 0.1918 (10⁻⁴, 1000 samples), 0.0820 (10⁻⁴, 10 000 samples).
- Inference 0.005 s against 2.2 s for the pseudo-spectral solver at 256×256 (Sec. 5.5).
- Cost per layer: $O(n\log n)$ with FFT; $O(n k_{\max})$ for a direct truncated transform on $n$ points (Sec. 4).

**Key results.**
- Error is flat in $s$ when training and testing at the same $s$ (measured, Tables 3–4); these are not cross-resolution tests (my reading).
- Zero-shot super-resolution is shown by one example: trained at 64×64×20, evaluated at 256×256×80, viscosity 10⁻⁴; no error is given (simulated, Fig. 1).
- Truncating a Navier–Stokes solution at 20 modes costs about 2 %; the network reaches ≤ 1 % with $k_{\max,j}=12$ because the nonlinearities regenerate higher modes (measured, Sec. 5.5).
- Non-periodic cases (Darcy, the time axis) work because "the linear transform $W$ (the bias term) … keeps the track of non-periodic boundary" (asserted).

**Assumptions and limitations.**
- Inputs are functions only (initial condition, coefficient field, first ten frames). Scalar parameters are not inputs: each viscosity has its own model. $R$ may depend on $\mathcal F a$; a linear dependence performed like the direct form, a neural-network one worse (asserted, Sec. 4).
- Aliasing, padding and cosine transforms are not discussed. Modes created by the nonlinearity above $k_{\max}$ are represented differently on different grids (my reading).

**Use for this project.**
- Adopt Eqs. (2) and (4) with $k_{\max}$ well below the 256-point Nyquist limit, so the same modes exist on both grids (my reading).
- No number here supports "256 to 512 with less than twofold error growth"; cross-resolution errors are in [45], Table 1.
- Non-periodic energy axis: only the $W$ argument here; cosine transform and padding are in [43], [45].
- Material parameters and lattice temperature must enter as constant input channels or through $R$; neither is tested here.
- Conflict: [43], Table 3, gives 0.0018–0.0020 for the same Burgers problem.

---
### [43] Kovachki et al., "Neural operator: learning maps between function spaces with applications to PDEs" — J. Mach. Learn. Res. 24, 89 (2023)

**Source read:** arXiv:2108.08481v6, text + TeX. Read: abstract, Secs. 1, 2.3, 3, 4.4–4.5, 5.1, 7.2–7.4, 8.2–8.3 (statements), 10. Skipped: Secs. 4.1–4.3, 5.2, 6, 7.1, 7.3.4, 9, proofs, appendices.
**Role in this project:** Theory behind the WP2 operator; source for the cosine-transform option.

**What it does.** Defines neural operators (pointwise maps composed with linear integral operators), proves discretisation invariance and universal approximation, and compares four kernel parameterisations on PDE benchmarks.

**Model / method.** One layer (Eq. 10): hidden field $v_t$, matrix $W_t$, kernel $\kappa^{(t)}$, measure $\nu_t$, bias function $b_t$, fixed map $\Pi_t$.

$$v_{t+1}(x) = \sigma_{t+1} \left ( W_t v_t( \Pi_{t} (x)) + \int_{D_t} \kappa^{(t)} (x,y) v_t(y) \: \text{d}\nu_t(y) + b_t(x) \right )$$

Kernel cost on $J$ points (Eqs. 28–31): graph (GNO), truncated Nyström sum, $O(JJ')$; low-rank (LNO), $O(J)$; multipole graph (MGNO), $O(J)$; Fourier (FNO), $u=\mathcal F^{-1}(R_\phi\cdot\mathcal F v)$, $O(J\log J)$.

**Numbers.**
- Burgers, relative L2, same train and test resolution $s$ (Table 3): FNO 0.0018 at $s$ = 256 and 512; DeepONet 0.0569, 0.0617, 0.0857 at $s$ = 256, 512, 8192. Darcy, $s=85$ (Table 2): 0.0108 and 0.0476.
- Suggested $k_{\max,j}$: one third to two thirds of the available modes (Sec. 4.4).

**Key results.**
- Theorem 8: for continuous functions and activations and fixed parameters, the discretised operator converges to the continuum one under mesh refinement (derived; no rate).
- Theorems 11–14: universal approximation on compact sets (derived); it "does not apply to the FNO formulation" (Sec. 4.4).
- FFT layers cover periodic grids "or … homogeneous Dirichlet (fast Fourier sine transform) or homogeneous Neumann (fast Fourier cosine transform) boundary conditions"; other domains by zero padding (asserted, Sec. 4.4).
- Training with input noise closes the clean/noisy gap (measured, Table 9).

**Assumptions and limitations.** Super-resolution is shown by pictures only (Figs. 9, 11).

**Use for this project.**
- Citable basis for the even-extension (cosine) layer; it implies zero slope of hidden fields at both ends of the mapped axis (my reading).
- No rate or number supports the 256-to-512 criterion.
- Conflicts: FNO Burgers error 0.0018 here, 0.0149 in [42]; the PDF header prints "JMLR 23 (2022) 1-97".

---
### [44] Lu et al., "Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators" — Nat. Mach. Intell. 3, 218 (2021)

**Source read:** arXiv:1910.03193v3 (15 Apr 2020), full text + TeX source. The preprint has three authors (Lu, Jin, Karniadakis) and the title "DeepONet: Learning nonlinear operators for identifying differential equations based on the universal approximation theorem of operators". The five-author journal version was not read.
**Role in this project:** Architecture of the DeepONet baseline that WP2 must report (its conservative form is [47]).

**What it does.** Turns the Chen–Chen universal approximation theorem for operators into a trainable architecture: a branch net encodes the input function at fixed sensors, a trunk net the output coordinate. On four ODE and PDE examples it generalises better than fully connected networks; error is characterised against sensor number, horizon and data-set size.

**Model / method.** Theorem 1 (Eq. 1): for any $\epsilon>0$ there exist sizes and weights with

$$\left|G(u)(y) - \sum_{k=1}^p \underbrace{\sum_{i=1}^n c_i^k \sigma\left(\sum_{j=1}^m \xi_{ij}^ku(x_j)+\theta_i^k\right)}_{branch} \underbrace{\sigma(w_k \cdot y+\zeta_k)}_{trunk} \right|<\epsilon$$

for all $u$ in a compact set of continuous functions and all $y$ in a compact set. The network (Eq. 2) uses deep sub-networks:

$$G(u)(y) \approx \sum_{k=1}^p b_k t_k + b_0$$

$G$ is the operator, $u$ the input function sampled at sensors $x_1,\dots,x_m$, $y$ the output location, $b_k$ the branch outputs, $t_k$ the trunk outputs, $b_0$ a bias. "Stacked" uses $p$ branch nets, "unstacked" one branch net with $p$ outputs. All inputs must use the same sensors; output locations are free. Sensor requirement (Theorem 2, Appendix C): $c(b-a)\kappa(m,V)e^{c(b-a)}<\varepsilon$, with $c$ a Lipschitz constant of the ODE right-hand side on $[a,b]$ and $\kappa(m,V)$ the piecewise-linear interpolation error of the input space; $\kappa\sim 1/(m^2l^2)$ for Gaussian random fields of correlation length $l$.

**Numbers.**
- Defaults (Tables 1–2): $m=100$; Gaussian random field inputs, $l=0.2$; 10⁴ training triplets $(u,y,G(u)(y))$; unstacked; trunk depth 3, branch depth 2, width 40 (ODEs) or 100 (PDE); Adam, learning rate 10⁻³, 5×10⁴ to 5×10⁵ iterations.
- Sensors, pendulum, horizon $T=1$ (Sec. 4.2.1): MSE ∝ 4.6^(−m) up to about 10 sensors, then flat; $T=5$ needs about 25. Fits: $m\propto\sqrt T$, $m\propto l^{-1}$.
- Horizon (Sec. 4.2.2): MSE ∝ 8^T at width 40 and 10⁴ points.
- Data-set size $x$ (Sec. 4.2.2): test MSE ∝ e^(−x/2000) below 10⁴ points, then $x^{-0.5}$.
- Diffusion–reaction PDE (Sec. 4.3): test MSE about 10⁻⁵ from 100 input functions with 1000 output points each.

**Key results.**
- The unstacked form with bias has the smallest test error and beats the best fully connected baseline (test MSE about 10⁻⁴) on the antiderivative (measured, Figs. 3–4).

**Assumptions and limitations.**
- Smooth, low-dimensional inputs; mean-squared errors only; no bounds on network size (authors).
- The branch net is tied to its sensors, so the input resolution cannot change after training (my reading; stated in [43], Sec. 5.1).

**Use for this project.**
- Baseline layout: branch input $f$ (or $f-f_{FD}$) at fixed energy sensors plus material parameters; trunk input $E$ (my reading).
- Sensors must resolve the narrowest input feature ($m\propto l^{-1}$); use the same mapped grid as the Fourier operator.
- A DeepONet trained with 256 sensors needs resampling to accept 512 points; the 256-to-512 test applies to its output side only.
- Comparison with Fourier operators (none in this paper; from others here): same data and resolution, relative L2 error. [43], Tables 2–3: Darcy 0.0476 against 0.0108; Burgers 0.0569 against 0.0018, the DeepONet error growing with resolution. [45], Table 2: Darcy 6.97 % against 1.22 % (Fourier operator with physics loss).
- The 8^T growth concerns a learned solution map; WP2 learns Q, where the horizon enters only through the time stepper (my reading).

---
### [45] Li et al., "Physics-informed neural operator for learning partial differential equations" — ACM/IMS J. Data Sci. 1 (2024)

**Source read:** arXiv:2111.03794v4 (29 Jul 2023), full text with appendices. The journal version was not compared.
**Role in this project:** Reference for adding physics losses to the WP2 Fourier operator; the only source in this group with cross-resolution error numbers.

**What it does.** Trains a Fourier neural operator with a data loss plus a PDE-residual loss imposed at a higher resolution than the data, and optionally fine-tunes it on one instance with the residual alone. The residual loss removes most error growth under super-resolution and reduces the data needed.

**Model / method.** Operator-learning losses (Eqs. 6, 7) for inputs $a_j\sim\mu$, solutions $u_j$ and operator $\mathcal G_\theta$; $\mathcal L_{pde}$ is the squared $L^2$ norm of the PDE residual plus weighted boundary and initial terms (Eqs. 3–4):

$$\mathcal J_{data}(\mathcal G_\theta) \approx \frac{1}{N}\sum_{j=1}^N \int_D |u_j(x) - \mathcal G_\theta(a_j)(x)|^2\,\mathrm dx, \qquad \mathcal J_{pde}(\mathcal G_\theta) = \mathbb E_{a\sim\mu}\big[\mathcal L_{pde}(a,\mathcal G_\theta(a))\big]$$

$\mathcal J_{pde}$ needs no labels, so fresh inputs can be drawn every iteration. Instance-wise fine-tuning (Sec. 3.2) minimises $\mathcal L_{pde}+\alpha\mathcal L_{op}$, with an anchor to the pre-trained parameters $\theta_0$:

$$\mathcal L_{op}\big(\mathcal G_{\theta_i}(a),\mathcal G_{\theta_0}(a)\big) := \Vert \mathcal G_{\theta_i}(a)-\mathcal G_{\theta_0}(a)\Vert^2_{\mathcal U}$$

Non-periodic outputs are handled by zero padding (Appendix C).

**Numbers.**
- Backbone: width 64, 8 or 12 modes, 4 layers, GeLU (Appendix A).
- Relative L2 error at the data resolution and two finer test grids (Table 1):

| Problem (data grid; residual grid) | Data only | Data + residual |
| --- | --- | --- |
| Burgers, viscosity 0.01 (32×25; 128×100) | 0.32, 3.32, 3.76 % | 0.17, 0.28, 0.38 % |
| Darcy (11×11; 61×61) | 5.41, 9.01, 9.46 % | 5.23, 1.56, 1.58 % |
| Kolmogorov flow (64×64×33; 256×256×65) | 8.28, 8.27, 8.30 % | 6.04, 6.02, 6.01 % |

  Test grids: Burgers 1×, 2×, 4× the data grid; Darcy 11², 61², 211²; Kolmogorov 64², 128², 256² in space (horizon 0.125). 800–1000 training samples.
- Data efficiency, Kolmogorov flow, horizon 0.125, 2200 residual instances (Table 3): 6.22 % with no data, 6.01 % with 800 samples, 5.04 % with 2200.
- Horizon 0.5 (Table 7): 4000 samples 25.2 %; plus 100 000 residual instances 24.2 %; residual only 74.4 %.
- Fine-tuning one instance (Table 4): 0.9 % in 473–608 s; a physics-informed network reaches 18.7 % in 4577 s.
- Transient flow, 4800 samples (Sec. 4.2): 3.04 % data only, 2.87 % with residual, 1.84 % after fine-tuning; 400 times faster than a GPU pseudo-spectral solver before fine-tuning.

**Key results.**
- A data-only Fourier operator lost a factor of ten at twice the training resolution on Burgers and nothing on Kolmogorov flow; the Darcy growth is attributed to an "unresolved" training grid (measured, Table 1).
- A label-free loss at the finer resolution removes most of that growth (measured, Table 1).
- Residual-only training works at the short horizon and fails at the longer one (measured, Tables 3, 7).
- Pre-training shortens fine-tuning but leaves its final error unchanged (measured, Table 4).

**Assumptions and limitations.**
- All physics terms are penalties (the anchor loss is called "a hard constraint" but is one too, my reading); the relative weight of $\mathcal J_{data}$ and $\mathcal J_{pde}$ is not stated.
- With all Fourier modes active an aliasing error "is inevitably present"; said to matter rarely because the mode count is far below the grid size (Sec. 2.4).

**Use for this project.**
- The 256-to-512 criterion is realistic only if 256 mapped points resolve the Fermi edge and the photon-energy steps (my reading); test it.
- Transferable recipe: label-free losses on unlimited sampled states at 512 points. For a collision integral these are the number and energy moments, the sign of entropy production and $Q[f_{FD}]=0$, not a PDE residual (my reading).
- Penalties buy robustness to resolution and some data efficiency, not exactness; exact conservation needs a construction as in [47] (my reading).
- Instance-wise fine-tuning takes minutes per instance: unusable inside a Maxwell run.

---
### [46] Xiao and Frank, "RelaxNet: a structure-preserving neural network to approximate the Boltzmann collision operator" — J. Comput. Phys. 490, 112317 (2023)

**Source read:** arXiv:2211.08149v1 (15 Nov 2022), full text + TeX source; Figs. 2, 7–15, 18, 19 inspected in the PDF. The published JCP version was not compared.
**Role in this project:** Direct prior art for WP2: a learned collision term Q[f] with built-in physical structure, embedded in a kinetic solver. It fixes what "soft conservation, architectural positivity and equilibrium" means in the novelty claim.

**What it does.** Replaces the fivefold classical Boltzmann collision integral by a learned BGK-type relaxation term: two fully connected networks ("RelaxNet") output a corrected target equilibrium and a velocity-dependent relaxation time. Advection and time stepping stay conventional (finite volume, Runge–Kutta); the result is called the "universal Boltzmann equation" (UBE). Trained on sampled (f, Q) pairs from the Shakhov model, a ν-BGK model and the full hard-sphere Boltzmann operator, it reproduces homogeneous relaxation, Mach 2–3 shock profiles and a lid-driven cavity at 50–108 times lower cost than the fast spectral method.

**Model / method.** The network learns the collision operator, not the flow map. Relaxation form (Eq. 13) inserted in the kinetic equation (Eq. 19):

$$\mathcal Q(f) = \nu\,(\mathcal E(f) - f), \qquad \partial_t f + \mathbf v \cdot \nabla_{\mathbf x} f = \mathrm{NN}_\theta (f,\bar\tau)$$

with (Eqs. 28, 16)

$$\mathcal E=\mathrm{ReLU}(\mathcal M \mathcal M^+), \quad \mathcal M^+=\exp \left(\alpha_1 + \boldsymbol\alpha_2 \cdot \mathbf v+\frac{1}{2}\alpha_3 \mathbf v^2\right), \quad \nu = \frac{1}{\tau}, \quad \tau = \bar\tau\, (1 + \mathrm{ELU} (\boldsymbol\beta_{\theta_2}))$$

- $f$ is the full distribution on a fixed grid of $N_v$ velocity points; $\langle\cdot\rangle$ is quadrature over velocity; $\boldsymbol\psi=(1,\mathbf v,\mathbf v^2/2)^T$ are the collision invariants.
- $\mathcal M(f)$ is the local Maxwellian built from the density, bulk velocity and temperature of $f$ (Eqs. 8, 15). $\bar\tau=\mu/p$ is the mean relaxation time (Eq. 14), with viscosity $\mu = \mu_\mathrm{ref}(T/T_\mathrm{ref})^\omega$ (Eq. 75) and pressure $p$.
- $\mathcal E$-net: $\mathbb R^{N_v} \to \mathbb R^{D+2}$, input $\mathcal M - f$ (the deviation from local equilibrium), output $\boldsymbol\alpha$. $\tau$-net (parameters $\theta_2$): $\mathbb R^{N_v+1}\to\mathbb R^{N_v}$, input $(\mathcal M-f,\bar\tau)$, output $\boldsymbol\beta_{\theta_2}$. Both use tanh hidden layers and weights without biases (Sec. 3.2).

Training cost and the condition on its weights (Eqs. 20, 21):

$$C(\theta) = \sum_{n,i,j} \Vert \mathrm{NN}_\theta - \mathcal Q_\mathrm{ref} \Vert (t^n,\mathbf x_i,\mathbf v_j) + \lambda \sum_{n,i} \Vert \langle \boldsymbol\psi\, \mathrm{NN}_\theta \rangle \Vert (t^n,\mathbf x_i) + \xi \sum_l^L \Vert \theta^l \Vert, \qquad \lambda N_t N_x > 1,\quad \xi N_t N_x N_v > 1$$

Guarantees proved (Theorems 2 and 4; Eqs. 30, 42), with $\varepsilon$ the minimum of $C$ and $\tilde{\boldsymbol\alpha}$ the multipliers of $\mathcal E$:

$$\Vert \langle \boldsymbol\psi\, \mathrm{NN}_\theta (g,\bar\tau) \rangle \Vert \leq \varepsilon, \qquad \partial_t \langle H(f) \rangle + \nabla_\mathbf x \cdot \langle \mathbf v H(f) \rangle \leq \Big( 1 + \sum_i \tilde\alpha_i \Big) \varepsilon$$

Training-data sampler (Eqs. 59, 63). Distributions are drawn from the entropy-closure family and joined across a random half-space:

$$f_u(\mathbf v)=\exp(\boldsymbol\gamma_u \cdot \boldsymbol\phi), \qquad f_\mathrm{ref} = f_L\, \mathscr H(\mathbf n \cdot \mathbf v) + f_R\, (1- \mathscr H(\mathbf n \cdot \mathbf v))$$

$\boldsymbol\phi$ holds $\boldsymbol\psi$ plus monomials up to degree $N_m$; the first three multipliers encode a Maxwellian with sampled $(\rho,\mathbf V,T)$, the others are random; samples are filtered by the condition number of the dual Hessian (Algorithm 1). Labels $\mathcal Q_\mathrm{ref}$ come from the chosen collision model.

How each structure is obtained:

| Structure | Mechanism | Kind |
| --- | --- | --- |
| $f \ge 0$ | $\mathcal E \ge 0$ (ReLU) and $\nu \ge 0$ (ELU $> -1$) for any $\theta$; relaxation form (Theorem 1) | hard, by architecture |
| Mass, momentum, energy | penalty $\lambda\Vert\langle\boldsymbol\psi \mathrm{NN}_\theta\rangle\Vert$ in Eq. (20); bound $\varepsilon$ "at any solution point in the dataset" (Theorem 2) | soft penalty |
| H-theorem | relaxation towards an exponential-family $\mathcal E$ plus Theorem 2; "holds asymptotically as $\varepsilon \to 0$" (Theorem 4) | up to training error |
| Maxwellian fixed point, Euler limit | bias-free nets fed $\mathcal M - f$: as $f\to\mathcal M$, $\boldsymbol\alpha \to 0$, $\mathcal M^+ \to 1$, $\mathcal E \to \mathcal M$ (proof of Theorem 6); exact at $f=\mathcal M$ because bias-free tanh layers map zero to zero (my reading) | hard, by architecture |
| Galilean invariance | asserted because the inputs are "scalar and Galilean invariant" (Theorem 5) | asserted |

**Numbers.**
- Velocity grids (dimensionless): 1-D $u\in[-8,8]$, 80 points (Table 6); 3-D $[-8,8]^3$, 80×28×28 (Table 7); shock $[-10,10]^3$, 80×28×28, 100 cells, $t\in(0,250]$, Ma 2 and 3 (Table 9); cavity $[-5,5]^3$, 28³, 45×45 cells, Kn = 0.075 (Table 11).
- Networks: 2 hidden layers, width equal to the input dimension (Sec. 6.2.1, 1-D case); 3-D sizes not stated.
- Training set: 10 000 sampled (f, Q) pairs per case (Tables 6, 7, 9, 11). Homogeneous cases: $\rho\sim\mathcal U(0.1,10)$, $U\sim\mathcal U(-1,1)$, $T\sim\mathcal U(0.1,5)$, $\gamma_4\sim\mathcal N(0,0.005)$, $\gamma_5$ truncated normal $(0,0.005)$ on $(-\infty,0]$, Kn $\sim\mathcal U(0.001,1)$, $N_m=4$.
- Optimiser Adam. Learning rate, epochs, batch size, $\lambda$, $\xi$: not stated. Final cost about $10^{-4}$ (Fig. 7; the axis appears logarithmic).
- Time integration: Tsitouras 5/4 Runge–Kutta, CFL 0.3 (homogeneous, $t\in(0,8]$, Kn = 1) or 0.5.
- Cost, homogeneous time series (Table 8): Boltzmann 2.167 s, 6.43 GB; BGK 1.076 ms; UBE 43.472 ms, 375.66 MB. Text: 1/50 the time, 1/17 the memory.
- Cost, one right-hand side, shock case (Table 10): Boltzmann 32.114 ms, 102.98 MB; BGK 5.893 µs; UBE 295.947 µs, 2.38 MB. Text: factor 108, 1/43 the memory. Hardware: not stated.
- Conservation in rollout, read from figures. Shakhov case (Fig. 10): density 0.85000 dips to about 0.84994 (7×10⁻⁵ relative) near $t=0.8$; energy rises from 1.1688 to about 1.1725 (+0.3 %), the reference to about 1.1711. ν-BGK case (Fig. 15): density 0.750 dips to about 0.7437 (−0.8 %), energy rises from 1.688 to about 1.735 (+2.8 %); the ν-BGK reference itself deviates more (0.7412, 1.769). The deviations are transient and have decayed by $t\approx4$–5.
- Error norms of Q or f: not stated anywhere; all accuracy comparisons are plots.

**Key results.**
- UBE tracks the full Boltzmann homogeneous relaxation where BGK deviates visibly (simulated, Figs. 16–17); the correction comes mainly from the velocity-dependent relaxation frequency, the learned equilibrium stays near Maxwellian (Fig. 18).
- The plotted entropy increases monotonically and follows the Boltzmann curve (simulated, Fig. 19).
- Shock profiles at Ma 2 and 3, including stress and heat flux, match the Boltzmann reference (simulated, Figs. 20–21); cavity velocity profiles match DSMC (simulated, Fig. 25).
- The bimodal initial conditions of the relaxation tests "can fall outside the training set" (Sec. 6.2.1). No rollout instability is reported.

**Assumptions and limitations.**
- Classical statistics, monatomic hard-sphere gas. $\mathcal E$ is restricted to exponentials of the collision invariants; the authors give this as the reason the Shakhov case is matched less well.
- Conservation is not exact, although Sec. 6.2.1 calls mass conservation "exactly satisfied" and the conclusion says "strictly preserves"; the introduction says conservation and the H-theorem hold "up to the training error".
- The shock case adds case-specific training data "on the fly" when residuals rise (Secs. 5.3, 6.3); it is not a pure out-of-distribution test and the speed-up excludes that cost (my reading).
- Dense networks tied to one velocity grid; no transfer between resolutions (my reading).
- "Same width as the input" would mean 62 720-wide layers in 3-D, inconsistent with the 375 MB reported; take the 3-D sizes from the code (my reading).
- As printed: $H(f)=-f\log f$ is stated with a dissipation inequality $\le 0$ (Eqs. 6–7), whereas the proof of Theorem 4 uses $f\log f$; the same proof sets $\tilde\alpha_i=\alpha_i\hat\alpha_i$ (a product) although $\mathcal M\mathcal M^+$ implies a sum; Fig. 2 draws the $\tau$-net input from $f$, the text says $\mathcal M - f$.
- Positivity is shown for the continuous equation, not for the explicit time step (my reading).
- Left open by the authors: multi-level relaxation models built from several RelaxNet blocks (Conclusion).

**Use for this project.**
- Novelty wording: [46] enforces conservation of mass, momentum and energy by a soft penalty; positivity and the equilibrium fixed point are architectural. An exact equilibrium fixed point is therefore not new. New relative to [46]: exact conservation, Fermi–Dirac statistics with Pauli blocking, a resolution-transferable operator.
- The plan's sentence that [46] is a structure-preserving variant of Fourier neural operators is inaccurate: RelaxNet is a fixed-grid fully connected relaxation model.
- Design to consider: feeding $f - f_{FD}$ to bias-free layers gives $Q[f_{FD}]=0$ without the second forward pass that $N[f]-N[f_{FD}]$ needs; this fails if other non-zero inputs (material parameters) enter additively (my reading).
- Numbers to beat: conservation drift of 10⁻⁴ to 10⁻² (WP2 target 10⁻⁶); 108-fold speed-up per evaluation (WP2 target 100).
- A learned relaxation model (target distribution plus energy-dependent rate) is a natural structured baseline next to the relaxation-time approximation [40]. For fermions the target must also stay below 1 and the entropy is the Fermi–Dirac one; neither is covered here (my reading).
- The entropy-closure sampler is a cheap source of the "perturbed off-trajectory states" of WP1.
- Pitfall: no error norms are reported, so there is no accuracy figure to quote.
- Code: "The open-source codes to reproduce the numerical results are available under the MIT license [1]" (reference [1] is Kinetic.jl).

---
### [47] Lee et al., "Structure-preserving operator learning: modeling the collision operator of kinetic equations" — arXiv:2402.16613 (2024)

**Source read:** arXiv v1 (26 Feb 2024), full text + TeX source; Figs. 1–3 inspected in the PDF. Authors: J. Y. Lee, S. Schotthöfer, T. Xiao (equal contribution), S. Krumscheid, M. Frank. No venue is stated in the paper (it is typeset in the ICML 2024 style).
**Role in this project:** Direct prior art for the hard conservation constraint of WP2, and the structure-preserving form of the DeepONet baseline that WP2 must report.

**What it does.** Learns the collision operator of the linear kinetic (radiation-transport) equation with a DeepONet whose trunk network is modified so that the output has zero mass moment for any input. Two constructions are compared with an unmodified DeepONet and a penalty-based ("soft") version. The hard versions conserve mass to 10⁻⁷–10⁻⁸ without losing accuracy; the penalty neither conserves better nor keeps accuracy. The surrogate is run inside a finite-volume solver on three test cases.

**Model / method.** The network learns the collision operator, not the flow map. Kinetic equation (Eq. 1; as printed, the absorption term lacks the factor $f$) and collision operator (Eq. 2):

$$\partial_t f+\mathbf{v} \cdot \nabla_{\mathbf{x}} f +\sigma_{\textup{a}} = \sigma_{\textup{s}}Q(f) + p, \qquad Q(f)(\mathbf{v})=\int_{\mathbf{S}^2} k(\mathbf{v}_*, \mathbf{v}) \left[ f(\mathbf{v}_*)-f(\mathbf{v})\right] d\mathbf{v}_*$$

$f(t,\mathbf x,\mathbf v)$ is the kinetic density, $\mathbf v$ lies on the unit sphere $\mathbf S^2$ (or in $[-1,1]$ for slab geometry), $k$ is the scattering kernel, $\sigma_\textup{a}$, $\sigma_\textup{s}$ are absorption and scattering coefficients, $p$ a source. The only collision invariant is $\varphi = 1$ (mass). Kernels: isotropic, and Henyey–Greenstein with anisotropy $g$ (Sec. 5.1):

$$k_g(\mathbf{v}_*,\mathbf{v}) = \frac{1-g^2}{4\pi(1-2g\,\mathbf{v}\cdot\mathbf{v}_* + g^2)^{3/2}}$$

Modified DeepONet (Eq. 8), with the free bias of the standard DeepONet replaced by an extra basis function:

$$Q_\theta(f)(\mathbf{v})=\sum_{k=1}^{p+1}  \beta_k(\bar{f})\tau_k(\mathbf{v})$$

$\bar f = [f(\mathbf v_1),\dots,f(\mathbf v_m)]$ are the values of the full distribution at $m$ sensor points, which coincide with the velocity quadrature nodes; the branch net gives coefficients $\beta_k$, the trunk net basis functions $\tau_k$; $\langle\cdot\rangle$ is the velocity integral.

- Method I, orthogonal trunk (Sec. 2.2.1, Eq. 9): set $\tau_{p+1}=\varphi$, $\beta_{p+1}=0$ and orthonormalise the extended set $[\varphi,\tau]$ before each forward pass, so that

$$\langle\tau_{k_1}\tau_{k_2}\rangle=0,\quad\forall\, 1\leq k_1\neq k_2\leq p+1$$

- Method II, bias adaption (Sec. 2.2.2, unnumbered): set $\tau_{p+1}(\mathbf v)=1$ and

$$\beta_{p+1}(\bar{f}) =  -\sum_{k=1}^p  \beta_k(\bar{f})\, \frac{\langle\varphi\tau_k\rangle}{\langle\varphi\rangle}$$

- Soft baseline (Eq. 11; as printed, the index $k_1$ appears twice in the penalty, presumably $\tau_{k_1}\tau_{k_2}$):

$$\mathcal{L}(Q_\theta,Q)= \frac{1}{T}\sum_{i=1}^T \lVert Q(\bar{f}_i)-Q_\theta(\bar{f}_i)\rVert_2^2 +\lambda \sum_{1\leq k_1\neq k_2\leq p+1}\lvert\langle\tau_{k_1}\tau_{k_1}\rangle\rvert$$

The hard variants are trained with the first term only (mean-squared error on Q).

Training data (Sec. 3, Eq. 13; Algorithm 3). Distributions are sampled from the minimal-entropy closure with entropy density $\eta(f)=f\log f-f$, for which $\eta_*' = \exp$:

$$f_\mathbf{u}(\mathbf{v}) = \eta'_*(\boldsymbol{\alpha}_\mathbf{u}\cdot \mathbf m(\mathbf{v}))$$

$\mathbf m$ is a spherical-harmonics basis; $\boldsymbol\alpha$ is drawn from a Gaussian centred on the multiplier of the equilibrium distribution, with $\alpha_1$ fixed so that the sample has unit mass; samples with entropy above a threshold $c$ are rejected. At inference the input is normalised and the output rescaled (Eq. 15), and $Q_\theta$ is set to zero by hand when $\bar f=0$:

$$Q_\theta(f) = \langle f\rangle\, Q_\theta(\tilde{f}), \qquad \tilde f = f/\langle f\rangle$$

How each structure is obtained:

| Structure | Mechanism | Kind |
| --- | --- | --- |
| Mass | Method I or II | hard, by construction |
| Momentum, energy | not invariants of this operator; extension of Method I claimed, not shown | none |
| Entropy dissipation | "subject to future work" (Sec. 6); entropy only selects training samples | none |
| Equilibrium fixed point | only $Q(f=0)=0$ is imposed | none |
| $f\ge0$ | "naturally preserved since $Q_\theta$ … contains a limited number of neurons" (Sec. 4.1) | asserted |

**Numbers.**
- Architectures, fully connected, tanh (Table 2): slab, branch 100–16–16–16, trunk 1–16–16–16; sphere, branch 800–100–100–16, trunk 3–100–100–16. Sensors $m$ = 100 (slab), 800 (sphere), 32 (lattice test).
- Training: 10 000 epochs, Adam, learning rate 10⁻³, full batch; best model chosen by lowest training loss (Sec. 5).
- Not stated: training-set size $T$, sampling width $\sigma$, entropy threshold $c$, moment order, $\lambda$, hardware, code availability.
- Mean relative L2 test error of Q (Table 1):

| Case | Vanilla, no bias | Vanilla, bias | Soft | Method I | Method II |
| --- | --- | --- | --- | --- | --- |
| Slab, isotropic | 0.0132 | 0.0146 | 0.0135 | 0.0084 | 0.0104 |
| Slab, HG $g=0.9$ | 0.1064 | 0.1074 | 0.1771 | 0.1056 | 0.0425 |
| Sphere, isotropic | 0.1086 | 0.1123 | 0.1435 | 0.0710 | 0.1063 |
| Sphere, HG $g=0.9$ | 0.2020 | 0.2051 | 0.3234 | 0.1453 | 0.1950 |

- Mass-moment error $\lvert\langle\varphi Q_\theta\rangle\rvert$ on unseen data, same column order (Table 3): slab HG 1.17×10⁻³, 1.03×10⁻³, 1.58×10⁻³, 4.16×10⁻⁷, 7.62×10⁻⁸; sphere HG 7.58×10⁻⁴, 6.75×10⁻⁴, 5.38×10⁻³, 1.51×10⁻⁸, 1.74×10⁻⁸.
- Solver speed, inflow test, iterations per second (Table 4): numerical quadrature 1.06, vanilla 41.82, soft 39.91, Method I 2.80, Method II 40.65. Text: speed-ups of 2.2 (Method I; the table gives 2.6) and "almost 40" (Method II).
- Homogeneous relaxation, 50 steps to $t=5$, initial mass 1 (Fig. 2, read from the plot): vanilla drifts to about 1.004, soft to about 0.9983, Method II stays at 1.000.
- Solver runs (Fig. 2 table; Secs. 5.4–5.5): inflow, 1000 cells, 375 steps, CFL 0.4, HG $g=0.9$; lattice, 1600 cells, 74 steps, isotropic kernel.
- Inflow test: pointwise absolute error at or below about 10⁻² for all variants, for a solution of order 0.5 (Fig. 3).

**Key results.**
- Both hard constructions reduce the mass-moment error by four to five orders of magnitude relative to vanilla and soft versions (measured, Table 3).
- The penalty does not conserve better than no constraint and is less accurate in three of four cases (measured, Tables 1, 3). Explanation given: the penalty is trivially minimised by $\tau = 0$, so the trunk outputs shrink (Fig. 1, Sec. 5.2).
- The hard variants are at least as accurate as the vanilla DeepONet in all four cases (measured, Table 1).
- Method I needs a stable orthogonalisation; Gram–Schmidt gives strongly oscillating basis functions (Fig. 1) and the solver runs about 15 times slower than with Method II (Table 4).
- In the inflow test the surrogate sees inputs "close to the boundary or even outside the convex hull of the training data" without reported instability; the lattice result is compared with the quadrature baseline only by plots (Fig. 4).

**Assumptions and limitations.**
- Linear operator, bounded velocity space, one invariant. The authors: "not directly applicable for the nonlinear Boltzmann equation with unbounded velocity domain".
- Accuracy is modest (1–20 % on Q) and rollouts are short (50–375 steps).
- Eq. (15) relies on Q being homogeneous of degree one in $f$ (my reading).
- The multi-invariant extension rests on the remark that invariants "are typically pairwise orthonormal"; $1$, $\mathbf v$, $\lvert\mathbf v\rvert^2$ are not orthogonal under a general quadrature weight (my reading).
- Conservation holds for the discrete quadrature on the fixed sensor set, to which the branch net is tied; residuals of 10⁻⁷–10⁻⁸ look like single-precision round-off (my reading).
- As printed: Fig. 1 refers to $p=8$, Table 2 has 16 outputs; the inflow text gives $x\in[0,1]$ with $\Delta x=0.01$, the settings table 1000 cells, and Fig. 3 plots $x$ up to 6.

**Use for this project.**
- Novelty wording: [47] is a hard constraint by construction, not a soft one. Method II is algebraically a projection, $Q_\theta=\sum_{k\le p}\beta_k\,(\tau_k-\langle\varphi\tau_k\rangle/\langle\varphi\rangle)$, which removes the mass-violating part of every basis function (my reading). A hard projection for one invariant of a linear operator therefore exists. Accurate statement: "[46] enforces conservation by a loss penalty; [47] enforces a single invariant (mass) exactly, by construction, for a linear operator". For [49] see its own digest.
- Not provided by [47]: two invariants with density-of-states weights (number, energy), a nonlinear Pauli-blocked operator, an equilibrium fixed point ([46] has one for the Maxwellian), an entropy constraint.
- The plan calls [47] a structure-preserving variant of Fourier operators for the classical Boltzmann operator; it is a DeepONet for linear transport.
- Baseline: implement "DeepONet [44, 47]" as Method II extended to two invariants; it is free at run time (40.65 against 41.82 iterations per second).
- The WP2 target of 10⁻⁶ conservation "by construction" matches what is measured here; use double precision (my reading).
- Pitfall for the planned entropy-production penalty: a penalty that zero output satisfies shrank the network output here.
- Do not copy the normalisation of Eq. (15): the fermionic operator is not homogeneous in $f$ (my reading).
- Sampling idea for off-trajectory states: Gaussian perturbations of the equilibrium Lagrange multipliers with an entropy cut-off.

---
### [48] Lee et al., "FourierSpecNet: neural collision operator approximation inspired by the Fourier spectral method for solving the Boltzmann equation" — Kinet. Relat. Models 21, 69 (2026); arXiv:2504.20408

**Source read:** arXiv v2 (6 Mar 2026) full text + TeX source, including Appendix B (ablations). Proofs of Prop. 3.1 and Appendix A (decay of kernel modes) skimmed only. The journal reference is not printed in the arXiv PDF (not verified). The project plan cites v1; v2 is the referee-revised version and its timing tables differ from v1 (see Numbers).
**Role in this project:** Closest prior art for WP2's central design choice (learn the collision operator Q[f], integrate with a standard time stepper) and the architecture/speed baseline the WP2 operator will be compared with.

**What it does.** The classical fast Fourier spectral method writes the Boltzmann collision operator as a sum of weighted convolutions of Fourier coefficients, with weights fixed by quadrature of the collision kernel. FourierSpecNet keeps that exact algebraic form and makes the weights trainable complex numbers on a fixed low-frequency band. It is trained on pairs (f, Q(f,f)) from the classical solver and then integrated with RK3. Because the parameters live on a fixed Fourier band, the same weights apply at any velocity-grid size. Results cover 2D Maxwellian, hard-sphere and inelastic kernels and one 3D case, all classical and spatially homogeneous.

**Model / method.** Collision operator (elastic or inelastic, restitution coefficient $e$; kernel $B=C|\boldsymbol v-\boldsymbol v_*|^{\alpha}$, $\alpha=0$ Maxwellian, $\alpha=1$ hard sphere):

$$Q(f,f)=\int_{\mathbb{R}^d}\int_{S^{d-1}} B(|\boldsymbol{v}-\boldsymbol{v}_*|,\cos\theta)\left(f(\boldsymbol{v}')f(\boldsymbol{v}_*')-f(\boldsymbol{v})f(\boldsymbol{v}_*)\right)d\omega\,d\boldsymbol{v}_* \qquad (1.2)$$

$f$ is assumed supported in a ball of radius $S$ and periodized on $D_T=[-T,T]^d$, $T=\frac{3+\sqrt2}{2}S$ (here $T$ is the half-width of the velocity box, not temperature). With $f_N=\sum_{\|k\|_\infty\le N}\hat f_k e^{i\frac{\pi}{T}k\cdot\boldsymbol v}$, Galerkin projection gives a weighted convolution in Fourier-index space:

$$\hat{Q}_k(f_N,f_N)=\sum_{\substack{\|l\|_\infty,\|m\|_\infty\le N\\ l+m=k}}\hat{B}(l,m)\,\hat f_l\hat f_m,\qquad \hat B(l,m)=\hat G(l,m)-\hat G(m,m) \qquad (2.7)$$

$$\hat{G}(l,m)=\int_{B_0(2S)}\int_{S^{d-1}}B(|q|,\hat q\cdot\omega)\,e^{\,i\frac{\pi}{T}q\cdot\left(\frac{1+e}{4}l-\frac{3-e}{4}m\right)}\,e^{-i\frac{\pi}{T}|q|\omega\cdot\frac{1+e}{4}(l+m)}\,d\omega\,dq \qquad (2.8)$$

Direct evaluation costs $O(N^{2d})$. The fast method (Gamba et al.) approximates the kernel modes by a separable sum from radial/angular quadrature ($M_{\rm quad}-1=N_rN_\sigma$ terms), giving FFT-evaluable convolutions:

$$\hat G^{\rm fast}(l,m)=\sum_{t=1}^{M_{\rm quad}-1}\gamma_t(l+m)\,\alpha_t(l)\,\beta_t(m) \qquad (2.11)$$

$$\hat Q_k^{\rm fast}(f_N,f_N)=\sum_{t=1}^{M_{\rm quad}}\Big[\gamma_t(k)\cdot\big((\alpha_t\hat f)*(\beta_t\hat f)\big)_k\Big],\qquad (\alpha_t\hat f)*(\beta_t\hat f)=\mathrm{FFT}\big[\mathrm{FFT}^{-1}(\alpha_t\hat f)\cdot\mathrm{FFT}^{-1}(\beta_t\hat f)\big] \qquad (2.13),(2.14)$$

(the loss term is absorbed as the extra term $t=M_{\rm quad}$). Because the integrand oscillates, $N_r\gtrsim O(N)$, so the fast method costs $O(MN^{d+1}\log N)$, $M=N_\sigma$. FourierSpecNet replaces the quadrature weights by trainable complex arrays, zero for $\|k\|_\infty>N_{\rm trun}$:

$$\hat Q^{nn}_k := \sum_{t=1}^{M}\Big[\gamma^{nn}_t(k)\big((\alpha^{nn}_t\hat f)*(\beta^{nn}_t\hat f)\big)_k\Big] \qquad (3.1)$$

with $M$ a fixed hyperparameter, so the cost is $O(MN^d\log N)$. There are no activation functions or hidden layers; the model is exactly bilinear in $f$. Training loss (supervised, on $Q$):

$$\mathrm{Loss}=\frac{1}{|\mathcal D|}\sum_{f\in\mathcal D}\frac{\|Q^{nn}(f,f)-Q(f,f)\|_{L^2(D_T)}}{\|Q(f,f)\|_{L^2(D_T)}} \qquad (3.2)$$

Consistency (Prop. 3.1, assuming $\hat G(l,m)\to0$ at large $\|l\|,\|m\|$, proved for VHS with rational $e$): for any $\varepsilon>0$ there exist $M$, $N_{\rm trun}$ and weights such that

$$\|Q(f,f)-Q^{nn}(f_N,f_N)\|_{L^2(D_T)}\le C\Big(\|f-f_N\|_{L^2(D_T)}+\frac{1}{N^r}\|Q(f_N,f_N)\|_{H^r(D_T)}\Big)+\varepsilon\|f\|^2_{L^\infty(D_T)} \qquad (3.5)$$

This is an existence statement, not a bound on the trained network.

**Numbers.**
- Default model: $N_{\rm trun}=8$, $M=2$, 3072 real parameters (App. B table); Adam, learning rate $10^{-2}$; training stopped when relative $L^2$ error on $Q$ < $10^{-2}$ (Sec. 4).
- Training data: 3000 input functions (1000 Gaussians with centre in $[-1,1]^d$, width 0.8–1.2; 1000 sums of two Gaussians; 1000 Gaussians times $1+$ random quadratic polynomial), each normalized to unit mass; targets $Q(f,f)$ from the fast spectral method of Gamba et al. (Sec. 3.1). No time trajectories are used. Training resolution $N=64$ per dimension (2D) (Sec. 4.1–4.3).
- Training cost on one A100 GPU: $2\times10^5$ epochs, 1 h 39 min (Maxwellian); $2\times10^5$, 1 h 42 min (hard sphere); $5\times10^5$, 4 h 19 min (inelastic) (Table 1).
- Time integration: RK3, $\Delta t=0.01$, to $t=5$ (Maxwellian), $t=2$ (hard sphere), $t=3$ (inelastic) (Sec. 4).
- Rollout error in $f$ (relative $L^2$, read from figures): Maxwellian against the exact BKW solution (text; the caption of Fig. 5 also mentions the spectral method) rises then plateaus at ≈$1\times10^{-2}$ for $N=32,64,128$ and ≈$7\times10^{-3}$ for $N=16$ (Fig. 5); hard sphere vs spectral solver plateaus between $10^{-3}$ and $10^{-2}$ (Fig. 8, only one tick labelled); inelastic $e=0.2$ grows monotonically to ≈$5\times10^{-3}$ at $t=3$ (Fig. 10) and stays below 3.2% at $t=5$ for $N_{\rm trun}\ge8$; $N_{\rm trun}=4$ gives 23.9% (App. B.1).
- Resolution transfer: trained at $N=64$, tested at 16, 32, 64, 128 with indistinguishable error curves (Figs. 5, 8, 10); ablations also run at $N=256$.
- Time to solution to $t=5$ (same A100 GPU, PyTorch, both methods), fast spectral / FourierSpecNet: $N=16$: 0.40 s / 0.51 s (×0.78); 32: 0.52 / 0.52 (×1.00); 64: 1.57 / 0.59 (×2.65); 128: 6.02 / 0.57 (×10.6); 256: 22.5 / 0.72 (×31.2); 512: 90.7 / 1.31 (×69.0) (Table 2, Maxwellian; Tables 3–4 give ×72 at $N=512$ for hard sphere and inelastic $e=0.5$).
- $M$ ablation ($e=0.2$): accuracy improves from $M=2$ to $M\ge5$, no gain beyond; training time 4 h 19 min ($M=2$) to 11 h 17 min ($M=16$) (App. B.2).
- 3D case: no error or timing numbers are given (Sec. 4.4, Fig. 11).

**Key results.**
- A 3072-parameter bilinear spectral ansatz reproduces the 2D collision operator to about 1% and gives stable RK3 rollouts over 200–500 steps (simulated).
- Zero-shot change of grid size leaves the rollout error unchanged within the trained frequency band (simulated).
- Speed-up over the classical fast spectral solver appears only for $N\ge64$ per dimension and grows with $N$ because the learned rank $M$ does not grow with $N$; at $N\le32$ the network is no faster (simulated).
- Mass, momentum and kinetic energy stay visually constant and entropy $\int f\ln f$ decreases monotonically in the plotted runs (Figs. 4, 7); no conservation error is quantified (simulated, graphical only).
- With a poor learning rate ($10^{-4}$) momentum conservation visibly deteriorates (App. B.3), i.e. conservation is learned, not guaranteed.

**Assumptions and limitations.**
- Classical Boltzmann statistics only. No quantum Boltzmann, Fermi–Dirac statistics or Pauli blocking is mentioned anywhere.
- Spatially homogeneous; no forcing or source term; one network per collision kernel (per $e$), not parametric.
- Conservation, entropy decay and positivity are not enforced. The moment statements (3.4) concern Fourier truncation of $f$, not the learned operator (my reading).
- Periodized velocity box with compactly supported $f$; authors name boundary conditions as future work.
- The $N=16$ error is lower only because of aliasing/smoothing (authors' statement).
- The timing baseline is the authors' own GPU implementation of the fast spectral method. A commented-out table in the TeX, labelled as the earlier arXiv version, lists 5.74 s to 251.5 s for $N=16$–128 against 0.6–0.8 s (×9 to ×310), so a speed-up quoted from v1 would overstate the gain now reported (my reading; v1 itself not read).

**Use for this project.**
- Adopt the workflow: supervised training on (f, Q) pairs built from equilibria plus structured perturbations, then integration by a standard stepper. The Fermi–Dirac analogue of the three training families is the natural starting dataset.
- Calibration point: an operator-level error of $10^{-2}$ gave solution-level error $10^{-2}$–$10^{-3}$; WP2 asks for 1% on Q and 1% on plasma frequency over 2 ps.
- Baseline to report: FourierSpecNet-style bilinear low-rank model. It cannot represent the Uehling–Uhlenbeck operator as is, because Pauli blocking adds terms cubic in $f$ (my reading).
- Pitfall for Gate G1: the gain comes from removing the $O(N)$ quadrature rank in $d\ge2$; on a 1D energy grid an FFT-evaluated integral is already $O(N\log N)$, so this paper gives no reason to expect a large speed-up over B0 (my reading).
- Conflict with the plan: the plan says energy conservation gives the convolution structure "that makes spectral Boltzmann solvers work [48]". In [48] the convolution is over Fourier indices ($l+m=k$) and follows from the product $f f_*$ and the dependence of the kernel on $\boldsymbol v-\boldsymbol v_*$ after periodization, not from an energy-conserving delta function in an energy variable (my reading). The analogy should be reworded.
- Supports the novelty claim: no degenerate statistics, electron–phonon coupling or driving.
- No code or data-availability statement.

---
### [49] Hu and Qi, "Solving the BGK model and Boltzmann equation by Fourier neural operator with conservative constraints" — arXiv:2510.13047 (2025; v2 2026)

**Source read:** arXiv v2 (12 Nov 2026) full text + TeX source, all sections. Figures inspected for the BGK error curves. No journal reference is given; the acknowledgments mention anonymous reviewers.
**Role in this project:** Reference for the WP2 ablation (autoregressive flow-map operator) and the only quantitative evidence in the library on how far soft conservation penalties get for a Boltzmann surrogate.

**What it does.** A standard Fourier neural operator (lifting, stacked Fourier layers, projection) is trained to map a velocity distribution at one time to the distribution at a later time for the spatially homogeneous BGK model (1V, 3V) and the Boltzmann equation with Maxwell molecules (2V, 3V). It learns the solution (flow) map, not the collision operator. A "conservative" variant, C-FNO, adds penalties on the change of mass, momentum and energy to the loss. C-FNO reduces moment drift relative to plain FNO but stays orders of magnitude less conservative than the classical spectral solver.

**Model / method.** Equations solved, $\partial_t f=\mathcal Q[f]$, with

$$\mathcal{Q}_{\text{Bol}}(f,f)(v) := \frac{1}{\text{Kn}}\int_{\mathbb{R}^d}\int_{\mathbb{S}^{d-1}}\mathcal{B}(v-v_*,\sigma)\left[f(v')f(v_*')-f(v)f(v_*)\right]d\sigma\,dv_* \qquad (2.3)$$

$$\mathcal{Q}_{\text{BGK}}(f)(v) := \frac{f_{\text{eq}}(v)-f(v)}{\text{Kn}},\qquad f_{\text{eq}}=\frac{\rho}{(2\pi T)^{d/2}}\exp\Big(-\frac{|v-u|^2}{2T}\Big) \qquad (2.7),(2.8)$$

Kn = 1 for Boltzmann; kernel $B=b(\cos\theta)|v-v_*|^{\gamma}$ with $\gamma=0$ only. Learned object: $f_t(v)=\mathcal G_t[f_0](v)\approx\mathcal G_\theta[f_0](v)$ (3.4). Network: input channels $[v, f_0(v)]$ lifted by $L$ (3.6); Fourier layer

$$\tilde F_t(k)=R_t(k)\hat F_{t-1}(k)\ \ (\|k\|_\infty\le k_{\max}),\ 0\ \text{otherwise};\qquad F_t(v)=\sigma\big(W_tF_{t-1}(v)+\mathcal K_\phi[F_{t-1}](v)\big) \qquad (3.8),(3.9)$$

with $R_t(k)\in\mathbb C^{N_c\times N_c}$ learnable (here $t$ indexes layers), then projection $P$. Loss:

$$\mathcal{L}_{\text{Cons}} := \mathcal{L}_{\text{MSE}} + \lambda_m \mathcal{C}_m + \lambda_p \mathcal{C}_p + \lambda_e \mathcal{C}_e \qquad (3.13)$$

$$\mathcal C_m=\frac1n\sum_{j=1}^n\Big|\int f_\theta^{(j)}dv-\int f_0^{(j)}dv\Big|^2,\quad \mathcal C_p=\frac1n\sum_j\Big\|\int v f_\theta^{(j)}dv-\int v f_0^{(j)}dv\Big\|_2^2,\quad \mathcal C_e=\frac1n\sum_j\Big|\int|v|^2f_\theta^{(j)}dv-\int|v|^2f_0^{(j)}dv\Big|^2$$

Two training modes: "point-to-point" (one endpoint map $f_0\mapsto f_t$) and "sequence-to-sequence" (one-step map applied autoregressively, loss summed over the rollout, Eq. 3.15). Link to the spectral method (Remark, Eq. 3.34): one Euler step of the fast spectral scheme is $f_{t+1}=f_t+\Delta t\,\mathcal F^{-1}[\alpha^p_k((\beta^p_m\hat f_{t,m})*\hat f_{t,l})]$, which has the form of the FNO kernel term (3.11); the difference is that the spectral operator is quadratic in $f$ whereas FNO layers are linear convolutions, the nonlinearity being left to the activation. No part of the architecture is derived from the collision kernel.

**Numbers.**
- BGK data: Kn = $\tau$ = 0.1; $v\in[-3,3]$, 64 points (1V) or $32^3$ (3V); forward Euler, $\Delta t=0.01$, 100 steps; 100 trajectories per case; initial data = Maxwellian ($\rho_0\in[0.8,1.2]$, $T_0\in[0.6,1.4]$, $u_0=0$) times $1+\varepsilon g(v)$, $g$ a sum of six cosines; $\varepsilon=0.2$ (near equilibrium) or 3.0 (far) (Sec. 3.3.1, 4.1.1).
- Boltzmann data: fast spectral solver (Hu–Qi), Maxwell molecules, domain $[-L,L]^d$, $L=\frac{3+\sqrt2}{2}S$, $S=2.5$, $N=32$ per direction, RK4 with $\Delta t=0.1$. Initial conditions: 300 (2V) or 50 (3V) from Gaussians, BKW-type states and perturbed Gaussians (Sec. 3.3.2).
- Networks: 1V BGK $N_F=32$ modes, $N_C=64$ channels, $N_L=4$ layers, 0.57 M parameters, 18.5 s per epoch on a MacBook M2, 200 epochs (Table 1, Sec. 4.1.1). 2V Boltzmann $N_F=16^2$, $N_C=64$, $N_L=4$, Adam $10^{-3}$, batch 8, 150 epochs; $\lambda_m=\lambda_p=\lambda_e$, value not stated (Sec. 4.2.1). 3V BGK $N_F=16$, $N_C=32$, $N_L=4$, $\lambda=5\times10^{-3}$ (Sec. 4.3.1). 3V Boltzmann $N_F=8^3$, $N_C=16$, $N_L=4$, 3000 optimizer steps, 417 s on CPU (Sec. 4.3.2, Table 15).
- 1V BGK autoregressive rollout, 100 steps, $L^2$ error (Figs. 2, 6): near equilibrium plateaus at ≈$9\times10^{-3}$ ($N_F=32$), ≈$1.3\times10^{-2}$ ($N_F=16$), ≈$1.8\times10^{-2}$ ($N_F=8$); far from equilibrium peaks at ≈0.10 near step 7 then relaxes to ≈$10^{-2}$ ($N_F=32$). $N_F=64$ overfits and oscillates; $N_L=3$ underfits.
- 3V BGK rollout, relative $L^2$ (Fig. 30): peaks ≈4.5% (FNO) and ≈3.7% (C-FNO) near step 15, settles at ≈3.5% and ≈2.6% by step 100.
- 2V BKW benchmark, $L^2$ error in $f$ at $t=1,4,7$ (Tables 5–7): C-FNO 1.09e-2, 1.45e-2, 2.18e-2; FNO 1.45e-2, 2.49e-2, 2.70e-2; spectral solver 8.40e-3, 7.20e-3, 7.14e-3.
- Conservation in the same test, $|E^{\rm Num}-E^{\rm True}|$ at $t=1,4,7$: FNO 7.12e-1, 8.87e-1, 5.64e-1; C-FNO 2.68e-3, 9.20e-2, 2.93e-2; spectral 2.92e-6, 7.53e-5, 3.09e-4. Mass error: FNO 3.40e-2, 2.75e-2, 3.67e-3; C-FNO 2.11e-3, 2.56e-2, 4.77e-4; spectral 2.26e-13.
- 3V BKW endpoint map $t=5.5\to6.5$ (Table 13): $L^2$ error C-FNO 6.78e-4, FNO 6.97e-4, spectral 6.49e-4; mass error 2.9e-5, 1.1e-3, 1.2e-9; energy error 5.5e-4, 7.4e-4, 1.9e-5.
- 2V unseen Gaussian initial data (width 0.2–1.0), $L^2$ error against the spectral solver at $t=1,4,7$ (Tables 10–12): C-FNO 1.93e-2, 2.58e-2, 2.77e-2; FNO 2.49e-2, 3.30e-2, 4.60e-2.
- Training-set size (2V C-FNO, $t=1$): 300 samples $L^2$ = 1.84e-2; 240, 4.30e-2; 180, 1.27e-1 (Table 9).
- Cost: 2V, spectral / FNO per evaluation: $N=16$ 0.0162 s / 0.0150 s; 32: 0.0595 / 0.0304; 64: 0.271 / 0.130 (Table 8). 3V, $32^3$, Apple M4 CPU: one RK4 step 1.749 s; ten steps ($t=5.5\to6.5$) 17.59 s; one FNO endpoint inference 0.0143 s, speed-up ≈1231× (Table 14).
- DeepONet comparison (1V BGK): FNO error consistently smaller; total training 30.9 min (FNO) vs 86.3 min (DeepONet) (Table 4).

**Key results.**
- FNO flow maps reproduce BGK relaxation under 100-step autoregressive rollout without blow-up; error peaks early and then decays as the solution approaches the Maxwellian (simulated).
- Soft penalties cut the energy error of the 2V Boltzmann surrogate by one to two orders of magnitude, but it remains $10^{2}$–$10^{4}$ times larger than the spectral solver's (simulated).
- A plain FNO with $L^2$ error ≈1.5% in $f$ can carry an energy error of order 0.6–0.9 (simulated). With $\rho=T=1$ in 2V the exact energy is 1 (my calculation, not stated), so this is an order-unity relative error.
- Accuracy degrades sharply with less data (factor 7 in error from 300 to 180 samples) (simulated).
- The large speed-up (≈$10^3$) is for replacing ten solver steps by one endpoint evaluation in 3V; per step in 2V the gain is about 2× (simulated).

**Assumptions and limitations.**
- Authors: soft constraints "remain insufficient to guarantee strict conservation, particularly when the training data are of limited quality"; hard constraints or structure-preserving corrections are future work; only Maxwell molecules and homogeneous problems.
- Fixed time increment; no driving term; one trained map per time increment.
- Resolution invariance is asserted ("potentially enabling zero-shot super-resolution") but never tested; Table 8 gives timings only.
- The H-theorem is stated (2.11) but entropy is neither monitored nor constrained; positivity is not enforced.
- No quantum Boltzmann, Fermi–Dirac or Pauli-blocking content.
- It is not stated whether the 2V results at $t=4$ and $t=7$ come from separately trained endpoint maps or from iterating a unit-time map (my reading).
- Table 8's caption says "one time step" while the text says "one evaluation"; what the spectral time covers is ambiguous (my reading).
- Errors in Tables 5–13 are printed as plain norms $\|f^{\rm Num}-f^{\rm True}\|$; whether they are normalized is not stated.

**Use for this project.**
- Template for the flow-map ablation: one-step FNO on $[E, f]$ channels, rollout loss (3.15), optional penalties (3.13).
- Quantitative argument for hard projection: penalties left energy errors of $3\times10^{-3}$–$9\times10^{-2}$, against the project's $10^{-6}$ target.
- Pitfall: monitor the moments that set the optical response, not only $L^2(f)$; energy weights the tails where a small pointwise error is amplified.
- Budget warning: 300 trajectories were marginal in 2V; the WP1 plan of $2\times10^4$ trajectories is well above this.
- DeepONet is already shown weaker than FNO on a relaxation problem, consistent with listing it as a baseline.
- Consistent with the plan's statement that [49] uses soft constraints. The plan's phrase "Fourier neural operators exist for the classical Boltzmann collision operator [48, 49]" is accurate for [48] only; [49] learns the evolution operator.
- No code or data-availability statement.

---
### [50] Zhao et al., "Data-driven construction of a generalized kinetic collision operator from molecular dynamics" — Phys. Rev. Lett. 135, 185101 (2025); arXiv:2503.24208

**Source read:** arXiv v2 (4 Apr 2025) full text + TeX source, with appendices. The PRL reference is not printed in the arXiv PDF (not verified).
**Role in this project:** WP2 design reference: a learned collision operator whose conservation laws and H-theorem hold by construction rather than by projection or penalty.

**What it does.** A Landau-type (metriplectic) collision operator with a neural-network kernel is fitted to molecular-dynamics (MD) data for a one-component plasma. The tensor form of the kernel guarantees conservation, entropy production and a Maxwellian equilibrium for any network weights. An anisotropic kernel reproduces MD relaxation at coupling parameter $\Gamma\approx2.3$, where the Landau operator relaxes too fast.

**Model / method.** Spatially homogeneous equation for the velocity PDF $f(\boldsymbol v,t)$:

$$\frac{\partial f(\boldsymbol v,t)}{\partial t}=\nabla\cdot\int\boldsymbol\omega\left[f(\boldsymbol v')\nabla f(\boldsymbol v)-f(\boldsymbol v)\nabla' f(\boldsymbol v')\right]d\boldsymbol v' \qquad (1)$$

If the $3\times3$ kernel $\boldsymbol\omega$ is symmetric positive semi-definite and, for unitary $\mathcal U$,

$$\boldsymbol\omega(\mathcal U\boldsymbol v,\mathcal U\boldsymbol v')=\mathcal U\boldsymbol\omega(\boldsymbol v,\boldsymbol v')\mathcal U^T,\qquad \boldsymbol\omega(\boldsymbol v,\boldsymbol v')=\boldsymbol\omega(\boldsymbol v',\boldsymbol v),\qquad \boldsymbol\omega(\boldsymbol v,\boldsymbol v')(\boldsymbol v-\boldsymbol v')=\boldsymbol 0 \qquad (2)$$

then density, momentum and kinetic energy are conserved, $S=-\int f\log f\,d\boldsymbol v$ does not decrease, positivity is preserved and the equilibrium is Maxwellian (Prop. A.1). Proposed kernel ("CM2"):

$$\boldsymbol\omega=\boldsymbol{\mathcal P}\left(g_r^2\,\widetilde{\boldsymbol r}\widetilde{\boldsymbol r}^T+g_s^2\,\widetilde{\boldsymbol s}\widetilde{\boldsymbol s}^T\right)\boldsymbol{\mathcal P} \qquad (3)$$

with $\boldsymbol u=\boldsymbol v-\boldsymbol v'$, $\boldsymbol r=\boldsymbol v+\boldsymbol v'$, $\boldsymbol s=\boldsymbol u\times\boldsymbol r$, $\boldsymbol{\mathcal P}=\boldsymbol I-\boldsymbol u\boldsymbol u^T/|\boldsymbol u|^2$, $\widetilde{\boldsymbol r}=\boldsymbol{\mathcal P}\boldsymbol r/|\boldsymbol{\mathcal P}\boldsymbol r|$, $\widetilde{\boldsymbol s}=\boldsymbol s/|\boldsymbol s|$. The scalars $g_r,g_s$ are neural networks of the invariants $(u,r,s)$ (Appendix A; the main text prints $g_r(u,s,t)$, apparently a typo). "CM1" is the isotropic case $\boldsymbol\omega=g_I^2(u)\boldsymbol{\mathcal P}$; Landau is $g_I^2\propto1/u$. Training uses a weak form, so no 3D density estimate is needed:

$$L=\sum_{m=1}^{N_m}\sum_{l=1}^{N_t}\sum_{k=1}^{N_\psi}\left(\left.\frac{\partial f^{(l,m)}}{\partial t}\right|_{\rm MD}-\left.\frac{\partial f^{(l,m)}}{\partial t}\right|_{\rm c},\psi_k\right)^2 \qquad (4)$$

The MD term is $\frac{1}{N\delta t}\sum_i[\psi_k(\boldsymbol v_i^{l+1})-\psi_k(\boldsymbol v_i^{l})]$ (5); the kinetic term, integrated by parts, is a particle-pair sum estimated from $P$ random pairs per step (7).

**Numbers.**
- MD: $10^6$ singly charged ions, periodic $(100\ \text{Å})^3$ box, neutralizing background; $k_BT=100$ eV ($\Gamma=0.232$) and 10 eV ($\Gamma=2.321$); velocity Verlet, $dt=10^{-4}$ fs; 0.4 fs equilibration, 20 fs production (Appendix C).
- Training: $N_m=3$ initial distributions (uniform, bi-Maxwellian, double-well); $N_\psi=4$ Gaussian-type test functions; $P=10^5$ (Appendix B).
- Networks: fully connected, 6 layers, 10 neurons per hidden layer; Adam, $6\times10^5$ steps, learning rate 0.01 decayed by 0.99 every 2000 steps (Appendix B).
- No error norms, timings or speed-ups are reported; comparisons are plotted PDFs.

**Key results.**
- At 100 eV, CM1 and Landau both match MD. At 10 eV, Landau (even with a refitted time constant) and CM1 deviate by $t=0.6$ fs, while CM2 agrees for initial conditions not used in training (simulated, Figs. 2–4).
- The learned kernel has $g_r^2<g_s^2$ (Fig. 1), attributed to energy transfer from collective motion of the colliding pair with its environment (asserted interpretation).
- The operator is a 3D integral, against 5D for Boltzmann and 7D for Balescu–Lenard (derived).

**Assumptions and limitations.**
- Classical statistics, one species, homogeneous, Markovian; no quantum or degenerate-plasma content.
- The solver used for Eq. (1) with the learned kernel is not described (my reading).
- Three training trajectories; accuracy judged visually.

**Use for this project.**
- Design pattern: encode conservation in the algebraic form (projector $\boldsymbol{\mathcal P}$, symmetry under $\boldsymbol v\leftrightarrow\boldsymbol v'$) and learn only sign-definite scalars. The planned WP2 projection gives number and energy conservation but no H-theorem; this form gives both.
- No direct transfer: this is a Fokker–Planck-type operator with Maxwellian fixed point, whereas the project's e–e integral is Boltzmann-type with a Fermi–Dirac fixed point (my reading).
- The weak-form loss suits noisy or moment-level data.
- The plan lists the second author as "W. Burby"; the paper gives Joshua W. Burby.
- No data or code-availability statement.

---
### [51] Burles et al., "A neural operator closure for Landau damping in electrostatic plasma" — J. Plasma Phys. 92, E149 (2026); arXiv:2607.29364

**Source read:** arXiv v1 (31 Jul 2026) full text + TeX source, all sections. The PDF is marked "Under consideration for publication in J. Plasma Phys."; the journal reference is not printed (not verified).
**Role in this project:** Closest prior art for WP3 (a Fourier neural operator closure embedded in a time-stepping solver) and for the WP2 training recipe (unrolled rollouts through a differentiable solver); source of the stability and extrapolation evidence the methods paper will be compared against.

**What it does.** A three-moment electron fluid model (density, velocity, pressure) with Poisson's equation is closed by an FNO that returns the heat-flux gradient $\partial_xq$. The FNO is trained "online": gradients are backpropagated through 100-step rollouts of a differentiable fluid solver, and the loss compares fluid moments with Vlasov moments, never the heat flux itself. The closure takes a 50-step history of the resolved fields (non-Markovian). One network reproduces linear and nonlinear (trapping) Landau damping, runs stably for 5000 solver steps, interpolates between trained amplitudes and fails above the trained amplitude range.

**Model / method.** 1D electrostatic electrons, static ions; moments $n=\int f\,dv$, $u=\frac1n\int vf\,dv$, $p=m\int(v-u)^2f\,dv$, $q=m\int(v-u)^3f\,dv$ (2.4–2.7). Fluid system:

$$\frac{\partial n}{\partial t}=-\frac{\partial}{\partial x}(nu),\qquad \frac{\partial u}{\partial t}=-u\frac{\partial u}{\partial x}-\frac{1}{mn}\frac{\partial p}{\partial x}+\frac{q_{\mathrm e}}{m}E,\qquad \frac{\partial p}{\partial t}=-u\frac{\partial p}{\partial x}-3p\frac{\partial u}{\partial x}-\frac{\partial q}{\partial x} \qquad (2.8)\text{–}(2.10)$$

$$\frac{\partial E}{\partial x}=\frac{q_{\mathrm e}}{\epsilon_0}\big(n_0-n(x,t)\big) \qquad (2.11)$$

Learned closure replacing the last term of (2.10):

$$\frac{\partial p}{\partial t}=-u\frac{\partial p}{\partial x}-3p\frac{\partial u}{\partial x}-\mathrm{FNO}\big(\mathbf M;\log A_{\rm scale}\big) \qquad (3.2)$$

$$\mathbf m_t\equiv\big(n_t,u_t,p_t,E_t,(\partial_xq)_t\big),\qquad \mathbf M_t\equiv\big(\mathbf m_{t-(K-1)\Delta t},\dots,\mathbf m_t\big),\qquad (\partial_xq)=\mathrm{FNO}\big(\mathbf M_t;\log A_{\rm scale}\big) \qquad (3.6),(3.7)$$

$$A_{\rm scale}=\Big(\frac1N\sum_x\big(n(x)-1\big)^2\Big)^{1/2} \qquad (3.23)$$

Inputs are perturbations ($n-1$, $p-1$, $u$, $E$, $\partial_xq$) divided by $A_{\rm scale}$ (floored at $10^{-6}$), then passed through arcsinh; $\log A_{\rm scale}$ is appended as a constant channel; the output is multiplied by $A_{\rm scale}$. The network's own past outputs are part of its input. Non-locality in space enters through the Fourier layers; memory through the $K$-step window. The network outputs $\partial_xq$ rather than $q$ so that its spatial mean can be non-zero and heat the plasma; a periodic $q$ would force $\int\partial_xq\,dx=0$. Time stepping is SSPRK3 (3.3–3.5); the FNO is evaluated at every substage, a single rolling slot of the history buffer holds the provisional substage state, $E$ is recomputed from Poisson at each substage, and the pressure update at a substage uses the closure from the previous substage (3.8–3.22). Space: second-order central differences, spectral Poisson solve, all in PyTorch. Analytic baseline (Hammett–Perkins, non-local in space, memoryless, linear):

$$\tilde q_k=-n_0\sqrt{\tfrac{8}{\pi}}\,v_{\rm th}\frac{ik}{|k|}\tilde T_k \qquad (2.12)$$

**Numbers.**
- Ground truth: Eulerian Vlasov–Poisson code (splitting, fourth-order Van Leer scheme), collisionless; $N_x=256$, $N_v=401$, $v\in[-6,6]v_{\rm th}$; $k=0.35$, $L=2\pi/k$, periodic; $n_0(x)=1+A\cos kx$; $\Delta t_{\rm kin}=10^{-3}\omega_{\rm p}^{-1}$, output every $0.1\,\omega_{\rm p}^{-1}$; runs of 150 (linear) and 500 (nonlinear) $\omega_{\rm p}^{-1}$ (Table 1).
- Fluid solver: 128 points (kinetic fields downsampled), $\Delta t=0.1\,\omega_{\rm p}^{-1}$ (Table 3).
- FNO (`neuraloperator` library): 8 Fourier layers, 32 modes, 64 hidden channels, GELU, $5K+1=251$ input channels, 1 output; $K=50$ steps = $5\,\omega_{\rm p}^{-1}$ of memory (Table 2).
- Training: rollout 100 steps ($10\,\omega_{\rm p}^{-1}$); relative $L^2$ loss on all moments except $\partial_xq$ over all rollout steps; AdamW, learning rate $10^{-3}$ with cosine annealing, 1000 epochs, batch 50, no gradient clipping; one A100 GPU with 40 GB (Table 3).
- Linear case $A=10^{-3}$ (train $t<75$, test $t\ge75$): damping rate $\gamma_{\rm kin}=-3.43\times10^{-2}$; FNO $-3.43\times10^{-2}$ in training (relative error $1.56\times10^{-3}$), $-3.91\times10^{-2}$ in test against kinetic $-3.30\times10^{-2}$ (error 18%); Hammett–Perkins error 2.7% (train) and 6.9% (test) (Table 4).
- Hammett–Perkins coefficient calibrated to $\sqrt{4.8/\pi}\approx1.24$, 22% below the canonical $\sqrt{8/\pi}\approx1.60$; it still gives $\omega_{\rm r}=1.16$ against kinetic 1.22 (Sec. 3.6, 4.1).
- Time-averaged relative $L^2$ errors, train / test (Table 5). Linear: $E$ 6.8e-3 / 8.8e-2; $n$ 1.8e-6 / 2.0e-6; $u$ 5.7e-3 / 1.0e-1; $p$ 8.3e-6 / 6.5e-6; $\partial_xq$ 6.1e-2 / 9.7e-2. Nonlinear $A=10^{-1}$ (split at $t=250$): $E$ 5.3e-3 / 1.4e-1; $n$ 1.4e-4 / 2.9e-3; $u$ 8.6e-3 / 1.5e-1; $p$ 7.6e-4 / 1.3e-2; $\partial_xq$ 1.6e-1 / 1.01.
- Time after the training window at which the $E$ error exceeds 10% / 50% / 100%: linear 1.4 / 45.1 / 55.4 $\omega_{\rm p}^{-1}$; nonlinear 1.4 / 207 / 270 (Table 6).
- Multi-amplitude model: trained on $A=10^{(j/4-3)}$, $j=0,\dots,8$ (Eq. 4.1); tested on 24 amplitudes (my count of the listed grid) from $7\times10^{-4}$ to $3\times10^{-1}$, one-digit mantissas (Sec. 4.3).
- Wall-clock cost, speed-up over the Vlasov code, and parameter count are not stated.

**Key results.**
- Online-trained closure is stable over the full 500 $\omega_{\rm p}^{-1}$ nonlinear run (5000 steps, three FNO calls each) although trained only on 100-step rollouts; it reproduces saturation level and bounce oscillation (simulated).
- Beyond the trained time window the linear case loses phase gradually rather than blowing up (simulated).
- The same architecture trained offline (regression on kinetic $\partial_xq$, same data and hyperparameters) fails in both regimes when deployed; the authors restrict this to their matched configuration and note that Huang et al. obtained a working offline FNO with finer time sampling (simulated).
- Interpolation between trained amplitudes shows no train/test gap; extrapolation below the trained range ($7$–$9\times10^{-4}$) works; extrapolation above ($A=0.2$, 0.3) fails, with over-predicted field energy (simulated, Figs. 5, 6, 9).
- In the nonlinear regime the learned $\partial_xq$ differs from the kinetic one by order unity while the resolved moments stay accurate; the authors interpret it as an effective, scheme-dependent flux (simulated; interpretation asserted).
- Total energy is not enforced; its drift stays bounded but is one to two orders of magnitude larger than the Vlasov drift at large amplitude (simulated, Fig. 8).
- Ablation and HSIC analysis: the closure depends on the physical fields (mainly pressure), not only on its own output history, and reads the memory window at lags spaced by the field-energy oscillation period (simulated, Figs. 10–14).

**Assumptions and limitations.**
- Authors: 1D, electrostatic, single wavenumber, collisionless; memory cost of backpropagation through rollouts limits rollout length and batch size; no architecture comparison; the closure is tied to this SSPRK3 scheme and time step; mode truncation may fail for multi-mode problems.
- No uncertainty estimate, ensemble or fallback mechanism.
- No test of resolution transfer; trained and deployed on 128 points at one time step.
- Train/test split for single-amplitude models is chronological within one trajectory, so "test" is temporal extrapolation of the same run.
- The amplitude-conditioning input is computed from the current state, so nothing flags that the model has left its trained range (my reading).
- Text quotes $\gamma_{\rm HP}=-3.52\times10^{-2}$; Table 4 prints $-3.53$ and $-3.54\times10^{-2}$ (as printed).

**Use for this project.**
- Adopt the training recipe for WP2/WP3: unrolled multi-step training through the differentiable solver, with the loss on the quantities the solver needs. Here an a-priori-accurate closure was unstable in deployment while a trajectory-trained one ran 50 times longer than its training rollout.
- Adopt the input scaling: divide by a state-derived amplitude, apply arcsinh, pass the log-amplitude as a channel. The authors call this one of the most consequential choices; the project's intensity range (0.1 GW/cm² to 10 TW/cm²) spans five decades.
- Validation targets to mirror: a rate (damping rate), a saturation level, a divergence-time table, and a conserved-energy drift plot.
- Warning for the WP2 criterion "below 3% for intensities 10% beyond the training range": this closure failed just above its training envelope, so out-of-range tests must be reported, and an out-of-distribution detector (the planned ensemble) fills a gap this paper leaves open.
- Contrast: this closure needs memory because the truncated moments do not determine the flux. The project's operator acts on the full $f(E)$, so an instantaneous map should suffice; the flow-map ablation with input history (as in [53, 54]) is the analogue.
- Pitfall: a closure trained through one integrator inherits its time step and stage structure; WP3's multirate stepping should be fixed before training or the operator kept solver-independent by learning $Q[f]$.
- No speed-up is reported, so [51] cannot be cited for cost.
- Collisionless classical plasma; nothing on collision integrals or degenerate electrons.
- No data or code-availability statement.

---
### [52] Luo et al., "Resolution-independent machine learning heat flux closure for ICF plasmas" — PRX Intelligence 1, 013017 (2026); arXiv:2604.03439

**Source read:** arXiv v1 (3 Apr 2026) full text + TeX source (main text, 7 pages) and its figures. The Supplemental Material, which holds the FNO architecture, training strategy and the large-time-step tests, is not in the library and was not read. The journal reference is not printed in the arXiv PDF (not verified).
**Role in this project:** Venue precedent for the methods paper (an FNO kinetic closure published in PRX Intelligence) and the main evidence in the library that a closure trained on coarse data can be embedded in a finer-resolution implicit solver (WP2 resolution criterion, WP3 embedding).

**What it does.** An FNO is trained on particle-in-cell (PIC) data to map the electron temperature profile $T_e(x)$ to the heat-flux divergence $\partial_xq(x)$ in a collisional plasma with non-local transport. The learned operator replaces the heat-flux term in the electron energy equation, which is then solved implicitly. The closed equation reproduces PIC temperature evolution to better than 1%, beats the standard SNB non-local model, extrapolates in time, and gives the same accuracy whether trained on full-resolution data or on data 6 times coarser in space and 10–20 times coarser in time.

**Model / method.** Learned map and closed equation (constant density):

$$\partial_xq=\mathcal F_{\rm FNO}(T_e),\qquad (3/2)\,\partial T_e/\partial t+\mathcal F_{\rm FNO}(T_e)=0 \qquad (1)$$

The closure is instantaneous (no memory) and non-local in space; non-locality is carried by the global Fourier layers, at cost $O(n\log n)$ in the number of grid points. Eq. (1) is solved "implicitly using an iterative scheme" with two iterations per time step for both the FNO and the SNB closure; the scheme is cited, not written out. Normalization: lengths by the mean free path $\lambda_0=4\pi\varepsilon_0^2m_e^2v_{\rm th}^4/(Zn_ee^4\ln\Lambda)$, times by the collision time $1/\nu_0$, temperature by $T_{e0}$, heat flux by $n_{e0}T_{e0}v_{th0}$, $v_{th}=\sqrt{T_e/m_e}$. Error measure:

$$\mathcal L^2(\mathcal E)=\Big(\sum_{x,t}[\mathcal E_g-\mathcal F_{\mathcal E}]^2\Big)^{1/2}\Big/\Big(\sum_{x,t}\mathcal E_g^2\Big)^{1/2},\qquad \mathcal E\in\{T_e,\partial_xq\}$$

Physics check (Epperlein–Short test): a perturbation decays as $\delta T\propto e^{-\gamma t}$ with $\gamma=2k^2\kappa/3n_{e0}$, so $\kappa/\kappa_{\rm SH}=\gamma/\gamma_{\rm SH}$, with $\gamma_{\rm SH}/\nu_0=(2/3)(k\lambda_0)^2(128/\sqrt{2\pi})(Z+0.24)/(Z+4.2)$ (unnumbered, as printed). FNO layer count, modes, width, loss and optimizer: not stated in the main text.

**Numbers.**
- Ground truth: OSIRIS PIC, 1D, periodic, length $L=100\lambda_0$, duration $\tau=30/\nu_0$; cell size half a Debye length; $\Delta t/\Delta x\approx0.9$; $Z=1$ hydrogen; $6.4\times10^4$ electron and $6.4\times10^4$ ion macroparticles per cell; constant Coulomb logarithm. Data saved at $\mathbf{dx}\approx0.06$, $\mathbf{dt}=0.01$ (about 1670 points per profile, my calculation).
- Hot-spot case: $T_e(x)=1+\exp[-x^2/(\alpha\ell)^2]$, $\ell=8.44$; training $\alpha\in\{1,1.25,1.5,1.75,2\}$; test $\alpha\in\{1.125,1.375,1.625,1.875\}$; Knudsen number $\lambda_0/L_T$ up to ≈0.17.
- Epperlein–Short (ES) case: $T_e=1+\delta T_0\cos(\beta k_0x)$, $\delta T_0=4\%$, $k_0=2\pi/L$, $\beta\in\{1,4,7,10,13,16\}$, i.e. $k\lambda_0\in\{0.063,0.251,0.440,0.628,0.817,1.000\}$.
- Training window $t\in(0,20]$; $t\in(20,30]$ is temporal extrapolation.
- Resolution study: 15 models $\mathcal F_{(n,m)}$ trained on data downsampled by $n=dx/\mathbf{dx}\in\{1,3,6\}$ and $m=dt/\mathbf{dt}\in\{1,2,3,5,10\}$ (Table I); $m=20,50,100$ also tried for $n=6$ (Fig. 2i). All are deployed in the solver at $dx=\mathbf{dx}$, $dt=\mathbf{dt}$.
- SNB baseline: 180 energy groups, uniform energy grid up to 25 times the highest temperature.
- Hot-spot errors (read from Fig. 2): $\mathcal L^2(T_e)\approx0.004$ for $\mathcal F_{(1,1)}$ against ≈0.023–0.026 for SNB; $\mathcal L^2(\partial_xq)\approx0.12$–0.14 against ≈0.52–0.54. Across all 15 models $\mathcal L^2(T_e)$ lies between ≈0.0035 and ≈0.0067. For $n=6$ the $\alpha$-averaged $\mathcal L^2(T_e)$ stays ≈0.004 up to $m=20$ and rises to ≈0.006 ($m=50$) and ≈0.008 ($m=100$).
- Per-time error $\mathcal L^2_t(T_e)$ stays below 1% over nearly all of $t\in(0,30]$ for both $\mathcal F_{(1,1)}$ and $\mathcal F_{(6,10)}$, peaking at ≈1.1% near $t\approx2$–3 for $\alpha=1.125$ (Fig. 3b).
- Rollout length: 3000 time steps, 6000 closure iterations to $t=30$. Training profiles are 2.3% ($m=10$) or 1.1% ($m=20$) of that number (Summary).
- Cost: solving Eq. (1) to $t=30$ at $dt=\mathbf{dt}$ takes ≈800 min with SNB and ≈20 min with the FNO on a single CPU, a factor ≈40 (Summary). Training cost and PIC cost: not stated.

**Key results.**
- Embedded in the implicit solver, the FNO closure tracks PIC temperature and heat-flux divergence where SNB under-predicts $\partial_xq$ (simulated, Figs. 1, 2, 4).
- Accuracy is essentially independent of the resolution of the training data up to $n=6$, $m=20$; a model trained on coarse data runs in a solver six times finer in space and ten times finer in time, at time indices never seen in training (simulated).
- Temporal extrapolation from $t\le20$ to $t=30$ holds at the 1% level (simulated, Fig. 3a,b).
- Generalization: unseen intermediate widths; sub-Gaussian, super-Gaussian, Lorentzian and sech initial profiles; Gaussian widths outside the trained range ($\alpha=0.75$ and 4) all agree with PIC at $t=5$ and $t=25$ (simulated, Fig. 3c–h; shown as profiles, no error numbers).
- ES test: all 15 models reproduce the PIC decay of $\delta T$; the inferred $\kappa/\kappa_{\rm SH}$ versus $k\lambda_0$ follows published Vlasov–Fokker–Planck results (OSHUN, KIPP) and a published fit, while SNB falls below at $k\lambda_0\gtrsim0.25$ (simulated, Fig. 5).
- A model trained on the hot-spot family cannot predict the ES case (stated by the authors as a limitation).

**Assumptions and limitations.**
- Authors: predictive capability is limited when initial conditions deviate substantially from the training set; robust generalization across widely separated conditions "remains challenging".
- 1D, constant density, $Z=1$, periodic boundaries, constant Coulomb logarithm, no magnetic field; closure depends on $T_e(x)$ only.
- The speed-up is against the SNB reduced model, not against the PIC code (my reading).
- No conservation constraint, uncertainty estimate or fallback is described in the main text. Writing the output as $\partial_xq$ rather than $q$ does not guarantee $\int\partial_xq\,dx=0$ on the periodic domain, so energy conservation is not built in (my reading).
- Larger solver steps ($10\,\mathbf{dt}$, $100\,\mathbf{dt}$) are said to work with more iterations per step; evidence is in the Supplemental Material (not read).

**Use for this project.**
- What made this publishable in PRX Intelligence, as a checklist for Paper B: (i) the incumbent physics model (SNB) run in the same solver as baseline; (ii) accuracy judged after embedding in the solver, not only a priori; (iii) a physics observable ($\kappa/\kappa_{\rm SH}$ versus $k\lambda_0$) compared with independent kinetic codes; (iv) interpolation, temporal extrapolation, changed functional shapes and out-of-range parameters; (v) a systematic train-resolution versus deploy-resolution table; (vi) an explicit failure case; (vii) a wall-clock comparison. The main text tests one learned model class only; no MLP, CNN or DeepONet baseline appears there (Supplemental Material not read).
- WP2 criterion "train at 256, evaluate at 512 with less than twofold error growth" is less demanding than what is shown here (6× in space with no measurable error growth), but here the operator acts on a smooth profile; the Fermi edge makes the project's case harder (my reading).
- Embedding pattern for WP3: treat the operator as the right-hand side inside an implicit iteration; stability came with two iterations per step over 3000 steps.
- The plan's statement that PRX Intelligence "published a Fourier-operator plasma closure [52]" matches the citation supplied; the arXiv version calls itself a Letter.
- Collisional classical plasma; nothing on collision operators, degenerate electrons or optical driving.
- Data availability (quoted): "The data are not publicly available. The data are available from the authors upon reasonable request." An open dataset and code would go beyond this precedent.

---
### [53] Shah and Cangi, "Accelerating electron dynamics simulations through machine learned time propagators" — ICML 2024 AI for Science Workshop; arXiv:2407.09628

**Source read:** arXiv v2 (25 Jul 2024) full text + TeX source (workshop paper, 9 pages).
**Role in this project:** Template, with [49], for the WP2 ablation (autoregressive flow-map operator) and the first of two papers showing what the Shah–Cangi group has demonstrated.

**What it does.** A Fourier neural operator replaces the time propagator of real-time TDDFT. It advances the electron density $n(x,t)$ directly, without Kohn–Sham orbitals, from a window of previous density snapshots, with a time step 8 times the reference solver's. The system is a 1D two-electron soft-Coulomb diatomic molecule under one fixed sinusoidal laser field. Over a 41-step rollout the FNO is about 4 times more accurate and 2 times faster than the Octopus code at the same coarse step.

**Model / method.** Reference: time-dependent Kohn–Sham equations, adiabatic LDA, Crank–Nicolson in Octopus, zero boundary conditions, with

$$v_{\rm ion}=-\frac{Z_1}{\sqrt{(x-\frac d2)^2+a^2}}-\frac{Z_2}{\sqrt{(x+\frac d2)^2+a^2}},\qquad v_{\rm las}(t)=A\sin\omega t \qquad (10),(11)$$

Autoregressive flow map on a sliding window of $T_{\rm in}$ snapshots:

$$\mathbf N_t=[\mathbf n_{t-T_{\rm in}+1},\dots,\mathbf n_t],\qquad \hat{\mathbf n}_{t+1}=\mathcal G_\theta(\mathbf N_t),\qquad \mathbf N_{t+1}=[\mathbf n_{t-T_{\rm in}+2},\dots,\mathbf n_t,\hat{\mathbf n}_{t+1}] \qquad (18)\text{–}(20)$$

Loss: MSE over all predicted steps and systems (21) plus a soft particle-number penalty

$$\mathcal L_{\rm norm}(\theta)=\lambda\Big(\sum_i\hat n^{(d)}_{t,i}\Delta x-2\Big)^2 \qquad (23)$$

The laser is not an input. The first $T_{\rm in}$ snapshots come from the reference solver.

**Numbers.**
- 729 molecules with varied $Z_1$, $Z_2$, $d$ ($a=1.0$); 600 train / 10 validation / 100 test (Sec. 3).
- Domain $[-6,6]$ a.u., $\Delta x=0.05$ a.u. (241 points); duration 10 a.u. (0.241 fs); reference step 0.025 a.u. (0.6 as); ML step 0.2 a.u. (4.83 as); 51 slices, $T_{\rm in}=10$, 41 predicted.
- Laser: $A=0.75$ a.u. ($1.97\times10^{16}$ W cm$^{-2}$), $\omega=2.0$ Ha (22.78 nm).
- Network: 3 Fourier layers, 16 modes, width 128 (baseline), 2000 epochs, Adam. $\lambda$ and hardware: not stated.
- Table 1 (errors against the fine-step reference; time per slice): baseline FNO MAE $5.118\times10^{-3}$, MSE $2.56\times10^{-4}$, 1.75 ms; Octopus at 0.2 a.u. MAE $19.6\times10^{-3}$, MSE $24.7\times10^{-4}$, 4.00 ms. MAE falls from $8.60$ to $5.65\times10^{-3}$ as $T_{\rm in}$ goes from 5 to 20 (width 64).
- Super-resolution (trained at 241 points, evaluated at 481): MAE $7.335\times10^{-3}$ at 1.75 ms, against 9 ms for the solver. Half-step time offset: MAE $5.779\times10^{-3}$.
- Particle number $1.99992\pm0.00039$ (target 2), deviation growing with time (Fig. 4).
- Dipole moment: mean absolute percentage error 8.34%, accumulating over the rollout (Fig. 5).
- Reverse-time rollout: MAE 0.0751.

**Key results.**
- A density-only flow map is stable over 41 steps at an 8× larger step and beats the coarse solver in accuracy and speed (simulated).
- Zero-shot doubling of spatial resolution raises MAE by about 40% at unchanged cost (simulated).
- Particle number is held to about $2\times10^{-4}$ relative by the penalty but drifts; time-reversal symmetry is broken (simulated).

**Assumptions and limitations.**
- 1D, two electrons, ALDA, dipole approximation, one laser setting, 0.24 fs horizon.
- Speed-up (about 2.3×) is against Octopus at the same coarse step, not the fine reference.
- Not unitary; no continuity or positivity constraint; authors call for stronger conservation and built-in time-reversal symmetry.

**Use for this project.**
- Recipe for the flow-map ablation: sliding window of past states, rollout MSE, soft conservation penalty, zero-shot grid change as a reported test.
- Benchmark: conservation at $10^{-4}$ relative with drift, against the project's $10^{-6}$ by construction.
- Pitfall: without the drive as input the map is tied to one pulse; chirped and double pulses need the field as input ([54]) or an operator for $Q[f]$ with analytic drive, as WP2 plans.
- Pitfall: the integrated observable (dipole) degrades faster than the state error; report plasma frequency and damping over the rollout.
- Competitor status: model molecules only; nothing on solids, kinetic equations or Maxwell coupling.
- No data-availability statement.

---
### [54] Shah and Cangi, "Machine learning time propagators for time-dependent density functional theory simulations" — Mach. Learn.: Sci. Technol. 7, 035019 (2026); arXiv:2508.16554

**Source read:** arXiv v3 (1 Dec 2025) full text + TeX source, with Supplementary Information (SI). The journal reference is not printed in the arXiv PDF (not verified).
**Role in this project:** Journal-length successor to [53]; fixes the state of the art of the Shah–Cangi group (named in the plan as possible competitors) and gives the design of a drive-conditioned flow map for the WP2 ablation.

**What it does.** The autoregressive FNO density propagator of [53] is extended to a range of laser wavelengths and intensities by feeding the laser history to the network with the density history. On 1D two-electron diatomic molecules the best model has about 10 times lower density MSE and is 8 times faster per step than Octopus at the same coarse time step.

**Model / method.** Reference: time-dependent Kohn–Sham equations, adiabatic LDA, soft-Coulomb potentials ($a=1$), Crank–Nicolson in Octopus, $v_{\rm las}(t)=A\sin\omega t$. Flow map:

$$\mathbf N_t=[\mathbf n_{t-T_{\rm in}+1},\dots,\mathbf n_t],\qquad \mathbf V_{l,t}=[v_{l,t-T_{\rm in,laser}+1},\dots,v_{l,t}],\qquad \hat{\mathbf n}_{t+1}=\mathcal G_\theta(\mathbf N_t,\mathbf V_{l,t}) \qquad (24)\text{–}(26)$$

Default input encoding, with $\mathbf g$ a normalized coordinate grid and $W_d$, $W_l$ learned linear lifts:

$$P_{\text{sep-add}}(\mathbf N_t,\mathbf V_{l,t},\mathbf g)=W_d[\mathbf N_t,\mathbf g]+W_l[\mathbf V_{l,t}] \qquad (32)$$

$$\mathcal L_{\rm total}=\mathcal L_{\rm data}+\lambda\mathcal L_{\rm int},\qquad \mathcal L_{\rm int}=\frac{1}{|\mathcal D|}\sum_{d\in\mathcal D}\frac{1}{|\mathcal T|}\sum_{t\in\mathcal T}\Big(\sum_i\hat n^{(d)}_{t,i}\Delta x-N^{(d)}_{\rm ref}\Big)^2 \qquad (42),(44)$$

$\mathcal L_{\rm data}$ is the rollout MSE (43); $N_{\rm ref}=2$. Training is unrolled: windows start from reference data and the model runs on its own predictions for the whole rollout (SI Algorithm 1). The density enters as $\log_{10}n$ scaled to $[-1,1]$.

**Numbers.**
- Data: 2048 simulations; intensity $10^{12}$–$10^{14}$ W cm$^{-2}$, wavelength 400–750 nm; after filtering 800 train / 150 validation / 200 test (Sec. 3.1, SI).
- Grid $[-9,9]$ a.u., $\Delta x=0.05$ a.u., 361 points; $t\in[0,5]$ fs; reference step 0.01 fs; ML step 0.1 fs; $T_{\rm in}=10$ inputs, 41 predicted steps.
- Network: 3 Fourier layers, width 128, 32 modes; AdamW, 800 epochs; $\lambda=0.1$ (SI).
- Encodings (Table 1), density MAE ($10^{-3}$) / dipole MAPE: density only 1.97 / 27.6%; concatenation 0.899 / 18.6%; separate + add 0.872 / 9.35%; space-time convolution 1.48 / 27.0%.
- Against Octopus at 0.1 fs, both judged against the 0.01 fs reference (Table 2): MAE $0.872\times10^{-3}$ vs $6.10\times10^{-3}$; MSE $3.20\times10^{-5}$ vs $31.2\times10^{-5}$; dipole SMAPE 6.08% vs 51.1%; particle number $1.9995\pm0.0014$ vs $2.0000\pm2\times10^{-7}$; time per step 3.07 ms vs 23.94 ms. Table 1 lists 1.35 ms for the same encoding (not reconciled).
- Table 3: 721-point grid MAE $0.942\times10^{-3}$; rollout to 10 fs MAE $0.896\times10^{-3}$, MSE $2.53\times10^{-5}$ (text prints $2.56\times10^{-5}$, as printed); time-reversed rollout with reversed field MAE $1.053\times10^{-3}$.

**Key results.**
- Adding the laser history roughly halves the density error and cuts the dipole percentage error threefold (simulated).
- Rollout to twice the training horizon (about 90 steps, my calculation) shows no error growth; a twice-finer grid raises MAE by 8% (simulated).
- The penalty holds particle number to about $10^{-3}$ absolute; the coarse solver holds it to $10^{-7}$ (simulated).
- Backward rollouts retrace the forward trajectory; [53] had found time-reversal symmetry broken for the density-only model (simulated).

**Assumptions and limitations.**
- Authors: 1D two-electron molecules, ALDA, sinusoidal driving; not unitary; only an integral penalty; 3D feasibility "needs yet to be shown"; continuity and non-negativity are future work.
- Speed-up is against the solver at the same coarse step. The first 10 snapshots (about 1 fs) come from the reference solver.

**Use for this project.**
- For the flow-map ablation adopt separate lifts for state and drive (Eq. 32), unrolled training, and a log-transformed state.
- The authors propose a predictor–corrector hybrid (FNO proposes large steps, a conservative integrator corrects drift), close to the plan's fallback to the exact integral.
- Competitor assessment: variable laser parameters are now handled, but the physics is orbital-free propagation of the density of 1D model molecules, with no collision integral, solid or Maxwell coupling. Stated future work (3D, many electrons, beyond ALDA) does not point toward Boltzmann kinetics.
- Data availability (quoted): "The data that support the findings of this study are available at RODARE record 3995, url: https://rodare.hzdr.de/record/3995."

---

## Neural operators in photonics; machine learning for ultrafast metrology
### [55] Gu et al., "NeurOLight: A Physics-Agnostic Neural Operator Enabling Parametric Photonic Device Simulation" — NeurIPS 2022 (Adv. Neural Inf. Process. Syst. 35, 14623); arXiv:2209.10098

**Source read:** arXiv v1 full text + TeX source (main text, appendices A1–A5, tables); figures checked in the PDF. The NeurIPS camera-ready was not compared.
**Role in this project:** Supports gap 3 of the plan (neural operators exist for Maxwell solvers) and gives WP3 a measure of how accurate a learned Maxwell solver is.

**What it does.** Learns the solution operator of the 2D frequency-domain (time-harmonic, steady-state) Maxwell equations for linear, isotropic silicon multimode-interference (MMI) devices. Input: permittivity map, wavelength, domain size and input-port source; output: the complex field $\mathbf{H}_y$. A modified Fourier neural operator is trained on finite-difference frequency-domain (FDFD) data only; "physics-agnostic" means no Maxwell residual in the loss. Headline: normalised mean absolute error 0.12 (tunable MMI) and 0.39 (etched MMI), 140–200× faster than a CPU FDFD solve.

**Model / method.**
PDE whose solution map is learned (Eq. 1; curl-of-curl form with absorbing boundaries; $\mathbf{J}_e,\mathbf{J}_m$ current sources, $\vec\epsilon_r$ relative permittivity):

$$\big((\mu_0^{-1}\nabla\times\nabla\times)-\omega^2\epsilon_0\vec{\epsilon}_r(\vec{r})\big)\mathbf{E}(\vec{r})=j\omega\mathbf{J}_e(\vec{r}),\qquad \big(\nabla\times(\vec{\epsilon}_r^{-1}(\vec{r})\nabla\times)-\omega^2\mu_0\epsilon_0\big)\mathbf{H}(\vec{r})=j\omega \mathbf{J}_m(\vec{r})$$

Layer (Eqs. 3–4): $v_{k+1}=\texttt{FFN}(\mathcal{K}v_k)+v_k$, where the $C$ latent channels are split in two halves that receive 1D Fourier-space multiplications along $z$ and $x$ ("cross-shaped"; $R^h, R^v$ learned complex weights), and FFN contains a 3×3 depth-wise convolution:

$$(\mathcal{K}^hv^h_k)(\vec{r})=\mathcal{F}_z^{-1}\big(R^h(z)\cdot\mathcal{F}_z(v_k^h(\vec{r}))\big),\qquad (\mathcal{K}^vv^v_k)(\vec{r})=\mathcal{F}_x^{-1}\big(R^v(x)\cdot\mathcal{F}_x(v_k^v(\vec{r}))\big)$$

Loss (Sec. 3.3): $\mathcal{L}=\|\Psi_\theta(\mathcal{E}(a))-\Psi^{\ast}(a)\|_1/\|\Psi^{\ast}(a)\|_1$.

Encoding of parameters and source (Sec. 3.2):
- Every device is resampled to one $M\times N$ grid (80 × 384 pixels, Fig. 3); physical size enters through the steps $\Delta l_x=l_x/M$, $\Delta l_z=l_z/N$.
- Domain, permittivity and frequency are merged into two "wave prior" images, $\mathcal{P}_z=e^{j\frac{2\pi\sqrt{\epsilon_r}}{\lambda}\mathbf{1}z^{T}\Delta l_z}$ and $\mathcal{P}_x=e^{j\frac{2\pi\sqrt{\epsilon_r}}{\lambda}x\mathbf{1}^{T}\Delta l_x}$: plane waves with the local in-material wavelength.
- Source $\mathbf{H}_y^J$: the known field in the input waveguide, zero elsewhere.
- Network input: $a^\dagger=(\epsilon_r,\mathbf{H}_y^J,\mathcal{P}_x,\mathcal{P}_z)\in\mathbb{C}^{4\times M\times N}$.
- Because the response is linear in the source (Eq. 5), training uses single-port simulations mixed on the fly with random complex weights.

**Numbers.**
- Devices (Table A4): length 20–30 µm, width 5.5–7 µm, λ = 1.53–1.565 µm, $\epsilon_r$ = 11.9–12.3 (tunable pads) or {2.07, 12.11} (etched); TM polarisation.
- Data and network (Sec. A2–A4): 7,680 tunable and 15,360 etched examples from the FDFD code `angler`; 12 layers (16 for etched), $C$ = 64, Fourier modes $(k_z,k_x)$ = (70, 40), 1.58 M parameters (2.11 M etched).
- Test error (Table 1): tunable MMI — UNet 0.801, FNO-2d 0.244, F-FNO 0.292, NeurOLight 0.122. Etched MMI — 0.792, 0.648, 0.525, 0.387.
- Speed: 12.1 ms per inference on one Quadro RTX 6000 (Table 2); FDFD nearly 4 s for 5 × 20 µm² on 8 CPU cores (Sec. 1); 140–200× (Fig. 1b); a 1550–1565 nm sweep in 2 nm steps takes 150 ms, quoted as 450× (Sec. 4.4).

**Key results.**
- The encoding matters: raw $(\epsilon_r,\lambda,\text{domain})$ inputs generalise poorly; permittivity plus wave prior is best (simulated, Fig. 10b).
- Many more Fourier modes are needed than in fluid benchmarks, yet keeping all modes overfits; removing the local depth-wise convolution raises the error from 0.122 to 0.144 (simulated, Fig. 10c, Table 2).
- Unseen wavelengths inside the training band work. A 3-port model applied to 4- and 5-port devices gives "unsatisfying" error until fine-tuned on 3.7 K and 4.6 K new examples (simulated, Sec. 4.4).

**Assumptions and limitations.**
- Linear, isotropic, 2D, one frequency at a time, one polarisation; device level only (authors).
- Field errors of 12–39 % are far from solver accuracy; nothing enforces boundary conditions or energy balance (my reading).
- The speed-up is measured against a 2D sparse linear solve on CPU (my reading).

**Use for this project.**
- Division of labour: supports the project's choice. A learned Maxwell solver reaches tens-of-percent field error and degrades off-distribution, whereas WP3 needs transfer-matrix-level agreement in the linear limit. In 1D a conventional FDTD is cheap and converges under grid refinement; the cost sits in the collision integral. The paper's only case for learning Maxwell is speed in 2D/3D, which WP3 defers (my reading).
- The authors state that building a Maxwell residual loss "is no easier than re-implementing a numerical solver".
- Transferable to WP2: encode parameters as physically meaningful fields, not raw scalars; pair a truncated Fourier kernel with a local convolution; use exact linearity for augmentation where it holds.
- Citation: the plan lists the arXiv preprint; the paper appeared at NeurIPS 2022. Code: https://github.com/JeremieMelo/NeurOLight.

---
### [56] Ma et al., "PIC²O-Sim: A Physics-Inspired Causality-Aware Dynamic Convolutional Neural Operator for Ultra-Fast Photonic Device FDTD Simulation" — APL Photonics 10, 036104 (2025); arXiv:2406.17810

**Source read:** arXiv v1 full text + TeX source (main text, appendix A.1–A.6, tables). The APL Photonics version was not available; all numbers are from v1.
**Role in this project:** Context for WP3 (the closest published learned time-domain Maxwell surrogate) and for the WP2 flow-map ablation: a documented case of autoregressive error accumulation and its mitigation.

**What it does.** Replaces 2D FDTD time-marching for linear dielectric photonic devices by a convolutional neural operator. One call maps the last $T_{in}$ field frames, the permittivity map and the next $T$ source frames to the next $T$ frames of $E_z$; calls are chained autoregressively. The receptive field is limited to the light cone and the convolution weights depend on local permittivity. Headline: 160 fs rollouts with 5–9 % normalised L2 error, 308–632× faster than the FDTD code MEEP.

**Model / method.**
- Learned map (Sec. 3.2): $\Psi_\theta:(E_{in}\in\mathbb{R}^{T_{in}\times H\times W},\ \epsilon_r,\ J_{1:T})\mapsto E_{1:T}$. Reference solver: the standard Yee update for TM$^z$ fields with a scalar $\epsilon[m,n]$ per cell (Eq. 1, App. A.1).
- Permittivity-adaptive convolution (Sec. 3.3, unnumbered); $v^l$ latent field at layer $l$, $\mathbf{W}$ learned static filter, $\mathbf{f}$ features of $1/\epsilon_r$ from a small convolutional encoder, $\Omega(i)$ a dilated $K\times K$ window:

$$z^{l}(i)=\sum_{j \in \Omega(i)} \big(\mathcal{K}\left(\mathbf{f}_i, \mathbf{f}_j\right)(\epsilon_r)\, \mathbf{W}^T(j)\big) \cdot v^l(j),\qquad \mathcal{K}\left(\mathbf{f}_i, \mathbf{f}_j\right)(\epsilon_r)=\exp \left(-\tfrac{1}{2}\left(\mathbf{f}_i-\mathbf{f}_j\right)^{\top}\left(\mathbf{f}_i-\mathbf{f}_j\right)\right)$$

- Receptive field (App. A.3); $\epsilon_r^{wg}$ waveguide permittivity, $L$ layers, dilation $s$:

$$R \approx \frac{Tc_0}{\sqrt{\epsilon_r^{wg}}}+30,\qquad K=\Big[\frac{R-5}{sL}\Big]$$

  (as printed; looks inconsistent because the first term is a length while 30 and 5 are pixel counts, so it must first be divided by the grid step). The extra 30 pixels are empirical.
- Source: $J_{1:T}$ enters as $T$ frames that are zero except on the source line. Loss and metric (Sec. 4.1): $\texttt{N-L2Norm}=\frac{1}{T}\sum_{t=1}^{T}\|\Psi_{\theta}(E_{in},\epsilon_r,J_{1:T})[t]-E^{\ast}[t]\|_2/\|E^{\ast}[t]\|_2$.
- Rollout (Sec. 3.4): the horizon is cut into $N$ stages of $T$ frames. Each stage has its own sub-model, trained on the imperfect output of the previous stage, and also receives that stage's hidden state.
- The paper never names the CFL condition. Its causality argument is the light cone: a field point depends only on fields within $c_0T/\sqrt{\epsilon_r}$, and a frame depends only on the immediately preceding state.

**Numbers.**
- Data (Sec. 4.1, App. A.2, Table 5): MEEP videos resampled to 140 nm pixels and 1 fs frames; FDTD itself typically steps at 0.167 fs with about λ/15 cells (Sec. 3.1). Training/test devices: MMI 20/5, micro-ring 6/5, metaline 32/8. Permittivity (cladding, core) {2.07, 12.11}, ring {1, 6}; source band 1.4–1.65 µm.
- Model (App. A.5): $T_{in}$ = 10, $T$ = 80 per stage, 2 stages, 8 layers, 72 channels, kernel 5 (ring 7), dilation 4; 2.4 M parameters (ring 4.4 M).
- Test error over 160 frames (Table 1): MMI 0.052, ring 0.085, metaline 0.086. Full-mode FNO (340 M parameters): 0.122, 0.423, 0.173. NeurOLight [55]: 0.140, 0.151, 0.185.
- Speed: 15,701 / 1,906 / 7,348 frames per second on an A6000 GPU (Table 1); 308–632× against MEEP on a 64-core CPU (Fig. 1a, App. A.4); MEEP needs about 1 min for 28 × 10.5 µm² over 1.3 ps (Sec. 1).
- Rollout ablations, MMI, plain-CNN backbone (Tables 2–3): stages × frames 1×160 → 0.0638; 2×80 → 0.0513; 4×40 → 0.0538; 8×20 → 0.0871. One network applied twice: single-call error 0.0213, rollout 0.0730. Two independent networks with hidden-state passing: 0.0486.

**Key results.**
- PIC²O-Sim has the lowest test error on all three device types. The margin over the best baseline is small on the ring (0.085 against 0.088, plain CNN) and the metaline (0.086 against 0.089, F-FNO) (simulated, Table 1).
- Rollout error, not single-call error, sets the accuracy; many short stages are worse than two long ones because iterations multiply (simulated, Tables 2–3, Fig. 10).
- Two input frames should suffice in theory; 8–10 were needed (simulated, Fig. 8).

**Assumptions and limitations.**
- "Our framework still observes rapidly accumulated errors with large-iteration rollout" (authors, Sec. 5).
- Linear, isotropic, non-dispersive, time-independent permittivity; 2D; one polarisation. No material model enters (my reading).
- Longest tested horizon is 160 fs. The abstract quotes 51.2 % lower error than baselines, the conclusion 49.1 % (my reading).

**Use for this project.**
- Division of labour: supports learning only the closure. A learned field propagator drifts to 5–9 % within 160 fs, while WP3 needs 2 ps at percent level. The authors grant that schemes embedding the FDTD update itself "have a theoretical guarantee on the solved fields" at fine resolution but limited speed-up; in 1D that cost is negligible. Time-varying and dispersive media are not covered.
- For the WP2 flow-map ablation and unrolled training: train later steps on the model's own earlier output; expect rollout error several times the one-step error.
- Maxwell is local in space and time, so a per-cell collision operator acting on the current state is the matching structure for the closure (my reading).
- Citation: the plan lists the arXiv preprint; published in APL Photonics (2025). Code: https://github.com/ScopeX-ASU/PIC2O-Sim.

---
### [57] Blackwood and Hawley (authors taken from related records), "Physics-informed Fourier neural operator for nanophotonic scattering in nonlinear media" — poster, APS Global Physics Summit 2026

**Source read:** abstracts only, and NOT from the APS page (https://meetings-archive.aps.org/smt/2026/mar-h00-01/159/), which could not be retrieved. Everything below comes from two Belmont University repository records by B. L. Blackwood and S. Hawley: a 2025 poster with nearly the same title (repository.belmont.edu/surs/320) and a 2026 talk, "Fast Prediction of Photonic Scattering via a Modified Fourier Neural Operator" (repository.belmont.edu/spark_presentations/622). No full text exists; that these records describe the APS poster is an assumption.
**Role in this project:** Context for gap 3 of the plan: with [55, 56], a third neural operator used as a Maxwell solver.

**What it does.** Per the Belmont abstracts, a physics-informed Fourier neural operator is trained on FDTD simulations of Maxwell's equations for scattering problems containing nonlinear materials. It maps field configurations forward in time, across several equivalent FDTD steps per call. The 2025 abstract also mentions a model that predicts optical spectra from initial conditions.

**Numbers.** (2026 talk abstract only)
- Up to 100× reduction in computation time against high-resolution FDTD.
- Normalised mean square error "as low as 0.24", "outperforming interpolation in all cases". Errors were measured on unseen permittivity distributions and source parameters, against training-set size and time-step size (asserted; no figures).
- Material model, type of nonlinearity, wavelengths, geometry, grid, network and training-set size: not stated.

**Assumptions and limitations.**
- Nothing to check; the "physics-informed" term is not described. A best-case error of 0.24 is far from solver accuracy (my reading).

**Use for this project.**
- Cite only as evidence that others apply Fourier operators to nonlinear nanophotonic scattering; take no equations or parameters from it.
- It learns the field propagator, not a material closure; the abstracts mention neither electron kinetics nor ENZ media (my reading).
- Check the APS page and author list before submission.

---
### [58] Zahavy et al., "Deep learning reconstruction of ultrashort pulses" — Optica 5, 666 (2018); arXiv:1803.06024

**Source read:** arXiv v1 full text + TeX source + figures. This is NOT the Optica article: v1 is the four-page NIPS 2017 workshop paper ("Deep Learning Reconstruction of Ultra-Short Pulses"; authors Zahavy, Dikopoltsev, Cohen, Mannor, Segev). The Optica version adds D. Moss and G. I. Haham and was not available. Verify against Optica before quoting numbers.
**Role in this project:** WP5 retrieval: precedent for neural-network FROG retrieval and for refining a retrieval by back-propagating through a differentiable model of the measurement.

**What it does.** Trains a convolutional network to invert a second-harmonic-generation (SHG) FROG trace into the complex pulse field and compares it with two iterative algorithms, PCGPA and ptychographic FROG. In simulation the network trained on noisy traces is more accurate than both below 20 dB signal-to-noise ratio. On one measured Ti:sapphire pulse it gives the lowest trace error, but only after a label-free refinement on the measured trace.

**Model / method.**
Measurement model (Sec. "Problem formulation", unnumbered); $E(t)$ unknown complex pulse, $\tau_j$ delays, $\omega_i$ spectrometer frequencies:

$$I_{measured}(\omega_i,\tau_j)=|\mathcal{F}\{E(t)E(t-\tau_j)\}(\omega_i)|^2$$

Training objective (Eq. 1). "FROG Net" is a fixed, differentiable implementation of the model above; CNN$(I;w)$ is the trainable network; $\lambda$ weights the trace-consistency term:

$$w^* = \underset{w}{\arg\min}\,\{\mathrm{loss}(I,E)\} = \underset{w}{\arg\min}\,\{ \| \mathrm{CNN}(I;w) - E \| _1 + \lambda \| \mathrm{FROG\ Net}(\mathrm{CNN}(I;w)) - I \| _1 \}$$

- Supervised stage: first term, on simulated (trace, pulse) pairs. "Unsupervised" stage: second term only, on the measured trace, with no pulse label.
- Network: of three architectures tried, only "Multires" (several filter sizes per convolution layer) is reported; sizes are not given.
- Ambiguities: an SHG-FROG trace fixes the pulse up to a constant phase, a time shift, and time reversal with conjugation. All equivalent labels are mapped to one representative before training, because ambiguous labels destabilise it.

**Numbers.**
- Training set: 60,000 simulated pulses, 10,000 for testing; random spectral phase (random complex vector times a Lorentzian envelope) on a Gaussian spectrum (Sec. "Sim2Sim").
- Noise: additive white Gaussian noise on the trace, 0–30 dB drawn uniformly during training.
- Error against noise (Fig. 3 left, read off the plot): at −10 dB, noise-trained network ≈ 0.04, both iterative methods ≈ 0.07, network trained without noise ≈ 0.24. At 0 dB: ≈ 0.02 against ≈ 0.05. The curves cross near 20 dB. At 30 dB: network ≈ 0.015, iterative ≈ 0.01.
- One 10 dB example (Fig. 3 right): pulse error PCGPA 0.0334, ptychography 0.0326, network 0.0145.
- Measured pulse, mean absolute trace error (Fig. 4): best refined networks 0.00498 (one refined variant, 0.00586, is worse than the iterative methods); PCGPA 0.00532; ptychography 0.00532; supervised-only networks 0.0060–0.0098.
- Trace grid, delay step, pulse duration, wavelength, network size, $\lambda$: not stated. The plotted delay axis spans about ±80 fs.

**Key results.**
- Noise-aware training makes the network more robust than PCGPA and ptychography at low signal-to-noise ratio; at high signal-to-noise ratio the iterative algorithms are more accurate (simulated; the second point is my reading of Fig. 3).
- A network trained only on simulations is worse than the iterative algorithms on real data; the trace-consistency refinement closes the gap, and the authors call it "crucial" (measured, Fig. 4).

**Assumptions and limitations.**
- SHG-FROG only: one unknown pulse and an instantaneous gate that is the pulse itself. No blind or double-blind case and no unknown gate (my reading).
- One experimental pulse, judged by trace error because the true pulse is unknown; the simulated pulse family bounds what the network can retrieve (my reading).

**Use for this project.**
- Adopt the two-stage pattern for ENZ-gated double-blind FROG: a fast first estimate, then gradient refinement through a differentiable trace model. In WP5 that model would be $|\mathcal{F}\{E_{probe}(t)\,G(t-\tau)\}|^2$ with the complex gate $G$ unknown; it fits the differentiable solver of WP3–WP4.
- Do not trust a purely supervised retrieval trained on simulated gates: the transfer gap appears here even for the simplest FROG geometry (my reading).
- Remove trivial ambiguities from training targets and from any loss that compares fields. The set quoted here holds for SHG-FROG. The set for a double-blind trace with a non-instantaneous gate differs, is not treated in this paper, and must be worked out separately.
- The paper says nothing on time resolution or ENZ materials. The plan groups it under "ENZ pulse metrology"; it is general FROG retrieval. The plan's author list matches arXiv v1, not Optica.

---
### [59] Zhang, Xu et al., "Hybrid deep reconstruction for vignetting-free upconversion imaging through scattering in ENZ materials" — Light Sci. Appl. 15, 327 (2026); arXiv:2508.13096

**Source read:** arXiv v1 full text (no TeX source) + Fig. 1 inspected in the PDF. The published version and the six supplementary videos were not read. H. Zhang and Y. Xu are equal-contribution first authors; 13 authors (UCLA/SLAC, Rochester, Ottawa, AFRL, Clemson, CREOL).
**Role in this project:** WP5 laboratory capability: shows which ENZ sample, beam geometry and nonlinear process the experimental side already operates. It offers no model or dataset for WP1–WP4.

**What it does.** An image-carrying infrared beam passes a scattering medium and overlaps with a femtosecond gate pulse on an ITO film at its ENZ wavelength. Four-wave mixing in the film upconverts only light arriving with the gate, which rejects multiply scattered photons. A two-stage network ("DeepTimeGate": supervised U-Net, then a small self-supervised deep-image-prior refinement) reconstructs the unscattered image and removes the vignetting of the gate. Reported gain over raw gated images: PSNR +124 %, SSIM +231 %, intersection-over-union about 10×.

**Model / method.**
- Nonlinear process (Sec. 2.2): third-order sum-frequency four-wave mixing with two pump (gate) photons and one signal photon,

$$\omega_{out} = 2\omega_{pump} + \omega_{signal}$$

  labelled "3ω" and "Sum FWM" in Fig. 1a. The visible output is recorded with a CMOS camera; an infrared camera records the unscattered mask as the training label.
- Geometry (Fig. 1a): gate and signal beams are non-collinear on the film; the upconverted beam leaves in a third direction.
- Reconstruction: a U-Net, then a three-layer convolutional network that learns a residual correction. No physical model of the ENZ response enters.

**Numbers.**
- Sample: ITO film of "subwavelength" thickness (value not stated); real permittivity crosses zero near 1510 nm (Sec. 2.2). Read off the Fig. 1a inset: Re ε ≈ 3.2 at 600 nm and ≈ −1.6 at 1800 nm; Im ε ≈ 0.35 at the crossing and ≈ 0.6 at 1800 nm. Substrate, deposition, carrier density, mobility: not stated. "M.Z.A. supported with ITO fabrication."
- Pulses: a "femtosecond" gate; "FWM at the ENZ wavelength". Wavelengths, durations, repetition rate, energies, intensities, angles, polarisations, conversion efficiency, gate window: not stated. The gate experiment itself is reported in the paper's Ref. 11 (Xu et al., arXiv:2503.21113; not in this library).
- Scatterers and data (Sec. 6.2, 6.4): polystyrene-bead suspensions (0.175–0.525 g/cm³) with USAF targets; 600-grit, 1500-grit and stacked diffusers with vortex-phase targets; 2,068 image pairs (Fig. 1 gives 1653/207/207 = 2,067).
- Averages (Table 1 and text): DeepTimeGate PSNR 31.8 dB, SSIM 0.83, IoU 0.25; U-Net alone 31.6 dB, 0.82, 0.24.

**Key results.**
- The ENZ four-wave-mixing gate works as an imaging time gate through bead suspensions and ground-glass diffusers (measured).
- The U-Net gives nearly all of the numerical gain; the refinement mainly removes speckle and broken edges, and the authors call its metric gains "incremental" (measured).

**Assumptions and limitations.**
- Demonstrated "in the weak scattering regime" (authors, Sec. 1).
- An imaging paper: no pump–probe transient, no time constant, no intensity dependence of the ITO response (my reading).

**Use for this project.**
- Capability indicated for WP5: an ITO film with ENZ wavelength near 1510 nm; femtosecond pulses at or near that wavelength (inferred from "3ω" and "FWM at the ENZ wavelength"; no numbers given); a non-collinear two-beam geometry on the film; detection of the upconverted signal on a silicon camera. Per the contribution statement, Y.X. and seven co-authors carried out the time-gated measurements. This answers part of the plan's open question on available ENZ wavelengths. Thickness, pulse duration and usable intensity must be asked for directly or taken from arXiv:2503.21113.
- The process shown is a parametric upconverting gate, $E_{out}\propto E_{pump}^2E_{signal}$, not the transient transmission or reflection gate that the plan calls "ENZ-gated double-blind FROG". The same film and geometry would support either; a spectrally resolved version of this signal would be a cross-correlation FROG (my reading; not discussed in the paper).
- Whether this gate is instantaneous on the pulse time scale is not addressed. WP1–WP3 could predict it.
- Data: "The code and datasets associated with this manuscript are publicly available at: https://doi.org/10.5281/zenodo.16551120." These are images, not pump–probe traces.
- The plan cites it as "related ENZ pulse metrology"; it is ENZ time-gated imaging.

---
### [60] Ball et al., "A space-time knife-edge in epsilon-near-zero films for ultrafast pulse characterization" — Laser Photonics Rev. 19, 2401462 (2025); arXiv:2408.00577

**Source read:** arXiv v1 full text including the supplement (no TeX source); equations and figures checked in the PDF. The published version was not compared. v1 has no "Methods" section although the text refers to one; experimental details are in "Supplemental 2".
**Role in this project:** WP5 and H3: the one paper of [58–60] that uses the transient index change of an ENZ film to characterise a pulse, and a quantitative warning about geometric time smearing in non-collinear pump–probe measurements.

**What it does.** A pump pulse strikes a thin Al-doped ZnO (AZO) film at an angle, so it reaches different positions on the film at different times and writes a refractive-index change that sweeps across the surface: the "space-time knife edge". A smaller, weaker probe is deflected by the local index gradient, and a quadrant detector records the deflection against pump–probe delay. Fitting two such traces with an analytic rate-equation model returns the pulse duration, the film relaxation time, the beam width and the angle of arrival. The method is an intensity cross-correlation with spatial sensitivity. It is not a FROG: nothing is spectrally resolved and there is no phase-retrieval algorithm.

**Model / method.**
Beam irradiances in the sample frame (Eqs. 1a, 1b); $n_g$ group index, $\tau_{pu},\tau_{pr}$ temporal widths, $w_x,w_y$ spatial widths, $(x_o,y_o)$ pump offset, and $a,b,c$ the coefficients of a Gaussian rotated by the angle of incidence $\theta$ (Eq. S1a):

$$I_{pump}(x,y,t) = n_g I_0\, e^{-\left(a(x-x_o)^2+2b(x-x_o)(y-y_o)+c(y-y_o)^2\right)}\, e^{\frac{-(t-t_0)^2}{2\tau_{pu}^2}},\qquad I_{probe}(x,y,t) = I_0\, e^{\left(-\frac{x^2}{2w_x^2}-\frac{y^2}{2w_y^2}\right)} e^{\frac{-t^2}{2\tau_{pr}^2}}$$

$$\begin{bmatrix} a & b\\ b & c\end{bmatrix}=\begin{bmatrix}\frac{\cos^2\theta}{2w_x^2}+\frac{\sin^2\theta}{2w_y^2} & \frac{\sin 2\theta}{4w_x^2}+\frac{\sin 2\theta}{4w_y^2}\\[4pt] \frac{\sin 2\theta}{4w_x^2}+\frac{\sin 2\theta}{4w_y^2} & \frac{\sin^2\theta}{2w_x^2}+\frac{\cos^2\theta}{2w_y^2}\end{bmatrix}$$

Material response (Eqs. S2a, S2b): the absorbed energy density $U$ follows the pump irradiance with absorption $\alpha$ and decays with one relaxation time $\tau_R$; the index change is taken proportional to $U$:

$$\frac{dU(x,y,t)}{dt}=\alpha I(x,y,t)-\frac{U(x,y,t)}{\tau_R},\qquad U(x,y,t)=e^{-t/\tau_R}\int_{-\infty}^{t}\alpha I(x,y,t')\,e^{-t'/\tau_R}\,dt'$$

(as printed; looks inconsistent because the integrating factor inside the integral should be $e^{+t'/\tau_R}$; also, the off-diagonal terms of the matrix above are printed as a sum, whereas a rotated Gaussian normally has their difference). Normal-incidence solution (Eq. S3), with $A(x,y)$ the spatial profile:

$$U(x,y,t)=\left(A(x,y)\frac{\sqrt{\pi}}{\sqrt{2}}\,\tau_p\, e^{\frac{\tau_p^2}{2\tau_R^2}}\left[\mathrm{erf}\!\left(\frac{\tau_R t+\tau_p^2}{\sqrt{2}\,\tau_p\tau_R}\right)+1\right]\right)e^{-\frac{t}{\tau_R}}$$

(as printed; with the corrected sign the erf argument would contain $\tau_R t-\tau_p^2$; check before implementing). Oblique-incidence solution used for the fits (Eq. 2), with $v$ the wavefront velocity outside the film:

$$U(x,y,t)=e^{\frac{-t_{diff}}{\tau_R}}\left[\frac{A\cdot\left(\mathrm{erf}\!\left(\frac{2v\tau_R C-1}{B}\right)-1\right)}{B}\right],\quad A=-\sqrt{\pi}\,e^{-\frac{4v\tau_R x\left(acv\tau_R x-b^2v\tau_R x+b+a\right)-1}{B^2\tau_R^2}},\quad B=2\sqrt{c+2b+a}\,v,\quad C=bx+ax-ctv-2btv-atv$$

(as printed; the supplement's Eq. S4 has $B\tau_R$ in the erf denominator and the prefactor $e^{-t(x)/\tau_R}$; S4 is the dimensionally consistent form). The pump arrival time varies across the spot by

$$t_{diff}=\frac{w_x}{v}\sin(\theta).$$

Measured signal (Eqs. 3a–3c): probe displacement at a detector a distance $L$ away, and normalised deflection on the quadrant cell ($n_{2,eff}$ effective nonlinear index, $t$ film thickness, $w_{pr}$ probe width):

$$\delta x=\frac{L\,n_{2,eff}\,I_{pu}\,t}{w_{pr}}\propto\frac{\partial U(x,y,t)}{\partial x},\qquad \Delta E_{n,x}(t)=\frac{E_{left}-E_{right}}{E_{total}}$$

The total detector signal also gives $\Delta T/T=\Delta E_{total}/E_{total}$.

Procedure: (i) probe offset perpendicular to the plane of incidence (position P4), deflection along $y$: the spatial gradient is fixed, so the trace is the pump–probe cross-correlation convolved with the relaxation, giving $\tau_p$, $\tau_R$ and $\Delta n$, assuming pump and probe are temporally identical. (ii) Probe at the pump centre (P1), deflection along $x$: a bipolar trace caused by the sweeping edge, giving $\theta$ or phase-front tilt and $w_x$. The fitting algorithm is not stated.

**Numbers.**
| Quantity | Value | Where |
| --- | --- | --- |
| Film | 240 nm Al:ZnO; substrate and deposition not stated | Suppl. 2 |
| ENZ crossover | ~1560 nm, ε″ = 0.46 (cursor in plot: ε′ = 0.003, ε″ = 0.457 at 1564.9 nm) | Suppl. 2 |
| Laser | 800 nm Ti:sapphire, 1 kHz, 100 fs nominal; 97 fs measured by FROG; two OPAs tunable 1200–1600 nm | Suppl. 2 |
| Pump | 1400 nm, peak irradiance 56 GW/cm² (incident or in-film not stated), ~600 µm (1/e² diameter), angles 10°, 20°, 30° | main text, Suppl. 2 |
| Probe | 1300 nm, ~170 µm (1/e² diameter), cross-polarised to pump, drawn at normal incidence in Fig. 1a | main text, Suppl. 2 |
| Bandwidth | 35 nm FWHM, both beams | main text |
| Fit results | $\tau_p$ ≈ 95 ± 5 fs; $\tau_R$ ≈ 400 ± 20 fs; $w_x$ ≈ 600 ± 30 µm; Δn = 0.15 | text, Fig. 3 |
| Sensitivity | $\Delta E_n$ ~ 10⁻³ (≈ 1 mrad) for Δn ~ 10⁻⁵ at ~500 MW/cm², with lock-in detection (asserted, citing Ferdinandus et al.; not shown here) | Discussion |
| Spatial sensitivity | 30 µm offset changes $\Delta E_{n,x}$ by 3.5 %; 0.1 % detectable, hence ~1 µm | text, Fig. 4 |
| Stated limits | electron–electron collision time assumed ~10 fs, thermalisation ~30 fs; nonlinear enhancement bandwidth ~300 nm, hence 8–10 fs at 1500 nm | Discussion |
| Cited, not measured here | damage threshold > 1 TW/cm² in the infrared; group index 4–10 | Introduction |

**Key results.**
- One parameter set fits all traces at 10°, 20° and 30°, and the fitted pulse duration, 95 ± 5 fs, agrees with the 97 fs FROG value (measured plus fit).
- The model has no intrinsic rise time: the index rises with the running integral of the pump. The data therefore do not require a finite material rise time at this pulse duration; the authors quote no bound (my reading).
- Recovery is a single exponential, 400 ± 20 fs, in this AZO film at 56 GW/cm² (fit).
- The $x$-deflection is bipolar and asymmetric in time for AZO, and symmetric for a 500 µm fused-silica control at 175 GW/cm², whose Kerr response is instantaneous (measured). The asymmetry is the signature of the finite relaxation.
- The angle-induced signal grows with angle (measured, 10–30°); the authors expect it to vanish for collinear beams.

**Assumptions and limitations.**
- Gaussian beams and pulses; pump and probe temporally identical; probe 3–5× smaller than pump (authors).
- "Because the technique employs a non-instantaneous nonlinearity, phase information (with the exception of the phase front tilt) is implicitly lost" (authors).
- For pulses shorter than the electron–electron collision time the distribution is athermal and "additional temporal signatures may arise"; they call the athermal response of ENZ films "still an open question" (authors).
- Index change proportional to absorbed energy with one relaxation time, with no intensity or depth dependence (my reading).
- $\tau_p$ is called a "1/e² temporal width" but enters as $e^{-t^2/2\tau_p^2}$; the conversion to FWHM is not given (my reading).
- AZO, not ITO; one film, one intensity, one wavelength pair (my reading).

**Use for this project.**
- Geometric smearing, relevant to H3 and to the WP5 design. From $t_{diff}$ with $v=c$, treating 600 µm as the full pump spot width: about 350 fs at 10°, 680 fs at 20°, 1 ps at 30°; across the 170 µm probe the pump arrival time still varies by about 100, 190 and 280 fs (my arithmetic). The angle that matters is the one between pump and probe. A spatially integrated transmission or reflection transient in such a geometry is broadened by this amount, independent of kinetics. The authors make the same point: in an angled Z-scan or reflection/transmission measurement "the space-time characteristics of the beam are not captured", and the resulting asymmetric focusing "can dominate the response at large angles (e.g. >30º)". The 1D FDTD of WP3 cannot represent this; WP5 needs a probe much smaller than the pump, a small pump–probe angle, and this term in the forward model (my reading).
- Gate shape for ENZ-gated FROG: a rise limited by the pump and a slow decay ($\tau_R$ ≈ 400 fs here), so the time resolution of a retrieved probe field would come from the rising edge only (my reading). Eq. S3 is a ready parametrisation of such a gate (a Gaussian convolved with an exponential) for a baseline or an initial guess.
- Feasibility: a 240 nm ENZ film gives deflection signals of 5–14 % at 56 GW/cm² without phase matching (Fig. 3). Deflection gives the index (phase) change and $\Delta T/T$ the amplitude change, so both parts of the gate are accessible even without spectral resolution (my reading).
- Uniqueness: no result on blind or double-blind retrieval. The fits are unique only under the Gaussian and identical-pulse assumptions. The statement that a non-instantaneous gate loses phase information applies to this intensity-type detection; [22] retrieves a complex, non-instantaneous gate from spectrally resolved traces and cites this paper as an intensity-profile method.
- Conflict: the authors take ~30 fs thermalisation as the speed limit of the framework, while [19], as summarised in the project plan, infers a 1–10 fs rise. This paper does not resolve the rise.
- Their Ref. 50 (Jaffray et al., Nat. Commun. 13, 3536 (2022)) is an ENZ-film FROG based on harmonic generation; it is not in this library. No data-availability statement in v1.
- The plan cites the arXiv preprint; the published reference is Laser Photonics Rev. 19, 2401462 (volume year 2025; [22] cites it with year 2024).

---

## Inference
### [61] Cranmer, Brehmer and Louppe, "The frontier of simulation-based inference" — Proc. Natl. Acad. Sci. 117, 30055 (2020); arXiv:1911.01429

**Source read:** arXiv v3 full text + TeX source (the whole review; nothing skipped) + Fig. 3 inspected. Journal pagination is from the reference manifest, not from the text read.
**Role in this project:** WP4 method reference; the plan cites it for posteriors "by Hamiltonian Monte Carlo or neural posterior estimation".

**What it does.** A short review of inference for models defined by a simulator whose likelihood $p(x|\theta)$ cannot be evaluated. It contrasts Approximate Bayesian Computation (ABC) and classical density estimation on summary statistics with methods built on neural density estimators, active learning, and access to simulator internals through automatic differentiation and probabilistic programming. It gives a taxonomy and qualitative recommendations, with no benchmarks and no simulation budgets.

**Model / method.**
Posterior (Eq. 1) and the source of intractability (Eq. 2); $\theta$ parameters, $x$ data, $z$ the simulator's latent variables:

$$p(\theta|x) = \frac{p(x|\theta) \, p(\theta)}{\int \mathrm{d} \theta' \, p(x|\theta') \, p(\theta')},\qquad p(x | \theta) = \int \mathrm{d} z\; p(x, z | \theta)$$

Quantities a simulator may expose beyond samples (Sec. 2C): I. $p(x|z,\theta)$; II. joint score $t(x,z|\theta) \equiv \nabla_\theta \log p(x, z | \theta)$; III. $\nabla_z \log p(x, z | \theta)$; IV. joint ratio $r(x,z| \theta, \theta') \equiv p(x, z | \theta) / p(x, z | \theta')$; V. $\nabla_\theta (x,z)$; VI. $\nabla_z x$.

Taxonomy (Sec. 3, Fig. 3 panels a–h):

| Method | Needs from the simulator | Amortised | Gives; remarks |
| --- | --- | --- | --- |
| ABC: rejection, MCMC, sequential Monte Carlo (a, b) | samples only | no | posterior samples. Needs summary statistics (prescribed or learned), a distance and a tolerance $\epsilon$; exact only as $\epsilon\to0$, where acceptance vanishes; rerun for each observation |
| Probabilistic programming (c, d) | control of every random draw; tractable $p(x\vert z,\theta)$ (I) or an ABC-like comparison | no | joint posterior of $\theta$ and $z$ |
| Differentiable simulator or tractable final stage (VI or I) | gradients or final-stage density | no | "asymptotically exact Bayesian inference" without tolerance or summaries (Sec. 3A) |
| Neural likelihood estimation, $\hat p(x\vert\theta)$ (e) | samples | yes | prior-independent; posterior needs a further MCMC or variational step |
| Neural posterior estimation, $\hat p(\theta\vert x)$ (f) | samples drawn using the prior | yes | posterior directly; "induces a prior dependence at every stage" |
| Neural likelihood-ratio estimation (g) | samples at two parameter points, or against the marginal | yes | a classifier, hence supervised and "often a simpler task"; like (e), "enables frequentist inference or model comparisons" |
| Surrogates trained with augmented data (h) | II and IV | yes | better sample efficiency |

Sequential (active-learning) variants update the proposal for $\theta$ toward the posterior of the observed data. This trades against amortisation, which needs surrogates "agnostic about the observed data" (Sec. 3D).

**Numbers.** None. No simulation counts, costs or accuracies are given.

**Key results.** (all asserted; this is a review)
- Recommendations (Sec. 3D): if any of I–VI is available, use it; "powerful algorithms are available for the case with differentiable simulators". Otherwise try a neural surrogate of the likelihood or the likelihood ratio "in most cases", the ratio when sampling from the surrogate is not needed.
- Diagnostics (Sec. 3C): calibration by parametric bootstrap, that is, simulate at known parameters and check coverage or credibility ("may require a large number of simulations"); a classifier trained to tell surrogate samples from simulator samples; checks of known expectation values of the estimators; varying reference parameters that should not change the result; ensembles; comparison with asymptotic distributions. "Passing these sanity checks does not guarantee that an estimator is correct, but failing them is an indication of problems."
- No diagnostic addresses misspecification; the usual remedy is more model flexibility and nuisance parameters.
- A time series is non-i.i.d. and must be treated as one high-dimensional observation (Sec. 1B); this applies to each pump–probe transient (my reading).

**Assumptions and limitations.**
- Written for stochastic simulators with intractable likelihood; qualitative throughout (authors).
- Bayes factors, evidence estimation, choice among discrete models and identifiability are not treated; "model comparisons" appears once. Hamiltonian Monte Carlo is not named, only MCMC and variational inference (my reading).
- The score is printed as $t(x|\theta) \equiv \nabla_\theta p(x|\theta)$ (as printed; looks inconsistent because the text defines it as the gradient of the log-likelihood).

**Use for this project.**
- WP4's forward model is deterministic and differentiable and its noise model is chosen by the analyst, so its likelihood is tractable: a "prescribed" model in the review's terms, not an implicit one. Hamiltonian Monte Carlo through the solver is then ordinary likelihood-based inference, and the review's first recommendation points the same way (my reading).
- Neural posterior estimation earns its place through amortisation: the synthetic identifiability study (gate G4) and fits to several datasets can reuse one trained network. It also marginalises nuisance inputs such as an uncertain pulse shape, since the network "implicitly marginalizes over all other (latent) variables". The costs are prior dependence and the calibration checks above (my reading).
- Bayes factors: the review gives no recipe. It notes that the evidence is an integral over $\theta$; a posterior estimator alone does not return it (my reading). Plan a separate evidence estimate and keep leave-one-condition-out prediction as the headline.
- Use the bootstrap coverage test on synthetic data from each of M0–M4 in the identifiability study. The warning that no diagnostic detects a wrong model family supports the plan's posterior-predictive caveat.
- Discrepancy: the review does not discuss Hamiltonian Monte Carlo, and among surrogates it recommends likelihood or ratio estimation, not posterior estimation. The plan marks this entry as cited from memory; journal, volume and first page agree with the manifest.

---

# Appendix. What was read for each paper

| Ref. | Key | Source read |
| --- | --- | --- |
| [1] | `alam2016` | abstract only (verbatim OpenAlex copy in `references/text/01_alam2016.abstract.txt`). Main text and Supplementary Materials NOT READ (paywalled, no PDF). Everything beyond the abstract is reconstructed from how papers in this l... |
| [2] | `caspani2016` | arXiv:1603.03581v2 full text + TeX source (accepted manuscript); figures checked in the PDF. "Read from figure" marks approximate values. |
| [3] | `reshef2019` | NOT READ — paywalled, no PDF and no abstract could be retrieved. Reconstructed from citing papers in this library and background knowledge; verify before use. |
| [4] | `kinsey2019` | NOT READ — paywalled, no PDF and no abstract could be retrieved. Reconstructed from citing papers in this library and background knowledge; verify before use. |
| [5] | `khurgin2021` | arXiv:2007.05613v1, full text (33-page Word-generated PDF; no TeX source, so equations were transcribed from the PDF pages and checked numerically where noted). The arXiv title is "Fast and slow nonlinearities in ENZ materials"... |
| [6] | `secondo2020` | arXiv:2004.00134v1, full text with Appendices A–D (Word-generated PDF, no TeX source; equations transcribed from the PDF pages). |
| [7] | `wang2019` | abstract only (verbatim OpenAlex copy in `references/text/07_wang2019.abstract.txt`). Main text NOT READ (paywalled, no PDF). Statements tagged **[lib: X]** are how library paper X describes this paper (Baxter = [9], Sarkar = [... |
| [8] | `secondo2022` | published PDF, full (cover page + 8 article pages; the file contains no supplementary material). |
| [9] | `baxter2022` | arXiv v1 full text with its Supplemental Information (SI) + TeX source; figures checked in the PDF. The published version was not available and may differ. Equation and figure numbers follow the arXiv v1 PDF (SI equations are 1... |
| [10] | `li2022` | arXiv v1 full text + TeX source; figures checked in the PDF. The Supplemental Material (set-up, Figs. S1–S4, Table S2 with fitted Drude and band parameters) is not in the library; the published version was not checked. Equation... |
| [11] | `scalora2020` | arXiv PDF, full (no TeX source; equations read from the page images). Published version not checked. Gaussian units. |
| [12] | `rodriguezsune2020` | arXiv PDF, full (no TeX source; equations read from the page images). Published version not checked. Gaussian units. |
| [13] | `tamashevich2023` | arXiv v2 full text + TeX source (appendices A–D, the path-integral derivations, skimmed). Published version not checked. Equation numbers follow the arXiv v2 PDF. |
| [14] | `sarkar2023` | arXiv:2207.14528v1 (29 Jul 2022), full text incl. appendices + TeX source; figure values read from the vector paths of the PDF. Published version not seen. arXiv title: "The electronic and thermal response of low electron densi... |
| [15] | `un2023` | arXiv:2210.08504v3 (30 Mar 2023), full text incl. appendices + TeX source; published version not seen. Figure values come from the vector paths of the authors' figure files; reading uncertainty is about a line width (±10 K on T... |
| [16] | `un2026` | arXiv:2601.10985v2 (26 Feb 2026), full text + TeX source; figure values from the vector paths of the PDF (±1 THz, ±0.01). The arXiv version lists three authors, Ieng-Wai Un, Subhajit Sarkar and Yonatan Sivan (the project plan c... |
| [17] | `zhou2020` | published PDF, full (main text, Methods, all four figures), cross-checked against the PMC full-text XML. Supplementary Notes 1-6 (setup, simulation model, higher-order effects, sample holder) were NOT available; anything that l... |
| [18] | `bohn2021` | arXiv:2107.05970v1, full text, Supplement 1, TeX source and figures. The preprint has no subtitle and three authors (Bohn, Luk, Hendry); the published version adds S. Horsley. Published version not compared. |
| [19] | `tirole2023` | arXiv:2206.04362v2 (PDF dated 16 June 2022), main text and Supplementary Materials in full, all figures. No TeX source. The published Nature Physics version was not compared; numbers are from the arXiv version. |
| [20] | `lustig2023` | PMC full-text XML (PMC11501756), main text and figure captions in full, plus the abstract file. No PDF, so figure panels were not seen; Supplementary Information not available. |
| [21] | `gurung2025` | published PDF, full (main text, Materials and Methods, Figs. 1-4), cross-checked against the PMC XML. The Supplementary Materials were NOT available; they hold the hydrodynamic equations (section S4), the fit function (eq. S23)... |
| [22] | `harwood2026` | arXiv v2 full text + TeX source, main text only. The Supplementary Information (SI) is not in the library; it holds the retrieval algorithm, the formula for $n_0(T)$, the Kane parameters, the depth-resolved transfer-matrix chec... |
| [23] | `ozlu2026` | arXiv v1 full text + TeX source, including the Supplementary Information (pump FROG, transmitted probe spectra, 100 fs simulation). |
| [24] | `choi2026` | arXiv v2 full text + TeX source, main text and Methods only. Supplementary Sections 1–13 are not in the library; they hold the two-temperature equations, the temperature-dependent heat capacity and electron–phonon coupling, the... |
| [25] | `segal2026` | arXiv v1 PDF, full (it has no methods section or supplementary material). |
| [26] | `galiffi2022` | arXiv:2111.08640v2 (18 Nov 2021) + TeX source; equation and figure numbers are those of the arXiv version. Read: Sec. 1 (introduction), Sec. 2 (temporal interfaces: boundary conditions, temporal slabs, anisotropic switching, di... |
| [27] | `hayran2022` | arXiv:2208.00459v2 (7 Sep 2022) full text + TeX source; figures checked in the PDF. The published OME version was not compared. |
| [28] | `saha2023` | published PDF, full (7 pages, 8267–8273), figures checked. |
| [29] | `lobet2023` | arXiv:2306.01314v2 PDF (42 pages), not the ACS version. Read: abstract, introduction, "Nonlinear properties of NZI media and their application to all-optical switching", "ENZ media for time-varying photonics", conclusions. Skip... |
| [30] | `solis2021` | arXiv:2103.06142v1 (10 Mar 2021) full text + TeX source. The published Photon. Res. version was not available; equation and figure numbers below are those of arXiv v1 and may differ in print. The arXiv title reads "in Presence ... |
| [31] | `horsley2023` | arXiv:2208.11778v1 (24 Aug 2022) full text + TeX source. The Supplementary Information is not in the arXiv files (only an uncaptioned SI figure), and the published PRL version was not available. The Drude-model derivation, all ... |
| [32] | `rethfeld2002` | abstract only (verbatim OpenAlex abstract in `references/text/32_rethfeld2002.abstract.txt`). Body NOT READ — paywalled, no PDF. Everything beyond the abstract is reconstructed from how [14], [36], [38], [39], [40] describe the... |
| [33] | `mueller2013` | abstract only (verbatim OpenAlex abstract in `references/text/33_mueller2013.abstract.txt`). Body NOT READ — paywalled, no PDF. Details beyond the abstract are reconstructed from [37], [38], [39], [40] and background knowledge ... |
| [34] | `carpene2006` | abstract only (verbatim OpenAlex abstract in `references/text/34_carpene2006.abstract.txt`). Body NOT READ — paywalled, no PDF. Details beyond the abstract are reconstructed from [9], [14], [15], [38], [39], [40] and background... |
| [35] | `dellavalle2012` | abstract only (verbatim OpenAlex abstract in `references/text/35_dellavalle2012.abstract.txt`). Body NOT READ — no PDF arrived in this session (the manifest lists an open copy at hal.science/hal-02394961). Details beyond the ab... |
| [36] | `dubi2019` | arXiv:1810.00565v4 (21 Jul 2019), full text with its Supplementary Information appendix, + TeX source. Authors: Y. Dubi, Y. Sivan. Published version not compared; equation numbers are those of the arXiv PDF. |
| [37] | `brown2016` | arXiv:1602.00625v2 (28 Jul 2016) full text + TeX source, including the supplementary permittivity tables shipped with the source. Published version not compared; equation numbers are those of the arXiv PDF. |
| [38] | `uehlein2022` | arXiv:2112.01851v1 (3 Dec 2021) full text + TeX source. The published Nanomaterials version was not available and was not compared; equation numbers are those of the arXiv PDF. |
| [39] | `uehlein2025` | arXiv:2503.09479v1 (12 Mar 2025) main text + Supplementary Information, PDF and TeX source. The published J. Appl. Phys. version was not available and was not compared. The SI carries a different title ("Non-equilibrium electro... |
| [40] | `roden2026` | arXiv:2601.11371v1 (16 Jan 2026) full text + TeX source. The published Phys. Plasmas version was not available and was not compared; equation numbers below are those of the arXiv PDF. |
| [41] | `maliyov2024` | arXiv:2311.07520v1 (13 Nov 2023) full text + TeX source. Published npj version not compared; equation numbers are those of the arXiv PDF. |
| [42] | `li2021fno` | arXiv:2010.08895v3 (ICLR version), full text + TeX source. |
| [43] | `kovachki2023` | arXiv:2108.08481v6, text + TeX. Read: abstract, Secs. 1, 2.3, 3, 4.4–4.5, 5.1, 7.2–7.4, 8.2–8.3 (statements), 10. Skipped: Secs. 4.1–4.3, 5.2, 6, 7.1, 7.3.4, 9, proofs, appendices. |
| [44] | `lu2021deeponet` | arXiv:1910.03193v3 (15 Apr 2020), full text + TeX source. The preprint has three authors (Lu, Jin, Karniadakis) and the title "DeepONet: Learning nonlinear operators for identifying differential equations based on the universal... |
| [45] | `li2021pino` | arXiv:2111.03794v4 (29 Jul 2023), full text with appendices. The journal version was not compared. |
| [46] | `xiao2022relaxnet` | arXiv:2211.08149v1 (15 Nov 2022), full text + TeX source; Figs. 2, 7–15, 18, 19 inspected in the PDF. The published JCP version was not compared. |
| [47] | `lee2024spol` | arXiv v1 (26 Feb 2024), full text + TeX source; Figs. 1–3 inspected in the PDF. Authors: J. Y. Lee, S. Schotthöfer, T. Xiao (equal contribution), S. Krumscheid, M. Frank. No venue is stated in the paper (it is typeset in the IC... |
| [48] | `lee2025fourierspecnet` | arXiv v2 (6 Mar 2026) full text + TeX source, including Appendix B (ablations). Proofs of Prop. 3.1 and Appendix A (decay of kernel modes) skimmed only. The journal reference is not printed in the arXiv PDF (not verified). The ... |
| [49] | `hu2025` | arXiv v2 (12 Nov 2026) full text + TeX source, all sections. Figures inspected for the BGK error curves. No journal reference is given; the acknowledgments mention anonymous reviewers. |
| [50] | `zhao2025` | arXiv v2 (4 Apr 2025) full text + TeX source, with appendices. The PRL reference is not printed in the arXiv PDF (not verified). |
| [51] | `burles2026` | arXiv v1 (31 Jul 2026) full text + TeX source, all sections. The PDF is marked "Under consideration for publication in J. Plasma Phys."; the journal reference is not printed (not verified). |
| [52] | `luo2026` | arXiv v1 (3 Apr 2026) full text + TeX source (main text, 7 pages) and its figures. The Supplemental Material, which holds the FNO architecture, training strategy and the large-time-step tests, is not in the library and was not ... |
| [53] | `shah2024` | arXiv v2 (25 Jul 2024) full text + TeX source (workshop paper, 9 pages). |
| [54] | `shah2025` | arXiv v3 (1 Dec 2025) full text + TeX source, with Supplementary Information (SI). The journal reference is not printed in the arXiv PDF (not verified). |
| [55] | `gu2022neurolight` | arXiv v1 full text + TeX source (main text, appendices A1–A5, tables); figures checked in the PDF. The NeurIPS camera-ready was not compared. |
| [56] | `ma2024pic2osim` | arXiv v1 full text + TeX source (main text, appendix A.1–A.6, tables). The APL Photonics version was not available; all numbers are from v1. |
| [57] | `aps2026pifno` | abstracts only, and NOT from the APS page (https://meetings-archive.aps.org/smt/2026/mar-h00-01/159/), which could not be retrieved. Everything below comes from two Belmont University repository records by B. L. Blackwood and S... |
| [58] | `zahavy2018` | arXiv v1 full text + TeX source + figures. This is NOT the Optica article: v1 is the four-page NIPS 2017 workshop paper ("Deep Learning Reconstruction of Ultra-Short Pulses"; authors Zahavy, Dikopoltsev, Cohen, Mannor, Segev). ... |
| [59] | `zhang2025` | arXiv v1 full text (no TeX source) + Fig. 1 inspected in the PDF. The published version and the six supplementary videos were not read. H. Zhang and Y. Xu are equal-contribution first authors; 13 authors (UCLA/SLAC, Rochester, ... |
| [60] | `ball2024` | arXiv v1 full text including the supplement (no TeX source); equations and figures checked in the PDF. The published version was not compared. v1 has no "Methods" section although the text refers to one; experimental details ar... |
| [61] | `cranmer2020` | arXiv v3 full text + TeX source (the whole review; nothing skipped) + Fig. 3 inspected. Journal pagination is from the reference manifest, not from the text read. |
| [62] | `un2026nanophot` | Published version (Nanophotonics 15, e70272, CC BY) all 11 pages + Supporting Information (9 pages), supplied by the project lead; figures read by eye. arXiv:2608.11756v1 read the day before. Added 2026-11-18. |
