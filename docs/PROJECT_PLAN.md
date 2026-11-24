# Maxwell-Kinetic Model of ENZ Electron Dynamics: Project Plan

Revision 2.2 · Nov 24, 2026 · @Yang

Revision 1 (Nov 16, 2026, "Learned Kinetic Closure for ENZ Electron Dynamics") is kept unchanged at `docs/archive/PROJECT_PLAN_v1_2026-11-16.md`. Reference numbers are the same in both; revision 2.1 adds \[62\].

## What changed in revision 2.2

Prompted by decision test GT2 (`docs/gates/GT2_fewcycle_closure_report.md`): is a learned few-cycle constitutive law worth building? The project lead asked for the test on 2026-11-24 and set the generative pulse-retrieval idea aside.

| Change | Reason | Evidence |
| --- | --- | --- |
| The learned operator stays deferred for few-cycle pumps as well. The condition "the density-matrix solver is too expensive in a film" is removed from the deferred section. | The density matrix of \[16, 62\] costs 11 s per cell for a 60 fs window. Replacing its populations by one electron temperature changes the field-resolved current by 1.3 % (median) and 5.8 % (worst) of its peak, less than the model's own uncertainty of 5-13 %. A fitted 17-pole law reaches 4.4 % at 6 ms per cell. | GT2 report |
| Tier 2b (density matrix of a particle) is built and validated: `enzkin.densmat`. | Needed as the truth for the test; reproduces the published results of \[16, 62\] to a few percent with one fitted number. The dipole element printed in \[16\] Eq. (16) is wrong and was reconstructed. | `docs/theory/tier2_density_matrix.md`, `scripts/validate_densmat.py` |
| Tier 2a uses a multi-pole response whose widths follow the electron temperature wherever the linear response is not a Drude pole; its validation target is the density matrix itself. | A single Drude pole misses the particle's current by 8-21 %; poles with temperature-independent widths miss it by up to 13 % although they fit the permittivity at every frozen temperature. | GT2 report |
| A nanocrystal route is added as a candidate: the particle model applies directly to field-resolved measurements on ITO nanocrystals. | It needs no film solver and runs in seconds. | to be assessed after reading Herbst et al. 2026 |

## What changed in revision 2.1

Prompted by Un, Sarkar & Sivan, Nanophotonics 15, e70272 (2026), now reference \[62\]: a density-matrix calculation, with the full collision integrals, of a 4 nm ITO sphere under a 5 fs pump.

| Change | Reason | Evidence |
| --- | --- | --- |
| Tier 2 (few-cycle pumps) is split. Tier 2a drives every kinetic cell of the WP3 film solver by the field-resolved J.E instead of the cycle-averaged absorbed power, and is validated against \[62\] in the nanosphere limit. Tier 2b, the density matrix of \[16, 62\], is kept only for coherent signatures. Flagship claims stay at pumps of 20 fs and longer until Tier 2a passes. | In \[62\] a thermal model driven by the instantaneous net absorbed power reproduces the density-matrix permittivity for a 5 fs pump. The density matrix absorbs through the finite size of the particle, which does not carry over to a film. | `references/LITERATURE_SUMMARY.md`, \[62\]; `docs/theory/un2026_notes.md` |
| The electron-phonon coupling (deformation potential D) and the impurity damping model become free material parameters in WP4, with a prior covering at least D = 8.6-17.2 eV. | The same group's two parameter sets differ by four in the coupling (\[14\]: 17.2 eV, no impurity term; \[62\]: 8.6 eV plus Brooks-Herring impurities). The linear permittivity fixes only the total damping; with the two sets a 100 fs transient recovers in 0.5 and 1.9 ps, against 60-400 fs measured. | `scripts/check_un2026.py` |
| WP1 has an optional Brooks-Herring impurity damping and the parameter set of \[62\] (done). | Needed to run the authors' current parameter set. | `src/enzkin/collisions/eimp.py`, `ito_un2026()` |
| M0 may take its electron-electron damping from the generalised Fermi-liquid rate of \[62\] instead of a free power law. | The formula reproduces the exact integral within 3 % near the Fermi level (checked with the corrections to its printed form). | `docs/theory/un2026_notes.md` |
| The observation table and H3 note that \[62\]'s theory predicts a larger change and a slower recovery than \[20\] measured, and that its authors point to film effects. | Supports testing H3 on \[20\]. | \[62\] Sec. 5 |

## What changed in revision 2

| Change | Reason | Evidence |
| --- | --- | --- |
| The neural-operator surrogate (old WP2) is dropped for now and moved to a deferred section. | Gate 1: the exact collision integral costs 1.7 ms per evaluation at 128 energy points, so a film run takes about half an hour on four CPU cores instead of the 1-10 days the plan assumed. | `docs/gates/G1_report.md` |
| The film solver (WP3) uses the exact integral and starts now. | It no longer waits for a surrogate. | same |
| The FFT random-k integral (B0) is a comparison model only. | It misses the exact permittivity transient by 7-221 %. | same |
| Hypothesis H2 is stated for Re ε and the plasma frequency, not for the sign of the transmittance. | A fixed-density thermal model already reverses the transmittance \[24\]. | `references/LITERATURE_SUMMARY.md`, Part I §5 |
| The two-photon model (M3) is tested only on 800 nm data. | With the 2.67 eV gap of \[23\] two-photon absorption is closed for pumps beyond 929 nm. | same |
| M0 is a two-temperature model with a temperature-dependent damping law. | \[38\]'s "extended" model adds non-thermal electrons; M0 needs the damping law to produce any reversal \[24\]. | Part I §2, §5 |
| The form of the current equation is derived from the kinetic equation, as a WP3 task. | \[30, 31\] give two mathematical forms but do not tie them to a mass change or a density change. | Part I §6 |
| The intensity range is limited to what a single band can hold, and the band width is a model parameter. | A third of the old sampling range saturates a 4 eV band; \[15\] is reproduced only with a wider band. | `docs/WORKLOG.md` §1.5, §1.7 |
| The observation table and the dataset table carry qualifiers; 26 citations are corrected. | Checked against the papers, Crossref and arXiv. | Part I §7, §10; `references/README.md` |
| The methods paper becomes a solver-and-benchmark paper without machine learning. | No operator to report. | this document |

## Summary

The project answers one physics question: what sets the speed and sign of the optical response of epsilon-near-zero (ENZ) indium tin oxide at high intensity and short pulse duration? Several groups give different answers to the same anomalies, and no existing model couples distribution-resolved electron kinetics to Maxwell's equations in a film.

The method is a Boltzmann solver for the conduction-electron distribution, with the full Pauli-blocked electron-electron collision integral, embedded in a time-domain Maxwell solver as the constitutive law of every cell. The coupled model is then used as a forward model to rank the competing mechanisms against data.

**What is new**

- The first film-scale Maxwell-kinetic simulation of ENZ pump-probe experiments, with no adiabatic approximation for the probe.
- A current equation for a heated non-parabolic band derived from the kinetic equation, which decides between the two forms used in the literature.
- A quantitative, Bayesian comparison of the published mechanisms on shared data, plus one prediction tested with double-blind FROG.
- An open, exactly conserving kinetic solver for degenerate electrons in ITO, with documented corrections to the published formulas.

**Outputs**

| Output | Target venue | Month |
| --- | --- | --- |
| Flagship letter: origin of the ENZ response time and sign reversal | PRL or Optica | 12 |
| Methods paper: solver, corrections to the published model, benchmark dataset and code | Phys. Rev. Applied or Comput. Phys. Commun. | 9 |

Gate 1 is closed. Three gates remain, plus a short confirmation gate (G1b) that the exact solver stays affordable across the whole condition matrix and on a GPU.

## Background and gap

ITO near its ENZ wavelength gives order-unity index changes in under a picosecond \[1, 2\], but theory only explains the low-intensity, long-pulse regime. Bracketed numbers refer to the reference list at the end.

**Established picture.** Sub-bandgap light heats conduction electrons. In a non-parabolic band the hot population is heavier on average, so the plasma frequency drops and the ENZ point redshifts \[5-8\]. Four levels of theory exist:

- Energy-balance and two-temperature models \[6-9\].
- Hydrodynamic Maxwell models with nonlocal and hot-electron terms \[11, 12\].
- Boltzmann-type kinetics in energy space \[14, 15\]. At 220 fs the kinetic and thermal permittivities coincide; non-thermal electrons matter for the rise under pumps of 30 fs and shorter.
- An energy-space density-matrix model for few-femtosecond pulses \[16\], applied to a 5 fs pump on a 4 nm sphere in \[62\]: the change is of order one and cumulative, electron-electron thermalisation takes about 12 fs at TW/cm2, and a thermal model driven by the instantaneous absorbed power reproduces the permittivity.

**Observations that these do not jointly explain**

| Observation | Source | Mechanism proposed |
| --- | --- | --- |
| Reflectivity rise of 1-10 fs fitted under a pump of 225 fs field FWHM (159 fs in intensity); the fit is one parameter of a non-dispersive r(t) of a saturated ITO/Au mirror, with no error bar | Tirole et al. \[19\] | None at the time; non-thermal electrons proposed later \[15\] |
| Fast rise and recovery (10-30 fs in transmission, hundreds of fs in reflection) under a 6 fs, 800 nm pump | Lustig et al. \[20\] | Open; the authors point to film thickness. The intraband theory of \[62\] predicts a larger change and a recovery set by slow electron-phonon cooling, and names film gradients, an interband Kerr term or multiphoton absorption as possible causes |
| Sign reversal of ΔR/R and ΔT/T under a 6.4 fs, 800 nm pump; a fit crossover at 4.9 TW/cm2, full reversal near 8-10 TW/cm2 | Ozlu et al. \[23\] | Interband two-photon absorption once intraband excitation lifts Pauli blocking |
| Non-monotonic, oscillatory complex Fresnel coefficients above 55 mJ/cm2 (136 fs, 1300 nm) | Harwood et al. \[22\] | Auger-like interband relaxation added to an energy-balance model |
| Oscillatory, sign-reversing transmittance over a few optical cycles (simulation, not a measurement) | Choi et al. \[24\] | Extreme electron temperatures at fixed carrier density |
| Relaxation time rising 2.7-fold (ITO) and 2.9-fold (AZO) with intensity, and by 23 % (ITO) and 42 % (AZO) with ENZ-mode excitation | Gurung et al. \[21\] | Hydrodynamic model |
| Larger time-refraction shift as the index change becomes sub-cycle (CdO, 4 um) | Segal et al. \[25\] | Two-photon interband absorption, not modelled quantitatively |

**Gaps**

1. **No distribution-resolved kinetic model is coupled to Maxwell in a film.** Temperature-level and hydrodynamic models have been run in time-domain Maxwell solvers \[9, 11, 12\], but they do not resolve the electron distribution. The kinetic model uses a time-dependent permittivity in the adiabatic approximation \[15\]. The density-matrix model treats a 4 nm sphere in the quasi-static limit and names layers as future work \[16, 62\]. The experiments use 20-620 nm films, often at oblique incidence, where field enhancement, depth-dependent absorption and probe dispersion all shape the measured transient.
2. **The mechanisms have never been tested against each other.** Each group fits its own data with its own model. No study runs several mechanisms through one forward model on a shared dataset.
3. **The current equation of a heated ENZ film has not been derived.** Models insert a time-varying plasma frequency into a Drude current equation in one of two ways \[9, 11, 31\]. The two give the same new frequency after a change but different amplitudes, and nobody has shown which one the kinetics imply.

**Why it matters.** Optical time refraction, time diffraction and proposed photonic time crystals all rely on these oxides \[17-19, 26, 28\], and the energy cost of fast index modulation is severe \[27\]. The reviews themselves call the ultimate response time an open question \[26, 29\].

## Central question and hypotheses

**Question.** What sets the rise time, recovery and sign of the ENZ response in ITO as pump intensity rises from 1 GW/cm2 toward the limit of the single-band model and pulse duration falls from 200 fs toward a few optical cycles?

**Hypotheses**

|  | Hypothesis | Falsified if |
| --- | --- | --- |
| H1 | Below the reversal threshold, intraband non-thermal kinetics at fixed carrier density explain rise and recovery. | One parameter set cannot fit transients across pulse durations. |
| H2 | A sign reversal of the change in Re ε (equivalently, a plasma frequency rising above its equilibrium value) requires a channel that changes carrier density. Two-photon absorption \[23\] and Auger-like excitation \[22\] are distinguishable by how they scale. | A fixed-density model reproduces a reversal of Re ε, or the two channels fit all data equally well. |
| H3 | Part of the apparent speed is electromagnetic. The probe crosses a time-varying, dispersive, depth-inhomogeneous film, so the measured rise is not the kinetic rise. | Maxwell-kinetic and local adiabatic models agree within experimental error in all conditions. |

H2 is deliberately not stated for the transmittance: a fixed-density thermal model with rising damping already flips its sign \[24\].

Discriminators for H2. Two-photon absorption follows the instantaneous intensity squared, gated by the depletion of the band bottom, so at fixed fluence it strengthens as the pulse shortens; \[23\]'s own simulation shows no reversal at 100 fs for the same fluence. Auger-like excitation follows the hot-electron population, so it tracks absorbed energy. Lag is a weak discriminator, because both published reversals peak 120-150 fs after the pump. The energetic constraint is the strong one: two 1300 nm photons carry 1.91 eV, below the 2.67 eV gap of \[23\].

Evidence so far for H1 and H3, from WP1: a thermal model built from the same operators misses the kinetic permittivity transient by 11-17 % for 20 fs pumps and 3-5 % at 220 fs; and several experiments already report film effects on the measured dynamics \[17, 18, 20, 60\]. \[62\] adds that intraband kinetics in a uniform particle recovers on the electron-phonon time, far slower than the 10-40 fs transmission recovery of \[20\]; if a film calculation confirms this, the fast recovery is outside H1 and a candidate for H3. Recovery times also depend on the electron-phonon coupling, which is not settled (WP4).

**Target claims for the flagship paper**

1. One parameter set reproduces published transients for pump durations of 20-225 fs and intensities from 1 GW/cm2 up to the single-band limit, within the stated experimental uncertainty.
2. The mechanisms are ranked by held-out predictive score and Bayes factor. The Auger-versus-two-photon comparison is restricted to 800 nm data.
3. The measured rise time is decomposed into a kinetic part and an electromagnetic part.
4. One prediction chosen by optimal experimental design is confirmed or refuted with double-blind FROG.

**Scope limits.**

- The Boltzmann description is cycle-averaged. One optical cycle at 1.2 um is 4 fs, so it is reliable for pumps longer than about 20 fs. For the sub-8 fs experiments \[20, 23\] the first step is Tier 2a (field-resolved absorption in the film solver); the coherent density-matrix level \[16, 62\] (Tier 2b) is needed only for coherent features. The field also displaces the Fermi sphere by 0.5-1.3 k_F in those experiments (`scripts/field_regime.py`), an intra-cycle effect neither level contains.
- A single conduction band of finite width holds a finite energy, about 5 x 10^8 J/m3 above the cold state for a 4 eV band. Conditions that deposit more are outside the model. As a rule of thumb this is an in-film intensity near 10^11 W/cm2 for a 100 fs pulse.

## WP1: Ground-truth kinetic solver (months 0-3)

**Status: first pass complete** (`src/enzkin`, `docs/theory/wp1_kinetic_model.md`, working log session 1).

**State and equation.** Each spatial cell carries an isotropic conduction-band distribution f(E, t) in a Kane band and a lattice temperature.

```latex
\frac{\partial f}{\partial t} = \left.\frac{\partial f}{\partial t}\right|_{\mathrm{abs}} + Q_{ee}[f] + Q_{e\text{-}ph}[f, T_l] + S_{\mathrm{inter}}[f, I], \qquad \frac{\hbar^2 k^2}{2 m_0^*} = E\,(1 + C E)
```

```latex
\omega_p^2 = \frac{e^2}{\varepsilon_0 m_0^*} \int g(E)\, f(E)\, w(E)\, dE, \qquad w(E) = \frac{1 + \tfrac{8}{3} C E (1 + C E)}{(1 + 2 C E)^3}
```

| Term | Treatment as built | Basis |
| --- | --- | --- |
| Intraband absorption | One-photon steps with a constant matrix element, normalised to the absorbed power density | \[14, 15\] |
| Electron-electron | Full Pauli-blocked integral with a statically screened Coulomb kernel, screening recomputed from f; conserves number and energy to round-off | \[14\], corrected |
| Electron-phonon | Acoustic deformation potential with momentum conservation, as a pair master equation with a symmetric cut-off | \[14\], corrected |
| Interband channels | Auger-like relaxation \[22\] and two-photon absorption with state filling \[23\], lifted from the published lumped models | \[22, 23\] |

Comparison models built from the same pieces: the FFT random-k integral of \[40\], a conserving relaxation-time model, and a two-temperature model.

**Done**

- Conservation of number and energy to 10^-13; relaxation to Fermi-Dirac; two-temperature limit for long pulses; 60 tests.
- Validation against \[15\]: absorptivity within 1-6 %, peak electron temperature within 2-7 % up to 22 GW/cm2 and within 5 % at 50-75 GW/cm2 with a 6 eV band.
- Gate 1 measured (see Timeline and gates).

**Remaining**

- [x] Confirm grid convergence over the full condition matrix, including weak excitation and the 6 eV band (gate G1b). Done 2026-11-17: an energy step of 31 meV holds 1 % everywhere (`docs/gates/G1b_report.md`).
- [ ] Port the exact integral to JAX for GPU and for gradients; benchmark it on a GPU (gate G1b). Port done and tested; the GPU benchmark is outstanding.
- [ ] Fix the band model for high intensity: choose the band width, state the validity limit. Rule found at G1b: widen the band (4, 6, 8 eV at a fixed 31 meV step) when the occupation of its top cell passes 4 × 10^-5. Still to do: rerun the validation against \[15\] and the Gate 1 conditions above about 4000 K with that rule, and state the limit above which the single-band model itself is not trustworthy (results above roughly 7000 K change by 1-3 % between 6, 8 and 10 eV).
- [ ] Close the validation gaps against \[15\]: the change of Re ε is 20-28 % low, and Im ε depends on how the screening is treated. This needs the authors.
- [ ] Write to the authors of \[14-16, 62\] for the parameters the papers do not print, and ask which parameter set (D = 17.2 eV, or 8.6 eV with impurities) applies to which sample (list in `references/LITERATURE_SUMMARY.md`, Part I §11).
- [x] Optional Brooks-Herring impurity damping and the parameter set of \[62\] (`ito_un2026()`); the rates of \[62\] checked against enzkin (`scripts/check_un2026.py`, `docs/theory/un2026_notes.md`). Done 2026-11-18.
- [ ] Benchmark dataset: about 10^3 trajectories with the exact integral over 20-250 fs and 1.1-1.6 um, inside the single-band limit, with material parameters varied by 25 %. It serves regression tests, the methods paper, and any later surrogate.
- [ ] Compare with the reduced athermal model of \[39\] as a second baseline (optional).

**Tier 2: few-cycle pumps.** The sign-reversal data \[20, 23\] and the Auger-versus-two-photon test both sit in this regime.

- *Tier 2b, built (2026-11-24).* The energy-space density matrix of \[16, 62\] for a small particle in a uniform field: `enzkin.densmat`. With the dipole element reconstructed (the printed Eq. (16) of \[16\] fails the oscillator-strength sum rule) and ε∞ fitted, it reproduces the cold permittivity, the absorbed energy and power, the local-field spectral shift and the permittivity change published in \[16, 62\] to a few percent. Cost: 11 s per run. It is a particle model: it absorbs through the particle's finite size, so it must not be used for a film cell as it stands.
- *Tier 2a.* In the WP3 film solver, drive each cell's absorption by the instantaneous J.E from the field solver instead of the cycle-averaged power (kinetic cells, or thermal cells for M0). Test GT2 showed that an electron temperature heated by the instantaneous work, feeding a linear response whose damping follows that temperature, tracks the density matrix to 1-6 %. It also showed that the damping must be updated with the temperature in time, not through a frozen-temperature fit. For a bulk-like film the linear response is a Drude pole; for anything with a non-Drude dispersion use fitted poles (`ReducedResponse.fitted`).
- *Nanocrystal route (candidate).* Field-resolved measurements on ITO nanocrystals under two-cycle pulses (Herbst et al., Adv. Sci. 13, e16818 (2026), cited by \[62\]; not yet read) could be compared with `enzkin.densmat` directly, without a film solver.
- *Still missing for films.* A coherent treatment of collision-assisted absorption, and the field-driven coupling between coherences, which carries the intra-cycle band nonlinearity when the field displaces the Fermi sphere by 0.5-1.3 k_F.

## WP3: Maxwell-kinetic closure (months 1-7)

WP3 embeds the exact kinetic solver in a one-dimensional time-domain Maxwell solver. Pump and probe both propagate through a film in which every cell evolves kinetically, with no adiabatic approximation. It starts now. (The work-package numbers of revision 1 are kept; WP2 is in the deferred section.)

**Solver**

- 1D FDTD through the film thickness. Normal incidence first, then TM oblique incidence, which is needed for ENZ-mode excitation \[21\].
- Films of 20-310 nm on a substrate, with 2-5 nm cells: roughly 100 kinetic cells.
- Multirate stepping: fields at about 5 as, kinetics at 0.25-1 fs with a rate-limited step.
- The cycle-averaged absorbed power J.E in each cell drives that cell's absorption term; for few-cycle pumps the instantaneous J.E does (Tier 2a).
- Written in JAX, so the same code runs on CPU and GPU and gives gradients for WP4.

**Constitutive law: a derivation task.** The literature uses two current equations for a time-varying plasma frequency:

```latex
\text{(a)}\quad \frac{\partial J}{\partial t} + \gamma J = \varepsilon_0\, \omega_p^2(t)\, E
\qquad\qquad
\text{(b)}\quad \frac{\partial J}{\partial t} + \left(\gamma - \frac{d \ln \omega_p^2}{dt}\right) J = \varepsilon_0\, \omega_p^2(t)\, E
```

Form (a) is used in \[9, 11, 12\] and is the switched-density model of \[30\]; form (b) is the main form of \[31\] and, at fixed density, equals revision 1's equation with the effective-mass term. Reading (a) as "density change" and (b) as "mass change" is this project's hypothesis, not a result of \[30, 31\]. Which form heated ITO obeys depends on whether the heating collisions conserve momentum: electron-electron scattering does, assisted absorption does not. The task is to derive the equation from the l = 1 moment of the kinetic equation and, if needed, to carry the current-carrying part of the distribution in each cell. When benchmarking against \[9\] the extra term is off.

**Cost.** One run needs about 100 cells x 4000 steps x 4 stages = 1.6 x 10^6 collision evaluations: about 1.8 core-hours at 128 energy points, 16 at 256 (measured single-thread timings). No surrogate is needed.

**Validation ladder**

1. Linear limit matches transfer-matrix reflection and transmission, including ENZ-mode absorption.
2. Intensity-dependent reflection and transmission of Alam et al. \[1\].
3. Time-refraction frequency shifts of Zhou et al. \[17\] and spatiotemporal refraction of Bohn et al. \[18\].
4. Relaxation-time trends with intensity and angle of Gurung et al. \[21\].

Parameters may be tuned only at step 2. Steps 3 and 4 are predictions. The FDTD two-temperature model of \[9\] is the comparison baseline throughout; note that it reproduces \[1\] with two fitted parameters.

**Test of H3.** Run identical kinetics two ways: full Maxwell-kinetic, and local adiabatic permittivity with thin-film Fresnel formulas (the quasi-static model of \[18\] is exactly this baseline). Compare the probe transients and report the difference in apparent rise time against film thickness, angle and probe bandwidth.

**Deferred.** Lateral structure and non-collinear geometry: in \[60\] the pump front sweeps a 600 um spot in 0.35-1 ps, which a 1D solver cannot represent. Nonlocal hydrodynamic terms \[11, 12\]: the 35-50 nm blueshift in \[12\] was computed with a band-edge mass of 0.033 and is about 30 times smaller for realistic ITO parameters, so it is not expected to matter.

**Deliverable.** The first film-scale Maxwell-kinetic simulation of an ENZ pump-probe experiment.

## WP4: Mechanism discrimination (months 5-10)

WP4 runs the published mechanisms through one forward model and lets the data rank them.

**Model family**

| Model | Physics | Extra free parameters | Source |
| --- | --- | --- | --- |
| M0 | Thermal electrons only: two-temperature model with full Fermi integrals, temperature-dependent chemical potential, and a damping law that rises with electron temperature. | Electron-phonon coupling; damping exponent and reference temperature, or the generalised Fermi-liquid electron-electron rate of \[62\] with an Umklapp fraction | \[9, 24, 62\] |
| M1 | Non-thermal intraband kinetics at fixed carrier density | Electron-electron prefactor, deformation potential D, impurity damping (constant rate or Brooks-Herring density), Umklapp fraction of electron-electron damping, band width | \[14, 15, 62\] |
| M2 | M1 plus Auger-like interband relaxation | Rate, gap, valence-band mass | \[22\] |
| M3 | M1 plus interband two-photon absorption with state filling. Open only for pumps below about 930 nm. | Two-photon coefficient, gap, valence-band mass | \[23\] |
| M4 | M2 and M3 together | Both |  |

\[22\] and \[23\] assume different gaps (3.9 and 2.67 eV) and print no valence-band parameters; both are fitted or taken from the authors.

**Damping and cooling are separate parameters.** The linear permittivity fixes only the total momentum-relaxation rate. How it splits between electron-phonon and impurity scattering sets the cooling rate, because only phonons take energy. \[14, 15\] use D = 17.2 eV with a constant impurity rate; \[62\] uses D = 8.6 eV with Brooks-Herring impurities (four times weaker coupling). The prior on D spans at least 8.6-17.2 eV, the impurity model is a discrete choice, and recovery times (\[1, 21\], WP5) constrain the split. These are sample parameters in the hierarchical fit.

**Pipeline**

1. **Identifiability on synthetic data.** Generate data from each model with realistic noise and confirm the pipeline recovers the generating model. This comes before any real data.
2. **Point estimates** by gradient descent through the differentiable solver.
3. **Posteriors.** The forward model is deterministic and differentiable with an explicit noise model, so the likelihood is tractable. A film run takes minutes, which is too slow for long sequential chains; use parallel simulation-based inference with 10^3-10^4 runs \[61\], or Hamiltonian Monte Carlo if the GPU benchmark makes single runs fast enough.
4. **Ranking** by two measures: Bayes factor, and leave-one-condition-out prediction (fit some pulse durations, predict the rest). The second is the headline because it does not depend on priors.
5. **Optimal experimental design.** Choose the pump duration, intensity, angle and probe wavelength that maximise expected disagreement between the two leading models. This defines the new measurement in WP5.

**Data**

| Dataset | Observable | Use | Raw data |
| --- | --- | --- | --- |
| Alam et al. \[1\] | Reflection and transmission against intensity at 30 degrees; n2 and β against angle and wavelength | Calibrates the intraband response | paper not yet read |
| Zhou et al. \[17\] | Probe frequency shift (620 nm stack) | Tests time refraction | on request |
| Tirole et al. \[19\] | Time-diffraction fringes | Sensitive to rise time | on request |
| Gurung et al. \[21\] | Relaxation time against intensity and angle | Tests recovery dynamics | in the paper and supplement |
| Harwood et al. \[22\] | Complex reflection and transmission coefficients in time, 1300 nm | Amplitude and phase at high fluence; tests M2 | on request |
| Lustig et al. \[20\], Ozlu et al. \[23\] | Shifts, sign reversal, 800 nm few-cycle pumps | Tier 2a once validated; the only data that can separate M2 from M3 | on request \[23\] |
| Un et al. \[62\] | Simulated ε(t) under a 5 fs pump, electron-electron and electron-phonon rates (4 nm sphere) | Reproduced by `enzkin.densmat` (Tier 2b); benchmark for Tier 2a; not a measurement | on request |
| Herbst et al., Adv. Sci. 13, e16818 (2026), cited by \[62\] (not yet read) | Transmitted field of ITO nanocrystals under two-cycle pumps, by fieldoscopy | Candidate few-cycle dataset; comparable with `enzkin.densmat` directly (nanocrystal route) | to be checked |
| Choi et al. \[24\] | Simulated transients, public code | Cross-check of M0 | Figshare |
| Own double-blind FROG data | Measured pump and probe fields plus complex gate | Removes assumed-pulse systematics | own |

Samples differ, so the fit is hierarchical: material parameters per sample, mechanism parameters shared.

**Caveat.** The verdict is conditional on the model family. Posterior predictive checks are reported for every model. If none fits, that is the finding, and it points to physics outside the family such as coherent effects.

## WP5: Experimental validation (months 8-11)

WP5 turns the model into a tested prediction, using ENZ-gated double-blind FROG to measure both fields and the complex gate.

**Why this measurement.** The film's transient response acts as a complex gate. The closest precedent is the retrieval of r(t) in \[22\]. Double-blind retrieval also returns the pump and probe fields at the sample, so the simulation is driven by measured pulses instead of assumed ones. \[58\] shows that a learned retrieval needs refinement through a differentiable forward model on the measured trace; \[60\] is an intensity cross-correlation in which phase is lost; \[59\] documents the group's ITO gating set-up. None of them treats double-blind retrieval with a complex gate, so its uniqueness and time resolution have to be established here.

**Measurement matrix**

| Variable | Range | Tests |
| --- | --- | --- |
| Pump peak intensity | Linear regime to just below damage | H1 and the onset of H2 |
| Pump duration at fixed fluence | As the source allows, by chirp or bandwidth | Two-photon against Auger-like (H2); needs a pump below about 930 nm |
| Angle of incidence, TM | Normal incidence to past the ENZ-mode angle | H3 and the trends of \[21\] |
| Film thickness | Two or three values | H3 |
| Probe wavelength | Both sides of the ENZ point | Sign of the response |

A collinear or near-collinear geometry is preferred: a crossing angle smears the gate by the transit of the pump front across the spot \[60\].

**Candidate predictions.** WP4 step 5 picks the final one.

- At fixed fluence, the reversal threshold shifts with pulse duration for two-photon absorption but stays nearly constant for Auger-like excitation.
- At fixed kinetics, the apparent rise time changes with film thickness and angle (H3).
- A second pump arriving inside the non-thermal window gives a non-additive response that a thermal model cannot reproduce.

**Protocol.** Freeze and timestamp the prediction before the measurement is taken.

**Open questions that set the scope**

- [ ] Which pump wavelengths and durations are available, and can duration be varied at fixed fluence? Is an 800 nm pump available?
- [ ] What in-film intensity do the samples survive?
- [ ] Which thicknesses and ENZ wavelengths are on hand, with Hall and ellipsometry data?
- [ ] Do existing FROG datasets already include intensity scans that WP4 can use?
- [ ] Who runs the measurements, and when is the setup free?

## Deferred: learned collision operator (old WP2)

The structure-preserving neural operator is not built now. It comes back if one of these is measured:

- WP4 inference is too slow with the exact solver even on a GPU. An emulator from model parameters to observables may then serve better than a collision-operator surrogate.
- A microscopic few-cycle model appears whose response does not reduce to an electron temperature plus a linear memory, for example one with coupling between coherences or with coherent collision-assisted absorption in a film. Test GT2 (2026-11-24) showed that the density matrix of \[16, 62\] does reduce that way and is cheap, so it gives a learned constitutive law nothing to learn.
- Films with lateral structure need thousands of cells.

What revision 1's literature review established, for that day: \[47\] already enforces one conservation law exactly and \[46\] has an exact equilibrium fixed point, so the novelty would be two weighted invariants, Pauli blocking, a Fermi-Dirac fixed point with exact conservation, and an entropy constraint. \[49\] learns a flow map, not a collision operator. Training through the solver is necessary \[51\]. The surrogate must deliver the collision rate as well as the collision term. The exact index identity that gives conservation needs a uniform energy grid. Details are in `references/LITERATURE_SUMMARY.md`, Part I §8.

## Timeline and gates

| Work package | Months | State |
| --- | --- | --- |
| WP1 Ground-truth kinetic solver | 0-3 | First pass complete; remaining items listed above |
| WP3 Maxwell-kinetic closure | 1-7 | Starts now |
| WP4 Mechanism discrimination | 5-10 | Not started |
| WP5 Experimental validation | 8-11 | Not started |
| Writing: methods paper month 9, flagship letter month 12 |  |  |

| Gate | Month | Pass criterion | If it fails |
| --- | --- | --- | --- |
| G1 | 2 | **Closed, not passed as written.** The exact integral is 6-46 times costlier than the FFT integral at 128-256 points (266 times at 512), and the FFT integral misses the transient by 7-221 %. | Taken: use the exact integral in WP3; drop the surrogate for now |
| G1b | 1 | 128 or 256 points reproduce the finest grid within 1 % of the transient amplitude over the whole condition matrix, including weak excitation and the 6 eV band; a film run on the exact integral fits in hours on a cluster node or minutes on a GPU. **Status 2026-11-17: grid and CPU parts met (128 points per 4 eV; 1.7-12 core-hours per film run for 4-8 eV bands); GPU part outstanding.** | Use the finer grid and recost WP3 and WP4; reopen the surrogate question if a film run exceeds a day |
| GT2 | 1 | Decision test, not a pass/fail gate: does the few-cycle density matrix contain anything a cheap closure misses, and is it too costly to run directly? **Closed 2026-11-24: no on both counts**; the learned operator stays deferred | n/a |
| G3 | 6 | Validation steps 1-3 pass with parameters fixed at step 2 | Revisit the absorption model and the current equation |
| G4 | 8 | Synthetic study recovers the generating model in at least 9 of 10 trials | Redesign the WP5 measurements and narrow the claims to H1 and H3 |

Gate G2 of revision 1 (acceptance of the neural operator) is suspended with WP2.

**Assumed effort.** One researcher at about half time and roughly two weeks of laser time in months 8-11.

**Compute.** A pilot of about 2,000 CPU core-hours and 100 GPU-hours for gate G1b, then about 50,000 core-hours, one GPU with at least 40 GB for 1,000-1,500 hours, and 500 GB of storage. The grid size decided at G1b sets these numbers (working log §1.8).

## Risks and mitigations

| Risk | Likelihood | Mitigation |
| --- | --- | --- |
| The exact solver is not as cheap as the first measurement suggests (finer grid needed, wider band, more cells) | Medium | Gate G1b measures it before WP3 depends on it; GPU port; the deferred operator |
| Ground-truth physics or parameters are wrong | High | Corrections to \[14\] are documented and tested; \[15\] is reproduced only in part; write to the authors; material parameters inferred per sample |
| The current equation cannot be closed at the level of an isotropic distribution | Medium | Carry the l = 1 part of the distribution per cell; report both forms and their observable difference |
| Few-cycle data lie outside Boltzmann validity | Certain for \[20, 23\] | Tier 2a (field-resolved absorption) in the film solver, benchmarked against the particle density matrix (Tier 2b, built); flagship claims limited to pumps of 20 fs and longer until Tier 2a exists. No microscopic few-cycle model for films exists yet. |
| The electron-phonon coupling and the impurity damping are not settled | Certain: the same group uses D = 17.2 eV \[14\] and 8.6 eV \[62\] | Infer them per sample from linear permittivity and recovery times; ask the authors; report how conclusions depend on the split |
| High-intensity data lie outside the single-band model | Certain above about 5 x 10^8 J/m3 absorbed | State the limit; widen the band; add higher bands only if the data demand it |
| Mechanisms cannot be separated with available data | High for Auger against two-photon at 1.1-1.6 um | Restrict that comparison to 800 nm data; synthetic identifiability study first; WP5 designed to break the degeneracy |
| Inference is too slow for sequential sampling | Medium | Parallel simulation-based inference; GPU port; parameter-to-observable emulator |
| Authors do not share raw data | Medium | Digitise figures; weight own FROG data more heavily |
| Samples damage below the reversal threshold | Unknown | Settle with the open questions in WP5; restrict claims to H1 and H3 if so |
| Another group publishes first | Medium | Conference abstract once Gate G3 passes; methods paper to arXiv at month 9 |

**Who could get there first**

- The Sivan group has the kinetic theory, its single-cycle extension \[62\], and names layer geometries as future work \[16, 62\].
- The Scalora group has time-domain hydrodynamic film models \[11, 12\] and reproduced \[20\] with one (Phys. Rev. A 112, 013502 (2025), cited by \[62\]).
- The Purdue groups have few-cycle experiments and their own models \[23, 24\].
- The Imperial group has complex-coefficient data at high fluence \[22\].

None combines kinetics, Maxwell coupling and cross-mechanism inference. All three are also natural collaborators.

## Publication strategy

**Paper A, flagship.** Working title: *What sets the response time of an epsilon-near-zero nonlinearity.*

| Figure | Content |
| --- | --- |
| 1 | Concept: kinetics in every cell of a film, pump and probe propagated together |
| 2 | One parameter set against published transients across pulse duration and intensity |
| 3 | Mechanism ranking, and the rise time split into kinetic and electromagnetic parts |
| 4 | The frozen prediction against the double-blind FROG measurement |

The case for PRL rests on three things: it settles a named disagreement among recent papers \[15, 22-24\]; H3 is conceptually new if confirmed, since it says the measured rise is not the kinetic rise; and the prediction is tested.

**Venue rule**

- A mechanism verdict plus a confirmed prediction goes to PRL.
- If the strongest result is experimental, such as a response time engineered through thickness or angle, it goes to Optica.
- If the mechanisms prove inseparable, the Maxwell-kinetic simulation and H3 go to Optica or Phys. Rev. Applied.
- Further fallbacks: Nature Communications, where \[24\] appeared, and ACS Photonics.

**Paper B, methods.** Working title: *An exactly conserving kinetic solver for hot electrons in transparent conducting oxides.* Content: the discretisation that conserves number and energy to round-off, the corrections to the formulas of \[14\], the reproduction of \[15\] and what it reveals about unprinted parameters, the finding that the random-k baseline \[40\] and the relaxation-time model miss the permittivity transient, and the open code and benchmark dataset. First choice is Phys. Rev. Applied, home of \[14, 15\]; second is Comput. Phys. Commun.

**Sequence.** Paper B to arXiv at month 9, Paper A at month 12. Submit a conference abstract (CLEO or Frontiers in Optics 2027) once Gate G3 passes, to establish priority on the Maxwell-kinetic result.

## References

Entries 1-61 were checked against Crossref and arXiv on Nov 17, 2026, and the numbering is that of revision 1. Entry 62 was added in revision 2.1. Download status, what was read for each paper, and the digests are in `references/README.md` and `references/LITERATURE_SUMMARY.md`.


**ENZ nonlinearity: foundations and models**

1. M. Zahirul Alam, Israel De Leon, Robert W. Boyd, "Large optical nonlinearity of indium tin oxide in its epsilon-near-zero region," Science 352, 795-797 (2016). [doi](https://doi.org/10.1126/science.aae0330). Not obtained (paywalled).
2. L. Caspani, R. P. M. Kaipurath, M. Clerici, et al., "Enhanced nonlinear refractive index in epsilon-near-zero materials," Physical Review Letters 116, 233901 (2016). [arXiv:1603.03581](https://arxiv.org/abs/1603.03581). [doi](https://doi.org/10.1103/physrevlett.116.233901).
3. Orad Reshef, Israel De Leon, M. Zahirul Alam, Robert W. Boyd, "Nonlinear optical effects in epsilon-near-zero media," Nature Reviews Materials 4, 535-551 (2019). [doi](https://doi.org/10.1038/s41578-019-0120-5). Not obtained (paywalled).
4. Nathaniel Kinsey, Clayton DeVault, Alexandra Boltasseva, Vladimir M. Shalaev, "Near-zero-index materials for photonics," Nature Reviews Materials 4, 742-760 (2019). [doi](https://doi.org/10.1038/s41578-019-0133-0). Not obtained (paywalled).
5. Jacob B. Khurgin, Matteo Clerici, Nathaniel Kinsey, "Fast and slow nonlinearities in epsilon-near-zero materials," Laser & Photonics Reviews 15, 2000291 (2020). [arXiv:2007.05613](https://arxiv.org/abs/2007.05613). [doi](https://doi.org/10.1002/lpor.202000291).
6. Ray Secondo, Jacob Khurgin, Nathaniel Kinsey, "Absorptive loss and band non-parabolicity as a physical origin of large nonlinearity in epsilon-near-zero materials," Optical Materials Express 10, 1545 (2020). [arXiv:2004.00134](https://arxiv.org/abs/2004.00134). [doi](https://doi.org/10.1364/ome.394111).
7. Heng Wang, Kang Du, Chuhao Jiang, et al., "Extended Drude model for intraband-transition-induced optical nonlinearity," Physical Review Applied 11, 064062 (2019). [doi](https://doi.org/10.1103/physrevapplied.11.064062). Not obtained (paywalled).
8. R. Secondo, A. Ball, B. Diroll, et al., "Deterministic modeling of hybrid nonlinear effects in epsilon-near-zero thin films," Applied Physics Letters 120, 031103 (2022). [doi](https://doi.org/10.1063/5.0077116).
9. Joshua Baxter, Adriana Pérez-Casanova, Luis Cortes-Herrera, et al., "Understanding the nonlinear optical response of epsilon near zero materials in the time-domain," Advanced Photonics Research 4, 2200280 (2023). [arXiv:2110.14806](https://arxiv.org/abs/2110.14806). [doi](https://doi.org/10.1002/adpr.202200280).
10. Chentao Li, Xinyu Tian, Guoce Yang, et al., "Invertible optical nonlinearity in epsilon-near-zero materials," Physical Review Research 5, 013198 (2023). [arXiv:2203.09600](https://arxiv.org/abs/2203.09600). [doi](https://doi.org/10.1103/physrevresearch.5.013198).
11. Michael Scalora, Jose Trull, Domenico de Ceglia, et al., "Electrodynamics of conductive oxides: intensity-dependent anisotropy, reconstruction of the effective dielectric constant, and harmonic generation," Physical Review A 101, 053828 (2020). [arXiv:2003.01059](https://arxiv.org/abs/2003.01059). [doi](https://doi.org/10.1103/physreva.101.053828).
12. L. Rodríguez-Suné, M. Scalora, A. S. Johnson, et al., "Study of second and third harmonic generation from an indium tin oxide nanolayer: influence of nonlocal effects and hot electrons," APL Photonics 5, 010801 (2020). [arXiv:1910.00313](https://arxiv.org/abs/1910.00313). [doi](https://doi.org/10.1063/1.5129627).
13. Yaraslau Tamashevich, Tornike Shubitidze, Luca Dal Negro, Marco Ornigotti, "Field theory description of the non-perturbative optical nonlinearity of epsilon-near-zero media," APL Photonics 9, 016105 (2024). [arXiv:2308.03181](https://arxiv.org/abs/2308.03181). [doi](https://doi.org/10.1063/5.0171708).

**Kinetic theory of low-electron-density Drude materials**

14. Subhajit Sarkar, Ieng Wai Un, Yonatan Sivan, "Electronic and thermal response of low-electron-density Drude materials to ultrafast optical illumination," Physical Review Applied 19, 014005 (2023). [arXiv:2207.14528](https://arxiv.org/abs/2207.14528). [doi](https://doi.org/10.1103/physrevapplied.19.014005).
15. Ieng-Wai Un, Subhajit Sarkar, Yonatan Sivan, "Electronic-based model of the optical nonlinearity of low-electron-density Drude materials," Physical Review Applied 19, 044043 (2023). [arXiv:2210.08504](https://arxiv.org/abs/2210.08504). [doi](https://doi.org/10.1103/physrevapplied.19.044043).
16. Ieng-Wai Un, Subhajit Sarkar, Yonatan Sivan, "Quantum-optical theory of the few-femtosecond nonlinear optical response of Drude metals with a nonparabolic conduction band," Physical Review B 113, 165411 (2026). [arXiv:2601.10985](https://arxiv.org/abs/2601.10985). [doi](https://doi.org/10.1103/9rs8-8zh8). Revision 1 listed I. W. Un as sole author.

**Time-varying ENZ experiments**

17. Yiyu Zhou, M. Zahirul Alam, Mohammad Karimi, et al., "Broadband frequency translation through time refraction in an epsilon-near-zero material," Nature Communications 11, 2180 (2020). [doi](https://doi.org/10.1038/s41467-020-15682-2).
18. Justus Bohn, Ting Shan Luk, Simon Horsley, Euan Hendry, "Spatiotemporal refraction of light in an epsilon-near-zero indium tin oxide layer: frequency shifting effects arising from interfaces," Optica 8, 1532 (2021). [arXiv:2107.05970](https://arxiv.org/abs/2107.05970). [doi](https://doi.org/10.1364/optica.436324). The arXiv version has three authors (no Horsley) and no subtitle.
19. Romain Tirole, Stefano Vezzoli, Emanuele Galiffi, et al., "Double-slit time diffraction at optical frequencies," Nature Physics 19, 999-1002 (2023). [arXiv:2206.04362](https://arxiv.org/abs/2206.04362). [doi](https://doi.org/10.1038/s41567-023-01993-w).
20. Eran Lustig, Ohad Segal, Soham Saha, et al., "Time-refraction optics with single cycle modulation," Nanophotonics 12, 2221-2230 (2023). [doi](https://doi.org/10.1515/nanoph-2023-0126). PDF not obtained; read from Europe PMC full text.
21. Sudip Gurung, Subhajit Bej, Quynh Dang, et al., "Control of ultrafast hot electron dynamics in epsilon-near-zero conductive oxide thin films," Science Advances 11, eadu8850 (2025). [doi](https://doi.org/10.1126/sciadv.adu8850).
22. Anthony C. Harwood, Sean Z. J. Lim, T. V. Raziman, et al., "Competing Auger-like and intraband excitations drive anomalous modulations in indium tin oxide," Optica 13, 1943 (2026). [arXiv:2605.21432](https://arxiv.org/abs/2605.21432). [doi](https://doi.org/10.1364/optica.607371). Revision 1 cited the arXiv v1 title, "Intraband and interband competition drives ultrafast modulations of indium tin oxide."
23. Mustafa Goksu Ozlu, Colton Fruhling, Ian Hoffman, et al., "Nonlinearity reversal in epsilon-near-zero indium tin oxide driven by few-cycle light pulse," arXiv preprint (2026). [arXiv:2606.11567](https://arxiv.org/abs/2606.11567).
24. Jae Ik Choi, Vahagn Mkhitaryan, Colton Fruhling, et al., "Pathway to optical-cycle dynamic photonics: extreme electron temperatures in transparent conducting oxides," Nature Communications 17, 11129 (2026). [arXiv:2512.24641](https://arxiv.org/abs/2512.24641). [doi](https://doi.org/10.1038/s41467-026-77671-1).
25. Ohad Segal, Noa Konforty, Oded Schiller, et al., "Sub-cycle time-refraction at optical frequencies," CLEO 2026, FM2E.2 (2026). [arXiv:2601.05566](https://arxiv.org/abs/2601.05566). [doi](https://doi.org/10.1364/cleo_fs.2026.fm2e.2).

**Time-varying media: context and constraints**

26. Emanuele Galiffi, Romain Tirole, Shixiong Yin, et al., "Photonics of time-varying media," Advanced Photonics 4, 1-32 (2022). [arXiv:2111.08640](https://arxiv.org/abs/2111.08640). [doi](https://doi.org/10.1117/1.ap.4.1.014002).
27. Zeki Hayran, Jacob B. Khurgin, Francesco Monticone, "hw versus hk: dispersion and energy constraints on time-varying photonic materials and time crystals," Optical Materials Express 12, 3904 (2022). [arXiv:2208.00459](https://arxiv.org/abs/2208.00459). [doi](https://doi.org/10.1364/ome.471672).
28. Soham Saha, Ohad Segal, Colton Fruhling, et al., "Photonic time crystals: a materials perspective," Optics Express 31, 8267 (2023). [doi](https://doi.org/10.1364/oe.479257).
29. Michaël Lobet, Nathaniel Kinsey, Iñigo Liberal, et al., "New horizons in near-zero refractive index photonics and hyperbolic metamaterials," ACS Photonics 10, 3805-3820 (2023). [arXiv:2306.01314](https://arxiv.org/abs/2306.01314). [doi](https://doi.org/10.1021/acsphotonics.3c00747).
30. Diego M. Solís, Raphael Kastner, Nader Engheta, "Time-varying materials in the presence of dispersion: plane-wave propagation in a Lorentzian medium with temporal discontinuity," Photonics Research 9, 1842 (2021). [arXiv:2103.06142](https://arxiv.org/abs/2103.06142). [doi](https://doi.org/10.1364/prj.427368).
31. S. A. R. Horsley, E. Galiffi, Y.-T. Wang, "Eigenpulses of dispersive time-varying media," Physical Review Letters 130, 203803 (2023). [arXiv:2208.11778](https://arxiv.org/abs/2208.11778). [doi](https://doi.org/10.1103/physrevlett.130.203803).

**Non-equilibrium electron kinetics in metals**

32. B. Rethfeld, A. Kaiser, M. Vicanek, G. Simon, "Ultrafast dynamics of nonequilibrium electrons in metals under femtosecond laser irradiation," Physical Review B 65, 214303 (2002). [doi](https://doi.org/10.1103/physrevb.65.214303). Not obtained (paywalled).
33. B. Y. Mueller, B. Rethfeld, "Relaxation dynamics in laser-excited metals under nonequilibrium conditions," Physical Review B 87, 035139 (2013). [doi](https://doi.org/10.1103/physrevb.87.035139). Not obtained (paywalled).
34. E. Carpene, "Ultrafast laser irradiation of metals: beyond the two-temperature model," Physical Review B 74, 024301 (2006). [doi](https://doi.org/10.1103/physrevb.74.024301). Not obtained (paywalled).
35. G. Della Valle, M. Conforti, S. Longhi, G. Cerullo, D. Brida, "Real-time optical mapping of the dynamics of nonthermal electrons in thin gold films," Physical Review B 86, 155139 (2012). [doi](https://doi.org/10.1103/physrevb.86.155139). PDF not obtained; abstract only.
36. Yonatan Dubi, Yonatan Sivan, "'Hot' electrons in metallic nanostructures: non-thermal carriers or heating?," Light: Science & Applications 8, 89 (2019). [arXiv:1810.00565](https://arxiv.org/abs/1810.00565). [doi](https://doi.org/10.1038/s41377-019-0199-x).
37. Ana M. Brown, Ravishankar Sundararaman, Prineha Narang, William A. Goddard, Harry A. Atwater, "Ab initio phonon coupling and optical response of hot electrons in plasmonic metals," Physical Review B 94, 075120 (2016). [arXiv:1602.00625](https://arxiv.org/abs/1602.00625). [doi](https://doi.org/10.1103/physrevb.94.075120).
38. Markus Uehlein, Sebastian T. Weber, Baerbel Rethfeld, "Implementation of the electronic non-equilibrium in the two-temperature model," Nanomaterials 12, 1655 (2022). [arXiv:2112.01851](https://arxiv.org/abs/2112.01851). [doi](https://doi.org/10.3390/nano12101655).
39. M. Uehlein, H. T. Snowden, C. Seibel, et al., "Capturing non-equilibrium electron dynamics in metals accurately and efficiently," Journal of Applied Physics 138, 063103 (2025). [arXiv:2503.09479](https://arxiv.org/abs/2503.09479). [doi](https://doi.org/10.1063/5.0276556).
40. Stephanie Roden, Christopher Seibel, Tobias Held, et al., "Thermalization of optically excited Fermi systems: electron-electron collisions in solid metals," Physics of Plasmas 33, 053906 (2026). [arXiv:2601.11371](https://arxiv.org/abs/2601.11371). [doi](https://doi.org/10.1063/5.0312335).
41. Ivan Maliyov, Jia Yin, Jia Yao, Chao Yang, Marco Bernardi, "Dynamic mode decomposition of nonequilibrium electron-phonon dynamics: accelerating the first-principles real-time Boltzmann equation," npj Computational Materials 10, 123 (2024). [arXiv:2311.07520](https://arxiv.org/abs/2311.07520). [doi](https://doi.org/10.1038/s41524-024-01308-4).

**Neural operators and learned kinetic closures**

42. Zongyi Li, Nikola Kovachki, Kamyar Azizzadenesheli, et al., "Fourier neural operator for parametric partial differential equations," ICLR 2021 (2021). [arXiv:2010.08895](https://arxiv.org/abs/2010.08895).
43. Nikola Kovachki, Zongyi Li, Burigede Liu, et al., "Neural operator: learning maps between function spaces with applications to PDEs," J. Mach. Learn. Res. 24, 89 (2023). [arXiv:2108.08481](https://arxiv.org/abs/2108.08481).
44. Lu Lu, Pengzhan Jin, Guofei Pang, Zhongqiang Zhang, George Em Karniadakis, "Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators," Nature Machine Intelligence 3, 218-229 (2021). [arXiv:1910.03193](https://arxiv.org/abs/1910.03193). [doi](https://doi.org/10.1038/s42256-021-00302-5). The arXiv version has three authors and a different title.
45. Zongyi Li, Hongkai Zheng, Nikola Kovachki, et al., "Physics-informed neural operator for learning partial differential equations," ACM / IMS Journal of Data Science 1, 1-27 (2024). [arXiv:2111.03794](https://arxiv.org/abs/2111.03794). [doi](https://doi.org/10.1145/3648506).
46. Tianbai Xiao, Martin Frank, "RelaxNet: a structure-preserving neural network to approximate the Boltzmann collision operator," Journal of Computational Physics 490, 112317 (2023). [arXiv:2211.08149](https://arxiv.org/abs/2211.08149). [doi](https://doi.org/10.1016/j.jcp.2023.112317).
47. Jae Yong Lee, Steffen Schotthöfer, Tianbai Xiao, Sebastian Krumscheid, Martin Frank, "Structure-preserving operator learning: modeling the collision operator of kinetic equations," arXiv preprint (2024). [arXiv:2402.16613](https://arxiv.org/abs/2402.16613).
48. Jae Yong Lee, Gwang Jae Jung, Byungchan Lim, Hyung Ju Hwang, "FourierSpecNet: neural collision operator approximation inspired by the Fourier spectral method for solving the Boltzmann equation," Kinetic and Related Models 21, 69-103 (2026). [arXiv:2504.20408](https://arxiv.org/abs/2504.20408). [doi](https://doi.org/10.3934/krm.2026004).
49. Boyun Hu, Kunlun Qi, "Solving the BGK model and Boltzmann equation by Fourier neural operator with conservative constraints," arXiv preprint (2025). [arXiv:2510.13047](https://arxiv.org/abs/2510.13047).
50. Yue Zhao, Joshua Burby, Andrew Christlieb, Huan Lei, "Data-driven construction of a generalized kinetic collision operator from molecular dynamics," Physical Review Letters 135, 185101 (2025). [arXiv:2503.24208](https://arxiv.org/abs/2503.24208). [doi](https://doi.org/10.1103/gml7-f2ng).
51. Samuel Burles, Enrico Camporeale, Oreste Pezzi, "A neural operator closure for Landau damping in electrostatic plasma," Journal of Plasma Physics 92, E149 (2026). [arXiv:2607.29364](https://arxiv.org/abs/2607.29364). [doi](https://doi.org/10.1017/s0022377826102517).
52. M. Luo, A. R. Bell, F. Miniati, S. M. Vinko, G. Gregori, "Resolution-independent machine learning heat flux closure for ICF plasmas," PRX Intelligence 1, 013017 (2026). [arXiv:2604.03439](https://arxiv.org/abs/2604.03439). [doi](https://doi.org/10.1103/9l4n-mnz6).
53. Karan Shah, Attila Cangi, "Accelerating electron dynamics simulations through machine learned time propagators," ICML 2024 AI4Science Workshop (2024). [arXiv:2407.09628](https://arxiv.org/abs/2407.09628).
54. Karan Shah, Attila Cangi, "Machine learning time propagators for time-dependent density functional theory simulations," Machine Learning: Science and Technology 7, 035019 (2026). [arXiv:2508.16554](https://arxiv.org/abs/2508.16554). [doi](https://doi.org/10.1088/2632-2153/ae68d4).

**Neural operators in photonics; machine learning for ultrafast metrology**

55. Jiaqi Gu, Zhengqi Gao, Chenghao Feng, et al., "NeurOLight: a physics-agnostic neural operator enabling parametric photonic device simulation," Advances in Neural Information Processing Systems 35, 14623-14636 (2022). [arXiv:2209.10098](https://arxiv.org/abs/2209.10098). [doi](https://doi.org/10.52202/068431-1063).
56. Pingchuan Ma, Haoyu Yang, Zhengqi Gao, Duane S. Boning, Jiaqi Gu, "PIC2O-Sim: a physics-inspired causality-aware dynamic convolutional neural operator for ultra-fast photonic device FDTD simulation," APL Photonics 10, 036104 (2025). [arXiv:2406.17810](https://arxiv.org/abs/2406.17810). [doi](https://doi.org/10.1063/5.0242728).
57. B. L. Blackwood, S. Hawley (per the related Belmont University records), "Physics-informed Fourier neural operator for nanophotonic scattering in nonlinear media," poster, APS Global Physics Summit 2026 (https://meetings-archive.aps.org/smt/2026/mar-h00-01/159/). No paper exists. Related Belmont University records by B. L. Blackwood and S. Hawley carry the abstract.
58. Tom Zahavy, Alex Dikopoltsev, Daniel Moss, et al., "Deep learning reconstruction of ultrashort pulses," Optica 5, 666 (2018). [arXiv:1803.06024](https://arxiv.org/abs/1803.06024). [doi](https://doi.org/10.1364/optica.5.000666). The arXiv version is the four-page workshop paper with five authors.
59. Hao Zhang, Yang Xu, Wenwen Zhang, et al., "Hybrid deep reconstruction for vignetting-free upconversion imaging through scattering in ENZ materials," Light: Science & Applications 15, 327 (2026). [arXiv:2508.13096](https://arxiv.org/abs/2508.13096). [doi](https://doi.org/10.1038/s41377-026-02375-6).
60. Adam Ball, Ray Secondo, Dhruv Fomra, et al., "A space-time knife-edge in epsilon-near-zero films for ultrafast pulse characterization," Laser & Photonics Reviews 19, 2401462 (2024). [arXiv:2408.00577](https://arxiv.org/abs/2408.00577). [doi](https://doi.org/10.1002/lpor.202401462).

**Inference**

61. Kyle Cranmer, Johann Brehmer, Gilles Louppe, "The frontier of simulation-based inference," Proceedings of the National Academy of Sciences 117, 30055-30062 (2020). [arXiv:1911.01429](https://arxiv.org/abs/1911.01429). [doi](https://doi.org/10.1073/pnas.1912789117).

**Added in revision 2.1**

62. Ieng-Wai Un, Subhajit Sarkar, Yonatan Sivan, "Origin, strength, and speed of the nonlinear optical response of transparent conducting oxides to single-cycle optical pulses," Nanophotonics 15, e70272 (2026). [arXiv:2608.11756](https://arxiv.org/abs/2608.11756). [doi](https://doi.org/10.1002/nap2.70272). Open access (CC BY 4.0).
