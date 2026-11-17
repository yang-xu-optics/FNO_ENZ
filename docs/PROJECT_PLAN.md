# Learned Kinetic Closure for ENZ Electron Dynamics: Project Plan

Nov 16, 2026 · @Yang

## Summary

The project answers one physics question: what sets the speed and sign of the optical response of epsilon-near-zero (ENZ) indium tin oxide at high intensity and short pulse duration? Four groups currently give four different answers to the same anomalies, and no existing model couples distribution-resolved electron kinetics to Maxwell's equations in a film.

The method is a structure-preserving neural operator that learns the electron collision integral from full Boltzmann solutions. It is embedded in a time-domain Maxwell solver as a constitutive closure, then used as a differentiable forward model to rank the competing mechanisms against data.

**What is new**

- The first neural operator for a degenerate-fermion (Uehling-Uhlenbeck) collision integral with Pauli blocking, electron-phonon coupling and optical driving.
- The first film-scale Maxwell-kinetic simulation of ENZ pump-probe experiments, with no adiabatic approximation for the probe.
- A quantitative, Bayesian verdict on four published mechanisms, plus one prediction tested with double-blind FROG.

**Outputs**

| Output | Target venue | Month |
| --- | --- | --- |
| Flagship letter: origin of the ENZ response time and sign reversal | PRL or Optica | 12 |
| Methods paper: operator, solver, open dataset and code | PRX Intelligence or Phys. Rev. Applied | 10 |

The plan has four go/no-go gates. The most important is at month 2: if the exact kinetic solver is already fast enough, the machine-learning component shrinks and the Maxwell-kinetic physics carries the paper alone.

## Background and gap

ITO near its ENZ wavelength gives order-unity index changes in under a picosecond \[1, 2\], but theory only explains the low-intensity, long-pulse regime. Bracketed numbers refer to the reference list at the end.

**Established picture.** Sub-bandgap light heats conduction electrons. In a non-parabolic band the hot population is heavier on average, so the plasma frequency drops and the ENZ point redshifts \[5-8\]. Four levels of theory exist:

- Two-temperature and extended two-temperature models \[6-9\].
- Hydrodynamic Maxwell models with nonlocal and hot-electron terms \[11, 12\].
- Boltzmann-type kinetics in energy space, which attributes the real-permittivity change to non-thermal electrons \[14, 15\].
- An energy-space density-matrix model for few-femtosecond pulses \[16\].

**Anomalies that these do not jointly explain**

| Observation | Source | Mechanism proposed |
| --- | --- | --- |
| Rise time of 1-10 fs inferred under a 225 fs pump | Tirole et al. \[19\] | None at the time; non-thermal electrons proposed later \[15\] |
| Fast rise and fast recovery under near-single-cycle pumping | Lustig et al. \[20\] | Open |
| Sign reversal of the response above about 5 TW/cm2, sub-8 fs pump | Ozlu et al. \[23\] | Interband two-photon absorption once intraband excitation lifts Pauli blocking |
| Non-monotonic, oscillatory complex Fresnel coefficients at high fluence | Harwood et al. \[22\] | Auger-type interband scattering added to the two-temperature model |
| Oscillatory, sign-reversing transmittance over a few optical cycles | Choi et al. \[24\] | Extreme electron temperatures |
| Relaxation time up by more than 200% with intensity and 40% with ENZ-mode excitation | Gurung et al. \[21\] | Hydrodynamic model |
| Larger time-refraction shift as the index change becomes sub-cycle | Segal et al. \[25\] | Not modelled quantitatively in the abstract |

**Three gaps**

1. **No distribution-resolved kinetic model is coupled to Maxwell in a film.** Temperature-level and hydrodynamic models have been run in time-domain Maxwell solvers \[9, 11, 12\], but they do not resolve the electron distribution. The kinetic models use a time-dependent permittivity in the adiabatic approximation \[15\]. The density-matrix model treats a 4 nm sphere in the quasi-static limit and leaves films and propagation to future work \[16\]. The experiments use 20-310 nm films, often at oblique incidence, where field enhancement, depth-dependent absorption and probe dispersion all shape the measured transient.
2. **The mechanisms have never been tested against each other.** Each group fits its own data with its own model. No study runs several mechanisms through one forward model on a shared dataset.
3. **Operator learning has not reached this physics.** Fourier neural operators exist for the classical Boltzmann collision operator \[48, 49\], with structure-preserving variants \[46, 47\]. Neural operators also exist for TDDFT density propagation \[53, 54\], Maxwell solvers \[55-57\] and plasma closures \[51, 52\]. None treats degenerate electrons with Pauli blocking under optical driving.

**Why it matters.** Optical time refraction, time diffraction and proposed photonic time crystals all rely on these oxides \[17-19, 26, 28\], and the energy cost of fast index modulation is severe \[27\]. Whether sub-cycle modulation is feasible depends on which time constant limits the response and what sets it.

## Central question and hypotheses

**Question.** What sets the rise time, recovery and sign of the ENZ response in ITO as pump intensity rises from 1 GW/cm2 to several TW/cm2 and pulse duration falls from 200 fs toward a few optical cycles?

**Hypotheses**

|  | Hypothesis | Falsified if |
| --- | --- | --- |
| H1 | Below the reversal threshold, intraband non-thermal kinetics at fixed carrier density explain rise and recovery. | One parameter set cannot fit transients across pulse durations. |
| H2 | Sign reversal requires a channel that changes carrier density. Two-photon absorption \[23\] and Auger-type excitation \[22\] are distinguishable by how they scale. | A fixed-density model reproduces the sign flip, or the two channels fit all data equally well. |
| H3 | Part of the apparent speed is electromagnetic. The probe crosses a time-varying, dispersive, depth-inhomogeneous film, so the measured rise is not the kinetic rise. | Maxwell-kinetic and local adiabatic models agree within experimental error in all conditions. |

Expected discriminators for H2, to be confirmed by the model: two-photon absorption follows the instantaneous intensity squared, so at fixed fluence it strengthens as the pulse shortens. Auger-type excitation follows the hot-electron population, so it tracks absorbed energy and lags the pulse.

H3 has not been tested. Temperature-level models have been propagated through films \[9\], but no model that resolves the electron distribution has.

**Target claims for the flagship paper**

1. One parameter set reproduces published transients for pump durations of 20-225 fs and intensities from 1 GW/cm2 to the reversal threshold, within the stated experimental uncertainty.
2. The four mechanisms are ranked by held-out predictive score and Bayes factor.
3. The measured rise time is decomposed into a kinetic part and an electromagnetic part.
4. One prediction chosen by optimal experimental design is confirmed or refuted with double-blind FROG.

**Scope limit.** The Boltzmann description is cycle-averaged. One optical cycle at 1.2 um is 4 fs, so it is reliable for pumps longer than about 20 fs. The sub-8 fs experiments \[20, 23\] need the coherent density-matrix level \[16\]. That is a stretch goal (Tier 2 below), and the flagship claims do not depend on it.

## WP1: Ground-truth kinetic solver (months 0-3)

WP1 builds the reference solver that generates the training data and anchors every later claim. Its first job is to measure whether a learned surrogate is needed at all.

**State and equation.** Each spatial cell carries an isotropic conduction-band distribution f(E, t) in a Kane band, a lattice temperature, and a carrier density n for the models with interband channels.

```latex
\frac{\partial f}{\partial t} = \left.\frac{\partial f}{\partial t}\right|_{\mathrm{abs}} + Q_{ee}[f] + Q_{e\text{-}ph}[f, T_l] + S_{\mathrm{inter}}[f, n, I], \qquad \frac{\hbar^2 k^2}{2 m_0^*} = E\,(1 + C E)
```

The optical response follows from f without assuming a temperature:

```latex
\omega_p^2(t) = \frac{e^2}{3\varepsilon_0} \int g(E)\, v^2(E) \left(-\frac{\partial f}{\partial E}\right) dE
```

| Term | Treatment | Basis |
| --- | --- | --- |
| Intraband absorption | Phonon- and impurity-assisted absorption in steps of the photon energy, driven by the local field from WP3 | \[14, 32\] |
| Electron-electron | Full collision integral with Pauli blocking and a screened Coulomb interaction | \[15, 32, 40\] |
| Electron-phonon | Energy exchange with the lattice, with momentum conservation and Umklapp treated as in \[14\] | \[14\] |
| Interband channels | Two-photon absorption with state filling \[23\]; Auger-type excitation and recombination \[22\]; switched per model | \[22, 23\] |

**Two fidelity levels for the electron-electron term**

- **B0, random-k with a constant matrix element** \[40\]. Energy conservation should reduce the integral to convolutions over energy, which FFTs evaluate in O(N log N) on a uniform grid. This needs checking, but if it holds B0 is a fast, exact, non-ML baseline.
- **B1, momentum-resolved screened interaction**, with screening recomputed from f. The kernel does not factorise and costs O(N^3) per evaluation. This is the surrogate's real target.

**Material parameters.** Start from the parameter set of \[14, 15\] so results compare directly. Approximate ranges to sample: electron density 0.5-1.5 x 10^21 cm^-3, band-edge mass 0.3-0.4 electron masses, non-parabolicity near 0.4 eV^-1. Fix each sample's values from ellipsometry and Hall data.

**Dataset**

- Pump wavelength 1.1-1.6 um, duration 20-250 fs, in-film peak intensity 0.1 GW/cm2 to 10 TW/cm2, with chirped and double pulses.
- Material parameters varied by 25% around nominal, so the operator is parametric and usable for inference.
- About 2 x 10^4 trajectories on 512 energy points over 2 ps, plus perturbed off-trajectory states. Roughly 40 GB at 2 fs sampling.

**Validation**

- Reproduce the published electron-temperature and permittivity transients of \[15\] for its 220 fs and 30 fs cases.
- Conserve particle number and total energy to 10^-6; relax to Fermi-Dirac; recover the extended two-temperature limit \[38\].
- Compare against the reduced athermal model of \[39\] as a second non-ML baseline.

**Tier 2 (stretch).** Implement the energy-space density matrix of \[16\] for few-cycle pumps. Its state is two-dimensional in energy, which suits a 2D Fourier operator.

**Deliverable.** An open solver and dataset. Contacting the authors of \[14-16\] early is the fastest route to a verified ground truth and lowers scooping risk.

## WP2: Structure-preserving neural operator (months 2-6)

WP2 learns the collision operator Q\[f\], not the solution map. That keeps the known excitation term analytic and lets conservation laws be enforced exactly.

A Fourier architecture has the right inductive bias here. Energy conservation gives collision integrals a convolution structure in energy, the same property that makes spectral Boltzmann solvers work \[48\].

| Choice | Decision | Reason |
| --- | --- | --- |
| What is learned | The right-hand sides for electron-electron and electron-phonon scattering, integrated by a standard time stepper | Free step size, exact drive term, interpretable output |
| Input | Deviation of f from a Fermi-Dirac function matched in density and energy, plus material parameters and lattice temperature | Removes the sharp Fermi edge from what the network must represent |
| Energy axis | Mapped grid clustered near the Fermi level, with even extension (cosine transform) | The axis is non-periodic; avoids Gibbs ringing |
| Hard constraints | A projection that removes any number- or energy-violating component; output defined as N\[f\] - N\[f\_FD\] so equilibrium is an exact fixed point | Builds on \[46, 47, 49\], which use soft or architectural constraints |
| Soft constraints | Penalty on negative entropy production; limiter keeping f between 0 and 1 | H-theorem and Pauli principle |
| Training | Supervised on (f, Q) pairs, then unrolled multi-step rollouts with noise injection | Prevents drift off the training manifold |
| Uncertainty | Five-member ensemble; disagreement triggers fallback to the exact integral | Makes the closure trustworthy inside Maxwell |

**Acceptance criteria**

- Relative L2 error on Q below 1% on held-out trajectories, and below 3% for intensities 10% beyond the training range.
- Over a 2 ps rollout, error on plasma frequency and damping below 1%.
- Conservation of number and energy at 10^-6, by construction.
- At least 100x speed-up per evaluation over B1 at 512 energy points, batched on one GPU.
- Trained at 256 points and evaluated at 512 with less than twofold error growth.

**Baselines the paper must report**

- Relaxation-time approximation with energy-dependent lifetime \[40\].
- Extended two-temperature model \[38\] and the athermal electron model \[39\].
- FFT-evaluated B0.
- DeepONet \[44, 47\] and a fixed-grid MLP.

If the operator does not beat the best non-ML baseline on accuracy per unit cost, the paper says so and uses that baseline in WP3.

**Ablation.** An autoregressive flow-map operator, f(t) to f(t + dt), in the style of \[49, 53\]. It is faster but conserves only approximately.

**Tooling.** JAX or PyTorch, so gradients flow end to end through WP3 and WP4. One GPU is enough for training.

## WP3: Maxwell-kinetic closure (months 4-8)

WP3 embeds the operator in a one-dimensional time-domain Maxwell solver. Pump and probe both propagate through a film in which every cell evolves kinetically, with no adiabatic approximation.

**Solver**

- 1D FDTD through the film thickness. Normal incidence first, then TM oblique incidence, which is needed for ENZ-mode excitation \[21\].
- Films of 20-310 nm on a substrate, with 2-5 nm cells. That gives roughly 100 kinetic cells.
- Multirate stepping: fields at about 5 as, kinetics at 0.1-0.25 fs.
- The cycle-averaged absorbed power J.E in each cell drives that cell's absorption term.

**Constitutive law.** The current in each cell follows from the momentum moment of the kinetic equation:

```latex
\frac{\partial J}{\partial t} + \left(\gamma[f] + \frac{\dot{m}^*[f]}{m^*[f]}\right) J = \varepsilon_0\, \omega_p^2[f]\, E
```

A changing effective mass and a changing carrier density give different current equations. A permittivity-only description cannot tell them apart, and the dispersive time-varying-media literature shows the distinction matters \[30, 31\]. Here the kinetics decide.

**Why the surrogate is decisive here.** One run needs about 100 cells x 10^4 kinetic steps = 10^6 collision evaluations. At an assumed 0.1-1 s per B1 evaluation that is 1-10 days per run; with the batched operator it should be minutes. WP4 needs 10^3-10^4 runs. These are order-of-magnitude estimates, to be measured at Gate 1.

**Validation ladder**

1. Linear limit matches transfer-matrix reflection and transmission, including ENZ-mode absorption.
2. Intensity-dependent reflection and transmission of Alam et al. \[1\].
3. Time-refraction frequency shifts of Zhou et al. \[17\] and spatiotemporal refraction of Bohn et al. \[18\].
4. Relaxation-time trends with intensity and angle of Gurung et al. \[21\].

Parameters may be tuned only at step 2. Steps 3 and 4 are predictions. The FDTD thermal model of \[9\] is the comparison baseline throughout.

**Test of H3.** Run identical kinetics two ways: full Maxwell-kinetic, and local adiabatic permittivity with thin-film Fresnel formulas. Compare the probe transients and report the difference in apparent rise time against film thickness, angle and probe bandwidth.

**Deferred.** Nonlocal hydrodynamic terms \[11, 12\] and lateral structure. Add the nonlocal correction if step 1 misplaces the ENZ resonance in the thinnest films, where \[12\] reports a blueshift of tens of nanometres.

**Deliverable.** The first film-scale Maxwell-kinetic simulation of an ENZ pump-probe experiment.

## WP4: Mechanism discrimination (months 6-10)

WP4 runs the published mechanisms through one forward model and lets the data rank them.

**Model family**

| Model | Physics | Extra free parameters | Source |
| --- | --- | --- | --- |
| M0 | Thermal electrons only: extended two-temperature model with full Fermi integrals and temperature-dependent chemical potential. Covers the extreme-temperature explanation. | Coupling constants | \[9, 24, 38\] |
| M1 | Non-thermal intraband kinetics at fixed carrier density | Scattering prefactors, Umklapp ratio | \[14, 15\] |
| M2 | M1 plus Auger-type interband excitation and recombination | Auger coefficient, gap | \[22\] |
| M3 | M1 plus interband two-photon absorption with state filling | Two-photon coefficient, gap | \[23\] |
| M4 | M2 and M3 together | Both |  |

The operator is conditioned on carrier density, so M1-M4 share one trained network. Only the cheap rate terms differ.

**Pipeline**

1. **Identifiability on synthetic data.** Generate data from each model with realistic noise and confirm the pipeline recovers the generating model. This comes before any real data.
2. **Point estimates** by gradient descent through the differentiable solver.
3. **Posteriors** by Hamiltonian Monte Carlo or neural posterior estimation \[61\].
4. **Ranking** by two measures: Bayes factor, and leave-one-condition-out prediction (fit some pulse durations, predict the rest). The second is the headline because it does not depend on priors.
5. **Optimal experimental design.** Choose the pump duration, intensity, angle and probe wavelength that maximise expected disagreement between the two leading models. This defines the new measurement in WP5.

**Data**

| Dataset | Observable | Use |
| --- | --- | --- |
| Alam et al. \[1\] | Reflection and transmission against intensity and angle | Calibrates the intraband response |
| Zhou et al. \[17\] | Probe frequency shift | Tests time refraction |
| Tirole et al. \[19\] | Time-diffraction fringes | Sensitive to rise time |
| Gurung et al. \[21\] | Relaxation time against intensity and angle | Tests recovery dynamics |
| Harwood et al. \[22\] | Complex reflection and transmission coefficients in time | Amplitude and phase at high fluence |
| Ozlu et al. \[23\] | Sign reversal | Tier 2 only (sub-8 fs pump) |
| Own double-blind FROG data | Measured pump and probe fields plus complex gate | Removes assumed-pulse systematics |

Request raw data from the authors and digitise figures where that fails. Samples differ, so the fit is hierarchical: material parameters per sample, mechanism parameters shared.

**Caveat.** The verdict is conditional on the model family. Posterior predictive checks are reported for every model. If none fits, that is the finding, and it points to physics outside the family such as coherent effects.

## WP5: Experimental validation (months 8-11)

WP5 turns the model into a tested prediction, using ENZ-gated double-blind FROG to measure both fields and the complex gate.

**Why this measurement.** The film's transient response acts as a complex gate, the quantity retrieved in \[22\]. Double-blind retrieval also returns the pump and probe fields at the sample, so the simulation is driven by measured pulses instead of assumed ones. Related ENZ pulse metrology is in \[58-60\].

**Measurement matrix**

| Variable | Range | Tests |
| --- | --- | --- |
| Pump peak intensity | Linear regime to just below damage | H1 and the onset of H2 |
| Pump duration at fixed fluence | As the source allows, by chirp or bandwidth | Two-photon against Auger-type (H2) |
| Angle of incidence, TM | Normal incidence to past the ENZ-mode angle | H3 and the trends of \[21\] |
| Film thickness | Two or three values | H3 |
| Probe wavelength | Both sides of the ENZ point | Sign of the response \[10\] |

**Candidate predictions.** WP4 step 5 picks the final one.

- At fixed fluence, the sign-reversal threshold shifts with pulse duration for two-photon absorption but stays nearly constant for Auger-type excitation.
- At fixed kinetics, the apparent rise time changes with film thickness and angle (H3).
- A second pump arriving inside the non-thermal window gives a non-additive response that a thermal model cannot reproduce.

**Protocol.** Freeze and timestamp the prediction before the measurement is taken.

**Open questions that set the scope**

- [ ] Which pump durations are available, and can duration be varied at fixed fluence?
- [ ] What in-film intensity do the samples survive?
- [ ] Which thicknesses and ENZ wavelengths are on hand, with Hall and ellipsometry data?
- [ ] Do existing FROG datasets already include intensity scans that WP4 can use?
- [ ] Who runs the measurements, and when is the setup free?

## Timeline and gates

The project runs 12 months, with work packages overlapping by one to two months and a go/no-go gate inside each of the first four.

&#91;embedded content: project timeline · 5 work packages, 4 gates, 2 submissions\]

Each gate sits where its work package first produces a number that can change the plan.

| Gate | Month | Pass criterion | If it fails |
| --- | --- | --- | --- |
| G1 | 2 | B1 costs at least 100x more per evaluation than B0, and B0 misses B1's permittivity transient by more than 2% | Use FFT-evaluated B0 inside WP3, drop the operator from the flagship, and keep ML only where inference needs speed |
| G2 | 5 | All WP2 acceptance criteria met | Hybrid solver with exact fallback, or the flow-map variant |
| G3 | 7 | Validation steps 1-3 pass with parameters fixed at step 2 | Add nonlocal terms and revisit the absorption model |
| G4 | 9 | Synthetic study recovers the generating model in at least 9 of 10 trials | Redesign the WP5 measurements and narrow the claims to H1 and H3 |

**Assumed effort.** One researcher at about half time, one GPU, and roughly two weeks of laser time in months 9-11. Adjust once the open questions in WP5 are answered.

## Risks and mitigations

The largest risk is that the exact solver proves fast enough to make the surrogate unnecessary. The plan is built so the physics paper survives that outcome.

| Risk | Likelihood | Mitigation |
| --- | --- | --- |
| Exact kinetics are cheap enough without ML | Medium | Gate 1 measures it. Fall back to the FFT baseline inside Maxwell and keep ML only where inference needs speed. |
| Ground-truth physics or parameters are wrong | High | Parametric operator; material parameters inferred per sample; several model families compared; collaboration with the authors of \[14-16\]. |
| Few-cycle data lie outside Boltzmann validity | Certain for \[20, 23\] | Flagship claims limited to pumps of 20 fs and longer. Tier 2 is a stretch goal. |
| Operator drifts in rollout or out of distribution | Medium | Unrolled training, ensemble disagreement check, fallback to the exact integral. |
| Mechanisms cannot be separated with available data | Medium | Synthetic identifiability study first; WP5 measurements designed to break the degeneracy. |
| Authors do not share raw data | Medium | Digitise figures; weight own FROG data more heavily. |
| Samples damage below the reversal threshold | Unknown | Settle with the open questions in WP5; restrict claims to H1 and H3 if so. |
| Referees discount an ML-led paper | Medium | Lead with the physics question. Put the operator in methods and a companion paper. Report every non-ML baseline. |
| Another group publishes first | Medium | Post the methods paper to arXiv at month 10, sooner if Gate 2 passes early. |

**Who could get there first**

- The Sivan group has the kinetic theory and names film geometries as future work \[16\].
- The Purdue groups have few-cycle experiments and their own models \[23, 24\].
- The Imperial group has complex-coefficient data at high fluence \[22\].
- Shah and Cangi have Fourier operators for electron dynamics, though for TDDFT in molecules \[53, 54\].

None combines kinetics, Maxwell coupling and cross-mechanism inference. The first three are also natural collaborators.

## Publication strategy

Two papers: a physics letter whose abstract never needs the word "surrogate", and a methods paper that carries the machine learning.

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

**Paper B, methods.** Working title: *A structure-preserving neural operator for the degenerate-electron collision integral as a Maxwell closure.* First choice is PRX Intelligence, which published a Fourier-operator plasma closure \[52\]. Second choice is Phys. Rev. Applied, home of \[14, 15\]. It releases the solver, dataset and trained weights.

**Sequence.** Paper B to arXiv at month 10, Paper A at month 12. Submit a conference abstract (CLEO or Frontiers in Optics 2027) once Gate 3 passes, to establish priority on the Maxwell-kinetic result.

## References

Entries marked † were cited from memory and not re-opened while drafting. Check them against the publisher record before use. All others were opened or confirmed in a source during drafting.

**ENZ nonlinearity: foundations and models**

1. M. Z. Alam, I. De Leon, R. W. Boyd, "Large optical nonlinearity of indium tin oxide in its epsilon-near-zero region," Science 352, 795 (2016).
2. L. Caspani et al., "Enhanced nonlinear refractive index in epsilon-near-zero materials," Phys. Rev. Lett. 116, 233901 (2016).
3. O. Reshef, I. De Leon, M. Z. Alam, R. W. Boyd, "Nonlinear optical effects in epsilon-near-zero media," Nat. Rev. Mater. 4, 535 (2019).
4. N. Kinsey, C. DeVault, A. Boltasseva, V. M. Shalaev, "Near-zero-index materials for photonics," Nat. Rev. Mater. 4, 742 (2019).
5. J. B. Khurgin, M. Clerici, N. Kinsey, "Fast and slow nonlinearities in epsilon-near-zero materials," Laser Photonics Rev. 15, 2000291 (2021).
6. R. Secondo, J. Khurgin, N. Kinsey, "Absorptive loss and band non-parabolicity as a physical origin of large nonlinearity in epsilon-near-zero materials," Opt. Mater. Express 10, 1545 (2020).
7. H. Wang et al., "Extended Drude model for intraband-transition-induced optical nonlinearity," Phys. Rev. Applied 11, 064062 (2019).
8. R. Secondo et al., "Deterministic modeling of hybrid nonlinear effects in epsilon-near-zero thin films," Appl. Phys. Lett. 120, 031103 (2022). [doi](https://doi.org/10.1063/5.0077116)
9. J. Baxter et al., "Understanding the nonlinear optical response of epsilon near zero materials in the time-domain," [arXiv:2110.14806](https://arxiv.org/abs/2110.14806); doi 10.1002/adpr.202200280.
10. C. Li et al., "Invertible optical nonlinearity in epsilon-near-zero materials," [arXiv:2203.09600](https://arxiv.org/abs/2203.09600) (2022).
11. M. Scalora et al., "Electrodynamics of conductive oxides: intensity-dependent anisotropy, reconstruction of the effective dielectric constant, and harmonic generation," Phys. Rev. A 101, 053828 (2020).
12. L. Rodriguez-Sune et al., "Study of second and third harmonic generation from an indium tin oxide nanolayer: influence of nonlocal effects and hot electrons," APL Photonics 5, 010801 (2020).
13. Y. Tamashevich, T. Shubitidze, L. Dal Negro, M. Ornigotti, "Field theory description of the non-perturbative optical nonlinearity of epsilon-near-zero media," [arXiv:2308.03181](https://arxiv.org/abs/2308.03181) (2023).

**Kinetic theory of low-electron-density Drude materials**

14. S. Sarkar, I. W. Un, Y. Sivan, "Electronic and thermal response of low-electron-density Drude materials to ultrafast optical illumination," Phys. Rev. Applied 19, 014005 (2023).
15. I. W. Un, S. Sarkar, Y. Sivan, "Electronic-based model of the optical nonlinearity of low-electron-density Drude materials," Phys. Rev. Applied 19 (2023); [arXiv:2210.08504](https://arxiv.org/abs/2210.08504v2).
16. I. W. Un, "Quantum-optical theory of the few-femtosecond nonlinear optical response of Drude metals with a nonparabolic conduction band," Phys. Rev. B 113, 165411 (2026). [pdf](https://tzin.bgu.ac.il/~sivanyon/Papers/LEDDs/single_cycle_article.pdf)

**Time-varying ENZ experiments**

17. Y. Zhou et al., "Broadband frequency translation through time refraction in an epsilon-near-zero material," Nat. Commun. 11, 2180 (2020). †
18. J. Bohn, T. S. Luk, S. Horsley, E. Hendry, "Spatiotemporal refraction of light in an epsilon-near-zero indium tin oxide layer," Optica 8, 1532 (2021). †
19. R. Tirole et al., "Double-slit time diffraction at optical frequencies," Nat. Phys. 19, 999 (2023). †
20. E. Lustig et al., "Time-refraction optics with single cycle modulation," Nanophotonics 12, 2221 (2023). [doi](https://www.degruyterbrill.com/document/doi/10.1515/nanoph-2023-0126/html)
21. S. Gurung et al., "Control of ultrafast hot electron dynamics in epsilon-near-zero conductive oxide thin films," Sci. Adv. 11, eadu8850 (2025). [doi](https://www.science.org/doi/10.1126/sciadv.adu8850)
22. A. C. Harwood et al., "Intraband and interband competition drives ultrafast modulations of indium tin oxide," [arXiv:2605.21432](https://arxiv.org/abs/2605.21432) (2026).
23. M. G. Ozlu et al., "Nonlinearity reversal in epsilon-near-zero indium tin oxide driven by few-cycle light pulse," [arXiv:2606.11567](https://arxiv.org/abs/2606.11567) (2026).
24. J. I. Choi et al., "Pathway to optical-cycle dynamic photonics: extreme electron temperatures in transparent conducting oxides," Nat. Commun. (2026), doi 10.1038/s41467-026-77671-1; [arXiv:2512.24641](https://arxiv.org/abs/2512.24641).
25. O. Segal et al., "Sub-cycle time-refraction at optical frequencies," [arXiv:2601.05566](https://arxiv.org/abs/2601.05566) (2026).

**Time-varying media: context and constraints**

26. E. Galiffi et al., "Photonics of time-varying media," Adv. Photonics 4, 014002 (2022); [arXiv:2111.08640](https://arxiv.org/pdf/2111.08640).
27. Z. Hayran, J. B. Khurgin, F. Monticone, "hw versus hk: dispersion and energy constraints on time-varying photonic materials and time crystals," Opt. Mater. Express 12, 3904 (2022); [arXiv:2208.00459](https://arxiv.org/pdf/2208.00459).
28. S. Saha et al., "Photonic time crystals: a materials perspective," Opt. Express 31, 8267 (2023).
29. M. Lobet et al., "New horizons in near-zero refractive index photonics and hyperbolic metamaterials," [arXiv:2306.01314](https://arxiv.org/pdf/2306.01314) (2023). † (first author)
30. D. M. Solis, R. Kastner, N. Engheta, "Time-varying materials in the presence of dispersion: plane-wave propagation in a Lorentzian medium with temporal discontinuity," Photon. Res. 9, 1842 (2021). †
31. S. A. R. Horsley, E. Galiffi, Y.-T. Wang, "Eigenpulses of dispersive time-varying media," Phys. Rev. Lett. 130, 203803 (2023). †

**Non-equilibrium electron kinetics in metals**

32. B. Rethfeld, A. Kaiser, M. Vicanek, G. Simon, "Ultrafast dynamics of nonequilibrium electrons in metals under femtosecond laser irradiation," Phys. Rev. B 65, 214303 (2002). †
33. B. Y. Mueller, B. Rethfeld, "Relaxation dynamics in laser-excited metals under nonequilibrium conditions," Phys. Rev. B 87, 035139 (2013). †
34. E. Carpene, "Ultrafast laser irradiation of metals: beyond the two-temperature model," Phys. Rev. B 74, 024301 (2006). †
35. G. Della Valle et al., "Real-time optical mapping of the dynamics of nonthermal electrons in thin gold films," Phys. Rev. B 86, 155139 (2012). †
36. Y. Dubi, Y. Sivan, "'Hot' electrons in metallic nanostructures: non-thermal carriers or heating?," Light Sci. Appl. 8, 89 (2019). †
37. A. M. Brown, R. Sundararaman, P. Narang, W. A. Goddard III, H. A. Atwater, "Ab initio phonon coupling and optical response of hot electrons in plasmonic metals," Phys. Rev. B 94, 075120 (2016).
38. M. Uehlein, S. T. Weber, B. Rethfeld, "Implementation of the electronic non-equilibrium in the two-temperature model," [arXiv:2112.01851](https://arxiv.org/abs/2112.01851); doi 10.3390/nano12101655.
39. M. Uehlein et al., "Capturing non-equilibrium electron dynamics in metals accurately and efficiently," [arXiv:2503.09479](https://arxiv.org/abs/2503.09479) (2025).
40. S. Roden et al., "Thermalization of optically excited Fermi systems: electron-electron collisions in solid metals," [arXiv:2601.11371](https://arxiv.org/abs/2601.11371) (2026).
41. I. Maliyov, J. Yin, J. Yao, C. Yang, M. Bernardi, "Dynamic mode decomposition of nonequilibrium electron-phonon dynamics: accelerating the first-principles real-time Boltzmann equation," npj Comput. Mater. 10, 123 (2024).

**Neural operators and learned kinetic closures**

42. Z. Li et al., "Fourier neural operator for parametric partial differential equations," ICLR 2021; arXiv:2010.08895. †
43. N. Kovachki et al., "Neural operator: learning maps between function spaces with applications to PDEs," J. Mach. Learn. Res. 24, 89 (2023). †
44. L. Lu, P. Jin, G. Pang, Z. Zhang, G. E. Karniadakis, "Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators," Nat. Mach. Intell. 3, 218 (2021). †
45. Z. Li et al., "Physics-informed neural operator for learning partial differential equations," arXiv:2111.03794. †
46. T. Xiao, M. Frank, "RelaxNet: a structure-preserving neural network to approximate the Boltzmann collision operator," [arXiv:2211.08149](https://arxiv.org/pdf/2211.08149) (2022).
47. J. Y. Lee, S. Schotthofer, M. Frank and co-authors, "Structure-preserving operator learning: modeling the collision operator of kinetic equations," [arXiv:2402.16613](https://arxiv.org/abs/2402.16613) (2024). † (author order)
48. J. Y. Lee, G. J. Jung, B. C. Lim, H. J. Hwang, "FourierSpecNet: neural collision operator approximation inspired by the Fourier spectral method for solving the Boltzmann equation," [arXiv:2504.20408](https://arxiv.org/abs/2504.20408v1) (2025).
49. B. Hu, K. Qi, "Solving the BGK model and Boltzmann equation by Fourier neural operator with conservative constraints," [arXiv:2510.13047](https://arxiv.org/pdf/2510.13047) (2025).
50. Y. Zhao, W. Burby, A. Christlieb, H. Lei, "Data-driven construction of a generalized kinetic collision operator from molecular dynamics," arXiv:2503.24208 (2025).
51. S. Burles, E. Camporeale, O. Pezzi, "A neural operator closure for Landau damping in electrostatic plasma," [arXiv:2607.29364](https://arxiv.org/abs/2607.29364) (2026).
52. M. Luo, A. R. Bell, F. Miniati, S. M. Vinko, G. Gregori, "Resolution-independent machine learning heat flux closure for ICF plasmas," [PRX Intelligence](https://journals.aps.org/prxintelligence/abstract/10.1103/9l4n-mnz6) (2026); [arXiv:2604.03439](https://arxiv.org/abs/2604.03439).
53. K. Shah, A. Cangi, "Accelerating electron dynamics simulations through machine learned time propagators," ICML 2024 AI4Science Workshop; [arXiv:2407.09628](https://arxiv.org/pdf/2407.09628).
54. K. Shah, A. Cangi, "Machine learning time propagators for time-dependent density functional theory simulations," [arXiv:2508.16554](https://arxiv.org/html/2508.16554v3) (2025).

**Neural operators in photonics; machine learning for ultrafast metrology**

55. J. Gu et al., "NeurOLight: a physics-agnostic neural operator enabling parametric photonic device simulation," [arXiv:2209.10098](https://arxiv.org/pdf/2209.10098) (2022).
56. P. Ma, H. Yang, Z. Gao, D. S. Boning, J. Gu, "PIC2O-Sim: a physics-inspired causality-aware dynamic convolutional neural operator for ultra-fast photonic device FDTD simulation," [arXiv:2406.17810](https://arxiv.org/abs/2406.17810) (2024).
57. "Physics-informed Fourier neural operator for nanophotonic scattering in nonlinear media," poster, [APS Global Physics Summit 2026](https://meetings-archive.aps.org/smt/2026/mar-h00-01/159/). A related Belmont University talk lists B. L. Blackwood and S. Hawley.
58. T. Zahavy, A. Dikopoltsev, O. Cohen, S. Mannor, M. Segev, "Deep learning reconstruction of ultrashort pulses," Optica 5, 666 (2018); [arXiv:1803.06024](https://arxiv.org/pdf/1803.06024).
59. H. Zhang, Y. Xu, et al., "Hybrid deep reconstruction for vignetting-free upconversion imaging through scattering in ENZ materials," [arXiv:2508.13096](https://arxiv.org/abs/2508.13096) (2025).
60. A. Ball et al., "A space-time knife-edge in epsilon-near-zero films for ultrafast pulse characterization," [arXiv:2408.00577](https://arxiv.org/abs/2408.00577) (2024).

**Inference**

61. K. Cranmer, J. Brehmer, G. Louppe, "The frontier of simulation-based inference," Proc. Natl. Acad. Sci. 117, 30055 (2020). †
