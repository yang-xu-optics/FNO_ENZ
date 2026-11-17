# Working log

Chronological record of what was done, what was found and what was decided. Newest session last. Reference numbers `[n]` are the project plan's. Add an entry every session; record negative results and dead ends too.

**Conventions.** Each session has: scope, what was done, results (with the script that reproduces them), deviations from the plan, open items. Decisions that change the plan also go in the decision table at the end.

---

## Session 1 — 2026-11-17 — repository set-up, reference library, literature ingestion, WP1 first pass

**Scope requested.** Create a folder holding all reference articles and download what can be downloaded; ingest all of them into one summary file that serves as the theory cornerstone; carry out the first step of the plan (WP1, ground-truth kinetic solver); keep this log; organise the project so later work packages slot in.

### 1.1 Environment and repository

- Repository `yang-xu-optics/FNO_ENZ` was empty. Work is on branch `claude/eager-euler-1aan36`.
- The session's default cloud environment blocked every literature host. A dedicated cloud environment **`FNO_ENZ`** was created (personal scope, saved for this repository) with package registries plus an allowlist of 19 hosts: arXiv (+ export API), OpenAlex, Crossref, Semantic Scholar, Nature (+ its login redirect), Optica, AIP, De Gruyter, the Sivan group site, Europe PMC (+ EBI API), APS meetings archive, OSTI, eScholarship, the Segev group site, HAL, Belmont repository. Future sessions for this repository should start in that environment.
- Machine used: 4 CPU cores, 15 GB RAM, no GPU. Python 3.13, numpy 2.5, scipy 1.18, numba 0.68, h5py 3.16, matplotlib 3.11.
- A setup script was saved on the `FNO_ENZ` environment: on a freshly built container it runs `pip install -e ".[dev]"` in the checkout. It has not yet been exercised by a fresh session.
- **The repository is public.** Consequence: the PDFs, TeX sources and extracted text are git-ignored and rebuilt by `references/fetch_references.py`; only the manifest, the index and the written summary are committed. Decision recorded as D1.
- Layout created: `docs/` (plan, this log, theory notes, gate reports), `references/`, `src/enzkin/`, `tests/`, `scripts/`, `results/wp1/`, `data/` (ignored). See `README.md`.

### 1.2 Reference library

- All 61 references were resolved against Crossref and the arXiv API (titles, authors, venue, DOI, arXiv id and version). Every entry the plan marked † ("cited from memory") is confirmed. 26 citations need updating (published since, author lists, a retitled paper); the list is in `references/README.md`. Notable: [16] has three authors, not one; [22] was retitled and published in Optica; [15] is Phys. Rev. Applied 19, 044043.
- **51 of 61 PDFs obtained**: 47 from arXiv (with TeX source for most), [8] via OSTI, [17] from Nature Communications, [21] via eScholarship, [28] from the Segev group site, and [16] from arXiv (also on the authors' site).
- **Not obtained (10):**
  - 7 paywalled with no legitimate open copy found: [1] Alam 2016, [3] Reshef 2019, [4] Kinsey 2019, [7] Wang 2019, [32] Rethfeld 2002, [33] Mueller 2013, [34] Carpene 2006. Abstracts retrieved from OpenAlex for all but [3], [4].
  - 2 open-access but behind sites that reject scripted clients: [20] Lustig 2023 (full text read from Europe PMC XML instead) and [35] Della Valle 2012 (abstract only).
  - [57] is a poster with no paper; abstracts of the two related Belmont University records were retrieved.
- Things that did not work: a helper agent's web-fetch tool could not reach any publisher host (it is outside the session proxy); repeated retries against challenge-protected publisher pages were blocked by the session's permission system and were not pursued.
- To complete the folder: download the nine missing PDFs in a browser with library access and save them as `references/pdfs/NN_key.pdf` (names in `references/README.md`); `python references/fetch_references.py --text` then extracts text. The digests of the seven paywalled papers should be redone from the real papers.

### 1.3 Literature ingestion

- Method: ten parallel reading passes, one per topical group, each given the same template (source read, role in project, model equations with the paper's numbering, numbers with units and locations, results, limitations, use for the project) and paper-specific questions. Equations were taken from TeX sources where available. Readers were told to state only what the paper states and to flag printed inconsistencies.
- Output: [`references/LITERATURE_SUMMARY.md`](../references/LITERATURE_SUMMARY.md), about 65 000 words. Part I is the synthesis (physics picture, equations adopted with corrections, parameter sets, time scales, mechanisms, electromagnetic side, experiments as datasets, machine learning, inference, **contradictions with the plan**, what must come from authors). Part II is one digest per paper. An appendix lists what was read for each.
- Verification: quantities marked "[verified in code]" in Part I were re-derived and checked numerically while building the solver. The reading pass on [14–16] independently reached the same three corrections to [14] that the solver work found (§1.4).
- Caveats: supplements were mostly unavailable ([22], [24], [17], [21], [10], [52]); arXiv versions were read, not journal versions; seven digests are reconstructions (see §1.2); figure read-offs are approximate and marked.

**Findings that change the plan** (details in the summary, Part I §10):

1. Two-photon absorption across a 2.67 eV gap needs λ < 929 nm, so it is closed over the plan's entire 1.1–1.6 µm dataset range and for [22]'s 1300 nm pump. M2-versus-M3 can only be tested on 800 nm-pumped data ([20], [23]) — the few-cycle regime the plan defers to Tier 2. [22] and [23] also assume different gaps (3.9 and 2.67 eV).
2. A fixed-density thermal model already reverses the sign of the *transmittance* [24]. H2 has to be stated for Re ε or ω_p. M0 needs a temperature-dependent damping law with free parameters.
3. The "density change vs mass change" reading of the two current equations is the project's extrapolation; [30, 31] give the two mathematical forms but do not make that identification. Which form heated ITO obeys depends on whether the heating collisions conserve momentum and must be derived in WP3.
4. [47] already enforces one conservation law exactly; [46] already has an exact equilibrium fixed point; [49] learns a flow map, not a collision operator. The novelty claim stands (no paper treats Pauli blocking or two weighted invariants) but must be reworded.
5. "Train at 256, evaluate at 512 with less than twofold error growth" is unsupported by the cited FNO papers.
6. [58–60] are not precedents for double-blind FROG with a complex gate; the closest is [22]'s retrieval.
7. Several anomaly-table entries need qualifiers ([19] pump is 159 fs in intensity; [21]'s ">200 %" and "40 %"; [24] is theory; [15] models a different experiment from [19]).

### 1.4 WP1: the solver (`src/enzkin`)

Specification: [`docs/theory/wp1_kinetic_model.md`](theory/wp1_kinetic_model.md). Summary of design decisions:

| Decision | Reason |
| --- | --- |
| Units eV–fs–nm | all rates O(10⁻³–1), no SI exponents |
| Uniform, cell-centred energy grid, E_max = 4 eV | makes two-body energy conservation an index identity; avoids k = 0 |
| e–e level B1 = screened-Coulomb integral of [14] with κ recomputed from f | the plan's ground truth and surrogate target |
| e–e level B0 = random-k integral of [40] by FFT | confirmed: with a constant matrix element the integral is one convolution plus one correlation; equals the direct sum to 10⁻¹⁴ |
| e–ph as a pair master equation with a symmetric momentum cut-off | the cut-off printed in [14] breaks number conservation and triples G; the symmetric one reproduces [14] Eq. 27 exactly |
| U_abs and U_pair carried in the state | makes the energy audit a linear invariant, exact under Runge–Kutta |
| RK4 with rate-limited step | fixed 0.5 fs steps blow up at high excitation (B0 scatters at several per fs) |
| Interband terms as moment-level lifts of [22] and [23] | the papers give lumped models only; valence-band parameters are not printed |
| Thermal model built from the same operators | kinetic minus thermal isolates non-thermal electrons |

Corrections to printed formulas (all in the theory note): spurious 4π³ in [14] Eq. 15 (confirmed by a footnote in [16]); extra ħ in [14] Eq. 16–17; prefactor convention c_ee = 3/16 (printed) vs 1 (direct-term derivation), kept as a parameter; e–ph cut-off; [14] Eq. 26 garbled.

Reproduced numbers from [14, 15] with the printed parameters: E_F = 0.8805 eV (0.88), k_F = 3.541 nm⁻¹ (3.54), γ_e = 12.67 J m⁻³ K⁻² (12.7), q_TF = 3.83 nm⁻¹, K = 0.0651 eV⁻² fs⁻¹, G(T→0) = 2.51×10¹⁶ W m⁻³ K⁻¹ (2.5×10¹⁶), zero crossing of Re ε at 224–225 THz (≈225).

Performance of B1 after restructuring the lookup tables along diagonals (4.4× faster than the first version): 1.2 ms (N = 128), 10 ms (256), 78 ms (512) per evaluation on 4 cores. B0: 0.24–0.30 ms.

Tests: `pytest` runs 60 tests in about 40 s — conservation (10⁻¹¹ relative), Fermi–Dirac fixed points, H-theorem, FFT vs direct sum, Fermi-liquid limit, relaxation to the matched Fermi–Dirac function, detailed balance and coupling constant of the e–ph term, exact audits in driven runs with and without interband terms, stiff-regime stability, thermal limit for long pulses, transfer-matrix energy balance.

Dead ends and fixes worth remembering:

- First e–e kernel indexed N×N tables along diagonals with stride N+1 (cache-hostile): 350 ms at N = 512. Fixed by storing tables in (d, i₂) order.
- Fixed-step RK4 at 0.5 fs diverged (NaN) for B0 at 150 GW/cm² in-film: B0's rate grows with ρ³ and reached several per fs. Fixed by the rate-limited step.
- The Auger-like term re-solved for (μ, T_e) by nested bracketing at every call (60 ms): replaced by a warm-started 2-D Newton solve (0.08 ms) and a tabulated n₀(T).
- **The e–ph collision rate entering the optical damping was under-resolved.** The first Gate 1 grid check (N = 128 against 256) showed T_e and ω_p² converged to 0.15 % but ε″ and the damping off by tens of per cent. Cause: the rate was taken as the functional derivative of the discrete master equation, which has no quasi-elastic (sub-cell) phonons and samples the Bose factor at grid nodes; at N = 128 it gave 0.024 fs⁻¹ at E_F against a continuum 0.059. Fixed by adding the intra-cell term and hat-integrated Bose factors (theory note §5); the optical average is now identical from N = 64 to 1024. The population dynamics and the energy transfer were never affected. All validation and Gate 1 numbers in this log are from runs after the fix; the calibrated impurity rate dropped from 0.063 to the value quoted in §1.5.
- Two numba-parallel jobs on the 4-core machine thrash (a 30 s run took 450 s). Run benchmarks, tests and dataset generation one at a time.
- The first Gate 1 case list included an absorbed energy of 10⁹ J/m³, more than the 4 eV band can hold; the excitation term then has nothing to excite. The excitation now reports the power actually deposited so the audit stays exact, and the benchmark stays within the band's capacity.

### 1.5 WP1 validation against [15]

Script: `scripts/validate_un2023.py` → `results/wp1/validation_un2023.json`, figure `results/wp1/fig_validation_un2023.png`. Kinetic model, B1, N = 128.

Set-up. [15] does not print ε∞, the e–impurity rate, the gold permittivity or the incidence medium. Linear optics were calibrated to its plotted 300 K permittivity (0.172 + 0.399i at 230 THz), giving ε∞ = 3.451 and γ_imp = 0.0352 fs⁻¹ (the reading pass on [15] independently estimated 0.01–0.02 fs⁻¹ from its curves). Geometry: a scan of candidate stacks showed that only "60° in air, refracted through the glass substrate to 35.3°, then 40 nm ITO and 100 nm gold" has an absorptance that falls monotonically as Re ε rises from its 300 K value, as [15]'s absorptivity does. My first assumption (60° inside the glass) gave absorptance *rising* with heating and peak temperatures up to 3× too high; that run is superseded.

Result, self-consistent thin-film drive (transfer matrix + ε[f], screening recomputed from f):

| Pump | T_e,max this work / [15] (K) | Absorptivity this work / [15] | ε′_max this work / [15] | ε″_max this work / [15] |
| --- | --- | --- | --- | --- |
| 220 fs, 2.5 GW/cm² | 2510 / 2464 (+1.9 %) | 0.765 / 0.777 | 0.423 / — | 0.479 / — |
| 220 fs, 5 | 3424 / 3313 (+3.4 %) | 0.673 / 0.668 | 0.600 / 0.765 | 0.528 / 0.495 |
| 220 fs, 10 | 4564 / 4348 (+5.0 %) | 0.558 / 0.545 | 0.835 / — | 0.590 / — |
| 220 fs, 22 | 6344 / 5902 (+7.5 %) | 0.426 / 0.418 | 1.180 / 1.476 | 0.676 / 0.547 |
| 220 fs, 50 | 9767 / 8393 (+16 %) | 0.311 / 0.318 | 1.669 / 2.013 | 0.780 / 0.594 |
| 220 fs, 75 | 13539 / 10262 (+32 %) | 0.264 / 0.281 | 1.989 / 2.315 | 0.833 / 0.630 |
| 30 fs, 161 | 8412 / 7470 (+13 %) | 0.406 / — | 1.429 / 1.717 | 0.710 / 0.546 |
| 30 fs, 22 | 3270 / 3187 (+2.6 %) | 0.728 / — | 0.539 / 0.682 | 0.506 / 0.482 |

Reading of the comparison:

- Absorptivity agrees within 1–6 % at all six intensities, and peak electron temperature within 2–7 % up to 22 GW/cm² (T_e ≤ 6000 K). The self-consistent loop, the thermodynamics and the e–ph coupling are therefore right at the level the paper can be read.
- **Peak temperature is 16–32 % too high at 50–75 GW/cm² with the stated 4 eV band, and within 5 % with a 6 eV band.** Rerunning the three highest 220 fs cases with E_max = 6 eV (N = 192, same ΔE) gives T_e,max = 6203, 8777, 10771 K against [15]'s 5902, 8393, 10262 K (+5.1, +4.6, +5.0 %), absorptivity 0.429, 0.315, 0.269 against 0.418, 0.318, 0.281 (these three runs predate the damping-rate fix below, which changes T_e,max by less than 0.2 %). Cause: at these intensities the absorbed energy (up to 5.8×10⁸ J/m³) approaches what a 4 eV band can hold (about 2.7 eV per electron at uniform occupation, i.e. ≈5×10⁸ J/m³ above the cold state), so the electron heat capacity collapses and the matched temperature runs away. [14] states a band "total size" of 4 eV, but a comment in its TeX source says the band width was being adjusted, and [15] evidently used a wider one. A uniform +5 % remains at all intensities. **E_max is therefore a model parameter that matters above ≈30 GW/cm²**; the default stays at the printed 4 eV and the wider band is recommended for high-intensity runs (open item 3).
- **Δε′ is 20–28 % smaller than in [15] at every intensity.** Consistent with the paper's plotted permittivity implying a plasma frequency 8 % above what its printed n, m₀, C give (ħω_p = 1.87 vs 1.73 eV) and ε∞ = 4; not resolvable without the authors.
- **ε″_max matches [15] to 1–6 % if the Thomas–Fermi wavevector is frozen at its 300 K value** (0.501, 0.578, 0.627, 0.657, 0.576, 0.485 against 0.495, 0.547, 0.594, 0.630, 0.546, 0.482 for 5, 22, 50, 75 GW/cm² at 220 fs and the two 30 fs cases), and is 5–32 % too high when κ follows f. This suggests [15] did not recompute the screening from the hot distribution, although its Eq. 15 defines q_TF through f. To ask the authors.
- With the paper's own absorptivity imposed as a Gaussian energy deposit ("energy" mode) T_e,max is +3 % at 2.5 GW/cm² but +26 % at 50: at high intensity the true absorption is front-loaded, so this mode is only a clean check at low intensity.

Status against the plan's WP1 validation list:

| Plan item | Status |
| --- | --- |
| Reproduce T_e and permittivity transients of [15], 220 fs and 30 fs | partly: absorptivity within 1–6 % everywhere; T_e within 2–7 % up to 22 GW/cm² and for 30 fs/22 GW/cm², and within 5 % at 50–75 GW/cm² once the band is widened to 6 eV; Δε′ 20–28 % low and ε″ sensitive to the screening treatment; closing these needs the authors' unprinted parameters |
| Conserve number and total energy to 10⁻⁶ | done: 10⁻¹³ |
| Relax to Fermi–Dirac | done (test) |
| Recover the two-temperature limit [38] | done for long pulses (test: T_e,max within 3 %, ω_p² within 1 %) |
| Compare with the reduced athermal model of [39] | not done: [39] is written for transition/noble-metal DOS and gives no timings; a conserving RTA baseline is in place instead |
| Tier 2 density-matrix model [16] | not started (stretch) |

### 1.6 Gate 1

See [`docs/gates/G1_report.md`](gates/G1_report.md). Script: `scripts/gate1_benchmark.py` → `results/wp1/gate1.json`.

**Verdict: Gate 1 does not pass as written.** Its two criteria point in opposite directions.

| Criterion | Measured |
| --- | --- |
| B1 costs ≥ 100× B0 per evaluation | ×1 at N = 64, ×6 at 128, ×46 at 256, ×266 at 512 (B1: 0.3, 1.7, 12.8, 78 ms; B0: 0.23–0.29 ms; 4 CPU threads) |
| B0 misses B1's permittivity transient by > 2 % | yes in all 16 conditions: 7–221 % on Re ε, 3–47 % on ω_p² |
| (grid check) N = 128 against 256, B1 | every observable within 0.22 % of the transient amplitude |

Consequences:

- The fallback the plan attaches to a failed gate, "use FFT-evaluated B0 inside WP3", should not be taken: B0 is too inaccurate. Its error comes mainly from its collision *rate* (which sets the optical damping and grows with the cube of the density of states), less from the distribution it produces.
- The exact integral B1 is affordable for WP3 at the resolution the transients need: a film run of 10⁶ evaluations is about 28 minutes at N = 128 on this machine (3.5 h at 256, 22 h at 512), against the plan's assumed 1–10 days.
- The surrogate's justification moves from WP3 to WP4 (10³–10⁴ film runs, i.e. 470–4700 h here). A GPU port of B1 has to be benchmarked before that case can be judged, and the WP2 criterion "≥100× over B1 at 512 points" restated against it. Any surrogate must also deliver the collision rate, not only Q[f].
- Side result for H1: the thermal two-temperature model built from the same operators misses the B1 permittivity transient by 11–17 % for 20 fs pumps and 3–5 % for 220 fs pumps.

This is the outcome the plan lists as its largest risk; the plan's response to it (Maxwell–kinetic physics carries the flagship paper, the operator goes to the methods paper and to inference speed) applies. The decision is recorded as pending in the decision table (D9).

### 1.7 Pilot dataset

Script: `scripts/make_dataset.py` → `data/pilot.h5` (git-ignored; regenerate with `python scripts/make_dataset.py --n 8 --N 128 --out data/pilot.h5`).

- 8 trajectories, B1, N = 128, 2 ps at 2 fs sampling, each with f, Q_ee, Q_eph, Q_exc, the e–e rate, scalars, and 32 perturbed off-trajectory states with their collision terms. 35–43 s and 4.2 MB per trajectory; audit 10⁻¹¹–10⁻¹⁵; f stays in [0, 1].
- Sampled conditions ranged from 30 fs to 246 fs, 0.1 to 58 GW/cm² in-film, 1132–1547 nm, peak T_e 326–7486 K.
- **A third of the plan's sampling range does not fit in a single band.** Drawing from the plan's ranges (20–250 fs, 0.1 GW/cm²–10 TW/cm² in-film), 34 % of cases would nominally deposit more than half of what the 4 eV band can hold and 28 % more than all of it. The first attempt, without a screen, stalled on such a case (148 fs at 1 TW/cm²). The generator now rejects and resamples them and records the count. The plan's upper intensity range needs a model with higher bands, or an explicit validity limit of roughly 10¹¹ W/cm² × (100 fs / pulse duration) in-film.
- Scaling to the plan's 2×10⁴ trajectories: about 220 hours on this machine at N = 128 (46× more at 512) and ≈80 GB at this sampling. A cluster or the GPU port is needed; runs are independent and seeded per index (`--start`, `--n`).

### 1.8 Compute resource estimate (asked: would a GPU cluster help, and how much to request)

Measured on this machine (Intel Xeon 2.1 GHz), one B1 evaluation: 4.0 / 37 / 231 ms on one thread and 1.2 / 9.5 / 86 ms on four threads at N = 128 / 256 / 512. Runs are independent, so cluster cost is counted in single-thread core-hours.

| Workload | Evaluations | N = 128 | N = 256 | N = 512 |
| --- | --- | --- | --- | --- |
| One 2 ps trajectory (0D) | 1.8×10⁴ | 0.02 core-h | 0.19 | 1.2 |
| WP1 dataset, 2×10⁴ trajectories | 3.6×10⁸ | 400 core-h | 3 700 | 23 000 |
| One WP3 film run (100 cells, 4000 RK4 steps) | 1.6×10⁶ | 1.8 core-h | 16.5 | 103 |
| WP3 development, ~300 film runs | | 530 core-h | 5 000 | 31 000 |
| WP4, 10³–10⁴ film runs | | 1 800–18 000 core-h | 16 500–165 000 | 10⁵–10⁶ |

Not measured: any GPU number (no GPU here). Working estimate, to be replaced by a benchmark: a batched GPU implementation of B1 is worth roughly one 32–64-core node. Memory of a batched B1 on GPU: 100 cells × N³ × 8 bytes per intermediate tensor, i.e. 1.7 GB at N = 128, 13 GB at 256, 107 GB at 512 (chunk over cells above 128).

Conclusions: (i) the grid size sets everything, so the convergence sweep of open item 4 comes before any large request; (ii) WP1 and WP3 need CPU cores, not GPUs, and the code runs on a cluster unchanged; (iii) a GPU is needed for WP2 training and for gradient-based inference through a differentiable (JAX/PyTorch) port in WP4; (iv) sequential samplers in WP4 are limited by wall time per film run (about 2 minutes on a 64-core node at N = 128), which is where a surrogate or a GPU port earns its place. Suggested request: about 50 000 CPU core-hours, one GPU with at least 40 GB for 1 000–1 500 GPU-hours, 500 GB storage, preceded by a small pilot allocation for the GPU benchmark and the convergence sweep.

### 1.9 Open items, in priority order

1. **Write to the authors of [14–16]** (the plan already recommends it): ε∞, e–impurity rate, gold permittivity, grid and time step of [15]; whether κ was recomputed from f; the origin of the 8 % plasma-frequency gap; the 3/16 prefactor; e–ph treatment above 6000 K. This closes the remaining validation gaps.
2. **Obtain the nine missing PDFs** and redo the seven reconstructed digests, [1] first (it is the WP4 calibration dataset and its two-temperature parameters were never read).
3. **Fix the band model for high intensity** (§1.5) before generating the full dataset: choose E_max (6 eV reproduces [15]; the plan's upper range of 10 TW/cm² in-film saturates any single-band model and needs either higher bands or a stated validity limit), rerun the validation at N ≥ 192, and repeat Gate 1 with it.
4. **Decide what Gate 1 means for WP2** (decision D9): port B1 to a GPU (a batched contraction over three indices) and benchmark it; rerun the grid-convergence check over the whole case matrix at N = 128, 256, 512, including weak excitation where the cold Fermi edge is sharpest; then restate the WP2 speed criterion against that baseline or shrink WP2 as the plan foresees.
5. **Derive the l = 1 (current) moment of the kinetic equation** to decide between the two current equations of WP3; this is also a candidate result in its own right.
6. **Revise the plan** for the seven findings of §1.3 (TPA wavelength constraint, H2 wording, M0 definition, novelty wording, 256→512 criterion, WP5 metrology precedents, anomaly-table qualifiers).
7. Microscopic interband kernels (impact ionisation with a valence band; two-photon matrix element) once [22]'s supplement and a valence-band model are in hand.
8. Add stimulated emission study: the `net` excitation variant exists but has not been compared with `absorption_only` at high intensity.
9. Environment setup script (install `requirements.txt`, optionally run `fetch_references.py`) can be saved to the `FNO_ENZ` environment once the dependency list settles in WP2.

---

## Decision table

| Id | Date | Decision | Why |
| --- | --- | --- | --- |
| D1 | 2026-11-17 | Reference PDFs, TeX sources and extracted text are not committed; manifest + fetch script + written summary are | public repository; most papers are not licensed for redistribution |
| D2 | 2026-11-17 | Default e–e prefactor c_ee = 3/16 (as printed in [14]); treated as a fit parameter | comparability with [14, 15]; literature prefactors differ by ×3 |
| D3 | 2026-11-17 | e–ph term uses the symmetric triangle cut-off, not the printed one | exact number conservation; reproduces the published coupling |
| D4 | 2026-11-17 | Screening recomputed from f by default; `dynamic_screening_from_f=False` available | plan says "screening recomputed from f"; frozen screening matches [15]'s ε″ better (open question) |
| D5 | 2026-11-17 | Linear optics of the reference material calibrated to [15]'s plotted ε(230 THz, 300 K): ε∞ = 3.451, γ_imp = 0.0352 fs⁻¹ (`enzkin.calibrate`) | not printed in the paper |
| D6 | 2026-11-17 | [15] validation geometry: 60° in air through the glass substrate | only geometry reproducing the published absorptivity trend |
| D7 | 2026-11-17 | Interband models M2, M3 implemented as moment-level lifts with a declared free valence-band mass | published models are lumped; parameters not printed |
| D8 | 2026-11-17 | Default band width stays at the printed E_max = 4 eV; E_max is exposed as a material parameter and 6 eV is recommended above ≈30 GW/cm² | 4 eV saturates and overshoots T_e by 16–32 % at 50–75 GW/cm²; 6 eV reproduces [15] within 5 % |
| D9 | 2026-11-17 | **Pending (project lead):** Gate 1 did not pass as written. Proposed: run WP3 on exact B1 at N = 128–256; do not use B0; decide on the surrogate after a GPU B1 benchmark | B1 is only 6–46× costlier than B0 at the resolution needed, B0 misses the transient by 7–221 % |
| D10 | 2026-11-17 | Dataset generator rejects cases whose nominal deposit exceeds half the band's capacity | a third of the plan's sampling range saturates a single 4 eV band |
