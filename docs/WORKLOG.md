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

(Wording note: an earlier version of the Gate 1 report said the surrogate is "not needed for the flagship physics". That was too broad. It is not needed to run the WP3 film simulations or the H3 test; the mechanism ranking of WP4 is also a flagship claim and is where the cost question remains.)

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

### 1.9 Progress report and repository visibility

- Progress report: `python scripts/make_progress_report.py` → `docs/reports/progress_2026-11-17.html`. One self-contained page (inline CSS, inline SVG charts with hover read-outs, light and dark themes); every chart and result table is read from `results/wp1/validation_un2023.json`, `results/wp1/gate1.json` and `results/wp1/gate1_series.npz`. Also published as a private page at https://claude.ai/artifact/TPWgtV8XsG8qHVSaf3eQ7J (republishing the same file keeps the link).
- Request to make the repository private: **not done, cannot be done from a working session.** The session's GitHub access has no operation that changes repository settings. The owner changes it on GitHub: repository → Settings → General → Danger Zone → Change repository visibility → Make private. The repository was still public when this entry was written.
- Once it is private: decide whether to track the paper PDFs (about 170 MB). If yes, remove the `references/pdfs/` line from `.gitignore` (and optionally `references/sources/`, `references/text/`), update the "never commit PDFs" rule in `CLAUDE.md` and the note in `references/README.md`, and record it in the decision table.

### 1.10 Open items, in priority order

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

## Session 2 — 2026-11-17 — plan revision 2, JAX port of the exact integral, grid-convergence sweep

**Scope requested.** The project lead decided to drop the machine-learning surrogate for now (decision D9), and asked for the plan to be revised and for the first step of the new sequence to start: confirm that the exact solver is fast enough (grid-convergence sweep including weak excitation and the 6 eV band; GPU port of the exact integral).

### 2.1 Plan revision 2

- `docs/PROJECT_PLAN.md` is now revision 2, "Maxwell-Kinetic Model of ENZ Electron Dynamics". Revision 1 is unchanged at `docs/archive/PROJECT_PLAN_v1_2026-11-16.md`. Reference numbers are the same.
- Changes, each with its reason, are tabulated at the top of the plan. In short: the neural operator (old WP2) moves to a deferred section with the conditions under which it returns; WP3 runs on the exact integral and starts now; the FFT integral is a comparison model only; H2 is stated for Re ε; the two-photon model is tested on 800 nm data only; M0 is a two-temperature model with a damping law; the current equation becomes a derivation task; the intensity range is limited to what one band can hold; the observation and dataset tables carry the qualifiers found in the literature pass; the methods paper becomes a solver-and-benchmark paper.
- Gates: G1 closed (not passed as written); new G1b (confirmation of grid convergence and of film-run cost on a cluster node or GPU); G2 suspended with WP2; G3 and G4 unchanged.
- The reference list was regenerated from `references/manifest.json`, so all 61 citations carry the verified authors, venue, arXiv id and DOI.
- `README.md` and `CLAUDE.md` updated to match. The repository name is unchanged.

### 2.2 JAX port of the exact integral (`src/enzkin/collisions/ee_jax.py`)

- Same discrete integral as the numba kernel, written as array operations: for a fixed index shift d the kernel factor is an N×N matrix built from two vectors, and the sums are two matrix–vector products; `jax.lax.map` runs over the 2N−1 shifts with `d_batch` of them in memory at a time. Batched over cells with `vmap`, differentiable, 64-bit.
- Tests (`tests/test_ee_jax.py`, skipped when JAX is absent): agreement with numba to 10⁻¹¹ for the collision term and the rate, single-state and batched calls agree, number and energy conserved, Fermi–Dirac fixed point, and the gradient of the entropy production matches a finite difference.
- Reference timing on **one CPU core** (`scripts/bench_b1_jax.py --numba`, `results/wp1/bench_b1_jax_cpu.json`): 8.0 ms per cell at N = 128 and 48 ms at N = 256, against 3.7 and 27 ms for numba on the same core. JAX on a CPU is about twice as slow as numba; the port exists for the GPU and for gradients.
- **Not done: the GPU benchmark.** This environment has no GPU. On a GPU node: `pip install -U "jax[cuda12]"`, then `python scripts/bench_b1_jax.py --N 128 192 256 --cells 1 16 100 400`; the script prints the time per call, per cell and per film run, and writes `results/wp1/bench_b1_jax_gpu.json`. This closes the GPU half of gate G1b.

### 2.3 Grid-convergence sweep and gate G1b (`scripts/convergence_sweep.py`, `docs/gates/G1b_report.md`)

**Design.** Exact integral, reference material calibrated once (on 256 points), 1300 nm, 24 conditions: absorbed energy 10⁶, 10⁷, 3×10⁷, 10⁸, 3×10⁸ J/m³ for pulses of 20, 50, 100, 220 fs, and four fixed in-film intensities (100 and 500 GW/cm² at 20 fs, 10 and 40 GW/cm² at 220 fs). Peak electron temperatures run from 431 K to 10 100 K. Stages: A, all conditions at 128 and 256 points over 4 eV; B, five strong conditions over 6 eV at 192 and 384 points; C, four conditions at 512 points; D, the five strong conditions over 8 eV at 256 points; E, the two hottest over 10 eV at 320 points. 72 runs, about 1 h 40 min on this machine. Results: `results/wp1/convergence_sweep.json` (+ `_series.npz`); summary rebuilt with `--summarise`; figures by `scripts/make_figures.py`; report by `scripts/make_g1b_report.py`.

**Measure.** Largest deviation of the transient y(t) − y(t₀) from the finer grid's, divided by the finer grid's transient amplitude.

**Dead end worth remembering.** The first summary compared y(t) itself and showed 13–18 % "miss" at weak excitation. That was not a dynamic error: the cold-state permittivity differs between grids by a constant (−1.0×10⁻³ in Re ε and +8.9×10⁻⁵ in Im ε at 128 points against 512), which is up to 24 % of the weakest transient. It is the same for every condition and disappears when ε∞ and the impurity rate are calibrated on the grid in use, as the validation and gate scripts do. The sweep uses one material on all grids, so the measure was changed to the transient and the static offset is reported separately. **Rule: calibrate the linear optics on the grid in use.**

**Energy step.**

- 128 against 256 points, all 24 conditions: at most 0.72 % (220 fs, 10⁶ J/m³). 128 against 512 (four conditions): at most 0.96 %. 256 against 512: at most 0.24 %.
- Observed order of convergence 1.3–1.7. Extrapolated to zero step, the 128-point error is 1.1 % and the 256-point error 0.4 % for the weakest long pulse.
- The error falls with excitation (0.4 % or less against the next grid above 1800 K) and is the same fraction in T_e, ω_p², Re ε and Im ε: a scale error from the quadrature of the cold Fermi edge (k_B T = 26 meV against a 31 meV cell), not a shape error.
- 6 eV band: 192 against 384 points, at most 0.21 %.

**Band width — the larger error, and not anticipated by the gate.** At equal energy step:

| Narrow → wide | within 1 % up to | beyond |
| --- | --- | --- |
| 4 → 6 eV | T_e,max ≈ 3500 K | 2.7 % at 5900 K, 12 % at 8800 K (T_e,max itself 16 % too high) |
| 6 → 8 eV | ≈ 6700 K | 3.1 % near 8400 K |
| 8 → 10 eV | not within 1 % at the hottest | 1.3 % near 8400 K |

- All three pairs fall on one curve against the largest occupation reached in the narrower band's top cell (`f_top_max`, recorded by every run): below 1 % while it stays under 4×10⁻⁵, above 1 % from 2×10⁻⁴. **Rule: widen the band when the top-cell occupation passes 4×10⁻⁵**, keeping the 31 meV step (128, 192, 256 points for 4, 6, 8 eV).
- A narrow band always errs upward (electrons pile up under the cut-off), so earlier 4 eV numbers above about 4000 K overstate T_e and the transient. This affects the Gate 1 conditions at 3×10⁸ J/m³ and is the same effect as the T_e overshoot against [15] (§1.5).
- For the hottest condition the dependence converges slowly: 12.6 %, 3.1 %, 1.2 % for 4→6, 6→8, 8→10 eV. Several eV above the Fermi level the Kane dispersion is an extrapolation, so this is also a statement about the model: results above roughly 7000 K depend at the percent level on states it does not describe.

**Cost** (`results/wp1/bench_b1_jax_cpu.json`, rerun with nothing else on the machine, one core): numba 3.7, 9.9, 27.0, 67.3 ms per evaluation at 128, 192, 256, 320 points; a film run of 1.6×10⁶ evaluations is 1.7, 4.4, 12.0, 29.9 core-hours. (The Gate 1 report counted 10⁶ evaluations per film run; the count here is 100 cells × 4000 RK4 steps × 4 stages.)

**Gate G1b.** Grid part met at 128 points per 4 eV (at the line for the weakest excitation, comfortably elsewhere; 256 points meet it with a margin of 4). 6 eV grid part met. CPU cost part met. **GPU part outstanding**: no GPU here (§2.2). Consequences for WP3: exact integral, 31 meV step, band width set by the top-cell rule.

**Not covered.** Thin-film (self-consistent field) runs, other wavelengths, other material parameters, time-step error.

### 2.4 Housekeeping

- `pytest` on the bare command line picks an interpreter without numpy in this container; `python3 -m pytest` works (65 tests pass, 53 s). `CLAUDE.md` updated.
- The CPU benchmark file was overwritten by the rerun; §2.2 quotes the new numbers.
- Commits are pushed to `claude/eager-euler-1aan36`; partial sweep results were committed once mid-run.

### 2.5 Open items, in priority order (supersedes §1.10)

1. **Run the GPU benchmark** on a GPU node (§2.2) and rerun `python scripts/make_g1b_report.py`; this closes G1b.
2. **Start WP3** on the exact integral: derive the l = 1 (current) moment of the kinetic equation to choose the current equation; one-dimensional Maxwell solver coupled to one kinetic cell per grid point; validation ladder of the plan.
3. **Rerun what used the 4 eV band above about 4000 K** with the top-cell rule: the [15] validation at 50–75 GW/cm² (6 or 8 eV) and the Gate 1 conditions at 3×10⁸ J/m³. Add an automatic warning in `enzkin.integrate.run` when `f_top_max` passes the threshold.
4. **Write to the authors of [14–16]**: ε∞, impurity rate, gold permittivity, grid, time step and band cut-off of [15]; whether κ follows f; the 8 % plasma-frequency gap; the 3/16 prefactor.
5. **Obtain the nine missing PDFs** and redo the seven reconstructed digests, [1] first.
6. Extend the convergence sweep to the thin-film drive and to a second wavelength once WP3 has a film solver.
7. Benchmark dataset (plan WP1): regenerate with the band rule in place of the fixed 4 eV screen.
8. Repository visibility is still public; the owner changes it on GitHub (§1.9).

---

## Session 3 — 2026-11-18 — scope questions: the adiabatic approximation and few-cycle pumping

**Scope requested.** Two questions from the project lead, no change to the plan: (1) does the project drop the adiabatic approximation of [15]; (2) can it predict ENZ dynamics driven by few-cycle pulses, and what exists on few-cycle-driven non-equilibrium electron dynamics in general.

### 3.1 Answers given

- (1) Yes for the optics, in WP3: no memory-free ε(t), no steady-state transfer-matrix field, no film average. WP1 keeps all three on purpose to match [15]. WP3 keeps cycle-averaged excitation and instantaneous rates. The adiabatic model stays as the H3 baseline.
- (2) Not as planned. The plan already limits flagship claims to pumps of 20 fs and longer (Background, risk table) and puts few-cycle pumps in the Tier-2 stretch goal. Removing the adiabatic optics is necessary but not sufficient: few-cycle pumping also breaks cycle-averaged absorption, the isotropic distribution and, marginally, Markovian collisions (§3.2).

### 3.2 Field regime of the pumping conditions (`scripts/field_regime.py` → `results/wp1/field_regime.json`)

Peak displacement of the Fermi sphere by the local field, dk = eE_loc/ħω, against k_F = 3.54 nm⁻¹, and pulse length in optical cycles. Local fields for [15] from the thin-film stack of the validation (cold film, so upper values); for [20], [23] estimated as vacuum field × normal-incidence Fresnel transmission.

| Condition | cycles (FWHM) | E_loc (V/nm) | dk/k_F | mass factor at leading edge / at E_F |
| --- | --- | --- | --- | --- |
| [15], 220 fs, 22 GW/cm² incident | 51 | 0.97 | 0.29 | 1.20 |
| [15], 220 fs, 75 GW/cm² incident | 51 | 1.80 | 0.53 | 1.38 |
| [15], 220 fs, 200 GW/cm² (their stated adiabatic limit, 2.7 V/nm in [15]) | 51 | 2.93 | 0.87 | 1.63 |
| G1b sweep, 20 fs, 500 GW/cm² in-film | 4.6 | 1.94 | 0.57 | 1.41 |
| [20], 6 fs, 800 nm, ~2 TW/cm² | 2.2 | 3.12 | 0.57 | 1.41 |
| [23], 6.4 fs, 800 nm, 11 TW/cm² | 2.4 | 7.32 | 1.33 | 1.99 |

- The local field computed for the [15] stack at 200 GW/cm² (2.9 V/nm) matches the 2.7 V/nm that [15] quotes, a check on the reconstructed geometry.
- **Finding (my reading, not tested):** the field-driven displacement is already about half of k_F at the top of the range validated against [15] with 220 fs pulses. In the Kane band this modulates the electrons' mass within each cycle (the "fast" nonlinearity of [5–7, 11, 12]), which an isotropic, cycle-averaged f(E) cannot carry. For long pulses the accumulated heating term is expected to dominate ([5]: χ_slow/χ_fast up to ~40), but the WP1 model cannot check this. It is a reason to keep the l = 1 (current-carrying) part of the distribution in the WP3 derivation task rather than only the current equation.
- For 800 nm few-cycle pumping ħω = 1.55 eV exceeds E_F = 0.88 eV, and at 1300 nm the two are comparable, so a semiclassical field-driven Boltzmann equation (force term) cannot represent the quantised absorption; the density-matrix level of [16] is the minimum for those data. Time scales near the Markov limit: collision duration ħ/E_F ≈ 0.75 fs, unscreened plasma period ≈ 2.3 fs.

### 3.3 Literature

- Read: the companion Letter of [16] (Un, Sarkar, Sivan, arXiv:2608.11756, Aug 2026; main text, not the SI), digested at the end of [16] in `references/LITERATURE_SUMMARY.md`; Part I §4 (time scales), §6 (H3) and §11 updated. Key points: a thermal model driven by the oscillating net absorption Ṗ·E reproduces their full density-matrix permittivity for a 5 fs pulse; e–e thermalisation ≈12 fs at 1 TW/cm², rate ∝ |E₀|²; their prediction disagrees with the 40 fs transmission recovery of [20], and they name film gradients and transport as a possible cause (H3-type). Also checked: the published version of [16] (Phys. Rev. B 113, 165411) has the same content as the arXiv version that was digested.
- Found, not read, not in the library (candidates if the few-cycle regime is taken up): Scalora et al., Phys. Rev. A 112, 013502 (2025), time-domain model reproducing [20]; Herbst et al., Adv. Sci. 13, e16818 (2026), fieldoscopy of ITO nanocrystals with 10.7 fs two-cycle pulses; Yang et al., Nat. Phys. 15, 1022 (2019), high harmonics from an ENZ CdO film. General few-cycle solid-state dynamics: Borsch et al., Nat. Rev. Mater. 8, 668 (2023) (review); Sommer et al., Nature 534, 86 (2016); Volkov et al., Nat. Phys. 15, 1145 (2019); Hafez et al., Nature 561, 507 (2018); Maxwell-TDDFT (Yamada and Yabana, Phys. Rev. B 99, 245103 (2019)); Maxwell + semiconductor Bloch equations (Kilen et al., Phys. Rev. Lett. 125, 083901 (2020); Floss et al., Phys. Rev. A 97, 011401 (2018)); quantum kinetics (Fürst et al., Phys. Rev. Lett. 78, 3733 (1997); Huber et al., Nature 414, 286 (2001)).
- Caveat for the plan's Tier 2 (my reading): the density matrix of [16] absorbs through finite-size (Landau-type) damping in a 4 nm box and explicitly neglects phonon- and electron-assisted absorption, which is what absorbs in a 40–300 nm film. Implementing [16] in a film cell needs a coherent treatment of collision-assisted absorption first; it is a formulation problem, not only a cost problem.

### 3.4 Open items added

- **For the project lead:** whether the few-cycle regime ([20, 23, 25], Herbst et al.) stays a stretch goal or becomes a work package. Plan unchanged until decided.
- WP3 derivation task: carry the l = 1 part of the distribution and quantify the intra-cycle (Kane-mass) nonlinearity at the [15] conditions (dk/k_F ≈ 0.5).
- Cheap intermediate few-cycle model suggested by the companion Letter: WP1 populations (or M0) driven by the instantaneous J·E from the Maxwell solver instead of the cycle-averaged absorption term; validate against the Letter's Fig. 4 if its parameters can be matched.

## Session 4 — 2026-11-18 — ingest of Un, Sarkar & Sivan (Nanophotonics 2026) and plan revision 2.1

**Scope requested.** The project lead supplied the published version and Supporting Information of Un, Sarkar & Sivan, "Origin, strength, and speed of the nonlinear optical response of transparent conducting oxides to single-cycle optical pulses", Nanophotonics 15, e70272 (2026), and asked to understand and ingest it, then revise and re-evaluate the project if necessary.

### 4.1 Ingest

- Read: all 11 pages and the 9-page SI (figures by eye). Added as reference **[62]**: manifest entry (`un2026nanophot`, arXiv:2608.11756, CC BY), full digest in `references/LITERATURE_SUMMARY.md` (after [16], whose short note on the arXiv version now points to [62]), Part I §1, §2, §3, §4, §6, §7, §10, §11 and the appendix updated, `references/README.md` updated. PDFs are in the local library (`references/pdfs/62_un2026nanophot*.pdf`, git-ignored); `fetch_references.py` now stores supplementary PDFs as `NN_key_si.pdf`.
- What the paper does: density matrix (populations + coherences) of a 4 nm ITO sphere, quasi-static, 5 fs pump at 175 THz, E₀ = 0.6–2.4 V/nm, with the e–e, e–ph and (new) Brooks–Herring e–impurity integrals. Findings: O(1), cumulative permittivity change (not Kerr-like) with a 2ω ripple from coherences; stimulated emission comparable to absorption at the highest field; a thermal model driven by the instantaneous net absorbed power Ṗ·E reproduces the density-matrix permittivity; e–e thermalisation ≈12 fs at "1 TW/cm²" (T̃_e ≈ 4500 K), ∝ E₀²; generalised Fermi-liquid formula; e–ph relaxation slow and distribution-independent; prediction at odds with [20] (larger change, slower decay), with film gradients/transport among the suggested causes.
- Printed-formula problems found (all in `docs/theory/un2026_notes.md`): (S9) sign in the bracket; (S22) missing 1/4; (S23) tanh argument; (S26) missing 1/(k_BT)²; (S25) missing 1/(2π)³; factor 2 between (S2) and (S7); ε∞ vs ε_b in the screening; ρ units. Unprinted: ε∞, ε_b, band top, grid, T̃_e(t).

### 4.2 Checks against enzkin (`scripts/check_un2026.py` → `results/wp1/check_un2026.json`)

- e–e: (S19) is the B1 integral with c_ee = 3/16. The generalised Fermi-liquid coefficients, implemented independently, match Fig. S1 to 2–7 % (A₀(E_F) = 9.46/8.51/6.06/4.65/3.73 vs ≈9.4/8.6/6.3/4.9/4.0 at 300–7000 K) once (S22)'s 1/4 is restored, and the formula reproduces the B1 thermal rate within 3 % near μ. Absolute rate: B1 gives 20 fs at 4500 K, 12 fs at ≈6500 K; the paper's "12 fs at 4500 K" is 1.5–2× faster; not resolvable without T̃_e(t).
- e–ph: Fig. 5c matches enzkin with **D = 8.6 eV** (within 1.2–1.5×), not with the 17.2 eV of [14] (6× too high). The authors changed the deformation potential and added Brooks–Herring impurities (n_imp = 0.05 nm⁻³, Z = 1).
- Consequences (bulk film, ε∞ = 3.45): Γ = 0.074 vs 0.175 fs⁻¹; ε(230 THz) = 0.14 + 0.17i vs 0.172 + 0.399i ([15]); G = 6.6×10¹⁵ vs 2.6×10¹⁶ W m⁻³ K⁻¹; Re ε recovery (1/e) of a 100 fs, 10⁸ J/m³ transient **1.9 ps vs 0.52 ps**, against 60–400 fs measured in ITO films. The e–ph coupling is not settled even within the group that defines the model.

### 4.3 Code

- `src/enzkin/collisions/eimp.py`: Brooks–Herring out-scattering rate (corrected S9), screened with ε_b and the κ of the current f; added to the damping in `KineticModel.damping` when `mat.n_imp > 0`. Off by default: no existing result changes.
- `Material` gains `n_imp` (nm⁻³, default 0) and `Z_imp` (default 1); preset `ito_un2026()` (D = 8.6 eV, n_imp = 0.05, Z = 1, gamma_imp = 0).
- Tests: closed form against the elastic-sphere angular integral; damping off by default, on and linear in n_imp. Theory: `docs/theory/wp1_kinetic_model.md` §5, §8.

### 4.4 Re-evaluation and plan revision 2.1 (`docs/PROJECT_PLAN.md`)

- **Unchanged:** the central question, H1–H3, the flagship scope (pumps ≥20 fs), WP3 as the next build. The paper supports the project's direction: the authors themselves cannot explain [20] with intraband kinetics in a uniform particle and point to film effects, which is what WP3/H3 test.
- **Tier 2 split** (D16): Tier 2a = field-resolved J·E drive of the WP3 cells, validated against [62] Fig. 4 in the nanosphere limit; Tier 2b = density matrix, only for coherent features (needs a coherent collision-assisted absorption term for films first).
- **Damping vs cooling** (D15): D (e–ph coupling) and the impurity damping model are free per-sample parameters in WP4 (M0, M1), prior D ≥ 8.6–17.2 eV; M1 also gets an Umklapp fraction for the e–e share of the damping; M0 may use the generalised Fermi-liquid e–e rate.
- Data table: [62] (simulated benchmark) and Herbst et al. 2026 (fieldoscopy of ITO nanocrystals under two-cycle pulses; not yet read) added. Risks: new row for the unsettled e–ph coupling; few-cycle row mitigated by Tier 2a. Competitors: Scalora group's 2025 hydrodynamic model of [20].
- Not changed but worth knowing: the 2026 parameter set gives an ITO film with less than half the loss of [15]'s; WP1's validation against [15] keeps the [14, 15] set, which [15] itself used.

### 4.5 Open items added (on top of §2.5)

1. Ask the authors (with [14–16]): T̃_e(t) of Figs. 4–5, ε∞ and the screening constant, band top and grid, the ½ of (S2) vs (S7), the intensity convention, and which parameter set applies to which sample.
2. Tier 2a when WP3 has a field solver: reproduce [62] Fig. 4a,b in the quasi-static sphere limit.
3. Read Herbst et al., Adv. Sci. 13, e16818 (2026) and Scalora et al., Phys. Rev. A 112, 013502 (2025) before WP3's validation ladder is fixed; both are cited by [62] and bear on [20] and on few-cycle data.
4. When fitting recovery data (WP3 step 4, WP4), report results for both damping splits until the data decide.

## Session 5 — 2026-11-24 — where machine learning could enter the few-cycle work (discussion only)

**Scope requested.** The project lead asked how deep learning or generative models could be brought into the study of ITO electron dynamics under few-cycle pulses. Advice only; the plan is unchanged.

**Recommendation given (my assessment, not a decision):**

- ML earns a place where the physics model is too expensive for the loop it has to run in, where the observable is high-dimensional and the inverse problem ill-posed, or where the constitutive law itself is unknown. Few-cycle pumping meets all three; the ≥20 fs regime met none (Gate 1).
- Ranked options: (1) a learned, causal few-cycle constitutive law (field history → polarisation/current and absorbed energy) trained on density-matrix runs, with causality, passivity, the Lindhard/Drude linear limit and the Tier 2a thermal limit built in, embedded in the WP3 field solver: the condition under which the deferred operator returns (plan, deferred section); (2) amortised simulation-based inference (normalising-flow or diffusion posteriors, neural ratio estimation for M0–M4) on field-resolved data; (3) a diffusion prior for pump, probe and complex gate in double-blind FROG (WP5), building on [58, 59]; (4) equation discovery (sparse or symbolic regression, universal differential equations) for the current equation where the Fermi-sphere displacement reaches 0.5–1.3 k_F.
- Not recommended: a surrogate of the collision integral (Gate 1), black-box pulse-to-ΔR/R regression on the few available datasets, PINNs for the stiff integro-differential kinetic equation.
- Deciding experiment before any of (1): run the density matrix of [62] in its own uniform-field setting over field strength, duration and CEP, and check whether the Tier 2a thermal model misses anything a field-resolved measurement could see. If it does not, (1) shrinks to the coherent features.
- Literature checked in a quick search (not added to the library): Borthakur, Hirschman & Carbajo, arXiv:2511.08764 (diffusion pulse retrieval from partial FROG traces); Herrera Rodríguez & Kananenka, arXiv:2409.11320 (transformer for non-Markovian dissipative dynamics); simulation-based inference with neural posterior estimation in X-ray and galaxy spectroscopy (arXiv:2506.05911, 2211.09126); no simulation-based inference on ultrafast pump-probe data and no learned constitutive law for a hot-electron ENZ medium turned up.

**Open item.** The project lead decides whether to add an ML work package (gated on Tier 2a and the deciding experiment above) in a revision 2.2.

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
| D9 | 2026-11-17 | **Decided by the project lead: drop the machine-learning surrogate for now.** WP3 runs on the exact integral B1; B0 is a comparison model only; the surrogate returns only under the conditions listed in the plan's deferred section | Gate 1: B1 is only 6–46× costlier than B0 at the resolution needed, B0 misses the transient by 7–221 % |
| D11 | 2026-11-17 | Plan revised to revision 2; revision 1 archived | D9 and the findings of the literature pass |
| D12 | 2026-11-17 | The exact integral is ported to JAX (64-bit) for GPU use and gradients; numba stays the CPU reference | WP3 and WP4 need GPU batching over cells and gradients |
| D10 | 2026-11-17 | Dataset generator rejects cases whose nominal deposit exceeds half the band's capacity | a third of the plan's sampling range saturates a single 4 eV band |
| D14 | 2026-11-18 | Un, Sarkar & Sivan, Nanophotonics 2026 added to the reference list as [62]; plan revision 2.1 | supplied by the project lead; the only collision-resolved single-cycle calculation for a TCO |
| D15 | 2026-11-18 | Deformation potential D (e–ph coupling) and the impurity damping model (constant or Brooks–Herring) are inferred per sample in WP4, prior D ≥ 8.6–17.2 eV; Brooks–Herring damping and the [62] preset added to WP1 (off by default) | the authors' two parameter sets differ by 4 in G; recovery 0.5 vs 1.9 ps for the same pulse (`scripts/check_un2026.py`) |
| D16 | 2026-11-18 | Tier 2 split: 2a = field-resolved J·E drive of the WP3 cells, validated against [62] Fig. 4; 2b = density matrix only for coherent features. Flagship scope stays ≥20 fs until 2a passes | [62]: thermal model driven by instantaneous Ṗ·E reproduces the density matrix; its finite-size absorption does not transfer to films |
| D13 | 2026-11-17 | Grid rule for WP3 onward: energy step 31 meV (32 points per eV); band width 4, 6 or 8 eV chosen so that the top-cell occupation stays below 4×10⁻⁵; linear optics calibrated on the grid in use. Replaces the intensity threshold of D8 | convergence sweep, §2.3: step error ≤ 1 %, band cut-off error up to 12 % |
