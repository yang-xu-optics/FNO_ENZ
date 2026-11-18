# FNO_ENZ: notes for working sessions

Research project: what sets the speed and sign of the optical response of ENZ indium tin oxide under intense ultrashort pulses. Plan: `docs/PROJECT_PLAN.md` (revision 2.1: kinetic solver, Maxwell-kinetic film model, mechanism inference, experiment; the neural-operator surrogate is deferred; few-cycle pumps are Tier 2a/2b). Revision 1 is in `docs/archive/`. Start every session by reading the last entry and the open items of `docs/WORKLOG.md`.

## Standing rules

- **Keep the working log.** Append a dated session entry to `docs/WORKLOG.md` for everything done: scope, results with the script that reproduces them, dead ends, deviations from the plan, open items. Decisions that change the plan also go in its decision table.
- **Literature.** `references/LITERATURE_SUMMARY.md` is the theory cornerstone; reference numbers `[n]` are the plan's. When a paper is (re)read, update its digest there and, if it changes a conclusion, Part I. Never commit PDFs, TeX sources or extracted text (public repository); they are rebuilt with `python references/fetch_references.py --sources --text`.
- **Equations.** Any departure from a published formula goes in `docs/theory/` with the reason, and in the log.
- **Results.** A number enters a gate report (`docs/gates/`) only if a script under `scripts/` regenerates it into `results/`. Progress reports (`docs/reports/`) are built by `scripts/make_progress_report.py` from the same files.
- **Structure.** Solver code in `src/enzkin/`; later work packages add `src/enzkin/maxwell/` (WP3), `inference/` (WP4), `experiment/` (WP5). Large data go to `data/` (git-ignored).
- **No surrogate for now.** The film solver uses the exact collision integral (`ElectronElectronB1`, or `ElectronElectronB1Jax` on a GPU). The FFT integral B0 is a comparison model only. The conditions under which a learned operator returns are in the plan's deferred section.

## Commands

```bash
pip install -e ".[dev]"
```

```bash
python -m pytest
```

The full test suite takes about 35 s on 4 cores when nothing else is running; the numba kernels use all cores, so do not run benchmarks and tests at the same time.

## Things that are easy to get wrong

- Units inside `enzkin` are eV, fs, nm, K (`enzkin.constants`).
- The energy grid must stay uniform: exact energy conservation relies on i₁ + i₂ = i₃ + i₄.
- Default band width is 4 eV as printed in [14]. Keep the energy step at 31 meV (32 points per eV) and widen the band (6 eV with N = 192, 8 eV with N = 256) when the top-cell occupation `f_top_max` passes 4×10⁻⁵, roughly T_e above 3500 K for 4 eV and 6700 K for 6 eV; a narrow band overstates T_e and the transient (log §2.3, `docs/gates/G1b_report.md`).
- Two published ITO parameter sets exist and differ by four in the electron–phonon coupling: `ito_sarkar2023()` ([14, 15]: D = 17.2 eV, constant impurity rate, used for the validation against [15]) and `ito_un2026()` ([62]: D = 8.6 eV plus Brooks–Herring impurities). Say which one a result uses; the split between e–ph and impurity damping sets the cooling time (log §4.2).
- Calibrate ε∞ and γ_imp (`enzkin.calibrate`) on the grid in use; the cold permittivity shifts by about 10⁻³ between grids, which is a large fraction of a weak transient.
- `ElectronElectronB0` gets stiff (rates of several per fs) at high excitation; always integrate with the adaptive `enzkin.integrate.run`.
- Cloud sessions should use the `FNO_ENZ` environment (network allowlist for arXiv and the repositories).
- JAX is an optional dependency (`pip install -e ".[jax]"`); its tests are skipped when it is absent. Importing `enzkin.collisions.ee_jax` switches JAX to 64-bit floats.
