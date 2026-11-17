# FNO_ENZ: notes for working sessions

Research project: what sets the speed and sign of the optical response of ENZ indium tin oxide under intense ultrashort pulses. Plan: `docs/PROJECT_PLAN.md` (five work packages, four gates). Start every session by reading the last entry and the open items of `docs/WORKLOG.md`.

## Standing rules

- **Keep the working log.** Append a dated session entry to `docs/WORKLOG.md` for everything done: scope, results with the script that reproduces them, dead ends, deviations from the plan, open items. Decisions that change the plan also go in its decision table.
- **Literature.** `references/LITERATURE_SUMMARY.md` is the theory cornerstone; reference numbers `[n]` are the plan's. When a paper is (re)read, update its digest there and, if it changes a conclusion, Part I. Never commit PDFs, TeX sources or extracted text (public repository); they are rebuilt with `python references/fetch_references.py --sources --text`.
- **Equations.** Any departure from a published formula goes in `docs/theory/` with the reason, and in the log.
- **Results.** A number enters a gate report (`docs/gates/`) only if a script under `scripts/` regenerates it into `results/`. Progress reports (`docs/reports/`) are built by `scripts/make_progress_report.py` from the same files.
- **Structure.** Solver code in `src/enzkin/`; later work packages add `src/enzkin/operator/` (WP2) and `src/enzkin/maxwell/` (WP3), `inference/` (WP4), `experiment/` (WP5). Large data go to `data/` (git-ignored).

## Commands

```bash
pip install -e ".[dev]"
```

```bash
pytest
```

The full test suite takes about 35 s on 4 cores when nothing else is running; the numba kernels use all cores, so do not run benchmarks and tests at the same time.

## Things that are easy to get wrong

- Units inside `enzkin` are eV, fs, nm, K (`enzkin.constants`).
- The energy grid must stay uniform: exact energy conservation relies on i₁ + i₂ = i₃ + i₄.
- Default band width is 4 eV as printed in [14]; it saturates above about 30 GW/cm² in the geometry of [15]. Use `E_max=6.0` with N scaled to keep ΔE for high-intensity runs (log §1.5).
- `ElectronElectronB0` gets stiff (rates of several per fs) at high excitation; always integrate with the adaptive `enzkin.integrate.run`.
- Cloud sessions should use the `FNO_ENZ` environment (network allowlist for arXiv and the repositories).
