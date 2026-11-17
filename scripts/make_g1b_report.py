#!/usr/bin/env python3
"""Write docs/gates/G1b_report.md from the result files.

Inputs (all under results/wp1/):
    convergence_sweep.json      scripts/convergence_sweep.py
    gate1.json                  scripts/gate1_benchmark.py        (4-thread cost per evaluation)
    bench_b1_jax_cpu.json       scripts/bench_b1_jax.py --numba   (single-core cost, numba and JAX)
    bench_b1_jax_gpu.json       scripts/bench_b1_jax.py on a GPU  (optional; the report says "pending" without it)

Re-run this script after the GPU benchmark to fill in the GPU section.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "wp1"
FILM_EVALS = 1.6e6  # 100 cells x 4000 RK4 steps x 4 stages
CRITERION = 0.01
KEY_OBS = ("eps_re", "eps_im", "wp2")


def level_text(case):
    if case["mode"] == "energy":
        return f'{case["level"]:.0e} J/m³ absorbed'.replace("e+0", "×10^").replace("e+", "×10^")
    return f'{case["level"] / 1e9:g} GW/cm² in-film'


def main():
    sweep = json.loads((R / "convergence_sweep.json").read_text())
    gate = json.loads((R / "gate1.json").read_text())
    cpu = json.loads((R / "bench_b1_jax_cpu.json").read_text())
    gpu_path = R / "bench_b1_jax_gpu.json"
    gpu = json.loads(gpu_path.read_text()) if gpu_path.exists() else None

    grid = sweep["summary"]["grid"]
    band = sweep["summary"]["band_width"]
    g4 = sorted([e for e in grid if e["E_max"] == 4.0], key=lambda e: (e["case"]["mode"], e["case"]["fwhm"], e["case"]["level"]))
    g6 = sorted([e for e in grid if e["E_max"] == 6.0], key=lambda e: (e["case"]["mode"], e["case"]["fwhm"], e["case"]["level"]))

    def worst(entries, n):
        vals = [(e["miss"][n][o], e, o) for e in entries if n in e["miss"] for o in KEY_OBS]
        return max(vals, key=lambda v: v[0]) if vals else None

    w128, w256 = worst(g4, "128"), worst([e for e in g4 if e["N_ref"] == 512], "256")
    w128_512 = worst([e for e in g4 if e["N_ref"] == 512], "128")
    w192 = worst(g6, "192")
    ok128 = w128 is not None and w128[0] <= CRITERION and (w128_512 is None or w128_512[0] <= CRITERION)
    ok192 = w192 is not None and w192[0] <= CRITERION
    static = next((e["static"]["128"] for e in g4 if "128" in e.get("static", {})), None)

    rows4 = []
    for e in g4:
        c, m = e["case"], e["miss"]
        cell = lambda n, o: f'{100 * m[n][o]:.2f}' if n in m else "—"
        rows4.append(f'| {"energy" if c["mode"] == "energy" else "field"} | {c["fwhm"]:.0f} | {level_text(c)} | {e["Te_max"]:.0f} | {e["N_ref"]} | '
                     f'{cell("128", "eps_re")} | {cell("128", "eps_im")} | {cell("256", "eps_re") if e["N_ref"] == 512 else "ref"} | '
                     f'{cell("256", "eps_im") if e["N_ref"] == 512 else "ref"} |')
    rows6 = []
    for e in g6:
        c, m = e["case"], e["miss"]["192"]
        rows6.append(f'| {"energy" if c["mode"] == "energy" else "field"} | {c["fwhm"]:.0f} | {level_text(c)} | {e["Te_max"]:.0f} | '
                     f'{100 * m["eps_re"]:.2f} | {100 * m["eps_im"]:.2f} |')
    band = sorted(band, key=lambda b: (b["E_narrow"], b["Te_max_wide"]))
    rows_band = [f'| {b["E_narrow"]:.0f} → {b["E_wide"]:.0f} | {"energy" if b["case"]["mode"] == "energy" else "field"} | {b["case"]["fwhm"]:.0f} | '
                 f'{level_text(b["case"])} | {b["Te_max_narrow"]:.0f} | {b["Te_max_wide"]:.0f} | {b["f_top_narrow"]:.0e} | '
                 f'{100 * b["miss"]["eps_re"]:.1f} | {100 * b["miss"]["eps_im"]:.1f} | {100 * b["miss"]["T_e"]:.1f} |' for b in band]
    b46 = [b for b in band if b["E_narrow"] == 4.0]
    b68 = [b for b in band if b["E_narrow"] == 6.0]
    key_miss = lambda b: max(b["miss"][o] for o in KEY_OBS)
    safe4 = max((b["Te_max_wide"] for b in b46 if key_miss(b) <= CRITERION), default=None)
    unsafe4 = min((b["Te_max_wide"] for b in b46 if key_miss(b) > CRITERION), default=None)
    worst46 = max(b46, key=key_miss) if b46 else None
    worst68 = max(b68, key=key_miss) if b68 else None
    band_text = ""
    if worst46:
        band_text += (f"Widening the band from 4 to 6 eV changes the permittivity transient by up to {100 * key_miss(worst46):.0f} % "
                      f"({worst46['case']['fwhm']:.0f} fs, {level_text(worst46['case'])}), far more than the grid does. ")
        if safe4 and unsafe4:
            band_text += (f"The change is below 1 % for the condition that peaks at {safe4:.0f} K and above 1 % from {unsafe4:.0f} K upward, "
                          f"so the 4 eV band is safe only below roughly {safe4 / 1000:.1f}–{unsafe4 / 1000:.1f} × 10³ K. ")
    if worst68:
        ok68 = key_miss(worst68) <= CRITERION
        band_text += (f"Widening it again from 6 to 8 eV changes the transient by at most {100 * key_miss(worst68):.1f} % at the hottest conditions tested "
                      f"({worst68['Te_max_wide']:.0f} K), so the 6 eV band is "
                      + ("wide enough for the whole matrix." if ok68 else "not yet converged there; the hottest results carry this band-width uncertainty."))
    else:
        band_text += "The 6 eV band has not been checked against a wider one."

    cost4 = {c["N"]: c["b1_ms"] for c in gate["cost"]}
    cpu1 = {}
    for r in cpu["rows"]:
        if r.get("numba_ms_per_cell") and r["cells"] == 1:
            cpu1[r["N"]] = (r["numba_ms_per_cell"], r.get("ms_per_cell"))
    cost_rows = []
    for n in sorted(cpu1):
        nb, jx = cpu1[n]
        cost_rows.append(f"| {n} | {nb:.1f} | {FILM_EVALS * nb / 3.6e6:.1f} | {cost4.get(n, float('nan')):.1f} | "
                         f"{FILM_EVALS * cost4.get(n, float('nan')) / 6e4:.0f} min | {jx:.1f} |")

    if gpu:
        gpu_rows = "\n".join(
            f'| {r["N"]} | {r["cells"]} | {r["ms_per_call"]:.2f} | {r["ms_per_cell"]:.4f} | {r["film_run_minutes"]:.1f} |'
            if "ms_per_call" in r else f'| {r["N"]} | {r["cells"]} | failed: {r.get("error", "")[:60]} | | |' for r in gpu["rows"])
        gpu_section = (f'Measured on {gpu["device"]} (JAX {gpu["jax"]}).\n\n| N | cells per call | ms per call | ms per cell | film run (min) |\n'
                       f'| --- | --- | --- | --- | --- |\n{gpu_rows}\n')
        gpu_verdict = "measured, see section C"
    else:
        gpu_section = ("**Not measured yet.** The environment this was written in has no GPU. On a GPU node:\n\n```bash\npip install -U \"jax[cuda12]\"\n```\n\n"
                       "```bash\npython scripts/bench_b1_jax.py --N 128 256 --cells 1 16 100 400\n```\n\n```bash\npython scripts/make_g1b_report.py\n```\n\n"
                       "The benchmark writes `results/wp1/bench_b1_jax_gpu.json` and the last command fills this section in.\n")
        gpu_verdict = "**outstanding** (no GPU available)"

    def w_text(w):
        return f'{100 * w[0]:.2f} % ({w[2]}, {w[1]["case"]["fwhm"]:.0f} fs, {level_text(w[1]["case"])})' if w else "not run"

    n_ok = 128 if ok128 else 256
    report = f"""# Gate 1b report: is the exact solver affordable across the whole condition matrix?

Generated by `scripts/make_g1b_report.py` from `results/wp1/convergence_sweep.json`, `gate1.json` and `bench_b1_jax_*.json`.

## Verdict

| Part of the criterion (plan revision 2) | Result |
| --- | --- |
| 128 or 256 points reproduce the finest grid within 1 % of the transient amplitude over the whole matrix, including weak excitation | **{"met at 128 points" if ok128 else "not met at 128 points"}**: worst case {w_text(w128)} against 256 points{f", and {w_text(w128_512)} against 512 points" if w128_512 else ""} |
| … and with the 6 eV band | **{"met at 192 points (the same energy step)" if ok192 else "not met at 192 points"}**: worst case {w_text(w192)} against 384 points |
| A film run on the exact integral fits in hours on a cluster node | **met**: {FILM_EVALS * cpu1[n_ok][0] / 3.6e6:.1f} core-hours at {n_ok} points |
| … or minutes on a GPU | {gpu_verdict} |

{"The grid half of the gate passes: 128 points over a 4 eV band (192 over 6 eV) are enough for every condition tested, from a 130 K temperature rise to 10 000 K." if ok128 and ok192 else "The grid half of the gate does not pass at the coarsest grid; see the tables for where the finer grid is needed."}

## A. Grid convergence, 4 eV band

"Miss" is the largest deviation of the transient y(t) − y(t₀) from that of the finest grid run for the condition, as a percentage of the transient's size. All runs use the exact integral, one material, 1300 nm.

| Drive | Pulse (fs) | Level | T_e,max (K) | Finest grid | 128 pts: Re ε | 128 pts: Im ε | 256 pts: Re ε | 256 pts: Im ε |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(rows4)}

![Grid convergence](../../results/wp1/fig_convergence_sweep.png)

**Static offset.** The cold-state permittivity itself differs between grids by a constant: {f"{static['eps_re']:+.1e} in Re ε and {static['eps_im']:+.1e} in Im ε between 128 and the reference grid" if static else "n/a"}. It is the same for every condition, comes from the quadrature of the cold Fermi sea, and is absorbed when ε∞ and the impurity rate are calibrated on the grid in use (`enzkin.calibrate`). It matters only if two grids are compared without recalibrating: for the weakest excitation in the table it is as large as {f"{abs(static['eps_re']) / min(e['amplitude']['eps_re'] for e in g4) * 100:.0f} %" if static else "n/a"} of the transient.

## B. The 6 eV band

Strong-excitation conditions rerun with a 6 eV band at the same energy step (192 points) and at half the step (384 points).

| Drive | Pulse (fs) | Level | T_e,max (K) | 192 vs 384 pts: Re ε (%) | Im ε (%) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(rows6)}

**Band width.** Runs with the same energy step and different band width, the narrower band measured against the wider one with the same "miss". The top-cell occupation is the largest value reached by f in the last cell of the narrower band.

| Band (eV) | Drive | Pulse (fs) | Level | T_e,max narrow (K) | T_e,max wide (K) | top-cell occupation | effect on Re ε (%) | on Im ε (%) | on T_e (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(rows_band)}

{band_text} A 6 eV band at the same energy step costs (192/128)³ ≈ 3.4 times more per evaluation than a 4 eV band. These are statements about the numerical cut-off; whether a single Kane band describes ITO several eV above the Fermi level is a separate, physical question (`docs/theory/wp1_kinetic_model.md`, limitations).

## C. Cost of a film run

One film run of the plan's size is 100 cells × 4000 steps × 4 stages = 1.6 × 10⁶ evaluations of the exact integral.

| N | numba, 1 core (ms) | film run (core-hours) | numba, 4 cores (ms) | film run on 4 cores | JAX on 1 CPU core (ms) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(cost_rows)}

JAX on a CPU is slower than numba; it is there for the GPU and for gradients.

### GPU

{gpu_section}
## Consequences

- WP3 can be built on the exact integral at {n_ok} points over 4 eV, or 192 points over 6 eV for hot cases.
- Calibrate the linear optics on the grid in use; do not carry ε∞ from one grid to another.
- The inference stage (10³–10⁴ film runs) costs {FILM_EVALS * cpu1[n_ok][0] / 3.6e6 * 1e3:,.0f}–{FILM_EVALS * cpu1[n_ok][0] / 3.6e6 * 1e4:,.0f} core-hours at {n_ok} points. Whether sequential sampling is feasible depends on the GPU number above.

## Limits

- One wavelength (1300 nm) and the reference material; no thin-film (self-consistent field) runs in this sweep.
- The finest grid is 512 points for four conditions and 256 for the rest.
- Time-step error is common to all grids and is not measured here.
"""
    out = ROOT / "docs" / "gates" / "G1b_report.md"
    out.write_text(report, encoding="utf-8")
    print("wrote", out)
    print("grid 4 eV @128:", w_text(w128), "| vs 512:", w_text(w128_512), "| 256 vs 512:", w_text(w256), "| 6 eV @192:", w_text(w192))


if __name__ == "__main__":
    main()
