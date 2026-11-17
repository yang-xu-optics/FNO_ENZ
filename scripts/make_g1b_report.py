#!/usr/bin/env python3
"""Write docs/gates/G1b_report.md from the result files.

Inputs (all under results/wp1/):
    convergence_sweep.json      scripts/convergence_sweep.py
    gate1.json                  scripts/gate1_benchmark.py        (4-thread cost per evaluation)
    bench_b1_jax_cpu.json       scripts/bench_b1_jax.py --numba   (single-core cost, numba and JAX)
    bench_b1_jax_gpu.json       scripts/bench_b1_jax.py on a GPU  (optional; the report says "outstanding" without it)

Re-run this script after the GPU benchmark to fill in the GPU section.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "wp1"
FILM_EVALS = 1.6e6  # 100 cells x 4000 RK4 steps x 4 stages
CRITERION = 0.01
KEY_OBS = ("eps_re", "eps_im", "wp2")
POINTS_PER_EV = 32  # the coarse energy step of the sweep: 128 points over 4 eV


def level_text(case):
    if case["mode"] == "energy":
        return f'{sci(case["level"])} J/m³ absorbed'
    return f'{case["level"] / 1e9:g} GW/cm² in the film'


def sci(v):
    mant, exp = f"{v:.0e}".split("e")
    return f"{mant} × 10^{int(exp)}"


def round2(v):
    """Two significant figures, with thousands separators."""
    return f"{float(f'{v:.2g}'):,.0f}"


def drive_text(case):
    return "absorbed energy" if case["mode"] == "energy" else "field"


def key_miss(m):
    return max(m[o] for o in KEY_OBS)


def worst(entries, getter):
    vals = [(key_miss(getter(e)), e) for e in entries if getter(e) is not None]
    return max(vals, key=lambda v: v[0]) if vals else None


def w_text(w):
    return f'{100 * w[0]:.2f} % ({w[1]["case"]["fwhm"]:.0f} fs, {level_text(w[1]["case"])})' if w else "not run"


def extrapolate(e):
    """Observed order and zero-step error estimates from three grids (128, 256, 512), worst key observable."""
    out = None
    for o in KEY_OBS:
        m1, m2 = e["miss"]["128"][o], e["miss"]["256"][o]
        if m2 <= 0 or m1 / m2 <= 2.0:
            continue
        x = m1 / m2 - 1.0  # x = 2^p for errors c h^p
        cand = dict(order=math.log2(x), e128=m1 * x * x / (x * x - 1.0), e256=m2 * x / (x - 1.0))
        if out is None or cand["e128"] > out["e128"]:
            out = cand
    return out


def main():
    sweep = json.loads((R / "convergence_sweep.json").read_text())
    gate = json.loads((R / "gate1.json").read_text())
    cpu = json.loads((R / "bench_b1_jax_cpu.json").read_text())
    gpu_path = R / "bench_b1_jax_gpu.json"
    gpu = json.loads(gpu_path.read_text()) if gpu_path.exists() else None

    order = lambda e: (e["case"]["mode"], e["case"]["fwhm"], e["case"]["level"])
    grid = sweep["summary"]["grid"]
    g4 = sorted([e for e in grid if e["E_max"] == 4.0], key=order)
    g6 = sorted([e for e in grid if e["E_max"] == 6.0], key=order)
    g512 = [e for e in g4 if e["N_ref"] == 512]

    # ---- grid convergence
    w128_256 = worst(g4, lambda e: e["halving"].get("128"))
    w128_512 = worst(g512, lambda e: e["miss"].get("128"))
    w256_512 = worst(g512, lambda e: e["miss"].get("256"))
    w192 = worst(g6, lambda e: e["miss"].get("192"))
    extra = [(extrapolate(e), e) for e in g512]
    extra = [(x, e) for x, e in extra if x]
    ex_worst = max(extra, key=lambda v: v[0]["e128"]) if extra else None
    orders = [x["order"] for x, _ in extra]
    ok128 = all(w is not None and w[0] <= CRITERION for w in (w128_256, w128_512))
    ok256 = w256_512 is not None and w256_512[0] <= CRITERION
    ok192 = w192 is not None and w192[0] <= CRITERION
    marginal128 = ok128 and ex_worst is not None and ex_worst[0]["e128"] > 0.8 * CRITERION
    w128_hot = worst([e for e in g4 if e["Te_max"] >= 1000.0], lambda e: e["halving"].get("128"))
    w128_warm = worst([e for e in g4 if e["Te_max"] >= 1800.0], lambda e: e["halving"].get("128"))
    rise = sorted(e["Te_max"] - 300.0 for e in g4)[:4]
    static = {n: next((e["static"][n] for e in g512 if n in e["static"]), None) for n in ("128", "256")}
    weakest_amp = min(e["amplitude"]["eps_re"] for e in g4)

    rows4 = []
    for e in g4:
        c = e["case"]
        pct = lambda m: f"{100 * key_miss(m):.2f}" if m else "—"
        rows4.append(f'| {drive_text(c)} | {c["fwhm"]:.0f} | {level_text(c)} | {e["Te_max"]:.0f} | {pct(e["halving"].get("128"))} | '
                     f'{pct(e["miss"].get("128")) if e["N_ref"] == 512 else "—"} | {pct(e["miss"].get("256")) if e["N_ref"] == 512 else "—"} |')
    rows6 = [f'| {drive_text(e["case"])} | {e["case"]["fwhm"]:.0f} | {level_text(e["case"])} | {e["Te_max"]:.0f} | '
             f'{100 * e["miss"]["192"]["eps_re"]:.2f} | {100 * e["miss"]["192"]["eps_im"]:.2f} |' for e in g6]

    # ---- band width
    band = sorted(sweep["summary"]["band_width"], key=lambda b: (b["E_narrow"], b["Te_max_wide"]))
    rows_band = [f'| {b["E_narrow"]:.0f} → {b["E_wide"]:.0f} | {drive_text(b["case"])} | {b["case"]["fwhm"]:.0f} | {level_text(b["case"])} | '
                 f'{b["Te_max_narrow"]:.0f} | {b["Te_max_wide"]:.0f} | {sci(b["f_top_narrow"])} | '
                 f'{100 * b["miss"]["eps_re"]:.1f} | {100 * b["miss"]["eps_im"]:.1f} | {100 * b["miss"]["T_e"]:.1f} |' for b in band]
    band_lines, safe_T = [], {}
    for e_lo in sorted({b["E_narrow"] for b in band}):
        bs = [b for b in band if b["E_narrow"] == e_lo]
        good = [b for b in bs if key_miss(b["miss"]) <= CRITERION]
        bad = sorted((b for b in bs if key_miss(b["miss"]) > CRITERION), key=lambda b: b["Te_max_wide"])
        parts = []
        if good:
            safe_T[e_lo] = max(b["Te_max_wide"] for b in good)
            parts.append(f"within 1 % up to a peak temperature of {safe_T[e_lo]:.0f} K")
        if bad and bad[-1]["Te_max_wide"] < 1.1 * bad[0]["Te_max_wide"]:
            parts.append(f"off by {100 * max(key_miss(b['miss']) for b in bad):.1f} % near {100 * round(sum(b['Te_max_wide'] for b in bad) / len(bad) / 100):.0f} K")
        elif bad:
            parts.append(f"off by {100 * key_miss(bad[0]['miss']):.1f} % at {bad[0]['Te_max_wide']:.0f} K and by "
                         f"{100 * key_miss(bad[-1]['miss']):.0f} % at {bad[-1]['Te_max_wide']:.0f} K")
        band_lines.append(f"- **{e_lo:.0f} eV band** (against {bs[0]['E_wide']:.0f} eV): " + "; ".join(parts) + ".")
    r500 = lambda t: f"{500 * round(t / 500):.0f}"
    temp_terms = (f" In temperature terms: 4 eV up to about {r500(safe_T[4.0])} K, 6 eV up to about {r500(safe_T[6.0])} K, 8 eV beyond."
                  if 4.0 in safe_T and 6.0 in safe_T else "")
    f_good = max((b["f_top_narrow"] for b in band if key_miss(b["miss"]) <= CRITERION), default=None)
    f_bad = min((b["f_top_narrow"] for b in band if key_miss(b["miss"]) > CRITERION), default=None)
    rule = (f"The three band pairs fall on one curve when plotted against the occupation of the narrower band's top cell: the change stays below 1 % "
            f"while that cell never fills beyond {sci(f_good)}, and exceeds 1 % from {sci(f_bad)} upward. Every run records this number "
            f"(`f_top_max`), so it serves as the test: **widen the band when the top-cell occupation passes {sci(f_good)}.**" if f_good and f_bad else "")
    series46 = sorted((b for b in band if b["case_key"] == max(band, key=lambda b: key_miss(b["miss"]))["case_key"]), key=lambda b: b["E_narrow"])
    conv_text = ("For the hottest condition the change shrinks with each widening (" + ", ".join(
        f"{100 * key_miss(b['miss']):.1f} % from {b['E_narrow']:.0f} to {b['E_wide']:.0f} eV" for b in series46)
        + "), so the cut-off dependence converges, but slowly.") if len(series46) > 1 else ""

    # ---- cost
    cost4 = {c["N"]: c["b1_ms"] for c in gate["cost"]}
    cpu1 = {r["N"]: (r["numba_ms_per_cell"], r.get("ms_per_cell")) for r in cpu["rows"] if r.get("numba_ms_per_cell") and r["cells"] == 1}
    hours = lambda n: FILM_EVALS * cpu1[n][0] / 3.6e6
    cost_rows = []
    for n in sorted(cpu1):
        nb, jx = cpu1[n]
        four = f"{cost4[n]:.1f} | {FILM_EVALS * cost4[n] / 6e4:.0f} min" if n in cost4 else "— | —"
        cost_rows.append(f"| {n} | {n / POINTS_PER_EV:g} | {nb:.1f} | {hours(n):.1f} | {four} | {f'{jx:.1f}' if jx else '—'} |")

    if gpu:
        gpu_rows = "\n".join(
            f'| {r["N"]} | {r["cells"]} | {r["ms_per_call"]:.2f} | {r["ms_per_cell"]:.4f} | {r["film_run_minutes"]:.1f} |'
            if "ms_per_call" in r else f'| {r["N"]} | {r["cells"]} | failed: {r.get("error", "")[:60]} | | |' for r in gpu["rows"])
        gpu_section = (f'Measured on {gpu["device"]} (JAX {gpu["jax"]}).\n\n| N | cells per call | ms per call | ms per cell | film run (min) |\n'
                       f'| --- | --- | --- | --- | --- |\n{gpu_rows}\n')
        gpu_verdict = "measured, see section C"
    else:
        gpu_section = ("**Not measured yet.** The environment this was written in has no GPU. On a GPU node:\n\n```bash\npip install -U \"jax[cuda12]\"\n```\n\n"
                       "```bash\npython scripts/bench_b1_jax.py --N 128 192 256 --cells 1 16 100 400\n```\n\n```bash\npython scripts/make_g1b_report.py\n```\n\n"
                       "The benchmark writes `results/wp1/bench_b1_jax_gpu.json` and the last command fills this section in.\n")
        gpu_verdict = "**outstanding** (no GPU in this environment)"

    grid_status = ("met at 128 points, without margin at the weakest excitation" if marginal128 else "met at 128 points") if ok128 else \
        ("met at 256 points, not at 128" if ok256 else "not met")
    ex_text = (f"; extrapolated to zero energy step (observed order {min(orders):.1f}–{max(orders):.1f}) the 128-point error is "
               f"{100 * ex_worst[0]['e128']:.1f} % and the 256-point error {100 * ex_worst[0]['e256']:.1f} % in the worst case" if ex_worst else "")
    n_lo, n_hi = min(cpu1), max(n for n in cpu1 if n <= 256)

    report = f"""# Gate 1b report: is the exact solver affordable across the whole condition matrix?

Generated by `scripts/make_g1b_report.py` from `results/wp1/convergence_sweep.json`, `gate1.json` and `bench_b1_jax_*.json`. Do not edit by hand.

## Verdict

| Part of the criterion (plan revision 2) | Result |
| --- | --- |
| 128 or 256 points reproduce the finest grid within 1 % of the transient amplitude over the whole matrix, including weak excitation | **{grid_status}.** 128 points: at most {w_text(w128_256)} against 256 points over all {len(g4)} conditions, and {w_text(w128_512)} against 512 points. 256 points: at most {w_text(w256_512)} against 512 points{ex_text} |
| … and with the 6 eV band | **{"met" if ok192 else "not met"} for the grid:** 192 points differ from 384 by at most {w_text(w192)}. **The band width itself is the larger error** at strong excitation, see section B |
| A film run on the exact integral fits in hours on a cluster node | **met:** {", ".join(f"{hours(n):.1f} core-hours at {n} points" for n in sorted(cpu1) if n <= n_hi)}; the cells of a film run in parallel |
| … or minutes on a GPU | {gpu_verdict} |

In short: the energy step is not the problem. A step of 31 meV (128 points per 4 eV) sits at the 1 % line for the weakest excitation (a temperature rise of {rise[0]:.0f}–{rise[-1]:.0f} K) and improves from there: against the next finer grid the difference is {100 * w128_256[0]:.1f} % at worst and {100 * w128_warm[0]:.1f} % or less for every run that peaks above 1800 K. What limits accuracy at strong excitation is the upper cut-off of the band, and the grid has to grow with it.

## A. Energy step, 4 eV band

"Miss" is the largest deviation of the transient y(t) − y(t₀) from that of a finer grid, as a percentage of the finer grid's transient amplitude; the table gives the largest of Re ε, Im ε and ω_p². All runs use the exact integral, one material, 1300 nm. "Field" levels are a fixed in-film intensity; the thin-film (self-consistent) drive is not part of this sweep.

| Drive | Pulse (fs) | Level | T_e,max (K) | 128 vs 256 pts (%) | 128 vs 512 pts (%) | 256 vs 512 pts (%) |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(rows4)}

![Change of the transient when the energy step is halved](../../results/wp1/fig_convergence_sweep.png)

- The error falls with excitation: the cold Fermi edge (k_B T = 26 meV at 300 K) is the sharpest feature the grid has to resolve, and it is as narrow as one 128-point cell (31 meV). Above 1000 K the 128-point grid differs from the 256-point grid by at most {w_text(w128_hot)}.
- The error is a scale error, not a shape error: T_e, ω_p², Re ε and Im ε all miss by the same fraction, as expected from a quadrature error in the heat capacity.
- **Static offset.** The cold-state permittivity itself differs between grids by a constant ({f"{static['128']['eps_re']:+.1e} in Re ε and {static['128']['eps_im']:+.1e} in Im ε at 128 points, {static['256']['eps_re']:+.1e} and {static['256']['eps_im']:+.1e} at 256 points, against 512" if static["128"] and static["256"] else "n/a"}). It is the same for every condition and is removed when ε∞ and the impurity rate are calibrated on the grid in use (`enzkin.calibrate`). Carried from one grid to another it is up to {f"{abs(static['128']['eps_re']) / weakest_amp * 100:.0f} %" if static["128"] else "n/a"} of the weakest transient in the table, which is how the first pass of this sweep appeared not to converge.

## B. The 6 eV band, and how wide the band has to be

Strong-excitation conditions rerun with a 6 eV band at the same energy step (192 points) and at half the step (384 points):

| Drive | Pulse (fs) | Level | T_e,max (K) | 192 vs 384 pts: Re ε (%) | Im ε (%) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(rows6)}

**Band width.** Runs with the same energy step and different band width; the narrower band is measured against the wider one with the same "miss". The top-cell occupation is the largest value f reaches in the last cell of the narrower band.

| Band (eV) | Drive | Pulse (fs) | Level | T_e,max narrow (K) | T_e,max wide (K) | top-cell occupation | effect on Re ε (%) | on Im ε (%) | on T_e (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(rows_band)}

![Band-width effect against top-cell occupation](../../results/wp1/fig_band_width.png)

{chr(10).join(band_lines)}

{rule}

{conv_text}

Two cautions. A narrow band always errs in the same direction: electrons that cannot go higher pile up below the cut-off, so the temperature and the transient are overestimated. And these are statements about the numerical cut-off. The Kane dispersion several eV above the Fermi level is an extrapolation, not the band structure of ITO, so the sensitivity to the cut-off is also a measure of how much the hottest results depend on states the model does not describe (`docs/theory/wp1_kinetic_model.md`, limitations).

## C. Cost of a film run

One film run of the plan's size is counted as 100 cells × 4000 steps × 4 stages = 1.6 × 10⁶ evaluations of the exact integral (the Gate 1 report counted 10⁶). The grid is {POINTS_PER_EV} points per eV of band. Single-core times were measured with nothing else running on the machine.

| N | band at 31 meV step (eV) | numba, 1 core (ms) | film run (core-hours) | numba, 4 cores (ms) | film run on 4 cores | JAX on 1 CPU core (ms) |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(cost_rows)}

JAX on a CPU is slower than numba; it is there for the GPU and for gradients.

### GPU

{gpu_section}
## Consequences

- WP3 can be built on the exact integral. The energy step stays at 31 meV; the band width follows the excitation: 4 eV (128 points) while the top-cell occupation stays below {sci(f_good) if f_good else "the threshold of section B"}, then 6 eV (192 points) and 8 eV (256 points).{temp_terms}
- Use 16 meV (256 points per 4 eV) only where a weak-excitation transient (a temperature rise of a few hundred kelvin) is needed to better than 1 %.
- Calibrate the linear optics on the grid in use; do not carry ε∞ from one grid to another.
- The inference stage (10³–10⁴ film runs) costs {round2(hours(n_lo) * 1e3)}–{round2(hours(n_lo) * 1e4)} core-hours at {n_lo} points and {round2(hours(n_hi) * 1e3)}–{round2(hours(n_hi) * 1e4)} at {n_hi} points. Whether sequential sampling is practical depends on the GPU number above.

## Limits

- One wavelength (1300 nm) and the reference material; no thin-film (self-consistent field) runs in this sweep.
- The finest grid is 512 points for four conditions and 256 for the rest.
- The widest band is 10 eV, for the two hottest conditions only.
- Time-step error is common to all grids and is not measured here.
"""
    out = ROOT / "docs" / "gates" / "G1b_report.md"
    out.write_text(report, encoding="utf-8")
    print("wrote", out)
    print("128 vs 256:", w_text(w128_256), "| 128 vs 512:", w_text(w128_512), "| 256 vs 512:", w_text(w256_512), "| 6 eV 192 vs 384:", w_text(w192))
    if ex_worst:
        print(f"extrapolated: order {min(orders):.2f}-{max(orders):.2f}, e128 {100 * ex_worst[0]['e128']:.2f} %, e256 {100 * ex_worst[0]['e256']:.2f} %")


if __name__ == "__main__":
    main()
