#!/usr/bin/env python3
"""Grid-convergence sweep for the exact collision integral B1 (follow-up to Gate 1, plan revision 2, G1b).

Gate 1 found that 128 energy points reproduce the 256-point permittivity transient to 0.2 %, but for one
condition only.  This sweep asks the same question over the whole condition matrix, adds the weak
excitations where the cold Fermi edge is sharpest (k_B T = 26 meV against dE = 31 meV at N = 128), and
repeats the strong cases with the 6 eV band that the validation against ref. [15] calls for.

Stages (results are written after every run, so an interrupted sweep keeps what it has):

  A   all conditions, 4 eV band, N = 128 and 256
  B   strong conditions, 6 eV band, N = 192 and 384 (same dE as 128 and 256 at 4 eV)
  C   four conditions, 4 eV band, N = 512

One material is used throughout (linear optics calibrated once, on the 256-point grid), so differences
between grids include any static offset as well as the dynamics.

"Miss" of grid N against the finest grid available for that condition, for an observable y(t):
    max_t |y_N - y_ref| / max_t |y_ref - y_ref(t_0)|

Usage:  python scripts/convergence_sweep.py [--stages A B C] [--quick]
        python scripts/convergence_sweep.py --summarise      # recompute the summary from the stored series
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from enzkin.calibrate import calibrate_linear_optics
from enzkin.drive import FieldDrive, PowerDrive, Pulse
from enzkin.integrate import run
from enzkin.material import ito_sarkar2023
from enzkin.model import KineticModel

OBS = ("eps_re", "eps_im", "wp2", "gamma", "T_e")
OUT = Path("results/wp1/convergence_sweep.json")
SERIES = Path("results/wp1/convergence_sweep_series.npz")
SAMPLE_DT = 4.0  # fs


def all_cases(quick=False):
    cases = []
    fwhms = [20.0, 220.0] if quick else [20.0, 50.0, 100.0, 220.0]
    levels = [1e6, 1e8] if quick else [1e6, 1e7, 3e7, 1e8, 3e8]
    for fwhm in fwhms:
        for u in levels:
            cases.append(dict(mode="energy", fwhm=fwhm, level=u))
    for fwhm, inten in ([(20.0, 1e11)] if quick else [(20.0, 1e11), (20.0, 5e11), (220.0, 1e10), (220.0, 4e10)]):
        cases.append(dict(mode="field", fwhm=fwhm, level=inten))
    return cases


def key(case):
    return f'{case["mode"]}_{case["fwhm"]:.0f}fs_{case["level"]:.0e}'


STRONG = [("energy", 20.0, 3e8), ("energy", 220.0, 1e8), ("energy", 220.0, 3e8), ("field", 220.0, 4e10), ("field", 20.0, 5e11)]
FINEST = [("energy", 20.0, 1e6), ("energy", 20.0, 1e8), ("energy", 220.0, 1e6), ("energy", 220.0, 3e8)]


def pick(cases, wanted):
    return [c for c in cases if (c["mode"], c["fwhm"], c["level"]) in wanted]


def one_run(mat, N, case):
    pulse = Pulse(1300.0, case["fwhm"])
    drive = PowerDrive.from_absorbed_energy(pulse, case["level"]) if case["mode"] == "energy" else FieldDrive(pulse, case["level"])
    model = KineticModel(mat, N=N, drive=drive)
    dt = min(1.0, case["fwhm"] / 40.0)
    every = int(round(SAMPLE_DT / dt))
    t_lo = -2.5 * case["fwhm"]
    span = math.ceil((max(2.5 * case["fwhm"], 150.0) + 250.0 - t_lo) / SAMPLE_DT) * SAMPLE_DT
    t0 = time.perf_counter()
    tr = run(model, t_lo, t_lo + span, dt, sample_every=every, omega_probe=pulse.omega)
    return tr, time.perf_counter() - t0


def miss(a, ref):
    amp = float(np.max(np.abs(ref - ref[0])))
    return float(np.max(np.abs(a - ref)) / amp) if amp > 0 else float("nan")


def summarise(series: dict, runs: list) -> dict:
    """Per condition and band: miss of every coarser grid against the finest one that was run."""
    out = []
    groups = {}
    for r in runs:
        groups.setdefault((r["case_key"], r["E_max"]), []).append(r)
    for (ck, emax), rs in groups.items():
        rs = sorted(rs, key=lambda r: r["N"])
        ref = rs[-1]
        if len(rs) < 2:
            continue
        entry = dict(case_key=ck, case=ref["case"], E_max=emax, N_ref=ref["N"], Te_max=ref["Te_max"],
                     amplitude={o: float(np.max(np.abs(series[f'{ref["id"]}_{o}'] - series[f'{ref["id"]}_{o}'][0]))) for o in OBS},
                     miss={})
        for r in rs[:-1]:
            entry["miss"][str(r["N"])] = {o: miss(series[f'{r["id"]}_{o}'], series[f'{ref["id"]}_{o}']) for o in OBS}
        out.append(entry)
    # band-width effect: 6 eV against 4 eV at the same dE (N = 192 vs 128), in units of the 6 eV transient
    band = []
    for (ck, emax), rs in groups.items():
        if emax != 6.0:
            continue
        r6 = next((r for r in rs if r["N"] == 192), None)
        r4 = next((r for r in groups.get((ck, 4.0), []) if r["N"] == 128), None)
        if r6 and r4:
            n = min(len(series[f'{r6["id"]}_T_e']), len(series[f'{r4["id"]}_T_e']))
            band.append(dict(case_key=ck, case=r6["case"], Te_max_4eV=r4["Te_max"], Te_max_6eV=r6["Te_max"],
                             miss={o: miss(series[f'{r4["id"]}_{o}'][:n], series[f'{r6["id"]}_{o}'][:n]) for o in OBS}))
    return dict(grid=out, band_width=band)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", nargs="*", default=["A", "B", "C"])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--summarise", action="store_true")
    args = ap.parse_args()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    state = json.loads(OUT.read_text()) if OUT.exists() else dict(runs=[])
    series = dict(np.load(SERIES)) if SERIES.exists() else {}
    if args.summarise:
        state["summary"] = summarise(series, state["runs"])
        OUT.write_text(json.dumps(state, indent=1))
        print("summary rewritten")
        return

    mat4 = calibrate_linear_optics(ito_sarkar2023(), 256)
    mat6 = mat4.with_(E_max=6.0)
    state["material"] = dict(eps_inf=mat4.eps_inf, gamma_imp=mat4.gamma_imp, calibrated_on_N=256)
    cases = all_cases(args.quick)
    plan = []
    if "A" in args.stages:
        plan += [(mat4, N, c) for N in (128, 256) for c in cases]
    if "B" in args.stages:
        plan += [(mat6, N, c) for N in (192, 384) for c in pick(cases, STRONG)]
    if "C" in args.stages:
        plan += [(mat4, 512, c) for c in pick(cases, FINEST)]

    done = {(r["case_key"], r["E_max"], r["N"]) for r in state["runs"]}
    for mat, N, case in plan:
        ck = key(case)
        if (ck, mat.E_max, N) in done:
            continue
        tr, wall = one_run(mat, N, case)
        rid = f"{ck}_E{mat.E_max:.0f}_N{N}"
        for o in OBS:
            series[f"{rid}_{o}"] = tr[o]
        series[f"{rid}_t"] = tr.t
        state["runs"].append(dict(id=rid, case_key=ck, case=case, E_max=mat.E_max, N=N, wall=wall,
                                  Te_max=float(tr["T_e"].max()), n_steps=tr.meta["n_steps"], h_min=tr.meta["h_min"],
                                  f_top_max=float(tr.f[:, -1].max()),
                                  audit=float(np.abs(tr["audit"] - tr["audit"][0]).max() / max(tr["u_abs"][-1], 1e-30))))
        state["summary"] = summarise(series, state["runs"])
        OUT.write_text(json.dumps(state, indent=1))
        np.savez_compressed(SERIES, **series)
        print(f"{rid:<34} T_e,max {tr['T_e'].max():7.0f} K  steps {tr.meta['n_steps']:5d}  f_top {tr.f[:, -1].max():.1e}  {wall:6.0f} s", flush=True)
    print("done")


if __name__ == "__main__":
    main()
