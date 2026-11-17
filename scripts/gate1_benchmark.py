#!/usr/bin/env python3
"""Gate 1 (project plan, month 2): is a learned surrogate of the collision integral needed at all?

Pass criterion in the plan:
    "B1 costs at least 100x more per evaluation than B0, and B0 misses B1's permittivity transient
     by more than 2 %."

This script measures both numbers.

Part A - cost.   Wall time per evaluation of B1 (screened Coulomb, O(N^3)) and B0 (random-k FFT, O(N log N))
                 at several grid sizes, on this machine.
Part B - accuracy.  The same optical excitation is run with B1, B0, a conserving relaxation-time model and
                 the thermal (two-temperature) model.  Two drive modes:
                   energy  identical absorbed power p_abs(t) for all models (isolates the kinetics);
                   field   identical in-film field, absorption through each model's own eps''(t).
                 "Miss" of model X against B1 for an observable y(t):
                     miss = max_t |y_X(t) - y_B1(t)| / max_t |y_B1(t) - y_B1(t_0)|
                 i.e. the largest deviation in units of the size of the transient itself.
Part C - grid check.  One B1 case repeated at twice the resolution.

Usage:  python scripts/gate1_benchmark.py [--N 128] [--quick] [--out results/wp1/gate1.json]
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import time
from pathlib import Path

import numba
import numpy as np

from enzkin import thermo
from enzkin.collisions.ee import ElectronElectronB0, ElectronElectronB1, ElectronElectronRTA
from enzkin.drive import FieldDrive, PowerDrive, Pulse
from enzkin.grid import EnergyGrid
from enzkin.integrate import run
from enzkin.material import ito_sarkar2023
from enzkin.model import KineticModel
from enzkin.ttm import ThermalModel

OBS = ("eps_re", "eps_im", "wp2", "gamma", "T_e")


def cost_table(mat, sizes, reps=7):
    rows = []
    for N in sizes:
        g = EnergyGrid(mat.band, mat.E_max, N)
        f0 = thermo.fd_of_density(g, mat.n_e, 300.0)
        sh = int(round(0.95 / g.dE))
        f = f0.copy()
        bump = 0.08 * f0[:-sh] * (1 - f0[sh:])
        f[sh:] += bump
        f[:-sh] -= bump * g.rho[sh:] / g.rho[:-sh]
        b1, b0 = ElectronElectronB1(g, mat), ElectronElectronB0(g, mat)
        b1.rhs(f0), b0.rhs(f0)  # JIT / plan warm-up
        t1 = []
        for i in range(reps):
            fp = f * (1 + 1e-9 * (i + 1))  # new screening wavevector -> lookup tables rebuilt, as in a real run
            t = time.perf_counter()
            b1.rhs(fp)
            t1.append(time.perf_counter() - t)
        n0 = 300
        t = time.perf_counter()
        for _ in range(n0):
            b0.rhs(f)
        t0 = (time.perf_counter() - t) / n0
        rows.append(dict(N=N, b1_ms=1e3 * float(np.median(t1)), b0_ms=1e3 * t0, ratio=float(np.median(t1) / t0)))
        print(f"  N={N:4d}:  B1 {rows[-1]['b1_ms']:8.2f} ms   B0 {rows[-1]['b0_ms']:6.3f} ms   ratio {rows[-1]['ratio']:6.0f}", flush=True)
    return rows


def build(mat, N, kind, drive):
    if kind == "thermal":
        return ThermalModel(mat, N=N, ee=ElectronElectronB1, drive=drive)
    ee = dict(B1=ElectronElectronB1, B0=ElectronElectronB0, RTA=ElectronElectronRTA)[kind]
    return KineticModel(mat, N=N, ee=ee, drive=drive)


def transient(mat, N, kind, case):
    pulse = Pulse(case["wavelength"], case["fwhm"])
    if case["mode"] == "energy":
        drive = PowerDrive.from_absorbed_energy(pulse, case["level"])
    else:
        drive = FieldDrive(pulse, case["level"])
    model = build(mat, N, kind, drive)
    fwhm = case["fwhm"]
    t_lo, t_hi = -2.5 * fwhm, max(2.5 * fwhm, 300.0) + 300.0
    dt = min(0.5, fwhm / 60.0)
    n_samp = max(int(round(2.0 / dt)), 1)  # sample every ~2 fs
    t0 = time.perf_counter()
    tr = run(model, t_lo, t_hi, dt, sample_every=n_samp, omega_probe=pulse.omega)
    wall = time.perf_counter() - t0
    return tr, wall


def miss(a, ref):
    amp = np.max(np.abs(ref - ref[0]))
    return float(np.max(np.abs(a - ref)) / amp) if amp > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--quick", action="store_true", help="fewer cases (smoke test)")
    ap.add_argument("--out", default="results/wp1/gate1.json")
    ap.add_argument("--series", default="results/wp1/gate1_series.npz")
    args = ap.parse_args()

    # linear optics calibrated to [15]: eps(230 THz, 300 K) = 0.172 + 0.399i  (scripts/validate_un2023.py)
    mat = ito_sarkar2023(eps_inf=3.4509, gamma_imp=0.0631)
    print("Part A - cost per collision-integral evaluation")
    cost = cost_table(mat, [128, 256] if args.quick else [64, 128, 256, 512])

    # absorbed energy densities: 6e7 J/m^3 ~ T_e 3000 K, 3e8 ~ 6500 K, 1e9 ~ 12000 K
    cases = []
    for fwhm in ([30.0, 220.0] if args.quick else [20.0, 50.0, 100.0, 220.0]):
        for u in ([3e8] if args.quick else [6e7, 3e8, 1e9]):
            cases.append(dict(mode="energy", wavelength=1300.0, fwhm=fwhm, level=u))
    for fwhm in ([30.0] if args.quick else [20.0, 220.0]):
        for inten in ([1e11] if args.quick else [2e10, 2e11]):
            cases.append(dict(mode="field", wavelength=1300.0, fwhm=fwhm, level=inten))

    print(f"Part B - transients at N = {args.N}")
    rows, series = [], {}
    for ic, case in enumerate(cases):
        res = {}
        for kind in ("B1", "B0", "RTA", "thermal"):
            res[kind], wall = transient(mat, args.N, kind, case)
            res[kind].meta["wall"] = wall
        ref = res["B1"]
        row = dict(case=case, Te_max_B1=float(ref["T_e"].max()), nonthermal_max_B1=float(ref["nonthermal"].max()),
                   eps_re_amp_B1=float(np.max(np.abs(ref["eps_re"] - ref["eps_re"][0]))),
                   wall={k: v.meta["wall"] for k, v in res.items()}, steps={k: v.meta["n_steps"] for k, v in res.items()},
                   miss={k: {o: miss(res[k][o], ref[o]) for o in OBS} for k in ("B0", "RTA", "thermal")})
        rows.append(row)
        for k, v in res.items():
            for o in ("t",) + OBS + ("nonthermal",):
                series[f"c{ic}_{k}_{o}"] = v.t if o == "t" else v[o]
        m = row["miss"]
        print(f"  {case['mode']:6} fwhm {case['fwhm']:5.0f} fs  level {case['level']:.1e}  T_e,max {row['Te_max_B1']:6.0f} K | "
              f"miss eps' : B0 {100 * m['B0']['eps_re']:5.1f} %  RTA {100 * m['RTA']['eps_re']:5.1f} %  thermal {100 * m['thermal']['eps_re']:5.1f} % | "
              f"miss wp2: B0 {100 * m['B0']['wp2']:5.1f} %  thermal {100 * m['thermal']['wp2']:5.1f} % | "
              f"wall B1 {row['wall']['B1']:.0f} s  B0 {row['wall']['B0']:.0f} s", flush=True)

    print("Part C - grid check (B1, 50 fs, 3e8 J/m^3)")
    chk = dict(mode="energy", wavelength=1300.0, fwhm=50.0, level=3e8)
    a, _ = transient(mat, args.N, "B1", chk)
    b, wall_fine = transient(mat, 2 * args.N, "B1", chk)
    grid = dict(N=args.N, N_fine=2 * args.N, wall_fine=wall_fine,
                miss={o: miss(np.interp(b.t, a.t, a[o]), b[o]) for o in OBS})
    print("  ", {k: f"{100 * v:.2f} %" for k, v in grid["miss"].items()})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(
        machine=dict(platform=platform.platform(), cpu_count=os.cpu_count(), numba_threads=numba.get_num_threads(),
                     numba=numba.__version__, numpy=np.__version__),
        material=dict(eps_inf=mat.eps_inf, gamma_imp=mat.gamma_imp, c_ee=mat.c_ee, n_e=mat.n_e),
        N=args.N, cost=cost, cases=rows, grid_check=grid), indent=1))
    np.savez_compressed(args.series, **series)
    print("wrote", out, "and", args.series)


if __name__ == "__main__":
    main()
