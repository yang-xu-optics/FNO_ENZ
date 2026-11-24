#!/usr/bin/env python3
"""Decision test for a learned few-cycle constitutive law (plan, deferred operator; working log session 6).

Question.  Under few-cycle pumping, does the coherent density-matrix model of [16, 62] contain anything in
its field-resolved response that a cheap physics closure misses?  If not, a learned operator trained on it
has nothing to learn.

Truth:     L0  density matrix: populations with collisions + all pair oscillators        (enzkin.densmat)
Closures:  L1  Fermi-Dirac populations at T_e(t) + all pair oscillators                  ("thermal model" of [62])
           L2  T_e(t) + 17 or 9 fixed poles, strengths fitted to chi(omega; T_e), widths rising with the
               collision rate                                                             (few-pole law)
           L2f the same 17 poles with temperature-independent widths                      (adiabatic fit)
           L3  T_e(t) + one Drude pole, bulk plasma frequency, fitted size damping       (plan Tier 2a)
All closures are driven by the instantaneous work E_loc dP/dt.

Setting: the 4 nm ITO sphere of [16, 62] (uniform quasi-static field), Gaussian pulses at 175 THz (and two
at 230 THz), scanned over peak field, duration and carrier-envelope phase.

Observable: the sheet current J = dP/dt, which is what a field-resolved transmission measurement of a thin
layer of particles sees.  For each model the response to a weak pulse, scaled, is subtracted to obtain the
nonlinear part dJ.  Reported per condition:

    nl          max|dJ_0| / max|J_0|                 how nonlinear the truth is
    err_tot     max|J_k - J_0| / max|J_0|            what a measurement of the total field would see
    err_nl      max|dJ_k - dJ_0| / max|J_0|          the same with each model's linear response removed
    U ratio     absorbed energy, model / truth
    h3 ratio    third-harmonic amplitude of J, model / truth

Output: results/tier2/closure_test.json, results/tier2/closure_test_series.npz
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from enzkin import constants as c
from enzkin.densmat import ParticleDM, ReducedResponse, gaussian_pulse
from enzkin.material import ito_sarkar2023, ito_un2026

OUT = Path("results/tier2/closure_test.json")
SERIES = Path("results/tier2/closure_test_series.npz")
E_WEAK = 1e-3
KEEP = {("175", 5.0, 1.2, 0.0), ("175", 5.0, 2.4, 0.0), ("175", 20.0, 1.2, 0.0), ("175", 5.0, 2.4, 1.57)}


def conditions(quick=False):
    out = []
    grid = {5.0: (0.3, 0.6, 1.2, 1.8, 2.4), 10.0: (0.3, 0.6, 1.2, 1.8), 20.0: (0.3, 0.6, 1.2)}
    if quick:
        grid = {5.0: (0.6, 2.4), 20.0: (1.2,)}
    for fwhm, fields in grid.items():
        for E0 in fields:
            for cep in ((0.0,) if quick else (0.0, 1.57)):
                out.append(dict(f_thz=175.0, fwhm=fwhm, E0=E0, cep=cep))
    if not quick:
        out += [dict(f_thz=230.0, fwhm=5.0, E0=E0, cep=0.0) for E0 in (1.2, 2.4)]
    return out


def window(fwhm):
    half = max(2.5 * fwhm, 12.0)
    return -half, half + 25.0


def harmonic(t, y, f_thz, n):
    """Spectral amplitude of y in a band around n times the carrier."""
    dt = t[1] - t[0]
    Y = np.abs(np.fft.rfft(y * np.hanning(y.size), 1 << 15))
    fr = np.fft.rfftfreq(1 << 15, dt) * 1e3
    sel = (fr > (n - 0.5) * f_thz) & (fr < (n + 0.5) * f_thz)
    return float(Y[sel].max())


def model_uncertainty(mat, N):
    """How much the density matrix itself moves under choices the papers do not fix.

    Two variations of the truth, each recalibrated to Re eps(175 THz) = -1.6 as the reference is:
      c_reg 0.9   regularisation of the dipole element at the level spacing, 10 % smaller
      Umklapp 1/2 only half of the e-e collision rate dephases the polarisation ([16] footnote: the Umklapp share
                  "is likely to be of order unity" but unknown)
    Reported like the closure errors: largest deviation of J from the reference, in units of its peak.
    """
    ref = ParticleDM(mat, N=N).calibrate_eps_inf(-1.6, 175.0)
    alt_c = ParticleDM(mat, N=N, c_reg=0.9).calibrate_eps_inf(-1.6, 175.0)
    alt_u = ParticleDM(mat, N=N)
    ee_rate = alt_u.km.ee.rate
    damping_full = alt_u.km.damping
    alt_u.km.damping = lambda f, T_ph, rate_ee=None: damping_full(f, T_ph, 0.5 * (ee_rate(f) if rate_ee is None else rate_ee))
    alt_u._tables(alt_u.T_tab.size, float(alt_u.T_tab[-1]))  # rebuild the damping tables with the halved e-e share
    alt_u.calibrate_eps_inf(-1.6, 175.0)
    rows = []
    for fwhm, E0 in ((5.0, 0.6), (5.0, 1.2), (5.0, 2.4), (20.0, 1.2)):
        t_lo, t_hi = window(fwhm)
        p = gaussian_pulse(E0, fwhm, 175.0)
        r0, rc, ru = ref.run(p, t_lo, t_hi), alt_c.run(p, t_lo, t_hi), alt_u.run(p, t_lo, t_hi)
        jmax = float(np.abs(r0.J).max())
        rows.append(dict(fwhm=fwhm, E0=E0, c_reg=float(np.abs(rc.J - r0.J).max() / jmax), umklapp_half=float(np.abs(ru.J - r0.J).max() / jmax)))
        print(f"model uncertainty {fwhm:.0f} fs {E0} V/nm: c_reg 0.9 -> {100 * rows[-1]['c_reg']:.1f} %, Umklapp 1/2 -> {100 * rows[-1]['umklapp_half']:.1f} %", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--un2026", action="store_true", help="parameter set of [62] Table S1 instead of [14, 15]")
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--E-max", type=float, default=8.0)
    ap.add_argument("--uncertainty-only", action="store_true", help="add the model-uncertainty block to an existing result file")
    args = ap.parse_args()
    if args.uncertainty_only:
        out = OUT if not args.un2026 else OUT.with_name("closure_test_un2026.json")
        res = json.loads(out.read_text())
        res["model_uncertainty"] = model_uncertainty((ito_un2026 if args.un2026 else ito_sarkar2023)(E_max=args.E_max), args.N)
        out.write_text(json.dumps(res, indent=1))
        print("updated", out)
        return

    mat = (ito_un2026 if args.un2026 else ito_sarkar2023)(E_max=args.E_max)
    t0 = time.perf_counter()
    dm = ParticleDM(mat, N=args.N).calibrate_eps_inf(-1.6, 175.0)
    closures = {"L2 17 poles": ReducedResponse.fitted(dm, 16), "L2 9 poles": ReducedResponse.fitted(dm, 8),
                "L2f 17 poles, fixed widths": ReducedResponse.fitted(dm, 16, collisional=False),
                "L3 Drude": ReducedResponse.drude(dm, 175.0)}
    build = time.perf_counter() - t0
    w175 = c.frequency_thz_to_omega(175.0)
    lin = {"L0 density matrix": dm.permittivity(w175)}
    lin.update({k: v.permittivity(w175) for k, v in closures.items()})
    print(f"built in {build:.1f} s; eps_inf {dm.mat.eps_inf:.3f}; eps(175 THz) at 300 K: " +
          "; ".join(f"{k} {v.real:+.2f}{v.imag:+.2f}i" for k, v in lin.items()), flush=True)

    def run_all(pulse, t_lo, t_hi):
        runs = {"L0 density matrix": dm.run(pulse, t_lo, t_hi), "L1 thermal": dm.run(pulse, t_lo, t_hi, populations="thermal")}
        runs.update({k: v.run(pulse, t_lo, t_hi) for k, v in closures.items()})
        return runs

    weak, rows, series = {}, [], {}
    for cond in conditions(args.quick):
        f_thz, fwhm, E0, cep = cond["f_thz"], cond["fwhm"], cond["E0"], cond["cep"]
        t_lo, t_hi = window(fwhm)
        wk = (f_thz, fwhm, cep)
        if wk not in weak:
            weak[wk] = run_all(gaussian_pulse(E_WEAK, fwhm, f_thz, cep), t_lo, t_hi)
        runs = run_all(gaussian_pulse(E0, fwhm, f_thz, cep), t_lo, t_hi)
        ref = runs["L0 density matrix"]
        scale = E0 / E_WEAK
        dJ0 = ref.J - scale * weak[wk]["L0 density matrix"].J
        jmax = float(np.abs(ref.J).max())
        h3_0 = harmonic(ref.t, ref.J, f_thz, 3)
        row = dict(cond, nl=float(np.abs(dJ0).max() / jmax), U=ref.meta["W_field"], T_e_max=float(ref.T_e.max()),
                   f_top_max=ref.meta["f_top_max"], h3_over_h1=h3_0 / harmonic(ref.t, ref.J, f_thz, 1), models={})
        for name, r in runs.items():
            if name == "L0 density matrix":
                row["wall_L0"] = r.meta["wall"]
                continue
            dJ = r.J - scale * weak[wk][name].J
            U = r.meta.get("W_field", float(np.trapezoid(r.E_loc * r.J, r.t) / c.E2_OVER_EPS0))
            row["models"][name] = dict(err_tot=float(np.abs(r.J - ref.J).max() / jmax), err_nl=float(np.abs(dJ - dJ0).max() / jmax),
                                       err_nl_rel=float(np.abs(dJ - dJ0).max() / max(np.abs(dJ0).max(), 1e-300)),
                                       U_ratio=float(U / ref.meta["W_field"]), T_e_max=float(r.T_e.max()),
                                       h3_ratio=harmonic(r.t, r.J, f_thz, 3) / h3_0, wall=r.meta["wall"])
        rows.append(row)
        key = (f"{f_thz:.0f}", fwhm, E0, cep)
        if key in KEEP:
            tag = f"f{f_thz:.0f}_tau{fwhm:.0f}_E{E0:.1f}_cep{cep:.2f}"
            series[f"{tag}_t"] = ref.t
            series[f"{tag}_Einc"] = ref.E_inc
            for name, r in runs.items():
                series[f"{tag}_J_{name}"] = r.J
                series[f"{tag}_dJ_{name}"] = r.J - scale * weak[wk][name].J
                series[f"{tag}_Te_{name}"] = r.T_e
        m = row["models"]
        print(f"{f_thz:.0f} THz {fwhm:4.0f} fs {E0:3.1f} V/nm cep {cep:.2f}: nl {100 * row['nl']:5.1f} %  T_e {row['T_e_max']:6.0f} K  " +
              "  ".join(f"{k.split()[0]}{'-' + k.split()[1] if k.startswith('L2') else ''}: {100 * v['err_tot']:4.1f}/{100 * v['err_nl']:4.1f}/{v['U_ratio']:.2f}"
                        for k, v in m.items()) + f"  [{row['wall_L0']:.0f} s]", flush=True)

    # cost per run of the truth and the closures on a common window (60 fs)
    cost = {}
    p = gaussian_pulse(1.2, 5.0, 175.0)
    for name, fn in (("L0 density matrix", lambda: dm.run(p, -20.0, 40.0)), ("L1 thermal", lambda: dm.run(p, -20.0, 40.0, populations="thermal")),
                     *((k, (lambda v=v: v.run(p, -20.0, 40.0))) for k, v in closures.items())):
        r = fn()
        cost[name] = dict(wall_s=r.meta["wall"], oscillators=r.meta["n_osc"], collision_calls=r.meta.get("n_coll", 0))
    res = dict(parameter_set=mat.name, N=args.N, E_max=args.E_max, eps_inf=dm.mat.eps_inf, sum_rule=dm.sum_rule(),
               eps_175={k: [v.real, v.imag] for k, v in lin.items()}, build_s=build, cost_60fs=cost, rows=rows,
               fit_residual={k: getattr(v, "fit_residual", None) for k, v in closures.items()})
    out = OUT if not args.un2026 else OUT.with_name("closure_test_un2026.json")
    if args.quick:
        out = out.with_name(out.stem + "_quick.json")
    if not args.quick:
        res["model_uncertainty"] = model_uncertainty(mat, args.N)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    if not args.quick and not args.un2026:
        np.savez_compressed(SERIES, **series)
    print("cost for a 60 fs window:", {k: f"{v['wall_s']:.2f} s" for k, v in cost.items()})
    print("wrote", out)


if __name__ == "__main__":
    main()
