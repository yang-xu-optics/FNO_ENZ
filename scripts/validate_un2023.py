#!/usr/bin/env python3
"""WP1 validation against Un, Sarkar & Sivan, Phys. Rev. Applied 19, 044043 (2023)  [ref. 15].

The paper models a 40 nm ITO film on glass under 100 nm of gold, p-polarised 230 THz pump at 60 degrees,
220 fs (2.5-75 GW/cm^2) and 30 fs (161 and 22 GW/cm^2).  Its figures give the targets tabulated below
(values extracted from the vector graphics of the arXiv source; see references/LITERATURE_SUMMARY.md [15]).

The paper does not print eps_inf, the electron-impurity rate, the gold permittivity or the incidence
medium, so an exact reproduction is not possible from the paper alone.  The stack assumed here is
air -> glass substrate -> 40 nm ITO -> 100 nm Au, with the 60 degrees measured in air (35.3 degrees inside
the glass): of the geometries tried, it is the one whose absorptance falls monotonically as Re eps rises
from its 300 K value, as the paper's absorptivity does.  Three comparisons are run, from least to most
dependent on the unknowns:

  energy   absorbed energy imposed from the paper's own absorptivity (Gaussian in time).  Tests the
           thermodynamics and kinetics only: C_e(T_e), mu(T_e), e-ph coupling, thermalisation.
  film     fully self-consistent thin-film drive (transfer matrix + eps[f]), linear optics calibrated to the
           paper's 300 K permittivity at 230 THz, screening recomputed from f at every step.
  film-static  the same with the Thomas-Fermi wavevector frozen at its 300 K value.

Usage:  python scripts/validate_un2023.py [--N 128] [--model kinetic|thermal] [--out results/wp1/validation_un2023.json]
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from enzkin import constants as c
from enzkin.calibrate import calibrate_linear_optics
from enzkin.collisions.ee import ElectronElectronB1
from enzkin.drive import PowerDrive, Pulse, ThinFilmDrive
from enzkin.integrate import run
from enzkin.material import ito_sarkar2023
from enzkin.model import KineticModel
from enzkin.ttm import ThermalModel

# --- targets from [15] (BE/"ANTHEM" columns) ----------------------------------------------------------
#        I0 [GW/cm2], fwhm [fs], T_e,max [K], t(T_e,max) [fs], T_ph(990 fs) [K], eps'_max, eps''_max, absorptivity
TARGETS = [
    (2.5, 220.0, 2464, 124, 319.2, None, None, 0.777),
    (5.0, 220.0, 3313, 117, 332.5, 0.765, 0.495, 0.668),
    (10.0, 220.0, 4348, 109, 352.7, None, None, 0.545),
    (22.0, 220.0, 5902, 105, 389.6, 1.476, 0.547, 0.418),
    (50.0, 220.0, 8393, 108, 456.3, 2.013, 0.594, 0.318),
    (75.0, 220.0, 10262, 115, 507.2, 2.315, 0.630, 0.281),
    (161.0, 30.0, 7470, 26, None, 1.717, 0.546, None),
    (22.0, 30.0, 3187, 27, None, 0.682, 0.482, None),
]
EPS_300K = 0.172 + 0.399j  # [15] Fig. 1 inset at 230 THz
F_PUMP_THZ = 230.0
THICKNESS, ANGLE = 40.0, 60.0  # angle of incidence in air
EPS_GLASS, EPS_AU = 2.25, -80.0 + 8.0j  # assumed (not given in [15]); gold ~ Johnson-Christy at 1.3 um
ANGLE_GLASS = math.degrees(math.asin(math.sin(math.radians(ANGLE)) / math.sqrt(EPS_GLASS)))  # 35.3 degrees
# intensity inside the glass for an intensity I0 in air (flux through the film plane is conserved; the
# p-polarised Fresnel loss at the air-glass face, 0.2 % at 60 degrees, is neglected)
GLASS_OVER_AIR = math.cos(math.radians(ANGLE)) / math.cos(math.radians(ANGLE_GLASS))


def calibrated_material(N: int):
    """eps_inf and gamma_imp such that eps(230 THz, 300 K) equals the paper's plotted value."""
    return calibrate_linear_optics(ito_sarkar2023(), N, F_PUMP_THZ, EPS_300K)


def incident_fluence_on_film(i0_gw: float, fwhm: float) -> float:
    """Incident energy per unit film area, J/m^2 (projection by cos(theta))."""
    return i0_gw * 1e13 * fwhm * 1e-15 * math.sqrt(math.pi / (4 * math.log(2))) * math.cos(math.radians(ANGLE))


def one_run(mat, N, model_cls, mode, i0, fwhm, absorptivity):
    lam = c.C_LIGHT / (F_PUMP_THZ * 1e-3)
    pulse = Pulse(lam, fwhm)
    if mode == "energy":
        u_abs = absorptivity * incident_fluence_on_film(i0, fwhm) / (THICKNESS * 1e-9)  # J/m^3
        drive = PowerDrive.from_absorbed_energy(pulse, u_abs)
    else:
        drive = ThinFilmDrive(pulse, i0 * 1e9 * GLASS_OVER_AIR, thickness=THICKNESS, angle_deg=ANGLE_GLASS,
                              polarization="p", eps_in=EPS_GLASS, eps_out=1.0, layers_after=[(EPS_AU, 100.0)])
    model = model_cls(mat, N=N, drive=drive)
    if mode == "film-static":
        model.ee.dynamic_screening_from_f = False
    t_lo, t_hi = (-450.0, 1000.0) if fwhm > 100 else (-80.0, 460.0)
    dt = 0.5 if fwhm > 100 else 0.2
    tr = run(model, t_lo, t_hi, dt, sample_every=10 if fwhm > 100 else 5)
    i = int(tr["T_e"].argmax())
    out = dict(mode=mode, model=model_cls.__name__, I0=i0, fwhm=fwhm,
               Te_max=float(tr["T_e"][i]), t_Te_max=float(tr.t[i]),
               Tph_990=float(np.interp(990.0, tr.t, tr.T_ph)) if t_hi >= 990 else None,
               eps_re_max=float(tr["eps_re"].max()), eps_im_max=float(tr["eps_im"].max()),
               eps_re_0=float(tr["eps_re"][0]), eps_im_0=float(tr["eps_im"][0]),
               absorptivity=float(tr["u_abs"][-1] / c.J_PER_M3_IN_EV_PER_NM3 * THICKNESS * 1e-9
                                  / incident_fluence_on_film(i0, fwhm)),
               audit=float(np.abs(tr["audit"] - tr["audit"][0]).max() / max(tr["u_abs"][-1], 1e-30)),
               n_steps=tr.meta["n_steps"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--model", choices=["kinetic", "thermal"], default="kinetic")
    ap.add_argument("--modes", nargs="*", default=["energy", "film", "film-static"])
    ap.add_argument("--out", default="results/wp1/validation_un2023.json")
    args = ap.parse_args()

    mat = calibrated_material(args.N)
    cls = KineticModel if args.model == "kinetic" else ThermalModel
    print(f"calibrated linear optics: eps_inf = {mat.eps_inf:.4f}, gamma_imp = {mat.gamma_imp:.4f} /fs")
    rows = []
    for (i0, fwhm, te, t_te, tph, e1, e2, a) in TARGETS:
        for mode in args.modes:
            if mode == "energy" and a is None:
                continue
            t0 = time.perf_counter()
            r = one_run(mat, args.N, cls, mode, i0, fwhm, a)
            r["target"] = dict(Te_max=te, t_Te_max=t_te, Tph_990=tph, eps_re_max=e1, eps_im_max=e2, absorptivity=a)
            rows.append(r)
            print(f"I0={i0:6.1f} GW/cm2 {fwhm:5.0f} fs  {mode:<11}  T_e,max {r['Te_max']:7.0f} K (paper {te})  "
                  f"A {r['absorptivity']:.3f} (paper {a})  eps' {r['eps_re_max']:.3f} (paper {e1})  "
                  f"eps'' {r['eps_im_max']:.3f} (paper {e2})  [{time.perf_counter() - t0:.0f} s]", flush=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(N=args.N, model=args.model, eps_inf=mat.eps_inf, gamma_imp=mat.gamma_imp,
                                   assumptions=dict(eps_glass=EPS_GLASS, eps_au=[EPS_AU.real, EPS_AU.imag],
                                                    thickness_nm=THICKNESS, angle_in_air_deg=ANGLE,
                                                    angle_in_glass_deg=ANGLE_GLASS),
                                   rows=rows), indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
