#!/usr/bin/env python3
"""Validation of the density-matrix particle model (enzkin.densmat) against refs. [16] and [62].

The model is that of [16, 62] with the dipole strength reconstructed from first principles (Eq. 16 of [16]
fails the oscillator-strength sum rule as printed; docs/theory/tier2_density_matrix.md).  One number is
fitted: eps_inf, to Re eps(175 THz) = -1.6 of [16] Fig. 2b.  Everything else is compared, not tuned.

Targets (from the digests in references/LITERATURE_SUMMARY.md):
  [16] Fig. 1b   |D|^2 at E' = E_F, L = 4 nm (by eye from the log plot)
  [16] Fig. 2b   eps(175 THz) = -1.6 + 1.4i, Re eps = 0 at 215 THz           (vector paths of the PDF)
  [16] Fig. 6    E0 = 1.2 V/nm: absorbed power between -0.14 and +0.33 eV/nm^3/fs, U settles at 0.63 eV/nm^3
  [16] Fig. 8c   local-field spectral peak 179, 178, 177, 176, 173, 170, 168 THz for E0 = 0 ... 2.4 V/nm
  [62] Fig. 4    E0 = 1.2 V/nm: eps(160 THz) from -2.6 + 1.85i to about -1.0 + 2.25i
  [62] Sec. 4.2  T_e about 4500 K and 1/eta_ee about 12 fs at E0 = 1.8 V/nm

Parameter set: ito_sarkar2023 ([16] states the parameters of [14, 15]); --un2026 uses Table S1 of [62].

Output: results/tier2/densmat_validation.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from enzkin import constants as c
from enzkin.densmat import ParticleDM, box_dipole_weight, gaussian_pulse
from enzkin.material import ito_sarkar2023, ito_un2026

OUT = Path("results/tier2/densmat_validation.json")
W = c.frequency_thz_to_omega
PEAKS_16 = {0.0: 179, 0.4: 178, 0.6: 177, 0.8: 176, 1.2: 173, 1.8: 170, 2.4: 168}


def weight_as_printed(k, kp, L):
    """V |D|^2 / e^2 from Eq. (16) of [16] exactly as printed (nm^5)."""
    x, y = k * L, kp * L
    d2, dm = x * x - y * y, (x - y) ** 2
    curly = ((d2**2 - 16) / 4 * np.log((4 + np.abs(d2)) / (4 + dm)) - (d2**2 + 1) * np.log((1 + np.abs(d2)) / (1 + dm))
             - 0.75 * dm**2 * np.log(np.abs((x - y) / (x + y))))
    return 8 * np.pi**2 * L**5 * curly / (x * y * d2**2)


def sum_rule_at(grid, band, weight, E):
    """int rho(E') (E' - E) W(E, E') dE' over hbar^2 w_opt(E) / m0 (1 if the oscillator strength is complete)."""
    i = int(np.argmin(np.abs(grid.E - E)))
    with np.errstate(divide="ignore", invalid="ignore"):
        w = weight(grid.k[i], grid.k)
    w[i] = 0.0
    return float(np.sum(w * grid.rho * (grid.E - grid.E[i])) * grid.dE / (2 * band.h2m * grid.w_opt[i]))


def spectrum_peak(t, y):
    """Peak frequency (THz) and FWHM (THz) of |Y(omega)|, from a zero-padded FFT with parabolic refinement."""
    dt = t[1] - t[0]
    n = 1 << 16
    Y = np.abs(np.fft.rfft(y, n))
    fr = np.fft.rfftfreq(n, dt) * 1e3
    sel = (fr > 60) & (fr < 400)
    i = np.argmax(np.where(sel, Y, 0.0))
    half = Y[i] / 2
    lo = i - np.argmax(Y[i::-1] < half)
    hi = i + np.argmax(Y[i:] < half)
    return float(fr[i]), float(fr[hi] - fr[lo])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--un2026", action="store_true")
    ap.add_argument("--N", type=int, default=192)
    ap.add_argument("--E-max", type=float, default=6.0)
    args = ap.parse_args()
    mat = (ito_un2026 if args.un2026 else ito_sarkar2023)(E_max=args.E_max)
    dm = ParticleDM(mat, N=args.N).calibrate_eps_inf(-1.6, 175.0)
    band = mat.band
    kF = band.k_fermi(mat.n_e)
    EF = band.fermi_energy(mat.n_e)
    res = dict(parameter_set=mat.name, N=args.N, E_max=args.E_max, L=dm.L, c_reg=dm.c_reg, eps_inf=dm.mat.eps_inf,
               sum_rule=dm.sum_rule())

    # dipole element against [16] Fig. 1b
    V, e2 = 4 / 3 * np.pi * 2.0**3, 1.602176634e-19**2
    res["dipole_log10_C2m2"] = {f"{x} E_F": float(np.log10(box_dipole_weight(band.k(x * EF), kF, dm.L) * 1e-45 / (V * 1e-27) * e2))
                                for x in (0.5, 1.5, 2.0)}
    res["dipole_paper_by_eye"] = {"0.5 E_F": -58.6, "1.5 E_F": -58.3, "2.0 E_F": -60.0}
    res["sum_rule_at_EF"] = dict(reconstructed=sum_rule_at(dm.grid, band, lambda k, kp: box_dipole_weight(k, kp, dm.L), EF),
                                 eq16_as_printed=sum_rule_at(dm.grid, band, lambda k, kp: weight_as_printed(k, kp, dm.L), EF))
    print(f"oscillator-strength sum rule at E_F: reconstructed {res['sum_rule_at_EF']['reconstructed']:.2f}, "
          f"Eq. (16) as printed {res['sum_rule_at_EF']['eq16_as_printed']:.3g}")

    # linear permittivity
    fs = np.linspace(150.0, 320.0, 341)
    eps = dm.permittivity(W(fs))
    res["eps_175"] = [float(dm.permittivity(W(175.0)).real), float(dm.permittivity(W(175.0)).imag)]
    res["eps_160"] = [float(dm.permittivity(W(160.0)).real), float(dm.permittivity(W(160.0)).imag)]
    res["zero_crossing_THz"] = float(fs[np.argmin(np.abs(eps.real))])
    res["eps_curve"] = dict(f_thz=fs[::10].tolist(), re=eps.real[::10].tolist(), im=eps.imag[::10].tolist())

    # nonlinear runs
    rows = []
    for E0 in (1e-3, 0.4, 0.6, 0.8, 1.2, 1.8, 2.4):
        r = dm.run(gaussian_pulse(E0, 5.0, 175.0), -20.0, 40.0, probe_thz=160.0)
        pk, width = spectrum_peak(r.t, r.E_loc)
        power = r.E_loc * r.J / c.E2_OVER_EPS0
        i5 = int(np.argmin(np.abs(r.t - 5.0)))
        row = dict(E0=E0, peak_THz=pk, fwhm_THz=width, U_field=r.meta["W_field"], power_min=float(power.min()),
                   power_max=float(power.max()), T_e_max=float(r.T_e.max()), T_e_5fs=float(r.T_e[i5]),
                   eps160_start=[float(r.meta["eps"][0].real), float(r.meta["eps"][0].imag)],
                   eps160_end=[float(r.meta["eps"][-1].real), float(r.meta["eps"][-1].imag)],
                   f_top_max=r.meta["f_top_max"], wall=r.meta["wall"])
        # e-e rate at E_F, 5 fs after the peak, from the populations at that time
        rk = dm.run(gaussian_pulse(E0, 5.0, 175.0), -20.0, 5.0, keep_f=True)
        rate = dm.km.ee.rate(np.clip(rk.f[-1], 0.0, 1.0))
        row["eta_ee_EF_5fs"] = float(np.interp(EF, dm.grid.E, rate))
        row["eta_ee_EF+1eV_5fs"] = float(np.interp(EF + 1.0, dm.grid.E, rate))
        rows.append(row)
        print(f"E0 {E0:5.3f} V/nm: peak {pk:6.1f} THz (paper {PEAKS_16.get(round(E0, 1) if E0 > 0.01 else 0.0)}), FWHM {width:5.1f}, "
              f"U {row['U_field']:.3f} eV/nm^3, power {row['power_min']:+.3f}..{row['power_max']:+.3f}, T_e(5 fs) {row['T_e_5fs']:.0f} K, "
              f"eps160 {row['eps160_start'][0]:+.2f}{row['eps160_start'][1]:+.2f}i -> {row['eps160_end'][0]:+.2f}{row['eps160_end'][1]:+.2f}i, "
              f"1/eta_ee(E_F, 5 fs) {1 / row['eta_ee_EF_5fs']:.0f} fs, f_top {row['f_top_max']:.1e}, {row['wall']:.1f} s", flush=True)
    res["runs"] = rows
    res["paper"] = dict(eps_175=[-1.6, 1.4], zero_crossing_THz=215, U_1p2=0.63, power_1p2=[-0.14, 0.33],
                        peaks_THz=PEAKS_16, eps160_1p2=[[-2.6, 1.85], [-1.0, 2.25]], T_e_1p8=4500, tau_ee_1p8_fs=12)
    out = OUT if not args.un2026 else OUT.with_name("densmat_validation_un2026.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    print(f"eps_inf {res['eps_inf']:.3f}; sum rule {res['sum_rule']:.3f}; eps(175) {res['eps_175'][0]:+.2f}{res['eps_175'][1]:+.2f}i (paper -1.6+1.4i); "
          f"zero crossing {res['zero_crossing_THz']:.0f} THz (paper 215); log10|D|^2 {res['dipole_log10_C2m2']}")
    print("wrote", out)


if __name__ == "__main__":
    main()
