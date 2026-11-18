#!/usr/bin/env python3
"""Check the rates of Un, Sarkar & Sivan, Nanophotonics 15, e70272 (2026) against enzkin.

The paper (a 5 fs pulse on a 4 nm ITO sphere, density-matrix model) prints, in its Supporting
Information, the collision terms and a parameter set (Table S1).  Four things are checked here.

1. e-e rate.  Their exact energy-space rate (S19) has the same kernel and the same prefactor as the
   enzkin B1 integral with c_ee = 3/16.  For Fermi-Dirac distributions at several temperatures the B1
   rate is compared with their reading "about 12 fs, nearly energy independent, at T_e ~ 4500 K"
   (main text, Fig. 5a, Fig. S1a).
2. Generalised Fermi-liquid formula (S24)-(S26).  Implemented independently on a fine grid, with three
   corrections to the printed forms (see docs/theory/un2026_notes.md): the Fermi-Dirac product (S22) lacks
   a factor 1/4; the coefficients A_n carry a factor 1/(k_B T)^2 that (S26) omits; and K of (S25) is the
   T -> 0 coefficient of the B1 kernel, i.e. it carries the 1/(2 pi)^3 of (S19) that (S25) omits.  The A_n are compared with Fig. S1(b)-(d) and the formula
   with the B1 rate.
3. e-ph rate for their deformation potential D = 8.6 eV (Table S1) and for D = 17.2 eV of ref. [14], against
   Fig. 5c (about 0.01 1/fs at E_F, 0.06 1/fs at 3 E_F).
5. Consequence for a transient: 100 fs pulse at 1300 nm depositing 1e8 J/m^3, both parameter sets,
   1/e recovery time of Re eps and of T_e after their peaks.
4. Damping budget at 300 K: e-e + e-ph + impurity, with their Brooks-Herring impurity rate (S9) (sign in
   the bracket corrected, see notes) for n_imp = 0.05 nm^-3, Z = 1, against the enzkin calibration to [15]
   (D = 17.2 eV, constant gamma_imp = 0.0352 1/fs); resulting permittivity at 230 THz and the e-ph
   coupling G.

Values read from the paper's figures by eye are marked as such in the output; they carry +-10 % at best.

Output: results/wp1/check_un2026.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from enzkin import constants as c
from enzkin.calibrate import calibrate_linear_optics
from enzkin.collisions.ee import ElectronElectronB1, _antiderivative, b1_prefactor, thomas_fermi_kappa
from enzkin.collisions.eimp import brooks_herring_rate
from enzkin.collisions.eph import ElectronPhonon
from enzkin.grid import EnergyGrid
from enzkin.material import ito_sarkar2023
from enzkin import optics
from enzkin.drive import PowerDrive, Pulse
from enzkin.integrate import run
from enzkin.material import ito_un2026
from enzkin.model import KineticModel
from enzkin.thermo import fd_of_density, heat_capacity, mu_of_density

OUT = Path("results/wp1/check_un2026.json")
N = 256
KB = c.KB

# read by eye from the paper (main Fig. 5, SI Fig. S1)
PAPER = dict(
    ee_rate_4500K="~0.083 1/fs (12 fs), nearly flat in E (main text; Fig. 5a, S1a green, E0 = 1.8 V/nm, t = 5 fs)",
    A0_at_EF={300: 9.4, 1000: 8.6, 3000: 6.3, 5000: 4.9, 7000: 4.0},
    A2_at_EF={300: 0.95, 1000: 0.85, 3000: 0.55, 5000: 0.42, 7000: 0.35},
    A1_at_EF={300: -0.1, 1000: -0.3, 3000: -0.75, 5000: -0.85, 7000: -0.85},
    eph_rate={"0.1 E_F": 1e-3, "1 E_F": 0.01, "2 E_F": 0.03, "3 E_F": 0.06},
)


def generalised_fl(grid, mat, T, E_eval, N_int=600):
    """A_0, A_1, A_2 of (S26) (with 1/(k_B T)^2) and K of (S25) (with the B1 prefactor) at temperature T."""
    band = mat.band
    n = mat.n_e
    mu = mu_of_density(grid, n, T)
    kappa = thomas_fermi_kappa(grid, mat, fd_of_density(grid, n, T))
    kT = KB * T
    k_mu = float(band.k(max(mu, 1e-6)))
    phi_mu = float(_antiderivative(2 * k_mu, kappa))
    P = b1_prefactor(mat, kappa)
    K = P * (1 + 2 * mat.C * mu) ** 3 * phi_mu / (2 * k_mu)
    dE = mat.E_max / N_int
    e = (np.arange(N_int) + 0.5) * dE
    E1, E2 = np.meshgrid(e, e, indexing="ij")
    k1, k2 = band.k(E1), band.k(E2)
    m12 = band.mass_factor(E1) * band.mass_factor(E2)
    u = (E1 - E2) / kT
    # (S22) with the factor 1/4 that the printed form omits: the Fermi-Dirac product equals
    # cosh(x/2) / [4 cosh(x1/2) cosh(x2/2) cosh(x3/2)], x = (E - mu)/k_B T
    denom = 4.0 * np.cosh((E1 - mu) / (2 * kT)) * np.cosh((E2 - mu) / (2 * kT)) * np.cosh(u / 2)
    t = np.tanh(u / 2)
    out = []
    for E in E_eval:
        E3 = E + E1 - E2
        ok = (E3 > 0) & (E3 < mat.E_max)
        k, k3 = float(band.k(E)), band.k(np.clip(E3, 1e-9, None))
        qmax = np.minimum(k1 + k3, k + k2)
        qmin = np.maximum(np.abs(k1 - k3), np.abs(k - k2))
        phi = np.where(ok & (qmax > qmin), _antiderivative(qmax, kappa) - _antiderivative(qmin, kappa), 0.0)
        base = phi * m12 * band.mass_factor(np.clip(E3, 0, None)) / denom
        norm = (2 * k_mu / k) / phi_mu * (1 + 2 * mat.C * mu) ** -3 / kT**2 * dE * dE
        A = [float(norm * np.sum(base * (-0.5 * t) ** p)) for p in range(3)]
        x = E - mu
        out.append(dict(E=float(E), A0=A[0], A1=A[1], A2=A[2],
                        rate_fl=float(K * (A[0] * kT**2 + A[1] * kT * x + A[2] * x * x)),
                        rate_fl_standard=float(K * ((math.pi * kT) ** 2 + x * x))))
    return dict(T=T, mu=mu, kappa=kappa, K=K, points=out)


def main():
    mat = calibrate_linear_optics(ito_sarkar2023(), N)  # eps_inf, gamma_imp on this grid, D = 17.2 eV
    grid = EnergyGrid(mat.band, mat.E_max, N)
    ee = ElectronElectronB1(grid, mat)
    EF = mat.band.fermi_energy(mat.n_e)
    idx = {name: int(np.argmin(np.abs(grid.E - x * EF))) for name, x in (("0.5 E_F", 0.5), ("1 E_F", 1.0), ("1.5 E_F", 1.5))}
    idx["E_F + 1 eV"] = int(np.argmin(np.abs(grid.E - (EF + 1.0))))
    res = dict(paper_readings=PAPER, material=dict(eps_inf=mat.eps_inf, gamma_imp=mat.gamma_imp, D_def=mat.D_def, N=N))

    # 1. B1 rate for thermal distributions
    ee_rows = []
    for T in (300.0, 1000.0, 2000.0, 3000.0, 4500.0, 6000.0):
        r = ee.rate(fd_of_density(grid, mat.n_e, T))
        ee_rows.append(dict(T=T, **{k: float(r[i]) for k, i in idx.items()}))
    res["ee_B1_thermal"] = ee_rows

    # 2. generalised Fermi liquid
    E_eval = [0.5 * EF, EF, 1.5 * EF]
    fl = []
    for T in (300.0, 1000.0, 3000.0, 5000.0, 7000.0):
        g = generalised_fl(grid, mat, T, E_eval)
        r = ee.rate(fd_of_density(grid, mat.n_e, T))
        for p in g["points"]:
            p["rate_B1"] = float(np.interp(p["E"], grid.E, r))
        fl.append(g)
    res["generalised_fl"] = fl

    # 3. e-ph rate, two deformation potentials
    f300 = fd_of_density(grid, mat.n_e, 300.0)
    eph_rows = {}
    for D in (8.6, 17.2):
        m = mat.with_(D_def=D)
        r = ElectronPhonon(grid, m).rate(f300, 300.0)
        eph_rows[f"D={D}"] = {name: float(np.interp(x * EF, grid.E, r)) for name, x in
                              (("0.1 E_F", 0.1), ("1 E_F", 1.0), ("2 E_F", 2.0), ("3 E_F", 3.0))}
        G = ElectronPhonon(grid, m).coupling(mat.n_e, 1000.0, 300.0) / c.W_PER_M3_IN_WORKING
        eph_rows[f"D={D}"]["G_W_m3_K_at_T_e_1000K"] = float(G)
    res["eph_rate"] = eph_rows

    # 4. damping budget and permittivity at 230 THz, 300 K
    kappa300 = thomas_fermi_kappa(grid, mat, f300)
    w = c.frequency_thz_to_omega(230.0)
    r_ee = ee.rate(f300)
    budgets = {}
    for label, D, imp in (("enzkin calibration to [15] (D 17.2 eV, constant gamma_imp)", 17.2, None),
                          ("Table S1 of Un 2026 (D 8.6 eV, Brooks-Herring n_imp 0.05, Z 1)", 8.6, "BH")):
        m = mat.with_(D_def=D) if imp is None else mat.with_(D_def=D, n_imp=0.05, Z_imp=1.0, gamma_imp=0.0)
        r_eph = ElectronPhonon(grid, m).rate(f300, 300.0)
        r_imp = brooks_herring_rate(grid, m, kappa300) if imp == "BH" else np.full(N, mat.gamma_imp)
        eta = r_ee + r_eph + r_imp
        eps = complex(optics.permittivity(grid, m, f300, w, eta))
        wp2, gam = optics.drude_equivalent(grid, m, f300, eta)
        i = idx["1 E_F"]
        budgets[label] = dict(eta_ee_EF=float(r_ee[i]), eta_eph_EF=float(r_eph[i]), eta_imp_EF=float(r_imp[i]),
                              Gamma_drude=gam, eps_230THz=[eps.real, eps.imag], eps_inf_used=m.eps_inf)
    res["damping_300K"] = budgets
    # cooling time scale C_e / G at T_e = 1000 K for both deformation potentials
    ce = heat_capacity(grid, mat.n_e, 1000.0)  # eV / nm^3 / K
    res["heat_capacity_1000K_J_m3_K"] = float(ce / c.J_PER_M3_IN_EV_PER_NM3)
    for D in (8.6, 17.2):
        G = eph_rows[f"D={D}"]["G_W_m3_K_at_T_e_1000K"] * c.W_PER_M3_IN_WORKING
        eph_rows[f"D={D}"]["C_e_over_G_fs_at_1000K"] = float(ce / G)

    # 5. recovery of a standard transient with the two parameter sets (N = 128, 4 eV band)
    tr_rows = {}
    pulse = Pulse(1300.0, 100.0)
    for label, m in (("[14, 15] set (D 17.2 eV)", calibrate_linear_optics(ito_sarkar2023(), 128)),
                     ("[62] Table S1 set (D 8.6 eV, Brooks-Herring)", ito_un2026())):
        model = KineticModel(m, N=128, drive=PowerDrive.from_absorbed_energy(pulse, 1e8))
        tr = run(model, -250.0, 4000.0, 1.0, sample_every=4, omega_probe=pulse.omega)
        row = dict(Te_max=float(tr["T_e"].max()))
        for o in ("eps_re", "T_e"):
            y = tr[o] - tr[o][0]
            i = int(np.argmax(np.abs(y)))
            after = np.nonzero(np.abs(y[i:]) < np.abs(y[i]) / math.e)[0]
            row[f"{o}_peak_change"] = float(y[i])
            row[f"{o}_recovery_1e_fs"] = float(tr.t[i + after[0]] - tr.t[i]) if after.size else None
        tr_rows[label] = row
    res["transient_100fs_1e8"] = tr_rows

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1))

    print("1. B1 e-e rate (1/fs) for Fermi-Dirac f:")
    for r in ee_rows:
        print(f"   T = {r['T']:5.0f} K  " + "  ".join(f"{k}: {v:.4f}" for k, v in r.items() if k != "T"))
    print("   paper:", PAPER["ee_rate_4500K"])
    print("2. generalised Fermi liquid at E_F (paper Fig. S1 by eye in brackets) and rates at 0.5/1/1.5 E_F:")
    for g in fl:
        p = g["points"][1]
        T = int(g["T"])
        print(f"   T = {T:5d} K  A0 {p['A0']:5.2f} [{PAPER['A0_at_EF'][T]}]  A1 {p['A1']:+5.2f} [{PAPER['A1_at_EF'][T]}]  "
              f"A2 {p['A2']:4.2f} [{PAPER['A2_at_EF'][T]}]   " +
              "  ".join(f"{q['E'] / EF:.1f}E_F: FL {q['rate_fl']:.4f} / B1 {q['rate_B1']:.4f} / std {q['rate_fl_standard']:.4f}"
                        for q in g["points"]))
    print("3. e-ph rate at 300 K (1/fs); paper Fig. 5c by eye:", PAPER["eph_rate"])
    for k, v in eph_rows.items():
        print("  ", k, {a: (f"{b:.4f}" if b < 1 else f"{b:.3e}") for a, b in v.items()})
    print("4. damping budget at E_F, 300 K:")
    for k, v in budgets.items():
        print(f"   {k}: ee {v['eta_ee_EF']:.4f}, eph {v['eta_eph_EF']:.4f}, imp {v['eta_imp_EF']:.4f} 1/fs; "
              f"Drude Gamma {v['Gamma_drude']:.3f} 1/fs; eps(230 THz) = {v['eps_230THz'][0]:.3f} + {v['eps_230THz'][1]:.3f}i")
    print("5. 100 fs pulse, 1e8 J/m^3:")
    for k, v in tr_rows.items():
        print("  ", k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items()})
    print("wrote", OUT)


if __name__ == "__main__":
    main()
