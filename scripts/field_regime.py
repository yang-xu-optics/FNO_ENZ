#!/usr/bin/env python3
"""Where each pumping condition sits relative to the cycle-averaged, isotropic kinetic model.

The WP1 Boltzmann model (and WP3 as planned) evolves an isotropic distribution f(E) driven by a
cycle-averaged absorption term.  Two dimensionless numbers say how far a pump is from that picture:

    cycles      pulse FWHM / optical period.  Cycle averaging and a golden-rule absorption rate need
                many cycles; the plan's limit is about 20 fs (5 cycles at 1.3 um).
    dk / k_F    peak displacement of the Fermi sphere by the field, dk = e E_loc / (hbar omega), in units of
                the Fermi wave number.  In a parabolic band a rigid displacement leaves the Drude current
                linear; in the Kane band it changes the electrons' mass within the cycle, which is the
                "fast" intraband nonlinearity of refs. [5-7, 11, 12] that an isotropic f(E) cannot carry.
                Collisions of a displaced sphere also absorb differently from the weak-field rate when
                dk is comparable with k_F.

Also reported: the energy and the mass factor (1 + 2CE) of the electron at the leading edge of the
displaced sphere (k = k_F + dk along the field), against (1 + 2C E_F).

Local fields: the [15] geometry is computed with the same thin-film stack as scripts/validate_un2023.py,
at the cold permittivity (the enhancement falls as the film heats, so these are upper values at the pulse
peak).  For [20] and [23] only the incident intensity is published; the local field is estimated as the
vacuum field times the normal-incidence Fresnel transmission into ITO with the reference material's cold
permittivity at the pump wavelength.  These are order-of-magnitude estimates, not inputs to any gate.

Output: results/wp1/field_regime.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from enzkin import constants as c
from enzkin.calibrate import calibrate_linear_optics
from enzkin.drive import Pulse, ThinFilmDrive
from enzkin.material import ito_sarkar2023
from enzkin.model import KineticModel

OUT = Path("results/wp1/field_regime.json")
N = 128

# [15] geometry (as in scripts/validate_un2023.py)
THICKNESS, ANGLE = 40.0, 60.0
EPS_GLASS, EPS_AU = 2.25, -80.0 + 8.0j
ANGLE_GLASS = math.degrees(math.asin(math.sin(math.radians(ANGLE)) / math.sqrt(EPS_GLASS)))
GLASS_OVER_AIR = math.cos(math.radians(ANGLE)) / math.cos(math.radians(ANGLE_GLASS))  # as in validate_un2023.py


def cold_eps(mat, wavelength_nm):
    model = KineticModel(mat, N=N)
    f0, _ = model.split(model.initial_state(300.0))
    return complex(model.permittivity(f0, 300.0, 2.0 * math.pi * c.C_LIGHT / wavelength_nm))


def kane_energy(band, k):
    """E from hbar^2 k^2 / 2 m0 = E (1 + C E)."""
    a = c.HBAR ** 2 * k ** 2 / (2.0 * band.m0_rel * c.M_E)
    return (-1.0 + math.sqrt(1.0 + 4.0 * band.C * a)) / (2.0 * band.C)


def main():
    mat = calibrate_linear_optics(ito_sarkar2023(), N)
    band = mat.band
    kF = band.k_fermi(mat.n_e)
    EF = band.fermi_energy(mat.n_e)
    rows = []

    def add(label, source, wavelength, fwhm, e_loc, how):
        hw = c.wavelength_to_photon_energy(wavelength)
        period = wavelength / c.C_LIGHT
        dk = e_loc / hw  # (eV/nm) / eV = 1/nm
        e_lead = kane_energy(band, kF + dk)
        rows.append(dict(label=label, source=source, wavelength_nm=wavelength, fwhm_fs=fwhm, period_fs=period,
                         cycles=fwhm / period, photon_eV=hw, E_loc_V_per_nm=e_loc, local_field=how,
                         dk_per_nm=dk, dk_over_kF=dk / kF, E_lead_eV=e_lead,
                         mass_factor_ratio=(1 + 2 * band.C * e_lead) / (1 + 2 * band.C * EF)))

    # [15]: 220 fs at 1300 nm, incident intensities in GW/cm^2 (air), thin-film stack
    eps13 = cold_eps(mat, 1300.0)
    for i0 in (22.0, 75.0, 200.0):
        drive = ThinFilmDrive(Pulse(1300.0, 220.0), i0 * 1e9 * GLASS_OVER_AIR, thickness=THICKNESS,
                              angle_deg=ANGLE_GLASS, polarization="p", eps_in=EPS_GLASS, eps_out=1.0,
                              layers_after=[(EPS_AU, 100.0)])
        e_loc = c.intensity_to_field(drive.local_intensity(0.0, eps13))  # film-averaged |E_loc| at the peak
        add(f"[15] 220 fs, {i0:g} GW/cm2 incident" + (" (their adiabatic limit)" if i0 == 200.0 else ""),
            "Un et al. 2023 [15]", 1300.0, 220.0, e_loc, "thin-film stack of validate_un2023.py, cold film")

    # sweep 'field' cases: fixed in-film intensity, field from I = |E|^2 / (2 Z0) as in FieldDrive
    for fwhm, inten in ((20.0, 1e11), (20.0, 5e11)):
        add(f"G1b sweep, {fwhm:g} fs, {inten / 1e9:g} GW/cm2 in the film", "convergence_sweep.py", 1300.0, fwhm,
            c.intensity_to_field(inten), "FieldDrive convention")

    # 800 nm few-cycle experiments: vacuum field x Fresnel transmission at normal incidence
    eps08 = cold_eps(mat, 800.0)
    t08 = abs(2.0 / (1.0 + np.sqrt(eps08)))
    for label, src, fwhm, i_tw in (("[20] 6 fs, 800 nm, ~2 TW/cm2", "Lustig et al. 2023 [20]", 6.0, 2.0),
                                   ("[23] 6.4 fs, 800 nm, 5 TW/cm2", "Ozlu et al. 2026 [23]", 6.4, 5.0),
                                   ("[23] 6.4 fs, 800 nm, 11 TW/cm2", "Ozlu et al. 2026 [23]", 6.4, 11.0)):
        add(label, src, 800.0, fwhm, t08 * c.intensity_to_field(i_tw * 1e12), "vacuum field x |t| at normal incidence")

    out = dict(material=dict(n_e=mat.n_e, k_F=kF, E_F=EF, m0_rel=band.m0_rel, C=band.C,
                             eps_1300=[eps13.real, eps13.imag], eps_800=[eps08.real, eps08.imag]),
               rows=rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print(f"k_F = {kF:.3f} /nm, E_F = {EF:.3f} eV, (1+2C E_F) = {1 + 2 * band.C * EF:.2f}")
    print(f"{'condition':<52} {'cycles':>6} {'E_loc':>6} {'dk/kF':>6} {'E_lead':>6} {'m ratio':>7}")
    for r in rows:
        print(f"{r['label']:<52} {r['cycles']:6.1f} {r['E_loc_V_per_nm']:6.2f} {r['dk_over_kF']:6.2f} "
              f"{r['E_lead_eV']:6.2f} {r['mass_factor_ratio']:7.2f}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
