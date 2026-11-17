"""Physical constants in the package's working units.

Working units ("eV-fs-nm"):

    energy   eV          time    fs          length  nm
    mass     eV fs^2 / nm^2                  temperature  K
    density  1 / nm^3    energy density  eV / nm^3    power density  eV / (nm^3 fs)
    field    V / nm      (so that e * E has units eV / nm)

Everything inside :mod:`enzkin` uses these units; the helpers at the bottom convert to and
from SI at the boundaries (pulse intensities, heat capacities, coupling constants).
"""
from __future__ import annotations

import math

# --- CODATA values -------------------------------------------------------------------
HBAR = 0.6582119569  # eV fs
KB = 8.617333262e-5  # eV / K
C_LIGHT = 299.792458  # nm / fs
M_E = 510998.95 / C_LIGHT**2  # electron mass, eV fs^2 / nm^2  (= 5.68563 eV fs^2/nm^2)
E2_OVER_EPS0 = 4.0 * math.pi * 1.439964548  # e^2 / eps0, eV nm  (= 18.0951 eV nm)
EPS0 = 1.0 / E2_OVER_EPS0  # eps0 in e^2 / (eV nm); eps0 * E^2 with E in V/nm gives eV / nm^3
H2M_E = HBAR**2 / (2.0 * M_E)  # hbar^2 / (2 m_e), eV nm^2  (= 0.0380998)

# --- SI conversion factors -----------------------------------------------------------
EV_IN_J = 1.602176634e-19
J_PER_M3_IN_EV_PER_NM3 = 1.0 / (EV_IN_J * 1e27)  # 1 J/m^3  -> eV/nm^3
KG_PER_M3_IN_WORKING = 1.0 / EV_IN_J * 1e30 / 1e18 / 1e27  # 1 kg/m^3 -> eV fs^2 / nm^5
W_PER_M3_IN_WORKING = J_PER_M3_IN_EV_PER_NM3 * 1e-15  # 1 W/m^3 -> eV / (nm^3 fs)
Z0_OHM = 376.730313668  # vacuum impedance


def intensity_to_field(intensity_w_cm2: float, n_refr: float = 1.0) -> float:
    """Peak field amplitude |E0| in V/nm of a plane wave of intensity I = n |E0|^2 / (2 Z0)."""
    i_si = intensity_w_cm2 * 1e4  # W/m^2
    return math.sqrt(2.0 * Z0_OHM * i_si / n_refr) * 1e-9


def field_to_intensity(field_v_nm: float, n_refr: float = 1.0) -> float:
    """Inverse of :func:`intensity_to_field`; returns W/cm^2."""
    return n_refr * (field_v_nm * 1e9) ** 2 / (2.0 * Z0_OHM) * 1e-4


def wavelength_to_photon_energy(lambda_nm: float) -> float:
    """Photon energy in eV for a vacuum wavelength in nm."""
    return 2.0 * math.pi * HBAR * C_LIGHT / lambda_nm


def frequency_thz_to_omega(f_thz: float) -> float:
    """Angular frequency in rad/fs from a frequency in THz."""
    return 2.0 * math.pi * f_thz * 1e-3
