"""Isotropic non-parabolic (Kane) conduction band.

    hbar^2 k^2 / (2 m0) = gamma(E) = E (1 + C E)            [14, Eq. (1)], [15, Eq. (7)]

All derived quantities follow from this single relation:

    k(E)        = sqrt(gamma / h2m)                          h2m = hbar^2 / (2 m0)
    rho(E)      = (1 + 2 C E) k / (2 pi^2 h2m)               density of states incl. spin  [14, Eq. (2)]
    v(E)        = 2 h2m k / (hbar (1 + 2 C E))               group velocity
    m_opt(E)    = m0 (1 + 2 C E)                             momentum (optical) mass  hbar k / v
    w_opt(E)    = (1 + 8 C gamma / 3) / (1 + 2 C E)^3        weight in omega_p^2 = e^2/(eps0 m0) int rho f w_opt dE
                                                             [15, Eq. (9) and its footnote]
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import constants as c


@dataclass(frozen=True)
class KaneBand:
    """Kane band with band-edge mass ``m0_rel`` (units of m_e) and non-parabolicity ``C`` (1/eV)."""

    m0_rel: float = 0.3964
    C: float = 0.4191

    @property
    def h2m(self) -> float:
        """hbar^2 / (2 m0) in eV nm^2."""
        return c.H2M_E / self.m0_rel

    @property
    def m0(self) -> float:
        """Band-edge mass in working units (eV fs^2 / nm^2)."""
        return self.m0_rel * c.M_E

    def gamma(self, E):
        E = np.asarray(E, dtype=float)
        return E * (1.0 + self.C * E)

    def mass_factor(self, E):
        """1 + 2 C E  (= d gamma / dE = m_opt(E) / m0)."""
        return 1.0 + 2.0 * self.C * np.asarray(E, dtype=float)

    def k(self, E):
        """Wavevector magnitude in 1/nm."""
        return np.sqrt(np.maximum(self.gamma(E), 0.0) / self.h2m)

    def dos(self, E):
        """Density of states per volume and energy, both spins, 1/(eV nm^3)."""
        return self.mass_factor(E) * self.k(E) / (2.0 * np.pi**2 * self.h2m)

    def velocity(self, E):
        """Group velocity in nm/fs."""
        return 2.0 * self.h2m * self.k(E) / (c.HBAR * self.mass_factor(E))

    def optical_weight(self, E):
        """w_opt(E): contribution of an electron at E to omega_p^2, relative to a band-edge electron."""
        return (1.0 + 8.0 * self.C * self.gamma(E) / 3.0) / self.mass_factor(E) ** 3

    # --- zero-temperature relations ------------------------------------------------
    def k_fermi(self, n: float) -> float:
        """Fermi wavevector (1/nm) for electron density n (1/nm^3)."""
        return float((3.0 * np.pi**2 * n) ** (1.0 / 3.0))

    def energy_of_k(self, k):
        """Invert the dispersion: E(k)."""
        g = self.h2m * np.asarray(k, dtype=float) ** 2
        if self.C == 0.0:
            return g
        return (np.sqrt(1.0 + 4.0 * self.C * g) - 1.0) / (2.0 * self.C)

    def fermi_energy(self, n: float) -> float:
        """Zero-temperature Fermi energy (eV) for electron density n (1/nm^3)."""
        return float(self.energy_of_k(self.k_fermi(n)))

    def density_of_fermi_energy(self, e_f: float) -> float:
        return float(self.k(e_f) ** 3 / (3.0 * np.pi**2))
