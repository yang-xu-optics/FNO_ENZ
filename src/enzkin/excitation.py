"""Intraband photo-excitation term (assisted absorption in steps of the photon energy).

Refs. [14] Eq. (4)-(6) and [15] Eq. (2)-(3):

    (df/dt)_exc(E) = B(t) [ D rho(E - hw) f(E - hw) (1 - f(E))  -  D rho(E + hw) f(E) (1 - f(E + hw)) ]

with a constant squared matrix element D and B(t) fixed by energy conservation,
int E rho (df/dt)_exc dE = p_abs(t), the absorbed power density (Poynting's theorem).  The term moves
electrons up by one photon energy; it does not say which momentum-conserving partner (phonon,
impurity, surface) assisted the transition - only their cumulative rate enters, through p_abs.

Discrete form
-------------
A transition takes an electron from cell j to energy E_j + hw.  Writing hw = (m + w) dE, 0 <= w < 1, the
electron is deposited in cells j+m and j+m+1 with weights (1-w, w).  Number is conserved exactly for any
photon energy, whether or not it is commensurate with the grid.

Two variants:

``absorption_only``  the form printed in [14, 15]: upward transitions only.
``net``              upward minus stimulated downward transitions with the same matrix element, i.e. the
                     Pauli factors combine to f(E - hw) - f(E).  It is the physical form once the
                     distribution is hot enough for emission to matter ([15] Sec. III notes that the
                     absorption-only form piles electrons up at the band top at high intensity).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .grid import EnergyGrid


@dataclass
class Excitation:
    grid: EnergyGrid
    photon_energy: float  # eV
    mode: str = "absorption_only"  # or "net"
    _m: int = field(init=False, repr=False)
    _w: float = field(init=False, repr=False)

    def __post_init__(self):
        if self.mode not in ("absorption_only", "net"):
            raise ValueError(f"unknown excitation mode {self.mode!r}")
        x = self.photon_energy / self.grid.dE
        self._m = int(np.floor(x))
        self._w = float(x - self._m)
        if self._m + 1 >= self.grid.N:
            raise ValueError("photon energy exceeds the energy grid")

    def shape(self, f: np.ndarray):
        """Return (g, e_rate): df/dt and energy rate per unit B, so that B = p_abs / e_rate."""
        g_ = self.grid
        rho = g_.rho
        dN = np.zeros_like(f)  # d(rho f)/dt per unit B
        for shift, weight in ((self._m, 1.0 - self._w), (self._m + 1, self._w)):
            if weight == 0.0 or shift == 0:
                continue
            lo, hi = slice(None, -shift), slice(shift, None)  # lower state j, upper state j + shift
            up = rho[lo] * f[lo] * rho[hi] * (1.0 - f[hi])
            if self.mode == "net":
                up = up - rho[hi] * f[hi] * rho[lo] * (1.0 - f[lo])
            t = weight * up
            dN[lo] -= t
            dN[hi] += t
        e_rate = float(np.sum(g_.E * dN) * g_.dE)
        return dN / rho, e_rate

    def rhs(self, f: np.ndarray, p_abs: float) -> np.ndarray:
        """df/dt in 1/fs for an absorbed power density p_abs in eV/(nm^3 fs)."""
        if p_abs == 0.0:
            return np.zeros_like(f)
        g, e_rate = self.shape(f)
        if e_rate <= 0.0:  # population inverted on the photon scale: no net absorption possible
            return np.zeros_like(f)
        return g * (p_abs / e_rate)
