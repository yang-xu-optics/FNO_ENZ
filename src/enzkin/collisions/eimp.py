"""Electron-impurity scattering (Brooks-Herring, screened Coulomb), damping only.

Elastic scattering on ionised donors leaves an isotropic f(E) unchanged, so it enters the model only
through the damping rate eta(E) of the optical response.  The rate is the out-scattering rate (the
functional derivative of the collision term), as for the e-e and e-ph contributions:

    eta_imp(E) = (2 pi / hbar) n_imp (Z e^2 / eps0 eps_b)^2 / (2 pi)^3 * (pi / k^2) (dk/dE)
                 * [ (1 + kappa^2 / 2k^2)^2 - 1 ]^-1

which is Eq. (S9) of Un, Sarkar & Sivan, Nanophotonics 15, e70272 (2026) (ref. [62]) with the sign inside
the bracket corrected (printed as 1 - kappa^2/2k^2, which makes the rate negative for kappa < 2k).
The corrected form follows from integrating |V(q)|^2 = (Z e^2/eps0 eps_b)^2 / (q^2 + kappa^2)^2 over the
elastic sphere; see docs/theory/un2026_notes.md.  kappa is the Thomas-Fermi wavevector of the current f,
screened by eps_b as in the B1 e-e integral.
"""
from __future__ import annotations

import math

import numpy as np

from .. import constants as c
from ..grid import EnergyGrid
from ..material import Material


def brooks_herring_rate(grid: EnergyGrid, mat: Material, kappa: float) -> np.ndarray:
    """eta_imp(E) on the grid, 1/fs.  Zero when ``mat.n_imp`` is zero."""
    if mat.n_imp <= 0.0:
        return np.zeros(grid.N)
    band = mat.band
    k = grid.k
    dkdE = band.mass_factor(grid.E) / (2.0 * band.h2m * k)  # from hbar^2 k^2 / 2 m0 = E (1 + C E)
    v0 = mat.Z_imp * c.E2_OVER_EPS0 / mat.eps_b  # eV nm
    pref = 2.0 * math.pi * mat.n_imp / c.HBAR * v0**2 / (2.0 * math.pi) ** 3
    b = kappa**2 / (2.0 * k**2)
    return pref * math.pi / k**2 * dkdE / ((1.0 + b) ** 2 - 1.0)
