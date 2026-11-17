"""Optical response of the conduction electrons computed directly from f(E) (no temperature assumed).

Plasma frequency, ref. [15] Eq. (10) and the footnote before it (equivalent forms, related by parts):

    omega_p^2[f] = e^2/(eps0 m0) * int rho(E) f(E) w_opt(E) dE
                 = e^2/(3 eps0)  * int rho(E) v(E)^2 (-df/dE) dE          (the form in the project plan)

Permittivity with energy-dependent damping, ref. [15] Eq. (9):

    eps(omega) = eps_inf - e^2/(eps0 m0) * int rho(E) f(E) w_opt(E) / (omega + i eta(E))^2 dE
    eta(E)     = 1/tau_ee(E) + 1/tau_eph(E) + gamma_imp                    (Matthiessen)

Note the damped Lindhard form (omega + i eta)^2: to first order in eta/omega it equals a Drude term
omega_p^2 / (omega (omega + i Gamma)) with Gamma = 2 <eta>, where <.> is the average weighted by
rho f w_opt.  :func:`drude_equivalent` returns that (omega_p^2, Gamma) pair, which is what the current
equation of WP3,  dJ/dt + Gamma J = eps0 omega_p^2 E,  needs.
"""
from __future__ import annotations

import numpy as np

from . import constants as c
from .grid import EnergyGrid
from .material import Material


def plasma_frequency_sq(grid: EnergyGrid, mat: Material, f: np.ndarray) -> float:
    """omega_p^2 in 1/fs^2 (unscreened: eps = eps_inf - omega_p^2/omega^2 in the collisionless limit)."""
    return float(c.E2_OVER_EPS0 / mat.band.m0 * np.sum(grid.rho * f * grid.w_opt) * grid.dE)


def plasma_frequency_sq_gradient_form(grid: EnergyGrid, mat: Material, f: np.ndarray) -> float:
    """Same quantity from the (-df/dE) form; differs from the above only by the boundary term at E_max."""
    v2 = mat.band.velocity(grid.edges) ** 2
    return float(c.E2_OVER_EPS0 / 3.0 * grid.minus_dfdE_weighted(f, grid.rho_edges * v2))


def permittivity(grid: EnergyGrid, mat: Material, f: np.ndarray, omega, eta: np.ndarray):
    """Complex relative permittivity at angular frequency ``omega`` (rad/fs; scalar or array).

    ``eta`` is the total damping rate eta(E) in 1/fs on the grid.
    """
    omega = np.asarray(omega, dtype=float)
    w = c.E2_OVER_EPS0 / mat.band.m0 * grid.rho * f * grid.w_opt * grid.dE
    denom = (omega[..., None] + 1j * eta) ** 2
    return mat.eps_inf - np.sum(w / denom, axis=-1)


def drude_equivalent(grid: EnergyGrid, mat: Material, f: np.ndarray, eta: np.ndarray):
    """Return (omega_p^2 [1/fs^2], Gamma [1/fs]) of the Drude term matching eps to first order in eta/omega."""
    w = grid.rho * f * grid.w_opt
    wp2 = float(c.E2_OVER_EPS0 / mat.band.m0 * np.sum(w) * grid.dE)
    gamma = float(2.0 * np.sum(w * eta) / np.sum(w))
    return wp2, gamma


def optical_mass(grid: EnergyGrid, mat: Material, f: np.ndarray) -> float:
    """Density-averaged optical mass m* = n e^2 / (eps0 omega_p^2), in units of m_e."""
    n = grid.density(f)
    return float(c.E2_OVER_EPS0 * n / plasma_frequency_sq(grid, mat, f) / c.M_E)


def enz_frequency(grid: EnergyGrid, mat: Material, f: np.ndarray, eta: np.ndarray, bracket=(0.3, 6.0)) -> float:
    """Angular frequency (rad/fs) at which Re eps = 0."""
    from scipy.optimize import brentq

    return float(brentq(lambda w: permittivity(grid, mat, f, w, eta).real, *bracket))
