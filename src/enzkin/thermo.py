"""Fermi-Dirac thermodynamics on an :class:`~enzkin.grid.EnergyGrid`.

All functions use the grid's own quadrature, so a Fermi-Dirac function "matched in density and energy"
to a non-equilibrium f(E) matches the *discrete* moments exactly.  That exactness matters twice:
(i) the extracted electron temperature of [14, Eq. (19)] / [15] is well defined, and (ii) WP2 feeds the
network the deviation f - f_FD, which must carry zero number and zero energy.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import constants as c
from .grid import EnergyGrid


def fermi(E, mu: float, T: float):
    """Fermi-Dirac occupation; T in kelvin, overflow-safe."""
    x = (np.asarray(E, dtype=float) - mu) / (c.KB * max(T, 1e-9))
    return 0.5 * (1.0 - np.tanh(0.5 * x))


def mu_of_density(grid: EnergyGrid, n: float, T: float) -> float:
    """Chemical potential (eV) of the Fermi-Dirac distribution with density n at temperature T."""
    span = 60.0 * c.KB * T + 1.0

    def g(mu):
        return grid.density(fermi(grid.E, mu, T)) - n

    lo, hi = -span, grid.E_max + span
    while g(lo) > 0:
        lo -= span
    while g(hi) < 0:
        hi += span
    return float(brentq(g, lo, hi, xtol=1e-14, rtol=1e-13))


def fd_of_density(grid: EnergyGrid, n: float, T: float):
    return fermi(grid.E, mu_of_density(grid, n, T), T)


def energy_of_T(grid: EnergyGrid, n: float, T: float) -> float:
    return grid.energy(fd_of_density(grid, n, T))


def _match_brentq(grid: EnergyGrid, n: float, u: float, T_bounds=(1.0, 5e6)):
    def g(T):
        return energy_of_T(grid, n, T) - u

    lo, hi = T_bounds
    if g(lo) >= 0:
        T = lo
    elif g(hi) <= 0:
        T = hi
    else:
        T = float(brentq(g, lo, hi, xtol=1e-9, rtol=1e-12))
    return mu_of_density(grid, n, T), T


def _match_newton(grid: EnergyGrid, n: float, u: float, mu: float, T: float, tol: float = 1e-12, max_iter: int = 30):
    """Two-dimensional Newton solve for (mu, T); returns None if it does not converge."""
    E, w = grid.E, grid.rho * grid.dE
    for _ in range(max_iter):
        kT = c.KB * T
        f = fermi(E, mu, T)
        h = f * (1.0 - f)
        r1 = np.sum(w * f) - n
        r2 = np.sum(w * E * f) - u
        if abs(r1) <= tol * n and abs(r2) <= tol * abs(u):
            return mu, T
        a11 = np.sum(w * h) / kT
        a12 = np.sum(w * h * (E - mu)) / (kT * T)
        a21 = np.sum(w * E * h) / kT
        a22 = np.sum(w * E * h * (E - mu)) / (kT * T)
        det = a11 * a22 - a12 * a21
        if not np.isfinite(det) or det == 0.0:
            return None
        d_mu = (-r1 * a22 + r2 * a12) / det
        d_T = (-a11 * r2 + a21 * r1) / det
        lam = 1.0
        while T + lam * d_T <= 0.2 * T or abs(lam * d_mu) > 2.0:  # stay in the basin
            lam *= 0.5
            if lam < 1e-6:
                return None
        mu, T = mu + lam * d_mu, T + lam * d_T
    return None


def match_fermi_dirac(grid: EnergyGrid, f, guess=None, T_bounds=(1.0, 5e6)):
    """Return (mu, T, f_FD) of the Fermi-Dirac function with the same density and energy as ``f``.

    This is the "extracted electron temperature" of refs. [14, 15].  ``guess = (mu, T)`` warm-starts a
    Newton iteration (a few O(N) passes); without it, or if Newton fails, a bracketing solve is used.  On a
    finite band the energy of a Fermi-Dirac function saturates as T -> infinity; if ``f`` carries more
    energy than that the upper bound is returned.
    """
    n, u = grid.density(f), grid.energy(f)
    sol = _match_newton(grid, n, u, *guess) if guess is not None else None
    if sol is None:
        mu, T = _match_brentq(grid, n, u, T_bounds)
        if guess is None:  # polish to full precision
            sol = _match_newton(grid, n, u, mu, T)
    mu, T = sol if sol is not None else (mu, T)
    return mu, T, fermi(grid.E, mu, T)


def heat_capacity(grid: EnergyGrid, n: float, T: float, rel_step: float = 1e-3) -> float:
    """Electronic heat capacity C_e = dU/dT at fixed density, eV/(nm^3 K)."""
    h = T * rel_step
    return (energy_of_T(grid, n, T + h) - energy_of_T(grid, n, T - h)) / (2.0 * h)


def entropy_density(grid: EnergyGrid, f) -> float:
    """Fermionic entropy density  s/k_B = -int rho [f ln f + (1-f) ln(1-f)] dE,  1/nm^3."""
    fc = np.clip(f, 1e-300, 1.0)
    gc = np.clip(1.0 - f, 1e-300, 1.0)
    return float(-np.sum(grid.rho * (f * np.log(fc) + (1.0 - f) * np.log(gc))) * grid.dE)
