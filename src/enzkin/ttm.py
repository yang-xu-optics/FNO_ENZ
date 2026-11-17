"""Thermal-electron baseline (model M0): a two-temperature model built from the *same* operators.

    dU_e/dt        = p_abs(t) - P_e-ph(T_e, T_ph)
    C_ph dT_ph/dt  =            P_e-ph(T_e, T_ph)

The electrons are a Fermi-Dirac distribution at all times (instantaneous thermalisation).  Nothing is
linearised or fitted: U_e(T_e) and mu(T_e) are the exact Fermi integrals on the kinetic grid (so C_e(T_e)
and the temperature-dependent chemical potential are built in), and P_e-ph is the kinetic electron-phonon
operator applied to that Fermi-Dirac function (so G(T_e, T_ph) with momentum conservation is built in).
A kinetic run and a thermal run therefore differ *only* by non-thermal electrons, which is the comparison
refs. [15] (TTM vs "ANTHEM") and [38] make and what hypothesis H1 needs.

The class mirrors :class:`enzkin.model.KineticModel` (``initial_state``, ``rhs``, ``observe``, ``split``)
so :func:`enzkin.integrate.run` drives both.
State  y = [ U_e, T_ph, U_abs ].
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import optics, thermo
from .collisions.ee import ElectronElectronB1
from .collisions.eph import ElectronPhonon
from .grid import EnergyGrid
from .material import Material


@dataclass
class ThermalModel:
    mat: Material
    N: int = 256
    ee: Any = None  # only its collision *rate* is used, for the optical damping
    drive: Any = None
    grid: EnergyGrid = field(init=False)
    eph: ElectronPhonon = field(init=False)
    n_rhs: int = field(init=False, default=0)
    _guess: tuple | None = field(init=False, default=None, repr=False)

    def __post_init__(self):
        self.grid = EnergyGrid(self.mat.band, self.mat.E_max, self.N)
        if self.ee is None:
            self.ee = ElectronElectronB1(self.grid, self.mat)
        elif isinstance(self.ee, type):
            self.ee = self.ee(self.grid, self.mat)
        self.eph = ElectronPhonon(self.grid, self.mat)

    # --- state ---------------------------------------------------------------------------
    def initial_state(self, T0: float = 300.0) -> np.ndarray:
        f0 = thermo.fd_of_density(self.grid, self.mat.n_e, T0)
        self._guess = (thermo.mu_of_density(self.grid, self.mat.n_e, T0), T0)
        return np.array([self.grid.energy(f0), T0, 0.0])

    def distribution(self, u_e: float):
        """(mu, T_e, f_FD) for electron energy density u_e at the material's carrier density."""
        g, n = self.grid, self.mat.n_e
        sol = thermo._match_newton(g, n, u_e, *self._guess) if self._guess is not None else None
        if sol is None:
            sol = thermo._match_brentq(g, n, u_e)
        self._guess = sol
        return sol[0], sol[1], thermo.fermi(g.E, *sol)

    def split(self, y):
        """(f, T_ph), so that integrate.run can sample a thermal run like a kinetic one."""
        return self.distribution(float(y[0]))[2], float(y[1])

    # --- dynamics ------------------------------------------------------------------------
    def damping(self, f, T_ph):
        return self.ee.rate(f) + self.eph.rate(f, T_ph) + self.mat.gamma_imp

    def rhs(self, t: float, y: np.ndarray, want_stiffness: bool = False):
        self.n_rhs += 1
        u_e, T_ph = float(y[0]), float(y[1])
        _, T_e, f = self.distribution(u_e)
        _, p_lat = self.eph.rhs(f, T_ph)
        p_abs = 0.0
        if self.drive is not None:
            eps = None
            if getattr(self.drive, "needs_eps", False):
                eps = complex(optics.permittivity(self.grid, self.mat, f, self.drive.pulse.omega, self.damping(f, T_ph)))
            p_abs = self.drive.power(t, eps)
        dy = np.array([p_abs - p_lat, p_lat / self.mat.C_ph_working, p_abs])
        if want_stiffness:
            return dy, abs(dy[0]) / (0.02 * max(u_e, 1e-12))  # <= 2 % change of U_e per step
        return dy

    def observe(self, t: float, y: np.ndarray, omega_probe: float | None = None) -> dict:
        g = self.grid
        u_e, T_ph = float(y[0]), float(y[1])
        mu, T_e, f = self.distribution(u_e)
        eta = self.damping(f, T_ph)
        wp2, gamma = optics.drude_equivalent(g, self.mat, f, eta)
        out = dict(t=t, n=g.density(f), u=u_e, T_ph=T_ph, u_abs=float(y[2]), T_e=T_e, mu=mu, wp2=wp2, gamma=gamma,
                   m_opt=optics.optical_mass(g, self.mat, f), nonthermal=0.0,
                   audit=u_e + self.mat.C_ph_working * T_ph - float(y[2]))
        omega = omega_probe if omega_probe is not None else (self.drive.pulse.omega if self.drive is not None else None)
        if omega is not None:
            eps = complex(optics.permittivity(g, self.mat, f, omega, eta))
            out.update(eps_re=eps.real, eps_im=eps.imag, eps_th_re=eps.real, eps_th_im=eps.imag)
            if self.drive is not None:
                out["p_abs"] = self.drive.power(t, eps)
        return out
