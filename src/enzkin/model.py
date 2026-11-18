"""The kinetic model: one spatial cell of ENZ material.

State   y = [ f(E_0) ... f(E_{N-1}),  T_ph,  U_abs,  U_pair ]

    df/dt    = (df/dt)_exc[f; p_abs]  +  Q_ee[f]  +  Q_eph[f, T_ph]  +  S_inter[f, T_ph, I]
    C_ph dT_ph/dt = (power leaving the electrons through phonons)  (+ interband energy bookkeeping)
    dU_abs/dt  = p_abs (+ interband field power)   cumulative energy taken from the field (bookkeeping)
    dU_pair/dt = interband pair power              energy stored in electron-hole pairs (bookkeeping)

U_abs and U_pair are carried in the state so that the energy audit
    U_e + C_ph T_ph + U_pair - U_abs = const
is a linear invariant of the ODE system and therefore holds to round-off under Runge-Kutta stepping.

This is the equation of the project plan (WP1) and of refs. [14, 15].  The pieces are swappable:

* ``ee``      :class:`~enzkin.collisions.ee.ElectronElectronB1` (ground truth), ``...B0`` (FFT baseline),
              ``...RTA`` (relaxation-time baseline), or - from WP2 on - a learned operator with the same
              ``rhs(f, with_rate=...)`` interface.
* ``drive``   any object from :mod:`enzkin.drive`.
* ``interband``  optional list of source terms from :mod:`enzkin.interband` (models M2-M4).

Everything the optical side needs is derived from f in :meth:`KineticModel.observe`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

from . import constants as c
from . import optics, thermo
from .collisions.ee import ElectronElectronB1, thomas_fermi_kappa
from .collisions.eimp import brooks_herring_rate
from .collisions.eph import ElectronPhonon
from .excitation import Excitation
from .grid import EnergyGrid
from .material import Material


@dataclass
class KineticModel:
    mat: Material
    N: int = 256
    ee: Any = None  # defaults to B1
    drive: Any = None
    excitation_mode: str = "absorption_only"
    interband: Sequence = ()
    eph_momentum_conservation: bool = True
    heat_lattice: bool = True
    grid: EnergyGrid = field(init=False)
    eph: ElectronPhonon = field(init=False)
    exc: Excitation | None = field(init=False, default=None)
    n_rhs: int = field(init=False, default=0)

    def __post_init__(self):
        self.grid = EnergyGrid(self.mat.band, self.mat.E_max, self.N)
        if self.ee is None:
            self.ee = ElectronElectronB1(self.grid, self.mat)
        elif isinstance(self.ee, type):
            self.ee = self.ee(self.grid, self.mat)
        self.eph = ElectronPhonon(self.grid, self.mat, self.eph_momentum_conservation)
        if self.drive is not None:
            self.exc = Excitation(self.grid, self.drive.pulse.photon_energy, self.excitation_mode)

    # --- state helpers -----------------------------------------------------------------
    def initial_state(self, T0: float = 300.0) -> np.ndarray:
        f0 = thermo.fd_of_density(self.grid, self.mat.n_e, T0)
        return np.concatenate([f0, [T0, 0.0, 0.0]])

    def split(self, y):
        """(f, T_ph) views of a state vector."""
        return y[: self.grid.N], float(y[self.grid.N])

    def absorbed_energy(self, y) -> float:
        return float(y[self.grid.N + 1])

    def pair_energy(self, y) -> float:
        return float(y[self.grid.N + 2])

    def energy_audit(self, y) -> float:
        """U_e + C_ph T_ph + U_pair - U_abs  (eV/nm^3); constant along any trajectory."""
        f, T_ph = self.split(y)
        return self.grid.energy(f) + self.mat.C_ph_working * T_ph + self.pair_energy(y) - self.absorbed_energy(y)

    # --- damping and permittivity -------------------------------------------------------
    def damping(self, f, T_ph, rate_ee=None):
        """eta(E) = 1/tau_ee + 1/tau_eph + gamma_imp (+ Brooks-Herring impurities if mat.n_imp > 0), 1/fs."""
        if rate_ee is None:
            rate_ee = self.ee.rate(f)
        extra = sum(term.damping_extra(f) for term in self.interband if hasattr(term, "damping_extra"))
        if self.mat.n_imp > 0.0:
            extra = extra + brooks_herring_rate(self.grid, self.mat, thomas_fermi_kappa(self.grid, self.mat, f))
        return rate_ee + self.eph.rate(f, T_ph) + self.mat.gamma_imp + extra

    def permittivity(self, f, T_ph, omega, rate_ee=None):
        return optics.permittivity(self.grid, self.mat, f, omega, self.damping(f, T_ph, rate_ee))

    # --- right-hand side -----------------------------------------------------------------
    def rhs(self, t: float, y: np.ndarray, return_parts: bool = False, want_stiffness: bool = False):
        """dy/dt.  With ``want_stiffness`` also return an upper estimate (1/fs) of the fastest rate in the
        system, which the integrator uses to choose a stable and accurate step."""
        f, T_ph = self.split(y)
        self.n_rhs += 1
        need_eps = self.drive is not None and getattr(self.drive, "needs_eps", False)
        if need_eps or want_stiffness:
            d_ee, r_ee = self.ee.rhs(f, with_rate=True)
        else:
            d_ee, r_ee = self.ee.rhs(f), None
        d_eph, p_lat = self.eph.rhs(f, T_ph)

        p_abs, d_exc, eps = 0.0, 0.0, None
        if self.drive is not None:
            if need_eps:
                eps = complex(self.permittivity(f, T_ph, self.drive.pulse.omega, r_ee))
            d_exc, p_abs = self.exc.rhs_and_power(f, self.drive.power(t, eps))

        d_int, p_lat_int, p_field_int, p_pair = 0.0, 0.0, 0.0, 0.0
        if self.interband:
            has_field = self.drive is not None and hasattr(self.drive, "local_intensity")
            intensity = self.drive.local_intensity(t, eps) if has_field else 0.0
            for term in self.interband:
                out = term.rhs(f, T_ph, t=t, intensity=intensity, p_abs=p_abs)
                d_int = d_int + out["df"]
                p_lat_int += out["p_lattice"]
                p_field_int += out["p_field"]
                p_pair += out["pair_power"]

        dT = (p_lat + p_lat_int) / self.mat.C_ph_working if self.heat_lattice else 0.0
        if not self.heat_lattice:  # keep the audit exact: energy leaving to a frozen lattice counts as stored
            p_pair += p_lat + p_lat_int
        dy = np.concatenate([d_exc + d_ee + d_eph + d_int, [dT, p_abs + p_field_int, p_pair]])
        if want_stiffness:
            df = dy[: self.grid.N]
            lam = float(np.max(r_ee + self.eph.rate(f, T_ph)))  # collisional relaxation (stability)
            lam = max(lam, float(np.max(np.abs(df))) / 0.02)  # at most ~2 % change of any f(E) per step (accuracy)
            return dy, lam
        if return_parts:
            return dy, dict(exc=d_exc, ee=d_ee, eph=d_eph, inter=d_int, p_abs=p_abs, p_lat=p_lat, eps=eps, rate_ee=r_ee)
        return dy

    # --- diagnostics ------------------------------------------------------------------------
    def observe(self, t: float, y: np.ndarray, omega_probe: float | None = None) -> dict:
        """Scalars derived from the state: moments, extracted temperature, optical response."""
        g = self.grid
        f, T_ph = self.split(y)
        mu, T_e, f_fd = thermo.match_fermi_dirac(g, f)
        r_ee = self.ee.rate(f)
        eta = self.damping(f, T_ph, r_ee)
        wp2, gamma = optics.drude_equivalent(g, self.mat, f, eta)
        out = dict(
            t=t,
            n=g.density(f),
            u=g.energy(f),
            T_ph=T_ph,
            u_abs=self.absorbed_energy(y),
            u_pair=self.pair_energy(y),
            audit=self.energy_audit(y),
            T_e=T_e,
            mu=mu,
            wp2=wp2,
            gamma=gamma,
            m_opt=optics.optical_mass(g, self.mat, f),
            nonthermal=float(np.sum(g.rho * np.abs(f - f_fd)) * g.dE / g.density(f)),
            entropy=thermo.entropy_density(g, f),
            f_min=float(f.min()),
            f_max=float(f.max()),
        )
        omega = omega_probe if omega_probe is not None else (self.drive.pulse.omega if self.drive is not None else None)
        if omega is not None:
            eps = complex(optics.permittivity(g, self.mat, f, omega, eta))
            eps_th = complex(self.permittivity(f_fd, T_ph, omega))
            out.update(eps_re=eps.real, eps_im=eps.imag, eps_th_re=eps_th.real, eps_th_im=eps_th.imag)
            if self.drive is not None:
                out["p_abs"] = self.drive.power(t, eps)
        return out
