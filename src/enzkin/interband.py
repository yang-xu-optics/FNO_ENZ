"""Interband channels S_inter[f, T_ph, I]: the terms that change the conduction-band carrier density.

These are the switchable source terms behind the mechanism models of WP4:

    M1  intraband kinetics only                       -> no interband term
    M2  M1 + Auger-type excitation / recombination    -> :class:`AugerLikeRelaxation`      (ref. [22])
    M3  M1 + two-photon absorption with state filling -> :class:`TwoPhotonAbsorption`      (ref. [23])
    M4  both

Both published models are *lumped* (a density and an energy in [22]; three level populations in [23]).
Here they are lifted to the distribution level in the most conservative way that reproduces the published
rate equations when f is thermal - see each class.  Each term returns a dict

    df            contribution to df/dt, 1/fs
    p_lattice     power handed to the lattice, eV/(nm^3 fs)
    p_field       power drawn from the optical field *in addition to* the intraband p_abs
    pair_power    rate of change of the energy stored in electron-hole pairs (gap + hole kinetic energy)

so that the model's energy audit  U_e + C_ph T_ph + U_pair - U_abs = const  stays exact.

Valence band: a single parabolic band of mass ``m_v_rel`` whose maximum lies ``E_gap`` below the conduction
band minimum.  Neither [22] nor [23] prints a valence-band mass or density of states, so ``m_v_rel`` is a
declared free parameter (default 1.0 m_e), not a literature value.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from . import constants as c
from . import thermo
from .grid import EnergyGrid
from .material import Material


def _moment_source(grid: EnergyGrid, f_ref: np.ndarray, dn: float, du: float) -> np.ndarray:
    """Smooth df/dt with prescribed number rate ``dn`` and energy rate ``du``.

    Spanned by h(1-h) and E h(1-h) with h = f_ref: the directions in which a Fermi-Dirac function moves
    when its chemical potential and temperature change.  Applied to a thermal distribution it therefore
    moves it along the Fermi-Dirac family, which is exactly what a lumped (n, U) model such as [22] assumes.
    """
    w = np.clip(f_ref * (1.0 - f_ref), 1e-300, None)
    b1, b2 = w, grid.E * w
    r = grid.rho * grid.dE
    A = np.array([[np.sum(r * b1), np.sum(r * b2)], [np.sum(r * grid.E * b1), np.sum(r * grid.E * b2)]])
    x = np.linalg.solve(A, np.array([dn, du]))
    return x[0] * b1 + x[1] * b2


@dataclass
class ValenceBand:
    """Parabolic valence band, maximum at E = -E_gap on the conduction-band energy axis."""

    E_gap: float  # eV, conduction-band minimum to valence-band maximum
    m_v_rel: float = 1.0
    n_points: int = 400
    depth: float = 6.0  # eV of valence band included in hole integrals

    def __post_init__(self):
        self._eps = (np.arange(self.n_points) + 0.5) * self.depth / self.n_points  # hole kinetic energy
        h2m = c.H2M_E / self.m_v_rel
        self._rho = np.sqrt(self._eps / h2m) / (2.0 * math.pi**2 * h2m)
        self._d = self.depth / self.n_points

    def hole_density(self, mu: float, T: float) -> float:
        """Holes per nm^3 for chemical potential mu (measured from the conduction-band minimum)."""
        x = (mu + self.E_gap + self._eps) / (c.KB * T)
        return float(np.sum(self._rho * 0.5 * (1.0 - np.tanh(0.5 * x))) * self._d)

    def hole_energy(self, mu: float, T: float) -> float:
        x = (mu + self.E_gap + self._eps) / (c.KB * T)
        return float(np.sum(self._eps * self._rho * 0.5 * (1.0 - np.tanh(0.5 * x))) * self._d)


@dataclass
class AugerLikeRelaxation:
    """Auger-type interband excitation and recombination, ref. [22] Eqs. (3) and (5).

        dn_c/dt = eta [ n_0(T_e) - n_c ]            dU_c/dt |_inter = -E_gap dn_c/dt

    n_0(T) is the conduction-electron density the film would have if conduction and valence bands shared
    one chemical potential at temperature T (charge neutrality: n_0 - n_dopant = p).  Hot electrons
    therefore pull carriers across the gap at the cost of E_gap each, and the same term recombines them
    when the electrons cool.  [22] assigns eta = 1/120 fs^-1 (the observed delay; it is not fitted) and
    E_gap = 3.9 eV.

    Kinetic lift: the number and energy rates are injected through :func:`_moment_source` built on the
    Fermi-Dirac function matched to the current f, and T_e is that function's temperature.  With
    ``thermal_generation=False`` the target is the dopant density (recombination only), which serves as the
    recombination channel of model M3.
    """

    grid: EnergyGrid
    mat: Material
    eta: float = 1.0 / 120.0  # 1/fs
    E_gap: float = 3.9  # eV
    m_v_rel: float = 1.0
    thermal_generation: bool = True
    gamma_hole: float = 0.0  # extra optical damping per unit hole density, nm^3/fs  ([22]: damping ~ holes)
    valence: ValenceBand = field(init=False, repr=False)
    _guess: tuple | None = field(init=False, repr=False, default=None)
    _table: tuple | None = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.valence = ValenceBand(self.E_gap, self.m_v_rel)

    def _equilibrium_density_exact(self, T: float) -> float:
        g, n_d = self.grid, self.mat.n_e

        def neutrality(mu):
            return g.density(thermo.fermi(g.E, mu, T)) - n_d - self.valence.hole_density(mu, T)

        span = 60.0 * c.KB * T
        mu = brentq(neutrality, -self.E_gap - self.valence.depth - span, g.E_max + span, xtol=1e-12)
        return g.density(thermo.fermi(g.E, mu, T))

    def equilibrium_density(self, T: float) -> float:
        """n_0(T): conduction density with both bands at a common chemical potential.

        Tabulated once on a logarithmic temperature grid (it depends on T only) and interpolated.
        """
        if self._table is None:
            Ts = np.geomspace(100.0, 3.0e5, 160)
            self._table = (np.log(Ts), np.array([self._equilibrium_density_exact(t) for t in Ts]))
        lt, n0 = self._table
        return float(np.interp(math.log(min(max(T, 100.0), 3.0e5)), lt, n0))

    def hole_density(self, f) -> float:
        return max(self.grid.density(f) - self.mat.n_e, 0.0)

    def damping_extra(self, f) -> float:
        return self.gamma_hole * self.hole_density(f)

    def rhs(self, f, T_ph, **ctx) -> dict:
        g = self.grid
        mu, T_e, f_fd = thermo.match_fermi_dirac(g, f, guess=self._guess)
        self._guess = (mu, T_e)
        n_c = g.density(f)
        n_0 = self.equilibrium_density(T_e) if self.thermal_generation else self.mat.n_e
        rate = self.eta * (n_0 - n_c)
        df = _moment_source(g, f_fd, rate, -self.E_gap * rate)
        return dict(df=df, p_lattice=0.0, p_field=0.0, pair_power=self.E_gap * rate)


@dataclass
class TwoPhotonAbsorption:
    """Interband two-photon absorption with state filling, the mechanism of ref. [23].

    Pair generation rate  G = beta I^2 / (2 hw) * (1 - f(E*))          [pairs / (nm^3 fs)]

    Two pump photons lift a valence electron to the conduction-band energy E* fixed by a vertical
    transition, E* + (m0/m_v) E*(1 + C E*) = 2 hw - E_gap.  In degenerate ITO E* lies below the Fermi level,
    so the channel is Pauli-blocked (f(E*) = 1) until one-photon intraband absorption empties the bottom of
    the band - which the intraband excitation term of this package does on its own.  [23] expresses the
    same physics with lumped populations, rate beta I^2 N_1 (1 - N_2/N_2max), and takes the fundamental gap
    as 2.67 eV.  The channel is closed whenever 2 hw < E_gap (true for pumps beyond ~930 nm at that gap).

    ``beta_cm_per_GW`` is the usual two-photon coefficient defined by dI/dz = -beta I^2 with unblocked
    final states.  Recombination is not part of this term; pair it with
    ``AugerLikeRelaxation(thermal_generation=False)`` (model M3) or the full Auger term (model M4).
    """

    grid: EnergyGrid
    mat: Material
    photon_energy: float  # eV
    beta_cm_per_GW: float = 1.0
    E_gap: float = 2.67
    m_v_rel: float = 1.0
    E_star: float = field(init=False)
    _i: int = field(init=False, repr=False)
    _w: float = field(init=False, repr=False)

    def __post_init__(self):
        delta = 2.0 * self.photon_energy - self.E_gap
        self.E_star = float("nan")
        self._i, self._w = -1, 0.0
        if delta <= 0.0:
            return  # energetically closed
        r = self.mat.m0_rel / self.m_v_rel
        C = self.mat.C
        self.E_star = brentq(lambda e: e + r * e * (1.0 + C * e) - delta, 0.0, delta)
        x = self.E_star / self.grid.dE - 0.5
        self._i = int(np.clip(np.floor(x), 0, self.grid.N - 2))
        self._w = float(np.clip(x - self._i, 0.0, 1.0))

    @property
    def open(self) -> bool:
        return self._i >= 0

    def rhs(self, f, T_ph, intensity: float = 0.0, **ctx) -> dict:
        g = self.grid
        df = np.zeros_like(f)
        if not self.open or intensity <= 0.0:
            return dict(df=df, p_lattice=0.0, p_field=0.0, pair_power=0.0)
        i, w = self._i, self._w
        f_star = (1.0 - w) * f[i] + w * f[i + 1]
        # beta [cm/W] I^2 [W^2/cm^4] / (2 hw [J])  -> pairs / (cm^3 s)  -> / (nm^3 fs)
        gen = (self.beta_cm_per_GW * 1e-9) * intensity**2 / (2.0 * self.photon_energy * c.EV_IN_J) * 1e-21 * 1e-15
        gen *= max(1.0 - f_star, 0.0)
        # deposit in the two cells bracketing E*, weights chosen so the mean deposited energy is E*
        df[i] += gen * (1.0 - w) / (g.rho[i] * g.dE)
        df[i + 1] += gen * w / (g.rho[i + 1] * g.dE)
        e_dep = (1.0 - w) * g.E[i] + w * g.E[i + 1]
        p_field = 2.0 * self.photon_energy * gen
        return dict(df=df, p_lattice=0.0, p_field=p_field, pair_power=p_field - e_dep * gen)
