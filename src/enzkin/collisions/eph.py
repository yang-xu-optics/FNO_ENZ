"""Electron-phonon collision term: acoustic deformation potential, Debye spectrum, momentum conserved.

Continuum form, ref. [14] Eq. (9) (phonon emission and absorption at phonon energy eps, lattice at T_ph):

    (df/dt)(E) = D^2/(4 pi rho (hbar v_s)^4) * m0/(hbar k(E)) * int_0^{eps_max(E)} eps^2 d eps
                 { (n_B+1) [ (1-f(E)) f(E+eps) m(E+eps) - (1-f(E-eps)) f(E) m(E-eps) ]
                 +   n_B   [ (1-f(E)) f(E-eps) m(E-eps) - (1-f(E+eps)) f(E) m(E+eps) ] }

with m(E) = 1 + 2 C E, n_B = n_B(eps, T_ph) and eps_max(E) = min(E_D, 2 hbar v_s k(E)): momentum
conservation forbids phonons with q > 2k, which in ITO (k_F ~ q_D / 4) removes most of the Debye
sphere and lowers the energy-transfer rate ~30-fold [14].

Discrete form used here
-----------------------
Multiplying by rho(E) shows the kernel is symmetric in the two electron states:

    d(rho f)/dt |_E = kappa m(E) int d eps eps^2 m(E') [...] ,   kappa = m0^2 D^2 / (4 pi^3 hbar^3 rho (hbar v_s)^4)

so the term is a master equation between pairs of states (E_i, E_j = E_i + s dE) exchanging a phonon of
energy eps_s = s dE.  For each pair

    T_ij = W_ij [ (n_s + 1) f_j (1 - f_i) - n_s f_i (1 - f_j) ]     (net rate j -> i, per eV per nm^3 per fs)
    W_ij = kappa m_i m_j * int_{eps_lo}^{eps_hi} eps^2 hat_s(eps) d eps

where hat_s is the linear hat function centred on eps_s (so phonon energies between grid nodes are
represented, and the first energy moment int eps^3 is integrated exactly), and the limits are the
*exact* triangle condition  hbar v_s |k_j - k_i| <= eps <= min(E_D, hbar v_s (k_i + k_j)).
Using the symmetric bound k_i + k_j instead of 2 k(E) is the only departure from [14]; it differs at
O(eps/E) and is what makes the following hold to round-off rather than approximately:

* electron number is conserved exactly;
* a Fermi-Dirac distribution at T = T_ph is an exact fixed point (detailed balance);
* the energy leaving the electrons, sum_pairs T_ij eps_s, is known exactly and heats the lattice.

Cost: O(N * E_D/dE), negligible next to the electron-electron term.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import constants as c
from ..grid import EnergyGrid
from ..material import Material


def _hat_eps2_integral(s: int, dE: float, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """int_{lo}^{hi} eps^2 hat_s(eps) d eps with hat_s(eps) = max(0, 1 - |eps - s dE|/dE)."""
    c0 = s * dE

    def F_rise(x):  # antiderivative of eps^2 (eps - (c0 - dE))/dE on the rising side
        return (x**4 / 4.0 - (c0 - dE) * x**3 / 3.0) / dE

    def F_fall(x):  # antiderivative of eps^2 ((c0 + dE) - eps)/dE on the falling side
        return ((c0 + dE) * x**3 / 3.0 - x**4 / 4.0) / dE

    a1, b1 = np.clip(lo, c0 - dE, c0), np.clip(hi, c0 - dE, c0)
    a2, b2 = np.clip(lo, c0, c0 + dE), np.clip(hi, c0, c0 + dE)
    out = np.maximum(F_rise(b1) - F_rise(a1), 0.0) + np.maximum(F_fall(b2) - F_fall(a2), 0.0)
    return np.where(hi > lo, out, 0.0)


@dataclass
class ElectronPhonon:
    grid: EnergyGrid
    mat: Material
    momentum_conservation: bool = True  # False reproduces the "no m.c." curves of [14] Fig. 1(b), 5(c)
    W: list = field(init=False, repr=False)  # W[s-1][i] couples cells i and i+s
    eps: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        g, m = self.grid, self.mat
        dE = g.dE
        S = int(np.floor(m.E_debye / dE)) + 1
        S = min(S, g.N - 1)
        hv = m.hbar_v_sound
        self.W, eps = [], []
        for s in range(1, S + 1):
            ki, kj = g.k[:-s], g.k[s:]
            if self.momentum_conservation:
                lo = hv * np.abs(kj - ki)
                hi = np.minimum(m.E_debye, hv * (ki + kj))
            else:
                lo = np.zeros_like(ki)
                hi = np.full_like(ki, m.E_debye)
            w = m.kappa_eph * g.mfac[:-s] * g.mfac[s:] * _hat_eps2_integral(s, dE, lo, hi)
            self.W.append(w)
            eps.append(s * dE)
        self.eps = np.array(eps)

    def _bose(self, T_ph: float) -> np.ndarray:
        return 1.0 / np.expm1(self.eps / (c.KB * T_ph))

    def rhs(self, f: np.ndarray, T_ph: float):
        """Return (df/dt, power to lattice).

        df/dt in 1/fs; power in eV/(nm^3 fs), positive when electrons heat the lattice.
        """
        g = self.grid
        nb = self._bose(T_ph)
        dN = np.zeros_like(f)  # d(rho f)/dt
        power = 0.0
        for s, w in enumerate(self.W, start=1):
            fi, fj = f[:-s], f[s:]
            T = w * ((nb[s - 1] + 1.0) * fj * (1.0 - fi) - nb[s - 1] * fi * (1.0 - fj))
            dN[:-s] += T
            dN[s:] -= T
            power += self.eps[s - 1] * np.sum(T)
        return dN / g.rho, power * g.dE

    def rate(self, f: np.ndarray, T_ph: float) -> np.ndarray:
        """Collision rate 1/tau_e-ph(E) = -delta(df/dt)/delta f(E), ref. [14] Eq. (11); 1/fs."""
        g = self.grid
        nb = self._bose(T_ph)
        r = np.zeros_like(f)
        for s, w in enumerate(self.W, start=1):
            r[:-s] += w * (nb[s - 1] + f[s:])  # partner above: absorb to it / it emits into E
            r[s:] += w * (nb[s - 1] + 1.0 - f[:-s])  # partner below: emit to it
        return r / g.rho

    def coupling(self, n: float, T_e: float, T_ph: float) -> float:
        """Two-temperature coupling G = power / (T_e - T_ph) for a Fermi-Dirac f, eV/(nm^3 fs K)."""
        from ..thermo import fd_of_density

        _, p = self.rhs(fd_of_density(self.grid, n, T_e), T_ph)
        return p / (T_e - T_ph)
