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

Collision rate for the optical damping
--------------------------------------
The rate 1/tau_e-ph(E) of [14] Eq. (11) counts every phonon scattering event, including quasi-elastic
ones whose phonon energy is smaller than a grid cell.  Those events do not move population between cells,
so the master equation above rightly omits them, but they relax momentum and dominate the rate (the
integrand eps^2 n_B(eps) ~ eps k_B T is weighted towards small eps).  :meth:`ElectronPhonon.rate`
therefore (i) adds the intra-cell (s = 0) contribution and (ii) integrates the Bose factor over each hat
function instead of sampling it at the node.  With this the rate at E_F is within a few per cent of the
continuum value already at N = 128; sampled at the nodes it is three times too small there.
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
        self.W, eps, self._lim = [], [], []
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
            self._lim.append((lo, hi))
            eps.append(s * dE)
        self.eps = np.array(eps)
        # intra-cell (s = 0) quasi-elastic scattering: enters the collision rate only
        hi0 = np.minimum(m.E_debye, 2.0 * hv * g.k) if self.momentum_conservation else np.full(g.N, m.E_debye)
        self._lim0 = (np.zeros(g.N), hi0)
        self._W0 = m.kappa_eph * g.mfac**2 * _hat_eps2_integral(0, dE, *self._lim0)
        self._bose_cache = None

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

    def _bose_weights(self, T_ph: float):
        """Hat-integrated Bose factors  kappa m_i m_j int eps^2 n_B(eps, T_ph) hat_s(eps) d eps  for s = 0..S."""
        if self._bose_cache is not None and abs(self._bose_cache[0] - T_ph) < 1e-3:
            return self._bose_cache[1]
        g, dE, kT = self.grid, self.grid.dE, c.KB * T_ph
        x, wq = np.polynomial.legendre.leggauss(10)

        def integral(s, lo, hi):
            tot = np.zeros_like(lo)
            for a0, b0 in (((s - 1) * dE, s * dE), (s * dE, (s + 1) * dE)):  # rising and falling half of hat_s
                a = np.clip(lo, max(a0, 0.0), max(b0, 0.0))
                b = np.clip(hi, max(a0, 0.0), max(b0, 0.0))
                half = 0.5 * (b - a)
                e = a[:, None] + half[:, None] * (x[None, :] + 1.0)
                hat = np.clip(1.0 - np.abs(e - s * dE) / dE, 0.0, None)
                with np.errstate(divide="ignore", invalid="ignore"):
                    nb = np.where(e > 0, 1.0 / np.expm1(np.maximum(e, 1e-300) / kT), 0.0)
                tot += half * np.sum(wq[None, :] * e * e * nb * hat, axis=1)
            return tot

        out = [self.mat.kappa_eph * g.mfac**2 * integral(0, *self._lim0)]
        for s in range(1, len(self.W) + 1):
            out.append(self.mat.kappa_eph * g.mfac[:-s] * g.mfac[s:] * integral(s, *self._lim[s - 1]))
        self._bose_cache = (T_ph, out)
        return out

    def rate(self, f: np.ndarray, T_ph: float) -> np.ndarray:
        """Collision rate 1/tau_e-ph(E) of ref. [14] Eq. (11), 1/fs, including quasi-elastic scattering.

            1/tau_i = (1/rho_i) sum_s [ B_s(i, i+s) + W_s(i, i+s) f_{i+s}  +  B_s(i-s, i) + W_s(i-s, i) (1 - f_{i-s}) ]

        with W_s = kappa m m' int eps^2 hat_s and B_s the same integral weighted by n_B(eps, T_ph).  The s = 0
        term (partner in the same cell) contributes 2 B_0 + W_0.
        """
        g = self.grid
        B = self._bose_weights(T_ph)
        r = 2.0 * B[0] + self._W0
        for s, w in enumerate(self.W, start=1):
            r[:-s] += B[s] + w * f[s:]  # partner above: absorb a phonon into it, or it emits into E
            r[s:] += B[s] + w * (1.0 - f[:-s])  # partner below: emit a phonon into it
        return r / g.rho

    def coupling(self, n: float, T_e: float, T_ph: float) -> float:
        """Two-temperature coupling G = power / (T_e - T_ph) for a Fermi-Dirac f, eV/(nm^3 fs K)."""
        from ..thermo import fd_of_density

        _, p = self.rhs(fd_of_density(self.grid, n, T_e), T_ph)
        return p / (T_e - T_ph)
