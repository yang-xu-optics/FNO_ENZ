"""Electron-electron collision integrals in energy space, at two fidelity levels (plan, WP1).

Both are Uehling-Uhlenbeck (Pauli-blocked) two-body integrals for an isotropic distribution f(E):

    (df/dt)(E1) = sum over E2, E3, E4 with E1 + E2 = E3 + E4 of
                  W(1,2;3,4) [ f3 f4 (1-f1)(1-f2) - f1 f2 (1-f3)(1-f4) ]

On the uniform grid E1 + E2 = E3 + E4 is the index identity i1 + i2 = i3 + i4, and both kernels are
symmetric under (1,2) <-> (3,4) and (1 <-> 2, 3 <-> 4) once multiplied by rho(E1).  Number and energy
are therefore conserved to round-off, and any Fermi-Dirac function is an exact fixed point.

B1 - momentum-resolved, statically screened Coulomb interaction  (ref. [14] Eq. (15); refs. [32, 33])
----------------------------------------------------------------------------------------------------
Golden rule with |V(q)|^2 = [e^2 / (eps0 eps_b (q^2 + kappa^2))]^2, angular integrals done analytically
for isotropic f, k dk = m0 (1+2CE) dE / hbar^2:

    (df/dt)(E1) = P / k1 * sum_{E2,E3} dE^2 m2 m3 m4 Phi(1,2;3,4) [ ... ]
    P     = c_ee * (1/(hbar (2 pi)^3)) * (e^2/(eps0 eps_b))^2 * (m0 / (hbar^2 kappa))^3
    Phi   = [ q kappa/(q^2+kappa^2) + atan(q/kappa) ] evaluated from q_lo to q_hi   (0 if q_hi <= q_lo)
    q_lo  = max(|k1-k3|, |k2-k4|),   q_hi = min(k1+k3, k2+k4)
    kappa^2 = (e^2/(eps0 eps_b)) int rho(E) (-df/dE) dE           (Thomas-Fermi, recomputed from f)

``c_ee = 1`` is what the spin-summed direct term gives; ``c_ee = 3/16`` reproduces the prefactor of
[14] Eq. (15) exactly as printed (the two differ only by this number).  The kernel does not factorise,
so one evaluation costs O(N^3).  It is evaluated here with lookup tables of the antiderivative at
k_a + k_b and |k_a - k_b| stored along diagonals, which reduces the inner loop to a contiguous,
vectorisable min/max/multiply reduction.

B0 - "random-k" approximation with a constant matrix element  (ref. [40])
----------------------------------------------------------------------------------------------------
Momentum conservation is dropped and |M|^2 is a constant, leaving only energy conservation:

    (df/dt)(E1) = K0 * sum_{E2,E3} dE^2 rho2 rho3 rho4 [ ... ]

With a = rho f and b = rho (1-f) the double sum is a convolution followed by a correlation,

    gain(E1) = (1-f1) * sum_2 b_2 (a*a)_{1+2},      loss(E1) = f1 * sum_2 a_2 (b*b)_{1+2},

which FFTs evaluate in O(N log N).  ``K0`` is chosen so that B0 and B1 have the same Fermi-liquid
coefficient at T = 0, so a B0-B1 comparison tests the *shape* of the kernel, not its overall strength.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange
from scipy.signal import fftconvolve

from .. import constants as c
from ..grid import EnergyGrid
from ..material import Material


# ======================================================================================
# shared pieces
# ======================================================================================
def thomas_fermi_kappa(grid: EnergyGrid, mat: Material, f: np.ndarray) -> float:
    """Thomas-Fermi screening wavevector (1/nm) for an arbitrary isotropic f(E)."""
    s = grid.minus_dfdE_weighted(f, grid.rho_edges)
    return math.sqrt(max(c.E2_OVER_EPS0 / mat.eps_b * s, 1e-12))


def _antiderivative(q, kappa):
    return q * kappa / (q * q + kappa * kappa) + np.arctan(q / kappa)


def b1_prefactor(mat: Material, kappa: float) -> float:
    """P of the B1 integral, 1/(eV^2 fs nm)."""
    m0 = mat.band.m0
    return (
        mat.c_ee
        / (c.HBAR * (2.0 * math.pi) ** 3)
        * (c.E2_OVER_EPS0 / mat.eps_b) ** 2
        * (m0 / (c.HBAR**2 * kappa)) ** 3
    )


def fermi_liquid_K(mat: Material, kappa: float | None = None) -> float:
    """K in 1/tau_ee = K [(pi k_B T)^2 + (E - E_F)^2], 1/(eV^2 fs); cf. [14] Eqs. (17)-(18).

    Evaluated from the B1 kernel at T = 0: K = P/(2 k_F) * (1+2 C E_F)^3 * Phi(0 -> 2 k_F).
    """
    band = mat.band
    k_f = band.k_fermi(mat.n_e)
    e_f = band.fermi_energy(mat.n_e)
    if kappa is None:
        kappa = math.sqrt(c.E2_OVER_EPS0 / mat.eps_b * float(band.dos(e_f)))
    phi = float(_antiderivative(2.0 * k_f, kappa))
    return b1_prefactor(mat, kappa) / (2.0 * k_f) * (1.0 + 2.0 * mat.C * e_f) ** 3 * phi


# ======================================================================================
# B1: screened Coulomb, momentum resolved, O(N^3)
# ======================================================================================
@njit(parallel=True, fastmath=True, cache=True)
def _b1_tables(k, mfac, kappa):
    """Diagonal-ordered lookup tables so the innermost loop of the kernel runs over contiguous memory.

    For d = i2 - i4 (row index d + N - 1) and column i2:
        AP = A(k[i2] + k[i4]),  AM = A(|k[i2] - k[i4]|),  MM = mfac[i2] * mfac[i4],
    with A(q) = q kappa/(q^2 + kappa^2) + atan(q/kappa) the antiderivative of 2 kappa^3/(q^2+kappa^2)^2.
    """
    n = k.shape[0]
    ap = np.zeros((2 * n - 1, n))
    am = np.zeros((2 * n - 1, n))
    mm = np.zeros((2 * n - 1, n))
    k2 = kappa * kappa
    for r in prange(2 * n - 1):
        d = r - (n - 1)
        lo = d if d > 0 else 0
        hi = n + d if d < 0 else n
        for i2 in range(lo, hi):
            i4 = i2 - d
            qp = k[i2] + k[i4]
            qm = abs(k[i2] - k[i4])
            ap[r, i2] = qp * kappa / (qp * qp + k2) + math.atan(qp / kappa)
            am[r, i2] = qm * kappa / (qm * qm + k2) + math.atan(qm / kappa)
            mm[r, i2] = mfac[i2] * mfac[i4]
    return ap, am, mm


@njit(parallel=True, fastmath=True, cache=True)
def _b1_kernel(f, mfac, ap, am, mm, want_rate):
    """Return (S, R): S_i = sum m2 m3 m4 Phi [gain - loss],  R_i = sum m2 m3 m4 Phi [out + in coefficient]."""
    n = f.shape[0]
    out = np.zeros(n)
    rate = np.zeros(n)
    g = 1.0 - f
    for i1 in prange(n):
        f1 = f[i1]
        g1 = g[i1]
        acc = 0.0
        acc_r = 0.0
        for i3 in range(n):
            d = i3 - i1  # equals i2 - i4
            r = d + n - 1
            ap13 = ap[n - 1 - d, i1]  # pair (i1, i3): row index (i1 - i3) + n - 1
            am13 = am[n - 1 - d, i1]
            lo = d if d > 0 else 0
            hi = n + d if d < 0 else n
            s_in = 0.0  # sum_2 m2 m4 Phi (1-f2) f4
            s_out = 0.0  # sum_2 m2 m4 Phi f2 (1-f4)
            for i2 in range(lo, hi):
                phi = min(ap13, ap[r, i2]) - max(am13, am[r, i2])
                w = mm[r, i2] * max(phi, 0.0)
                s_in += w * g[i2] * f[i2 - d]
                s_out += w * f[i2] * g[i2 - d]
            m3 = mfac[i3]
            acc += m3 * (f[i3] * g1 * s_in - f1 * g[i3] * s_out)
            if want_rate:
                acc_r += m3 * (g[i3] * s_out + f[i3] * s_in)
        out[i1] = acc
        rate[i1] = acc_r
    return out, rate


@dataclass
class ElectronElectronB1:
    """Momentum-resolved screened-Coulomb collision integral (ground truth; target of the WP2 surrogate)."""

    grid: EnergyGrid
    mat: Material
    dynamic_screening_from_f: bool = True  # recompute kappa from f at every call (plan: "screening recomputed from f")
    kappa_fixed: float | None = None
    _cache: tuple | None = field(init=False, repr=False, default=None)

    def kappa(self, f) -> float:
        if not self.dynamic_screening_from_f:
            if self.kappa_fixed is None:
                e_f = self.mat.band.fermi_energy(self.mat.n_e)
                self.kappa_fixed = math.sqrt(c.E2_OVER_EPS0 / self.mat.eps_b * float(self.mat.band.dos(e_f)))
            return self.kappa_fixed
        return thomas_fermi_kappa(self.grid, self.mat, f)

    def _tables(self, kappa):
        if self._cache is not None and self._cache[0] == kappa:
            return self._cache[1]
        tabs = _b1_tables(self.grid.k, self.grid.mfac, kappa)
        self._cache = (kappa, tabs)
        return tabs

    def rhs(self, f: np.ndarray, with_rate: bool = False):
        """df/dt in 1/fs; optionally also the collision rate 1/tau_ee(E) of [14] Eq. (16)."""
        g = self.grid
        kappa = self.kappa(f)
        ap, am, mm = self._tables(kappa)
        s, r = _b1_kernel(np.ascontiguousarray(f, dtype=np.float64), g.mfac, ap, am, mm, with_rate)
        pref = b1_prefactor(self.mat, kappa) * g.dE**2 / g.k
        if with_rate:
            return pref * s, pref * r
        return pref * s

    def rate(self, f: np.ndarray) -> np.ndarray:
        return self.rhs(f, with_rate=True)[1]


# ======================================================================================
# B0: random-k, constant matrix element, O(N log N)
# ======================================================================================
@dataclass
class ElectronElectronB0:
    """Random-k / constant-matrix-element collision integral evaluated with FFT convolutions."""

    grid: EnergyGrid
    mat: Material
    K0: float | None = None  # (eV nm^3)^3 / (eV^2 fs); None -> matched to B1's Fermi-liquid coefficient
    dynamic_screening_from_f: bool = True  # rescale K0 with kappa[f] (one scalar per call, as in [40] Eq. (17))
    _rho_f: float = field(init=False, repr=False, default=0.0)

    def __post_init__(self):
        e_f = self.mat.band.fermi_energy(self.mat.n_e)
        self._rho_f = float(self.mat.band.dos(e_f))

    def strength(self, f=None) -> float:
        """K0 actually used: 1/tau = K0 rho_F^3 (E-E_F)^2 / 2 is made equal to K_FL(kappa) (E-E_F)^2."""
        if self.K0 is not None:
            return self.K0
        kappa = thomas_fermi_kappa(self.grid, self.mat, f) if (self.dynamic_screening_from_f and f is not None) else None
        return 2.0 * fermi_liquid_K(self.mat, kappa) / self._rho_f**3

    def _sums(self, f):
        g = self.grid
        n = g.N
        a = g.rho * f
        b = g.rho * (1.0 - f)
        aa = fftconvolve(a, a)  # (a*a)_m = sum_3 a_3 a_{m-3},  m = 0 .. 2N-2
        bb = fftconvolve(b, b)
        # c_1 = sum_2 x_2 y_{1+2}  ==  correlation; computed as convolution with reversed x
        gain_c = fftconvolve(b[::-1], aa)[n - 1 : 2 * n - 1]
        loss_c = fftconvolve(a[::-1], bb)[n - 1 : 2 * n - 1]
        return gain_c, loss_c

    def rhs(self, f: np.ndarray, with_rate: bool = False):
        g = self.grid
        gain_c, loss_c = self._sums(f)
        pref = self.strength(f) * g.dE**2
        d = pref * ((1.0 - f) * gain_c - f * loss_c)
        if with_rate:
            return d, pref * (gain_c + loss_c)
        return d

    def rate(self, f: np.ndarray) -> np.ndarray:
        return self.rhs(f, with_rate=True)[1]

    def rhs_direct(self, f: np.ndarray) -> np.ndarray:
        """O(N^3) direct summation of the same integral (reference for testing the FFT evaluation)."""
        g = self.grid
        a = g.rho * f
        b = g.rho * (1.0 - f)
        return self.strength(f) * g.dE**2 * _b0_direct(f, a, b)


@njit(parallel=True, cache=True)
def _b0_direct(f, a, b):
    n = f.shape[0]
    out = np.zeros(n)
    for i1 in prange(n):
        gain = 0.0
        loss = 0.0
        for i3 in range(n):
            d = i3 - i1
            lo = d if d > 0 else 0
            hi = n + d if d < 0 else n
            s_in = 0.0
            s_out = 0.0
            for i2 in range(lo, hi):
                i4 = i2 - d
                s_in += b[i2] * a[i4]
                s_out += a[i2] * b[i4]
            gain += a[i3] * s_in
            loss += b[i3] * s_out
        out[i1] = (1.0 - f[i1]) * gain - f[i1] * loss
    return out


# ======================================================================================
# relaxation-time approximation (baseline, refs. [36, 38, 40])
# ======================================================================================
@dataclass
class ElectronElectronRTA:
    """df/dt = -(f - f_FD[f]) / tau(E) with f_FD matched in density and energy and a Fermi-liquid tau(E).

    The plain energy-dependent RTA does not conserve number or energy; ``conserving=True`` projects the
    violation out (the "Lorentz term" device of [36]) so the baseline can be integrated over picoseconds.
    """

    grid: EnergyGrid
    mat: Material
    conserving: bool = True
    _guess: tuple | None = field(init=False, repr=False, default=None)

    def rhs(self, f: np.ndarray, with_rate: bool = False):
        from ..thermo import match_fermi_dirac

        g = self.grid
        mu, T, f_fd = match_fermi_dirac(g, f, guess=self._guess)
        self._guess = (mu, T)
        K = fermi_liquid_K(self.mat)
        inv_tau = K * ((math.pi * c.KB * T) ** 2 + (g.E - mu) ** 2)
        d = -(f - f_fd) * inv_tau
        if self.conserving:
            d = project_conserving(g, d, f_fd)
        return (d, inv_tau) if with_rate else d

    def rate(self, f):
        return self.rhs(f, with_rate=True)[1]


def project_conserving(grid: EnergyGrid, d: np.ndarray, f_ref: np.ndarray) -> np.ndarray:
    """Remove the number- and energy-violating part of a collision update ``d`` (1/fs).

    Subtracts a combination of the two basis functions h(1-h) and E h(1-h), h = f_ref, i.e. the directions
    in which a Fermi-Dirac function moves when mu and T change.  The result has zero number and energy
    rate.  The same projection is the hard constraint planned for the WP2 network output.
    """
    w = np.clip(f_ref * (1.0 - f_ref), 0.0, None)
    b1, b2 = w, grid.E * w
    r = grid.rho * grid.dE
    A = np.array(
        [
            [np.sum(r * b1), np.sum(r * b2)],
            [np.sum(r * grid.E * b1), np.sum(r * grid.E * b2)],
        ]
    )
    rhs = np.array([np.sum(r * d), np.sum(r * grid.E * d)])
    if abs(np.linalg.det(A)) < 1e-300:
        return d
    x = np.linalg.solve(A, rhs)
    return d - x[0] * b1 - x[1] * b2
