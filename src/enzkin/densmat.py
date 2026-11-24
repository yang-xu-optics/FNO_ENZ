"""Energy-space density matrix of a small conducting particle under a few-cycle pulse (plan Tier 2b).

Model of refs. [16, 62] (Un, Sarkar & Sivan): populations f(E, t) and coherences between conduction-band
states of different energy, dipole elements of a one-dimensional box of width L, the collision integrals of
WP1 on the populations, and the quasi-static field of a sphere.  Derivation, the reconstruction of the dipole
strength (Eq. 16 of [16] is unusable as printed) and the validation are in docs/theory/tier2_density_matrix.md.

Oscillator form.  With an isotropic f(E) every pair of energies (E_lo < E_hi) is a damped oscillator,

    Q'' + 2 eta_p Q' + Delta_p^2 Q = S_p (f_lo - f_hi) E_loc(t),        P / eps0 = sum_p Q_p,

    Delta_p = (E_hi - E_lo)/hbar,   eta_p = (eta(E_lo) + eta(E_hi))/2,
    S_p = (e^2/eps0) (Delta_p/hbar) W(E_lo, E_hi) rho_lo rho_hi dE^2         [1/fs^2],

and the work E_loc dQ_p/dt done on it moves electrons from E_lo to E_hi (or back: stimulated emission):

    d(rho f dE)/dt |_hi = - d(rho f dE)/dt |_lo = E_loc Q_p' / ((e/eps0) (E_hi - E_lo)).

This is Eqs. (13)-(15) of [16] exactly (a = Re, b = Im of the dipole-weighted coherence; Q' = Delta b).
E_loc = [3 eps_h E_inc - P/eps0] / (eps_inf + 2 eps_h) for the sphere, Eq. (17).

Three levels share one oscillator kernel:

    kinetic   populations evolve (photon term + e-e + e-ph collisions)               the model of [16, 62]
    thermal   f = Fermi-Dirac at T_e(t), dU/dt = E_loc dP/dt - P_eph(T_e)            "thermal model" of [62]
    reduced   a few oscillators whose strength and damping are tabulated against T_e (ReducedResponse)

Units: eV, fs, nm; fields in V/nm; P/eps0 in V/nm.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from numba import njit

from . import constants as c
from . import optics, thermo
from .material import Material
from .model import KineticModel

E_OVER_EPS0 = c.E2_OVER_EPS0  # e/eps0 in V nm (numerically e^2/eps0 in eV nm)
_XG, _WG = np.polynomial.legendre.leggauss(96)


# ======================================================================================
# dipole strength of a box
# ======================================================================================
def box_dipole_weight(k, kp, L: float, c_reg: float = 1.0) -> np.ndarray:
    """W(k, k') = V |D(k, k', L)|^2 / e^2 in nm^5: angular average of the squared box dipole element.

    States are plane waves across the box and standing waves along the field, so a transition keeps the
    transverse momentum and changes k_z -> k_z' with k_z'^2 = k'^2 - k^2 + k_z^2.  In the quasi-continuum
    limit the squared element is (32 / pi L) k_z^2 k_z'^2 / (k_z^2 - k_z'^2)^4 dk_z'; its (k_z - k_z')^-4
    singularity stands for the discrete level spacing and is replaced by [(k_z - k_z')^2 + (c_reg/L)^2]^-2.
    c_reg = 1 keeps the oscillator-strength sum over odd level differences (sum 1/s^2 = pi^2/8).
    """
    k = np.asarray(k, dtype=float)[..., None]
    kp = np.asarray(kp, dtype=float)[..., None]
    kappa2 = (c_reg / L) ** 2
    lo = np.sqrt(np.maximum(k * k - kp * kp, 0.0))
    u = 0.5 * (_XG + 1.0)
    kz = lo + (k - lo) * u**2  # nodes cluster where k_z or k_z' vanishes
    jac = (k - lo) * 2.0 * u * 0.5 * _WG
    kzp = np.sqrt(np.maximum(kp * kp - k * k + kz * kz, 0.0))
    integrand = kz**2 * kzp / (((kz - kzp) ** 2 + kappa2) ** 2 * (kz + kzp) ** 4 + 1e-300)
    return 64.0 * math.pi / (L * k[..., 0] * kp[..., 0]) * np.sum(integrand * jac, axis=-1)


# ======================================================================================
# time stepping kernels
# ======================================================================================
@njit(cache=True, fastmath=True)
def _fermi(E, mu, kT, out):
    for i in range(E.size):
        x = (E[i] - mu) / kT
        out[i] = 1.0 / (math.exp(x) + 1.0) if x < 500.0 else 0.0


@njit(cache=True, fastmath=True)
def _interp_row(tab, pos, out):
    """Linear interpolation between rows of ``tab`` at fractional row index ``pos``."""
    j = int(pos)
    if j < 0:
        j = 0
    if j > tab.shape[0] - 2:
        j = tab.shape[0] - 2
    w = pos - j
    for i in range(tab.shape[1]):
        out[i] = (1.0 - w) * tab[j, i] + w * tab[j + 1, i]


@njit(cache=True, fastmath=True)
def _temperature(U, U_tab):
    """Fractional index in the temperature table for electron energy U (U_tab is increasing)."""
    n = U_tab.size
    if U <= U_tab[0]:
        return 0.0
    if U >= U_tab[n - 1]:
        return float(n - 1)
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if U_tab[mid] <= U:
            lo = mid
        else:
            hi = mid
    return lo + (U - U_tab[lo]) / (U_tab[lo + 1] - U_tab[lo])


@njit(cache=True, fastmath=True)
def _deriv_pairs(Q, R, f, U, Einc, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc, thermal,
                 E, U_tab, mu_tab, T_tab, eta_tab, peph_tab, dQ, dR, df, wf, we):
    """Right-hand side for the pair-resolved models.  Returns (dU/dt, work done by the field per time)."""
    P = 0.0
    Rsum = 0.0
    for p in range(Q.size):
        P += Q[p]
        Rsum += R[p]
    Eloc = a_loc * Einc - b_loc * P
    work = Eloc * Rsum / E_OVER_EPS0
    dU = work
    if thermal:
        pos = _temperature(U, U_tab)
        j = int(pos)
        if j > T_tab.size - 2:
            j = T_tab.size - 2
        w = pos - j
        T = (1.0 - w) * T_tab[j] + w * T_tab[j + 1]
        mu = (1.0 - w) * mu_tab[j] + w * mu_tab[j + 1]
        _fermi(E, mu, c.KB * T, wf)
        _interp_row(eta_tab, pos, we)
        dU -= (1.0 - w) * peph_tab[j] + w * peph_tab[j + 1]
        fuse, euse = wf, we
    else:
        fuse, euse = f, eta
        for i in range(df.size):
            df[i] = dfc[i]
    for p in range(Q.size):
        i, j2 = hi[p], lo[p]
        dQ[p] = R[p]
        dR[p] = -(euse[i] + euse[j2]) * R[p] - D2[p] * Q[p] + S[p] * (fuse[j2] - fuse[i]) * Eloc
        if not thermal:
            n = Eloc * R[p]
            df[i] += n * xhi[p]
            df[j2] -= n * xlo[p]
    return dU, work


@njit(cache=True, fastmath=True)
def _advance_pairs(Q, R, f, U, Einc_half, n0, nsteps, dt, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc,
                   thermal, E, U_tab, mu_tab, T_tab, eta_tab, peph_tab):
    """RK4 over ``nsteps`` steps.  ``Einc_half`` holds the incident field on the half-step grid."""
    npair, N = Q.size, f.size
    k1q, k2q, k3q, k4q = np.empty(npair), np.empty(npair), np.empty(npair), np.empty(npair)
    k1r, k2r, k3r, k4r = np.empty(npair), np.empty(npair), np.empty(npair), np.empty(npair)
    k1f, k2f, k3f, k4f = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    tq, tr, tf, wf, we = np.empty(npair), np.empty(npair), np.empty(N), np.empty(N), np.empty(N)
    W_abs = 0.0
    for s in range(nsteps):
        e0, e1, e2 = Einc_half[2 * (n0 + s)], Einc_half[2 * (n0 + s) + 1], Einc_half[2 * (n0 + s) + 2]
        u1, w1 = _deriv_pairs(Q, R, f, U, e0, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc, thermal,
                             E, U_tab, mu_tab, T_tab, eta_tab, peph_tab, k1q, k1r, k1f, wf, we)
        tq[:] = Q + 0.5 * dt * k1q
        tr[:] = R + 0.5 * dt * k1r
        tf[:] = f + 0.5 * dt * k1f
        u2, w2 = _deriv_pairs(tq, tr, tf, U + 0.5 * dt * u1, e1, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc,
                             thermal, E, U_tab, mu_tab, T_tab, eta_tab, peph_tab, k2q, k2r, k2f, wf, we)
        tq[:] = Q + 0.5 * dt * k2q
        tr[:] = R + 0.5 * dt * k2r
        tf[:] = f + 0.5 * dt * k2f
        u3, w3 = _deriv_pairs(tq, tr, tf, U + 0.5 * dt * u2, e1, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc,
                             thermal, E, U_tab, mu_tab, T_tab, eta_tab, peph_tab, k3q, k3r, k3f, wf, we)
        tq[:] = Q + dt * k3q
        tr[:] = R + dt * k3r
        tf[:] = f + dt * k3f
        u4, w4 = _deriv_pairs(tq, tr, tf, U + dt * u3, e2, a_loc, b_loc, lo, hi, D2, S, xhi, xlo, eta, dfc,
                             thermal, E, U_tab, mu_tab, T_tab, eta_tab, peph_tab, k4q, k4r, k4f, wf, we)
        Q += dt / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)
        R += dt / 6.0 * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
        if not thermal:
            f += dt / 6.0 * (k1f + 2.0 * k2f + 2.0 * k3f + k4f)
        U += dt / 6.0 * (u1 + 2.0 * u2 + 2.0 * u3 + u4)
        W_abs += dt / 6.0 * (w1 + 2.0 * w2 + 2.0 * w3 + w4)
    return U, W_abs


@njit(cache=True, fastmath=True)
def _deriv_red(Q, R, U, Einc, a_loc, b_loc, D2, U_tab, S_tab, eta_tab, peph_tab, dQ, dR, ws, we):
    P = 0.0
    Rsum = 0.0
    for m in range(Q.size):
        P += Q[m]
        Rsum += R[m]
    Eloc = a_loc * Einc - b_loc * P
    pos = _temperature(U, U_tab)
    _interp_row(S_tab, pos, ws)
    _interp_row(eta_tab, pos, we)
    j = int(pos)
    if j > U_tab.size - 2:
        j = U_tab.size - 2
    w = pos - j
    for m in range(Q.size):
        dQ[m] = R[m]
        dR[m] = -2.0 * we[m] * R[m] - D2[m] * Q[m] + ws[m] * Eloc
    return Eloc * Rsum / E_OVER_EPS0 - ((1.0 - w) * peph_tab[j] + w * peph_tab[j + 1])


@njit(cache=True, fastmath=True)
def _advance_red(Q, R, U, Einc_half, n0, nsteps, dt, a_loc, b_loc, D2, U_tab, S_tab, eta_tab, peph_tab):
    M = Q.size
    k1q, k2q, k3q, k4q = np.empty(M), np.empty(M), np.empty(M), np.empty(M)
    k1r, k2r, k3r, k4r = np.empty(M), np.empty(M), np.empty(M), np.empty(M)
    ws, we = np.empty(M), np.empty(M)
    for s in range(nsteps):
        e0, e1, e2 = Einc_half[2 * (n0 + s)], Einc_half[2 * (n0 + s) + 1], Einc_half[2 * (n0 + s) + 2]
        u1 = _deriv_red(Q, R, U, e0, a_loc, b_loc, D2, U_tab, S_tab, eta_tab, peph_tab, k1q, k1r, ws, we)
        u2 = _deriv_red(Q + 0.5 * dt * k1q, R + 0.5 * dt * k1r, U + 0.5 * dt * u1, e1, a_loc, b_loc, D2, U_tab,
                        S_tab, eta_tab, peph_tab, k2q, k2r, ws, we)
        u3 = _deriv_red(Q + 0.5 * dt * k2q, R + 0.5 * dt * k2r, U + 0.5 * dt * u2, e1, a_loc, b_loc, D2, U_tab,
                        S_tab, eta_tab, peph_tab, k3q, k3r, ws, we)
        u4 = _deriv_red(Q + dt * k3q, R + dt * k3r, U + dt * u3, e2, a_loc, b_loc, D2, U_tab, S_tab, eta_tab,
                        peph_tab, k4q, k4r, ws, we)
        Q += dt / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)
        R += dt / 6.0 * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
        U += dt / 6.0 * (u1 + 2.0 * u2 + 2.0 * u3 + u4)
    return U


# ======================================================================================
# results
# ======================================================================================
@dataclass
class Response:
    """Sampled response of one run (fields in V/nm, P means P/eps0, U in eV/nm^3 above the cold state)."""

    t: np.ndarray
    E_inc: np.ndarray
    E_loc: np.ndarray
    P: np.ndarray
    J: np.ndarray  # d(P/eps0)/dt
    U: np.ndarray
    T_e: np.ndarray
    meta: dict = field(default_factory=dict)
    f: np.ndarray | None = None


def gaussian_pulse(E0: float, fwhm: float, f_thz: float, cep: float = 0.0):
    """E_inc(t) = E0 exp(-2 ln2 (t/fwhm)^2) cos(omega0 t + cep): the pulse of refs. [16, 62] (fwhm in intensity)."""
    w0 = c.frequency_thz_to_omega(f_thz)
    return lambda t: E0 * np.exp(-2.0 * math.log(2.0) * (np.asarray(t) / fwhm) ** 2) * np.cos(w0 * np.asarray(t) + cep)


# ======================================================================================
# the particle
# ======================================================================================
class ParticleDM:
    """Density-matrix model of a small sphere (diameter L) of the material ``mat`` in a host ``eps_h``."""

    def __init__(self, mat: Material, N: int = 192, L: float = 4.0, eps_h: float = 1.0, c_reg: float = 1.0,
                 T0: float = 300.0, n_T: int = 140, T_max: float = 4.0e4):
        self.mat, self.L, self.eps_h, self.c_reg, self.T0 = mat, L, eps_h, c_reg, T0
        self.km = KineticModel(mat, N=N)
        g = self.grid = self.km.grid
        lo, hi = np.triu_indices(N, 1)
        self.lo, self.hi = lo.astype(np.int64), hi.astype(np.int64)
        dEp = g.E[hi] - g.E[lo]
        W = np.concatenate([box_dipole_weight(g.k[hi[s]], g.k[lo[s]], L, c_reg)
                            for s in np.array_split(np.arange(lo.size), max(1, lo.size // 4000))])
        self.Delta = dEp / c.HBAR
        self.S = E_OVER_EPS0 * self.Delta / c.HBAR * W * g.rho[hi] * g.rho[lo] * g.dE**2
        self.xhi = 1.0 / (E_OVER_EPS0 * dEp * g.rho[hi] * g.dE)
        self.xlo = 1.0 / (E_OVER_EPS0 * dEp * g.rho[lo] * g.dE)
        self.f0 = thermo.fd_of_density(g, mat.n_e, T0)
        self.U0 = g.energy(self.f0)
        self._tables(n_T, T_max)

    # ---- tables against electron temperature (thermal and reduced levels, diagnostics) ----------
    def _tables(self, n_T, T_max):
        g, mat = self.grid, self.mat
        self.T_tab = np.geomspace(self.T0, T_max, n_T)
        self.mu_tab = np.empty(n_T)
        self.U_tab = np.empty(n_T)
        self.peph_tab = np.empty(n_T)
        self.eta_tab = np.empty((n_T, g.N))
        self.f_tab = np.empty((n_T, g.N))
        for i, T in enumerate(self.T_tab):
            mu = thermo.mu_of_density(g, mat.n_e, T)
            f = thermo.fermi(g.E, mu, T)
            self.mu_tab[i], self.f_tab[i] = mu, f
            self.U_tab[i] = g.energy(f) - self.U0
            self.eta_tab[i] = self.km.damping(f, self.T0)
            self.peph_tab[i] = self.km.eph.rhs(f, self.T0)[1]
        self.U_tab[0], self.peph_tab[0] = 0.0, 0.0

    def temperature_of(self, U):
        return np.interp(U, self.U_tab, self.T_tab)

    # ---- linear response -----------------------------------------------------------------------
    def susceptibility(self, f, eta, omega):
        """chi(omega) of the conduction electrons for populations f and damping eta(E); eps = eps_inf + chi."""
        omega = np.atleast_1d(np.asarray(omega, dtype=float))
        num = self.S * (f[self.lo] - f[self.hi])
        etap = 0.5 * (eta[self.lo] + eta[self.hi])
        out = np.array([np.sum(num / (self.Delta**2 - w * w - 2j * etap * w)) for w in omega])
        return out if out.size > 1 else complex(out[0])

    def permittivity(self, omega, T: float | None = None, f=None, eta=None):
        if f is None:
            T = self.T0 if T is None else T
            f = thermo.fd_of_density(self.grid, self.mat.n_e, T)
        if eta is None:
            eta = self.km.damping(np.clip(f, 0.0, 1.0), self.T0)
        return self.mat.eps_inf + self.susceptibility(f, eta, omega)

    def calibrate_eps_inf(self, target_re: float, f_thz: float) -> "ParticleDM":
        """Set eps_inf so that Re eps(f_thz) of the cold particle equals ``target_re``."""
        chi = self.susceptibility(self.f0, self.eta_tab[0], c.frequency_thz_to_omega(f_thz))
        self.mat = self.mat.with_(eps_inf=float(target_re - chi.real))
        return self

    def sum_rule(self):
        """Total oscillator strength of the cold particle over the bulk plasma frequency squared."""
        return float(np.sum(self.S * (self.f0[self.lo] - self.f0[self.hi])) /
                     optics.plasma_frequency_sq(self.grid, self.mat, self.f0))

    # ---- time domain ---------------------------------------------------------------------------
    def _field_factors(self, eps_inf=None):
        e = self.mat.eps_inf if eps_inf is None else eps_inf
        return 3.0 * self.eps_h / (e + 2.0 * self.eps_h), 1.0 / (e + 2.0 * self.eps_h)

    def run(self, e_inc, t0: float, t1: float, dt: float = 0.02, block: int = 5, populations: str = "kinetic",
            collisions: bool = True, keep_f: bool = False, probe_thz: float | None = None) -> Response:
        """Integrate from t0 to t1.  ``populations`` is "kinetic" or "thermal"; samples every ``block`` steps,
        which is also how often the collision term and the damping are refreshed in the kinetic model."""
        g, km = self.grid, self.km
        thermal = populations == "thermal"
        nblocks = int(math.ceil((t1 - t0) / (dt * block)))
        nsteps = nblocks * block
        th = t0 + 0.5 * dt * np.arange(2 * nsteps + 1)
        Einc_half = np.asarray(e_inc(th), dtype=float)
        a_loc, b_loc = self._field_factors()
        Q, R = np.zeros(self.lo.size), np.zeros(self.lo.size)
        f = self.f0.copy()
        U = 0.0
        D2 = self.Delta**2
        eta = self.eta_tab[0].copy()
        dfc = np.zeros(g.N)
        out = {k: np.empty(nblocks + 1) for k in ("t", "E_inc", "E_loc", "P", "J", "U", "T_e", "f_top", "dN")}
        if probe_thz is not None:
            out["eps"] = np.empty(nblocks + 1, dtype=complex)
            w_pr = c.frequency_thz_to_omega(probe_thz)
        fs = [] if keep_f else None
        n_coll = 0
        wall = time.perf_counter()

        def sample(i, U):
            P, J = float(Q.sum()), float(R.sum())
            out["t"][i] = t0 + i * block * dt
            out["E_inc"][i] = Einc_half[2 * i * block]
            out["E_loc"][i] = a_loc * out["E_inc"][i] - b_loc * P
            out["P"][i], out["J"][i] = P, J
            Ue = U if thermal else g.energy(f) - self.U0
            out["U"][i] = Ue
            out["T_e"][i] = self.temperature_of(Ue)
            out["f_top"][i] = f[-1] if not thermal else 0.0
            out["dN"][i] = 0.0 if thermal else g.density(f) / self.mat.n_e - 1.0
            if probe_thz is not None:
                if thermal:
                    pos = np.interp(Ue, self.U_tab, np.arange(self.U_tab.size))
                    j = min(int(pos), self.U_tab.size - 2)
                    w = pos - j
                    fx = thermo.fermi(g.E, (1 - w) * self.mu_tab[j] + w * self.mu_tab[j + 1],
                                      (1 - w) * self.T_tab[j] + w * self.T_tab[j + 1])
                    ex = (1 - w) * self.eta_tab[j] + w * self.eta_tab[j + 1]
                else:
                    fx, ex = f, eta
                out["eps"][i] = self.mat.eps_inf + self.susceptibility(fx, ex, w_pr)
            if keep_f:
                fs.append(f.copy())

        W_abs = 0.0
        for i in range(nblocks):
            if not thermal and collisions:
                fc = np.clip(f, 0.0, 1.0)
                d_ee, r_ee = km.ee.rhs(fc, with_rate=True)
                d_ph, _ = km.eph.rhs(fc, self.T0)
                dfc = d_ee + d_ph
                eta = km.damping(fc, self.T0, r_ee)
                n_coll += 1
            sample(i, U)
            U, w = _advance_pairs(Q, R, f, U, Einc_half, i * block, block, dt, a_loc, b_loc, self.lo, self.hi, D2,
                                  self.S, self.xhi, self.xlo, eta, dfc, thermal, g.E, self.U_tab, self.mu_tab,
                                  self.T_tab, self.eta_tab, self.peph_tab)
            W_abs += w
        sample(nblocks, U)
        meta = dict(populations=populations, wall=time.perf_counter() - wall, n_steps=nsteps, n_coll=n_coll,
                    n_osc=int(self.lo.size), W_field=float(W_abs), f_top_max=float(out["f_top"].max()),
                    dN_max=float(np.abs(out["dN"]).max()))
        r = Response(out["t"], out["E_inc"], out["E_loc"], out["P"], out["J"], out["U"], out["T_e"], meta,
                     np.array(fs) if keep_f else None)
        if probe_thz is not None:
            r.meta["eps"] = out["eps"]
        return r


# ======================================================================================
# reduced response: a few oscillators with temperature-dependent strength and damping
# ======================================================================================
class ReducedResponse:
    """A few-pole constitutive law with one state variable (the electron energy / temperature).

        Q_m'' + 2 eta_m(T_e) Q_m' + Delta_m^2 Q_m = S_m(T_e) E_loc,   P/eps0 = sum_m Q_m,
        dU/dt = E_loc dP/dt - P_eph(T_e),   T_e = T(U).

    ``fitted(dm, K)`` fits K fixed poles with temperature-dependent strengths to the particle's susceptibility.
    ``drude(dm, f_thz)`` is a single pole at zero frequency with the bulk plasma frequency of the heated
    Fermi-Dirac distribution and a size-broadened damping fitted to the particle's cold permittivity at
    one frequency: the Tier 2a closure of the plan.
    """

    def __init__(self, dm: ParticleDM, Delta, S_tab, eta_tab, eps_inf, label):
        self.dm, self.Delta, self.S_tab, self.eta_tab, self.eps_inf, self.label = dm, Delta, S_tab, eta_tab, eps_inf, label

    @classmethod
    def fitted(cls, dm: ParticleDM, K: int, f_lo: float = 20.0, f_hi: float = 1000.0, n_w: int = 160, width: float = 2.0,
               collisional: bool = True):
        """K fixed, overlapping Lorentz poles (plus one Drude pole) with non-negative strengths S_k(T_e) fitted to
        the particle's susceptibility chi(omega; T_e) between f_lo and f_hi (THz) at every temperature of the table.

        The poles are geometrically spaced in transition frequency and as wide as their spacing, so that together
        they can represent the continuum of pair oscillators (resonant, Landau-type absorption included).
        Non-negative strengths keep every pole passive.  With ``collisional`` the damping of every pole also
        rises with temperature by the strength-weighted mean change of the collision rate, eta(T_e) - eta(T_0):
        the pole width is then inhomogeneous (fixed) plus homogeneous (collisions), as in the pair model.
        """
        from scipy.optimize import nnls

        w = c.frequency_thz_to_omega(np.geomspace(f_lo, f_hi, n_w))
        Dk = np.geomspace(0.5 * w[0], 1.5 * w[-1], K)
        ratio = Dk[1] / Dk[0]
        gk = 0.5 * width * Dk * (ratio - 1.0 / ratio) / 2.0  # half the spacing to the neighbours
        Dk = np.concatenate([[0.0], Dk])
        gk = np.concatenate([[float(np.mean(dm.eta_tab[0]))], gk])
        S_tab = np.empty((dm.T_tab.size, K + 1))
        g_tab = np.empty((dm.T_tab.size, K + 1))
        resid = np.empty(dm.T_tab.size)
        num_all = dm.S[None, :] * (dm.f_tab[:, dm.lo] - dm.f_tab[:, dm.hi])
        D2 = dm.Delta**2
        eta_mean0 = None
        for i in range(dm.T_tab.size):
            etap = 0.5 * (dm.eta_tab[i, dm.lo] + dm.eta_tab[i, dm.hi])
            eta_mean = float(np.sum(num_all[i] * etap) / np.sum(num_all[i]))
            eta_mean0 = eta_mean if eta_mean0 is None else eta_mean0
            g_tab[i] = gk + (eta_mean - eta_mean0 if collisional else 0.0)
            basis = 1.0 / (Dk[None, :] ** 2 - w[:, None] ** 2 - 2j * g_tab[i][None, :] * w[:, None])
            A = np.vstack([basis.real, basis.imag])
            chi = np.array([np.sum(num_all[i] / (D2 - x * x - 2j * etap * x)) for x in w])
            b = np.concatenate([chi.real, chi.imag])
            scale = np.concatenate([1.0 / np.abs(chi), 1.0 / np.abs(chi)])  # relative error
            S_tab[i], _ = nnls(A * scale[:, None], b * scale)
            resid[i] = float(np.max(np.abs(basis @ S_tab[i] - chi) / np.abs(chi)))
        out = cls(dm, Dk, S_tab, g_tab, dm.mat.eps_inf, f"{K + 1} fitted poles")
        out.fit_residual = float(resid.max())
        return out

    @classmethod
    def drude(cls, dm: ParticleDM, f_thz: float):
        g, mat = dm.grid, dm.mat
        w = c.frequency_thz_to_omega(f_thz)
        eps_dm = dm.permittivity(w)
        wp2 = np.array([optics.plasma_frequency_sq(g, mat, f) for f in dm.f_tab])
        gam_bulk = np.array([optics.drude_equivalent(g, mat, f, e)[1] for f, e in zip(dm.f_tab, dm.eta_tab)])
        # eps = eps_inf' - wp2 / (w (w + i Gamma)) matched to the particle at T0:  Gamma from the phase of chi
        # solve Im: wp2 G / (w (w^2 + G^2)) = eps''  for the smaller root
        a, b, cc = eps_dm.imag * w, -wp2[0], eps_dm.imag * w**3
        gam0 = (-b - math.sqrt(b * b - 4 * a * cc)) / (2 * a)
        eps_inf = eps_dm.real + wp2[0] / (w * w + gam0 * gam0)
        gam = gam0 + (gam_bulk - gam_bulk[0])
        return cls(dm, np.zeros(1), wp2[:, None], 0.5 * gam[:, None], float(eps_inf), "Drude (Tier 2a)")

    def permittivity(self, omega, T: float | None = None):
        pos = 0.0 if T is None else float(np.interp(T, self.dm.T_tab, np.arange(self.dm.T_tab.size)))
        j = min(int(pos), self.dm.T_tab.size - 2)
        w = pos - j
        S = (1 - w) * self.S_tab[j] + w * self.S_tab[j + 1]
        eta = (1 - w) * self.eta_tab[j] + w * self.eta_tab[j + 1]
        return self.eps_inf + complex(np.sum(S / (self.Delta**2 - omega**2 - 2j * eta * omega)))

    def run(self, e_inc, t0, t1, dt: float = 0.02, block: int = 5) -> Response:
        dm = self.dm
        nblocks = int(math.ceil((t1 - t0) / (dt * block)))
        th = t0 + 0.5 * dt * np.arange(2 * nblocks * block + 1)
        Einc_half = np.asarray(e_inc(th), dtype=float)
        a_loc, b_loc = dm._field_factors(self.eps_inf)
        M = self.Delta.size
        Q, R, U = np.zeros(M), np.zeros(M), 0.0
        D2 = self.Delta**2
        out = {k: np.empty(nblocks + 1) for k in ("t", "E_inc", "E_loc", "P", "J", "U")}
        wall = time.perf_counter()
        for i in range(nblocks + 1):
            out["t"][i] = t0 + i * block * dt
            out["E_inc"][i] = Einc_half[2 * i * block]
            out["P"][i], out["J"][i], out["U"][i] = Q.sum(), R.sum(), U
            out["E_loc"][i] = a_loc * out["E_inc"][i] - b_loc * out["P"][i]
            if i < nblocks:
                U = _advance_red(Q, R, U, Einc_half, i * block, block, dt, a_loc, b_loc, D2, dm.U_tab, self.S_tab,
                                 self.eta_tab, dm.peph_tab)
        meta = dict(populations=self.label, wall=time.perf_counter() - wall, n_osc=M)
        return Response(out["t"], out["E_inc"], out["E_loc"], out["P"], out["J"], out["U"], dm.temperature_of(out["U"]), meta)
