"""Structural properties of the collision operators: conservation, fixed points, H-theorem, limits."""
import numpy as np
import pytest

from enzkin import constants as c
from enzkin import thermo
from enzkin.collisions.ee import (ElectronElectronB0, ElectronElectronB1, ElectronElectronRTA, fermi_liquid_K,
                                  project_conserving, thomas_fermi_kappa)
from enzkin.collisions.eph import ElectronPhonon

EE = [ElectronElectronB0, ElectronElectronB1]


def _scale(grid, d):
    return np.sum(grid.rho * np.abs(d)) * grid.dE


@pytest.mark.parametrize("cls", EE)
def test_ee_conserves_number_and_energy(cls, grid, mat, f_noneq):
    d = cls(grid, mat).rhs(f_noneq)
    s = _scale(grid, d)
    assert s > 0
    assert abs(grid.number_rate(d)) < 1e-11 * s
    assert abs(grid.energy_rate(d)) < 1e-11 * s


@pytest.mark.parametrize("cls", EE)
@pytest.mark.parametrize("T", [300.0, 4000.0])
def test_ee_fermi_dirac_is_fixed_point(cls, grid, mat, f_noneq, T):
    op = cls(grid, mat)
    f_fd = thermo.fd_of_density(grid, mat.n_e, T)
    assert np.abs(op.rhs(f_fd)).max() < 1e-11 * np.abs(op.rhs(f_noneq)).max()


@pytest.mark.parametrize("cls", EE)
def test_ee_h_theorem(cls, grid, mat, f_noneq):
    """Entropy production  -int rho ln(f/(1-f)) (df/dt) dE  is non-negative."""
    d = cls(grid, mat).rhs(f_noneq)
    sdot = -np.sum(grid.rho * np.log(f_noneq / (1 - f_noneq)) * d) * grid.dE
    assert sdot > 0


def test_b0_fft_equals_direct_sum(grid, mat, f_noneq):
    op = ElectronElectronB0(grid, mat)
    a, b = op.rhs(f_noneq), op.rhs_direct(f_noneq)
    assert np.abs(a - b).max() < 1e-11 * np.abs(b).max()


def test_thomas_fermi_wavevector(grid, mat, f_cold):
    # [16] footnote: Eq. (15) of [14] without the spurious 4 pi^3 gives q_TF = 3.83 nm^-1
    assert thomas_fermi_kappa(grid, mat, f_cold) == pytest.approx(3.83, abs=0.02)


def test_fermi_liquid_coefficient(mat):
    assert fermi_liquid_K(mat) == pytest.approx(0.065, rel=0.03)  # eV^-2 fs^-1, [14] Eq. (19) evaluated


@pytest.mark.parametrize("cls", EE)
def test_ee_rate_approaches_fermi_liquid_limit(cls, mat):
    from enzkin.grid import EnergyGrid

    g = EnergyGrid(mat.band, mat.E_max, 256)
    mu = thermo.mu_of_density(g, mat.n_e, 300.0)
    f = thermo.fermi(g.E, mu, 300.0)
    rate = cls(g, mat).rate(f)
    i = int(np.argmin(np.abs(g.E - (mu + 0.1))))
    fl = fermi_liquid_K(mat) * ((np.pi * c.KB * 300.0) ** 2 + (g.E[i] - mu) ** 2)
    assert rate[i] == pytest.approx(fl, rel=0.2)


def test_ee_relaxes_to_matched_fermi_dirac(grid, mat, f_noneq):
    """Free relaxation under B1 ends at the Fermi-Dirac function with the same number and energy."""
    op = ElectronElectronB1(grid, mat)
    _, _, f_target = thermo.match_fermi_dirac(grid, f_noneq)
    f = f_noneq.copy()
    dist0 = _scale(grid, f - f_target)
    t, dt = 0.0, 0.5
    while t < 3000.0:
        k1 = op.rhs(f)
        k2 = op.rhs(f + 0.5 * dt * k1)
        k3 = op.rhs(f + 0.5 * dt * k2)
        k4 = op.rhs(f + dt * k3)
        f = f + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    assert _scale(grid, f - f_target) < 2e-3 * dist0
    assert grid.density(f) == pytest.approx(grid.density(f_noneq), rel=1e-11)
    assert grid.energy(f) == pytest.approx(grid.energy(f_noneq), rel=1e-11)
    assert f.min() > -1e-12 and f.max() < 1 + 1e-12


def test_rta_projection_conserves(grid, mat, f_noneq):
    d = ElectronElectronRTA(grid, mat).rhs(f_noneq)
    s = _scale(grid, d)
    assert abs(grid.number_rate(d)) < 1e-10 * s and abs(grid.energy_rate(d)) < 1e-10 * s
    raw = ElectronElectronRTA(grid, mat, conserving=False).rhs(f_noneq)
    assert abs(grid.energy_rate(raw)) > 1e-4 * _scale(grid, raw)  # the plain tau(E) RTA does violate conservation


def test_project_conserving_removes_both_moments(grid, f_noneq):
    rng = np.random.default_rng(1)
    d = rng.standard_normal(grid.N) * f_noneq * (1 - f_noneq)
    p = project_conserving(grid, d, f_noneq)
    assert abs(grid.number_rate(p)) < 1e-12 * _scale(grid, d)
    assert abs(grid.energy_rate(p)) < 1e-12 * _scale(grid, d)


# ------------------------------------------------------------------ electron-phonon
def test_eph_conserves_number_and_reports_energy(grid, mat, f_noneq):
    d, power = ElectronPhonon(grid, mat).rhs(f_noneq, 300.0)
    assert abs(grid.number_rate(d)) < 1e-12 * _scale(grid, d)
    assert power == pytest.approx(-grid.energy_rate(d), rel=1e-10)  # energy lost by electrons = lattice gain


def test_eph_detailed_balance(grid, mat):
    op = ElectronPhonon(grid, mat)
    f = thermo.fd_of_density(grid, mat.n_e, 700.0)
    d, power = op.rhs(f, 700.0)
    d_hot, p_hot = op.rhs(thermo.fd_of_density(grid, mat.n_e, 2000.0), 700.0)
    assert np.abs(d).max() < 1e-12 * np.abs(d_hot).max()
    assert abs(power) < 1e-12 * abs(p_hot)
    assert p_hot > 0  # hot electrons heat the lattice
    assert op.rhs(thermo.fd_of_density(grid, mat.n_e, 300.0), 700.0)[1] < 0  # and a hot lattice heats electrons


def test_eph_coupling_matches_zero_temperature_formula(mat):
    """G(T -> 0) = D^2 m*^2 k_B k_F^4 (1 + 2 C E_F)^2 / (pi^3 rho hbar^3)   [14, Eq. (27)] ~ 2.5e16 W/m^3/K."""
    from enzkin.grid import EnergyGrid

    g = EnergyGrid(mat.band, mat.E_max, 512)
    k_f, e_f = mat.band.k_fermi(mat.n_e), mat.band.fermi_energy(mat.n_e)
    analytic = mat.kappa_eph * (1 + 2 * mat.C * e_f) ** 2 * c.KB * (2 * mat.hbar_v_sound * k_f) ** 4 / 4
    numeric = ElectronPhonon(g, mat).coupling(mat.n_e, 400.0, 300.0)
    assert analytic / c.W_PER_M3_IN_WORKING == pytest.approx(2.5e16, rel=0.02)
    assert numeric == pytest.approx(analytic, rel=0.05)


def test_eph_coupling_rises_with_electron_temperature_and_momentum_conservation_matters(grid, mat):
    op = ElectronPhonon(grid, mat)
    g_lo, g_hi = op.coupling(mat.n_e, 500.0, 300.0), op.coupling(mat.n_e, 5000.0, 300.0)
    assert g_hi > 2 * g_lo  # [14] Fig. 3(c)
    no_mc = ElectronPhonon(grid, mat, momentum_conservation=False).coupling(mat.n_e, 500.0, 300.0)
    assert no_mc > 10 * g_lo  # [14]: an order-of-magnitude overestimate without momentum conservation


def test_eph_rate_formula_is_minus_functional_derivative(grid, mat, f_noneq):
    op = ElectronPhonon(grid, mat)
    i = grid.N // 3
    h = 1e-7
    fp, fm = f_noneq.copy(), f_noneq.copy()
    fp[i] += h
    fm[i] -= h
    num = -(op.rhs(fp, 300.0)[0][i] - op.rhs(fm, 300.0)[0][i]) / (2 * h)
    assert op.rate(f_noneq, 300.0)[i] == pytest.approx(num, rel=1e-6)
