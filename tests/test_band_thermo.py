"""Band structure, thermodynamics and optics against the numbers quoted in refs. [14, 15]."""
import numpy as np
import pytest

from enzkin import constants as c
from enzkin import optics, thermo
from enzkin.grid import EnergyGrid


def test_fermi_energy_and_wavevector(mat):
    band = mat.band
    assert band.fermi_energy(mat.n_e) == pytest.approx(0.88, abs=0.005)  # [14] Sec. II
    assert band.k_fermi(mat.n_e) == pytest.approx(3.54, abs=0.01)  # [14] Sec. II B, nm^-1
    assert band.fermi_energy(mat.n_e) / c.KB == pytest.approx(10204, rel=0.005)  # T_F in [15] Fig. 2


def test_debye_wavevector(mat):
    assert mat.E_debye / mat.hbar_v_sound == pytest.approx(14.3, abs=0.05)  # [14] Sec. II B, nm^-1


def test_dos_integrates_to_density(mat):
    g = EnergyGrid(mat.band, mat.E_max, 4096)
    e_f = mat.band.fermi_energy(mat.n_e)
    n = np.sum(g.rho[g.E < e_f]) * g.dE
    assert n == pytest.approx(mat.n_e, rel=2e-3)


def test_velocity_is_dE_dk(mat):
    band = mat.band
    E = np.linspace(0.05, 3.5, 50)
    h = 1e-6
    dk = band.k(E + h) - band.k(E - h)
    assert np.allclose(band.velocity(E), 2 * h / dk / c.HBAR, rtol=1e-6)


def test_sommerfeld_coefficient(grid, mat):
    gamma_e = thermo.heat_capacity(grid, mat.n_e, 300.0) / 300.0 / c.J_PER_M3_IN_EV_PER_NM3
    assert gamma_e == pytest.approx(12.7, rel=0.02)  # J m^-3 K^-2, [14] Sec. III B


def test_chemical_potential_drops_with_temperature(grid, mat):
    mus = [thermo.mu_of_density(grid, mat.n_e, T) for T in (300.0, 3000.0, 10000.0, 20000.0)]
    assert all(a > b for a, b in zip(mus, mus[1:]))
    assert mus[-1] < 0.0  # "optical de-doping": mu goes negative at very high T_e ([15] Fig. 2b)


@pytest.mark.parametrize("T", [300.0, 2500.0, 12000.0])
def test_match_fermi_dirac_roundtrip(grid, mat, T):
    f = thermo.fd_of_density(grid, mat.n_e, T)
    mu, T_fit, f_fit = thermo.match_fermi_dirac(grid, f)
    assert T_fit == pytest.approx(T, rel=1e-8)
    assert grid.density(f_fit) == pytest.approx(grid.density(f), rel=1e-11)
    assert grid.energy(f_fit) == pytest.approx(grid.energy(f), rel=1e-11)
    # warm-started Newton lands on the same answer
    mu2, T2, _ = thermo.match_fermi_dirac(grid, f, guess=(mu + 0.02, T * 1.05))
    assert T2 == pytest.approx(T, rel=1e-8)


def test_match_of_nonthermal_state_has_same_moments(grid, f_noneq):
    _, _, f_fd = thermo.match_fermi_dirac(grid, f_noneq)
    assert grid.density(f_fd) == pytest.approx(grid.density(f_noneq), rel=1e-10)
    assert grid.energy(f_fd) == pytest.approx(grid.energy(f_noneq), rel=1e-10)


def test_plasma_frequency_two_forms_agree(grid, mat, f_cold):
    a = optics.plasma_frequency_sq(grid, mat, f_cold)
    b = optics.plasma_frequency_sq_gradient_form(grid, mat, f_cold)
    assert a == pytest.approx(b, rel=2e-3)
    # zero-temperature closed form: omega_p^2 = e^2 n / (eps0 m0 (1 + 2 C E_F))
    e_f = mat.band.fermi_energy(mat.n_e)
    closed = c.E2_OVER_EPS0 * mat.n_e / (mat.band.m0 * (1 + 2 * mat.C * e_f))
    assert a == pytest.approx(closed, rel=5e-3)


def test_plasma_frequency_redshifts_when_heated(grid, mat):
    """Hot electrons sit higher in a non-parabolic band, are heavier, and lower omega_p (refs. [5-8])."""
    wp2 = [optics.plasma_frequency_sq(grid, mat, thermo.fd_of_density(grid, mat.n_e, T)) for T in (300.0, 3000.0, 8000.0)]
    assert wp2[0] > wp2[1] > wp2[2]


def test_enz_point_near_225_thz(mat):
    from enzkin.model import KineticModel

    model = KineticModel(mat, N=96)
    f0, T = model.split(model.initial_state(300.0))
    w_enz = optics.enz_frequency(model.grid, mat, f0, model.damping(f0, T))
    assert w_enz / (2 * np.pi) * 1e3 == pytest.approx(225.0, abs=3.0)  # THz, [15] Sec. III
