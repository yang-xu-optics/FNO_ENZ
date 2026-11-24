"""Density-matrix particle model (enzkin.densmat): structure, conservation and the published benchmark."""
import numpy as np
import pytest

from enzkin import constants as c
from enzkin.densmat import ParticleDM, ReducedResponse, box_dipole_weight, gaussian_pulse
from enzkin.material import ito_sarkar2023


@pytest.fixture(scope="module")
def dm():
    return ParticleDM(ito_sarkar2023(), N=96).calibrate_eps_inf(-1.6, 175.0)


def test_dipole_weight_is_symmetric_positive_and_matches_closed_form_off_diagonal():
    k, kp, L = 2.328, 3.541, 4.0
    w = float(box_dipole_weight(k, kp, L))
    assert w == pytest.approx(float(box_dipole_weight(kp, k, L)), rel=1e-10)
    # unregularised closed form: (64 pi / L) J / (k k' (k^2 - k'^2)^4), J = int_0^k x^2 sqrt(a + x^2) dx, a = k'^2 - k^2
    a = kp * kp - k * k
    J = k * (2 * k * k + a) * kp / 8 - a * a / 8 * np.log((k + kp) / np.sqrt(a))
    closed = 64 * np.pi / L * J / (k * kp * a**4)
    assert float(box_dipole_weight(k, kp, L, c_reg=1e-3)) == pytest.approx(closed, rel=1e-6)
    assert 0.9 * closed < w < closed  # the regularisation lowers it, strongly only near the diagonal


def test_oscillator_strength_sum_rule(dm):
    # a 4 nm box keeps most of the bulk oscillator strength, a larger one more (0.84, 0.91, 0.95 at 4, 8, 16 nm
    # on grids that resolve the kernel, whose width in energy shrinks as 1/L)
    assert 0.78 < dm.sum_rule() < 0.90
    bigger = ParticleDM(ito_sarkar2023(), N=192, L=8.0, n_T=4)
    assert dm.sum_rule() < bigger.sum_rule() < 1.0


def test_cold_permittivity_matches_un2026(dm):
    # [16] Fig. 2b: eps(175 THz) = -1.6 + 1.4i (real part fitted through eps_inf), zero crossing at 215 THz
    eps = dm.permittivity(c.frequency_thz_to_omega(175.0))
    assert eps.real == pytest.approx(-1.6, abs=1e-9)
    assert eps.imag == pytest.approx(1.4, rel=0.15)
    fs = np.linspace(190.0, 240.0, 51)
    zero = fs[np.argmin(np.abs(dm.permittivity(c.frequency_thz_to_omega(fs)).real))]
    assert zero == pytest.approx(215.0, abs=6.0)


@pytest.mark.parametrize("populations", ["kinetic", "thermal"])
def test_weak_pulse_reproduces_frequency_domain_susceptibility(dm, populations):
    r = dm.run(gaussian_pulse(1e-4, 5.0, 175.0), -40.0, 120.0, populations=populations)
    fr = np.fft.rfftfreq(r.t.size, r.t[1] - r.t[0]) * 2 * np.pi
    i = int(np.argmin(np.abs(fr - c.frequency_thz_to_omega(175.0))))
    chi_td = np.conj(np.fft.rfft(r.P)[i] / np.fft.rfft(r.E_loc)[i])
    chi_fd = dm.susceptibility(dm.f0, dm.eta_tab[0], fr[i])
    assert chi_td == pytest.approx(chi_fd, rel=2e-3)


def test_energy_and_number_bookkeeping_without_collisions(dm):
    r = dm.run(gaussian_pulse(1.2, 5.0, 175.0), -15.0, 25.0, collisions=False)
    assert r.meta["dN_max"] < 1e-11  # every photon-driven transfer moves electrons between two levels
    assert r.U[-1] == pytest.approx(r.meta["W_field"], rel=1e-10)  # electron energy = work done by the field
    assert r.U[-1] > 0


def test_absorbed_energy_matches_un2026_fig6(dm):
    # [16] Fig. 6: E0 = 1.2 V/nm, 5 fs, 175 THz: the energy given to the electrons settles at 0.63 eV/nm^3 and the
    # instantaneous power swings between about -0.14 and +0.33 eV/nm^3/fs
    r = dm.run(gaussian_pulse(1.2, 5.0, 175.0), -15.0, 25.0)
    power = r.E_loc * r.J / c.E2_OVER_EPS0
    assert r.meta["W_field"] == pytest.approx(0.63, rel=0.08)
    assert power.max() == pytest.approx(0.33, rel=0.12)
    assert power.min() == pytest.approx(-0.14, rel=0.15)


def test_thermal_closure_tracks_density_matrix(dm):
    p = gaussian_pulse(1.2, 5.0, 175.0)
    full, th = dm.run(p, -15.0, 25.0), dm.run(p, -15.0, 25.0, populations="thermal")
    assert np.abs(th.J - full.J).max() / np.abs(full.J).max() < 0.05
    assert th.meta["W_field"] == pytest.approx(full.meta["W_field"], rel=0.05)


def test_fitted_poles_reproduce_cold_permittivity_and_are_passive(dm):
    red = ReducedResponse.fitted(dm, 16)
    w = c.frequency_thz_to_omega(175.0)
    assert red.permittivity(w) == pytest.approx(dm.permittivity(w), rel=0.05)
    assert np.all(red.S_tab >= 0.0)
    drude = ReducedResponse.drude(dm, 175.0)
    assert drude.permittivity(w) == pytest.approx(dm.permittivity(w), rel=1e-6)
