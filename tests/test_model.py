"""Whole-model checks: exact audits under driving, baselines, interband channels, transfer matrix."""
import numpy as np
import pytest

from enzkin import constants as c
from enzkin.collisions.ee import ElectronElectronB0, ElectronElectronB1
from enzkin.drive import FieldDrive, PowerDrive, Pulse, ThinFilmDrive
from enzkin.excitation import Excitation
from enzkin.integrate import run
from enzkin.interband import AugerLikeRelaxation, TwoPhotonAbsorption
from enzkin.model import KineticModel
from enzkin.tmm import absorptance_of_layer, reflection_transmission
from enzkin.ttm import ThermalModel


@pytest.mark.parametrize("mode", ["absorption_only", "net"])
@pytest.mark.parametrize("hw", [0.95, 0.777, 1.55])  # commensurate or not with the grid
def test_excitation_conserves_number_and_delivers_power(grid, f_noneq, mode, hw):
    exc = Excitation(grid, hw, mode)
    p_abs = 3.3e-4
    d = exc.rhs(f_noneq, p_abs)
    assert abs(grid.number_rate(d)) < 1e-12 * np.sum(grid.rho * np.abs(d)) * grid.dE
    assert grid.energy_rate(d) == pytest.approx(p_abs, rel=1e-12)


@pytest.mark.parametrize("ee", [ElectronElectronB0, ElectronElectronB1])
def test_driven_run_audits(mat, ee):
    model = KineticModel(mat, N=64, ee=ee, drive=FieldDrive(Pulse(1300.0, 60.0), 1.5e11))
    tr = run(model, -120.0, 300.0, 0.5, sample_every=20)
    assert np.abs(tr["n"] - tr["n"][0]).max() < 1e-10 * mat.n_e
    assert np.abs(tr["audit"] - tr["audit"][0]).max() < 1e-10 * tr["u_abs"][-1]
    assert tr["T_e"].max() > 1500 and tr.T_ph[-1] > 300.0
    assert tr["f_min"].min() > -1e-9 and tr["f_max"].max() < 1 + 1e-9
    # heated electrons in a non-parabolic band: plasma frequency drops, Re eps rises
    assert tr["wp2"].min() < 0.97 * tr["wp2"][0]
    assert tr["eps_re"].max() > tr["eps_re"][0] + 0.1


def test_adaptive_step_survives_stiff_regime(mat):
    """B0 at high excitation scatters at several per fs; a fixed 0.5 fs RK4 step would blow up."""
    model = KineticModel(mat, N=64, ee=ElectronElectronB0, drive=FieldDrive(Pulse(1300.0, 100.0), 2e11))
    tr = run(model, -200.0, 300.0, 0.5, sample_every=50)
    assert np.all(np.isfinite(tr.f))
    assert tr.meta["h_min"] < 0.25
    assert tr["f_min"].min() > -1e-9 and tr["f_max"].max() < 1 + 1e-9
    assert np.abs(tr["audit"] - tr["audit"][0]).max() < 1e-9 * tr["u_abs"][-1]


def test_thermal_model_matches_kinetics_for_long_pulse(mat):
    """For a pulse much longer than the e-e time the two-temperature limit is recovered ([15] Fig. 2)."""
    pulse = Pulse(1300.0, 220.0)
    out = []
    for cls in (KineticModel, ThermalModel):
        m = cls(mat, N=64, ee=ElectronElectronB1, drive=PowerDrive.from_absorbed_energy(pulse, 6e7))
        out.append(run(m, -450.0, 700.0, 0.5, sample_every=20))
    kin, th = out
    assert kin["T_e"].max() == pytest.approx(th["T_e"].max(), rel=0.03)
    assert kin.T_ph[-1] == pytest.approx(th.T_ph[-1], abs=1.0)
    assert kin["wp2"].min() == pytest.approx(th["wp2"].min(), rel=0.01)


def test_nonthermal_fraction_larger_for_short_pulse(mat):
    def peak_nonthermal(fwhm):
        m = KineticModel(mat, N=64, ee=ElectronElectronB1, drive=PowerDrive.from_absorbed_energy(Pulse(1300.0, fwhm), 6e7))
        return run(m, -2.5 * fwhm, 2.5 * fwhm, min(0.5, fwhm / 40), sample_every=10)["nonthermal"].max()

    assert peak_nonthermal(20.0) > 1.5 * peak_nonthermal(220.0)


# ------------------------------------------------------------------ interband
def test_two_photon_channel_closed_below_gap(grid, mat):
    assert not TwoPhotonAbsorption(grid, mat, Pulse(1300.0, 100.0).photon_energy, E_gap=2.67).open
    tpa = TwoPhotonAbsorption(grid, mat, Pulse(800.0, 10.0).photon_energy, E_gap=2.67)
    assert tpa.open and 0.0 < tpa.E_star < mat.band.fermi_energy(mat.n_e)  # final state is Pauli-blocked when cold


def test_two_photon_blocked_until_band_bottom_is_depleted(grid, mat, f_cold):
    tpa = TwoPhotonAbsorption(grid, mat, Pulse(800.0, 10.0).photon_energy, beta_cm_per_GW=5.0)
    cold = tpa.rhs(f_cold, 300.0, intensity=5e12)
    depleted = f_cold.copy()
    depleted[grid.E < 0.5] *= 0.5
    hot = tpa.rhs(depleted, 300.0, intensity=5e12)
    assert grid.number_rate(hot["df"]) > 1e6 * max(grid.number_rate(cold["df"]), 1e-300)
    # every pair costs two photons; the electron keeps E*, the rest is stored in the pair
    gen = grid.number_rate(hot["df"])
    assert hot["p_field"] == pytest.approx(2 * tpa.photon_energy * gen, rel=1e-12)
    assert hot["pair_power"] + grid.energy_rate(hot["df"]) == pytest.approx(hot["p_field"], rel=1e-10)


def test_auger_like_equilibrium_density(grid, mat):
    au = AugerLikeRelaxation(grid, mat)
    assert au.equilibrium_density(300.0) == pytest.approx(mat.n_e, rel=1e-6)
    assert au.equilibrium_density(20000.0) > 1.5 * mat.n_e  # hot electrons pull carriers across the gap


def test_interband_runs_keep_exact_energy_audit(mat):
    pulse = Pulse(1300.0, 136.0)
    model = KineticModel(mat, N=64, ee=ElectronElectronB0, drive=FieldDrive(pulse, 1.5e11))
    model.interband = [AugerLikeRelaxation(model.grid, mat)]
    tr = run(model, -300.0, 600.0, 0.5, sample_every=30)
    assert np.abs(tr["audit"] - tr["audit"][0]).max() < 1e-9 * tr["u_abs"][-1]
    assert tr["n"].max() > 1.02 * mat.n_e  # carrier density rises ...
    assert tr.t[tr["n"].argmax()] > 50.0  # ... and peaks after the pump, as in [22]
    assert tr["u_pair"].max() == pytest.approx(3.9 * (tr["n"].max() - mat.n_e), rel=1e-6)


# ------------------------------------------------------------------ transfer matrix
@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("angle", [0.0, 35.0, 60.0])
def test_tmm_energy_balance(pol, angle):
    layers = [(0.3 + 0.45j, 40.0)]
    R, T, _, _ = reflection_transmission(layers, 1300.0, angle, pol, eps_in=1.0, eps_out=2.25)
    A = absorptance_of_layer(layers, 0, 1300.0, angle, pol, eps_in=1.0, eps_out=2.25, n_sub=800)
    assert R + T + A == pytest.approx(1.0, abs=2e-5)


def test_tmm_fresnel_limit():
    import math

    R, _, _, _ = reflection_transmission([(2.25, 0.0)], 1300.0, 45.0, "p", 1.0, 2.25)
    th = math.radians(45.0)
    tt = math.asin(math.sin(th) / 1.5)
    rp = (1.5 * math.cos(th) - math.cos(tt)) / (1.5 * math.cos(th) + math.cos(tt))
    assert R == pytest.approx(rp**2, rel=1e-12)


def test_enz_film_absorbs_strongly_for_p_polarisation_only():
    """Berreman / ENZ-mode absorption: the thin film on a mirror absorbs p- but not s-polarised light."""
    layers = [(0.17 + 0.4j, 40.0), (-80 + 8j, 100.0)]
    a_p = absorptance_of_layer(layers, 0, 1303.0, 60.0, "p", eps_in=2.25, eps_out=1.0)
    a_s = absorptance_of_layer(layers, 0, 1303.0, 60.0, "s", eps_in=2.25, eps_out=1.0)
    assert a_p > 0.8 and a_s < 0.05


def test_thin_film_drive_absorptance_consistent_with_tmm(mat):
    pulse = Pulse(1303.0, 220.0)
    drive = ThinFilmDrive(pulse, 1e9, thickness=40.0, angle_deg=60.0, eps_in=2.25, layers_after=[(-80 + 8j, 100.0)])
    eps = 0.17 + 0.4j
    ref = absorptance_of_layer([(eps, 40.0), (-80 + 8j, 100.0)], 0, 1303.0, 60.0, "p", eps_in=2.25, eps_out=1.0, n_sub=16)
    assert drive.absorptance(eps) == pytest.approx(ref, rel=1e-9)
    # power density integrates to absorptance * incident flux / thickness
    flux = 1e9 * 1e4 * np.cos(np.radians(60.0))  # W/m^2 through the film plane
    p = drive.power(pulse.t0, eps) / c.W_PER_M3_IN_WORKING
    assert p * 40e-9 == pytest.approx(ref * flux, rel=1e-6)
