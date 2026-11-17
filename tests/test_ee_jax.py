"""The JAX version of the exact collision integral must reproduce the numba kernel and keep its invariants."""
import numpy as np
import pytest

jax = pytest.importorskip("jax")

from enzkin import thermo  # noqa: E402
from enzkin.collisions.ee import ElectronElectronB1  # noqa: E402
from enzkin.collisions.ee_jax import ElectronElectronB1Jax  # noqa: E402
from enzkin.grid import EnergyGrid  # noqa: E402


@pytest.fixture(scope="module")
def small(mat):
    g = EnergyGrid(mat.band, mat.E_max, 48)
    f0 = thermo.fd_of_density(g, mat.n_e, 300.0)
    rng = np.random.default_rng(3)
    sh = int(round(0.95 / g.dE))
    states = []
    for amp in (0.05, 0.2, 0.4):
        f = f0.copy()
        bump = amp * f0[:-sh] * (1.0 - f0[sh:])
        f[sh:] += bump
        f[:-sh] -= bump * g.rho[sh:] / g.rho[:-sh]
        states.append(np.clip(f * (1.0 + 0.05 * rng.standard_normal(g.N)), 1e-6, 1 - 1e-6))
    return g, np.array(states), ElectronElectronB1(g, mat), ElectronElectronB1Jax(g, mat, d_batch=8)


def test_matches_numba_kernel(small):
    g, states, ref, jx = small
    d_j, r_j = jx.rhs_and_rate(states)
    for i, f in enumerate(states):
        d_n, r_n = ref.rhs(f, with_rate=True)
        assert np.abs(np.asarray(d_j[i]) - d_n).max() < 1e-11 * np.abs(d_n).max()
        assert np.abs(np.asarray(r_j[i]) - r_n).max() < 1e-11 * np.abs(r_n).max()


def test_single_state_call_matches_batched(small):
    g, states, _, jx = small
    assert np.allclose(np.asarray(jx.rhs(states[1])), np.asarray(jx.rhs(states)[1]), rtol=1e-13, atol=0)


def test_conserves_number_and_energy(small):
    g, states, _, jx = small
    d = np.asarray(jx.rhs(states))
    for di in d:
        scale = np.sum(g.rho * np.abs(di)) * g.dE
        assert abs(g.number_rate(di)) < 1e-11 * scale
        assert abs(g.energy_rate(di)) < 1e-11 * scale


def test_fermi_dirac_is_fixed_point(small, mat):
    g, states, _, jx = small
    f_fd = thermo.fd_of_density(g, mat.n_e, 3000.0)
    assert np.abs(np.asarray(jx.rhs(f_fd))).max() < 1e-11 * np.abs(np.asarray(jx.rhs(states[1]))).max()


def test_differentiable(small):
    """Gradient of the entropy production with respect to f exists and is finite (needed by WP4)."""
    import jax.numpy as jnp

    g, states, _, jx = small
    w = jnp.asarray(g.rho * g.dE)

    def entropy_production(f):
        d, _ = jx._single(f)
        return -jnp.sum(w * jnp.log(f / (1.0 - f)) * d)

    f = jnp.asarray(states[1])
    val, grad = jax.value_and_grad(entropy_production)(f)
    assert float(val) > 0.0
    assert np.all(np.isfinite(np.asarray(grad)))
    # directional derivative agrees with a finite difference
    v = jnp.asarray(np.random.default_rng(0).standard_normal(g.N) * states[1] * (1 - states[1]))
    h = 1e-6
    fd = (entropy_production(f + h * v) - entropy_production(f - h * v)) / (2 * h)
    assert float(jnp.dot(grad, v)) == pytest.approx(float(fd), rel=1e-5)
