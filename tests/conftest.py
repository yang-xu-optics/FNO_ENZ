import numpy as np
import pytest

from enzkin import thermo
from enzkin.grid import EnergyGrid
from enzkin.material import ito_sarkar2023


@pytest.fixture(scope="session")
def mat():
    return ito_sarkar2023()


@pytest.fixture(scope="session")
def grid(mat):
    return EnergyGrid(mat.band, mat.E_max, 96)


@pytest.fixture(scope="session")
def f_cold(grid, mat):
    return thermo.fd_of_density(grid, mat.n_e, 300.0)


@pytest.fixture(scope="session")
def f_noneq(grid, mat, f_cold):
    """A strongly non-thermal state: cold Fermi sea plus one-photon shoulders and a random ripple."""
    rng = np.random.default_rng(7)
    sh = int(round(0.95 / grid.dE))
    f = f_cold.copy()
    bump = 0.2 * f_cold[:-sh] * (1.0 - f_cold[sh:])
    f[sh:] += bump
    f[:-sh] -= bump * grid.rho[sh:] / grid.rho[:-sh]
    f = np.clip(f * (1.0 + 0.05 * rng.standard_normal(grid.N)), 1e-6, 1.0 - 1e-6)
    return f
