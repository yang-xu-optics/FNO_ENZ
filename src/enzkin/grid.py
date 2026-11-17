"""Uniform, cell-centred energy grid.

    E_i = (i + 1/2) dE,   i = 0 .. N-1,   dE = E_max / N

Cell centres avoid the k = 0 singularity at the band bottom, and a uniform spacing makes energy
conservation exact on the grid: E_i + E_j = E_k + E_l  <=>  i + j = k + l.  Every collision operator in
this package relies on that identity, so number and energy are conserved to round-off, not to
discretisation error.

Moments are midpoint sums:  n = sum_i rho_i f_i dE,   u = sum_i E_i rho_i f_i dE.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .band import KaneBand


@dataclass(frozen=True)
class EnergyGrid:
    band: KaneBand
    E_max: float = 4.0
    N: int = 512
    E: np.ndarray = field(init=False, repr=False)
    edges: np.ndarray = field(init=False, repr=False)
    rho: np.ndarray = field(init=False, repr=False)
    k: np.ndarray = field(init=False, repr=False)
    mfac: np.ndarray = field(init=False, repr=False)
    w_opt: np.ndarray = field(init=False, repr=False)
    rho_edges: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        dE = self.E_max / self.N
        E = (np.arange(self.N) + 0.5) * dE
        edges = np.arange(self.N + 1) * dE
        object.__setattr__(self, "E", E)
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "rho", self.band.dos(E))
        object.__setattr__(self, "k", self.band.k(E))
        object.__setattr__(self, "mfac", self.band.mass_factor(E))
        object.__setattr__(self, "w_opt", self.band.optical_weight(E))
        object.__setattr__(self, "rho_edges", self.band.dos(edges))

    @property
    def dE(self) -> float:
        return self.E_max / self.N

    # --- moments ------------------------------------------------------------------------
    def density(self, f) -> float:
        """Electron density, 1/nm^3."""
        return float(np.sum(self.rho * f) * self.dE)

    def energy(self, f) -> float:
        """Kinetic energy density, eV/nm^3."""
        return float(np.sum(self.E * self.rho * f) * self.dE)

    def number_rate(self, dfdt) -> float:
        return float(np.sum(self.rho * dfdt) * self.dE)

    def energy_rate(self, dfdt) -> float:
        return float(np.sum(self.E * self.rho * dfdt) * self.dE)

    def minus_dfdE_weighted(self, f, weight_edges) -> float:
        """int weight(E) (-df/dE) dE evaluated on cell edges (exact summation by parts)."""
        return float(np.sum(weight_edges[1:-1] * (f[:-1] - f[1:])))
