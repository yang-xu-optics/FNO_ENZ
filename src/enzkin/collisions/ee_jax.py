"""JAX implementation of the exact (B1) electron-electron collision integral.

The same integral and the same discrete form as the numba kernel in :mod:`enzkin.collisions.ee`
(screened Coulomb, momentum resolved, screening recomputed from f), written as array operations so that it

* runs on a GPU without code changes (install a CUDA build of JAX),
* is batched over spatial cells: ``f`` may have any leading dimensions, the last one is energy,
* is differentiable with respect to ``f`` and to the physical parameters (needed for gradient-based
  fitting in WP4).

Structure.  With d = i3 - i1 = i2 - i4 fixed, the kernel factor

    Phi = [ min(A(k1+k3), A(k2+k4)) - max(A(|k1-k3|), A(|k2-k4|)) ]_+

is an N x N matrix built from one vector over i1 and one over i2, and the sums over i2 are two
matrix-vector products.  Looping d over its 2N-1 values gives the O(N^3) total; ``d_batch`` of those
matrices are held in memory at a time, so the peak memory is  cells x d_batch x N^2 x 8 bytes.

The collision integral conserves particle number and energy to round-off only in double precision, so
this module switches JAX to 64-bit floats when imported.
"""
from __future__ import annotations

import math

import jax
import jax.numpy as jnp

from .. import constants as c
from ..grid import EnergyGrid
from ..material import Material

jax.config.update("jax_enable_x64", True)

_TWO_PI_CUBED = (2.0 * math.pi) ** 3


def b1_single(f, k, mfac, rho_edges, dE, c_ee, eps_b, m0, d_batch=None):
    """(df/dt, collision rate) for one distribution ``f`` of shape (N,).  All arguments may be traced."""
    n = f.shape[0]
    g = 1.0 - f
    kappa = jnp.sqrt(jnp.maximum(c.E2_OVER_EPS0 / eps_b * jnp.sum(rho_edges[1:-1] * (f[:-1] - f[1:])), 1e-12))

    def antiderivative(q):
        return q * kappa / (q * q + kappa * kappa) + jnp.arctan(q / kappa)

    idx = jnp.arange(n)

    def one_d(d):
        i3 = idx + d  # final state of electron 1, as a function of i1
        i4 = idx - d  # final state of electron 2, as a function of i2
        ok13 = (i3 >= 0) & (i3 < n)
        ok24 = (i4 >= 0) & (i4 < n)
        i3c = jnp.clip(i3, 0, n - 1)
        i4c = jnp.clip(i4, 0, n - 1)
        k3, k4 = k[i3c], k[i4c]
        ap13, am13 = antiderivative(k + k3), antiderivative(jnp.abs(k - k3))
        ap24, am24 = antiderivative(k + k4), antiderivative(jnp.abs(k - k4))
        phi = jnp.maximum(jnp.minimum(ap13[:, None], ap24[None, :]) - jnp.maximum(am13[:, None], am24[None, :]), 0.0)
        w24 = jnp.where(ok24, mfac * mfac[i4c], 0.0)
        s_in = phi @ (w24 * g * f[i4c])  # sum_2 m2 m4 Phi (1 - f2) f4
        s_out = phi @ (w24 * f * g[i4c])  # sum_2 m2 m4 Phi f2 (1 - f4)
        m3 = jnp.where(ok13, mfac[i3c], 0.0)
        return m3 * (f[i3c] * g * s_in - f * g[i3c] * s_out), m3 * (g[i3c] * s_out + f[i3c] * s_in)

    acc, rate = jax.lax.map(one_d, jnp.arange(-(n - 1), n), batch_size=d_batch)
    pref = c_ee / (c.HBAR * _TWO_PI_CUBED) * (c.E2_OVER_EPS0 / eps_b) ** 2 * (m0 / (c.HBAR**2 * kappa)) ** 3 * dE**2 / k
    return pref * jnp.sum(acc, axis=0), pref * jnp.sum(rate, axis=0)


class ElectronElectronB1Jax:
    """Drop-in counterpart of :class:`enzkin.collisions.ee.ElectronElectronB1` for batched use.

    ``rhs(f)`` and ``rate(f)`` accept arrays of shape (N,) or (cells, N) and return JAX arrays.
    ``d_batch`` trades memory for speed; ``None`` lets JAX process one d at a time.
    """

    def __init__(self, grid: EnergyGrid, mat: Material, d_batch: int | None = 16):
        self.grid, self.mat = grid, mat
        k, mfac, rho_edges = jnp.asarray(grid.k), jnp.asarray(grid.mfac), jnp.asarray(grid.rho_edges)
        dE, c_ee, eps_b, m0 = grid.dE, mat.c_ee, mat.eps_b, mat.band.m0

        def single(f):
            return b1_single(f, k, mfac, rho_edges, dE, c_ee, eps_b, m0, d_batch)

        self._single = single
        self._batched = jax.jit(jax.vmap(single))

    def rhs_and_rate(self, f):
        f = jnp.asarray(f, dtype=jnp.float64)
        if f.ndim == 1:
            d, r = self._batched(f[None, :])
            return d[0], r[0]
        return self._batched(f)

    def rhs(self, f, with_rate: bool = False):
        d, r = self.rhs_and_rate(f)
        return (d, r) if with_rate else d

    def rate(self, f):
        return self.rhs_and_rate(f)[1]
