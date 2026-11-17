"""Time integration.

Classical fourth-order Runge-Kutta with a rate-limited step:

* every collision operator here has *linear* invariants (number, energy), and any Runge-Kutta method
  preserves linear invariants exactly, so conservation holds to round-off for any step size;
* the WP1 dataset needs uniform sampling in time, so outputs land exactly on a uniform sampling grid;
* the kinetic rates span orders of magnitude.  Near equilibrium they are < 0.1/fs, but a strongly
  excited distribution scatters at several per fs (worst for the random-k integral B0, whose rate grows
  with the cube of the density of states).  A fixed step is then either wasteful or unstable, so each
  step is limited to  h <= cfl / lambda,  with lambda the fastest collision rate (and the rate of change
  of f) reported by :meth:`enzkin.model.KineticModel.rhs`.  RK4 is stable on the negative real axis up to
  |lambda h| = 2.78; the default cfl = 1 leaves a factor > 2 for off-diagonal coupling.

``run`` returns a :class:`Trajectory` holding sampled distributions and the scalar diagnostics of
:meth:`enzkin.model.KineticModel.observe`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def rk4_step(rhs, t, y, dt, k1=None):
    if k1 is None:
        k1 = rhs(t, y)
    k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = rhs(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


@dataclass
class Trajectory:
    t: np.ndarray
    f: np.ndarray  # (n_samples, N)
    T_ph: np.ndarray
    obs: dict  # name -> array over samples
    E: np.ndarray
    meta: dict = field(default_factory=dict)

    def __getitem__(self, name):
        return self.obs[name]


def run(model, t0: float, t1: float, dt: float, y0=None, sample_every: int = 1, cfl: float = 1.0,
        adaptive: bool = True, omega_probe: float | None = None, observe: bool = True, progress=None,
        max_steps: int = 5_000_000) -> Trajectory:
    """Integrate ``model`` from t0 to t1 (fs).

    ``dt`` is the largest step taken and ``dt * sample_every`` the uniform sampling interval.  With
    ``adaptive`` (default) each step is shortened to cfl / lambda when the dynamics is fast.
    """
    y = model.initial_state() if y0 is None else np.array(y0, dtype=float)
    n_out = int(round((t1 - t0) / (dt * sample_every)))
    ts, fs, tph, obs = [], [], [], []

    def sample(t, y):
        ts.append(t)
        f, T_ph = model.split(y)
        fs.append(np.array(f))
        tph.append(T_ph)
        if observe:
            obs.append(model.observe(t, y, omega_probe))

    t, n_steps, h_min = t0, 0, dt
    sample(t, y)
    for i_out in range(1, n_out + 1):
        t_target = t0 + i_out * dt * sample_every
        while t_target - t > 1e-12 * max(1.0, abs(t_target)):
            if adaptive:
                k1, lam = model.rhs(t, y, want_stiffness=True)
                h = min(dt, cfl / lam if lam > 0 else dt, t_target - t)
            else:
                k1, h = model.rhs(t, y), min(dt, t_target - t)
            y = rk4_step(model.rhs, t, y, h, k1)
            t += h
            n_steps += 1
            h_min = min(h_min, h)
            if n_steps > max_steps:
                raise RuntimeError(f"exceeded {max_steps} steps at t = {t:.3f} fs (h = {h:.2e} fs)")
        t = t_target
        sample(t, y)
        if progress is not None:
            progress(i_out, n_out, t)
    obs_arr = {k: np.array([o[k] for o in obs]) for k in obs[0]} if obs else {}
    return Trajectory(np.array(ts), np.array(fs), np.array(tph), obs_arr, model.grid.E.copy(),
                      meta=dict(dt=dt, sample_dt=dt * sample_every, cfl=cfl, adaptive=adaptive,
                                n_steps=n_steps, h_min=h_min, n_rhs=model.n_rhs))
