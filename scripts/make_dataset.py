#!/usr/bin/env python3
"""WP1 dataset generator: trajectories of f(E, t) with their collision terms, for training the WP2 operator.

Sampling follows the project plan (WP1, "Dataset"):

    pump wavelength      1.1 - 1.6 um          (uniform)
    pump duration        20 - 250 fs FWHM      (log-uniform)
    in-film intensity    0.1 GW/cm^2 - 10 TW/cm^2   (log-uniform; |E_loc|^2 / 2 Z0)
    pulse shape          single | chirp-free double pulse (second pulse 0.2-1x, delay 0.5-4 FWHM)
    material             n_e, m0, C, c_ee, D_def varied by +-25 % around the [14, 15] values;
                         eps_inf and gamma_imp varied by +-10 % and +-50 %

Each trajectory stores, on a uniform time grid (default 2 fs):

    f(E, t)                         the distribution
    Q_ee(E, t), Q_eph(E, t)         the collision terms the operator must learn (computed with B1)
    Q_exc(E, t)                     the excitation term (analytic; kept for checks)
    rate_ee(E, t)                   collision rate entering the optical damping
    T_ph, T_e, mu, wp2, gamma, eps  scalars

plus the sampled parameters.  Along each trajectory ``--perturb`` extra off-trajectory states are stored:
the trajectory state plus a smooth random perturbation with zero number and energy, with their collision
terms, so the operator also sees states the solver would never visit (plan: "perturbed off-trajectory
states").

The full plan calls for ~2e4 trajectories on 512 points over 2 ps (about 40 GB).  With B1 at 512 points
one trajectory costs ~1.5 CPU-hours on 4 cores (docs/gates/G1_report.md), so the full set needs a GPU port
of B1 or a cluster.  Runs are independent and seeded per index: use --start/--n to split work.

Usage:  python scripts/make_dataset.py --n 64 --N 128 --out data/pilot.h5
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import numpy as np

from enzkin import thermo
from enzkin.collisions.ee import project_conserving
from enzkin.drive import FieldDrive, Pulse, gaussian_envelope
from enzkin.integrate import run
from enzkin.material import ito_sarkar2023
from enzkin.model import KineticModel

EPS_INF_REF, GAMMA_IMP_REF = 3.4509, 0.0631  # linear optics calibrated to [15] (scripts/validate_un2023.py)


class DoublePulse(Pulse):
    """Two Gaussian pulses of equal duration; ``ratio`` is the second-to-first intensity ratio."""

    def __init__(self, wavelength, fwhm, delay, ratio):
        super().__init__(wavelength, fwhm, 0.0)
        self.delay, self.ratio = delay, ratio

    def envelope(self, t):
        return gaussian_envelope(t, 0.0, self.fwhm) + self.ratio * gaussian_envelope(t, self.delay, self.fwhm)


def sample_case(rng: np.random.Generator) -> dict:
    u = rng.uniform
    case = dict(
        wavelength=u(1100.0, 1600.0),
        fwhm=float(np.exp(u(np.log(20.0), np.log(250.0)))),
        intensity=float(np.exp(u(np.log(1e8), np.log(1e13)))),
        double=bool(rng.random() < 0.25),
        n_e=1.5 * u(0.75, 1.25),
        m0_rel=0.3964 * u(0.75, 1.25),
        C=0.4191 * u(0.75, 1.25),
        c_ee=(3.0 / 16.0) * u(0.75, 1.25),
        D_def=17.2 * u(0.75, 1.25),
        eps_inf=EPS_INF_REF * u(0.9, 1.1),
        gamma_imp=GAMMA_IMP_REF * u(0.5, 1.5),
    )
    case["delay"] = case["fwhm"] * u(0.5, 4.0) if case["double"] else 0.0
    case["ratio"] = u(0.2, 1.0) if case["double"] else 0.0
    return case


def perturbed_states(model, f, T_ph, rng, k: int, amplitude: float = 0.05):
    """k states f + delta with delta smooth, confined to where f(1-f) > 0, carrying zero number and energy."""
    g = model.grid
    out_f, out_qee, out_qeph = [], [], []
    for _ in range(k):
        modes = rng.integers(1, 12)
        x = g.E / g.E_max
        bump = sum(rng.standard_normal() * np.cos(np.pi * m * x + rng.uniform(0, 2 * np.pi)) for m in range(1, modes + 1))
        delta = amplitude * rng.uniform(0.2, 1.0) * bump / max(np.abs(bump).max(), 1e-12) * np.sqrt(f * (1 - f))
        delta = project_conserving(g, delta, f)
        fp = np.clip(f + delta, 0.0, 1.0)
        out_f.append(fp)
        out_qee.append(model.ee.rhs(fp))
        out_qeph.append(model.eph.rhs(fp, T_ph)[0])
    return np.array(out_f), np.array(out_qee), np.array(out_qeph)


def one_trajectory(case: dict, N: int, t_span: float, dt_sample: float, n_perturb: int, rng):
    mat = ito_sarkar2023(**{k: case[k] for k in ("n_e", "m0_rel", "C", "c_ee", "D_def", "eps_inf", "gamma_imp")})
    pulse = (DoublePulse(case["wavelength"], case["fwhm"], case["delay"], case["ratio"]) if case["double"]
             else Pulse(case["wavelength"], case["fwhm"]))
    model = KineticModel(mat, N=N, drive=FieldDrive(pulse, case["intensity"]))
    t0 = -2.5 * case["fwhm"]
    dt = min(0.5, case["fwhm"] / 60.0)
    every = max(int(round(dt_sample / dt)), 1)
    tr = run(model, t0, t0 + t_span, dt, sample_every=every)

    n_t = len(tr.t)
    q_ee, q_eph, q_exc, r_ee = (np.empty((n_t, N)) for _ in range(4))
    for i in range(n_t):
        y = np.concatenate([tr.f[i], [tr.T_ph[i], 0.0, 0.0]])
        _, parts = model.rhs(tr.t[i], y, return_parts=True)
        q_ee[i], q_eph[i], q_exc[i], r_ee[i] = parts["ee"], parts["eph"], parts["exc"], parts["rate_ee"]
    pert = None
    if n_perturb > 0:
        idx = rng.choice(n_t, size=min(n_perturb, n_t), replace=False)
        fs, qe, qp = zip(*(perturbed_states(model, tr.f[i], tr.T_ph[i], rng, 1) for i in idx))
        pert = dict(f=np.concatenate(fs), q_ee=np.concatenate(qe), q_eph=np.concatenate(qp), T_ph=tr.T_ph[idx], t=tr.t[idx])
    return model, tr, dict(q_ee=q_ee, q_eph=q_eph, q_exc=q_exc, rate_ee=r_ee), pert


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=16, help="number of trajectories")
    ap.add_argument("--start", type=int, default=0, help="index of the first trajectory (sets the seed)")
    ap.add_argument("--N", type=int, default=128, help="energy points")
    ap.add_argument("--span", type=float, default=2000.0, help="trajectory length, fs")
    ap.add_argument("--dt-sample", type=float, default=2.0, help="sampling interval, fs")
    ap.add_argument("--perturb", type=int, default=32, help="off-trajectory states per trajectory")
    ap.add_argument("--seed", type=int, default=20261117)
    ap.add_argument("--out", default="data/pilot.h5")
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out, "a") as h5:
        h5.attrs.update(N=args.N, span_fs=args.span, dt_sample_fs=args.dt_sample, seed=args.seed,
                        e_max_ev=4.0, units="eV, fs, nm, K; f dimensionless; Q in 1/fs",
                        ee_level="B1 (screened Coulomb, kappa recomputed from f)")
        for i in range(args.start, args.start + args.n):
            name = f"traj_{i:06d}"
            if name in h5:
                continue
            rng = np.random.default_rng([args.seed, i])
            case = sample_case(rng)
            t0 = time.perf_counter()
            model, tr, q, pert = one_trajectory(case, args.N, args.span, args.dt_sample, args.perturb, rng)
            g = h5.create_group(name)
            g.attrs.update({k: v for k, v in case.items()})
            g.attrs.update(wall_s=time.perf_counter() - t0, n_steps=tr.meta["n_steps"], h_min_fs=tr.meta["h_min"])
            kw = dict(compression="gzip", compression_opts=4, shuffle=True)
            g.create_dataset("t", data=tr.t)
            g.create_dataset("E", data=tr.E)
            g.create_dataset("f", data=tr.f, **kw)
            for k, v in q.items():
                g.create_dataset(k, data=v, **kw)
            for k in ("T_e", "mu", "wp2", "gamma", "eps_re", "eps_im", "n", "u", "u_abs", "nonthermal", "audit"):
                g.create_dataset(k, data=tr[k])
            g.create_dataset("T_ph", data=tr.T_ph)
            if pert is not None:
                p = g.create_group("perturbed")
                for k, v in pert.items():
                    p.create_dataset(k, data=v, **(kw if v.ndim == 2 else {}))
            h5.flush()
            audit = np.abs(tr["audit"] - tr["audit"][0]).max() / max(tr["u_abs"][-1], 1e-30)
            print(f"{name}: {case['fwhm']:6.1f} fs  {case['intensity']:.2e} W/cm2  {case['wavelength']:.0f} nm  "
                  f"T_e,max {tr['T_e'].max():7.0f} K  steps {tr.meta['n_steps']:6d}  audit {audit:.1e}  "
                  f"f in [{tr['f_min'].min():.1e}, {tr['f_max'].max():.6f}]  {time.perf_counter() - t0:.0f} s", flush=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
