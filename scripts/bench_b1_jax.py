#!/usr/bin/env python3
"""Benchmark the JAX implementation of the exact collision integral (gate G1b, GPU part).

Run this on the machine whose speed you want to know.  On a GPU node install a CUDA build of JAX first:

    pip install -U "jax[cuda12]"
    python scripts/bench_b1_jax.py --N 128 256 --cells 1 16 100 400

It reports, for every grid size N and number of cells evaluated in one call, the wall time per call and
per cell, and what that means for one film run of the plan's size (100 cells x 4000 RK4 steps = 1.6e6
evaluations).  With ``--numba`` it also times the CPU numba kernel on the same machine for reference.
Results go to results/wp1/bench_b1_jax_<device>.json.

Memory.  The kernel holds  cells x d_batch x N^2  doubles per intermediate array.  If a case runs out of
memory, lower ``--d-batch`` (default 16) or the number of cells; out-of-memory cases are recorded and
skipped.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import jax
import numpy as np

from enzkin import thermo
from enzkin.collisions.ee import ElectronElectronB1
from enzkin.collisions.ee_jax import ElectronElectronB1Jax
from enzkin.grid import EnergyGrid
from enzkin.material import ito_sarkar2023

FILM_EVALS = 100 * 4000 * 4  # cells x steps x RK4 stages


def states(grid, mat, cells, seed=0):
    """A batch of hot, slightly non-thermal distributions with different temperatures (as in a film)."""
    rng = np.random.default_rng(seed)
    out = []
    for T in rng.uniform(500.0, 8000.0, size=cells):
        f = thermo.fd_of_density(grid, mat.n_e, float(T))
        out.append(np.clip(f * (1.0 + 0.02 * rng.standard_normal(grid.N)), 0.0, 1.0))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="*", default=[128, 256])
    ap.add_argument("--cells", type=int, nargs="*", default=[1, 16, 100])
    ap.add_argument("--d-batch", type=int, default=16)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--numba", action="store_true", help="also time the CPU numba kernel")
    args = ap.parse_args()

    dev = jax.devices()[0]
    print(f"JAX {jax.__version__} on {dev} ({dev.platform}); x64 = {jax.config.jax_enable_x64}")
    mat = ito_sarkar2023()
    rows = []
    for N in args.N:
        grid = EnergyGrid(mat.band, mat.E_max, N)
        op = ElectronElectronB1Jax(grid, mat, d_batch=args.d_batch)
        ref = ElectronElectronB1(grid, mat)
        numba_ms = None
        if args.numba:
            f1 = states(grid, mat, 1)[0]
            ref.rhs(f1)
            ts = []
            for i in range(args.reps):
                t = time.perf_counter()
                ref.rhs(f1 * (1 + 1e-9 * (i + 1)))
                ts.append(time.perf_counter() - t)
            numba_ms = 1e3 * float(np.median(ts))
        for cells in args.cells:
            f = states(grid, mat, cells)
            row = dict(N=N, cells=cells, d_batch=args.d_batch, numba_ms_per_cell=numba_ms)
            try:
                t = time.perf_counter()
                d, r = op.rhs_and_rate(f)
                d.block_until_ready()
                row["compile_s"] = time.perf_counter() - t
                ts = []
                for _ in range(args.reps):
                    t = time.perf_counter()
                    d, r = op.rhs_and_rate(f)
                    d.block_until_ready()
                    ts.append(time.perf_counter() - t)
                ms = 1e3 * float(np.median(ts))
                d_ref = ref.rhs(f[0])
                err = float(np.abs(np.asarray(d[0]) - d_ref).max() / np.abs(d_ref).max())
                row.update(ms_per_call=ms, ms_per_cell=ms / cells, max_rel_error_vs_numba=err,
                           film_run_minutes=FILM_EVALS / 100 * (ms / cells * 100) / 1e3 / 60 if cells else None)
                print(f"N={N:4d} cells={cells:4d}: {ms:9.2f} ms/call  {ms / cells:8.3f} ms/cell  "
                      f"film run {row['film_run_minutes']:7.1f} min  (compile {row['compile_s']:.1f} s, "
                      f"error vs numba {err:.1e}" + (f", numba {numba_ms:.2f} ms/cell)" if numba_ms else ")"), flush=True)
            except Exception as exc:  # out of memory and similar: record and continue
                row["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
                print(f"N={N:4d} cells={cells:4d}: FAILED {row['error']}", flush=True)
            rows.append(row)
    out = Path(f"results/wp1/bench_b1_jax_{dev.platform}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(device=str(dev), platform=dev.platform, jax=jax.__version__, host=platform.platform(),
                                   film_evaluations=FILM_EVALS, rows=rows), indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
