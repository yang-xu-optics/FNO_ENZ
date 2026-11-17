#!/usr/bin/env python3
"""Figures for the WP1 reports, drawn from the JSON/NPZ written by the other scripts.

    results/wp1/fig_validation_un2023.png   peak electron temperature and absorptivity vs pump intensity, this work vs [15]
    results/wp1/fig_gate1_cost.png          cost per collision-integral evaluation vs grid size, B1 vs B0
    results/wp1/fig_gate1_transients.png    change of Re(eps) for three pulse durations: B1, B0, thermal model
    results/wp1/fig_convergence_sweep.png   change of the permittivity transient when the energy step is halved
    results/wp1/fig_band_width.png          change of the transient when the band is widened, vs top-cell occupation

Style: one y-axis per panel, small multiples instead of dual axes, 2 px lines, recessive grid, a legend
plus direct labels, three categorical hues (blue, orange, aqua) that pass colour-vision checks on the
light surface used here.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

R = Path("results/wp1")
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0, "lines.markersize": 7,
    "legend.frameon": False, "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
})


def fig_validation():
    d = json.loads((R / "validation_un2023.json").read_text())
    rows = [r for r in d["rows"] if r["mode"] == "film" and r["fwhm"] > 100]
    i0 = np.array([r["I0"] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), constrained_layout=True)
    for ax, key, title, ylab in ((axes[0], "Te_max", "Peak electron temperature", "T_e,max (K)"),
                                 (axes[1], "absorptivity", "Fraction of the pump absorbed in the film", "absorptivity")):
        ours = np.array([r[key] for r in rows])
        paper = np.array([r["target"][key] for r in rows], dtype=float)
        ax.plot(i0, paper, "-o", color=BLUE, label="Un et al. 2023 [15]")
        ax.plot(i0, ours, "-s", color=ORANGE, label="this work (enzkin, B1)")
        ax.set_xscale("log")
        ax.set_xticks(i0)
        ax.set_xticklabels([f"{x:g}" for x in i0])
        ax.minorticks_off()
        ax.set_xlabel("pump peak intensity in air (GW/cm²)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        up = key == "Te_max"
        ax.annotate("[15]", (i0[-1], paper[-1]), xytext=(6, -3 if up else 6), textcoords="offset points", color=INK2, fontsize=9)
        ax.annotate("this work", (i0[-1], ours[-1]), xytext=(-52, 8 if up else -14), textcoords="offset points", color=INK2, fontsize=9)
    axes[0].legend(loc="upper left")
    axes[1].set_ylim(0, 1)
    fig.suptitle("220 fs pump on 40 nm ITO under gold: self-consistent thin-film run vs the published model",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.savefig(R / "fig_validation_un2023.png", dpi=200)
    plt.close(fig)


def fig_cost():
    d = json.loads((R / "gate1.json").read_text())
    n = np.array([c["N"] for c in d["cost"]])
    b1 = np.array([c["b1_ms"] for c in d["cost"]])
    b0 = np.array([c["b0_ms"] for c in d["cost"]])
    fig, ax = plt.subplots(figsize=(5.6, 4.0), constrained_layout=True)
    ax.plot(n, b1, "-o", color=BLUE, label="B1: screened Coulomb, O(N³)")
    ax.plot(n, b0, "-s", color=ORANGE, label="B0: random-k by FFT, O(N log N)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(n)
    ax.set_xticklabels([str(x) for x in n])
    ax.set_xlabel("energy grid points N")
    ax.set_ylabel("wall time per evaluation (ms)")
    ax.set_title(f"Cost of one collision-integral evaluation ({d['machine']['numba_threads']} CPU threads)")
    for x, y1, y0 in zip(n, b1, b0):
        if y1 / y0 >= 3:  # label the ratio where the two curves are far enough apart to hold it
            ax.annotate(f"×{y1 / y0:.0f}", (x, np.sqrt(y1 * y0)), ha="center", va="center", color=INK2, fontsize=9)
    ax.annotate("B1", (n[-1], b1[-1]), xytext=(6, -3), textcoords="offset points", color=INK2)
    ax.annotate("B0", (n[-1], b0[-1]), xytext=(6, -3), textcoords="offset points", color=INK2)
    ax.legend(loc="upper left")
    fig.savefig(R / "fig_gate1_cost.png", dpi=200)
    plt.close(fig)


def fig_transients():
    d = json.loads((R / "gate1.json").read_text())
    z = np.load(R / "gate1_series.npz")
    picks = [(i, c["case"]) for i, c in enumerate(d["cases"])
             if c["case"]["mode"] == "energy" and abs(c["case"]["level"] - 1e8) < 1 and c["case"]["fwhm"] in (20.0, 50.0, 220.0)]
    fig, axes = plt.subplots(1, len(picks), figsize=(3.5 * len(picks), 3.7), constrained_layout=True, sharey=True)
    for ax, (ic, case) in zip(np.atleast_1d(axes), picks):
        for kind, col, lab in (("B1", BLUE, "B1 (reference)"), ("B0", ORANGE, "B0"), ("thermal", AQUA, "thermal (two-temperature)")):
            t, y = z[f"c{ic}_{kind}_t"], z[f"c{ic}_{kind}_eps_re"]
            ax.plot(t, y - y[0], color=col, label=lab)
        m = d["cases"][ic]["miss"]
        ax.set_title(f"{case['fwhm']:.0f} fs pump")
        ax.text(0.98, 0.05, f"largest deviation from B1\nB0 {100 * m['B0']['eps_re']:.1f} %   thermal {100 * m['thermal']['eps_re']:.1f} %",
                transform=ax.transAxes, ha="right", va="bottom", color=INK2, fontsize=8.5)
        ax.set_xlabel("time relative to pump peak (fs)")
    np.atleast_1d(axes)[0].set_ylabel("change of Re ε at the pump frequency")
    np.atleast_1d(axes)[-1].legend(loc="upper right", fontsize=8.5)  # the long-pulse panel has free space there
    fig.suptitle("Same absorbed energy (1×10⁸ J/m³), three pulse durations: permittivity transient by collision model",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.savefig(R / "fig_gate1_transients.png", dpi=200)
    plt.close(fig)


def fig_convergence():
    """Grid-convergence sweep: change of the permittivity transient when the energy step is halved."""
    d = json.loads((R / "convergence_sweep.json").read_text())
    rows = d["summary"]["grid"]
    worst = lambda m: 100 * max(m["eps_re"], m["eps_im"])
    series = (("128", 4.0, BLUE, "o", "128 points against 256 (4 eV band)"),
              ("256", 4.0, ORANGE, "s", "256 points against 512 (4 eV band)"),
              ("192", 6.0, AQUA, "^", "192 points against 384 (6 eV band)"))
    fig, ax = plt.subplots(figsize=(7.4, 4.2), constrained_layout=True)
    for n, emax, col, mk, lab in series:
        pts = sorted((e["Te_max"], worst(e["halving"][n])) for e in rows if e["E_max"] == emax and n in e["halving"])
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, mk, color=col, label=lab, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 2.0)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate("1 % criterion", (ax.get_xlim()[0], 1.0), xytext=(6, 5), textcoords="offset points", color=INK2, fontsize=9)
    ax.set_xlabel("peak electron temperature of the run (K)")
    ax.set_ylabel("change of the permittivity transient\n(% of its amplitude, larger of Re ε and Im ε)")
    ax.set_title("Halving the energy step changes the transient by less than 1 % everywhere")
    ax.legend(loc="lower left")
    fig.savefig(R / "fig_convergence_sweep.png", dpi=200)
    plt.close(fig)


def fig_band_width():
    """Band-width effect at equal energy step, against the occupation of the narrower band's top cell."""
    d = json.loads((R / "convergence_sweep.json").read_text())
    band = d["summary"]["band_width"]
    worst = lambda m: 100 * max(m["eps_re"], m["eps_im"])
    series = ((4.0, BLUE, "o", "4 eV against 6 eV"), (6.0, ORANGE, "s", "6 eV against 8 eV"), (8.0, AQUA, "^", "8 eV against 10 eV"))
    fig, ax = plt.subplots(figsize=(7.4, 4.2), constrained_layout=True)
    for e_lo, col, mk, lab in series:
        pts = sorted((b["f_top_narrow"], worst(b["miss"])) for b in band if b["E_narrow"] == e_lo)
        if pts:
            x, y = zip(*pts)
            ax.plot(x, y, mk, color=col, label=lab, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate("1 % criterion", (ax.get_xlim()[0], 1.0), xytext=(6, 5), textcoords="offset points", color=INK2, fontsize=9)
    ax.set_xlabel("largest occupation reached in the top cell of the narrower band")
    ax.set_ylabel("change of the permittivity transient\n(% of its amplitude, larger of Re ε and Im ε)")
    ax.set_title("A band is wide enough while its top cell stays nearly empty")
    ax.legend(loc="upper left")
    fig.savefig(R / "fig_band_width.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    for fn in (fig_validation, fig_cost, fig_transients, fig_convergence, fig_band_width):
        try:
            fn()
            print("wrote", fn.__name__)
        except FileNotFoundError as exc:
            print("skipped", fn.__name__, "-", exc)
