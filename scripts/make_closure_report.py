#!/usr/bin/env python3
"""Write docs/gates/GT2_fewcycle_closure_report.md and its two figures from the result files.

Inputs (results/tier2/):
    densmat_validation.json        scripts/validate_densmat.py
    closure_test.json              scripts/fewcycle_closure_test.py
    closure_test_un2026.json       scripts/fewcycle_closure_test.py --un2026   (optional robustness check)
    closure_test_series.npz        time series of four conditions
Outputs:
    docs/gates/GT2_fewcycle_closure_report.md
    results/tier2/fig_closure_errors.png, results/tier2/fig_closure_traces.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "tier2"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0, "lines.markersize": 7,
    "legend.frameon": False, "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
})
L1, L2, L2F, L2S, L3 = "L1 thermal", "L2 17 poles", "L2f 17 poles, fixed widths", "L2 9 poles", "L3 Drude"
FILM_CELLS = 100


def fig_errors(rows):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5), constrained_layout=True, sharey=True)
    series = ((L1, BLUE, "o", "one temperature, all pair oscillators (L1)"),
              (L2, ORANGE, "s", "one temperature, 17 fitted poles (L2)"),
              (L3, AQUA, "^", "one temperature, Drude pole (L3, Tier 2a)"))
    for ax, key, title in ((axes[0], "err_tot", "Error in the total current"),
                           (axes[1], "err_nl", "Error in its nonlinear part")):
        for name, col, mk, lab in series:
            x = [r["T_e_max"] for r in rows]
            y = [100 * r["models"][name][key] for r in rows]
            ax.plot(x, y, mk, color=col, label=lab, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if v in (2000, 5000, 20000) else ""))
        ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.annotate("1 % of the peak current", (ax.get_xlim()[0], 1.0), xytext=(6, 5), textcoords="offset points", color=INK2, fontsize=9)
        ax.set_xlabel("peak electron temperature of the run (K)")
        ax.set_title(title)
    axes[0].set_ylabel("largest deviation from the density matrix\n(% of the peak current)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=3, fontsize=9)
    fig.savefig(R / "fig_closure_errors.png", dpi=200)
    plt.close(fig)


def fig_traces(series, tag, caption):
    t = series[f"{tag}_t"]
    ref = series[f"{tag}_dJ_L0 density matrix"]
    jmax = np.abs(series[f"{tag}_J_L0 density matrix"]).max()
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 5.6), constrained_layout=True, sharex=True, height_ratios=(1.4, 1))
    ax = axes[0]
    ax.plot(t, 100 * ref / jmax, color=INK, linewidth=2.6, label="density matrix (L0)")
    for name, col, lab in ((L1, BLUE, "one temperature, all pairs (L1)"), (L2, ORANGE, "17 fitted poles (L2)"), (L3, AQUA, "Drude pole (L3)")):
        ax.plot(t, 100 * series[f"{tag}_dJ_{name}"] / jmax, color=col, linewidth=1.6, label=lab)
    ax.set_ylabel("nonlinear part of the current\n(% of the peak current)")
    ax.set_title(caption)
    ax.legend(loc="lower left", fontsize=9)
    ax = axes[1]
    for name, col in ((L1, BLUE), (L2, ORANGE), (L3, AQUA)):
        ax.plot(t, 100 * (series[f"{tag}_dJ_{name}"] - ref) / jmax, color=col, linewidth=1.6)
    ax.axhline(0.0, color=INK2, linewidth=0.8)
    ax.set_ylabel("closure minus density matrix\n(% of the peak current)")
    ax.set_xlabel("time (fs)")
    ax.set_xlim(t[0], min(t[-1], 30.0))
    fig.savefig(R / "fig_closure_traces.png", dpi=200)
    plt.close(fig)


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def main():
    val = json.loads((R / "densmat_validation.json").read_text())
    res = json.loads((R / "closure_test.json").read_text())
    alt_path = R / "closure_test_un2026.json"
    alt = json.loads(alt_path.read_text()) if alt_path.exists() else None
    series = dict(np.load(R / "closure_test_series.npz"))
    rows = res["rows"]
    fig_errors(rows)
    fig_traces(series, "f175_tau5_E2.4_cep0.00", "Strongest 5 fs pulse (2.4 V/nm, 175 THz): what each closure misses")

    def worst(name, key, sel=lambda r: True):
        rr = [r for r in rows if sel(r)]
        i = int(np.argmax([r["models"][name][key] for r in rr]))
        return rr[i]["models"][name][key], rr[i]

    def cond(r):
        return f'{r["fwhm"]:.0f} fs, {r["E0"]:g} V/nm'

    w = {n: worst(n, "err_tot") for n in (L1, L2, L2S, L2F, L3)}
    wn = {n: worst(n, "err_nl") for n in (L1, L2, L2S, L2F, L3)}
    med = {n: float(np.median([r["models"][n]["err_tot"] for r in rows])) for n in (L1, L2, L2S, L2F, L3)}
    n_over = {n: sum(r["models"][n]["err_tot"] > 0.01 for r in rows) for n in (L1, L2, L3)}
    n_over5 = {n: sum(r["models"][n]["err_tot"] > 0.05 for r in rows) for n in (L1, L2, L3)}
    cep_pairs = {}
    for r in rows:
        cep_pairs.setdefault((r["f_thz"], r["fwhm"], r["E0"]), []).append(r["models"][L1]["err_tot"])
    cep_spread = max(abs(v[0] - v[1]) for v in cep_pairs.values() if len(v) == 2)
    cost = res["cost_60fs"]
    film_min = {k: v["wall_s"] * FILM_CELLS / 60 for k, v in cost.items()}
    u_dev = {n: max(abs(r["models"][n]["U_ratio"] - 1) for r in rows) for n in (L1, L2, L3)}
    h3 = {n: (min(r["models"][n]["h3_ratio"] for r in rows), max(r["models"][n]["h3_ratio"] for r in rows)) for n in (L1, L2, L3)}
    nl_range = (min(r["nl"] for r in rows), max(r["nl"] for r in rows))
    te_range = (min(r["T_e_max"] for r in rows), max(r["T_e_max"] for r in rows))
    ftop = max(r["f_top_max"] for r in rows)

    table = []
    for r in sorted(rows, key=lambda r: (r["f_thz"], r["fwhm"], r["E0"], r["cep"])):
        m = r["models"]
        table.append(f'| {r["f_thz"]:.0f} | {r["fwhm"]:.0f} | {r["E0"]:g} | {"0" if r["cep"] == 0 else "π/2"} | {r["T_e_max"]:.0f} | {pct(r["nl"], 0)} | '
                     f'{pct(m[L1]["err_tot"])} | {pct(m[L2]["err_tot"])} | {pct(m[L2F]["err_tot"])} | {pct(m[L2S]["err_tot"])} | {pct(m[L3]["err_tot"])} | {pct(m[L3]["err_nl"])} |')

    vr = {round(x["E0"], 1): x for x in val["runs"]}
    pk = val["paper"]["peaks_THz"]
    val_rows = [
        f'| ε(175 THz), 300 K | −1.6 + 1.4i | {val["eps_175"][0]:.2f} (fitted) + {val["eps_175"][1]:.2f}i |',
        f'| Zero of Re ε | 215 THz | {val["zero_crossing_THz"]:.0f} THz |',
        f'| Energy given to the electrons, 1.2 V/nm | 0.63 eV/nm³ | {vr[1.2]["U_field"]:.3f} |',
        f'| Instantaneous absorbed power, 1.2 V/nm | −0.14 … +0.33 eV/nm³/fs | {vr[1.2]["power_min"]:.3f} … +{vr[1.2]["power_max"]:.3f} |',
        f'| Local-field spectral peak at 0.4, 1.2, 2.4 V/nm | {pk["0.4"]}, {pk["1.2"]}, {pk["2.4"]} THz | {vr[0.4]["peak_THz"]:.0f}, {vr[1.2]["peak_THz"]:.0f}, {vr[2.4]["peak_THz"]:.0f} THz |',
        f'| ε(160 THz) after the pump, 1.2 V/nm ([62] Fig. 4) | about −1.0 + 2.25i | {vr[1.2]["eps160_end"][0]:.2f} + {vr[1.2]["eps160_end"][1]:.2f}i |',
        f'| 1/η_ee at E_F, 5 fs after the peak, 1.8 V/nm | about 12 fs | {1 / vr[1.8]["eta_ee_EF_5fs"]:.0f} fs |',
    ]
    val_rows = [row.replace(" -", " −").replace("| -", "| −") for row in val_rows]
    alt_text = ""
    if alt:
        a1 = max(r["models"][L1]["err_tot"] for r in alt["rows"])
        a2 = max(r["models"][L2]["err_tot"] for r in alt["rows"])
        a3 = max(r["models"][L3]["err_tot"] for r in alt["rows"])
        alt_text = (f"\n**Other parameter set.** With Table S1 of [62] (D = 8.6 eV, Brooks–Herring impurities) the worst errors are "
                    f"{pct(a1)} % (L1), {pct(a2)} % (L2) and {pct(a3)} % (L3): the same picture.\n")

    unc = res["model_uncertainty"]
    l1_at = {(r["fwhm"], r["E0"]): r["models"][L1]["err_tot"] for r in rows if r["f_thz"] == 175.0 and r["cep"] == 0.0}
    unc_rows = [f'| {u["fwhm"]:.0f} | {u["E0"]:g} | {pct(l1_at[(u["fwhm"], u["E0"])])} | {pct(u["c_reg"])} | {pct(u["umklapp_half"])} |' for u in unc]
    unc_min = min(max(u["c_reg"], u["umklapp_half"]) for u in unc)
    unc_max = max(max(u["c_reg"], u["umklapp_half"]) for u in unc)
    below_unc = all(l1_at[(u["fwhm"], u["E0"])] < max(u["c_reg"], u["umklapp_half"]) for u in unc)
    # verdict: no cost reason, and the most a learned operator could add is below the model's own uncertainty
    not_worth = film_min["L0 density matrix"] < 60 and below_unc
    report = f"""# Test GT2: is a learned few-cycle constitutive law worth building?

Generated by `scripts/make_closure_report.py` from `results/tier2/`. Do not edit by hand. Decision test requested by the project lead on 2026-11-24 ("option 1" of working-log session 5).

## Verdict

**{"No, not on the density-matrix model of [16, 62]." if not_worth else "Not settled: see the table."}** A learned operator would have to earn its place in one of three ways, and none holds here.

| Reason to learn an operator | Test | Result |
| --- | --- | --- |
| The microscopic model is too expensive to run in every cell | Wall time of the full density matrix | **Not the case.** {cost["L0 density matrix"]["wall_s"]:.0f} s per cell for a 60 fs window ({cost["L0 density matrix"]["oscillators"]:,} oscillators, {cost["L0 density matrix"]["collision_calls"]} collision evaluations) on 4 cores; about {film_min["L0 density matrix"]:.0f} min for {FILM_CELLS} cells. Revision 2 of the plan assumed it was too expensive without measuring it |
| The response carries a memory that no simple state variable captures | Replace the populations by one number, the electron temperature (L1) | **Not the case.** The field-resolved current differs by at most {pct(w[L1][0])} % of its peak ({cond(w[L1][1])}), median {pct(med[L1])} %, while the response itself is nonlinear by {pct(nl_range[0], 0)}–{pct(nl_range[1], 0)} %. It exceeds 1 % in {n_over[L1]} of {len(rows)} conditions and 5 % in {n_over5[L1]} |
| A compact fast law is needed and physics does not supply one | 17 fixed poles + one temperature (L2) | **Physics supplies it.** At most {pct(w[L2][0])} % from the density matrix ({cond(w[L2][1])}), median {pct(med[L2])} %, at {cost[L2]["wall_s"] * 1e3:.0f} ms per cell |

What a trained operator could add at best is the gap between L1 and the density matrix: {pct(med[L1])} % of the peak current in the median and {pct(w[L1][0])} % in the worst case. The model itself is less certain than that. Two choices the papers leave open move its current by {pct(unc_min)}–{pct(unc_max)} % (section B), and it leaves out the coupling between coherences and every absorption channel of a film (`docs/theory/tier2_density_matrix.md` §5).

The test also produced three results that stand on their own:

1. **The plan's Tier 2a closure (one Drude pole) is not adequate for this particle**: {pct(med[L3])} % median and {pct(w[L3][0])} % worst error. Most of it is linear: a 4 nm particle absorbs through its size, which a Drude pole fitted at one frequency does not reproduce across a few-cycle bandwidth. With each model's own linear response removed, the Drude pole still misses up to {pct(wn[L3][0])} % of the peak current.
2. **A few-cycle law must let the damping follow the temperature in time.** The same 17 poles with temperature-independent widths fit χ(ω; T_e) as well at every frozen temperature, but in time they err by up to {pct(w[L2F][0])} % against {pct(w[L2][0])} %: the electron temperature changes within the memory time of the polarisation. This is the non-adiabatic effect of the constitutive law, isolated.
3. **Eq. (16) of [16] is wrong as printed** and was reconstructed; with it the model reproduces the published results of [16, 62] to a few percent (section A).

## What was run

- Particle of [16, 62]: 4 nm ITO sphere in air, quasi-static field, parameters of [14, 15] ({res["parameter_set"]}), ε∞ = {res["eps_inf"]:.2f} fitted to Re ε(175 THz) = −1.6. Band of {res["E_max"]:g} eV on {res["N"]} points; the top cell never fills beyond {ftop:.0e}.
- {len(rows)} pulses: 175 THz with 5, 10, 20 fs duration and peak incident fields 0.3–2.4 V/nm, carrier-envelope phase 0 and π/2; two pulses at 230 THz. Peak electron temperatures {te_range[0]:.0f}–{te_range[1]:.0f} K.
- Truth L0: populations with e–e and e–ph collisions plus all pair oscillators. Closures: L1 = Fermi–Dirac populations at T_e(t); L2 = 17 fitted poles whose widths rise with the collision rate (9 poles and fixed widths as variants); L3 = one Drude pole. All closures are heated by the instantaneous work E_loc·Ṗ.
- Observable: the current J = Ṗ, which sets the field radiated by a thin layer of particles. "Error" is the largest deviation from L0 over the run, in percent of the peak of J. "Nonlinear part" is what remains after subtracting the scaled response to a weak pulse, for each model separately.

## A. The density-matrix model reproduces [16, 62]

`scripts/validate_densmat.py`; only ε∞ is adjusted.

| Quantity | [16] / [62] | this implementation |
| --- | --- | --- |
{chr(10).join(val_rows)}

The dipole element of [16] Eq. (16) gives {val["sum_rule_at_EF"]["eq16_as_printed"]:.1e} times the required oscillator strength at E_F as printed; the reconstructed one gives {val["sum_rule_at_EF"]["reconstructed"]:.2f} and matches [16] Fig. 1b away from the diagonal (log₁₀|D|² = {val["dipole_log10_C2m2"]["0.5 E_F"]:.2f} at 0.5 E_F against about {val["dipole_paper_by_eye"]["0.5 E_F"]} read from the figure). One disagreement remains: [62] quotes "∼4500 K" where this model has {vr[1.8]["T_e_5fs"]:.0f} K for the same e–e time; the absorbed energies agree, and the e–e time fits the higher temperature.

## B. Closure errors, all conditions

Errors in percent of the peak current. "nonlin." is the size of the nonlinear part of the true response.

| f (THz) | τ (fs) | E₀ (V/nm) | CEP | T_e,max (K) | nonlin. | L1 | L2, 17 poles | L2, fixed widths | L2, 9 poles | L3 Drude | L3, nonlinear part only |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(table)}

![Closure errors against peak electron temperature](../../results/tier2/fig_closure_errors.png)

![Nonlinear current of the strongest 5 fs pulse](../../results/tier2/fig_closure_traces.png)

**How certain is the truth?** The same runs with two choices varied that [16, 62] do not fix, each recalibrated to the same cold permittivity: the regularisation constant of the dipole element lowered by 10 %, and only half of the e–e rate allowed to dephase the polarisation ([16] says the Umklapp share is unknown and of order one).

| τ (fs) | E₀ (V/nm) | L1 against the density matrix (%) | density matrix, c_reg 0.9 against 1.0 (%) | density matrix, half the e–e dephasing (%) |
| --- | --- | --- | --- | --- |
{chr(10).join(unc_rows)}

- The carrier-envelope phase changes the L1 error by at most {pct(cep_spread, 2)} percentage points.
- Absorbed energy: within {pct(u_dev[L1])} % (L1), {pct(u_dev[L2])} % (L2) and {pct(u_dev[L3])} % (L3) of the density matrix.
- Third-harmonic amplitude of the current relative to the density matrix: {h3[L1][0]:.2f}–{h3[L1][1]:.2f} (L1), {h3[L2][0]:.2f}–{h3[L2][1]:.2f} (L2), {h3[L3][0]:.2f}–{h3[L3][1]:.2f} (L3). In this model the harmonics come from the 2ω ripple of the absorbed energy, which a temperature state carries.
{alt_text}
## C. Cost

| Model | State per cell | Wall time, 60 fs window (s) | {FILM_CELLS} cells (min) |
| --- | --- | --- | --- |
| L0 density matrix | {res["N"]} populations + {cost["L0 density matrix"]["oscillators"]:,} oscillators | {cost["L0 density matrix"]["wall_s"]:.1f} | {film_min["L0 density matrix"]:.0f} |
| L1 one temperature, all pairs | 1 + {cost[L1]["oscillators"]:,} oscillators | {cost[L1]["wall_s"]:.1f} | {film_min[L1]:.1f} |
| L2 17 fitted poles | 1 + 17 oscillators | {cost[L2]["wall_s"]:.3f} | {film_min[L2]:.2f} |
| L3 Drude | 1 + 1 oscillator | {cost[L3]["wall_s"]:.3f} | {film_min[L3]:.2f} |

Measured on this 4-core machine with 0.02 fs steps; building the tables for the closures takes {res["build_s"]:.0f} s once per material.

## Consequences

- **Do not build the learned operator now.** The decision of Gate 1 carries over to the few-cycle regime: the microscopic model is cheap, and its reduction to one temperature plus a linear memory is accurate to a few percent and known in closed form.
- **Tier 2a should use the multi-pole form, not a single Drude pole, wherever the linear response is not Drude** (particles; films only if their dispersion requires it), and must make the pole widths follow the temperature. For bulk-like films the Drude pole is the right linear response, and the result here says the energy-balance part of Tier 2a is sound (absorbed energy within {pct(u_dev[L3])} %).
- **Nanocrystal experiments are now within reach without a film solver**: the validated particle model runs in seconds and applies directly to field-resolved measurements on ITO nanocrystals (Herbst et al. 2026, cited by [62]; not yet read).
- **Where learning could still enter:** not as a constitutive law, but as an emulator for inference over material parameters (option 2 of session 5), or later if the microscopic model gains physics that does not reduce to a temperature (coherence–coherence coupling, collision-assisted absorption in films). Neither exists yet.

## Limits

- One microscopic model, with the omissions listed in `docs/theory/tier2_density_matrix.md` §5. The test says a learned operator has nothing to learn from *this* model; it cannot speak for physics the model lacks.
- Uniform field in a 4 nm particle: no propagation, no film.
- The detection threshold of a field-resolved experiment is not known to me; errors are given in percent of the peak current so that any threshold can be applied.
- Pulses at 175 and 230 THz only; durations 5–20 fs; fields up to 2.4 V/nm (peak temperatures up to {te_range[1]:.0f} K).
"""
    out = ROOT / "docs" / "gates" / "GT2_fewcycle_closure_report.md"
    out.write_text(report, encoding="utf-8")
    print("wrote", out)
    print("worst err_tot:", {k: f"{pct(v[0])} % ({cond(v[1])})" for k, v in w.items()})
    print("median err_tot:", {k: f"{pct(v)} %" for k, v in med.items()})


if __name__ == "__main__":
    main()
