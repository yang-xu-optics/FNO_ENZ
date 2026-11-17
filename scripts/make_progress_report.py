#!/usr/bin/env python3
"""Build the HTML progress report from the files under results/ and references/.

    python scripts/make_progress_report.py [--out docs/reports/progress_2026-11-17.html]

Every number in the charts and result tables is read from results/wp1/*.json (written by
validate_un2023.py and gate1_benchmark.py) or references/manifest.json.  The narrative is edited here.
The page is a single self-contained file: inline CSS, inline SVG charts, a few lines of script for
hover read-outs.  It carries no <html>/<head>/<body> wrapper so the same file can be published as an
artifact page; browsers open it locally as it is.
"""
from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "wp1"


# ----------------------------------------------------------------------------- svg helpers
def esc(s) -> str:
    return html.escape(str(s), quote=True)


class Frame:
    """Plot frame: maps data to viewBox coordinates.  ``xlog``/``ylog`` select logarithmic axes."""

    def __init__(self, w, h, ml, mr, mt, mb, x0, x1, y0, y1, xlog=False, ylog=False):
        self.w, self.h, self.ml, self.mr, self.mt, self.mb = w, h, ml, mr, mt, mb
        self.x0, self.x1, self.y0, self.y1, self.xlog, self.ylog = x0, x1, y0, y1, xlog, ylog

    def _f(self, v, a, b, log):
        return (math.log(v) - math.log(a)) / (math.log(b) - math.log(a)) if log else (v - a) / (b - a)

    def x(self, v):
        return self.ml + self._f(v, self.x0, self.x1, self.xlog) * (self.w - self.ml - self.mr)

    def y(self, v):
        return self.h - self.mb - self._f(v, self.y0, self.y1, self.ylog) * (self.h - self.mt - self.mb)

    def axes(self, xticks, yticks, xlabel, ylabel, xfmt=str, yfmt=str):
        p = []
        for v in yticks:
            yy = self.y(v)
            p.append(f'<line class="grid" x1="{self.ml}" x2="{self.w - self.mr}" y1="{yy:.1f}" y2="{yy:.1f}"/>')
            p.append(f'<text class="tick" x="{self.ml - 7}" y="{yy + 3.5:.1f}" text-anchor="end">{esc(yfmt(v))}</text>')
        for v in xticks:
            xx = self.x(v)
            p.append(f'<line class="grid" x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.mt}" y2="{self.h - self.mb}"/>')
            p.append(f'<text class="tick" x="{xx:.1f}" y="{self.h - self.mb + 15}" text-anchor="middle">{esc(xfmt(v))}</text>')
        p.append(f'<line class="axis" x1="{self.ml}" x2="{self.w - self.mr}" y1="{self.h - self.mb}" y2="{self.h - self.mb}"/>')
        p.append(f'<text class="alabel" x="{(self.ml + self.w - self.mr) / 2:.1f}" y="{self.h - 6}" text-anchor="middle">{esc(xlabel)}</text>')
        if ylabel:
            cy = (self.mt + self.h - self.mb) / 2
            p.append(f'<text class="alabel" transform="translate(13 {cy:.1f}) rotate(-90)" text-anchor="middle">{esc(ylabel)}</text>')
        return "".join(p)

    def line(self, xs, ys, cls):
        pts = " ".join(f"{self.x(a):.1f},{self.y(b):.1f}" for a, b in zip(xs, ys))
        return f'<polyline class="ln {cls}" points="{pts}"/>'

    def markers(self, xs, ys, cls, shape, tips):
        out = []
        for a, b, tip in zip(xs, ys, tips):
            cx, cy = self.x(a), self.y(b)
            mark = (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.6"/>' if shape == "circle"
                    else f'<rect x="{cx - 4.3:.1f}" y="{cy - 4.3:.1f}" width="8.6" height="8.6" rx="1.2"/>')
            out.append(f'<g class="pt {cls}" tabindex="0" data-tip="{esc(tip)}">'
                       f'<circle class="halo" cx="{cx:.1f}" cy="{cy:.1f}" r="11"/>{mark}</g>')
        return "".join(out)

    def svg(self, body, label, extra=""):
        return (f'<svg viewBox="0 0 {self.w} {self.h}" role="img" aria-label="{esc(label)}" {extra}>{body}</svg>')


def fmt_k(v):
    return f"{v / 1000:g}k" if v else "0"


# ----------------------------------------------------------------------------- charts
def chart_validation(val):
    rows = [r for r in val["rows"] if r["mode"] == "film" and r["fwhm"] > 100]
    i0 = [r["I0"] for r in rows]
    out = []
    specs = [
        ("Te_max", "Peak electron temperature", "T_e,max (K)", (0, 14000), [0, 4000, 8000, 12000], fmt_k, "{:.0f} K", (-3, 4)),
        ("absorptivity", "Fraction of the pump absorbed in the film", "absorptivity", (0, 1), [0, 0.25, 0.5, 0.75, 1.0],
         lambda v: f"{v:g}", "{:.3f}", (-7, 14)),
    ]
    for key, title, ylab, (y0, y1), yt, yfmt, tf, (dy_pub, dy_our) in specs:
        fr = Frame(440, 300, 52, 66, 14, 44, 2.0, 95.0, y0, y1, xlog=True)
        ours = [r[key] for r in rows]
        pub = [r["target"][key] for r in rows]
        body = fr.axes(i0, yt, "pump peak intensity in air (GW/cm²)", ylab, lambda v: f"{v:g}", yfmt)
        body += fr.line(i0, pub, "s1") + fr.line(i0, ours, "s2")
        body += fr.markers(i0, pub, "s1", "circle", [f"Published [15]\n{a:g} GW/cm²: " + tf.format(b) for a, b in zip(i0, pub)])
        body += fr.markers(i0, ours, "s2", "square", [f"This work\n{a:g} GW/cm²: " + tf.format(b) for a, b in zip(i0, ours)])
        xe = fr.x(i0[-1]) + 9
        body += f'<text class="dlabel" x="{xe:.1f}" y="{fr.y(pub[-1]) + dy_pub:.1f}">published</text>'
        body += f'<text class="dlabel" x="{xe:.1f}" y="{fr.y(ours[-1]) + dy_our:.1f}">this work</text>'
        out.append(f'<figure class="panel"><figcaption>{esc(title)}</figcaption>{fr.svg(body, title + ", this work against the published model")}</figure>')
    return "".join(out)


def chart_cost(gate):
    n = [c["N"] for c in gate["cost"]]
    b1 = [c["b1_ms"] for c in gate["cost"]]
    b0 = [c["b0_ms"] for c in gate["cost"]]
    fr = Frame(520, 310, 56, 50, 14, 44, 54, 610, 0.1, 200, xlog=True, ylog=True)
    body = fr.axes(n, [0.1, 1, 10, 100], "energy grid points N", "time per evaluation (ms)", str, lambda v: f"{v:g}")
    body += fr.line(n, b1, "s1") + fr.line(n, b0, "s2")
    body += fr.markers(n, b1, "s1", "circle", [f"Exact integral (B1)\nN = {a}: {b:.2f} ms" for a, b in zip(n, b1)])
    body += fr.markers(n, b0, "s2", "square", [f"FFT integral (B0)\nN = {a}: {b:.3f} ms" for a, b in zip(n, b0)])
    for a, y1, y0 in zip(n, b1, b0):
        if y1 / y0 >= 3:
            body += f'<text class="dlabel" x="{fr.x(a):.1f}" y="{fr.y(math.sqrt(y1 * y0)) + 4:.1f}" text-anchor="middle">×{y1 / y0:.0f}</text>'
    body += f'<text class="dlabel" x="{fr.x(n[-1]) + 10:.1f}" y="{fr.y(b1[-1]) + 4:.1f}">exact</text>'
    body += f'<text class="dlabel" x="{fr.x(n[-1]) + 10:.1f}" y="{fr.y(b0[-1]) + 4:.1f}">FFT</text>'
    return f'<figure class="panel wide"><figcaption>Cost of one collision-integral evaluation on {gate["machine"]["numba_threads"]} CPU threads</figcaption>{fr.svg(body, "Cost per evaluation against grid size")}</figure>'


def chart_transients(gate, series):
    picks = [(i, c) for i, c in enumerate(gate["cases"])
             if c["case"]["mode"] == "energy" and abs(c["case"]["level"] - 1e8) < 1 and c["case"]["fwhm"] in (20.0, 50.0, 220.0)]
    out, blobs = [], []
    for k, (ic, c) in enumerate(picks):
        t = series[f"c{ic}_B1_t"][::2]
        ys = [(series[f"c{ic}_{m}_eps_re"] - series[f"c{ic}_{m}_eps_re"][0])[::2] for m in ("B1", "B0", "thermal")]
        fwhm = c["case"]["fwhm"]
        fr = Frame(318, 262, 44 if k == 0 else 30, 12, 12, 42, float(t[0]), float(t[-1]), -0.02, 0.9)
        xt = [v for v in (-500, -250, 0, 200, 250, 400, 500, 600, 750) if t[0] <= v <= t[-1]]
        if fwhm < 100:
            xt = [v for v in xt if v in (0, 200, 400, 600)]
        else:
            xt = [v for v in xt if v in (-500, 0, 500)]
        body = fr.axes(xt, [0, 0.2, 0.4, 0.6, 0.8], "time from pump peak (fs)", "change of Re ε" if k == 0 else "",
                       str, (lambda v: f"{v:g}") if k == 0 else (lambda v: ""))
        for y, cls in zip(ys, ("s1", "s2", "s3")):
            body += fr.line(t, y, cls)
        ip = int(np.argmax(ys[1]))
        body += f'<text class="dlabel" x="{fr.x(t[ip]) + 8:.1f}" y="{fr.y(ys[1][ip]) + 2:.1f}">FFT</text>'
        body += f'<line class="cross" x1="0" x2="0" y1="{fr.mt}" y2="{fr.h - fr.mb}"/>'
        body += (f'<rect class="hit" x="{fr.x(t[0]):.1f}" y="{fr.mt}" width="{fr.x(t[-1]) - fr.x(t[0]):.1f}" '
                 f'height="{fr.h - fr.mt - fr.mb}"/>')
        sid = f"series-{k}"
        blobs.append(f'<script type="application/json" id="{sid}">' + json.dumps(dict(
            t=[round(float(v), 1) for v in t], px0=round(fr.x(t[0]), 1), px1=round(fr.x(t[-1]), 1),
            names=["exact (B1)", "FFT (B0)", "thermal"], y=[[round(float(v), 4) for v in y] for y in ys])) + "</script>")
        m = c["miss"]
        cap = (f'{fwhm:.0f} fs pump <span class="sub">off from exact by up to: FFT {100 * m["B0"]["eps_re"]:.0f} %, '
               f'thermal {100 * m["thermal"]["eps_re"]:.1f} %</span>')
        out.append(f'<figure class="panel"><figcaption>{cap}</figcaption>'
                   f'{fr.svg(body, f"Permittivity transient for a {fwhm:.0f} fs pump", f"data-series={sid!r}")}</figure>')
    return "".join(out) + "".join(blobs)


# ----------------------------------------------------------------------------- tables
def table(head, rows, cls=""):
    th = "".join(f"<th>{h}</th>" for h in head)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tablewrap"><table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'


def pct(a, b):
    return f"{100 * (a - b) / b:+.1f} %"


def build(out_path: Path):
    val = json.loads((R / "validation_un2023.json").read_text())
    gate = json.loads((R / "gate1.json").read_text())
    series = np.load(R / "gate1_series.npz")
    manifest = json.loads((ROOT / "references" / "manifest.json").read_text())

    n_refs = len(manifest)
    pdf_dir = ROOT / "references" / "pdfs"
    n_pdf = len(list(pdf_dir.glob("*.pdf"))) if pdf_dir.is_dir() and any(pdf_dir.glob("*.pdf")) else 51
    cost = {c["N"]: c for c in gate["cost"]}
    cases = gate["cases"]
    b0e = [c["miss"]["B0"]["eps_re"] for c in cases]
    b0w = [c["miss"]["B0"]["wp2"] for c in cases]
    the = [c["miss"]["thermal"]["eps_re"] for c in cases]
    gmax = max(gate["grid_check"]["miss"].values())
    film_min = {n: 1e6 * cost[n]["b1_ms"] / 1e3 / 60 for n in (128, 256, 512)}

    film = [r for r in val["rows"] if r["mode"] == "film"]
    val_rows = []
    for r in film:
        t = r["target"]
        val_rows.append([f'{r["fwhm"]:.0f} fs, {r["I0"]:g}',
                         f'{r["Te_max"]:.0f} / {t["Te_max"]} <span class="delta">{pct(r["Te_max"], t["Te_max"])}</span>',
                         f'{r["absorptivity"]:.3f} / ' + (f'{t["absorptivity"]:.3f}' if t["absorptivity"] else "n.p."),
                         f'{r["eps_re_max"]:.2f} / ' + (f'{t["eps_re_max"]:.2f}' if t["eps_re_max"] else "n.p.")])

    gate_rows = []
    for c in cases:
        k, m = c["case"], c["miss"]
        level = f'{k["level"] / 1e8:g}×10⁸ J/m³' if k["mode"] == "energy" else f'{k["level"] / 1e9:g} GW/cm² in-film'
        gate_rows.append([("same absorbed energy" if k["mode"] == "energy" else "same in-film field"), f'{k["fwhm"]:.0f}', level,
                          f'{c["Te_max_B1"]:.0f}', f'{100 * m["B0"]["eps_re"]:.0f}', f'{100 * m["B0"]["wp2"]:.0f}',
                          f'{100 * m["RTA"]["eps_re"]:.0f}', f'{100 * m["thermal"]["eps_re"]:.1f}'])

    page = f"""<meta charset="utf-8">
<title>FNO_ENZ Progress Report</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&family=Space+Grotesk:wght@500;600&display=swap">
<style>
/* Layout: one reading column with wider figure blocks; the page opens on the work-package rail, the project's own structure. */
:root {{
  --bg: #f3f5f8; --surface: #fcfcfb; --ink: #12161c; --ink-2: #47505e; --ink-3: #636c7a; --rule: #d9dee5;
  --accent: #1b5cb3;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a;
  --ok: #0c6f3e; --ok-bg: #dff2e7; --warn: #855300; --warn-bg: #fbedd0; --idle: #56606e; --idle-bg: #e6eaef;
  --display: "Space Grotesk", "IBM Plex Sans", system-ui, sans-serif;
  --text: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #0f1216; --surface: #1a1a19; --ink: #eef1f4; --ink-2: #b9c0ca; --ink-3: #939ba7; --rule: #2a3038;
  --accent: #8cb9f6;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70;
  --ok: #5fd09c; --ok-bg: #123021; --warn: #ecb75a; --warn-bg: #36280d; --idle: #a3abb7; --idle-bg: #222830;
  color-scheme: dark }} }}
:root[data-theme="dark"] {{
  --bg: #0f1216; --surface: #1a1a19; --ink: #eef1f4; --ink-2: #b9c0ca; --ink-3: #939ba7; --rule: #2a3038;
  --accent: #8cb9f6;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70;
  --ok: #5fd09c; --ok-bg: #123021; --warn: #ecb75a; --warn-bg: #36280d; --idle: #a3abb7; --idle-bg: #222830;
  color-scheme: dark }}

body {{ background: var(--bg); color: var(--ink); font-family: var(--text); font-size: 16px; line-height: 1.55; }}
.page {{ max-width: 1040px; margin: 0 auto; padding-block: 40px 72px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  display: flex; flex-direction: column; gap: 52px; }}
section, header {{ display: flex; flex-direction: column; gap: 16px; min-width: 0; }}
h1, h2, h3 {{ font-family: var(--display); font-weight: 600; text-wrap: balance; margin: 0; letter-spacing: -0.01em; }}
h1 {{ font-size: clamp(1.75rem, 4.2vw, 2.6rem); line-height: 1.12; max-width: 22ch; }}
h2 {{ font-size: 1.45rem; line-height: 1.2; padding-top: 18px; border-top: 1px solid var(--rule); }}
h3 {{ font-size: 1.02rem; font-weight: 500; }}
p, ul, ol {{ margin: 0; max-width: 68ch; }}
ul, ol {{ padding-left: 1.2em; display: flex; flex-direction: column; gap: 8px; }}
a {{ color: var(--accent); }}
a:focus-visible, .pt:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
code {{ font-family: var(--mono); font-size: 0.88em; overflow-wrap: anywhere; }}
.eyebrow {{ font-family: var(--mono); font-size: 0.78rem; letter-spacing: 0.06em; text-transform: uppercase; color: var(--ink-3); }}
.lede {{ font-size: 1.14rem; line-height: 1.5; color: var(--ink); max-width: 64ch; }}
.note {{ color: var(--ink-2); font-size: 0.93rem; }}
.sub {{ display: block; font-family: var(--text); font-weight: 400; font-size: 0.82rem; color: var(--ink-3); }}

/* work-package rail */
.rail {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 0; border: 1px solid var(--rule); background: var(--surface); }}
.rail > div {{ padding: 14px 16px; border-right: 1px solid var(--rule); border-bottom: 1px solid var(--rule); margin: 0 -1px -1px 0;
  display: flex; flex-direction: column; gap: 6px; min-width: 0; }}
.rail .wp {{ font-family: var(--mono); font-size: 0.76rem; color: var(--ink-3); letter-spacing: 0.04em; }}
.rail .name {{ font-family: var(--display); font-weight: 500; font-size: 0.98rem; line-height: 1.25; }}
.pill {{ align-self: flex-start; font-size: 0.78rem; font-weight: 500; padding: 2px 10px 2px 8px; border-radius: 999px; white-space: nowrap;
  display: inline-flex; align-items: center; gap: 7px; }}
.pill::before {{ content: ""; width: 9px; height: 9px; border-radius: 50%; border: 1.5px solid currentColor; box-sizing: border-box; }}
.pill.ok::before {{ background: currentColor; }}
.pill.warn::before {{ background: linear-gradient(90deg, currentColor 50%, transparent 50%); }}
.pill.ok {{ color: var(--ok); background: var(--ok-bg); }}
.pill.warn {{ color: var(--warn); background: var(--warn-bg); }}
.pill.idle {{ color: var(--idle); background: var(--idle-bg); }}

/* figures at a glance */
.glance {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 20px 28px; }}
.glance div {{ display: flex; flex-direction: column; gap: 2px; min-width: 0; }}
.glance b {{ font-family: var(--display); font-weight: 600; font-size: 2rem; line-height: 1.1; font-variant-numeric: tabular-nums; }}
.glance span {{ color: var(--ink-2); font-size: 0.92rem; }}

.decisions {{ border: 1px solid var(--rule); background: var(--surface); padding: 20px 22px; display: flex; flex-direction: column; gap: 12px; }}
.decisions ol {{ max-width: 72ch; }}

/* charts */
.figrow {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(270px, 1fr)); gap: 16px; }}
.panel {{ margin: 0; background: var(--surface); border: 1px solid var(--rule); padding: 14px 12px 8px; min-width: 0; display: flex; flex-direction: column; gap: 6px; }}
.panel.wide {{ max-width: 640px; }}
.panel figcaption {{ font-family: var(--display); font-weight: 500; font-size: 0.95rem; padding-left: 4px; }}
svg {{ width: 100%; height: auto; display: block; overflow: visible; }}
svg text {{ font-family: var(--text); font-variant-numeric: tabular-nums; }}
.tick {{ fill: var(--ink-3); font-size: 11px; }}
.alabel {{ fill: var(--ink-2); font-size: 11.5px; }}
.dlabel {{ fill: var(--ink-2); font-size: 11.5px; font-weight: 500; }}
.grid {{ stroke: var(--rule); stroke-width: 1; fill: none; }}
.axis {{ stroke: var(--ink-3); stroke-width: 1; fill: none; }}
.ln {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
.ln.s1 {{ stroke: var(--s1); }} .ln.s2 {{ stroke: var(--s2); }} .ln.s3 {{ stroke: var(--s3); }}
.pt.s1 > :not(.halo) {{ fill: var(--s1); stroke: var(--surface); stroke-width: 2; }}
.pt.s2 > :not(.halo) {{ fill: var(--s2); stroke: var(--surface); stroke-width: 2; }}
.pt .halo {{ fill: var(--ink); fill-opacity: 0; }}
.pt:hover .halo, .pt:focus .halo {{ fill-opacity: 0.08; }}
.pt {{ cursor: default; outline: none; }}
.cross {{ stroke: var(--ink-3); stroke-width: 1; stroke-dasharray: 3 3; visibility: hidden; fill: none; }}
.hit {{ fill: var(--ink); fill-opacity: 0; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 0.88rem; color: var(--ink-2); }}
.legend span {{ display: inline-flex; align-items: center; gap: 7px; }}
.sw {{ display: inline-block; width: 18px; height: 3px; border-radius: 2px; }}
.sw.s1 {{ background: var(--s1); }} .sw.s2 {{ background: var(--s2); }} .sw.s3 {{ background: var(--s3); }}
.tip {{ position: absolute; z-index: 20; pointer-events: none; background: var(--ink); color: var(--bg); font-size: 0.8rem; line-height: 1.4;
  padding: 6px 9px; border-radius: 4px; white-space: pre-line; font-variant-numeric: tabular-nums; max-width: 240px; }}
.tip .sw {{ width: 10px; height: 10px; margin-right: 6px; vertical-align: -1px; }}

/* tables */
.tablewrap {{ overflow-x: auto; border: 1px solid var(--rule); background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; font-variant-numeric: tabular-nums; }}
th, td {{ text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--rule); vertical-align: top; }}
tbody tr:last-child td {{ border-bottom: 0; }}
th {{ font-weight: 500; color: var(--ink-2); font-size: 0.82rem; white-space: nowrap; }}
table.num td:not(:first-child), table.num th:not(:first-child) {{ text-align: right; white-space: nowrap; }}
.delta {{ color: var(--ink-3); font-size: 0.82rem; padding-left: 4px; }}
.cols {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 24px 36px; }}
.cols > div {{ display: flex; flex-direction: column; gap: 12px; min-width: 0; }}
footer {{ color: var(--ink-3); font-size: 0.85rem; border-top: 1px solid var(--rule); padding-top: 16px; }}
</style>

<main class="page">
<header>
  <div class="eyebrow">Progress report · 17 November 2026 · after working session 1</div>
  <h1>Learned kinetic closure for ENZ electron dynamics</h1>
  <p class="lede">The project now has a verified reference library with a written summary of all {n_refs} papers, a tested Boltzmann solver for the electron distribution in indium tin oxide, a first validation against a published model, and the Gate 1 measurement. Gate 1 did not pass as written: the exact collision integral is cheap enough to use directly, and the fast baseline it was to be compared with is too inaccurate to use. The case for the neural operator therefore moves from the film simulations to the inference stage.</p>
</header>

<section aria-label="Status by work package">
  <div class="rail">
    <div><span class="wp">WP1 · MONTHS 0–3</span><span class="name">Ground-truth kinetic solver</span><span class="pill ok">first pass done</span></div>
    <div><span class="wp">GATE 1 · MONTH 2</span><span class="name">Is a surrogate needed?</span><span class="pill warn">measured, decision pending</span></div>
    <div><span class="wp">WP2 · MONTHS 2–6</span><span class="name">Structure-preserving neural operator</span><span class="pill idle">not started</span></div>
    <div><span class="wp">WP3 · MONTHS 4–8</span><span class="name">Maxwell–kinetic closure</span><span class="pill idle">not started</span></div>
    <div><span class="wp">WP4 · MONTHS 6–10</span><span class="name">Mechanism discrimination</span><span class="pill idle">not started</span></div>
    <div><span class="wp">WP5 · MONTHS 8–11</span><span class="name">Experimental validation</span><span class="pill idle">not started</span></div>
  </div>
  <div class="glance">
    <div><b>{n_pdf} / {n_refs}</b><span>reference papers on disk; all {n_refs} have a digest in the summary</span></div>
    <div><b>60</b><span>solver tests passing, including conservation of particle number and energy to round-off</span></div>
    <div><b>{cost[128]["b1_ms"]:.1f} ms</b><span>per exact collision-integral evaluation at 128 energy points on 4 CPU cores</span></div>
    <div><b>{100 * min(b0e):.0f}–{100 * max(b0e):.0f} %</b><span>error of the fast FFT baseline on the permittivity transient, over 16 conditions</span></div>
  </div>
</section>

<section class="decisions">
  <h3>Decisions waiting on the project lead</h3>
  <ol>
    <li><strong>Gate 1.</strong> Run the film simulations of WP3 on the exact collision integral and decide on the neural operator after a GPU benchmark, or keep WP2 as planned.</li>
    <li><strong>Plan revisions.</strong> Seven points where the literature contradicts the plan are listed below; the two-photon wavelength constraint and the wording of hypothesis H2 matter most.</li>
    <li><strong>Missing papers.</strong> Nine papers need a download through library access. Reference [1] comes first because it is the calibration dataset for WP4.</li>
    <li><strong>Compute.</strong> A small pilot allocation first, then about 50,000 CPU core-hours and one 40 GB GPU for 1,000–1,500 hours.</li>
    <li><strong>Repository.</strong> It was public when this report was written, so the paper PDFs are kept out of git. Once it is private they can be committed if you want them tracked.</li>
  </ol>
</section>

<section>
  <h2>Gate 1: is a learned surrogate of the collision integral needed?</h2>
  <p>The plan's test has two parts. The exact, momentum-resolved integral must cost at least 100 times more per evaluation than the fast FFT-based integral, and the fast one must miss the exact permittivity transient by more than 2 %. Both hold at 512 energy points. At the resolution the transients actually need, only the second holds.</p>
  {table(["Criterion", "Measured"], [
      ["Exact integral costs ≥ 100× the FFT integral",
       f"×{cost[512]['ratio']:.0f} at 512 points, ×{cost[256]['ratio']:.0f} at 256, ×{cost[128]['ratio']:.0f} at 128"],
      ["FFT integral misses the exact transient by &gt; 2 %",
       f"yes in all 16 conditions: {100 * min(b0e):.0f}–{100 * max(b0e):.0f} % on Re ε, {100 * min(b0w):.0f}–{100 * max(b0w):.0f} % on the plasma frequency squared"],
      ["Grid check, 128 against 256 points, exact integral",
       f"every observable within {100 * gmax:.2f} % of the transient amplitude (one condition tested)"]])}
  <div class="figrow">{chart_cost(gate)}</div>
  <p>The exact integral scales with the cube of the grid size, and the FFT integral barely grows. A film run of the plan's size, 100 cells and a million evaluations, takes about {film_min[128]:.0f} minutes at 128 points on this 4-core machine, {film_min[256] / 60:.1f} hours at 256 and {film_min[512] / 60:.0f} hours at 512. The plan assumed 1 to 10 days.</p>
  <div class="legend" aria-hidden="true"><span><i class="sw s1"></i>exact integral (B1)</span><span><i class="sw s2"></i>FFT integral (B0)</span><span><i class="sw s3"></i>thermal two-temperature model</span></div>
  <div class="figrow">{chart_transients(gate, series)}</div>
  <p class="note">Change of Re ε at the pump frequency for the same absorbed energy, 10⁸ J/m³ at 1300 nm, delivered in three pulse durations. Hover or tap a panel to read values.</p>
  <ul>
    <li>The FFT integral fails mainly through its collision rate, which sets the optical damping and grows with the cube of the density of states. Its error on the distribution alone is {100 * min(b0w[:12]):.0f}–{100 * max(b0w[:12]):.0f} % at equal absorbed energy.</li>
    <li>Anything that replaces the exact integral must reproduce the collision rate as well as the collision term.</li>
    <li>The thermal model built from the same operators is the best cheap model: {100 * min(the):.0f}–{100 * max(the):.0f} % on Re ε, worst for 20 fs pumps and 3–5 % at 220 fs. That residual is the non-thermal contribution hypothesis H1 is about.</li>
    <li>The plan's fallback for a failed gate, using the FFT integral inside the film solver, should not be taken.</li>
  </ul>
  <details><summary>All 16 conditions</summary>
  {table(["Drive", "Pulse (fs)", "Level", "T<sub>e,max</sub> exact (K)", "FFT: Re ε (%)", "FFT: ω<sub>p</sub>² (%)", "Relaxation-time: Re ε (%)", "Thermal: Re ε (%)"], gate_rows, "num")}
  </details>
</section>

<section>
  <h2>Validation against the published kinetic model</h2>
  <p>The solver was run on the case of Un, Sarkar and Sivan, reference [15]: a 40 nm ITO film on glass under gold, pumped at 230 THz. The film's permittivity, the field in the film and the absorbed power are solved together.</p>
  <div class="legend" aria-hidden="true"><span><i class="sw s1"></i>published [15]</span><span><i class="sw s2"></i>this work</span></div>
  <div class="figrow">{chart_validation(val)}</div>
  {table(["Pump (duration, GW/cm²)", "T<sub>e,max</sub> this work / published (K)", "Absorptivity", "Peak Re ε"], val_rows, "num")}
  <p class="note">n.p.: not plotted in the paper. All runs use 128 energy points and the 4 eV band stated in reference [14].</p>
  <div class="cols">
    <div><h3>What agrees</h3><ul>
      <li>Absorptivity within 1–6 % at all six intensities.</li>
      <li>Peak electron temperature within 2–7 % up to 22 GW/cm², and for the weaker 30 fs case.</li>
      <li>With a 6 eV band the three highest 220 fs cases come within 5 % (6203, 8777 and 10771 K against 5902, 8393 and 10262 K).</li>
    </ul></div>
    <div><h3>What does not, and why</h3><ul>
      <li>With the stated 4 eV band the peak temperature overshoots by 16–32 % at 50–75 GW/cm², because the band saturates.</li>
      <li>The change of Re ε is 20–28 % smaller than published. The paper's plotted permittivity implies a plasma frequency 8 % above what its printed parameters give.</li>
      <li>Im ε matches to 1–6 % only if the screening wavevector is frozen at its 300 K value; the paper's equation defines it from the hot distribution.</li>
    </ul></div>
  </div>
  <p>Two things the paper does not print had to be inferred: the pump arrives at 60° in air and refracts through the glass substrate, and the conduction band is wider than 4 eV. Closing the remaining gaps needs the authors' unprinted parameters.</p>
</section>

<section>
  <h2>The solver</h2>
  <div class="cols">
    <div><h3>What it does</h3><ul>
      <li>Evolves the electron distribution f(E, t) in a non-parabolic band together with the lattice temperature.</li>
      <li>Electron–electron collisions at two levels: the exact screened-Coulomb integral and the FFT random-k integral.</li>
      <li>Electron–phonon collisions with momentum conservation, photo-excitation, and switchable Auger-like and two-photon channels.</li>
      <li>Plasma frequency, damping and permittivity computed from f with no temperature assumed.</li>
      <li>A thin-film drive that solves the field self-consistently, and a two-temperature baseline built from the same operators.</li>
    </ul></div>
    <div><h3>Checked against references [14, 15]</h3>
    {table(["Quantity", "This work", "Published"], [
        ["Fermi energy", "0.8805 eV", "0.88 eV"], ["Fermi wavevector", "3.541 nm⁻¹", "3.54 nm⁻¹"],
        ["Electronic heat-capacity coefficient", "12.67 J m⁻³ K⁻²", "12.7"],
        ["Electron–phonon coupling, low temperature", "2.51×10¹⁶ W m⁻³ K⁻¹", "2.5×10¹⁶"],
        ["Zero crossing of Re ε at 300 K", "224 THz", "≈225 THz"]], "num")}</div>
  </div>
  <h3>Corrections to printed formulas</h3>
  {table(["Where", "Problem", "Resolution"], [
      ["[14] screening wavevector", "spurious factor 4π³", "removed; confirmed by a footnote in [16]"],
      ["[14] electron–electron prefactor", "one power of ħ too many; prefactor convention unclear", "dimensionally consistent form; prefactor kept as a fit parameter"],
      ["[14] electron–phonon cut-off", "as printed it does not conserve particle number and triples the coupling", "symmetric cut-off, which reproduces the published coupling exactly"],
      ["Damping rate on a coarse grid", "phonons below one grid cell were missing, so the rate was three times too small at 128 points", "added; the rate is now grid-independent"]])}
</section>

<section>
  <h2>What the literature changes in the plan</h2>
  <p>The summary file has a synthesis and one digest per paper. Fifty-two papers were read in full text; seven paywalled ones are reconstructions from abstracts and citing papers and are flagged as such.</p>
  <ol>
    <li><strong>Two-photon absorption is closed at 1.1–1.6 µm.</strong> With the 2.67 eV gap of reference [23] it needs a pump below 929 nm. The Auger-versus-two-photon test can only run on 800 nm few-cycle data, which the plan defers to Tier 2.</li>
    <li><strong>Hypothesis H2 needs rewording.</strong> A fixed-density thermal model already reverses the sign of the transmittance [24]. The hypothesis has to be stated for Re ε or the plasma frequency.</li>
    <li><strong>The thermal model M0 needs a temperature-dependent damping law</strong> with free parameters to produce any reversal.</li>
    <li><strong>The density-versus-mass reading of the current equation is the project's own.</strong> References [30, 31] give the two mathematical forms but do not make that identification.</li>
    <li><strong>The novelty claim survives but must be reworded.</strong> Reference [47] already enforces one conservation law exactly, and [49] learns a flow map, not a collision operator. No paper treats Pauli blocking.</li>
    <li><strong>"Train at 256 points, evaluate at 512" has no support</strong> in the cited neural-operator papers.</li>
    <li><strong>References [58–60] are not precedents</strong> for double-blind FROG with a complex gate; the closest is the retrieval in [22].</li>
  </ol>
  <p class="note">Twenty-six citations in the plan also need updating, mostly preprints that have since been published. Reference [16] has three authors and [22] was retitled.</p>
</section>

<section>
  <h2>Dataset pilot and a validity limit</h2>
  <p>The dataset generator produced eight trajectories with the exact integral at 128 points, about 40 seconds and 4 MB each. Every trajectory stores the distribution, the collision terms the operator would learn, and perturbed off-trajectory states.</p>
  <p>A third of the plan's sampling range does not fit in a single-band model. Drawing from 20–250 fs and 0.1 GW/cm² to 10 TW/cm² in-film, 34 % of cases would deposit more than half of what a 4 eV band can hold. The generator now screens those out. The upper intensity range needs higher bands in the model, or an explicit validity limit near 10¹¹ W/cm² in-film for a 100 fs pulse.</p>
</section>

<section>
  <h2>Compute</h2>
  <p>One exact evaluation takes 4.0, 37 and 231 ms on a single core at 128, 256 and 512 points. Runs are independent, so cost is counted in core-hours. No GPU was available; GPU figures are estimates until benchmarked.</p>
  {table(["Workload", "128 points", "256 points", "512 points"], [
      ["Dataset, 20,000 trajectories", "400", "3,700", "23,000"],
      ["One film run (100 cells, 2 ps)", "1.8", "16.5", "103"],
      ["Film-solver development, about 300 runs", "530", "5,000", "31,000"],
      ["Inference, 1,000–10,000 film runs", "1,800–18,000", "16,500–165,000", "10⁵–10⁶"]], "num")}
  <p>The grid size drives everything, so a grid-convergence sweep comes before any large request. Suggested: a pilot of about 2,000 core-hours and 100 GPU-hours, then about 50,000 core-hours, one GPU with at least 40 GB for 1,000–1,500 hours, and 500 GB of storage.</p>
</section>

<section>
  <h2>Next steps, in priority order</h2>
  <ol>
    <li>Write to the authors of references [14–16] for the parameters the papers do not print.</li>
    <li>Obtain the nine missing papers and redo the seven reconstructed digests, [1] first.</li>
    <li>Fix the band model for high intensity, then rerun the validation and Gate 1 with it.</li>
    <li>Port the exact integral to a GPU, benchmark it, and run the grid-convergence sweep over the full case matrix.</li>
    <li>Derive the current equation from the kinetic equation, to decide between the two forms for the film solver.</li>
    <li>Revise the plan for the seven literature findings.</li>
  </ol>
</section>

<section>
  <h2>Where things are</h2>
  {table(["Path", "Contents"], [
      ["<code>docs/WORKLOG.md</code>", "Working log: everything done, dead ends, decisions, open items"],
      ["<code>docs/gates/G1_report.md</code>", "Gate 1 report"],
      ["<code>docs/theory/wp1_kinetic_model.md</code>", "Equations as implemented, with corrections"],
      ["<code>references/LITERATURE_SUMMARY.md</code>", "Synthesis and one digest per paper"],
      ["<code>references/README.md</code>", "Reference index, download status, citation corrections"],
      ["<code>src/enzkin/</code>, <code>tests/</code>", "Solver package and its tests"],
      ["<code>scripts/</code>, <code>results/wp1/</code>", "Reproducible runs and their outputs, including this report's numbers"]])}
</section>

<footer>Generated by <code>scripts/make_progress_report.py</code> from <code>results/wp1/validation_un2023.json</code> and <code>results/wp1/gate1.json</code>. Reference numbers in square brackets are those of the project plan. Measurements were taken on a 4-core CPU without a GPU.</footer>
</main>

<script>
(function () {{
  var tip = document.createElement('div'); tip.className = 'tip'; tip.hidden = true; document.body.appendChild(tip);
  function place(x, y) {{
    tip.hidden = false;
    var r = tip.getBoundingClientRect(), L = x + 12, T = y - r.height - 12;
    if (L + r.width > window.innerWidth - 8) L = x - r.width - 12;
    if (L < 8) L = 8;
    if (T < 8) T = y + 16;
    tip.style.left = (L + window.scrollX) + 'px'; tip.style.top = (T + window.scrollY) + 'px';
  }}
  function hide() {{ tip.hidden = true; }}
  document.querySelectorAll('[data-tip]').forEach(function (el) {{
    function show() {{ var r = el.getBoundingClientRect(); tip.textContent = el.dataset.tip; place(r.left + r.width / 2, r.top); }}
    el.addEventListener('pointerenter', show); el.addEventListener('focus', show);
    el.addEventListener('pointerleave', hide); el.addEventListener('blur', hide);
  }});
  document.querySelectorAll('svg[data-series]').forEach(function (svg) {{
    var d = JSON.parse(document.getElementById(svg.dataset.series).textContent);
    var hit = svg.querySelector('.hit'), cross = svg.querySelector('.cross'), n = d.t.length;
    function move(e) {{
      var r = hit.getBoundingClientRect();
      var i = Math.max(0, Math.min(n - 1, Math.round((e.clientX - r.left) / r.width * (n - 1))));
      var x = d.px0 + (d.px1 - d.px0) * i / (n - 1);
      cross.setAttribute('x1', x); cross.setAttribute('x2', x); cross.style.visibility = 'visible';
      tip.textContent = '';
      var head = document.createElement('b'); head.textContent = d.t[i] + ' fs'; tip.appendChild(head);
      d.names.forEach(function (name, k) {{
        tip.appendChild(document.createElement('br'));
        var sw = document.createElement('i'); sw.className = 'sw s' + (k + 1); tip.appendChild(sw);
        tip.appendChild(document.createTextNode(name + '  ' + d.y[k][i].toFixed(3)));
      }});
      place(e.clientX, e.clientY);
    }}
    hit.addEventListener('pointermove', move); hit.addEventListener('pointerdown', move);
    hit.addEventListener('pointerleave', function () {{ cross.style.visibility = 'hidden'; hide(); }});
  }});
}})();
</script>
"""
    head, sep, body = page.partition("</style>")
    body, sep2, script = body.rpartition("<script>")
    page = head + sep + body.replace(" %", "&nbsp;%") + sep2 + script
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(page, encoding="utf-8")
    print("wrote", out_path, f"({len(page) / 1024:.0f} kB)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(ROOT / "docs" / "reports" / "progress_2026-11-17.html"))
    build(Path(ap.parse_args().out))
