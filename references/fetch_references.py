#!/usr/bin/env python3
"""Download the project's reference library from references/manifest.json.

The PDFs themselves are NOT tracked in git (the repository is public and most papers are not
licensed for redistribution).  This script rebuilds the local library in one command:

    python references/fetch_references.py            # PDFs  -> references/pdfs/NN_key.pdf
    python references/fetch_references.py --sources  # + arXiv TeX sources -> references/sources/NN_key/
    python references/fetch_references.py --text     # + plain text (pdftotext) -> references/text/NN_key.txt
    python references/fetch_references.py --only 14 15 16

Rules
-----
* For each reference the `sources` list in the manifest is tried in order until a real PDF arrives.
* A file that already exists is never re-downloaded, so a PDF you drop in by hand under the
  expected name (e.g. references/pdfs/01_alam2016.pdf for a paywalled paper) is picked up as-is.
* arXiv is contacted at most once every 3 s, as its terms of use ask.
* Needs outbound access to arxiv.org (and tzin.bgu.ac.il for [16]); publisher sites that reject
  scripted clients are reported as "manual" - fetch those in a browser.
"""
from __future__ import annotations

import argparse
import io
import json
import shutil
import subprocess
import sys
import tarfile
import gzip
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PDFS, SRCS, TEXT = HERE / "pdfs", HERE / "sources", HERE / "text"
UA = "Mozilla/5.0 (X11; Linux x86_64) fetch-references/1.0 (research library sync)"
_last_arxiv = 0.0


def _get(url: str, timeout: int = 90) -> bytes | None:
    global _last_arxiv
    if "arxiv.org" in url:
        wait = 3.1 - (time.time() - _last_arxiv)
        if wait > 0:
            time.sleep(wait)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
    except Exception as exc:  # network errors are reported, not fatal
        print(f"      ! {type(exc).__name__}: {str(exc)[:80]}  <- {url}")
        data = None
    if "arxiv.org" in url:
        _last_arxiv = time.time()
    return data


def stem(ref: dict) -> str:
    return f"{ref['n']:02d}_{ref['key']}"


def fetch_pdf(ref: dict) -> str:
    out = PDFS / f"{stem(ref)}.pdf"
    if out.exists() and out.stat().st_size > 10_000:
        return "present"
    for src in ref.get("sources", []):
        if src["kind"] == "abstract-page":
            continue
        data = _get(src["url"])
        if data and data[:5] == b"%PDF-":
            out.write_bytes(data)
            return f"downloaded ({src['kind']})"
    return "manual" if ref.get("sources") else "no-open-source"


def _safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    dest = dest.resolve()
    for m in tar.getmembers():
        target = (dest / m.name).resolve()
        if not (m.isfile() or m.isdir()) or dest not in target.parents and target != dest:
            continue  # skip links, devices and anything escaping dest
        tar.extract(m, dest, filter="data") if hasattr(tarfile, "data_filter") else tar.extract(m, dest)


def fetch_source(ref: dict) -> str:
    ax = ref.get("arxiv")
    if not ax:
        return "-"
    dest = SRCS / stem(ref)
    if dest.exists() and any(dest.iterdir()):
        return "present"
    data = _get(f"https://arxiv.org/e-print/{ax['version']}")
    if not data:
        return "failed"
    if data[:5] == b"%PDF-":
        return "pdf-only submission"
    dest.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as tar:
            _safe_extract(tar, dest)
        return "extracted"
    except tarfile.TarError:
        try:  # single gzipped .tex file
            (dest / "main.tex").write_bytes(gzip.decompress(data))
            return "extracted (single file)"
        except OSError:
            shutil.rmtree(dest, ignore_errors=True)
            return "unknown format"


def make_text(ref: dict) -> str:
    pdf, txt = PDFS / f"{stem(ref)}.pdf", TEXT / f"{stem(ref)}.txt"
    if not pdf.exists():
        return "-"
    if txt.exists():
        return "present"
    if not shutil.which("pdftotext"):
        return "pdftotext missing"
    r = subprocess.run(["pdftotext", "-layout", str(pdf), str(txt)], capture_output=True)
    return "ok" if r.returncode == 0 else "failed"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", action="store_true", help="also fetch arXiv TeX sources")
    ap.add_argument("--text", action="store_true", help="also extract plain text with pdftotext")
    ap.add_argument("--only", type=int, nargs="*", help="reference numbers to process")
    args = ap.parse_args()

    refs = json.loads((HERE / "manifest.json").read_text())
    for d in (PDFS, SRCS, TEXT):
        d.mkdir(exist_ok=True)
    tally: dict[str, int] = {}
    for ref in refs:
        if args.only and ref["n"] not in args.only:
            continue
        status = fetch_pdf(ref)
        line = f"[{ref['n']:2d}] {ref['key']:<24} pdf: {status}"
        if args.sources:
            line += f" | src: {fetch_source(ref)}"
        if args.text:
            line += f" | text: {make_text(ref)}"
        print(line, flush=True)
        tally[status.split(" ")[0]] = tally.get(status.split(" ")[0], 0) + 1
    print("\nsummary:", ", ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
