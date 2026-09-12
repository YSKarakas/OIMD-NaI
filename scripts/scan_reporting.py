#!/usr/bin/env python3
"""Corpus-wide reporting practice, and who already corrects by a simulated LCE.

Two questions, both answered by keyword over the full corpus rather than by
reading, so that they cover all 132 papers instead of the handful extracted by
hand.

  1. How much of what a light yield depends on does the literature mention at
     all -- surface treatment, reflector, light-collection efficiency?
  2. Which papers already do what this project proposed as its method: divide a
     measured photoelectron count by a light-collection efficiency obtained from
     an optical simulation?

READ EVERY NUMBER HERE AS AN UPPER BOUND.  A paper that contains the word
"polished" anywhere -- in its introduction, about someone else's crystal, in a
reference title -- counts as mentioning a surface treatment. Stating a surface
treatment FOR THE SAMPLE MEASURED is a stricter thing, and the hand-extracted
records in ly_records.csv are what measure that. The two are reported together
precisely because the gap between them is informative: the keyword scan says
what the literature could be doing, the extraction says what it does.

    python3 scripts/scan_reporting.py
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "literature" / "cache"
CLASSES = ROOT / "data" / "literature" / "fulltext_classes.csv"
OUT = ROOT / "data" / "literature" / "reporting_scan.csv"

PATTERNS = {
    "surface": re.compile(r"\bpolish|\bground\b|lapp?ed|saw[- ]cut|etch", re.I),
    "wrapping": re.compile(
        r"teflon|PTFE|ESR|VM2000|tyvek|lumirror|BaSO4|TiO2|Gore|millipore", re.I),
    "lce": re.compile(r"light[- ]collection efficiency|\bLCE\b", re.I),
    "optical_sim": re.compile(
        r"\bGeant4\b|\bGATE\b|\bLitrani\b|\bZemax\b|Monte[- ]Carlo", re.I),
}


def main() -> None:
    meta = {p["arxiv_id"].split("v")[0]: p
            for p in json.loads((CACHE / "candidates.json").read_text())["papers"]}
    with CLASSES.open() as fh:
        klass = {r["arxiv_id"]: r["klass"]
                 for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))}

    rows = []
    for path in sorted((CACHE / "txt").glob("*.txt")):
        text = path.read_text(errors="replace")
        flags = {name: bool(rx.search(text)) for name, rx in PATTERNS.items()}
        rows.append({
            "arxiv_id": path.stem,
            "klass": klass.get(path.stem, "?"),
            "title": meta.get(path.stem.split("v")[0], {}).get("title", ""),
            **{k: str(v).lower() for k, v in flags.items()},
            "lce_corrected": str(flags["lce"] and flags["optical_sim"]).lower(),
        })

    with OUT.open("w", newline="") as fh:
        fh.write("# Keyword scan of reporting practice. UPPER BOUNDS: a paper counts as\n")
        fh.write("# mentioning something if the word appears anywhere, including about\n")
        fh.write("# other people's crystals. Produced by scripts/scan_reporting.py.\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    absolute = [r for r in rows if r["klass"] == "absolute"]
    print(f"full texts: {len(rows)}   reporting an absolute yield: {len(absolute)}\n")
    labels = [
        ("surface", "mentions any surface treatment"),
        ("wrapping", "mentions any reflector or wrapping"),
        ("lce", "mentions light-collection efficiency"),
        ("optical_sim", "mentions an optical simulation tool"),
        ("lce_corrected", "both: plausibly corrects by a simulated LCE"),
    ]
    for name, subset in (("whole corpus", rows), ("absolute-yield papers", absolute)):
        print(f"--- {name} (n={len(subset)}) ---")
        for key, description in labels:
            n = sum(1 for r in subset if r[key] == "true")
            print(f"   {description:<44} {n:>4} / {len(subset)}   "
                  f"{100 * n / len(subset):.0f} %")
        print()

    print("Absolute-yield papers that already correct by a simulated LCE")
    print("(this is the prior art for the correction protocol):")
    for r in absolute:
        if r["lce_corrected"] == "true":
            print(f"   {r['arxiv_id']:16s} {r['title'][:62]}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
