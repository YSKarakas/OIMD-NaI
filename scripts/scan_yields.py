#!/usr/bin/env python3
"""Classify each cached paper by what KIND of light yield it reports.

This is the full-text screening step, and it is mechanical on purpose: it sorts
papers by the form of the number they publish, which is a property of the text,
not a judgement about the physics. The judgement -- what the number means and
under what conditions it was taken -- happens afterwards, by reading.

The taxonomy is itself one of the survey's results, because these are not
interchangeable quantities and the literature routinely presents them as if they
were:

  absolute    photons/MeV          a property of the material (plus geometry)
  detected    photoelectrons/MeV   a property of the material AND the
                                   photodetector's quantum efficiency AND the
                                   light collection of that particular setup.
                                   Comparing one of these to a photons/MeV
                                   number is comparing different quantities.
  relative    "% of NaI(Tl)",      meaningful only if the reference is stated
              "relative to ..."    together with the conditions it was taken in
  none        no light yield in numerical form

A paper can fall in more than one class; it is filed under the strongest it
reports (absolute > detected > relative > none).

    python3 scripts/scan_yields.py                # summary table
    python3 scripts/scan_yields.py --class absolute --sentences
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "literature" / "cache"
CLASSES = ROOT / "data" / "literature" / "fulltext_classes.csv"

# A number, then within a short distance a unit. Requiring the NUMBER is what
# separates "we measured 41000 photons/MeV" from "the light yield is important".
_NUM = r"(\d[\d\s,.]{0,12}\d|\d)"
_GAP = r"[^.\n]{0,25}?"

PHOTONS = re.compile(
    rf"{_NUM}\s*(?:[x×]\s*10\^?\d+\s*)?{_GAP}"
    r"(photons?\s*(?:/|per)\s*MeV|ph\s*/\s*MeV|photons?\s*MeV\s*-?\s*1)", re.I)
PHOTOELECTRONS = re.compile(
    rf"{_NUM}{_GAP}"
    r"(photo-?electrons?\s*(?:/|per)\s*MeV|phe\s*/\s*MeV|p\.?e\.?\s*/\s*MeV"
    r"|photo-?electrons?\s*(?:/|per)\s*keV|phe\s*/\s*keV)", re.I)
PER_KEV = re.compile(
    rf"{_NUM}{_GAP}(photons?\s*(?:/|per)\s*keV|ph\s*/\s*keV)", re.I)
RELATIVE = re.compile(
    r"(relative light (yield|output)|light (yield|output) relative to"
    r"|\d{1,3}\s*%\s*(of|relative to)\s*(the\s*)?(NaI|CsI|BGO|LYSO|LSO|YAP|reference)"
    r"|normali[sz]ed to (the )?(NaI|CsI|BGO|YAP|reference))", re.I)


def sentences(text: str, pattern: re.Pattern, limit: int = 6) -> list[str]:
    out = []
    for match in pattern.finditer(text):
        start = max(0, match.start() - 160)
        end = min(len(text), match.end() + 80)
        out.append(" ".join(text[start:end].split()))
        if len(out) >= limit:
            break
    return out


def classify(text: str) -> tuple[str, dict[str, int]]:
    counts = {
        "photons_per_mev": len(PHOTONS.findall(text)),
        "photons_per_kev": len(PER_KEV.findall(text)),
        "photoelectrons": len(PHOTOELECTRONS.findall(text)),
        "relative": len(RELATIVE.findall(text)),
    }
    if counts["photons_per_mev"] or counts["photons_per_kev"]:
        return "absolute", counts
    if counts["photoelectrons"]:
        return "detected", counts
    if counts["relative"]:
        return "relative", counts
    return "none", counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--class", dest="want", default=None,
                    choices=["absolute", "detected", "relative", "none"])
    ap.add_argument("--sentences", action="store_true")
    args = ap.parse_args()

    meta = {p["arxiv_id"].split("v")[0]: p
            for p in json.loads((CACHE / "candidates.json").read_text())["papers"]}

    rows, tally = [], Counter()
    for path in sorted((CACHE / "txt").glob("*.txt")):
        text = path.read_text(errors="replace")
        verdict, counts = classify(text)
        tally[verdict] += 1
        info = meta.get(path.stem.split("v")[0], {})
        rows.append({
            "arxiv_id": path.stem,
            "year": info.get("published", "")[:4],
            "klass": verdict,
            "title": info.get("title", ""),
            "journal_ref": info.get("journal_ref") or "",
            "doi": info.get("doi") or "",
            **counts,
        })

    CLASSES.parent.mkdir(parents=True, exist_ok=True)
    with CLASSES.open("w", newline="") as fh:
        fh.write("# Full-text screening: what KIND of light yield each paper reports.\n")
        fh.write("# Produced by scripts/scan_yields.py -- mechanical, rebuildable.\n")
        fh.write("# photons/MeV, photoelectrons/MeV and relative yields are different\n")
        fh.write("# quantities; the counts are how many numeric mentions of each appear.\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    if args.want:
        for row in rows:
            if row["klass"] != args.want:
                continue
            print(f"\n{'=' * 74}\n{row['arxiv_id']}  {row['year']}  {row['title'][:70]}")
            if args.sentences:
                text = (CACHE / "txt" / f"{row['arxiv_id']}.txt").read_text(errors="replace")
                for pattern in (PHOTONS, PER_KEV, PHOTOELECTRONS):
                    for sentence in sentences(text, pattern, 4):
                        print(f"   > {sentence[:230]}")
        return

    print(f"{len(rows)} papers with full text\n")
    for name in ("absolute", "detected", "relative", "none"):
        print(f"  {name:10s} {tally[name]:4d}")
    print(f"\nwrote {CLASSES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
