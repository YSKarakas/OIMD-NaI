#!/usr/bin/env python3
"""Show the parts of a cached paper that carry the survey's metadata.

Reading 149 full texts end to end is not feasible, and skimming is how surveys
acquire errors. This pulls out the lines that mention any of the quantities the
Paper 0 schema asks for, with context, so that the reading is targeted but the
evidence for each extracted field is still the paper's own words.

It deliberately does NOT decide anything. If a field is absent from the output
below, that is a prompt to search the text directly, not a licence to record
`not_stated` -- a paper may phrase something in a way these patterns miss, and
recording a false `not_stated` would corrupt the one statistic being measured.

    python3 scripts/read_paper.py 2403.02668
    python3 scripts/read_paper.py 2403.02668 --full        # whole text
    python3 scripts/read_paper.py --next 5                 # not yet extracted
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "literature" / "cache"
RECORDS = ROOT / "data" / "literature" / "ly_records.csv"

PATTERNS: dict[str, str] = {
    "yield": r"(light yield|light output|photons?\s*/\s*MeV|ph\s*/\s*MeV|phe\s*/\s*MeV"
             r"|photoelectrons?\s*(per|/)\s*MeV|photo-?electron yield)",
    "size": r"(\d+\s*[x×]\s*\d+\s*[x×]\s*\d+|\bmm\s*[x×]|\bcm\s*[x×]|diameter|thickness"
            r"|\bcubic\b|\bcube\b|dimensions?)",
    "surface": r"(polish|ground|lapp?ed|saw[- ]cut|etch|roughen|as[- ]grown|surface (treatment|finish))",
    "wrap": r"(teflon|PTFE|ESR|VM2000|tyvek|lumirror|millipore|BaSO4|TiO2|white paint|wrapp?|reflector)",
    "couple": r"(optical (grease|cement|glue)|silicone|coupl|index matching|RTV|EJ-550|BC-630)",
    "detector": r"(PMT|photomultiplier|SiPM|avalanche|APD|Hamamatsu|Photonis|ET Enterprises|R\d{4}|S\d{4}-)",
    "standard": r"(reference (crystal|standard)|normali[sz]ed to|calibrat|absolute light|YAP|NaI\(?Tl\)?\s*(standard|reference))",
    "source": r"(\^?\{?\d{2,3}\}?\s*(Cs|Co|Na|Am|Ba|Eu|Cd)|662\s*keV|511\s*keV|59\.5\s*keV|122\s*keV|1332\s*keV)",
    "shaping": r"(shaping time|integration (time|gate)|\bµs\b|\bus\b gate|micro-?second|time constant)",
    "temp": r"(room temperature|\bRT\b|\d+\s*K\b|temperature dependence|cooled|cryogenic)",
}


def _text_path(arxiv_id: str) -> Path:
    stem = re.sub(r"[^0-9A-Za-z.\-]", "_", arxiv_id)
    return CACHE / "txt" / f"{stem}.txt"


def _metadata(arxiv_id: str) -> dict | None:
    data = json.loads((CACHE / "candidates.json").read_text())
    for paper in data["papers"]:
        if paper["arxiv_id"].split("v")[0] == arxiv_id.split("v")[0]:
            return paper
    return None


def already_extracted() -> set[str]:
    if not RECORDS.exists():
        return set()
    with RECORDS.open() as fh:
        rows = csv.DictReader(ln for ln in fh if not ln.startswith("#"))
        return {r["arxiv_id"].split("v")[0] for r in rows if r.get("arxiv_id")}


def show(arxiv_id: str, *, full: bool = False, context: int = 0) -> None:
    path = _text_path(arxiv_id)
    if not path.exists():
        matches = sorted(CACHE.glob(f"txt/{arxiv_id.split('v')[0]}*.txt"))
        if not matches:
            print(f"no cached text for {arxiv_id}")
            return
        path = matches[0]

    meta = _metadata(arxiv_id) or {}
    print("=" * 78)
    print(f"{arxiv_id}   {meta.get('published', '')}")
    print(meta.get("title", ""))
    if meta.get("journal_ref"):
        print(f"journal: {meta['journal_ref']}")
    if meta.get("doi"):
        print(f"doi: {meta['doi']}")
    print("=" * 78)

    text = path.read_text(errors="replace")
    if full:
        print(text)
        return

    print("\n--- ABSTRACT ---")
    print(meta.get("summary", "(not in metadata)")[:900])

    lines = text.splitlines()
    for label, pattern in PATTERNS.items():
        rx = re.compile(pattern, re.I)
        hits = [i for i, ln in enumerate(lines) if rx.search(ln)]
        print(f"\n--- {label.upper()}  ({len(hits)} lines) ---")
        if not hits:
            print("  (no match -- check the text directly before recording not_stated)")
            continue
        shown, last = 0, -10
        for i in hits:
            if shown >= 6:
                print(f"  ... {len(hits) - shown} more")
                break
            if i - last > 1 or context:
                for j in range(max(0, i - context), min(len(lines), i + context + 1)):
                    stripped = " ".join(lines[j].split())
                    if stripped:
                        print(f"  {j:5d}| {stripped[:150]}")
                shown += 1
            last = i


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("arxiv_id", nargs="?")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--context", type=int, default=0)
    ap.add_argument("--next", type=int, default=0,
                    help="list this many cached papers not yet in ly_records.csv")
    args = ap.parse_args()

    if args.next:
        done = already_extracted()
        pending = [p.stem for p in sorted((CACHE / "txt").glob("*.txt"))
                   if p.stem.split("v")[0] not in done]
        print(f"{len(pending)} cached papers not yet extracted; next {args.next}:")
        for stem in pending[: args.next]:
            meta = _metadata(stem) or {}
            print(f"  {stem:16s} {meta.get('title', '')[:70]}")
        return

    if not args.arxiv_id:
        ap.error("give an arXiv id, or --next N")
    show(args.arxiv_id, full=args.full, context=args.context)


if __name__ == "__main__":
    main()
