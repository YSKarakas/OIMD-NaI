#!/usr/bin/env python3
"""Assemble the Paper 0 literature corpus from arXiv, reproducibly.

WHAT THIS DOES AND DOES NOT DO.  It gathers candidate papers and their full text.
It does NOT extract the measurement metadata: that is done by reading, and the
results go into data/literature/ly_records.csv by hand.  Automating the
extraction would mean inventing values for fields that papers do not report,
which is the precise thing this survey exists to measure.

THE SEARCH PROTOCOL is fixed here rather than described in prose, so that the
corpus can be rebuilt and the counts checked:

  source      arXiv API (export.arxiv.org/api/query)
  categories  physics.ins-det, cond-mat.mtrl-sci, nucl-ex, physics.med-ph
  queries     see QUERIES below
  inclusion   the paper reports an ORIGINAL light-yield measurement, in
              photons/MeV or photoelectrons/MeV, for an inorganic scintillator
  exclusion   compilations and reviews that only restate other people's numbers;
              organic scintillators and liquid scintillators; simulation-only
              papers with no measurement of their own

KNOWN BIAS, to be stated in the paper.  arXiv is open and scriptable, which is
why it is used, but its coverage is not neutral: instrumentation and
high-energy-physics work is over-represented and materials-science journals
(J. Lumin., Opt. Mater., J. Alloys Compd.) are under-represented, and those are
where much of the crystal-growth literature lives.  The resulting statistic is
therefore "of the arXiv-accessible papers reporting a light yield, X % state the
surface treatment", not a claim about the whole field.

    python3 scripts/collect_literature.py --fetch      # query arXiv, cache metadata
    python3 scripts/collect_literature.py --download   # fetch PDFs and extract text
    python3 scripts/collect_literature.py --status
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "literature" / "cache"
METADATA = CACHE / "candidates.json"
SCREENING = ROOT / "data" / "literature" / "screening.csv"

ARXIV_API = "https://export.arxiv.org/api/query"
# arXiv asks for no more than one request every three seconds.
POLITE_DELAY_S = 3.0
USER_AGENT = "scint-survey/0.1 (academic literature survey; contact via repository)"

ATOM = "{http://www.w3.org/2005/Atom}"

QUERIES: list[str] = [
    'abs:"light yield" AND abs:scintillator',
    'abs:"light output" AND abs:scintillator',
    'abs:"photons/MeV"',
    'abs:"photons per MeV" AND abs:crystal',
    'abs:scintillation AND abs:"energy resolution" AND abs:crystal',
    'abs:"light yield" AND abs:garnet',
    'abs:"light yield" AND abs:"single crystal"',
    'abs:scintillator AND abs:"Teflon"',
    'abs:scintillator AND abs:"wrapping"',
    'abs:"light collection efficiency" AND abs:scintillator',
]

CATEGORIES = ["physics.ins-det", "cond-mat.mtrl-sci", "nucl-ex", "physics.med-ph"]


def _get(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def search(query: str, *, max_results: int = 100) -> list[dict]:
    """One arXiv query, restricted to the categories above."""
    cats = " OR ".join(f"cat:{c}" for c in CATEGORIES)
    full = f"({query}) AND ({cats})"
    params = urllib.parse.urlencode({
        "search_query": full,
        "start": 0,
        "max_results": max_results,
        "sortBy": "relevance",
        "sortOrder": "descending",
    })
    root = ET.fromstring(_get(f"{ARXIV_API}?{params}"))

    out = []
    for entry in root.findall(f"{ATOM}entry"):
        arxiv_id = entry.findtext(f"{ATOM}id", "").rsplit("/", 1)[-1]
        out.append({
            "arxiv_id": arxiv_id,
            "title": " ".join((entry.findtext(f"{ATOM}title") or "").split()),
            "published": entry.findtext(f"{ATOM}published", "")[:10],
            "summary": " ".join((entry.findtext(f"{ATOM}summary") or "").split()),
            "authors": [a.findtext(f"{ATOM}name") for a in entry.findall(f"{ATOM}author")],
            "doi": entry.findtext("{http://arxiv.org/schemas/atom}doi"),
            "journal_ref": entry.findtext("{http://arxiv.org/schemas/atom}journal_ref"),
            "found_by": [query],
        })
    return out


def fetch() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    known: dict[str, dict] = {}
    if METADATA.exists():
        known = {r["arxiv_id"]: r for r in json.loads(METADATA.read_text())["papers"]}

    for i, query in enumerate(QUERIES):
        if i:
            time.sleep(POLITE_DELAY_S)
        try:
            results = search(query)
        except Exception as exc:  # a failed query must not lose the rest
            print(f"  QUERY FAILED  {query}: {exc}")
            continue
        new = 0
        for record in results:
            existing = known.get(record["arxiv_id"])
            if existing is None:
                known[record["arxiv_id"]] = record
                new += 1
            elif query not in existing["found_by"]:
                existing["found_by"].append(query)
        print(f"  {len(results):3d} hits, {new:3d} new   {query}")

    METADATA.write_text(json.dumps({
        "protocol": {
            "source": ARXIV_API,
            "categories": CATEGORIES,
            "queries": QUERIES,
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
        "papers": sorted(known.values(), key=lambda r: r["arxiv_id"]),
    }, indent=2) + "\n")
    print(f"\n{len(known)} unique candidates -> {METADATA.relative_to(ROOT)}")


# --- title/abstract screening ---------------------------------------------- #
#
# A systematic review screens on title and abstract before reading full text.
# The patterns below do the mechanical part of that -- they are deliberately
# crude, and their only job is to sort candidates into buckets that a human then
# confirms. Nothing is excluded on a machine's say-so alone: the screening file
# records the automatic verdict AND a column for the reviewed one, and the
# extraction step reads the reviewed column.

ORGANIC = re.compile(
    r"\b(plastic scintillat|liquid scintillat|organic scintillat|stilbene|anthracene"
    r"|EJ-2|EJ-3|BC-4|polystyrene|PVT\b|water[- ]based)", re.I)
YIELD = re.compile(
    r"(light yield|light output|photons?\s*/\s*MeV|ph/MeV"
    r"|photoelectrons?\s*/\s*MeV|phe/MeV)", re.I)


def auto_bucket(paper: dict) -> tuple[str, str]:
    text = f"{paper['title']} {paper['summary']}"
    if ORGANIC.search(text):
        return "exclude", "organic or liquid scintillator"
    if not YIELD.search(text):
        return "exclude", "abstract reports no light yield or output"
    return "screen", "abstract mentions a light yield -- read the full text"


def screen() -> None:
    """Write the screening sheet. Existing human decisions are preserved."""
    papers = json.loads(METADATA.read_text())["papers"]

    reviewed: dict[str, tuple[str, str]] = {}
    if SCREENING.exists():
        import csv as _csv
        with SCREENING.open() as fh:
            for row in _csv.DictReader(ln for ln in fh if not ln.startswith("#")):
                if row.get("reviewed_decision"):
                    reviewed[row["arxiv_id"]] = (row["reviewed_decision"],
                                                 row.get("reviewed_reason", ""))

    import csv
    SCREENING.parent.mkdir(parents=True, exist_ok=True)
    with SCREENING.open("w", newline="") as fh:
        fh.write("# Title/abstract screening for the Paper 0 light-yield survey.\n")
        fh.write("# auto_* is a crude keyword verdict; reviewed_* is what a human decided\n")
        fh.write("# after reading, and is what the extraction step uses. A blank\n")
        fh.write("# reviewed_decision means not yet reviewed.\n")
        writer = csv.writer(fh)
        writer.writerow(["arxiv_id", "published", "title", "auto_decision", "auto_reason",
                         "reviewed_decision", "reviewed_reason"])
        for paper in papers:
            decision, reason = auto_bucket(paper)
            rev_d, rev_r = reviewed.get(paper["arxiv_id"], ("", ""))
            writer.writerow([paper["arxiv_id"], paper["published"], paper["title"],
                             decision, reason, rev_d, rev_r])

    from collections import Counter
    counts = Counter(auto_bucket(p)[0] for p in papers)
    print(f"screened {len(papers)} candidates -> {SCREENING.relative_to(ROOT)}")
    print(f"  auto-excluded  {counts['exclude']}")
    print(f"  to read        {counts['screen']}")
    print(f"  human-reviewed {len(reviewed)}")


def _to_read() -> list[dict]:
    """Papers whose full text is wanted: auto-kept unless a human excluded them."""
    papers = {p["arxiv_id"]: p for p in json.loads(METADATA.read_text())["papers"]}
    if not SCREENING.exists():
        return [p for p in papers.values() if auto_bucket(p)[0] == "screen"]
    import csv
    out = []
    with SCREENING.open() as fh:
        for row in csv.DictReader(ln for ln in fh if not ln.startswith("#")):
            decision = row["reviewed_decision"] or row["auto_decision"]
            if decision != "exclude" and row["arxiv_id"] in papers:
                out.append(papers[row["arxiv_id"]])
    return out


def download(limit: int | None = None) -> None:
    """Fetch PDFs and extract text. Cached, so re-running is cheap and offline."""
    papers = _to_read()
    print(f"{len(papers)} papers pass screening")
    pdf_dir, txt_dir = CACHE / "pdf", CACHE / "txt"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    txt_dir.mkdir(parents=True, exist_ok=True)

    done = 0
    for paper in papers:
        if limit is not None and done >= limit:
            break
        stem = re.sub(r"[^0-9A-Za-z.\-]", "_", paper["arxiv_id"])
        pdf, txt = pdf_dir / f"{stem}.pdf", txt_dir / f"{stem}.txt"
        if txt.exists():
            continue
        try:
            if not pdf.exists():
                pdf.write_bytes(_get(f"https://arxiv.org/pdf/{paper['arxiv_id']}"))
                time.sleep(POLITE_DELAY_S)
            subprocess.run(["pdftotext", "-layout", str(pdf), str(txt)], check=True)
            print(f"  ok    {paper['arxiv_id']}  {paper['title'][:60]}")
        except Exception as exc:
            print(f"  FAIL  {paper['arxiv_id']}: {exc}")
            continue
        done += 1
    print(f"\n{len(list(txt_dir.glob('*.txt')))} papers with extracted text")


def status() -> None:
    if not METADATA.exists():
        print("no corpus yet -- run with --fetch")
        return
    data = json.loads(METADATA.read_text())
    txt = list((CACHE / "txt").glob("*.txt"))
    print(f"candidates          {len(data['papers'])}")
    print(f"text extracted      {len(txt)}")
    print(f"retrieved           {data['protocol']['retrieved_utc']}")
    records = ROOT / "data" / "literature" / "ly_records.csv"
    if records.exists():
        rows = [ln for ln in records.read_text().splitlines()
                if ln and not ln.startswith("#")]
        print(f"records extracted   {max(0, len(rows) - 1)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--screen", action="store_true")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()

    if args.fetch:
        fetch()
    if args.screen:
        screen()
    if args.download:
        download(args.limit)
    if args.status or not (args.fetch or args.screen or args.download):
        status()


if __name__ == "__main__":
    main()
