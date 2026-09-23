#!/usr/bin/env python3
"""List the retrieved full texts that carry no recorded disposition.

Of the 132 full texts the survey retrieved, each should either have been
extracted (data/literature/ly_records.csv) or read and deferred with a reason
(data/literature/ly_records.csv.pending). The remainder carry no disposition at
all. Until this file existed they could only be found by difference; the
manuscript says the archive lists them, so it now does, with the keyword class
scripts/scan_yields.py gave each:

  absolute   a yield in photons/MeV
  detected   a yield in photoelectrons/MeV
  relative   a yield relative to another crystal
  none       no light yield in numerical form

No reason was recorded for leaving any of them unread; the output says so
rather than supplying one.

    python3 scripts/list_undisposed.py        # writes data/literature/undisposed.csv
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
OUT = LIT / "undisposed.csv"


def _rows(path: Path) -> list[dict]:
    return list(csv.DictReader(ln for ln in path.read_text().splitlines() if not ln.startswith("#")))


def main() -> None:
    texts = _rows(LIT / "fulltext_classes.csv")
    extracted = {r["arxiv_id"] for r in _rows(LIT / "ly_records.csv")}
    # The pending file is prose: one entry per identifier at the start of a line.
    pending = {m.group(1) for m in re.finditer(r"^(\d{4}\.\d{4,5}v\d+)\b",
                                                (LIT / "ly_records.csv.pending").read_text(), re.M)}
    left = [t for t in texts if t["arxiv_id"] not in extracted and t["arxiv_id"] not in pending]
    left.sort(key=lambda t: (t["klass"] == "none", t["klass"], t["arxiv_id"]))
    with OUT.open("w", newline="") as fh:
        fh.write(
            "# Retrieved full texts with no recorded disposition: neither extracted\n"
            "# (ly_records.csv) nor read and deferred with a reason\n"
            "# (ly_records.csv.pending). Written by scripts/list_undisposed.py.\n"
            "# klass is the keyword classification of scripts/scan_yields.py. No\n"
            "# reason was recorded for leaving any of these unread.\n"
        )
        w = csv.writer(fh)
        w.writerow(["arxiv_id", "klass", "year", "title"])
        for t in left:
            w.writerow([t["arxiv_id"], t["klass"], t["year"], t["title"]])
    counts = Counter(t["klass"] for t in left)
    print(f"{len(texts)} full texts, {len(extracted)} extracted papers, {len(pending)} deferred")
    print(f"{len(left)} without a disposition: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
