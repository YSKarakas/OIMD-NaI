#!/usr/bin/env python3
"""How do published optical simulations obtain their bulk attenuation input?

Why this exists
---------------
The paper claims that feeding a *measured*, wavelength-dependent attenuation
curve to a Geant4 optical simulation of NaI(Tl) is not current practice. A
claim of that shape is easy to make and hard to trust, so it is measured here
rather than asserted, against the same 132 full texts the light-yield survey
reads, and the classification is mechanical and rebuildable.

What it measures, and what it does not
--------------------------------------
This is a keyword classification of full text, so it is an UPPER BOUND on good
practice in both directions: a paper counts as using a measured curve if the
words appear near each other, whether or not it really did, and a paper that
did so without saying it in words this scan knows is missed. It reads the
arXiv-accessible corpus of the light-yield survey, which is a sample of the
literature and not a census of optical simulations. What it supports is a
statement about a defined sample, which is the only kind of statement the
evidence allows.

    python3 scripts/scan_optical_practice.py
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TXT = ROOT / "data" / "literature" / "cache" / "txt"
OUT = ROOT / "data" / "literature" / "optical_practice.csv"

ATTEN = re.compile(r"(absorption|attenuation)\s+lengths?", re.I)
SIMULATION = re.compile(r"geant\s*-?4|\bGATE\b", re.I)
# "measured" in the sense that matters: the paper says the curve came from an
# instrument, or from published data it processed, not from a datasheet number.
MEASURED = re.compile(
    r"spectrophotometer|spectro-photometer|transmittance measurement|"
    r"transmission measurement|measured (?:the )?(?:optical )?"
    r"(?:transmission|transmittance)|digitis|digitiz", re.I)
# a curve rather than a single number
DISPERSIVE = re.compile(
    r"wavelength[- ]depend\w*|spectral depend\w*|as a function of wavelength|"
    r"wavelength dependence of the (?:absorption|attenuation)", re.I)
NAI = re.compile(r"NaI\s*[:(]?\s*Tl|sodium\s+iodide", re.I)

WINDOW = 400   # characters either side of an attenuation-length mention


def classify(text: str) -> dict:
    flat = re.sub(r"\s+", " ", text)
    hits = list(ATTEN.finditer(flat))
    near = " ".join(flat[max(0, m.start() - WINDOW):m.end() + WINDOW]
                    for m in hits)
    return {
        "mentions_attenuation": bool(hits),
        "runs_optical_sim": bool(SIMULATION.search(flat)),
        "dispersive_near_mention": bool(DISPERSIVE.search(near)),
        "measured_near_mention": bool(MEASURED.search(near)),
        "mentions_nai": bool(NAI.search(flat)),
    }


def main() -> int:
    if not TXT.is_dir():
        raise SystemExit(f"{TXT} not found -- run the survey fetch step first")
    rows = []
    for f in sorted(TXT.glob("*.txt")):
        c = classify(f.read_text(errors="ignore"))
        c["arxiv_id"] = f.stem
        rows.append(c)
    if not rows:
        raise SystemExit(f"no texts in {TXT}")

    cols = ["arxiv_id", "mentions_attenuation", "runs_optical_sim",
            "dispersive_near_mention", "measured_near_mention", "mentions_nai"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        fh.write(
            "# How published optical simulations obtain their bulk attenuation.\n"
            "# Produced by scripts/scan_optical_practice.py -- mechanical and\n"
            "# rebuildable. A KEYWORD classification of full text, so every\n"
            "# column is an upper bound on the practice it names: the words\n"
            "# appearing near each other is not proof the paper did the thing.\n"
            "# The corpus is the arXiv-accessible light-yield survey sample,\n"
            "# not a census of optical simulations.\n"
            "#\n")
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})

    sim = [r for r in rows if r["runs_optical_sim"] and r["mentions_attenuation"]]
    disp = [r for r in sim if r["dispersive_near_mention"]]
    meas = [r for r in disp if r["measured_near_mention"]]
    nai = [r for r in meas if r["mentions_nai"]]
    print(f"corpus                                            {len(rows):4d}")
    print(f"  run an optical simulation AND name an           {len(sim):4d}")
    print(f"  attenuation length")
    print(f"    ... as a wavelength-dependent quantity        {len(disp):4d}")
    print(f"        ... obtained from a measurement           {len(meas):4d}")
    for r in meas:
        print(f"            {r['arxiv_id']}")
    print(f"            ... for NaI(Tl)                       {len(nai):4d}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
