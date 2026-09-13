#!/usr/bin/env python3
"""Turn the extracted records into the survey's two results.

The first is reporting completeness: for each field needed to interpret a
published light yield, how often the paper actually gives it. The second is the
reference-value problem: absolute yields are usually measured against another
crystal, and the assumed value of that crystal varies between papers, so the
ruler itself is uncertain.

    python3 scripts/analyse_literature.py
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "data" / "literature" / "ly_records.csv"
CLASSES = ROOT / "data" / "literature" / "fulltext_classes.csv"

# The fields a reader needs in order to know what a reported number means, in
# the order they matter for correcting it back to an intrinsic yield.
REQUIRED = [
    ("dimensions_mm", "sample dimensions"),
    ("surface_treatment", "surface treatment"),
    ("wrapping", "wrapping / reflector"),
    ("coupling", "optical coupling"),
    ("coupling_index", "coupling refractive index"),
    ("photodetector", "photodetector model"),
    ("pde_stated", "photodetector QE or PDE"),
    ("reference_standard", "reference standard used"),
    ("shaping_time_us", "shaping / integration time"),
    ("temperature_K", "temperature"),
    ("ly_uncertainty", "uncertainty on the yield"),
]

MISSING = {"not_stated", "unclear", "", "no"}


def load() -> list[dict]:
    if not RECORDS.exists():
        return []
    with RECORDS.open() as fh:
        return list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))


def main() -> None:
    rows = load()
    if not rows:
        raise SystemExit("no records yet")

    papers = {r["arxiv_id"] for r in rows}
    print(f"{len(rows)} measurements extracted from {len(papers)} papers\n")

    print("=" * 74)
    print("REPORTING COMPLETENESS -- does the paper say what its number means?")
    print("=" * 74)
    print(f"  {'field':<32} {'stated':>8} {'of':>4} {'':>4} {'percent':>8}")
    for field, label in REQUIRED:
        stated = sum(1 for r in rows if r[field].strip().lower() not in MISSING)
        pct = 100 * stated / len(rows)
        bar = "#" * int(pct / 5)
        print(f"  {label:<32} {stated:>8} {len(rows):>4} {'':>4} {pct:>7.1f} % {bar}")

    complete = sum(
        1 for r in rows
        if all(r[f].strip().lower() not in MISSING for f, _ in REQUIRED)
    )
    print(f"\n  records stating ALL of the above: {complete} of {len(rows)} "
          f"({100 * complete / len(rows):.1f} %)")

    print()
    print("=" * 74)
    print("THE RULER PROBLEM -- what absolute yields are measured against")
    print("=" * 74)
    absolute = [r for r in rows if r["reference_standard"].lower().startswith("absolute")]
    relative = [r for r in rows if not r["reference_standard"].lower().startswith("absolute")
                and r["reference_standard"].strip().lower() not in MISSING]
    print(f"  measured absolutely (own QE calibration)  {len(absolute)}")
    print(f"  measured against another scintillator     {len(relative)}")
    for r in relative:
        print(f"    {r['material']:<14} {r['arxiv_id']:<14} -> {r['reference_standard']}")

    # The sharpest form of the problem, read from the structured columns rather
    # than parsed out of prose: papers that calibrate against the same reference
    # crystal do not agree on what that crystal is worth, so two yields quoted in
    # the same unit were measured with different rulers.
    #
    # PUBLISHED is what the compilations and other papers in this project's own
    # evidence base say for the same material. It is kept next to the assumed
    # values because the comparison is the point.
    PUBLISHED = {
        "BGO": ["8000 (arXiv:1607.05486 Table I)", "8200 (arXiv:1308.3908, measured)",
                "8500 (Miller et al., IEEE TNS 72 197 (2025))", "10000 (Bonesini arXiv:2505.06929)"],
        "NaI:Tl": ["38000 (Bonesini; also arXiv:1308.3908, measured)",
                   "40000 (ANAIS arXiv:1703.01262, cited)",
                   "41000 (Miller et al., IEEE TNS 72 197 (2025))"],
        "LYSO:Ce": ["26000 (arXiv:1309.3736, simulation input)",
                    "32000 (arXiv:1607.05486 Table I, same paper that assumes 33000)"],
        "GSO:Ce": ["11000 (arXiv:2104.07568, assumed)",
                   "12500 (arXiv:1607.05486 Table I)"],
    }
    assumed: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        material = r.get("reference_material", "")
        value = r.get("reference_value_ph_per_MeV", "")
        # "other" is a bucket of DIFFERENT reference materials, so a spread
        # computed across it would compare a plastic scintillator with a
        # silicate and report a meaningless number. Those rows are listed
        # individually below instead.
        if material and material not in {"absolute", "not_stated", "other"} and value.isdigit():
            assumed[material].add(value)

    print()
    print("=" * 74)
    print("THE RULER ITSELF -- what the reference crystals were assumed to be worth")
    print("=" * 74)
    if not assumed:
        print("  no relative measurements recorded yet")
    for material in sorted(set(assumed) | set(PUBLISHED)):
        values = sorted(assumed.get(material, set()), key=float)
        if values:
            lo, hi = float(values[0]), float(values[-1])
            spread = f"   ({100 * (hi - lo) / lo:.0f} % spread)" if len(values) > 1 else ""
            print(f"  {material}  assumed here as: {', '.join(values)} ph/MeV{spread}")
        else:
            print(f"  {material}")
        for other in PUBLISHED.get(material, []):
            print(f"      elsewhere published as {other}")

    others = [(r["material"], r["reference_standard"]) for r in rows
              if r.get("reference_material") == "other"]
    if others:
        print("\n  references outside the common set (listed, not pooled --")
        print("  a spread across different materials would mean nothing):")
        for material, standard in sorted(set(others)):
            print(f"    {material:<26} -> {standard}")

    print()
    print("=" * 74)
    print("SAME MATERIAL, DIFFERENT ANSWERS")
    print("=" * 74)
    by_material: dict[str, list[tuple[float, str, str]]] = defaultdict(list)
    for r in rows:
        if r["ly_unit"] != "photons_per_MeV":
            continue
        try:
            value = float(r["ly_value"])
        except ValueError:
            continue
        # Group by material AND dopant: an undoped and a doped crystal of the
        # same host are different scintillators, and comparing them would
        # manufacture a disagreement that is not there.
        key = f"{r['material'].strip()} : {r['dopant'].strip()}"
        by_material[key].append((value, r["arxiv_id"], r["dopant"]))
    for material, values in sorted(by_material.items()):
        # Only a disagreement BETWEEN papers is evidence about the literature.
        if len({paper for _, paper, _ in values}) < 2:
            continue
        lo, hi = min(v for v, _, _ in values), max(v for v, _, _ in values)
        print(f"  {material}: {lo:.0f} - {hi:.0f} ph/MeV  ({100 * (hi - lo) / lo:.0f} % spread)")
        for value, paper, dopant in sorted(values):
            print(f"    {value:>9.0f}  {paper:<14} {dopant}")

    if CLASSES.exists():
        print()
        print("=" * 74)
        print("FULL-TEXT SCREENING -- what form the literature reports yields in")
        print("=" * 74)
        with CLASSES.open() as fh:
            klasses = Counter(row["klass"]
                              for row in csv.DictReader(ln for ln in fh
                                                        if not ln.startswith("#")))
        total = sum(klasses.values())
        for name in ("absolute", "detected", "relative", "none"):
            print(f"  {name:<10} {klasses[name]:>4}  ({100 * klasses[name] / total:.0f} %)")
        print(f"  {'total':<10} {total:>4}")


if __name__ == "__main__":
    main()
