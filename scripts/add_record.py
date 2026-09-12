#!/usr/bin/env python3
"""Append one extracted measurement to data/literature/ly_records.csv, validated.

Hand-editing a 24-column CSV is how a survey acquires silent errors: a shifted
column, a typo in a category name, a field quietly left blank when it should say
`not_stated`. This refuses all three.

    python3 scripts/add_record.py --arxiv_id 1607.05486v2 --material "Tl2LaCl5" \\
        --dopant "Ce 10%" --ly_value 51000 --ly_unit photons_per_MeV \\
        --ly_uncertainty 5000 --sample_form single_crystal ...

Every column must be given. `not_stated` is a legitimate and required value: it
is the observation the survey is built to count, and the one thing that must
never be replaced by a plausible guess.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "data" / "literature" / "ly_records.csv"

COLUMNS = [
    "arxiv_id", "doi", "journal_ref", "year", "material", "dopant",
    "ly_value", "ly_unit", "ly_uncertainty", "sample_form", "dimensions_mm",
    "surface_treatment", "wrapping", "coupling", "coupling_index",
    "photodetector", "pde_stated", "reference_standard", "reference_material",
    "reference_value_ph_per_MeV", "source_isotope",
    "source_energy_keV", "shaping_time_us", "temperature_K", "notes",
    "extracted_on",
]

# Controlled vocabularies. A value outside these is a typo or a category the
# schema has not thought about; either way it should stop the script rather than
# land in the file and be counted as something it is not.
VOCAB: dict[str, set[str]] = {
    "ly_unit": {"photons_per_MeV", "photons_per_keV", "photoelectrons_per_MeV",
                "photoelectrons_per_keV", "relative", "not_stated"},
    "sample_form": {"single_crystal", "ceramic", "powder", "thin_film",
                    "composite", "glass", "unclear", "not_stated"},
    "surface_treatment": {"polished", "ground", "saw_cut", "etched", "as_grown",
                          "unclear", "not_stated"},
    "wrapping": {"teflon", "esr", "tyvek", "lumirror", "painted", "aluminium",
                 "millipore", "bas04", "tio2", "none", "multiple", "unclear",
                 "not_stated"},
    "coupling": {"grease", "oil", "air", "glue", "gel", "unclear", "not_stated"},
    "pde_stated": {"yes", "no"},
    # "absolute" means the paper calibrated its own photon counting rather than
    # normalising to another crystal. Anything else names the crystal it was
    # measured against, and reference_value_ph_per_MeV is what that crystal was
    # assumed to be worth -- which is the number this survey exists to compare.
    "reference_material": {"absolute", "NaI:Tl", "CsI:Tl", "BGO", "LYSO:Ce", "LSO:Ce",
                           "YAP:Ce", "other", "not_stated"},
}


# Words that mean the extractor is not sure. A row carrying any of them is
# refused: uncertainty belongs in ly_records.csv.pending, where nothing counts it.
HEDGES = ("provisional", "placeholder", "do not use", "needs re-read",
          "needs re-reading", "before this row is used", "not established")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    for column in COLUMNS:
        if column == "extracted_on":
            continue
        ap.add_argument(f"--{column}", required=True)
    args = vars(ap.parse_args())
    args["extracted_on"] = dt.date.today().isoformat()

    # A value the extractor is unsure of does not belong here at all. Three times
    # in one session a row was entered with the uncertainty written into the
    # notes -- "placeholder for the order of magnitude", "do not use until
    # re-read" -- and each had to be taken out again. A caveat in a notes field
    # is still a number in a data file, and the analysis does not read notes.
    # Put it in ly_records.csv.pending instead.
    problems = []
    lowered = f"{args['notes']} {args['ly_value']}".lower()
    for hedge in HEDGES:
        if hedge in lowered:
            problems.append(
                f"notes contain {hedge!r}: a value you are not sure of belongs in "
                "data/literature/ly_records.csv.pending, not in the records file"
            )
    for field, allowed in VOCAB.items():
        if args[field] not in allowed:
            problems.append(f"{field}={args[field]!r} not one of {sorted(allowed)}")
    for field, value in args.items():
        if value == "":
            problems.append(f"{field} is empty -- use not_stated if the paper does not say")
    if problems:
        raise SystemExit("refused:\n  " + "\n  ".join(problems))

    new = not RECORDS.exists()
    with RECORDS.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        if new:
            writer.writeheader()
        writer.writerow({c: args[c] for c in COLUMNS})
    print(f"recorded {args['material']} from {args['arxiv_id']}")


if __name__ == "__main__":
    main()
