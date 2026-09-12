"""The survey's data file must stay machine-readable and self-consistent.

A survey is only worth the integrity of its spreadsheet. These tests exist
because one hand-made correction to a citation introduced an unquoted comma,
which shifted every column after it in two rows -- silently, and in a way that
made the affected records report a year as a material name. Nothing in the
analysis would have complained.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "data" / "literature" / "ly_records.csv"
SCREENING = ROOT / "data" / "literature" / "screening.csv"


def _module(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def rows() -> list[dict]:
    with RECORDS.open() as fh:
        return list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))


def test_records_file_exists_and_parses():
    assert RECORDS.exists()
    assert rows(), "no records extracted yet"


def test_every_row_has_exactly_the_declared_columns():
    """A shifted column is the failure this whole file is here to catch."""
    add_record = _module("add_record")
    for i, row in enumerate(rows(), start=1):
        assert None not in row, f"row {i} has more fields than the header"
        missing = [c for c in add_record.COLUMNS if c not in row or row[c] is None]
        assert not missing, f"row {i} is missing {missing}"


def test_no_field_is_blank():
    """Blank is ambiguous: it could mean 'not stated' or 'not yet looked at'.

    The schema requires the distinction to be explicit, so blank is an error.
    """
    add_record = _module("add_record")
    for i, row in enumerate(rows(), start=1):
        for column in add_record.COLUMNS:
            assert row[column].strip() != "", f"row {i}: {column} is blank"


def test_controlled_vocabularies_hold():
    add_record = _module("add_record")
    for i, row in enumerate(rows(), start=1):
        for field, allowed in add_record.VOCAB.items():
            assert row[field] in allowed, (
                f"row {i}: {field}={row[field]!r} is outside the vocabulary"
            )


def test_light_yield_values_are_numeric_and_positive():
    for i, row in enumerate(rows(), start=1):
        value = float(row["ly_value"])
        assert value > 0, f"row {i}: non-positive yield"


def test_material_column_never_holds_a_year():
    """The exact symptom of the column shift that prompted these tests."""
    for i, row in enumerate(rows(), start=1):
        material = row["material"].strip()
        assert not (material.isdigit() and 1900 < int(material) < 2100), (
            f"row {i}: material is {material!r}, which is a year -- columns have shifted"
        )


def test_screening_sheet_is_consistent_with_the_corpus():
    if not SCREENING.exists():
        pytest.skip("screening sheet not built yet")
    with SCREENING.open() as fh:
        screened = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    assert screened
    ids = [r["arxiv_id"] for r in screened]
    assert len(ids) == len(set(ids)), "the same paper is screened twice"
    for row in screened:
        assert row["auto_decision"] in {"screen", "exclude"}
        assert row["auto_reason"], f"{row['arxiv_id']} screened without a reason"


def test_extracted_papers_were_all_in_the_screened_corpus():
    """A record must come from the declared sample, or the sample is not the sample."""
    if not SCREENING.exists():
        pytest.skip("screening sheet not built yet")
    with SCREENING.open() as fh:
        corpus = {r["arxiv_id"].split("v")[0]
                  for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))}
    for row in rows():
        base = row["arxiv_id"].split("v")[0]
        assert base in corpus, (
            f"{row['arxiv_id']} was extracted but is not in the screened corpus"
        )


def test_no_record_hedges_its_own_value():
    """A caveat in a notes field is still a number in a data file.

    Three rows were entered in one session with the doubt written into the notes
    -- "placeholder for the order of magnitude", "do not use until re-read" --
    and each had to be removed again. The analysis does not read notes, so a
    hedged row counts exactly like a confident one. Uncertain values go to
    ly_records.csv.pending; add_record.py now refuses them, and this checks the
    file directly in case a row arrives by some other route.
    """
    add_record = _module("add_record")
    for i, row in enumerate(rows(), start=1):
        lowered = row["notes"].lower()
        for hedge in add_record.HEDGES if hasattr(add_record, "HEDGES") else (
            "provisional", "placeholder", "do not use", "needs re-read",
        ):
            assert hedge not in lowered, (
                f"row {i} ({row['arxiv_id']}, {row['material']}) hedges its own "
                f"value with {hedge!r} -- it belongs in the pending file"
            )


def test_pending_file_exists_and_explains_itself():
    """Deferrals must be recorded, or "not extracted" is indistinguishable from
    "not noticed"."""
    pending = ROOT / "data" / "literature" / "ly_records.csv.pending"
    assert pending.exists()
    text = pending.read_text()
    assert "NOT recorded" in text or "EXCLUDED" in text
