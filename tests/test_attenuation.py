"""Photon-attenuation validation against Geant4.

These tests invoke the compiled ``g4data`` helper, so they are slower than the
pure-Python suite and are skipped automatically when it has not been built.
Run them with ``pytest -m geant4``.

Structure matters here. The suite separates two different kinds of check:

*Gates* compare against **primary** sources -- the PDG tables and the 511 keV
attenuation length, which is the most widely reproduced number in the PET
literature. A failure there means our code is wrong.

*Characterisation* tests compare against a **secondary** compilation
(arXiv:2505.06929 Table 1) and record which of its rows reconcile with Geant4
and which do not. A failure there does not necessarily mean our code is wrong;
it means the relationship to that compilation has changed and needs a human
decision. Two rows are known not to reconcile, and they are pinned explicitly
rather than quietly excluded -- the whole premise of Paper 0 is that published
light-yield and attenuation compilations disagree, so silently dropping the
disagreements would be the one unforgivable move.
"""

from __future__ import annotations

import pytest

from scint.attenuation import Geant4NotBuiltError, attenuation_lengths_cm
from scint.materials import Material, attenuation_depth_cm

pytestmark = pytest.mark.geant4

TOLERANCE_PERCENT = 5.0


def _lengths(material, energies):
    try:
        return attenuation_lengths_cm(material, energies)
    except Geant4NotBuiltError as exc:
        pytest.skip(str(exc))


# --------------------------------------------------------------------------- #
# Gates: primary sources
# --------------------------------------------------------------------------- #

# Attenuation length at 511 keV: the sharpest available cross-check on the whole
# formula -> mass fraction -> Geant4 chain.
ATTENUATION_511_KEV_CM = {
    "BGO": ("Bi4Ge3O12", 7.13, 1.05),
    "LYSO": ("Lu1.8Y0.2SiO5", 7.20, 1.20),
    "NaI": ("NaI", 3.67, 2.93),
}


@pytest.mark.parametrize("name", sorted(ATTENUATION_511_KEV_CM))
def test_attenuation_length_at_511_kev(name):
    formula, density, expected = ATTENUATION_511_KEV_CM[name]
    material = Material(name, formula, density)
    got = _lengths(material, [0.511])[0.511]
    deviation = abs(got - expected) / expected * 100
    assert deviation < TOLERANCE_PERCENT, (
        f"{name}: computed {got:.3f} cm, expected ~{expected} cm ({deviation:.1f} % off)"
    )


def test_denser_lutetium_compound_attenuates_more_than_luag():
    """LYSO carries more lutetium by mass than LuAG and is 7 % denser.

    Its attenuation length must therefore be shorter. The benchmark compilation
    lists identical 88 % depths for the two, which is what first exposed the
    problem with its LYSO row.
    """
    lyso = Material("LYSO", "Lu1.8Y0.2SiO5", 7.20)
    luag = Material("LuAG", "Lu3Al5O12", 6.73)
    assert lyso.mass_fractions["Lu"] > luag.mass_fractions["Lu"]
    assert _lengths(lyso, [0.2])[0.2] < _lengths(luag, [0.2])[0.2]


# --------------------------------------------------------------------------- #
# Characterisation: secondary compilation
# --------------------------------------------------------------------------- #

# 88 % attenuation depths from arXiv:2505.06929 Table 1 (src_A in
# data/benchmark_scintillators.csv), as (formula, density, d88@100keV, d88@200keV).
BENCHMARK_88_PERCENT_DEPTHS = {
    "Pr:LuAG": ("Lu3Al5O12", 6.73, 0.12, 0.65),
    "Ce:GAGG": ("Gd3Al2Ga3O12", 6.63, 0.18, 0.91),
    "LaBr3:Ce": ("LaBr3", 5.08, 0.33, 1.54),
    "CeBr3": ("CeBr3", 5.18, 0.31, 1.48),
    "CsI:Tl": ("CsI", 4.51, 0.23, 1.25),
    "NaI:Tl": ("NaI", 3.67, 0.35, 1.76),
    "BGO": ("Bi4Ge3O12", 7.13, 0.08, 0.43),
}

# Rows that do not reconcile with Geant4, with the measured discrepancy at the
# time of writing. Pinned so that a future data change surfaces as a test failure
# demanding a human decision, rather than passing unnoticed.
KNOWN_IRRECONCILABLE = {
    ("LYSO", 0.2): dict(formula="Lu1.8Y0.2SiO5", density=7.20, reference=0.67, min_deviation=10.0),
    ("BGO", 0.2): dict(formula="Bi4Ge3O12", density=7.13, reference=0.43, min_deviation=5.0),
}


@pytest.mark.parametrize("name", sorted(BENCHMARK_88_PERCENT_DEPTHS))
@pytest.mark.parametrize("energy", [0.1, 0.2])
def test_benchmark_rows_reconcile_with_geant4(name, energy):
    if (name, energy) in KNOWN_IRRECONCILABLE:
        pytest.skip(f"{name} at {energy * 1000:.0f} keV is a pinned discrepancy")
    formula, density, ref100, ref200 = BENCHMARK_88_PERCENT_DEPTHS[name]
    reference = ref100 if energy == 0.1 else ref200
    computed = attenuation_depth_cm(_lengths(Material(name, formula, density), [energy])[energy], 0.88)
    # The table quotes two decimals, so at 0.08 cm a single digit is already 6 %.
    # Compare against the larger of the percentage band and the quoting precision.
    tolerance = max(TOLERANCE_PERCENT / 100 * reference, 0.005)
    assert abs(computed - reference) <= tolerance, (
        f"{name} at {energy * 1000:.0f} keV: computed {computed:.4f} cm, "
        f"reference {reference} cm, tolerance {tolerance:.4f} cm"
    )


@pytest.mark.parametrize("key", sorted(KNOWN_IRRECONCILABLE))
def test_pinned_discrepancies_persist(key):
    """Regression guard on a *finding*, not on our code.

    If one of these starts agreeing, the note recorded against that row in
    data/benchmark_scintillators.csv must be revisited.
    """
    name, energy = key
    spec = KNOWN_IRRECONCILABLE[key]
    material = Material(name, spec["formula"], spec["density"])
    computed = attenuation_depth_cm(_lengths(material, [energy])[energy], 0.88)
    deviation = abs(computed - spec["reference"]) / spec["reference"] * 100
    assert deviation > spec["min_deviation"], (
        f"{name} at {energy * 1000:.0f} keV now agrees to {deviation:.1f} % "
        f"(was > {spec['min_deviation']} %). Revisit the note in the benchmark CSV."
    )
