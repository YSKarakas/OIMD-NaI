"""Layer 0 validation.

Reference values come from the PDG *Atomic and Nuclear Properties* tables
(pdg.lbl.gov/2024/AtomicNuclearProperties, retrieved 2026-09-11) rather than from
a single secondary compilation. That choice is deliberate: during development the
LYSO row of the benchmark compilation in ``data/benchmark_scintillators.csv``
turned out to disagree with a Geant4 computation by 17-20 %, while every other
row agreed to within 2 %. Validating against a primary tabulation prevents an
error in a secondary source from being baked into the pipeline.
"""

from __future__ import annotations

import pytest

from scint.materials import (
    FormulaError,
    Material,
    attenuation_depth_cm,
    parse_formula,
)

pytest_plugins: list[str] = []

# PDG 2024 Atomic and Nuclear Properties tables.
PDG_REFERENCE = {
    "BGO": {
        "formula": "Bi4Ge3O12",
        "density_g_cm3": 7.130,
        "radiation_length_g_cm2": 7.97,
        "radiation_length_cm": 1.118,
        "moliere_radius_g_cm2": 16.10,
        "moliere_radius_cm": 2.259,
        "nuclear_interaction_length_cm": 22.32,
    },
    "NaI": {
        "formula": "NaI",
        "density_g_cm3": 3.667,
        "radiation_length_g_cm2": 9.49,
        "radiation_length_cm": 2.588,
        "moliere_radius_g_cm2": 15.05,
        "moliere_radius_cm": 4.105,
        "nuclear_interaction_length_cm": 42.16,
    },
}

PDG_TOLERANCE_PERCENT = 5.0


# --------------------------------------------------------------------------- #
# Formula parsing
# --------------------------------------------------------------------------- #

def test_parses_simple_formula():
    assert parse_formula("NaI") == {"Na": 1.0, "I": 1.0}


def test_parses_garnet_with_multiple_subscripts():
    assert parse_formula("Gd3Al2Ga3O12") == {"Gd": 3.0, "Al": 2.0, "Ga": 3.0, "O": 12.0}


def test_parses_nested_groups_and_fractional_subscripts():
    """Multicomponent garnets must write naturally, including partial occupancy."""
    assert parse_formula("(Lu0.5Gd2.5)(Al2Ga3)O12") == {
        "Lu": 0.5, "Gd": 2.5, "Al": 2.0, "Ga": 3.0, "O": 12.0,
    }


def test_group_multiplier_applies_to_every_element_inside():
    assert parse_formula("(AlO2)3") == {"Al": 3.0, "O": 6.0}


def test_repeated_element_accumulates():
    assert parse_formula("CH3CH3") == {"C": 2.0, "H": 6.0}


@pytest.mark.parametrize("bad", ["", "(NaI", "NaI)", "Xx2", "3Na"])
def test_rejects_malformed_formulae(bad):
    with pytest.raises(FormulaError):
        parse_formula(bad)


# --------------------------------------------------------------------------- #
# Stoichiometry
# --------------------------------------------------------------------------- #

def test_mass_fractions_sum_to_one():
    m = Material("GAGG", "Gd3Al2Ga3O12", 6.63)
    assert sum(m.mass_fractions.values()) == pytest.approx(1.0, abs=1e-12)


def test_electron_fractions_sum_to_one():
    m = Material("GAGG", "Gd3Al2Ga3O12", 6.63)
    assert sum(m.electron_fractions.values()) == pytest.approx(1.0, abs=1e-12)


def test_massfrac_arg_round_trips_to_unity():
    """g4data rejects mass fractions that miss unity by more than 1e-6."""
    for formula, rho in [("Gd3Al2Ga3O12", 6.63), ("Lu1.8Y0.2SiO5", 7.20), ("Bi4Ge3O12", 7.13)]:
        m = Material("x", formula, rho)
        total = sum(float(p.split(":")[1]) for p in m.massfrac_arg().split(","))
        assert total == pytest.approx(1.0, abs=1e-9)


def test_rejects_non_positive_density():
    with pytest.raises(ValueError):
        Material("bad", "NaI", 0.0)


def test_zeff_exponent_is_reported_with_the_value():
    """Z_eff has no single accepted convention, so the exponent must travel with it."""
    m = Material("GAGG", "Gd3Al2Ga3O12", 6.63, zeff_exponent=2.94)
    assert m.summary()["z_eff_exponent"] == 2.94
    assert m.z_eff != Material("GAGG", "Gd3Al2Ga3O12", 6.63, zeff_exponent=3.5).z_eff


def test_zeff_lies_between_min_and_max_constituent_z():
    m = Material("GAGG", "Gd3Al2Ga3O12", 6.63)
    assert 8 < m.z_eff < 64  # O and Gd bound it


# --------------------------------------------------------------------------- #
# PDG validation gate
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", sorted(PDG_REFERENCE))
@pytest.mark.parametrize(
    "quantity",
    [
        "radiation_length_g_cm2",
        "radiation_length_cm",
        "moliere_radius_g_cm2",
        "moliere_radius_cm",
        "nuclear_interaction_length_cm",
    ],
)
def test_matches_pdg_tables(name, quantity):
    ref = PDG_REFERENCE[name]
    material = Material(name, ref["formula"], ref["density_g_cm3"])
    computed = getattr(material, quantity)
    deviation = abs(computed - ref[quantity]) / ref[quantity] * 100
    assert deviation < PDG_TOLERANCE_PERCENT, (
        f"{name}.{quantity}: computed {computed:.4f}, PDG {ref[quantity]:.4f} "
        f"({deviation:.1f} % off, tolerance {PDG_TOLERANCE_PERCENT} %)"
    )


# --------------------------------------------------------------------------- #
# Attenuation depth helper
# --------------------------------------------------------------------------- #

def test_attenuation_depth_inverts_the_exponential():
    from math import exp
    depth = attenuation_depth_cm(2.0, 0.88)
    assert 1.0 - exp(-depth / 2.0) == pytest.approx(0.88)


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.1, 1.5])
def test_attenuation_depth_rejects_out_of_range_fractions(bad):
    with pytest.raises(ValueError):
        attenuation_depth_cm(1.0, bad)
