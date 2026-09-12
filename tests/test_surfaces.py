"""The surface grid that defines the light-collection study's parameter space."""

from __future__ import annotations

import pytest

from scint.geant4 import Geant4NotBuiltError
from scint.surfaces import COUPLINGS, TREATMENTS, WRAPPINGS, catalogue, supported_options

pytestmark = pytest.mark.geant4


@pytest.fixture(scope="module")
def options():
    try:
        return catalogue()
    except Geant4NotBuiltError as exc:
        pytest.skip(str(exc))


def test_catalogue_covers_the_full_cross_product(options):
    assert len(options) == len(TREATMENTS) * len(WRAPPINGS) * len(COUPLINGS) == 36


def test_geant4_supports_twenty_four_combinations(options):
    """3 treatments x 6 air-coupled wrappings, plus the 3 x 2 glue-coupled cases."""
    # 21, not 24. The three bare finishes (polishedair, etchedair, groundair)
    # exist as Geant4 enumerators but have no look-up table behind them:
    # G4OpticalSurface::ReadLUTFile has no case for them, so nothing is loaded
    # and G4OpBoundaryProcess::DielectricLUT() then loops forever on an all-zero
    # angular distribution. This catalogue counted them as supported until that
    # hang was traced; the number is pinned here so the mistake cannot return.
    assert len(supported_options()) == 21


def test_every_treatment_works_with_every_air_coupled_reflector(options):
    """Every treatment pairs with every air-coupled REFLECTOR -- but not with
    'none', which has no look-up table and is handled analytically instead."""
    air = {(o.treatment, o.wrapping) for o in options if o.coupling == "air" and o.supported}
    reflectors = {w for w in WRAPPINGS if w != "none"}
    assert air == {(t, w) for t in TREATMENTS for w in reflectors}


def test_glue_coupling_is_available_only_for_specular_films(options):
    """Geant4 ships glue tables for lumirror and ESR only; the diffuse
    reflectors have none, and the catalogue must not invent them."""
    glued = {o.wrapping for o in options if o.coupling == "glue" and o.supported}
    assert glued == {"lumirror", "esr"}


def test_unsupported_combinations_carry_an_explanation(options):
    unsupported = [o for o in options if not o.supported]
    # 15: the 12 glue combinations Geant4 has no table for, plus the 3 bare
    # ones whose absence is what made the simulation hang.
    assert len(unsupported) == 15
    for option in unsupported:
        assert option.reason, f"{option.label} was refused without saying why"
        assert not option.g4_finish, "a refused combination must not name a finish"


def test_supported_combinations_name_a_distinct_geant4_finish(options):
    finishes = [o.g4_finish for o in options if o.supported]
    assert all(finishes)
    assert len(set(finishes)) == len(finishes), "two combinations map to the same finish"


def test_esr_and_vm2000_are_the_same_material(options):
    """Geant4 calls the 3M film vm2000; the rest of the field calls it ESR."""
    esr = {o.g4_finish for o in options if o.wrapping == "esr" and o.supported}
    assert all("vm2000" in finish for finish in esr)


def test_bare_surfaces_are_refused_with_the_real_reason():
    """A bare surface must be refused, and refused for the right reason.

    Not because of glue coupling, which is the generic message and was what this
    case used to return. The actual cause is that the look-up table does not
    exist, and a caller who reads the reason should be sent to the analytic model
    rather than to a different wrapping.
    """
    from scint.surfaces import resolve

    for treatment in ("polished", "etched", "ground"):
        result = resolve(treatment, "none", "air")
        assert not result.supported
        assert "no look-up table for a bare surface" in result.reason
        assert "UNIFIED" in result.reason


def test_no_bare_combination_is_listed_as_supported():
    """The failure mode this catalogue exists to prevent, checked directly."""
    assert not any(option.wrapping == "none" for option in supported_options())
