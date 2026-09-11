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
    assert len(supported_options()) == 24


def test_every_treatment_works_with_every_air_coupled_wrapping(options):
    air = {(o.treatment, o.wrapping) for o in options if o.coupling == "air" and o.supported}
    assert air == {(t, w) for t in TREATMENTS for w in WRAPPINGS}


def test_glue_coupling_is_available_only_for_specular_films(options):
    """Geant4 ships glue tables for lumirror and ESR only; the diffuse
    reflectors have none, and the catalogue must not invent them."""
    glued = {o.wrapping for o in options if o.coupling == "glue" and o.supported}
    assert glued == {"lumirror", "esr"}


def test_unsupported_combinations_carry_an_explanation(options):
    unsupported = [o for o in options if not o.supported]
    assert len(unsupported) == 12
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
