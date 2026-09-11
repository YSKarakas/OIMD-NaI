"""Guards on Geant4 version identification.

These tests encode a lesson learned the hard way: ``geant4-config --version``
reported ``11.4.0`` on a machine that was in fact running
``geant4-11-04-beta-01``, because the beta of a series already carries the
target ``G4VERSION_NUMBER``. Provenance records and cache keys must therefore be
built on the version *tag*, never on the config version string.
"""

from __future__ import annotations

import pytest

from scint.geant4 import Geant4NotBuiltError, geant4_env, is_prerelease, version_info, version_key


@pytest.mark.parametrize(
    "tag",
    [
        "geant4-11-04-beta-01",
        "geant4-11-02-beta-01",
        "geant4-11-03-rc-01",
        "geant4-11-04-ref-05",
        "GEANT4-11-04-BETA-01",
    ],
)
def test_prerelease_tags_are_detected(tag):
    assert is_prerelease(tag) is True


@pytest.mark.parametrize("tag", ["geant4-11-04", "geant4-11-04-patch-02", "geant4-11-03-patch-01"])
def test_release_tags_are_not_flagged(tag):
    assert is_prerelease(tag) is False


def test_unknown_tag_is_unknown_not_false():
    """A missing tag must not be reported as 'this is a release'."""
    assert is_prerelease(None) is None
    assert is_prerelease("") is None


def test_geant4_env_does_not_override_variables_already_set():
    """Inside the container the dataset paths are baked in and must survive."""
    env = geant4_env({"G4LEDATA": "/container/path"})
    assert env["G4LEDATA"] == "/container/path"
    assert env["G4REALSURFACEDATA"].endswith("RealSurface2.2")


def test_geant4_env_fills_in_missing_variables():
    env = geant4_env({})
    assert "G4LEDATA" in env and "G4ENSDFSTATEDATA" in env


@pytest.mark.geant4
def test_version_info_reports_the_library_tag_not_the_config_string():
    try:
        info = version_info()
    except Geant4NotBuiltError:
        pytest.skip("g4data not built")
    assert info["source"] == "g4data"
    assert info["tag"], "the version tag is what distinguishes a beta from a release"
    assert isinstance(info["is_prerelease"], bool)


@pytest.mark.geant4
def test_version_key_is_usable_in_a_cache_key():
    try:
        key = version_key()
    except Geant4NotBuiltError:
        pytest.skip("g4data not built")
    assert key and key != "unknown-geant4"
    assert "/" not in key and " " not in key
