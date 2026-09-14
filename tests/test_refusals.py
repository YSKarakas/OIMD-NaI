"""The three misconfigurations the application must refuse rather than run.

Each was once accepted silently: a wrapped crystal with no reflectance ran as a
perfect mirror; the 'davis' model handed an LBNL finish to the DAVIS reader and
hung; a bare 'ground' surface with sigma_alpha = 0 was a polished one under
another name. The refusal messages are the ones a user will read, so the tests
pin them. Needs a built binary; the host build is fine, the guards are ours.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build" / "sim" / "scint_optical"

BASE = """/random/setSeeds 1 1
/scint/crystal/shape cylinder
/scint/crystal/diameter 76.2 mm
/scint/crystal/length 76.2 mm
/scint/crystal/materialSpec materials/NaI_Tl.dat
/scint/surface/treatment {treatment}
/scint/surface/wrapping {wrapping}
/scint/surface/coupling air
/scint/surface/model {model}
{extra}
/scint/coupling/type grease
/scint/coupling/rindex 1.465
/scint/coupling/thickness 0.1 mm
/scint/readout/efficiency 1.0
/scint/source/energy 662 keV
/run/initialize
/run/beamOn 1
"""


def _run(tmp_path, **kw):
    macro = tmp_path / "m.mac"
    macro.write_text(BASE.format(**kw))
    proc = subprocess.run([str(BINARY), str(macro)], cwd=ROOT, capture_output=True,
                          text=True, timeout=120)
    return proc.returncode, proc.stdout + proc.stderr


@pytest.fixture(autouse=True)
def _need_binary():
    if not BINARY.exists():
        pytest.skip(f"{BINARY} not built")


def test_a_wrapped_crystal_without_a_reflectance_is_refused(tmp_path):
    rc, out = _run(tmp_path, treatment="polished", wrapping="teflon", model="lut", extra="")
    assert rc != 0
    assert "silent default of 1.0 is a lossless mirror" in out


def test_the_davis_model_is_refused_rather_than_hung(tmp_path):
    rc, out = _run(tmp_path, treatment="polished", wrapping="teflon", model="davis",
                   extra="/scint/surface/reflectivity 0.99")
    assert rc != 0
    assert "ReadLUTDAVISFile" in out


def test_a_bare_ground_surface_needs_a_roughness(tmp_path):
    rc, out = _run(tmp_path, treatment="ground", wrapping="none", model="lut", extra="")
    assert rc != 0
    assert "sigmaAlpha" in out
