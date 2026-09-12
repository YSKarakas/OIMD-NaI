"""The optical input models are pinned to the published numbers they came from.

Every assertion here is a value that appears in a paper or in a public database,
so a change to scint/optical.py that silently moves a refractive index or an
absorption length away from its source fails the suite rather than quietly
propagating into a result.
"""

from __future__ import annotations

import math

import pytest

from scint.optical import (
    gaussian_emission_sampled_in_wavelength,
    JELLISON_2012_ENDPOINTS,
    LI_1976,
    MAO_2008_ANCHOR_MM,
    MAO_2008_ANCHOR_NM,
    AbsorptionEdge,
    constant_dispersion,
    energy_to_wavelength_nm,
    gaussian_emission,
    wavelength_to_energy_eV,
)


# --- refractive index ------------------------------------------------------ #


def test_li_1976_reproduces_the_sodium_d_line_value():
    """n(589.3 nm) = 1.7745 is the classical handbook value for NaI.

    It is the check that the Sellmeier coefficients transcribed from
    refractiveindex.info really are Li's, and not a corrupted copy.
    """
    assert LI_1976(589.3) == pytest.approx(1.7745, abs=5e-5)


def test_jellison_reproduces_both_published_measurements():
    """The two values quoted verbatim in the abstract of Jellison et al. (2012).

    The coefficients in scint/optical.py are constructed to pass through these,
    so this test is what defines them; if it fails the coefficients have drifted.
    """
    assert JELLISON_2012_ENDPOINTS(436.0) == pytest.approx(1.839, abs=5e-4)
    assert JELLISON_2012_ENDPOINTS(633.0) == pytest.approx(1.786, abs=5e-4)


def test_refractiveindex_info_transcription_of_jellison_is_inconsistent():
    """Pin the discrepancy that made us refit rather than copy.

    refractiveindex.info encodes Jellison as a single-term Sellmeier with
    coefficients (0, 1.994, 0.176). That form reproduces the published 436 nm
    value but falls about 0.008 below the published 633 nm value -- four times the
    quoted uncertainty of 0.002, which a fit with a reported chi-squared of 1.02
    cannot produce. One of the two published items is wrong; this test records
    the inconsistency rather than choosing a winner, and fails if a future
    version of the database changes its numbers.
    """
    def transcribed(lam_nm: float) -> float:
        lam2 = (lam_nm / 1000.0) ** 2
        return math.sqrt(1 + 1.994 * lam2 / (lam2 - 0.176**2))

    assert transcribed(436.0) == pytest.approx(1.839, abs=5e-4)
    shortfall = 1.786 - transcribed(633.0)
    assert shortfall == pytest.approx(0.008, abs=1e-3)
    assert shortfall > 4 * 0.002 - 1e-3


def test_the_two_dispersion_models_disagree_at_the_emission_peak():
    """The finding that motivates carrying both: 1.6 % apart at 415 nm.

    That is not a rounding difference. It propagates into the critical angle and
    therefore into every light-collection number computed from it.
    """
    li, jellison = LI_1976(415.0), JELLISON_2012_ENDPOINTS(415.0)
    assert li == pytest.approx(1.8218, abs=1e-3)
    assert jellison == pytest.approx(1.8504, abs=1e-3)
    assert (jellison - li) / li == pytest.approx(0.0157, abs=2e-3)


def test_the_emission_band_reaches_below_the_validated_range_of_both_models():
    """Neither published measurement covers the blue half of the emission band.

    Jellison's measurement starts at 436 nm -- above the 415 nm emission peak --
    so the peak itself is extrapolated. This is why the material files carry an
    explicit extrapolation warning.
    """
    assert JELLISON_2012_ENDPOINTS.extrapolating(415.0)
    assert JELLISON_2012_ENDPOINTS.extrapolating(350.0)
    assert not JELLISON_2012_ENDPOINTS.extrapolating(500.0)
    # Li's range starts at 250 nm, so it spans the band without extrapolating.
    assert not LI_1976.extrapolating(340.0)


def test_constant_dispersion_is_flat_and_never_extrapolates():
    flat = constant_dispersion(1.85, source="test")
    assert flat(300.0) == flat(700.0) == 1.85
    assert not flat.extrapolating(100.0)


# --- absorption ------------------------------------------------------------ #


def test_absorption_edge_honours_the_measured_anchor():
    """Whatever the slope, the model must pass through the one measured point."""
    for e_u in (0.10, 0.175, 0.25):
        edge = AbsorptionEdge(e_u)
        assert edge(MAO_2008_ANCHOR_NM) == pytest.approx(MAO_2008_ANCHOR_MM, rel=1e-9)


def test_absorption_edge_saturates_at_the_bulk_length():
    edge = AbsorptionEdge(0.175, bulk_length_mm=2000.0)
    assert edge(1000.0) == pytest.approx(2000.0, rel=1e-3)
    assert edge(700.0) < 2000.0


def test_absorption_edge_is_monotonic_in_wavelength():
    edge = AbsorptionEdge(0.175)
    lengths = [edge(lam) for lam in range(300, 701, 10)]
    assert all(b > a for a, b in zip(lengths, lengths[1:]))


def test_scanned_slope_brackets_the_published_flat_values():
    """The scan range is chosen so that it spans what the literature asserts.

    Published Geant4 setups use flat absorption lengths of 500 mm and 1000 mm at
    the emission peak. The scanned family must contain both, otherwise the
    systematic band would not cover current practice.
    """
    shallow, steep = AbsorptionEdge(0.10)(415.0), AbsorptionEdge(0.25)(415.0)
    assert steep < 500.0 < 1000.0 < shallow


def test_a_bulk_length_shorter_than_the_anchor_is_refused():
    """Refuse the contradiction rather than silently producing a wrong curve."""
    with pytest.raises(ValueError, match="shorter than the measured anchor"):
        AbsorptionEdge(0.175, bulk_length_mm=30.0)(415.0)


def test_the_crystal_is_strongly_absorbing_inside_its_own_emission_band():
    """The physical point: 59 mm at 365 nm, against a 76.2 mm long crystal.

    A flat absorption length of 500-2000 mm cannot represent this, which is why
    the flat variants exist only as a comparison.
    """
    assert AbsorptionEdge(0.175)(365.0) < 76.2
    assert gaussian_emission(365.0, 415.0, 65.0) > 0.02


# --- emission and units ---------------------------------------------------- #


def test_gaussian_emission_peaks_at_one_and_has_the_requested_fwhm():
    assert gaussian_emission(415.0, 415.0, 65.0) == 1.0
    assert gaussian_emission(415.0 - 32.5, 415.0, 65.0) == pytest.approx(0.5, abs=1e-12)
    assert gaussian_emission(415.0 + 32.5, 415.0, 65.0) == pytest.approx(0.5, abs=1e-12)


def test_wavelength_energy_round_trip_matches_the_cpp_constant():
    """The Python and C++ sides must use the same hc, or the tables shift."""
    assert wavelength_to_energy_eV(415.0) == pytest.approx(2.98757, abs=1e-5)
    assert energy_to_wavelength_nm(wavelength_to_energy_eV(415.0)) == pytest.approx(415.0)


# --- the generated material files ------------------------------------------ #


def test_material_files_are_in_sync_with_their_generator():
    """The .dat files are generated; a hand edit or a stale file must be caught.

    Without this, someone can edit a number in a material file, get a result from
    it, and have the provenance comment beside it silently describe a different
    number.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "make_material_variants", root / "scripts" / "make_material_variants.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    stale = []
    for path, variant, note, dispersion, absorption, fwhm, jacobian in module.VARIANTS:
        expected = module.build(
            variant=variant, variant_note=note,
            dispersion=dispersion, absorption=absorption, fwhm_nm=fwhm,
            jacobian_corrected=jacobian,
        )
        if (root / path).read_text() != expected:
            stale.append(path)
    assert not stale, (
        "these material files differ from what the generator produces: "
        f"{stale}. Run scripts/make_material_variants.py."
    )


def test_no_provisional_values_remain_in_the_baseline_material():
    """Every optical input in the baseline must now carry a source.

    PROVISIONAL marked a value whose source had not been secured. The baseline
    file is not allowed to contain one again without this test being changed
    deliberately.
    """
    from pathlib import Path

    text = (Path(__file__).resolve().parents[1] / "materials" / "NaI_Tl.dat").read_text()
    assert "PROVISIONAL" not in text
    # The two inputs that genuinely are assumptions must still say so out loud.
    assert "ASSUMED" in text
    for section in ("[RINDEX]", "[ABSLENGTH_mm]", "[SCINTILLATIONCOMPONENT1]",
                    "[SCINTILLATIONCOMPONENT2]"):
        assert section in text


def test_baseline_absorption_table_dips_below_the_crystal_length():
    """The physical claim the baseline file now makes, checked in the file itself."""
    from pathlib import Path

    text = (Path(__file__).resolve().parents[1] / "materials" / "NaI_Tl.dat").read_text()
    block = text.split("[ABSLENGTH_mm]")[1].split("[")[0]
    rows = [line.split() for line in block.strip().splitlines() if line and line[0].isdigit()]
    table = {float(a): float(b) for a, b in rows}
    assert table[360.0] < 76.2, "the crystal must be opaque to its own blue tail"
    assert table[500.0] > 1000.0, "and transparent on the red side"


def test_the_absorption_anchor_matches_its_own_derivation():
    """MAO_2008_ANCHOR_MM must be what scripts/derive_nai_abslength.py computes.

    Without this the constant in scint/optical.py and the derivation that
    justifies it could drift apart, and the material files would carry a
    provenance note describing a number they no longer contain.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "derive_nai_abslength", root / "scripts" / "derive_nai_abslength.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    n365 = module.n_li1976(MAO_2008_ANCHOR_NM)
    t365 = 0.5 * module.t_theoretical(module.n_li1976(800.0))
    derived = module.abslength_from_transmittance(t365, n365, module.LENGTH_MM)
    assert derived == pytest.approx(MAO_2008_ANCHOR_MM, abs=1.0)


def test_the_absorption_anchor_is_robust_to_the_assumptions_behind_it():
    """The derivation's stated robustness, 55-60 mm, checked rather than asserted.

    The assumptions varied are the ones the paper does not pin down: which
    dispersion model to use, how close the 800 nm transmittance sits to its
    theoretical limit, and which radiation length the authors used for the
    sample size.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "derive_nai_abslength", root / "scripts" / "derive_nai_abslength.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    results = []
    for n_of in (module.n_li1976, module.n_jellison_endpoints):
        for t800_factor in (1.00, 0.98, 0.95):
            for x0 in (9.49, 9.67):
                length = 10.0 * 1.5 * x0 / module.DENSITY
                t365 = 0.5 * t800_factor * module.t_theoretical(n_of(800.0))
                results.append(
                    module.abslength_from_transmittance(t365, n_of(365.0), length)
                )
    assert 54.0 < min(results) < max(results) < 62.0


def test_jacobian_correction_undoes_geant4s_energy_space_sampling():
    """Geant4 integrates the emission table over photon ENERGY, not wavelength.

    G4Scintillation::BuildInverseCdfTable (Geant4 11.4.2) builds the CDF with the
    trapezium weight Energy(i) - Energy(i-1), so a table written down as a curve
    in wavelength -- which is what digitising a published emission figure gives --
    is sampled as p(lambda) proportional to w(lambda)/lambda^2, and comes out
    blue-shifted. This test reproduces Geant4's own construction on both tables
    and checks that the correction puts the mean back where it was written.

    The measured shift in simulation was 411.25 nm against 415 nm written; the
    numbers below are that same effect computed directly.
    """
    import numpy as np

    from scint.optical import PLANCK_C_EV_NM

    lam = np.arange(300.0, 600.001, 5.0)

    def sampled_mean(values: np.ndarray) -> float:
        energy = PLANCK_C_EV_NM / lam
        order = np.argsort(energy)
        e, w = energy[order], values[order]
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * np.diff(e) * (w[1:] + w[:-1]))])
        cdf /= cdf[-1]
        u = np.linspace(0.0, 1.0, 400_001)[1:-1]
        return float((PLANCK_C_EV_NM / np.interp(u, cdf, e)).mean())

    plain = np.array([gaussian_emission(x, 415.0, 65.0) for x in lam])
    fixed = np.array([gaussian_emission_sampled_in_wavelength(x, 415.0, 65.0) for x in lam])

    assert sampled_mean(plain) == pytest.approx(411.24, abs=0.15)
    assert sampled_mean(fixed) == pytest.approx(415.0, abs=0.15)
    # Leading order the shift is -2 sigma^2 / lambda_0.
    sigma = 65.0 / (2 * math.sqrt(2 * math.log(2)))
    assert sampled_mean(plain) - 415.0 == pytest.approx(-2 * sigma**2 / 415.0, abs=0.3)
