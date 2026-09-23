#!/usr/bin/env python3
"""Generate NaI:Tl material specifications from the sourced optical models.

One baseline file (``materials/NaI_Tl.dat``) and a family of variants
(``materials/variants/``), each varying exactly one optical input.  The point of
the family is that three of the inputs are not pinned down by any single
measurement, so the honest output is not a number but a spread: run the family,
report the range.

Run from the repository root:

    python3 scripts/make_material_variants.py

The files are generated, not hand-edited, so that the provenance text in them
cannot drift away from the models in ``scint/optical.py`` that produced them.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scint.optical import (  # noqa: E402
    JELLISON_2012_ENDPOINTS,
    LI_1976,
    AbsorptionEdge,
    MeasuredAttenuation,
    ScatteringSplit,
    TabulatedAttenuation,
    Dispersion,
    constant_absorption,
    constant_dispersion,
    gaussian_emission,
    gaussian_emission_sampled_in_wavelength,
    sample_wavelengths,
    NullAttenuation,
    TabulatedEmission,
)

ROOT = Path(__file__).resolve().parents[1]

# The band is sampled well past the point where the Gaussian has died, because
# Geant4 interpolates linearly between the endpoints of a property table and
# clamps outside it: a table that stops while the value is still appreciable
# produces a hard spectral edge that is an artefact of the table, not physics.
BAND_LOW_NM, BAND_HIGH_NM, BAND_STEP_NM = 300.0, 600.0, 5.0

# The refractive index and absorption tables must span at least the emission
# band; Geant4 refuses to transport a photon whose energy falls outside RINDEX.
OPTICAL_LOW_NM, OPTICAL_HIGH_NM, OPTICAL_STEP_NM = 280.0, 700.0, 10.0

EMISSION_PEAK_NM = 415.0
EMISSION_FWHM_NM = 65.0

HEADER = """\
# NaI:Tl -- material specification for the optical simulation.
#
# GENERATED FILE.  Produced by scripts/make_material_variants.py from the models
# in scint/optical.py.  Do not edit by hand: regenerate instead, so that the
# provenance notes below cannot drift away from the code that made the numbers.
#
# VARIANT: {variant}
{variant_note}
#
# PROVENANCE OF THE NON-OPTICAL PARAMETERS
#
#   src_A = Bonesini, arXiv:2505.06929, Table 1
#   src_B = Miller et al., IEEE TNS 72 197 (2025), arXiv:2403.02668, Table I (Geant4 parameter compilation)
#
# NOTE -- the two sources disagree about this material:
#   light yield   src_A 38,000 ph/MeV   src_B 41,000 ph/MeV      (7.9 %)
#   decay         src_A 250 ns single   src_B 220 ns (96 %) + 1500 ns (4 %)
# The structure of the decay differs, not merely its value. src_B is used here
# because it is the compilation assembled specifically for Geant4 and benchmarked
# against measurement; the disagreement is itself a subject of this study.

name          = {name}
density_g_cm3 = 3.67                    # src_A; PDG gives 3.667 for pure NaI
composition   = Na:0.153373922,I:0.846626078   # exact, from stoichiometry

scintillation_yield_per_MeV = 41000     # src_B
yield_ratio_1               = 0.96      # src_B
yield_ratio_2               = 0.04      # src_B
time_constant_1_ns          = 220       # src_B
time_constant_2_ns          = 1500      # src_B

# RESOLUTIONSCALE inflates the variance of the photon number above Poisson.
# src_B uses 3.50 for this material, stating that the resolution scale "was
# calculated from the energy resolution" -- that is, fitted to the observable it
# would then be used to reproduce. Using it would make simulated energy
# resolution a fit rather than a prediction, so this study keeps it at 1 and
# reports the statistical resolution as a lower bound. The gap to measurement is
# the non-proportionality and transfer term, which Geant4 does not model.
resolution_scale = 1.0
"""


def _wrap(prefix: str, text: str, width: int = 78) -> str:
    """Comment-wrap a provenance string without breaking words."""
    return textwrap.fill(
        " ".join(text.split()),
        width=width,
        initial_indent=prefix,
        subsequent_indent=prefix,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _section(title: str, note: str, rows: list[tuple[float, float]], header: str) -> str:
    body = "\n".join(f"{lam:<6.0f} {value:.6g}" for lam, value in rows)
    return f"\n# --- {title} {'-' * max(0, 72 - len(title))}\n{_wrap('# ', note)}\n[{header}]\n{body}\n"


def build(
    *,
    variant: str,
    variant_note: str,
    dispersion: Dispersion,
    absorption,
    fwhm_nm: float,
    name: str = "NaI_Tl",
    jacobian_corrected: bool = False,
    emission=None,
) -> str:
    optical_lams = sample_wavelengths(OPTICAL_LOW_NM, OPTICAL_HIGH_NM, OPTICAL_STEP_NM)
    # A tabulated band is sampled over its own support, on the same 5 nm grid:
    # SSLG4's table runs to 650 nm and 1.2 % of it (integrated the way Geant4
    # samples, in photon energy) lies above the 600 nm the Gaussian variants
    # stop at. That tail is the least absorbed part of the band, so truncating
    # it would bias the comparison in the direction the variant is testing.
    if emission is not None:
        band_lams = sample_wavelengths(
            min(BAND_LOW_NM, emission._lam[0]), max(BAND_HIGH_NM, emission._lam[-1]),
            BAND_STEP_NM,
        )
    else:
        band_lams = sample_wavelengths(BAND_LOW_NM, BAND_HIGH_NM, BAND_STEP_NM)

    extrapolated = [lam for lam in optical_lams if dispersion.extrapolating(lam)]
    rindex_note = dispersion.source
    if extrapolated:
        rindex_note += (
            f" EXTRAPOLATION: the source claims validity over "
            f"{dispersion.valid_nm[0]:.0f}-{dispersion.valid_nm[1]:.0f} nm, but "
            f"{len(extrapolated)} of the {len(optical_lams)} tabulated points lie "
            f"outside it ({min(extrapolated):.0f}-{max(extrapolated):.0f} nm). "
            f"For this material that is unavoidable -- the emission band reaches "
            f"below every published measurement of its own refractive index."
        )

    emission_note = (
        f"Gaussian band, peak {EMISSION_PEAK_NM:.0f} nm (src_A), "
        f"FWHM {fwhm_nm:.0f} nm. The peak position is well established; the SHAPE "
        f"AND WIDTH ARE ASSUMED. The real band is asymmetric, and no open source "
        f"publishes it as a table -- src_B plots it in Fig. 8b but the underlying "
        f"data sits behind an IEEE DataPort subscription, and the LBNL scintillator "
        f"library lists the peak only. The width is therefore scanned (55/65/75 nm) "
        f"rather than claimed. It is entangled with the absorption edge: the blue "
        f"side of this band runs straight into it. "
        + (
            "JACOBIAN-CORRECTED: the tabulated values are multiplied by "
            "(lambda/lambda_0)^2 so that the distribution Geant4 actually samples is "
            "the Gaussian in wavelength written above. Geant4 integrates the table "
            "over photon energy, so without this factor the sampled band is "
            "blue-shifted by 2 sigma^2 / lambda_0 = 3.8 nm."
            if jacobian_corrected else
            "AS TABULATED: the values are written straight down as a Gaussian in "
            "wavelength, which is what digitising a published emission curve "
            "produces. Geant4 integrates the table over photon ENERGY, so the band "
            "it samples is blue-shifted by 3.8 nm relative to the curve written "
            "here (measured: mean 411.25 nm against the 415 nm written). This is "
            "the baseline because it is what current practice produces; the "
            "jacobian variant shows what the choice costs."
        )
    )

    text = HEADER.format(variant=variant, variant_note=_wrap("# ", variant_note), name=name)
    text += _section(
        "Refractive index", rindex_note,
        [(lam, dispersion(lam)) for lam in optical_lams], "RINDEX",
    )
    text += _section(
        "Bulk absorption length", absorption.source,
        [(lam, absorption(lam)) for lam in optical_lams], "ABSLENGTH_mm",
    )
    scattering_mm = getattr(absorption, "scattering_mm", None)
    if scattering_mm:
        text += _section(
            "Rayleigh scattering length",
            "Grey (wavelength-independent) scattering at the lower bound the "
            "transmittance plateau allows. Values in mm: Geant4's internal length "
            "unit is the millimetre and G4OpRayleigh reads RAYLEIGH unscaled.",
            [(lam, scattering_mm) for lam in optical_lams], "RAYLEIGH",
        )
    if emission is not None:
        emission_note = emission.source
        band = [(lam, emission(lam)) for lam in band_lams]
    else:
        shape = (gaussian_emission_sampled_in_wavelength if jacobian_corrected
                 else gaussian_emission)
        band = [(lam, shape(lam, EMISSION_PEAK_NM, fwhm_nm)) for lam in band_lams]
    for component in ("SCINTILLATIONCOMPONENT1", "SCINTILLATIONCOMPONENT2"):
        note = emission_note if component.endswith("1") else (
            "The two decay components are given the same emission spectrum: no "
            "source reports a different band for the slow component of NaI:Tl."
        )
        text += _section(f"Emission spectrum ({component[-1]})", note, band, component)
    return text


SSLG4_BAND = TabulatedEmission(key="emission_sslg4", source=(
    "PUBLISHED TABLE, NOT A MEASUREMENT AND NOT THIS STUDY'S GAUSSIAN. The "
    "NaI(Tl) emission component distributed for Geant4 by SSLG4 (Kandemir et "
    "al., Comput. Phys. Commun. 306 (2025) 109385), read from "
    "data/optical/sslg4_nai_emission.csv, which records the upstream commit and "
    "the file's sha256. The library states that the data for its commercially "
    "available scintillators come from the manufacturers' published data "
    "sheets and attributes this entry to Luxium Solutions; no temperature, "
    "excitation, sample or uncertainty travels with it. It is carried as a "
    "variant because it is the only openly distributed machine-readable band "
    "we located and because its shape lies outside what this study assumes: "
    "peak 424.3 nm against the 415 nm every compilation tabulates, full width "
    "at half maximum 52.8 nm against the 55-75 nm scanned here. Tabulated as "
    "distributed, so Geant4 samples it in photon energy exactly as it would a "
    "digitised published curve -- the convention of Appendix B.2 applies to it "
    "unchanged."))

TRANSPARENT = NullAttenuation(source=(
    "NO BULK ATTENUATION AT ALL, which is what one published Geant4 material "
    "library supplies for NaI(Tl): SSLG4 (Comput. Phys. Commun. 306 (2025) "
    "109385) binds RINDEX and SCINTILLATIONCOMPONENT1 for its isc-2004 entry "
    "and leaves the ABSLENGTH line commented out in the material macro, with "
    "the file it names absent from the distribution. The value written is "
    "1e7 mm, longer than any path a photon takes in a 76.2 mm crystal, so "
    "nothing is absorbed in bulk. Not an assumption anyone defends: the top of "
    "what the published record actually hands to Geant4, arrived at the same "
    "way the missing REFLECTIVITY is -- not by a decision but by a blank."))

BASELINE_EDGE = AbsorptionEdge(urbach_energy_eV=0.175)

# The measured curve, and its digitisation-error bounds. These are not a scan:
# the slope that the AbsorptionEdge family had to guess at is measured here, so
# what the hi/lo pair spans is the error on a measurement rather than the range
# of models somebody might choose.
FLAT_185 = constant_dispersion(1.85, source=(
    "The single value quoted for NaI:Tl on essentially every manufacturer "
    "datasheet and in most simulation papers. Kept as a variant because it is "
    "what the field does, not because it is a measurement."))

MEASURED = MeasuredAttenuation()
# The other reading of the same curve: as much of the attenuation as the
# plateau allows assigned to scattering, which redirects light instead of
# removing it. 4800 mm is the bound Section 2.2 of the paper states (4795 mm
# before rounding); the paper quotes the rounded figure, so the run uses it.
MEASURED_SCATTER_MAX = ScatteringSplit(MEASURED, 4800.0, source=(
    "The measured curve (Mao et al. Fig. 2, digitised) with grey scattering at "
    "the 4.8 m lower bound the 750-800 nm plateau allows, and absorption "
    "carrying the remainder: 703 mm at 415 nm against 614 mm total. The "
    "extreme of the absorption/scattering ambiguity a single-beam "
    "transmittance leaves; not a measurement of either."))
MEASURED_HI = MeasuredAttenuation(sigma=+1.0)
MEASURED_LO = MeasuredAttenuation(sigma=-1.0)
# The wavelength axis is the other half of the digitisation error, and below
# 390 nm it is the larger half. Kept as a separate pair because the two axes
# are independent: they are combined in quadrature at the end, not stacked.
MEASURED_LAMHI = MeasuredAttenuation(sigma_lam=+1.0)
MEASURED_LAMLO = MeasuredAttenuation(sigma_lam=-1.0)
# The index variants invert the transmittance with their OWN index, so that a
# file's ABSLENGTH and RINDEX blocks are one description rather than two.
MEASURED_JELLISON = MeasuredAttenuation(dispersion=JELLISON_2012_ENDPOINTS)
MEASURED_FLAT185 = MeasuredAttenuation(dispersion=FLAT_185)

FLAT_2000 = constant_absorption(
    2000.0,
    source=(
        "Flat 2000 mm -- effectively transparent. This is the assumption this "
        "study started from, kept as a variant to measure what it costs. It "
        "cannot be right: Mao et al. (2008) measured the crystal to transmit only "
        "50 % at 365 nm over 39 mm, which is inside the emission band."
    ),
)
FLAT_1000 = constant_absorption(
    1000.0,
    source=(
        "Flat 1000 mm. A value that appears in published Geant4 scintillator "
        "setups. Flat absorption cannot represent an absorption edge that cuts "
        "into the emission band; included to compare against the edge model."
    ),
)
# The one structured NaI(Tl) attenuation model found implemented in a published
# Geant4 simulation, read from the corrected figure of its corrigendum. Not a
# measurement and not a bound: the comparison of Section 2 of the paper, run.
BROWN_CORRECTED = TabulatedAttenuation(
    key="brown2023",
    source=(
        "As implemented by J. M. C. Brown, Appl. Radiat. Isot. 168 (2021) 109368, "
        "read from the corrected Fig. 1 of the corrigendum, Appl. Radiat. Isot. 194 "
        "(2023) 110721 (data/optical/brown2021_nai_attenuation.csv). Brown's Table "
        "A3 credits the curve to Mao et al. (2008) without stating the derivation; "
        "it is about 1.8 times shorter than the inversion of that measurement at "
        "415 nm. Carried so that the difference can be stated in light collection."
    ),
)
FLAT_500 = constant_absorption(
    500.0,
    source=(
        "Flat 500 mm. The other value that appears in published Geant4 "
        "scintillator setups. Same caveat as the 1000 mm variant."
    ),
)

VARIANTS: list[tuple[str, str, str, Dispersion, object, float, bool]] = [
    # (path, variant key, note, dispersion, absorption, emission FWHM)
    (
        "materials/NaI_Tl.dat", "baseline",
        "Li 1976 dispersion + measured-anchor absorption edge (E_U = 0.175 eV, the "
        "midpoint of the scanned range) + 65 nm Gaussian band. Li rather than "
        "Jellison is the baseline because it is what the field's reference "
        "compilation actually uses, so this baseline is comparable with published "
        "work; Jellison is carried as the systematic variation.",
        LI_1976, BASELINE_EDGE, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_rindex_jellison.dat", "rindex = Jellison 2012",
        "Refractive index from the 2012 direct measurement instead of the 1976 fit.",
        JELLISON_2012_ENDPOINTS, BASELINE_EDGE, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_rindex_flat185.dat", "rindex = flat 1.85",
        "The datasheet value, with no dispersion at all.",
        constant_dispersion(1.85, source=(
            "The single value quoted for NaI:Tl on essentially every manufacturer "
            "datasheet and in most simulation papers. Kept as a variant because it "
            "is what the field does, not because it is a measurement.")),
        BASELINE_EDGE, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_edge010.dat", "absorption edge E_U = 0.10 eV",
        "Shallow edge: most transparent member of the scanned family.",
        LI_1976, AbsorptionEdge(0.10), EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_edge025.dat", "absorption edge E_U = 0.25 eV",
        "Steep edge: most absorbing member of the scanned family.",
        LI_1976, AbsorptionEdge(0.25), EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_rindex_jellison.dat",
        "measured absorption + Jellison 2012 refractive index",
        "The refractive-index alternative, on top of the measured attenuation.",
        JELLISON_2012_ENDPOINTS, MEASURED_JELLISON, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_rindex_flat185.dat",
        "measured absorption + flat refractive index 1.85",
        "The flat-index practice, on top of the measured attenuation, "
        "inverted with the same flat index.",
        FLAT_185, MEASURED_FLAT185, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_fwhm55.dat",
        "measured absorption + emission FWHM 55 nm",
        "The narrow end of the emission-band scan, on the measured attenuation.",
        LI_1976, MEASURED, 55.0, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_fwhm75.dat",
        "measured absorption + emission FWHM 75 nm",
        "The wide end of the emission-band scan, on the measured attenuation.",
        LI_1976, MEASURED, 75.0, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_jacobian.dat",
        "measured absorption + Jacobian-corrected emission table",
        "The tabulation convention of Appendix B.2, on the measured attenuation, "
        "so that the measurement-consistent family counts the same term the "
        "scanned family does.",
        LI_1976, MEASURED, EMISSION_FWHM_NM, True,
    ),
    (
        "materials/variants/NaI_Tl_abs_measured.dat",
        "absorption = measured (Mao et al. Fig. 2, digitised)",
        "The published transmittance curve, inverted. Replaces the scanned "
        "Urbach slope with a measurement.",
        LI_1976, MEASURED, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_scatter_max.dat",
        "absorption = measured, with the maximal scattering share the plateau allows",
        "Same total attenuation as the measured curve, but the largest part of "
        "it that the 800 nm plateau permits (L_s = 4.8 m, grey) is scattering "
        "rather than absorption. Bounds what the absorption-only implementation "
        "of the curve can overstate.",
        LI_1976, MEASURED_SCATTER_MAX, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_meas_emission_sslg4.dat",
        "measured absorption + the SSLG4 tabulated emission band",
        "The one openly distributed machine-readable NaI(Tl) emission band, on "
        "the measured attenuation. Replaces this study's assumed Gaussian with "
        "a published table whose peak and width both lie outside the assumed "
        "band, so what it measures is the cost of the shape assumption rather "
        "than of the width scan.",
        LI_1976, MEASURED, EMISSION_FWHM_NM, False, SSLG4_BAND,
    ),
    (
        "materials/variants/NaI_Tl_abs_none.dat",
        "no bulk attenuation (SSLG4 supplies none)",
        "A crystal transparent by omission, which is what a user of the SSLG4 "
        "NaI(Tl) entry gets. The top of the range of absorption models the "
        "published record hands to Geant4, and unlike the flat 2000 mm at the "
        "old top of that range it is somebody else's.",
        LI_1976, TRANSPARENT, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_measured_hi.dat",
        "absorption = measured, +1 sigma of the transmittance digitisation error",
        "Transmittance-axis upper bound of the measured curve.",
        LI_1976, MEASURED_HI, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_measured_lo.dat",
        "absorption = measured, -1 sigma of the transmittance digitisation error",
        "Transmittance-axis lower bound of the measured curve.",
        LI_1976, MEASURED_LO, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_measured_lamhi.dat",
        "absorption = measured, wavelength axis +1 sigma",
        "Wavelength-axis bound: the curve moved 3.72 nm to the red, which is "
        "the residual left by the third calibration handle. At 365 nm this "
        "reproduces the analytic inversion of the stated cut-off.",
        LI_1976, MEASURED_LAMHI, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_measured_lamlo.dat",
        "absorption = measured, wavelength axis -1 sigma",
        "Wavelength-axis bound: the curve moved 3.72 nm to the blue.",
        LI_1976, MEASURED_LAMLO, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_brown.dat",
        "absorption = as implemented by Brown (2021), corrected figure (2023)",
        "The structured model the Geant4 optical-property compilation of Miller "
        "et al. redistributes, with everything else at the baseline, so that its "
        "distance from the measured curve is a number in light collection.",
        LI_1976, BROWN_CORRECTED, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_flat2000.dat", "absorption = flat 2000 mm",
        "No absorption edge -- the assumption this study started from.",
        LI_1976, FLAT_2000, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_flat1000.dat", "absorption = flat 1000 mm",
        "No absorption edge -- a published flat value.",
        LI_1976, FLAT_1000, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_abs_flat500.dat", "absorption = flat 500 mm",
        "No absorption edge -- the other published flat value.",
        LI_1976, FLAT_500, EMISSION_FWHM_NM, False,
    ),
    (
        "materials/variants/NaI_Tl_emission_fwhm55.dat", "emission FWHM = 55 nm",
        "Narrower band: less emission overlapping the absorption edge.",
        LI_1976, BASELINE_EDGE, 55.0, False,
    ),
    (
        "materials/variants/NaI_Tl_emission_fwhm75.dat", "emission FWHM = 75 nm",
        "Wider band: more emission overlapping the absorption edge.",
        LI_1976, BASELINE_EDGE, 75.0, False,
    ),
    (
        "materials/variants/NaI_Tl_emission_jacobian.dat",
        "emission band Jacobian-corrected",
        "Same 415 nm / 65 nm Gaussian, but tabulated so that the distribution Geant4 "
        "SAMPLES is that Gaussian in wavelength. Geant4 builds the emission CDF by "
        "integrating the table over photon energy, so writing a wavelength-space curve "
        "down directly -- what digitising a published figure produces -- makes the "
        "sampled band 3.8 nm bluer than the curve drawn. Verified in simulation: the "
        "uncorrected table gives a mean emitted wavelength of 411.25 nm against the "
        "415 nm written, matching the analytic expectation of 411.24 nm.",
        LI_1976, BASELINE_EDGE, EMISSION_FWHM_NM, True,
    ),
]


def main() -> None:
    for entry in VARIANTS:
        path, variant, note, dispersion, absorption, fwhm, jacobian = entry[:7]
        emission = entry[7] if len(entry) > 7 else None
        out = ROOT / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            build(
                variant=variant, variant_note=note,
                dispersion=dispersion, absorption=absorption, fwhm_nm=fwhm,
                jacobian_corrected=jacobian, emission=emission,
            )
        )
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
