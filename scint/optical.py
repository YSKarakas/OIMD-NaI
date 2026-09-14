"""Optical input models for scintillator material specifications.

Geant4 takes refractive index, bulk absorption length and emission spectrum as
*inputs*.  They are not derived from chemistry, so every one of them is a piece
of literature that has to be traced, and the traceability is the point of this
module: each model below carries its source in its docstring, and the generated
material files repeat it, so that a number in a result can always be walked back
to a published measurement.

Three findings drove the design, all recorded in ``docs/05_OPTIK_GIRDILER.md``:

1.  The refractive index of NaI used throughout the scintillator-simulation
    literature descends from a *single measurement made in 1923*.  Jellison et
    al. (2012) open their abstract by saying so.  Their own measurement disagrees
    with it by 1.7 % at the NaI:Tl emission peak.
2.  Nobody publishes a machine-readable NaI:Tl emission spectrum.  The band shape
    is therefore an assumption, and is treated as one.
3.  The published bulk absorption lengths (500 mm, 1000 mm, "transparent") are
    flat numbers, but NaI:Tl has an absorption edge that cuts *into* its own
    emission band -- Mao et al. (2008) measured the crystal to be 50 % opaque at
    365 nm.  A flat absorption length cannot represent that.

Because none of these three inputs is pinned down by a single unambiguous
measurement, every model here is one member of a family, and the intended use is
to run the family and report the spread.  See ``scripts/make_material_variants.py``.

Wavelengths are in nm, energies in eV, lengths in mm throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

__all__ = [
    "PLANCK_C_EV_NM",
    "wavelength_to_energy_eV",
    "energy_to_wavelength_nm",
    "Dispersion",
    "LI_1976",
    "JELLISON_2012_ENDPOINTS",
    "constant_dispersion",
    "DISPERSIONS",
    "AbsorptionEdge",
    "constant_absorption",
    "MAO_2008_ANCHOR_NM",
    "MAO_2008_ANCHOR_MM",
    "gaussian_emission",
    "gaussian_emission_sampled_in_wavelength",
    "sample_wavelengths",
]

# CODATA 2018 hc in eV nm, as used by the C++ side (MaterialLibrary.cc).
PLANCK_C_EV_NM = 1239.841984


def wavelength_to_energy_eV(wavelength_nm: float) -> float:
    return PLANCK_C_EV_NM / wavelength_nm


def energy_to_wavelength_nm(energy_eV: float) -> float:
    return PLANCK_C_EV_NM / energy_eV


# --------------------------------------------------------------------------- #
# Refractive index
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Dispersion:
    """A refractive index model n(lambda), with its provenance attached.

    ``valid_nm`` is the range over which the source claims validity.  Calling
    outside it is allowed -- the emission band of NaI:Tl reaches below every
    published measurement of its own refractive index, so extrapolation is
    unavoidable -- but :meth:`extrapolating` reports when it happens so that the
    material file can say so.
    """

    key: str
    label: str
    source: str
    valid_nm: tuple[float, float]
    _n: Callable[[float], float]

    def __call__(self, wavelength_nm: float) -> float:
        return self._n(wavelength_nm)

    def extrapolating(self, wavelength_nm: float) -> bool:
        low, high = self.valid_nm
        return not (low <= wavelength_nm <= high)


def _sellmeier(wavelength_nm: float, constant: float, terms: Sequence[tuple[float, float]]) -> float:
    """n^2 - 1 = constant + sum_i B_i lambda^2 / (lambda^2 - C_i^2), lambda in um."""
    lam2 = (wavelength_nm / 1000.0) ** 2
    n2 = 1.0 + constant
    for b, c in terms:
        n2 += b * lam2 / (lam2 - c * c)
    if n2 <= 0:
        raise ValueError(f"Sellmeier gave n^2 = {n2} at {wavelength_nm} nm")
    return math.sqrt(n2)


LI_1976 = Dispersion(
    key="li1976",
    label="Li 1976 Sellmeier",
    source=(
        "H. H. Li, 'Refractive index of alkali halides and its wavelength and "
        "temperature derivatives', J. Phys. Chem. Ref. Data 5, 329-528 (1976), "
        "doi:10.1063/1.555536; coefficients as tabulated by refractiveindex.info "
        "(database/data/main/NaI/nk/Li.yml, CC0), verified by reproducing the "
        "classical sodium-D-line value n(589.3 nm) = 1.7745. "
        "The database entry itself records that the fit 'is based on a single "
        "data point'; Jellison et al. (2012) identify that point as Spangenberg, "
        "Z. Kristallogr. 57, 494-534 (1923). "
        "This is the model used by Mao, Zhang and Zhu, IEEE TNS 55, 2425 (2008), "
        "which in turn is the source cited for NaI:Tl by the Geant4 parameter "
        "compilation of Miller et al., IEEE TNS 72 197 (2025), arXiv:2403.02668 (Fig. 8b)."
    ),
    valid_nm=(250.0, 40000.0),
    _n=lambda lam: _sellmeier(lam, 0.478, ((1.532, 0.170), (4.27, 86.21))),
)


JELLISON_2012_ENDPOINTS = Dispersion(
    key="jellison2012",
    label="Jellison 2012 (published Sellmeier fit)",
    source=(
        "G. E. Jellison Jr., L. A. Boatner, J. O. Ramey, J. A. Kolopus, "
        "L. A. Ramey, D. J. Singh, 'Refractive index of sodium iodide', "
        "J. Appl. Phys. 111, 043521 (2012), doi:10.1063/1.3689746. "
        "Minimum-deviation measurement at six wavelengths from 436 to 633 nm, "
        "fitted by the authors to a one-term Sellmeier form with chi^2 = 1.02. "
        "The coefficients used here, n^2 = 1 + 1.994 lambda^2/(lambda^2 - 0.176^2) "
        "with lambda in um, are the paper's fit as transcribed by "
        "refractiveindex.info (database/data/main/NaI/nk/Jellison.yml). "
        "Table I of the paper (consulted 14 September 2026) lists 435.8 nm 1.839, "
        "488.0 1.814, 514.5 1.804, 546.1 1.799, 578.0 1.786, 633.0 1.778, all "
        "+- 0.002; the fit reproduces them with residuals of at most 0.003 (546.1 nm), "
        "four within +- 0.002, and chi^2/ndf = 1.01 against the published 1.02. The "
        "abstract's '633 nm (n = 1.786)' is a typographical error -- it is the "
        "578 nm value. An earlier version of this file, working from the abstract "
        "alone, rejected the fit and pinned a one-term form to the two abstract "
        "values, giving 1.8504 at 415 nm against the fit's 1.8524; corrected."
    ),
    valid_nm=(436.0, 633.0),
    _n=lambda lam: _sellmeier(lam, 0.0, ((1.994, 0.176),)),
)


def constant_dispersion(value: float, *, source: str) -> Dispersion:
    """A wavelength-independent refractive index -- the usual simplification."""
    return Dispersion(
        key=f"flat{value:g}".replace(".", "p"),
        label=f"constant n = {value:g}",
        source=source,
        valid_nm=(0.0, math.inf),
        _n=lambda lam, v=value: v,
    )


FLAT_185 = constant_dispersion(
    1.85,
    source=(
        "The single value quoted for NaI:Tl on essentially every manufacturer "
        "datasheet and in most simulation papers. Kept as a variant because it "
        "is what the field does, not because it is a measurement: it is a "
        "rounded reading of one of the dispersion curves at the emission peak."
    ),
)

DISPERSIONS: dict[str, Dispersion] = {
    d.key: d for d in (LI_1976, JELLISON_2012_ENDPOINTS, FLAT_185)
}


# --------------------------------------------------------------------------- #
# Bulk absorption
# --------------------------------------------------------------------------- #

# The one measured constraint we have on NaI:Tl self-absorption.  Derived in
# scripts/derive_nai_abslength.py from the cut-off wavelength quoted by
# Mao, Zhang and Zhu, IEEE Trans. Nucl. Sci. 55, 2425-2431 (2008), Section III,
# for a sample 1.5 radiation lengths long.  The derivation gives 55-60 mm across
# the full range of assumptions tested (refractive index model, how close the
# 800 nm transmittance sits to its theoretical limit, which X0 the authors used).
MAO_2008_ANCHOR_NM = 365.0
MAO_2008_ANCHOR_MM = 59.0


@dataclass(frozen=True)
class AbsorptionEdge:
    """Bulk absorption length with an exponential edge, anchored to measurement.

    The absorption coefficient is taken to rise exponentially with photon energy,

        alpha(E) = 1/L_bulk + alpha_anchor * exp((E - E_anchor) / E_U)

    which is the Urbach form.  Two things should be said about it plainly.

    *What is measured*: ``alpha_anchor``, from the Mao et al. (2008) cut-off
    wavelength.  The crystal really is half-opaque at 365 nm over 39 mm; a flat
    absorption length of 500 mm or 2000 mm cannot reproduce that.

    *What is not measured*: the slope ``E_U``.  A single anchor point fixes the
    height of the edge, not its steepness, and the absorbing species here is the
    Tl+ activator band rather than the intrinsic NaI band edge at ~5.9 eV, so the
    Urbach form is phenomenology, not theory.  ``E_U`` is therefore a scanned
    systematic: the range 0.10-0.25 eV brackets the flat absorption lengths that
    appear in the literature (roughly 3500 mm down to 300 mm at the 415 nm
    emission peak).  Mao et al.'s statement that their measured transmittance
    "approaches the theoretical limit" disfavours the strongly absorbing end.
    """

    urbach_energy_eV: float
    bulk_length_mm: float = 2000.0
    anchor_nm: float = MAO_2008_ANCHOR_NM
    anchor_length_mm: float = MAO_2008_ANCHOR_MM

    @property
    def key(self) -> str:
        return f"edge{self.urbach_energy_eV:.3f}".replace(".", "p")

    def __call__(self, wavelength_nm: float) -> float:
        e = wavelength_to_energy_eV(wavelength_nm)
        e_anchor = wavelength_to_energy_eV(self.anchor_nm)
        # The anchor is the *total* absorption length at 365 nm, so the edge term
        # carries only what the wavelength-independent bulk term does not.
        alpha_anchor = 1.0 / self.anchor_length_mm - 1.0 / self.bulk_length_mm
        if alpha_anchor <= 0.0:
            raise ValueError(
                f"bulk length {self.bulk_length_mm} mm is shorter than the "
                f"measured anchor {self.anchor_length_mm} mm at {self.anchor_nm} nm"
            )
        alpha_edge = alpha_anchor * math.exp((e - e_anchor) / self.urbach_energy_eV)
        return 1.0 / (1.0 / self.bulk_length_mm + alpha_edge)

    @property
    def source(self) -> str:
        return (
            f"Urbach-form absorption edge, E_U = {self.urbach_energy_eV:.3f} eV "
            f"(ASSUMED -- scanned systematic, not measured), anchored to "
            f"L_abs({self.anchor_nm:.0f} nm) = {self.anchor_length_mm:.0f} mm "
            f"derived in scripts/derive_nai_abslength.py from the 50 % cut-off "
            f"wavelength reported by Mao, Zhang and Zhu, IEEE TNS 55, 2425 (2008), "
            f"Section III; far from the edge it saturates at "
            f"L_bulk = {self.bulk_length_mm:.0f} mm."
        )


def constant_absorption(length_mm: float, *, source: str) -> Callable[[float], float]:
    fn: Callable[[float], float] = lambda lam, v=length_mm: v
    fn.source = source  # type: ignore[attr-defined]
    fn.key = f"flat{length_mm:g}mm"  # type: ignore[attr-defined]
    return fn


BROWN_CORRECTED_CURVE = "data/optical/brown2021_nai_attenuation.csv"


class TabulatedAttenuation:
    """Attenuation length read from a digitised published curve.

    Log-linear in wavelength between the tabulated points and held at the end
    values outside them. This is how the one structured NaI(Tl) attenuation
    model found implemented in a published Geant4 simulation -- Brown (2021),
    as corrected by the 2023 corrigendum -- is carried into a material variant,
    so that its distance from the measurement can be stated in light-collection
    efficiency and not only in millimetres.
    """

    def __init__(self, path: str | None = None, root: str | None = None, *,
                 key: str, source: str):
        import csv as _csv
        import os

        base = root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.path = path or os.path.join(base, BROWN_CORRECTED_CURVE)
        self.key = key
        self.source = source
        lams: list[float] = []
        lens: list[float] = []
        with open(self.path) as fh:
            for row in _csv.DictReader(ln for ln in fh if not ln.startswith("#")):
                lams.append(float(row["wavelength_nm"]))
                lens.append(float(row["attenuation_length_mm"]))
        if not lams:
            raise ValueError(f"no rows in {self.path}")
        order = sorted(range(len(lams)), key=lambda i: lams[i])
        self._lam = [lams[i] for i in order]
        self._log = [math.log10(lens[i]) for i in order]

    def __call__(self, lam_nm: float) -> float:
        import bisect

        if lam_nm <= self._lam[0]:
            return 10.0 ** self._log[0]
        if lam_nm >= self._lam[-1]:
            return 10.0 ** self._log[-1]
        i = bisect.bisect_left(self._lam, lam_nm)
        f = (lam_nm - self._lam[i - 1]) / (self._lam[i] - self._lam[i - 1])
        return 10.0 ** (self._log[i - 1] + f * (self._log[i] - self._log[i - 1]))


# --------------------------------------------------------------------------- #
# Emission spectrum
# --------------------------------------------------------------------------- #


def gaussian_emission(wavelength_nm: float, peak_nm: float, fwhm_nm: float) -> float:
    """Normalised Gaussian emission band, peak value 1.

    A Gaussian in wavelength is an ASSUMPTION, not a measurement.  The real
    NaI:Tl band is asymmetric, and a search of the open literature turns up no
    machine-readable tabulation of it -- every source publishes a picture.  The
    peak position (415 nm) is well established; the width is not, so it is
    scanned.  This matters more than it looks: the short-wavelength side of the
    band runs into the absorption edge at 365 nm, so band width and
    self-absorption are entangled.
    """
    sigma = fwhm_nm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return math.exp(-0.5 * ((wavelength_nm - peak_nm) / sigma) ** 2)


def gaussian_emission_sampled_in_wavelength(
    wavelength_nm: float, peak_nm: float, fwhm_nm: float
) -> float:
    """The table that makes Geant4 SAMPLE a Gaussian in wavelength.

    Geant4 (verified in the 11.4 beta and the 11.4.2 release) builds the emission CDF by integrating the tabulated values
    over photon ENERGY (``G4Scintillation::BuildInverseCdfTable``: the trapezium
    weight is ``Energy(i) - Energy(i-1)``).  The numbers in a
    SCINTILLATIONCOMPONENT table are therefore a probability density in energy,
    not in wavelength.

    That matters because published emission curves are plotted against
    wavelength.  Tabulating such a curve directly -- which is what digitising a
    figure naturally produces -- makes Geant4 sample

        p(lambda) proportional to  w(lambda) / lambda^2

    which is blue-shifted relative to the curve that was drawn.  For a 415 nm
    band of 65 nm FWHM the mean moves to 411.24 nm, a shift of -3.76 nm, or
    -2 sigma^2 / lambda_0 to leading order.  Measured in simulation: 411.25 nm.

    Multiplying by the Jacobian |d lambda / dE| ~ lambda^2 undoes it, so that the
    sampled distribution really is the Gaussian in wavelength that was intended.
    Whether that is what one wants is a physics question -- emission bands from a
    localised centre are often closer to Gaussian in energy -- but it should be a
    choice, not an accident, which is why both forms are carried as variants.
    """
    return gaussian_emission(wavelength_nm, peak_nm, fwhm_nm) * (wavelength_nm / peak_nm) ** 2


def sample_wavelengths(low_nm: float, high_nm: float, step_nm: float) -> list[float]:
    n = int(round((high_nm - low_nm) / step_nm))
    return [low_nm + i * step_nm for i in range(n + 1)]


# --------------------------------------------------------------------------- #
# Measured attenuation, digitised from a published transmittance curve.
#
# Everything above this line is a model. This is a measurement, and it replaces
# the scanned Urbach slope that used to dominate the model envelope: the slope
# is no longer a free parameter because the curve it describes is published,
# and scripts/digitise_mao_fig2.py turns that picture back into numbers.
#
# Three honest limits travel with it. The quantity is an EFFECTIVE ATTENUATION
# length -- a single-beam transmittance cannot separate absorption from
# scattering. The inversion is well conditioned only where the crystal is not
# yet transparent: above about 450 nm the measured transmittance sits close to
# its Fresnel limit, so a digitisation error of half a percent in T moves the
# attenuation length by tens of percent. And below about 400 nm the curve is so
# steep that the WAVELENGTH axis, not the transmittance axis, carries the error
# -- at 365 nm a 3.7 nm shift is worth a factor 1.8 in L_att, where half a
# percent in T is worth 2.6 %. Both axes are therefore propagated, coherently
# and independently; `uncertainty_mm` combines them and reports either
# component on request.
# --------------------------------------------------------------------------- #

MAO_2008_CURVE = "data/optical/mao2008_nai_transmittance.csv"

# The stroke of the plotted curve, measured by digitise_mao_fig2.py with the
# same colour mask that traces it, at the 400 dpi it renders: 8 px tall on the
# transparent plateau (5.45 px per % T) and 6 px wide on the steep rise
# (0.711 px per nm). An earlier version of this file assumed a 6 px stroke and
# carried 0.0055; a referee measured it. Both numbers below are HALF-WIDTHS of
# that stroke, taken as one-sigma coherent shifts of the whole curve. That is
# deliberately conservative -- the tracer reads the centre of the stroke, and
# the GUM standard uncertainty of a rectangular half-width would be a factor
# sqrt(3) smaller -- and it is stated as a half-width, not as a fitted sigma.
STROKE_PX_VERTICAL, STROKE_PX_HORIZONTAL = 8, 6
T_PX_PER_PCT, PX_PER_NM = 5.45, 0.711
DIGITISATION_SIGMA_T = STROKE_PX_VERTICAL / 2 / T_PX_PER_PCT / 100      # 0.0073

# Wavelength axis. Two independent terms, combined in quadrature: the reading
# half-width above (4.2 nm), and the calibration residual. Three handles fix
# the axis -- the emission peak printed in the panel, the excitation peak
# printed in the panel, and the 50 % cut-off the authors state in their text.
# About the applied offset (the mean of the two in-panel handles) they leave
# -0.7, +0.7 and +3.7 nm. Two of the three agree, so the 3.7 nm is not evidence
# of an axis error; it is a discrepancy between the plotted curve and the
# stated cut-off, of unknown origin, and it is carried in full rather than
# averaged away. Shifting the curve by it puts the half-transmittance point at
# the stated 365 nm by construction -- an identity, not a check.
CALIBRATION_RESIDUAL_NM = 3.7
DIGITISATION_SIGMA_LAM_NM = math.hypot(STROKE_PX_HORIZONTAL / 2 / PX_PER_NM,
                                       CALIBRATION_RESIDUAL_NM)             # 5.6 nm


def _running_median(values: list[float], window: int = 9) -> list[float]:
    half = window // 2
    out = []
    for i in range(len(values)):
        lo, hi = max(0, i - half), min(len(values), i + half + 1)
        chunk = sorted(values[lo:hi])
        out.append(chunk[len(chunk) // 2])
    return out


class MeasuredAttenuation:
    """Effective attenuation length interpolated from a digitised curve.

    `sigma` shifts the curve vertically by that many standard deviations of the
    transmittance digitisation error; `sigma_lam` shifts the wavelength axis by
    that many standard deviations of the axis calibration residual. Both are
    deliberately coherent shifts rather than random ones: a digitisation bias
    moves every point the same way, and that is the conservative assumption for
    an envelope. The two are independent of each other -- one is set by the
    stroke width, the other by the axis handles -- so the variants are generated
    separately and combined in quadrature, not stacked.
    """

    def __init__(self, path: str | None = None, sigma: float = 0.0,
                 root: str | None = None, sigma_lam: float = 0.0,
                 dispersion: "Dispersion | None" = None):
        import csv as _csv
        import os

        base = root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.path = path or os.path.join(base, MAO_2008_CURVE)
        self.sigma = sigma
        self.sigma_lam = sigma_lam
        # The inversion needs n(lambda) for the Fresnel terms. A material
        # variant that carries a different RINDEX must invert with that same
        # index, or its ABSLENGTH and RINDEX blocks disagree with each other
        # by about 2 % at the edge -- a referee's point. Default: Li 1976.
        self.dispersion = dispersion or LI_1976
        lams: list[float] = []
        trans: list[float] = []
        with open(self.path) as fh:
            for row in _csv.DictReader(ln for ln in fh if not ln.startswith("#")):
                lams.append(float(row["wavelength_nm"]))
                trans.append(float(row["transmittance"]))
        if not lams:
            raise ValueError(f"no rows in {self.path}")
        # A wavelength-axis error says the true wavelength of every traced point
        # is displaced by the same amount, so the axis moves and the inversion
        # follows it -- n(lambda) is evaluated where the point actually sits.
        lams = [lm + sigma_lam * DIGITISATION_SIGMA_LAM_NM for lm in lams]
        shifted = [min(0.999, max(1e-6, t + sigma * DIGITISATION_SIGMA_T))
                   for t in trans]
        lengths = [self._invert(lm, t, dispersion=self.dispersion) for lm, t in zip(lams, shifted)]
        smoothed = _running_median(lengths)
        # The curve is physically monotone in this range; pixel quantisation is
        # not. Enforcing it keeps the interpolation from wobbling.
        for i in range(1, len(smoothed)):
            smoothed[i] = max(smoothed[i], smoothed[i - 1])
        self._lams, self._lengths = lams, smoothed
        self.source = (
            "Effective attenuation length from the NaI(Tl) transmittance curve of "
            "Mao, Zhang & Zhu, IEEE Trans. Nucl. Sci. 55 (2008) 2425, Fig. 2, "
            "digitised by scripts/digitise_mao_fig2.py and inverted with the "
            "paper's own transmittance expression over its 38.8 mm sample. "
            "Calibration was checked against the emission and excitation peaks "
            "printed inside the panel, and the recovered 50 % cut-off against "
            "the 365 nm stated in the text. This is an EFFECTIVE attenuation: a "
            "single-beam transmittance does not separate absorption from "
            "scattering, so implementing it as ABSLENGTH overestimates the loss "
            "by whatever fraction is scattering."
        )
        if sigma:
            self.source += (
                f" Shifted by {sigma:+.0f} sigma of the transmittance digitisation "
                f"error ({DIGITISATION_SIGMA_T:.4f} in T) to bound the curve."
            )
        if self.dispersion is not LI_1976:
            self.source += f" Inverted with the {self.dispersion.key} dispersion, matching this variant's RINDEX."
        if sigma_lam:
            self.source += (
                f" Wavelength axis shifted by {sigma_lam:+.0f} sigma of the "
                f"calibration residual ({DIGITISATION_SIGMA_LAM_NM:.2f} nm) to "
                f"bound the curve."
            )

    @staticmethod
    def _invert(lam_nm: float, transmittance: float,
                length_mm: float = 38.8, dispersion=None) -> float:
        n = (dispersion or LI_1976)(lam_nm)
        refl = ((n - 1.0) / (n + 1.0)) ** 2
        qa = transmittance * refl ** 2
        qb = (1.0 - refl) ** 2
        qc = -transmittance
        disc = qb * qb - 4.0 * qa * qc
        a = (-qb + math.sqrt(disc)) / (2.0 * qa)
        a = min(1.0 - 1e-12, max(1e-12, a))
        # Above the Fresnel limit the inversion is unbounded: a shifted-up
        # transmittance can exceed what a lossless crystal transmits, and the
        # length then runs to 1e13 mm. Cap at 100 m, which is transparent for
        # every purpose here and says so in the material file.
        return min(-length_mm / math.log(a), 1.0e5)

    def uncertainty_components_mm(self, wavelength_nm: float) -> dict[str, float]:
        """Half-spreads from each axis separately, and their quadrature sum.

        Reported separately because which one dominates flips inside the
        emission band: the wavelength axis below about 390 nm, where the curve
        is steep, and the transmittance axis above about 430 nm, where it is
        flat and close to the Fresnel limit.
        """
        def half(**kw):
            hi = MeasuredAttenuation(self.path, dispersion=self.dispersion, **{k: +1.0 for k in kw})(wavelength_nm)
            lo = MeasuredAttenuation(self.path, dispersion=self.dispersion, **{k: -1.0 for k in kw})(wavelength_nm)
            return abs(hi - lo) / 2.0

        d_t = half(sigma=1.0)
        d_lam = half(sigma_lam=1.0)
        return {"transmittance": d_t, "wavelength": d_lam,
                "combined": math.hypot(d_t, d_lam)}

    def uncertainty_mm(self, wavelength_nm: float) -> float:
        """Both axes combined in quadrature, as a half-spread in mm."""
        return self.uncertainty_components_mm(wavelength_nm)["combined"]

    def bounds_mm(self, wavelength_nm: float) -> dict[str, float]:
        """Asymmetric bounds: the inversion is convex, so +1 sigma and -1 sigma
        move the length by different amounts, and on the steep edge by very
        different amounts. Each axis' up/down excursions are combined in
        quadrature on their own side; the symmetric half-spread hides this."""
        c = self(wavelength_nm)
        up, dn = 0.0, 0.0
        for kw in ({"sigma": +1.0}, {"sigma": -1.0}, {"sigma_lam": +1.0}, {"sigma_lam": -1.0}):
            v = MeasuredAttenuation(self.path, dispersion=self.dispersion, **kw)(wavelength_nm) - c
            if v > 0: up = math.hypot(up, v)
            else: dn = math.hypot(dn, -v)
        return {"central": c, "lower": c - dn, "upper": c + up}

    def __call__(self, wavelength_nm: float) -> float:
        lams, lengths = self._lams, self._lengths
        if wavelength_nm <= lams[0]:
            return lengths[0]
        if wavelength_nm >= lams[-1]:
            return lengths[-1]
        import bisect

        i = bisect.bisect_left(lams, wavelength_nm)
        x0, x1 = lams[i - 1], lams[i]
        y0, y1 = lengths[i - 1], lengths[i]
        f = (wavelength_nm - x0) / (x1 - x0)
        return y0 + f * (y1 - y0)
