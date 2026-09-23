#!/usr/bin/env python3
"""Digitise the NaI(Tl) photoluminescence band from Mao, Zhang & Zhu (2008), Fig. 2.

Why this exists
---------------
The emission band is the one optical input this study had treated as having
no measured shape: the only numerical spectrum located was a manufacturer's
table (SSLG4), and the baseline is an assumed Gaussian. But the figure the
attenuation curve is digitised from carries a measured band as well. Each
panel of Mao et al.'s Fig. 2 draws, beside the transmittance (green), the
excitation (red) and photo-luminescence (blue) spectra of the same sample,
measured with a Hitachi F4500 fluorescence spectrophotometer with the
excitation beam at 10 degrees to the sample normal "so that the
photo-luminescence spectra collected are not affected by internal absorption
in the sample" (Section III of the paper). That is a measured band shape with
its instrument and geometry stated. It is photoluminescence under ultraviolet
excitation, not scintillation under ionising radiation, and the paper does
not say whether the spectra were corrected for the instrument's spectral
response; both caveats travel with the output.

What it does
------------
The page is rendered and calibrated exactly as scripts/digitise_mao_fig2.py
does for the transmittance -- same panel frame, same wavelength calibration,
same offset from the two peak wavelengths the authors print inside the panel
-- so the two curves share one wavelength axis. The blue stroke is then traced
along its centreline: column by column where the stroke is shallow, row by row
where it is steep (both flanks), the same rule the transmittance trace uses.
The intensity is measured from the panel's lower frame rule, which is the
zero of the plot, and normalised to the band's maximum.

Checks printed, all of which must pass:
  * the traced maximum, after the calibration offset, lies within the
    +-0.7 nm the calibration leaves on the printed 410 nm;
  * every gap wider than 3 nm between traced points is one where another
    curve's stroke is drawn over the blue one. There is exactly one: the red
    excitation curve crosses the steep blue flank near 360 nm. The gap is
    written into the output header; Geant4, like the table reader, bridges it
    linearly.

    python3 scripts/digitise_mao_fig2_emission.py --pdf mao2008.pdf
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import digitise_mao_fig2 as tr  # noqa: E402  (the transmittance digitiser)

OUT = ROOT / "data" / "optical" / "mao2008_nai_photoluminescence.csv"

# A stroke run no longer than this locates the curve well along the direction
# it was measured in; the same threshold as the transmittance trace.
NARROW_PX = 14

# The calibration leaves the two printed peaks at -0.7 and +0.7 nm about the
# applied offset; the traced maximum must sit inside that.
PEAK_TOLERANCE_NM = 1.0
MAX_GAP_NM = 3.0
# The one occlusion allowed, and how wide it may be: the red excitation curve
# crosses the blue flank. Anything else is a tracing failure.
MAX_OCCLUDED_GAP_NM = 10.0


def runs(idx: np.ndarray) -> list[tuple[int, int]]:
    out, start, prev = [], int(idx[0]), int(idx[0])
    for v in idx[1:]:
        v = int(v)
        if v - prev > 3:
            out.append((start, prev))
            start = v
        prev = v
    out.append((start, prev))
    return out


def centreline(mask: np.ndarray, panel: dict) -> list[tuple[float, float]]:
    """(x, y) points on the stroke's centreline, in page pixels."""
    pts: list[tuple[float, float]] = []
    top, bottom = panel["top"] + 2, panel["bottom"] - 1
    left, right = panel["left"] + 3, panel["right"] - 3
    for x in range(left, right):
        yy = np.where(mask[top:bottom, x])[0]
        if len(yy):
            for s, e in runs(yy + top):
                if e - s <= NARROW_PX:
                    pts.append((float(x), (s + e) / 2.0))
    for y in range(top, bottom):
        xx = np.where(mask[y, left:right])[0]
        if len(xx):
            for s, e in runs(xx + left):
                if e - s <= NARROW_PX:
                    pts.append(((s + e) / 2.0, float(y)))
    return pts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", type=Path, required=True,
                    help="Mao, Zhang & Zhu, IEEE TNS 55 (2008) 2425")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory() as tmp:
        page = tr.render(args.pdf, Path(tmp))
    m = tr.masks(page)
    panel = tr.PANEL

    raw_nm = lambda x: tr.LAM_REF_NM + (x - tr.LAM_REF_PX) / tr.PX_PER_NM   # noqa: E731
    em = raw_nm(tr.peak_of(m["blue"], panel))
    ex = raw_nm(tr.peak_of(m["red"], panel))
    offset = ((tr.STATED_EMISSION_NM - em) + (tr.STATED_EXCITATION_NM - ex)) / 2
    to_nm = lambda x: raw_nm(x) + offset                                    # noqa: E731
    zero_px = float(panel["bottom"])

    pts = centreline(m["blue"], panel)
    # One value per wavelength bin of one pixel column: the mean of the
    # centreline points that fall in it.
    by_col: dict[int, list[float]] = {}
    for x, y in pts:
        by_col.setdefault(int(round(x)), []).append(y)
    cols = sorted(by_col)
    lam = np.array([to_nm(c) for c in cols])
    height = np.array([zero_px - float(np.mean(by_col[c])) for c in cols])
    peak = float(height.max())
    inten = height / peak

    lam_peak = float(lam[int(np.argmax(inten))])
    print("calibration (shared with the transmittance trace)")
    print(f"  applied offset {offset:+.2f} nm")
    print(f"  traced band maximum at {lam_peak:.1f} nm against the "
          f"{tr.STATED_EMISSION_NM:.0f} nm printed")
    if abs(lam_peak - tr.STATED_EMISSION_NM) > PEAK_TOLERANCE_NM:
        raise SystemExit("the traced maximum is not where the panel says the "
                         "emission peak is; do not use this output")
    gaps = np.diff(lam)
    occluded = []
    for i in np.where(gaps > MAX_GAP_NM)[0]:
        # Between the two traced columns, is the blue stroke's path covered by
        # another colour? Look for red or green pixels in the columns between
        # them, over the rows the flank spans there.
        c0, c1 = cols[i], cols[i + 1]
        y0, y1 = sorted((float(np.mean(by_col[c0])), float(np.mean(by_col[c1]))))
        rows = slice(int(y0) - 4, int(y1) + 5)
        covered = any((m["red"][rows, c] | m["green"][rows, c]).any() for c in range(c0 + 1, c1))
        if not covered or gaps[i] > MAX_OCCLUDED_GAP_NM:
            raise SystemExit(f"gap of {gaps[i]:.1f} nm in the trace at {lam[i]:.1f} nm "
                             "that no other curve covers; the trace has failed there")
        occluded.append((float(lam[i]), float(lam[i + 1])))
        print(f"  occluded by another curve's stroke: {lam[i]:.1f}-{lam[i + 1]:.1f} nm "
              f"({gaps[i]:.1f} nm), bridged linearly")
    if len(occluded) > 1:
        raise SystemExit("more than one occluded gap; expected only the red crossing")

    half = inten >= 0.5
    lo = float(np.interp(0.5, inten[:int(np.argmax(inten)) + 1], lam[:int(np.argmax(inten)) + 1]))
    hi_part = slice(int(np.argmax(inten)), None)
    hi = float(np.interp(0.5, inten[hi_part][::-1], lam[hi_part][::-1]))
    print(f"  traced range {lam.min():.1f}-{lam.max():.1f} nm, "
          f"FWHM {hi - lo:.1f} nm ({lo:.1f}-{hi:.1f})")
    print(f"  intensity at the ends: {inten[0]:.3f} at {lam[0]:.1f} nm, "
          f"{inten[-1]:.3f} at {lam[-1]:.1f} nm")
    del half

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fh.write(
            "# NaI(Tl) photo-luminescence band, digitised from Fig. 2 of\n"
            "#   R. Mao, L. Zhang, R.-Y. Zhu, IEEE Trans. Nucl. Sci. 55 (2008)\n"
            "#   2425-2431, doi:10.1109/TNS.2008.2000776\n"
            "#\n"
            "# Produced by scripts/digitise_mao_fig2_emission.py. The article PDF\n"
            "# is not redistributed here; it was obtained from the authors' page at\n"
            "#   https://www.its.caltech.edu/~rzhu/papers/08_tns_crystal.pdf\n"
            f"#   sha256 {digest}\n"
            "#\n"
            "# The blue curve of the NaI(Tl) panel: photo-luminescence measured\n"
            "# with a Hitachi F4500 fluorescence spectrophotometer, excitation beam\n"
            "# at 10 degrees to the sample normal so that internal absorption does\n"
            "# not shape the spectrum (Section III of the paper). Ultraviolet\n"
            "# excitation, not ionising radiation; the paper does not state whether\n"
            "# the spectrum was corrected for the instrument's spectral response.\n"
            "# Intensity is measured from the panel's lower frame rule and\n"
            "# normalised to the band maximum; it is an ordinate per unit\n"
            "# wavelength as plotted, whatever the instrument's own convention.\n"
            "#\n"
            f"# calibration offset applied: {offset:+.2f} nm, shared with\n"
            "#   data/optical/mao2008_nai_transmittance.csv\n"
            f"# traced maximum {lam_peak:.1f} nm (printed in the panel: "
            f"{tr.STATED_EMISSION_NM:.0f} nm); FWHM {hi - lo:.1f} nm\n"
            + "".join(f"# occluded by the red excitation curve between {a:.1f} and "
                      f"{b:.1f} nm; no points there, bridged linearly on reading\n"
                      for a, b in occluded)
            + "# the band as traced starts at the first blue pixel; the paper's\n"
            "#   spectrum is not drawn below it\n"
            "#\n"
        )
        w = csv.writer(fh)
        w.writerow(["wavelength_nm", "intensity"])
        for lm, v in zip(lam, inten):
            w.writerow([f"{lm:.2f}", f"{v:.5f}"])
    try:
        shown = args.out.relative_to(ROOT)
    except ValueError:
        shown = args.out
    print(f"\nwrote {shown}  ({len(lam)} points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
