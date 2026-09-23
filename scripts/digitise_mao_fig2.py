#!/usr/bin/env python3
"""Digitise the NaI(Tl) transmittance curve from Mao, Zhang & Zhu (2008), Fig. 2.

Why this exists
---------------
The absorption model in this study was anchored to a single number quoted in
that paper's text -- the wavelength at which transmittance falls to half its
800 nm value, 365 nm for NaI(Tl) -- and its *slope* across the emission band
was a scanned parameter, because no measurement of the slope was available.
The slope is what dominates the model envelope. The curve is published; it is
just published as a picture. This script turns the picture back into numbers.

What it does
------------
1. Renders page 3 of the paper (Fig. 2) at 400 dpi.
2. Locates the NaI(Tl) panel from the axis rules.
3. Calibrates the wavelength axis from the printed tick labels (250, 500,
   750 nm) and the transmittance axis from the green tick marks (20, 40, 60,
   80 %).
4. VALIDATES that calibration against two numbers printed inside the panel by
   the authors themselves -- the emission peak (410 nm) and the excitation
   peak (346 nm) -- and applies the residual offset.
5. Traces the green transmittance stroke. The rise is nearly vertical, so it
   is traced row-wise (x as a function of T) and the shallow parts column-wise
   (T as a function of x); tracing the rise column-wise flattens it and was
   the first thing this script got wrong.
6. Inverts T to an effective attenuation length using the same expression the
   paper uses for its own theoretical limit.

Three checks are printed, and all must pass for the output to mean anything:
the calibration check in step 4; the recovered cut-off wavelength against the
365 nm the paper states in its text, within a tolerance chosen so that the
rejected column-wise trace fails it; and the residual that check leaves against
the wavelength-axis systematic the analysis propagates. The last of these is
the point of the exercise. The residual is 3.7 nm, it is not negligible -- at
365 nm it is worth a factor 1.8 in attenuation length, against 2.6 % for the
transmittance error -- and it is carried as an uncertainty rather than
explained away.

The paper's PDF is not redistributed with this repository. It is available
from the authors' institutional page; the URL and checksum are recorded in
data/optical/mao2008_nai_transmittance.csv.

    python3 scripts/digitise_mao_fig2.py --pdf mao2008.pdf
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scint.optical import LI_1976  # noqa: E402

OUT = ROOT / "data" / "optical" / "mao2008_nai_transmittance.csv"

# Page 3 of the article carries Fig. 2; 400 dpi puts the NaI(Tl) panel at
# roughly 494 x 540 px, which is the resolution the numbers below refer to.
FIG_PAGE, DPI = 2, 400

# Panel frame, found from the long dark rules (see _find_frame).
PANEL = dict(left=2146, right=2640, top=423, bottom=963)

# Axis calibration, from the printed labels.
LAM_REF_PX, LAM_REF_NM = 2249.5, 250.0     # centre of the "250" label
PX_PER_NM = (2605.0 - 2249.5) / 500.0      # "250" to "750"
T_ZERO_PX, T_PX_PER_PCT = 968.0, 545.0/100.0   # green ticks at 20/40/60/80 %

# The authors print these two peak wavelengths inside the panel. They are the
# calibration check: whatever offset reproduces them is applied to everything.
STATED_EMISSION_NM, STATED_EXCITATION_NM = 410.0, 346.0

# Section III: "All samples, except the NaI:Tl sample, have a cubic shape with
# a dimension of 1.5 X0. The NaI:Tl sample is a cylinder of 1.5 X0 long".
SAMPLE_LENGTH_MM = 38.8

# The paper's own stated cut-off for NaI(Tl), used as the second check.
STATED_CUTOFF_NM = 365.0

# How far the recovered cut-off may sit from the stated one. This is not a
# comfort margin: it is the discriminator between this tracing method and the
# one that was rejected. Traced row-wise, the residual is 3.7 nm; traced
# column-wise throughout -- the first thing this script got wrong -- it was
# 8.0 nm. A tolerance of 5 nm passes the first and fails the second, which is
# the whole point of having the check. Whatever residual does survive is not
# waved away either: it is carried into the analysis as the wavelength-axis
# systematic DIGITISATION_SIGMA_LAM_NM in scint/optical.py, and that constant
# and this measurement are checked against each other below.
CUTOFF_TOLERANCE_NM = 5.0


def render(pdf: Path, workdir: Path) -> np.ndarray:
    from PIL import Image

    stem = workdir / "page"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-png",
                    "-f", str(FIG_PAGE), "-l", str(FIG_PAGE), str(pdf), str(stem)],
                   check=True, capture_output=True)
    hits = sorted(workdir.glob("page*.png"))
    if not hits:
        raise SystemExit("pdftoppm produced no image; is poppler installed?")
    return np.asarray(Image.open(hits[0]).convert("RGB")).astype(int)


def masks(a: np.ndarray):
    r, g, b = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    return dict(
        green=(g > 70) & (g - r > 35) & (g - b > 35),
        blue=(b > 90) & (b - r > 45) & (b - g > 45),
        red=(r > 90) & (r - g > 45) & (r - b > 45),
    )


def peak_of(mask: np.ndarray, panel: dict) -> float:
    """Sub-pixel x of a curve's apex, by parabola through its topmost points."""
    tops = {}
    for x in range(panel["left"] + 3, panel["right"] - 3):
        yy = np.where(mask[panel["top"] + 2:panel["bottom"] - 1, x])[0]
        if len(yy):
            tops[x] = yy.min()
    apex = min(tops, key=lambda k: tops[k])
    xs = sorted(k for k in tops if abs(k - apex) <= 12)
    coef = np.polyfit(np.array(xs, float), np.array([tops[k] for k in xs], float), 2)
    return -coef[1] / (2 * coef[0])


def longest_run(idx: np.ndarray) -> tuple[int, int]:
    runs, start, prev = [], idx[0], idx[0]
    for v in idx[1:]:
        if v - prev > 3:
            runs.append((start, prev))
            start = v
        prev = v
    runs.append((start, prev))
    return max(runs, key=lambda t: t[1] - t[0])


def trace(green: np.ndarray, panel: dict, to_nm, to_pct) -> list[tuple[float, float]]:
    """Trace the stroke: row-wise where it is steep, column-wise where it is not.

    The distinction matters more than it looks. Traced column-wise throughout,
    the near-vertical rise collapses into a handful of columns each holding a
    150-px-tall stroke, and the recovered cut-off comes out 8 nm too blue.
    """
    curve: dict[float, float] = {}
    for y in range(panel["top"] + 1, panel["bottom"]):
        xx = np.where(green[y, panel["left"] + 3:panel["right"] - 3])[0]
        if not len(xx):
            continue
        xx = xx + panel["left"] + 3
        s, e = longest_run(xx)
        if e - s <= 14:                       # steep here: the row locates x well
            curve.setdefault(round(to_nm((s + e) / 2.0), 1), to_pct(y))
    for x in range(panel["left"] + 3, panel["right"] - 3):
        yy = np.where(green[panel["top"] + 2:panel["bottom"] - 1, x])[0]
        if not len(yy):
            continue
        yy = yy + panel["top"] + 2
        s, e = longest_run(yy)
        if e - s <= 14:                       # shallow here: the column locates T well
            curve[round(to_nm(x), 1)] = to_pct((s + e) / 2.0)
    return sorted(curve.items())


def attenuation_length(lam_nm: float, transmittance: float, length_mm: float):
    """Invert the paper's own transmittance expression for the bulk term.

    T = (1-R)^2 a / (1 - R^2 a^2),  a = exp(-L/L_att),  R = ((n-1)/(n+1))^2.

    What comes back is an EFFECTIVE ATTENUATION length: a transmittance
    transmittance cannot separate absorption from scattering.
    """
    n = LI_1976(lam_nm)
    refl = ((n - 1) / (n + 1)) ** 2
    qa, qb, qc = transmittance * refl ** 2, (1 - refl) ** 2, -transmittance
    disc = qb * qb - 4 * qa * qc
    if disc < 0:
        return None
    a = (-qb + math.sqrt(disc)) / (2 * qa)
    if not 0.0 < a < 1.0:
        return None
    return -length_mm / math.log(a)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdf", type=Path, required=True,
                    help="Mao, Zhang & Zhu, IEEE TNS 55 (2008) 2425")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory() as tmp:
        page = render(args.pdf, Path(tmp))
    m = masks(page)

    raw_nm = lambda x: LAM_REF_NM + (x - LAM_REF_PX) / PX_PER_NM     # noqa: E731
    to_pct = lambda y: (T_ZERO_PX - y) / T_PX_PER_PCT                # noqa: E731

    # ---- check 1: the panel's own printed peak wavelengths ----------------
    em = raw_nm(peak_of(m["blue"], PANEL))
    ex = raw_nm(peak_of(m["red"], PANEL))
    offset = ((STATED_EMISSION_NM - em) + (STATED_EXCITATION_NM - ex)) / 2
    print("calibration check against the panel's own labels")
    print(f"  emission   : read {em:6.1f} nm   printed {STATED_EMISSION_NM:.0f} nm"
          f"   ({STATED_EMISSION_NM - em:+.1f})")
    print(f"  excitation : read {ex:6.1f} nm   printed {STATED_EXCITATION_NM:.0f} nm"
          f"   ({STATED_EXCITATION_NM - ex:+.1f})")
    print(f"  applied offset {offset:+.2f} nm")
    if abs(offset) > 6.0:
        raise SystemExit("calibration offset is implausibly large; check the frame")

    to_nm = lambda x: raw_nm(x) + offset                             # noqa: E731
    pts = trace(m["green"], PANEL, to_nm, to_pct)
    lam = np.array([p[0] for p in pts])
    tr = np.array([p[1] for p in pts]) / 100.0
    keep = (lam >= 300) & (lam <= 790)
    lam, tr = lam[keep], tr[keep]

    # ---- check 2: the cut-off the paper states in its text ----------------
    plateau = float(np.median(tr[lam > 750]))
    rise = (lam > 320) & (lam < 420)
    order = np.argsort(tr[rise])
    cutoff = float(np.interp(plateau / 2, tr[rise][order], lam[rise][order]))
    print("\ncut-off check (transmittance = 50 % of its 800 nm value)")
    print(f"  plateau    : {100*plateau:.1f} %")
    print(f"  recovered  : {cutoff:.1f} nm")
    print(f"  stated     : {STATED_CUTOFF_NM:.0f} nm      "
          f"difference {cutoff - STATED_CUTOFF_NM:+.1f} nm")
    residual = STATED_CUTOFF_NM - cutoff
    if abs(residual) > CUTOFF_TOLERANCE_NM:
        raise SystemExit(
            f"recovered cut-off is {residual:+.1f} nm from the stated value, "
            f"past the {CUTOFF_TOLERANCE_NM:.0f} nm tolerance. That is the size "
            "of error the rejected column-wise trace produced; do not use this "
            "output.")
    from scint.optical import (CALIBRATION_RESIDUAL_NM, STROKE_PX_VERTICAL,
                               STROKE_PX_HORIZONTAL)
    print(f"  residual   : {residual:+.2f} nm, carried as the calibration term of "
          f"the wavelength-axis systematic")
    print(f"  the three handles about the applied offset: emission "
          f"{STATED_EMISSION_NM - em - offset:+.1f}, excitation "
          f"{STATED_EXCITATION_NM - ex - offset:+.1f}, cut-off {residual:+.1f} nm")
    if abs(abs(residual) - CALIBRATION_RESIDUAL_NM) > 0.1:
        raise SystemExit(
            f"the residual measured here ({abs(residual):.2f} nm) no longer "
            f"matches CALIBRATION_RESIDUAL_NM ({CALIBRATION_RESIDUAL_NM:.2f} nm) "
            "in scint/optical.py, which the paper propagates. Update it.")

    # ---- check 3: the stroke width the error model rests on -------------------
    # Measured with the same mask that traces the curve: vertical thickness on
    # the transparent plateau, horizontal extent on the steep rise. Half of each
    # is the reading half-width scint/optical.py carries.
    vt, ht = [], []
    for x in range(PANEL["left"] + 3, PANEL["right"] - 3):
        if 450 <= raw_nm(x) <= 740:
            yy = np.where(m["green"][PANEL["top"] + 2:PANEL["bottom"] - 1, x])[0]
            if len(yy):
                s, e = longest_run(yy); vt.append(e - s + 1)
    for y in range(PANEL["top"] + 1, PANEL["bottom"]):
        if 15 <= to_pct(y) <= 60:
            xx = np.where(m["green"][y, PANEL["left"] + 3:PANEL["right"] - 3])[0]
            if len(xx):
                s, e = longest_run(xx); ht.append(e - s + 1)
    vt_med, ht_med = int(np.median(vt)), int(np.median(ht))
    print("\nstroke width, measured with the tracing mask")
    print(f"  vertical, on the plateau   : {vt_med} px  (half = {vt_med/2/T_PX_PER_PCT/100:.4f} in T)")
    print(f"  horizontal, on the rise    : {ht_med} px  (half = {ht_med/2/PX_PER_NM:.1f} nm)")
    if (vt_med, ht_med) != (STROKE_PX_VERTICAL, STROKE_PX_HORIZONTAL):
        raise SystemExit(
            f"stroke widths measured here ({vt_med}, {ht_med}) px differ from "
            f"({STROKE_PX_VERTICAL}, {STROKE_PX_HORIZONTAL}) in scint/optical.py, "
            "which the error model rests on. Update it.")

    # Every traced point is written. Where the transmittance sits above the
    # Fresnel limit of a lossless crystal the inversion has no solution and the
    # length column is left empty; dropping those rows, as an earlier version
    # did, silently biased the plateau low by excluding exactly the points
    # that say the crystal is transparent there.
    rows = []
    for lm, t in zip(lam, tr):
        la = attenuation_length(float(lm), float(t), SAMPLE_LENGTH_MM)
        rows.append((float(lm), float(t), la))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fh.write(
            "# NaI(Tl) optical transmittance, digitised from Fig. 2 of\n"
            "#   R. Mao, L. Zhang, R.-Y. Zhu, IEEE Trans. Nucl. Sci. 55 (2008)\n"
            "#   2425-2431, doi:10.1109/TNS.2008.2000776\n"
            "#\n"
            "# Produced by scripts/digitise_mao_fig2.py. The article PDF is not\n"
            "# redistributed here; it was obtained from the authors' page at\n"
            "#   https://www.its.caltech.edu/~rzhu/papers/08_tns_crystal.pdf\n"
            f"#   sha256 {digest}\n"
            "#\n"
            "# The sample is a cylinder 1.5 X0 = 38.8 mm long (Section III).\n"
            "# attenuation_length_mm inverts the paper's own transmittance\n"
            "# expression and is EFFECTIVE: a transmittance measurement cannot\n"
            "# separate absorption from scattering. It is empty where the\n"
            "# transmittance exceeds the Fresnel limit of a lossless crystal\n"
            "# and the inversion has no solution.\n"
            "#\n"
            f"# calibration offset applied: {offset:+.2f} nm, from the emission\n"
            f"#   ({em:.1f} vs {STATED_EMISSION_NM:.0f}) and excitation\n"
            f"#   ({ex:.1f} vs {STATED_EXCITATION_NM:.0f}) peaks printed in the panel\n"
            f"# recovered cut-off {cutoff:.1f} nm vs {STATED_CUTOFF_NM:.0f} nm stated\n"
            "#\n"
        )
        w = csv.writer(fh)
        w.writerow(["wavelength_nm", "transmittance", "attenuation_length_mm"])
        for lm, t, la in rows:
            w.writerow([f"{lm:.2f}", f"{t:.5f}", "" if la is None else f"{la:.3f}"])
    try:
        shown = args.out.relative_to(ROOT)
    except ValueError:
        shown = args.out
    print(f"\nwrote {shown}  ({len(rows)} points, "
          f"{lam.min():.0f}-{lam.max():.0f} nm)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
