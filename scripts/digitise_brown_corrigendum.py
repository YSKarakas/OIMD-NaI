#!/usr/bin/env python3
"""Digitise the NaI(Tl) attenuation length implemented by Brown (2021), from
the corrected figure in the 2023 corrigendum.

Why this exists
---------------
The Geant4 optical-property compilation of Miller et al. (IEEE TNS 72 (2025)
197) takes its NaI(Tl) data from J. M. C. Brown, Appl. Radiat. Isot. 168
(2021) 109368, whose Fig. A.11 shows the attenuation length that simulation
fed to Geant4. A corrigendum (Appl. Radiat. Isot. 194 (2023) 110721) states
that Fig. A11 as printed contained "typographic errors" and replaces it. The
replacement is not a cosmetic change for NaI(Tl): the curve as first printed
sat about a factor 2.8 lower across the visible than the corrected one.

An earlier revision of this study digitised the 2021 figure (from the arXiv
preprint) and quoted a factor 2.5 between that curve and the measurement of
Mao et al. at 415 nm. That number was read from a figure its own author has
withdrawn. This script reads the corrected figure instead; the superseded
digitisation is kept beside it as
data/optical/brown2021_fig_a11_as_printed_nai_attenuation.csv so that the
difference between the two is itself on record.

What it does
------------
1. Extracts the figure from page 1 of the corrigendum at its embedded
   resolution (1457 x 542 px, a JPEG at 300 ppi) -- no re-rendering.
2. Calibrates both axes from the plotted grid lines, which sit at 0.5 eV and
   0.5 decade intervals, by least squares, and reports the residuals. The left
   frame (1.0 eV) is not used in the fit and serves as the check.
3. Traces the NaI(Tl) attenuation stroke (blue, dashed) column-wise where it
   is shallow and row-wise where it is steep, exactly as the Mao digitiser
   does, after masking the flat refractive-index line and the legend.
4. Converts to wavelength and attenuation length: the right-hand axis is
   log10(L / mm) = the left-hand axis value, 0 to 4.

The reading error is the stroke's half-thickness on the log axis, reported in
percent. The corrigendum PDF is not redistributed with this repository; its
checksum is recorded in the output.

    python3 scripts/digitise_brown_corrigendum.py --pdf corrigendum.pdf
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
OUT = ROOT / "data" / "optical" / "brown2021_nai_attenuation.csv"
SUPERSEDED = ROOT / "data" / "optical" / "brown2021_fig_a11_as_printed_nai_attenuation.csv"

# The corrigendum as downloaded from the publisher on 14 Sept 2026.
CORRIGENDUM_SHA256 = "4a00635cb946554af61d1195e5a21f217bc587082dd09a54e4d0a038e5f2a7cd"
CORRIGENDUM_DOI = "10.1016/j.apradiso.2023.110721"

# Figure geometry. Both panels share one embedded image; NaI(Tl) is in the left
# one. The axes are energy 1-5 eV (linear) and 0-4 on the left, which the
# right-hand axis labels as 10^0-10^4 mm.
IMAGE_SIZE = (1457, 542)
E_GRID = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]      # vertical grid lines
Y_GRID = [3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5]            # horizontal grid lines, top down
E_LEFT_FRAME = 1.0
# The legend box occupies the top right of the panel down to about 3.3 on the
# left axis; the attenuation curve never rises above 3.0, so everything above
# this value is excluded from the trace and the trace is checked against it.
LEGEND_FLOOR_Y = 3.25
STEEP_RUN_PX = 8                 # a column run longer than this is traced row-wise
MAX_JUMP_PX = 12                 # continuity: how far the stroke may move per column
HC_EV_NM = 1239.841984

# Calibration tolerance: the grid lines must fit a linear axis to within a
# pixel, and the unused left frame must land on 1.0 eV to within a pixel.
CALIBRATION_TOLERANCE_PX = 1.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_figure(pdf: Path, workdir: Path) -> np.ndarray:
    from PIL import Image

    subprocess.run(["pdfimages", "-png", "-f", "1", "-l", "1", str(pdf),
                    str(workdir / "img")], check=True, capture_output=True)
    for p in sorted(workdir.glob("img*.png")):
        im = Image.open(p)
        if im.size == IMAGE_SIZE:
            return np.asarray(im.convert("RGB")).astype(int)
    raise SystemExit(f"no {IMAGE_SIZE} image on page 1; is this the corrigendum?")


def clusters(idx: np.ndarray) -> list[float]:
    """Mean position of each run of adjacent indices."""
    out, run = [], [idx[0]]
    for v in idx[1:]:
        if v - run[-1] <= 2:
            run.append(v)
        else:
            out.append(float(np.mean(run)))
            run = [v]
    out.append(float(np.mean(run)))
    return out


def calibrate(panel: np.ndarray):
    r, g, b = panel[:, :, 0], panel[:, :, 1], panel[:, :, 2]
    grey = (np.abs(r - g) < 12) & (np.abs(g - b) < 12) & (r > 150) & (r < 225)
    dark = panel.sum(axis=2) < 200
    cols = clusters(np.where(grey.mean(axis=0) > 0.4)[0])
    rows = clusters(np.where(grey.mean(axis=1) > 0.4)[0])
    if len(cols) != len(E_GRID) or len(rows) != len(Y_GRID):
        raise SystemExit(f"expected {len(E_GRID)} vertical and {len(Y_GRID)} "
                         f"horizontal grid lines, found {len(cols)} and {len(rows)}")
    bx, ax = np.polyfit(cols, E_GRID, 1)          # E = ax + bx * x
    by, ay = np.polyfit(rows, Y_GRID, 1)          # y_axis = ay + by * y_px
    res_x = np.array(E_GRID) - (ax + bx * np.array(cols))
    res_y = np.array(Y_GRID) - (ay + by * np.array(rows))
    left_frame = clusters(np.where(dark.mean(axis=0) > 0.5)[0])[0]
    frame_err_px = (E_LEFT_FRAME - (ax + bx * left_frame)) / bx
    return dict(ax=ax, bx=bx, ay=ay, by=by,
                rms_x_px=float(np.sqrt(np.mean((res_x / bx) ** 2))),
                rms_y_px=float(np.sqrt(np.mean((res_y / by) ** 2))),
                left_frame_px=left_frame, frame_err_px=float(frame_err_px),
                x_lo=int(round(cols[0] - (cols[1] - cols[0]))), x_hi=int(round(cols[-1])),
                y_top=int(round(rows[0] - (rows[1] - rows[0]))),
                y_bottom=int(round(rows[-1] + (rows[-1] - rows[-2]))))


def runs_of(idx: np.ndarray) -> list[tuple[int, int]]:
    out, s, p = [], idx[0], idx[0]
    for v in idx[1:]:
        if v - p > 2:
            out.append((s, p))
            s = v
        p = v
    out.append((s, p))
    return out


def trace(panel: np.ndarray, cal: dict):
    r, g, b = panel[:, :, 0], panel[:, :, 1], panel[:, :, 2]
    blue = (b > 180) & (r < 130) & (g < 170)
    to_y = lambda ypx: cal["ay"] + cal["by"] * ypx
    to_e = lambda xpx: cal["ax"] + cal["bx"] * xpx
    x0, x1, y0, y1 = cal["x_lo"] + 2, cal["x_hi"] - 1, cal["y_top"] + 2, cal["y_bottom"] - 1

    # The flat refractive-index line spans the whole panel width: mask it.
    full = blue[:, x0:x1].mean(axis=1) > 0.9
    solid_rows = np.where(full)[0]
    work = blue.copy()
    work[solid_rows, :] = False
    # ...and the legend, which is the only other blue above the plateau.
    legend_row = int((LEGEND_FLOOR_Y - cal["ay"]) / cal["by"])
    work[:legend_row, :] = False

    # The stroke is dashed, and in a dash gap the topmost blue in a column is
    # the dotted emission curve underneath. So the trace follows the stroke by
    # continuity: a column's run counts only if it sits within MAX_JUMP_PX of
    # where the previous column left the curve; a gap column is skipped.
    pts: dict[float, float] = {}
    steep_cols: list[tuple[int, int, int]] = []
    widths: list[int] = []
    prev_y: float | None = None
    for x in range(x0, x1):
        yy = np.where(work[y0:y1, x])[0]
        if not len(yy):
            continue
        runs = runs_of(yy + y0)
        if prev_y is None:
            s, e = runs[0]                      # first column: nothing else above the stroke
        else:
            near = [(min(abs(s - prev_y), abs(e - prev_y)), s, e) for s, e in runs]
            d, s, e = min(near)
            if d > MAX_JUMP_PX:
                continue                        # a dash gap: the stroke is absent here
        if e - s + 1 <= STEEP_RUN_PX:
            pts[to_e(x)] = to_y((s + e) / 2.0)
            widths.append(e - s + 1)
            prev_y = (s + e) / 2.0
        else:
            steep_cols.append((x, s, e))
            prev_y = e                          # the curve only descends: leave at the bottom
    if steep_cols:
        xs = [c[0] for c in steep_cols]
        lo, hi = min(xs) - 2, max(xs) + 3
        for y in range(min(c[1] for c in steep_cols), max(c[2] for c in steep_cols) + 1):
            xx = np.where(work[y, lo:hi])[0]
            if not len(xx):
                continue
            s, e = max(runs_of(xx + lo), key=lambda t: t[1] - t[0])
            if e - s + 1 <= STEEP_RUN_PX:
                pts.setdefault(to_e((s + e) / 2.0), to_y(y))
    stroke_px = float(np.median(widths))
    return sorted(pts.items()), stroke_px, solid_rows


def interp_log(curve, lam_nm: float) -> float:
    lam = np.array([c[0] for c in curve]); lg = np.log10([c[1] for c in curve])
    return float(10 ** np.interp(lam_nm, lam, lg))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdf", type=Path, required=True,
                    help="J. M. C. Brown, Appl. Radiat. Isot. 194 (2023) 110721")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    digest = sha256(args.pdf)
    if digest != CORRIGENDUM_SHA256:
        print(f"WARNING: PDF sha256 {digest[:16]} differs from the recorded "
              f"{CORRIGENDUM_SHA256[:16]}; the output records what was read")

    with tempfile.TemporaryDirectory() as tmp:
        img = extract_figure(args.pdf, Path(tmp))
    panel = img[:, : img.shape[1] // 2]
    cal = calibrate(panel)
    print(f"axis calibration from grid lines: rms residual {cal['rms_x_px']:.2f} px "
          f"(energy), {cal['rms_y_px']:.2f} px (log length); "
          f"left frame predicts 1.0 eV to {cal['frame_err_px']:+.2f} px")
    ok = max(cal["rms_x_px"], cal["rms_y_px"], abs(cal["frame_err_px"])) < CALIBRATION_TOLERANCE_PX
    print(f"  calibration check: {'PASS' if ok else 'FAIL'}")
    if not ok:
        return 1

    curve_e, stroke_px, solid_rows = trace(panel, cal)
    half_decades = 0.5 * stroke_px * abs(cal["by"])
    n_line = cal["ay"] + cal["by"] * float(np.mean(solid_rows))
    print(f"stroke width {stroke_px:.0f} px -> half-thickness {half_decades:.4f} decades "
          f"= {100 * (10 ** half_decades - 1):.1f} % in attenuation length")
    print(f"flat refractive-index line reads n = {n_line:.3f} "
          f"(the figure states no number; the value the earlier revision quoted was 1.85)")
    top = max(v for _, v in curve_e)
    print(f"highest traced point {top:.3f} on the left axis (legend floor {LEGEND_FLOOR_Y})")
    if top > LEGEND_FLOOR_Y - 0.1:
        print("  the trace reaches the legend mask: FAIL"); return 1

    rows = sorted(((HC_EV_NM / e, e, 10 ** v) for e, v in curve_e), key=lambda t: t[0])
    curve = [(lam, L) for lam, _, L in rows]
    l415 = interp_log(curve, 415.0)
    print(f"{len(rows)} points, {rows[0][0]:.0f}-{rows[-1][0]:.0f} nm; "
          f"attenuation length at 415 nm: {l415:.0f} mm")
    if SUPERSEDED.exists():
        with SUPERSEDED.open() as fh:
            old = [(float(r["wavelength_nm"]), float(r["attenuation_length_mm"]))
                   for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))]
        o415 = interp_log(old, 415.0)
        print(f"figure as first printed (superseded file): {o415:.0f} mm at 415 nm; "
              f"corrected / as printed = {l415 / o415:.2f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fh.write(
            "# NaI(Tl) attenuation length AS IMPLEMENTED IN A PUBLISHED GEANT4 MODEL,\n"
            "# digitised from the CORRECTED Fig. 1 (left panel, blue dashed) of the corrigendum\n"
            "#   J. M. C. Brown, Appl. Radiat. Isot. 194 (2023) 110721,\n"
            f"#   doi:{CORRIGENDUM_DOI}\n"
            f"#   PDF sha256 {digest}\n"
            "# which replaces Fig. A11 of\n"
            "#   J. M. C. Brown, Appl. Radiat. Isot. 168 (2021) 109368,\n"
            "#   doi:10.1016/j.apradiso.2020.109368; preprint arXiv:1908.04565.\n"
            "# The figure as first printed is digitised in\n"
            f"#   {SUPERSEDED.name}\n"
            "# and is kept only to document the size of the correction.\n"
            "#\n"
            "# This is NOT a measurement. It is the model one published simulation fed to\n"
            "# Geant4, recorded here so that it can be compared with the measurement in\n"
            "# mao2008_nai_transmittance.csv. It matters because the parameter compilation\n"
            "# that supplies NaI(Tl) to other Geant4 users (Miller et al., IEEE TNS 72\n"
            "# (2025) 197) takes its NaI(Tl) optical data from this work. Brown's Table A3\n"
            "# (Table 1 of the corrigendum) attributes the NaI(Tl) row, absorption length\n"
            "# included, to Mao et al. (2008); how the curve was derived from that\n"
            "# transmittance measurement is not stated.\n"
            "#\n"
            f"# Axes calibrated from the grid lines (rms residual {cal['rms_x_px']:.2f} px in\n"
            f"# energy, {cal['rms_y_px']:.2f} px in log length); the left frame lands on 1.0 eV\n"
            f"# to {cal['frame_err_px']:+.2f} px. Stroke {stroke_px:.0f} px on a log axis: half-\n"
            f"# thickness {100 * (10 ** half_decades - 1):.1f} %. Steep parts traced row-wise.\n"
            "# Generated by scripts/digitise_brown_corrigendum.py.\n"
            "wavelength_nm,energy_eV,attenuation_length_mm\n")
        for lam, e, L in rows:
            fh.write(f"{lam:.2f},{e:.4f},{L:.1f}\n")
    print(f"wrote {args.out.relative_to(ROOT)} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
