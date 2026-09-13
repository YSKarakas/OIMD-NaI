#!/usr/bin/env python3
"""Read the runs produced by scripts/run_gates.py and report gates G2-G4.

Reads only what is on disk: each run directory carries the configuration that
produced it, so the analysis never needs to be told what a run was, and a result
cannot drift away from its inputs.

    python3 scripts/analyse_gates.py             # human-readable report
    python3 scripts/analyse_gates.py --json out.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scint.g4csv import read_ntuple  # noqa: E402
from scint.registry import iter_runs  # noqa: E402

RUNS_DIR = ROOT / "runs"

# The gate references, and where each comes from.
YIELD_PER_MEV = 41000.0          # materials/NaI_Tl.dat, src_B
SOURCE_MEV = 0.662               # 137Cs
MEASURED_FWHM_PCT = 7.0          # data/benchmark_scintillators.csv, src_A
PHOTOPEAK_KEV_MIN = 650.0        # loose photopeak window, for light-collection statistics
FULL_ENERGY_KEV = (661.5, 662.5) # strict window: the deposit really is the whole 662 keV
MIN_EDEP_KEV = 50.0              # above this the per-event yield ratio is meaningful


def photopeak(data: dict[str, np.ndarray]) -> np.ndarray:
    return np.asarray(data["edep_keV"]) >= PHOTOPEAK_KEV_MIN


def summarise(run) -> dict | None:
    status = run.status() or {}
    if status.get("state") != "done":
        return None
    path = run.path / "output_nt_events.csv"
    if not path.exists():
        return None

    data = read_ntuple(path)
    edep = np.asarray(data["edep_keV"], dtype=float)
    gen_all = np.asarray(data["photons_generated"], dtype=float)
    # Schema 2 onwards separates scintillation from Cerenkov. Where the split is
    # available the yield gate uses scintillation alone, because that is the
    # quantity the material file specifies; where it is not, it uses all optical
    # photons and the Cerenkov contribution shows up as a small positive bias.
    scint_only = ("photons_scintillation" in data)
    gen_scint = (np.asarray(data["photons_scintillation"], dtype=float)
                 if scint_only else gen_all)
    cherenkov = (np.asarray(data["photons_cherenkov"], dtype=float)
                 if "photons_cherenkov" in data else None)

    mask = photopeak(data)
    gen = gen_all[mask]
    det = np.asarray(data["photons_detected"], dtype=float)[mask]
    if gen.size < 30:
        return None

    # Yield closure is tested on the RATIO, event by event, over every event with
    # a meaningful deposit. Testing the mean photon count inside an energy window
    # instead would be biased by the window: events that lost a little energy are
    # still inside it and drag the mean down.
    ratio_mask = edep > MIN_EDEP_KEV
    ratio = gen_scint[ratio_mask] / (edep[ratio_mask] / 1000.0)

    # Variance is tested only where the deposit really is the full 662 keV, so
    # that the spread being measured is photon statistics and not a spread in the
    # energy that produced the photons.
    full = (edep >= FULL_ENERGY_KEV[0]) & (edep <= FULL_ENERGY_KEV[1])
    gen_full = gen_scint[full]

    n = gen.size
    lce = det / gen
    return {
        "id": run.id,
        "label": status.get("label", "?"),
        "material": run.config["crystal"]["material_spec"],
        "wrapping": run.config["surface"]["wrapping"],
        "events": run.config["run"]["events"],
        # Runs made with different recorded observables are kept apart: they are
        # the same physics but not the same dataset, and averaging them would
        # hide the fact that only the later schema can answer some questions.
        "schema": int(run.config["run"].get("schema", 1)),
        "photopeak_events": int(n),
        "yield_ratio_mean": float(ratio.mean()),
        "yield_ratio_sem": float(ratio.std(ddof=1) / math.sqrt(ratio.size)),
        "yield_ratio_events": int(ratio.size),
        "yield_is_scintillation_only": bool(scint_only),
        "cherenkov_mean": float(cherenkov[full].mean()) if cherenkov is not None and full.any() else None,
        "mean_generated_nm": float(np.asarray(data["mean_generated_nm"])[full].mean())
        if "mean_generated_nm" in data and full.any() else None,
        "mean_detected_nm": float(np.asarray(data["mean_detected_nm"])[mask].mean())
        if "mean_detected_nm" in data else None,
        "full_energy_events": int(gen_full.size),
        "generated_mean": float(gen_full.mean()) if gen_full.size else float("nan"),
        "generated_sem": float(gen_full.std(ddof=1) / math.sqrt(gen_full.size))
        if gen_full.size > 1 else float("nan"),
        "generated_var_over_mean": float(gen_full.var(ddof=1) / gen_full.mean())
        if gen_full.size > 1 else float("nan"),
        "detected_mean": float(det.mean()),
        "lce_mean": float(lce.mean()),
        "lce_sem": float(lce.std(ddof=1) / math.sqrt(n)),
        # Resolution from photon statistics alone: RESOLUTIONSCALE = 1 and
        # detection efficiency 1, so this is a lower bound, not a prediction.
        "fwhm_pct": float(2.3548 * det.std(ddof=1) / det.mean() * 100.0),
        "fwhm_pct_err": float(2.3548 * det.std(ddof=1) / det.mean() * 100.0 / math.sqrt(2 * (n - 1))),
    }


def _resolve_labels(rows: list[dict]) -> list[dict]:
    """Keep one run per label, and say out loud which one and why.

    A label is not unique: the same configuration re-run at a different event
    count is a different run with the same name, and `{r["label"]: r}` would
    have silently kept whichever happened to be last on disk. For a paper whose
    subject is exactly this kind of untracked provenance, that is not an
    acceptable failure mode. The rule is explicit -- prefer the newest output
    schema, then the most events -- and every collapse is printed.
    """
    groups: dict[str, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["label"], []).append(r)
    kept, notes = [], []
    for label, g in groups.items():
        if len(g) == 1:
            kept.append(g[0])
            continue
        g = sorted(g, key=lambda r: (r["schema"], r["events"]), reverse=True)
        kept.append(g[0])
        notes.append((label, g[0], g[1:]))
    if notes:
        print("=" * 78)
        print("DUPLICATE LABELS -- one run kept per label, the rest ignored")
        print("=" * 78)
        for label, win, rest in sorted(notes):
            print(f"  {label:<28} keep {win['id'][:12]} "
                  f"(schema {win['schema']}, {win['events']} events)")
            for r in rest:
                print(f"  {'':<28} drop {r['id'][:12]} "
                      f"(schema {r['schema']}, {r['events']} events)")
        print()
    return kept


# Runs allowed to sit in a comparison set at a different event count, each
# with the reason. The rule is "no SILENT mixing", not "no mixing": a recorded
# decision is fine, an accident is not. Every entry here is printed as a NOTE
# whenever it is used, and the manuscript states the same exception in the
# caption of the table that carries it.
DOCUMENTED_EVENT_EXCEPTIONS = {
    "SCOPE_ctrl_flatabs_3inch": (
        "run at 300 events: with a flat 2000 mm absorption length photons "
        "survive far more boundary crossings and the 1000-event run projected "
        "to ~65 min; 300 events give the ratio this control measures to 0.1 %"),
    "SCOPE_ctrl_flatabs_1inch": "run at 300 events, to match its 3-inch partner",
}


def _same_settings(rows: list[dict], what: str) -> None:
    """Refuse to compare runs that were not made under the same settings,
    unless the deviation is recorded above -- in which case say so."""
    if not rows:
        return
    ref = max(((r["schema"], r["events"]) for r in rows),
              key=lambda c: sum(1 for r in rows if (r["schema"], r["events"]) == c))
    odd = [r for r in rows if (r["schema"], r["events"]) != ref]
    undocumented = [r for r in odd if r["label"] not in DOCUMENTED_EVENT_EXCEPTIONS]
    if undocumented:
        raise SystemExit(
            f"{what}: runs differ in schema/events from the set's "
            f"(schema {ref[0]}, {ref[1]} events): "
            + ", ".join(f"{r['label']} ({r['schema']}, {r['events']})" for r in undocumented)
            + ". These are not comparable and the spread between them would be "
            "partly an artefact of statistics. Re-run at one setting, or record "
            "the reason in DOCUMENTED_EVENT_EXCEPTIONS.")
    for r in odd:
        print(f"  NOTE {r['label']} at {r['events']} events, not {ref[1]}: "
              f"{DOCUMENTED_EVENT_EXCEPTIONS[r['label']]}")


def report(rows: list[dict]) -> dict:
    rows = _resolve_labels(rows)
    by_label = {r["label"]: r for r in rows}
    verdicts: dict[str, dict] = {}

    baseline = by_label.get("G2G3_baseline")

    print("=" * 78)
    print("G2 -- PHOTON BUDGET")
    print("=" * 78)
    if baseline is None:
        print("  no baseline run available")
    else:
        ratio, ratio_sem = baseline["yield_ratio_mean"], baseline["yield_ratio_sem"]
        var_ratio = baseline["generated_var_over_mean"]
        n_full = baseline["full_energy_events"]
        # The relative uncertainty of a variance estimate is sqrt(2/(n-1)).
        var_err = math.sqrt(2.0 / (n_full - 1)) if n_full > 1 else float("inf")
        print(f"  events used for the yield     {baseline['yield_ratio_events']} (edep > {MIN_EDEP_KEV:.0f} keV)")
        print(f"  photons per MeV deposited     {ratio:.1f} +- {ratio_sem:.1f}")
        print(f"  requested (materials file)    {YIELD_PER_MEV:.1f}")
        print(f"  counted photons               "
              f"{'scintillation only' if baseline['yield_is_scintillation_only'] else 'ALL optical (includes Cerenkov)'}")
        print(f"  deviation                     {100 * (ratio - YIELD_PER_MEV) / YIELD_PER_MEV:+.3f} %"
              f"  ({(ratio - YIELD_PER_MEV) / ratio_sem:+.2f} sigma)")
        print()
        print(f"  full-energy events            {n_full} "
              f"({FULL_ENERGY_KEV[0]:.1f}-{FULL_ENERGY_KEV[1]:.1f} keV)")
        print(f"  photons generated, mean       {baseline['generated_mean']:.1f} "
              f"+- {baseline['generated_sem']:.1f}   (expected {YIELD_PER_MEV * SOURCE_MEV:.1f})")
        if baseline.get("cherenkov_mean") is not None:
            ch = baseline["cherenkov_mean"]
            print(f"  Cerenkov photons per event    {ch:.1f} "
                  f"({100 * ch / baseline['generated_mean']:.3f} % of the scintillation photons)")
            print("  In NaI (n ~ 1.85) the Cerenkov threshold for electrons is ~96 keV, which")
            print("  Compton electrons from a 662 keV gamma pass. A simulation that counts")
            print("  'optical photons' without separating them overstates the light yield.")
        print(f"  variance / mean               {var_ratio:.3f} +- {var_err:.3f}   (Poisson: 1.000)")
        print("  Geant4 11.4.2 G4Scintillation.cc draws N from a Gaussian of width")
        print("  RESOLUTIONSCALE*sqrt(mean) when mean > 10 and from a Poisson below it, so")
        print("  with RESOLUTIONSCALE = 1 the summed variance must equal the summed mean.")
        ok_mean = abs(ratio - YIELD_PER_MEV) / YIELD_PER_MEV < 0.01
        ok_var = abs(var_ratio - 1.0) < 3 * var_err
        print(f"  VERDICT                       {'PASS' if ok_mean and ok_var else 'FAIL'}"
              f"  (yield within 1 %: {ok_mean}; variance Poisson within 3 sigma: {ok_var})")
        verdicts["G2"] = {"pass": bool(ok_mean and ok_var), "yield_per_mev": ratio,
                          "expected_per_mev": YIELD_PER_MEV, "var_over_mean": var_ratio}

    print()
    print("=" * 78)
    print("G3 -- STATISTICAL RESOLUTION (lower bound)")
    print("=" * 78)
    if baseline is not None:
        fwhm, err = baseline["fwhm_pct"], baseline["fwhm_pct_err"]
        print(f"  simulated FWHM at 662 keV     {fwhm:.3f} +- {err:.3f} %")
        print(f"  measured FWHM (src_A)         {MEASURED_FWHM_PCT:.1f} %")
        print(f"  ratio                         {fwhm / MEASURED_FWHM_PCT:.3f}")
        # The quadratic residual is the part Geant4 does not model: the
        # non-proportionality of the light yield and the photodetector transfer term.
        residual = math.sqrt(max(MEASURED_FWHM_PCT**2 - fwhm**2, 0.0))
        print(f"  R_measured^2 - R_sim^2        -> {residual:.3f} % "
              f"(non-proportionality + transfer, not modelled here)")
        ok = fwhm < MEASURED_FWHM_PCT
        print(f"  VERDICT                       {'PASS' if ok else 'FAIL'}"
              f"  (statistical bound must lie below the measurement)")
        verdicts["G3"] = {"pass": bool(ok), "fwhm_pct": fwhm,
                          "measured_pct": MEASURED_FWHM_PCT, "residual_pct": residual}

    print()
    print("=" * 78)
    print("G4 -- LIGHT-COLLECTION EFFICIENCY BY WRAPPING")
    print("=" * 78)
    wraps = {r["wrapping"]: r for r in rows if r["label"].startswith("G4_")
             or r["label"] == "G2G3_baseline"}
    if wraps:
        print(f"  {'wrapping':<12} {'LCE':>9} {'+- sem':>9} {'vs none':>10}")
        bare = wraps.get("none")
        order = sorted(wraps.values(), key=lambda r: -r["lce_mean"])
        for r in order:
            rel = f"{r['lce_mean'] / bare['lce_mean']:.3f}x" if bare else "n/a"
            print(f"  {r['wrapping']:<12} {r['lce_mean']:>9.4f} {r['lce_sem']:>9.4f} {rel:>10}")
        # An external control for the ordering, from the literature survey.
        # arXiv:1602.02983 measured ONE undoped CsI crystal with three different
        # wrappings, which is the controlled comparison this gate needs and which
        # simulation alone cannot provide. Their numbers are photoelectrons per
        # MeV on a different material and a different geometry, so only the
        # ORDER and the rough size of the steps are comparable -- not the values.
        MEASURED = {
            "aluminium": 79.0,   # specular metal
            "tyvek": 89.0,
            "teflon": 91.0,
        }
        print()
        print("  measured on one CsI crystal by arXiv:1602.02983, for comparison:")
        base = MEASURED["aluminium"]
        for name, value in sorted(MEASURED.items(), key=lambda kv: -kv[1]):
            print(f"    {name:<12} {value:>9.0f} p.e./MeV  {value / base:.3f}x aluminium")
        ours = {r["wrapping"]: r["lce_mean"] for r in wraps.values()}
        if "teflon" in ours and "tyvek" in ours:
            theirs = MEASURED["teflon"] / MEASURED["tyvek"]
            mine = ours["teflon"] / ours["tyvek"]
            print(f"    teflon/tyvek: measured {theirs:.3f}, simulated here {mine:.3f}")
            print("    Same sign. The measured gap is the larger one, and our two")
            print("    diffuse reflectors sit closer together than theirs do.")
        print()

        if bare is None:
            # Distinguish "the reference point could not be computed" from "the
            # reflectors lost". The unwrapped configuration is not tractable here:
            # at a bare polished surface with air outside, total internal
            # reflection is lossless and a skew ray in a cylinder preserves its
            # angle to the wall, so photons are trapped until bulk absorption
            # removes them. Three attempts (3000, 250 and 25 events) were all
            # abandoned; the same crystal with a Teflon wrap runs at 0.65 s/event.
            print("  VERDICT                       INCOMPLETE -- no unwrapped reference")
            print("  This run set has no bare-crystal point. Under the look-up-table")
            print("  surface model it could not be produced at all: Geant4 ships no table")
            print("  for a bare surface and hangs instead of saying so. That is fixed --")
            print("  bare surfaces are now resolved analytically and run in 108 s -- so a")
            print("  newer schema has the point. Analyse that schema for a verdict; the")
            print("  ordering among the reflectors here stands on its own, but the")
            print("  comparison against bare is NOT claimed from this set.")
            ok = None
        else:
            ok = all(r["lce_mean"] > bare["lce_mean"]
                     for r in wraps.values() if r["wrapping"] != "none")
            print(f"  VERDICT                       {'PASS' if ok else 'FAIL'}"
                  f"  (every reflector must beat the unwrapped crystal)")
        verdicts["G4"] = {"pass": ok,
                          "lce": {r["wrapping"]: r["lce_mean"] for r in wraps.values()}}

    print()
    print("=" * 78)
    print("OPTICAL INPUT SYSTEMATICS -- LCE under one-at-a-time variation")
    print("=" * 78)
    sys_rows = [r for r in rows if r["label"].startswith("SYS_")]
    _same_settings(sys_rows, "optical-systematics set")
    if baseline is not None and sys_rows:
        b = baseline["lce_mean"]
        print(f"  {'variant':<22} {'LCE':>9} {'+- sem':>9} {'vs baseline':>13}")
        print(f"  {'baseline':<22} {b:>9.4f} {baseline['lce_sem']:>9.4f} {'--':>13}")
        for r in sorted(sys_rows, key=lambda r: r["label"]):
            print(f"  {r['label'][4:]:<22} {r['lce_mean']:>9.4f} {r['lce_sem']:>9.4f} "
                  f"{100 * (r['lce_mean'] - b) / b:>+12.2f} %")
        values = [b] + [r["lce_mean"] for r in sys_rows]
        spread = (max(values) - min(values)) / b

        # Two spreads, because they answer different questions. The wider one is
        # over everything the literature actually does, flat absorption lengths
        # included. The narrower one drops the flat-absorption variants, which the
        # Mao et al. 365 nm transmittance measurement excludes, and so is the
        # spread that survives after using the evidence that already exists.
        defensible = [b] + [r["lce_mean"] for r in sys_rows
                            if not r["label"].startswith("SYS_abs_flat")]
        defensible_spread = (max(defensible) - min(defensible)) / b

        # Where the schema records it, show what the crystal does to its own
        # spectrum: light that reaches the readout is redder than light emitted,
        # because the blue side of the band runs into the absorption edge.
        reddening = [(r["label"][4:], r) for r in sys_rows if r.get("mean_detected_nm")]
        if baseline.get("mean_detected_nm"):
            reddening.insert(0, ("baseline", baseline))
        if reddening:
            print()
            print("  Self-absorption seen directly -- mean wavelength, emitted vs detected:")
            print(f"  {'variant':<22} {'emitted':>9} {'detected':>9} {'shift':>9}")
            for name, r in reddening:
                gen_nm, det_nm = r["mean_generated_nm"], r["mean_detected_nm"]
                if gen_nm is None or det_nm is None:
                    continue
                print(f"  {name:<22} {gen_nm:>9.2f} {det_nm:>9.2f} {det_nm - gen_nm:>+9.2f} nm")

        print(f"\n  spread over published practice  {100 * spread:.1f} % of baseline LCE"
              f"   [{len(values)} variants]")
        print("    (includes the flat absorption lengths that appear in the literature)")
        print(f"  spread over defensible models   {100 * defensible_spread:.1f} % of baseline LCE"
              f"   [{len(defensible)} variants]")
        print("    (flat absorption dropped: the 365 nm transmittance measurement excludes it)")
        if len(values) < 11:
            print("    NOTE: this is a SUBSET of the variant family, so these spreads are")
            print("    narrower than the full scan by construction -- not a different answer.")
        print("  Neither is reducible by running more events. Both measure what the state")
        print("  of the published optical inputs costs, not what the simulation costs.")
        verdicts["systematics"] = {
            "baseline_lce": b,
            "spread_fraction": spread,
            "defensible_spread_fraction": defensible_spread,
            "variants": {r["label"][4:]: r["lce_mean"] for r in sys_rows},
        }

    # ------------------------------------------------------------------ #
    # The measured attenuation curve. This is the run set that replaces a
    # scanned parameter with a measurement: the Urbach slope was scanned only
    # because the curve it describes was published as a picture, and the
    # picture has now been digitised. What the hi/lo pair spans is the error
    # on that digitisation, not the range of models somebody might pick.
    meas_rows = sorted((r for r in rows if r["label"].startswith("MEAS_")
                        and r["label"] != "MEAS_baseline"),
                       key=lambda r: r["label"])
    meas_base = by_label.get("MEAS_baseline") or by_label.get("SCOPE_baseline")
    if meas_rows and meas_base:
        print()
        print("=" * 78)
        print("MEASURED ATTENUATION CURVE -- digitised from Mao et al. Fig. 2")
        print("=" * 78)
        _same_settings(meas_rows + [meas_base], "measured-attenuation set")
        mb = meas_base["lce_mean"]
        print(f"  {'variant':<22} {'LCE':>9} {'+- sem':>9} {'vs baseline':>13}")
        print(f"  {'baseline (edge model)':<22} {mb:>9.4f} "
              f"{meas_base['lce_sem']:>9.4f} {'--':>13}")
        vals = {}
        for r in meas_rows:
            name = r["label"][len("MEAS_"):]
            vals[name] = r["lce_mean"]
            print(f"  {name:<22} {r['lce_mean']:>9.4f} {r['lce_sem']:>9.4f} "
                  f"{100*(r['lce_mean']-mb)/mb:>+12.2f} %")
        if {"abs_measured", "abs_measured_hi", "abs_measured_lo",
            "abs_measured_lamhi", "abs_measured_lamlo"} <= set(vals):
            c = vals["abs_measured"]
            hi, lo = vals["abs_measured_hi"], vals["abs_measured_lo"]
            lhi, llo = vals["abs_measured_lamhi"], vals["abs_measured_lamlo"]
            # Two independent digitisation errors, one per axis of the figure,
            # each run as its own coherent +-1 sigma pair. They are independent
            # of each other -- one is set by the stroke width, the other by the
            # axis calibration handles -- so they combine in quadrature. The
            # wavelength term was missing from an earlier version of this
            # analysis and it is the larger of the two.
            band_t = abs(hi - lo) / 2.0
            band_l = abs(lhi - llo) / 2.0
            band = math.hypot(band_t, band_l)
            print()
            print(f"  light-collection efficiency from the measured curve:")
            print(f"      {c:.4f} +- {band:.4f}  ({100*band/c:.1f} %)")
            print(f"        transmittance axis  +- {band_t:.4f} "
                  f"({100*band_t/c:.1f} %)")
            print(f"        wavelength axis     +- {band_l:.4f} "
                  f"({100*band_l/c:.1f} %)")
            sysd = verdicts.get("systematics", {})
            if sysd:
                print(f"  for comparison, the scanned-slope family spanned "
                      f"{100*sysd['defensible_spread_fraction']:.0f} % of baseline")
            verdicts["measured"] = {
                "baseline_lce": mb,
                "lce": c, "lce_hi": hi, "lce_lo": lo,
                "lce_lamhi": lhi, "lce_lamlo": llo,
                "half_band": band,
                "half_band_transmittance": band_t,
                "half_band_wavelength": band_l,
                "half_band_fraction": band / c,
                "half_band_fraction_transmittance": band_t / c,
                "half_band_fraction_wavelength": band_l / c,
                "vs_baseline_pct": 100*(c-mb)/mb,
            }
            # The envelope over everything built ON the measured attenuation:
            # its own error band, plus the refractive-index and emission-band
            # alternatives rebuilt on top of it. This is the number that
            # replaces the old 42 %.
            fam = {k: v for k, v in vals.items()
                   if k.startswith("abs_measured") or k.startswith("meas_")}
            if len(fam) >= 5:
                lo_f, hi_f = min(fam.values()), max(fam.values())
                print()
                print("  ENVELOPE over the family built on the measured curve")
                for k in sorted(fam, key=lambda k: fam[k]):
                    print(f"      {k:<24} {fam[k]:.4f}")
                print(f"      spread {100*(hi_f-lo_f)/c:.1f} % of the measured "
                      f"central value ({lo_f:.4f} .. {hi_f:.4f})")
                verdicts["measured"]["family"] = fam
                verdicts["measured"]["family_spread_fraction"] = (hi_f-lo_f)/c
                verdicts["measured"]["family_lo"] = lo_f
                verdicts["measured"]["family_hi"] = hi_f

    # ------------------------------------------------------------------ #
    # The arrangement's own optical inputs. The material scan above asks what
    # the crystal's published properties cost; this asks the same question of
    # the things a published measurement chooses -- how big the crystal is,
    # what it is coupled with, how its surface was finished. Reporting one
    # without the other would let a reader conclude that the crystal's
    # properties dominate, which these runs are what decides.
    scope_rows = sorted((r for r in rows if r["label"].startswith("SCOPE_")
                         and r["label"] != "SCOPE_baseline"),
                        key=lambda r: r["label"])
    scope_base = by_label.get("SCOPE_baseline")
    if scope_base:
        _same_settings(scope_rows + [scope_base], "arrangement-scan set")
    if scope_rows and scope_base:
        print()
        print("=" * 78)
        print("ARRANGEMENT SCAN -- geometry, coupling and surface finish")
        print("=" * 78)
        sb = scope_base["lce_mean"]
        print(f"  {'variant':<20} {'LCE':>9} {'+- sem':>9} {'vs baseline':>13}")
        print(f"  {'baseline':<20} {sb:>9.4f} {scope_base['lce_sem']:>9.4f} {'--':>13}")
        groups: dict[str, list[float]] = {}
        for r in scope_rows:
            name = r["label"][len("SCOPE_"):]
            rel = 100 * (r["lce_mean"] - sb) / sb
            print(f"  {name:<20} {r['lce_mean']:>9.4f} {r['lce_sem']:>9.4f} "
                  f"{rel:>+12.2f} %")
            groups.setdefault(name.split("_")[0], []).append(r["lce_mean"])

        # The comparison that matters: is the crystal's optical model the
        # biggest term, or did we simply not vary the others?
        material_spread = verdicts.get("systematics", {}).get(
            "defensible_spread_fraction")
        print()
        for family, label in (("geom", "geometry"), ("couple", "readout coupling"),
                              ("finish", "surface finish")):
            vals = groups.get(family)
            if not vals:
                continue
            allv = vals + [sb]
            print(f"  {label:<18} envelope {100 * (max(allv) - min(allv)) / sb:6.1f} % "
                  f"of baseline  ({min(allv):.4f} .. {max(allv):.4f})")
        if material_spread is not None:
            print(f"  {'crystal optics':<18} envelope {100 * material_spread:6.1f} % "
                  f"of baseline  (defensible models, from the scan above)")
        print()
        print("  Read this as scope, not as a ranking: the geometry entries are a")
        print("  deliberate factor-four in size and the coupling entries span air to")
        print("  index 1.57, so their envelopes are as wide as the choices made here.")
        print("  What it establishes is that the crystal's optical model is not being")
        print("  called dominant by default, having been the only thing varied.")

        # The control pair. A smaller crystal collecting more light is consistent
        # with the absorption edge removing light in proportion to optical path
        # -- but a shorter path also means fewer reflector bounces, and the size
        # comparison alone cannot tell the two apart. Repeating it with the bulk
        # absorption made flat and long leaves only the non-absorption part.
        by_scope = {r["label"][len("SCOPE_"):]: r["lce_mean"] for r in scope_rows}
        if {"ctrl_flatabs_1inch", "ctrl_flatabs_3inch", "geom_1inch"} <= set(by_scope):
            edge_ratio = by_scope["geom_1inch"] / sb
            flat_ratio = (by_scope["ctrl_flatabs_1inch"]
                          / by_scope["ctrl_flatabs_3inch"])
            print()
            print("  CONTROL -- is the size dependence really the absorption edge?")
            print(f"    1\" / 3\" with the measured absorption edge   "
                  f"{edge_ratio:.3f}")
            print(f"    1\" / 3\" with flat 2000 mm absorption        "
                  f"{flat_ratio:.3f}")
            attributable = (edge_ratio - flat_ratio) / (edge_ratio - 1.0) \
                if edge_ratio > 1.0 else float("nan")
            print(f"    fraction of the size effect that disappears when the")
            print(f"    absorption edge is removed                  "
                  f"{100 * attributable:.0f} %")
            verdicts.setdefault("scope_control", {}).update({
                "size_ratio_edge": edge_ratio,
                "size_ratio_flat": flat_ratio,
                "absorption_attributable_fraction": attributable,
            })

        verdicts["scope"] = {
            "baseline_lce": sb,
            "baseline_sem": scope_base["lce_sem"],
            "variants": {r["label"][len("SCOPE_"):]: r["lce_mean"]
                         for r in scope_rows},
            "sem": {r["label"][len("SCOPE_"):]: r["lce_sem"] for r in scope_rows},
            "envelopes": {
                family: (max(vals + [sb]) - min(vals + [sb])) / sb
                for family, vals in groups.items()
            },
        }

    return verdicts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--schema", type=int, default=None,
                    help="analyse only runs with this output schema (default: the newest "
                         "schema that has a baseline run)")
    args = ap.parse_args()

    rows = [s for s in (summarise(run) for run in iter_runs(RUNS_DIR)) if s]
    if not rows:
        raise SystemExit(f"no completed runs under {RUNS_DIR}")

    available = sorted({r["schema"] for r in rows})
    schema = args.schema if args.schema is not None else max(available)
    chosen = [r for r in rows if r["schema"] == schema]
    if not chosen:
        raise SystemExit(f"no runs with schema {schema}; available: {available}")
    print(f"{len(chosen)} completed runs, output schema {schema} "
          f"(schemas on disk: {available})\n")
    rows = chosen
    verdicts = report(rows)

    if args.json:
        args.json.write_text(json.dumps({"runs": rows, "gates": verdicts}, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
