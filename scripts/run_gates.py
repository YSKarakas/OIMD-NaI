#!/usr/bin/env python3
"""Run the optical validation gates G2-G4 and the optical-input systematics scan.

Each configuration becomes a registry run: the directory name is the hash of the
configuration, so a result can never be separated from the inputs that produced
it, and re-running an identical configuration is refused rather than silently
overwriting the earlier result.

    python3 scripts/run_gates.py --set gates        # G2/G3 baseline + G4 wrappings
    python3 scripts/run_gates.py --set systematics  # the optical-input family
    python3 scripts/run_gates.py --set scope        # geometry, coupling, finish
    python3 scripts/run_gates.py --set all

Gates, and what each is actually testing:

G2  PHOTON BUDGET.  Does Geant4 generate the number of photons the material file
    asks for, with the statistics it asks for?  For full-energy 662 keV events the
    expected mean is 41000 ph/MeV * 0.662 MeV = 27142, and with RESOLUTIONSCALE = 1
    the variance should equal the mean.  No free parameters; a failure here means
    the generation side is wrong and nothing downstream can be trusted.

G3  STATISTICAL RESOLUTION.  The photopeak width from photon statistics alone, with
    RESOLUTIONSCALE = 1 and detection efficiency 1.  This is a LOWER BOUND on the
    energy resolution, not a prediction of it: the gap to the measured 7.0 % FWHM
    is the non-proportionality and transfer term that Geant4 does not model.  The
    gate is that the simulated width must come out BELOW the measurement.  If it
    came out above, the model would be unphysical.

G4  LIGHT-COLLECTION ORDERING.  LCE across the supported wrappings.  Absolute LCE
    depends on the optical inputs being scanned, but the ORDERING should be robust:
    a specular reflector (ESR) and diffuse reflectors (Teflon, Lumirror, Tyvek,
    TiO2) must all collect more light than an unwrapped crystal.

    What "unwrapped" means here needs saying, because it is not obviously the
    worst case. The crystal is POLISHED and the world is air, so a photon outside
    the critical angle is totally internally reflected with no loss at all, and in
    a cylinder a skew ray preserves its angle to the wall -- such a photon is
    trapped until the bulk absorbs it or it reaches the readout face, where the
    grease coupling opens the critical angle and lets it out. A bare polished
    crystal is therefore a light pipe, not a sieve. Whether that beats a diffuse
    wrap is a question for the simulation, not an assumption; the gate states the
    conventional expectation so that a violation is visible rather than absorbed.

The systematics set varies one optical input at a time over the family generated
by scripts/make_material_variants.py, holding geometry and surface fixed.  Its
output is not a number but a spread.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scint.registry import (  # noqa: E402
    Run,
    RunExistsError,
    checksum,
    create_run,
    record_status,
    run_id as _peek_id,
)

BINARY = ROOT / "build" / "sim" / "scint_optical"
RUNS_DIR = ROOT / "runs"

# Fixed for every gate: a 3" x 3" NaI(Tl) cylinder read out over its full face,
# the geometry the published 7.0 % FWHM refers to.
BASE_GEOMETRY = {
    "crystal": {
        "shape": "cylinder",
        "diameter_mm": 76.2,
        "length_mm": 76.2,
    },
    "surface": {
        "treatment": "polished",
        "wrapping": "teflon",
        "coupling": "air",
        "model": "lut",
        "reflectivity": 0.99,   # REFLECTANCE_440NM["teflon"]; kept in step below
    },
    # The index is now explicit rather than left at the application's default,
    # because it is one of the inputs this study varies: a light-collection
    # correction that fixes the coupling while scanning the crystal is only
    # half an answer.
    "coupling": {"type": "grease", "rindex": 1.465, "thickness_mm": 0.1},
    # Detection efficiency is deliberately 1: this study separates light
    # COLLECTION from photodetector quantum efficiency, and mixing them is one of
    # the reasons published light yields are not comparable.
    "readout": {"efficiency": 1.0},
    "source": {"energy_keV": 662.0, "offset_z_mm": -200.0},
    # Geant4 11.4.0 removed the clause of G4Scintillation::IsApplicable that
    # excluded optical photons (11.3.0 and the 11.4 beta had it), and
    # G4OpticalPhysics attaches the process to every particle IsApplicable
    # accepts. An optical photon absorbed in the bulk deposits its ~3 eV
    # (G4OpAbsorption::PostStepDoIt), and at 41 000 photons/MeV that deposit
    # scintillates: 0.12 new photons per absorbed one, a cascade worth +8.6 %
    # in the photon budget, +20 % in its variance and +40 % in mean detection
    # time, measured here on the release. It is not a model of NaI(Tl) -- the
    # yield per eV of absorbed light is not the yield per eV of ionisation --
    # so it is switched off, as a recorded input rather than a silent one.
    "physics": {"scintillation_from_optical_photons": False},
}

WRAPPINGS = ["none", "teflon", "lumirror", "tyvek", "tio", "esr"]

# Reflectance of each wrapping, as a flat number applied through REFLECTIVITY.
#
# Without it Geant4 treats every reflector as a perfect mirror -- the property
# defaults to 1.0 when absent, and the look-up-table models carry angular
# distributions only -- which is how an earlier campaign of this study ran. A
# referee found it. The values are the reflection coefficients at 440 nm in
# Table I of M. Janecek, IEEE Trans. Nucl. Sci. 59 (2012) 490, read from the
# printed table (the LBNL look-up tables these wrappings use come from the same
# group's measurements). 440 nm is the nearest tabulated wavelength to the
# 415 nm emission peak, and the same paper says which of these fail inside the
# band: ESR's reflectivity drops sharply below 395 nm and TiO2's below 420 nm,
# and PTFE reaches its value in 380-500 nm only above 0.5 mm thickness. A flat
# number therefore overstates ESR and TiO2 in the blue tail, and the
# reflectance scan in the scope set bounds what that is worth.
REFLECTANCE_440NM = {
    "teflon": 0.99,      # ACE Teflon tape (matte), n x 0.06 mm
    "esr": 0.985,        # 3M ESR film, 0.065 mm
    "lumirror": 0.98,    # Toray Lumirror, 0.24 mm
    "tyvek": 0.97,       # DuPont Tyvek paper, n x 0.11 mm
    "tio": 0.955,        # Saint-Gobain titanium dioxide paint, 0.14-0.18 mm
    "none": None,        # bare: no wrapper, nothing to reflect
}

# The arrangement's own optical inputs: geometry, readout coupling and surface
# finish. The baseline for all of them is BASE_GEOMETRY, so each of these runs
# differs from the baseline in exactly one thing, the same discipline the
# material scan follows.
#
# Geometry is here because the absorption result depends on optical path
# length: an absorption edge that costs a 76.2 mm crystal a great deal should
# cost a 25.4 mm crystal much less, and if it does not, the model is wrong.
# Coupling is here because a published measurement of the same crystal with and
# without optical grease differs by about 80 %, so leaving it fixed while
# calling the crystal's properties dominant would not be an argument.
SCOPE_VARIANTS: dict[str, dict[str, dict]] = {
    # The baseline itself, so the scan is self-contained. Making the coupling
    # index explicit in BASE_GEOMETRY changed every configuration hash without
    # changing any physics -- 1.465 is what DetectorConstruction.hh already
    # defaulted to -- so this run also proves that: it must reproduce the
    # earlier teflon light-collection efficiency exactly.
    "baseline":      {},
    # --- geometry, at fixed aspect ratio ---------------------------------
    "geom_1inch":    {"crystal": {"diameter_mm": 25.4, "length_mm": 25.4}},
    "geom_2inch":    {"crystal": {"diameter_mm": 50.8, "length_mm": 50.8}},
    "geom_4inch":    {"crystal": {"diameter_mm": 101.6, "length_mm": 101.6}},
    # --- geometry, at fixed diameter: path length alone -------------------
    "geom_long":     {"crystal": {"length_mm": 152.4}},
    "geom_short":    {"crystal": {"length_mm": 25.4}},
    # --- readout coupling -------------------------------------------------
    "couple_air":    {"coupling": {"type": "air", "rindex": 1.0}},
    "couple_gel146": {"coupling": {"rindex": 1.46}},
    "couple_n15":    {"coupling": {"rindex": 1.50}},
    "couple_n157":   {"coupling": {"rindex": 1.57}},
    # --- surface finish ---------------------------------------------------
    "finish_ground": {"surface": {"treatment": "ground"}},
    # --- wrapping reflectance ---------------------------------------------
    # The input the earlier campaign left at Geant4's default. 1.00 is that
    # default, a lossless mirror; 0.95 and 0.90 bracket a worn or thin Teflon
    # tape and the in-band cut-offs Janecek reports for ESR and TiO2.
    "reflectance_100": {"surface": {"reflectivity": 1.00}},
    "reflectance_095": {"surface": {"reflectivity": 0.95}},
    "reflectance_090": {"surface": {"reflectivity": 0.90}},
}

# The control pair described above. Keyed separately because they need a
# different material file as well as an override, and because they are a
# diagnostic rather than part of the scan proper.
SCOPE_CONTROLS: dict[str, tuple[str, dict[str, dict]]] = {
    "ctrl_flatabs_3inch": ("materials/variants/NaI_Tl_abs_flat2000.dat", {}),
    "ctrl_flatabs_1inch": ("materials/variants/NaI_Tl_abs_flat2000.dat",
                           {"crystal": {"diameter_mm": 25.4, "length_mm": 25.4}}),
}

# Bumped when the simulation's OUTPUT SCHEMA changes, so that runs made with a
# different set of recorded observables get different identifiers instead of
# colliding with the earlier ones. It is not a physics parameter: schema 1 and
# schema 2 with the same seed must produce identical photon counts, and that
# equality is itself a check that the added instrumentation perturbs nothing.
#   1: event, edep, generated, detected, first/mean detection time
#   2: + scintillation/Cerenkov split of the generated photons,
#        mean generated wavelength, mean detected wavelength
#   3: + photons removed by the optical step cap; and bare surfaces are now
#        resolved analytically rather than through a look-up table that does
#        not exist, which is what made the unwrapped configuration hang
#   4: surface.reflectivity is a REQUIRED input and part of the hash. Every run
#      before schema 4 was made with Geant4's silent default of a lossless
#      wrapper, and with a pre-release toolkit; schema 4 runs are made in the
#      container against the 11.4.2 release. Bumping the schema is what keeps
#      the two campaigns from ever being averaged or compared by accident.
#   5: physics.scintillation_from_optical_photons is a recorded input, off.
#      Schema-4 runs were made on the release with Geant4's new default of
#      letting absorbed optical photons scintillate; they are kept as the
#      evidence for the size of that effect and are compared with nothing.
CONFIG_SCHEMA = 5


def config_for(
    *,
    label: str,
    material_spec: str,
    events: int,
    seed: int,
    wrapping: str | None = None,
    overrides: dict[str, dict] | None = None,
) -> dict:
    """Build one run configuration.

    `overrides` is a shallow per-section update -- {"crystal": {"length_mm":
    25.4}} -- so that a scan over geometry or coupling needs no new keyword
    argument here. It is applied to the deep copy, so it changes the hash and
    therefore the run identity, which is the point: a different configuration
    must never reuse another one's results.
    """
    cfg = json.loads(json.dumps(BASE_GEOMETRY))  # deep copy
    cfg["crystal"]["material_spec"] = material_spec
    if wrapping is not None:
        cfg["surface"]["wrapping"] = wrapping
        cfg["surface"]["reflectivity"] = REFLECTANCE_440NM[wrapping]
    for section, values in (overrides or {}).items():
        cfg[section].update(values)
    if cfg["surface"]["wrapping"] == "none":
        cfg["surface"]["reflectivity"] = None
    # The label is deliberately NOT part of the configuration: the run id is a
    # hash of the physics, so two gates that ask for the same physics resolve to
    # the same run and the second is reused rather than recomputed.
    cfg["run"] = {"events": events, "seed": seed, "schema": CONFIG_SCHEMA}
    return {"label": label, "config": cfg}


def macro_for(cfg: dict, output_stem: Path) -> str:
    c, s, k, r, src, run = (
        cfg["crystal"], cfg["surface"], cfg["coupling"],
        cfg["readout"], cfg["source"], cfg["run"],
    )
    return f"""\
# Generated by scripts/run_gates.py -- do not edit; the configuration of record
# is config.yaml in this run directory, and its hash is the run identifier.
/random/setSeeds {run['seed']} 1

/scint/crystal/shape {c['shape']}
/scint/crystal/diameter {c['diameter_mm']} mm
/scint/crystal/length {c['length_mm']} mm
/scint/crystal/materialSpec {c['material_spec']}

/scint/surface/treatment {s['treatment']}
/scint/surface/wrapping {s['wrapping']}
/scint/surface/coupling {s['coupling']}
/scint/surface/model {s['model']}
{f"/scint/surface/reflectivity {s['reflectivity']}" if s.get('reflectivity') is not None else "# bare crystal: no wrapping reflectance"}

/scint/coupling/type {k['type']}
/scint/coupling/rindex {k['rindex']}
/scint/coupling/thickness {k['thickness_mm']} mm

/scint/readout/efficiency {r['efficiency']}

/scint/source/energy {src['energy_keV']} keV
/scint/source/offsetZ {src['offset_z_mm']} mm

/scint/output/file {output_stem}
/scint/output/format csv

/run/initialize
{"" if cfg.get("physics", {}).get("scintillation_from_optical_photons", True) else "/process/inactivate Scintillation opticalphoton"}
/run/beamOn {run['events']}
"""


def execute(item: tuple[str, Run]) -> tuple[str, bool, float]:
    label, run = item
    stem = run.path / "output"
    macro = run.path / "run.mac"
    macro.write_text(macro_for(run.config, stem))
    log = run.path / "geant4.log"

    # Hash the material file into the record. The configuration names the file;
    # only the hash pins its contents, and a material file that changed between
    # two runs would otherwise make them look comparable when they are not.
    material = ROOT / run.config["crystal"]["material_spec"]
    started = time.time()
    record_status(
        run, label=label, state="running", started_unix=started,
        material_sha256=checksum(material) if material.exists() else None,
    )
    with log.open("w") as fh:
        proc = subprocess.run(
            [str(BINARY), str(macro)], cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT
        )
    elapsed = time.time() - started

    ntuple = run.path / "output_nt_events.csv"
    ok = proc.returncode == 0 and ntuple.exists()
    record_status(
        run,
        state="done" if ok else "failed",
        returncode=proc.returncode,
        elapsed_s=round(elapsed, 1),
        output=str(ntuple.relative_to(ROOT)) if ntuple.exists() else None,
        output_sha256=checksum(ntuple) if ntuple.exists() else None,
    )
    return label, ok, elapsed


def gather(set_name: str, events: int, seed: int) -> list[dict]:
    configs: list[dict] = []
    if set_name in ("gates", "all"):
        # G2 and G3 come from the same baseline run: the photon budget and the
        # photopeak width are two readings of one dataset.
        configs.append(config_for(
            label="G2G3_baseline", material_spec="materials/NaI_Tl.dat",
            events=events, seed=seed,
        ))
        # G4: the wrapping scan.
        for wrap in WRAPPINGS:
            configs.append(config_for(
                label=f"G4_wrap_{wrap}", material_spec="materials/NaI_Tl.dat",
                events=events, seed=seed, wrapping=wrap,
            ))
    if set_name in ("systematics", "all"):
        for spec in sorted((ROOT / "materials" / "variants").glob("*.dat")):
            configs.append(config_for(
                label=f"SYS_{spec.stem.replace('NaI_Tl_', '')}",
                material_spec=str(spec.relative_to(ROOT)),
                events=events, seed=seed,
            ))
    if set_name in ("measured", "all"):
        # The digitised transmittance curve and its four digitisation-error
        # bounds -- two per axis, transmittance and wavelength -- run against the
        # same baseline as the arrangement scan so all of them are directly
        # comparable: same schema, same seed, same events.
        for label in ("abs_measured", "abs_measured_hi", "abs_measured_lo",
                      "abs_measured_lamhi", "abs_measured_lamlo",
                      "meas_rindex_jellison", "meas_rindex_flat185",
                      "meas_fwhm55", "meas_fwhm75", "meas_jacobian"):
            configs.append(config_for(
                label=f"MEAS_{label}",
                material_spec=f"materials/variants/NaI_Tl_{label}.dat",
                events=events, seed=seed,
            ))
        configs.append(config_for(
            label="MEAS_baseline", material_spec="materials/NaI_Tl.dat",
            events=events, seed=seed,
        ))
    if set_name in ("scope", "all"):
        # Everything in SCOPE_VARIANTS is an optical input of the arrangement
        # rather than of the crystal. They are scanned for the same reason the
        # material properties are: a study that varies the crystal's optical
        # properties while holding the surface, the coupling and the geometry
        # fixed cannot say the crystal's properties are what matters.
        for label, overrides in SCOPE_VARIANTS.items():
            configs.append(config_for(
                label=f"SCOPE_{label}", material_spec="materials/NaI_Tl.dat",
                events=events, seed=seed, overrides=overrides,
            ))
        for label, (spec, overrides) in SCOPE_CONTROLS.items():
            configs.append(config_for(
                label=f"SCOPE_{label}", material_spec=spec,
                events=events, seed=seed, overrides=overrides,
            ))
    return configs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", default="all",
                    choices=["gates", "systematics", "scope", "measured", "all"])
    ap.add_argument("--events", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--only", default=None,
                    help="comma-separated substrings; run only labels containing one of them")
    ap.add_argument("--binary", type=Path, default=None,
                    help="override the simulation binary (e.g. a build with a newer output schema)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    

    if args.binary is not None:
        global BINARY
        BINARY = args.binary.resolve()

    if not BINARY.exists():
        raise SystemExit(f"{BINARY} not built -- run cmake --build build/sim first")

    configs = gather(args.set, args.events, args.seed)
    if args.only:
        wanted = [w.strip() for w in args.only.split(",") if w.strip()]
        configs = [c for c in configs if any(w in c["label"] for w in wanted)]
    runs: list[tuple[str, Run]] = []
    seen: dict[str, str] = {}
    for item in configs:
        label, cfg = item["label"], item["config"]
        if args.dry_run:
            print(f"would run {label:32s} {_peek_id(cfg)}")
            continue
        # The run identifier hashes the configuration, which names the material
        # file by path. If that file's CONTENT has changed since the run was
        # made -- a regenerated variant -- the run on disk no longer represents
        # the configuration, and reusing it would pair a label with physics it
        # did not simulate. Such a run is moved aside, not deleted, and redone.
        rid = _peek_id(cfg)
        existing = RUNS_DIR / rid / "status.json"
        material = ROOT / cfg["crystal"]["material_spec"]
        if existing.exists() and material.exists():
            recorded = json.loads(existing.read_text()).get("material_sha256")
            current = checksum(material)
            if recorded and recorded != current:
                stale_dir = RUNS_DIR / "_stale"
                stale_dir.mkdir(exist_ok=True)
                target = stale_dir / f"{rid}-material-{recorded[:8]}"
                (RUNS_DIR / rid).rename(target)
                print(f"stale {label:32s} {rid}  material file changed since the run; "
                      f"moved to {target.relative_to(ROOT)}")
        try:
            run = create_run(cfg, runs_dir=RUNS_DIR, seed=cfg["run"]["seed"])
        except RunExistsError:
            first = seen.get(rid)
            note = f"identical physics to {first}" if first else "already on disk"
            print(f"reuse {label:32s} {rid}  ({note})")
            # Record the alias on the run so the analysis does not have to
            # re-derive it from this script later.
            st = json.loads(existing.read_text()) if existing.exists() else {}
            aliases = st.get("aliases", [])
            if label not in aliases and label != st.get("label"):
                aliases.append(label); st["aliases"] = aliases
                existing.write_text(json.dumps(st, indent=2) + "\n")
            continue
        seen[run.id] = label
        runs.append((label, run))
        print(f"queue {label:32s} {run.id}")

    if args.dry_run or not runs:
        return

    print(f"\nrunning {len(runs)} configurations on {args.jobs} workers "
          f"({args.events} events each)\n", flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for label, ok, elapsed in pool.map(execute, runs):
            print(f"{'ok  ' if ok else 'FAIL'} {label:32s} {elapsed:7.1f} s", flush=True)


if __name__ == "__main__":
    main()
