#!/usr/bin/env python3
"""Run the optical validation gates G2-G4 and the optical-input systematics scan.

Each configuration becomes a registry run: the directory name is the hash of the
configuration, so a result can never be separated from the inputs that produced
it, and re-running an identical configuration is refused rather than silently
overwriting the earlier result.

    python3 scripts/run_gates.py --set gates        # G2/G3 baseline + G4 wrappings
    python3 scripts/run_gates.py --set systematics  # the optical-input family
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
    },
    "coupling": {"type": "grease", "thickness_mm": 0.1},
    # Detection efficiency is deliberately 1: this study separates light
    # COLLECTION from photodetector quantum efficiency, and mixing them is one of
    # the reasons published light yields are not comparable.
    "readout": {"efficiency": 1.0},
    "source": {"energy_keV": 662.0, "offset_z_mm": -200.0},
}

WRAPPINGS = ["none", "teflon", "lumirror", "tyvek", "tio", "esr"]

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
CONFIG_SCHEMA = 3


def config_for(
    *,
    label: str,
    material_spec: str,
    events: int,
    seed: int,
    wrapping: str | None = None,
) -> dict:
    cfg = json.loads(json.dumps(BASE_GEOMETRY))  # deep copy
    cfg["crystal"]["material_spec"] = material_spec
    if wrapping is not None:
        cfg["surface"]["wrapping"] = wrapping
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

/scint/coupling/type {k['type']}
/scint/coupling/thickness {k['thickness_mm']} mm

/scint/readout/efficiency {r['efficiency']}

/scint/source/energy {src['energy_keV']} keV
/scint/source/offsetZ {src['offset_z_mm']} mm

/scint/output/file {output_stem}
/scint/output/format csv

/run/initialize
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
    return configs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", default="all", choices=["gates", "systematics", "all"])
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
        try:
            run = create_run(cfg, runs_dir=RUNS_DIR, seed=cfg["run"]["seed"])
        except RunExistsError:
            first = seen.get(_peek_id(cfg))
            note = f"identical physics to {first}" if first else "already on disk"
            print(f"reuse {label:32s} {_peek_id(cfg)}  ({note})")
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
