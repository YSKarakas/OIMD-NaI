# Cost-Efficient Scintillator Discovery Platform

A reproducible computational pipeline for screening scintillator compositions under
**cost** and **critical-raw-material (CRM) supply-risk** constraints, with detector-level
validation in Geant4.

## Why

Reported scintillator light yield is **not an intrinsic material property**. It is the
joint output of *material × geometry × surface finish × readout chain*. Two peer-reviewed
compilations disagree by 13–18 % on the same materials (CsI:Tl 54,000 vs 61,000 ph/MeV;
BGO 10,000 vs 8,500 ph/MeV), and in the perovskite literature powdered samples show up to
a 3× apparent enhancement over polished single crystals purely from scattering disrupting
waveguiding. Consequently **no screening study currently has trustworthy input data.**

> You cannot search for a cheaper scintillator until you can compare scintillators fairly.

## Scope

| Stage | Claim |
|---|---|
| **Paper 0** | A Geant4-based geometry-correction protocol for reported light yields |
| **Paper 1** | An open, benchmarked, cost- and CRM-aware screening framework |
| **Paper 2** | Ga-reduced garnet candidates (LHCb PicoCal / CERN DRD6 context) |
| **Paper 3** | Synthesis and characterisation |

This repository currently implements the infrastructure for Papers 0 and 1.

## What Geant4 can and cannot do here

**It cannot predict light yield from chemistry.** Geant4 is a radiation-transport code; the
scintillation yield is an *input*. This project uses Geant4 for what it is uniquely good at:
optical photon transport, light-collection efficiency, self-absorption, surface and wrapping
effects, and detector-level energy resolution. Intrinsic properties come from electronic-structure
databases and published measurements, always with stated uncertainty.

## Repository layout

```
docs/      Strategy and literature documents (Turkish)
scint/     Python package: materials, optical models, cost/CRM, run registry, runner, reporting
sim/       Geant4 C++ application and helper tools
scripts/   Derivations and drivers — each one re-runnable and self-documenting
configs/   Run and sweep configurations (YAML)
materials/ Geant4 material specifications (GENERATED — see below)
data/      Sourced input data — every numeric value carries a source and a date
runs/      Run outputs (git-ignored; regenerable from configs)
site/      Generated static report (git-ignored)
tests/     pytest suite
```

## Reproducibility contract

Every simulation run is a directory under `runs/<run_id>/` containing:

| File | Contents |
|---|---|
| `config.yaml` | The complete input. Nothing outside this file influences the result. |
| `env.json` | git SHA (+dirty flag), Geant4/ROOT/Python versions, OS, RNG seed, timestamps |
| `status.json` | Exit state, wall time, output checksums |
| outputs | ROOT / CSV |

`run_id = sha256(canonical(config))[:16]`, so the same configuration always maps to the same
run and the input↔output correspondence can never be lost. `runs/` is **not** committed —
it is regenerable by construction.

**Data provenance rule:** no unsourced number enters `data/`. A value without a source and a
quotation date is recorded as `null` and surfaces as missing in reports.

**Material files are generated, not written.** `materials/*.dat` comes out of
`scripts/make_material_variants.py`, which builds each optical property from a named model in
`scint/optical.py` and writes that model's full provenance into the file beside the numbers.
Editing a `.dat` by hand is caught by the test suite, because the provenance comment and the
value it describes must not be able to drift apart.

**Where a value has no single trustworthy source, it is scanned, not chosen.** Three of the
optical inputs for NaI:Tl are in that position: the refractive index descends from a single
1923 measurement and a 2012 direct measurement disagrees with it by 1.6 % at the emission
peak; the bulk absorption length has exactly one usable published anchor; and no
machine-readable emission spectrum exists in the open literature at all. Rather than pick a
number, the study runs a family of variants that brackets what the literature asserts and
reports the resulting spread as a systematic. See `docs/05_OPTIK_GIRDILER.md` and
`data/literature/optical_provenance.csv`.

### Running the validation gates

```bash
python3 scripts/derive_nai_abslength.py       # the absorption-length derivation, from the paper's own numbers
python3 scripts/make_material_variants.py     # regenerate materials/ from scint/optical.py
python3 scripts/run_gates.py --set all        # G2-G4 plus the optical systematics family
python3 scripts/analyse_gates.py              # verdicts
```

## Setup

Requires Geant4 11.4.0, ROOT 6.36, CMake ≥ 3.20, a C++17 compiler, Python 3.12.

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e . --config-settings editable_mode=compat
cmake -S sim -B build/sim -DCMAKE_BUILD_TYPE=Release && cmake --build build/sim -j8
./build/sim/g4data elements --out data/derived/atomic_data.csv
```

`editable_mode=compat` is not optional: setuptools' default editable install uses an
import-hook `.pth` that did not load at interpreter startup on this setup, so scripts
failed to import `scint` while `python -c` appeared to work (it was picking the package
up from the working directory). The compat mode writes a plain path entry instead.

A pinned Docker image (`containers/Dockerfile`) reproduces the full toolchain; see
`containers/README.md`.

### Check your Geant4 before trusting a number

```bash
./build/sim/g4data version
```

`geant4-config --version` **cannot tell a beta from a release**: the beta of a series
already carries the target `G4VERSION_NUMBER`, so both report e.g. `11.4.0`. Only the
version *tag* discriminates. Every run records the tag and an `is_prerelease` flag in
`env.json`, and `scripts/validate_layer0.py` prints a warning banner when the toolchain
is a pre-release build. Attenuation caches are keyed by the tag, so results from one
Geant4 build are never served to another.

## Licence

TBD — to be chosen in consultation with the host institution before any public release.
