#!/usr/bin/env python3
"""The analytic expectation behind gate G1, recorded with its inputs.

G1 compares the simulated probability that a 662 keV photon interacts in the
76.2 mm crystal with 1 - exp(-L/lambda). The attenuation length is Geant4's own
(G4EmCalculator, G4EmStandardPhysics_option4), and it includes coherent
(Rayleigh) scattering, which deposits no energy and so cannot count as an
interaction in the simulated tally. This script asks the container's Geant4
11.4.2 for the four terms of the attenuation coefficient separately and writes
both expectations -- with and without Rayleigh -- beside the simulated value
and its binomial uncertainty, to runs_scratch/diag/g1_expectation.json.

The Geant4 query is made by sim/tools/g4data.cc, built in the container into
build/container-tools/ (not build/container, whose simulation binary is the
archived one):

    docker run --rm -v "$PWD":/work -w /work scint:11.4.2 bash -c \\
        '. /opt/geant4/bin/geant4.sh && cmake -S sim -B build/container-tools \\
         -DCMAKE_BUILD_TYPE=Release && cmake --build build/container-tools --target g4data'
    python3 scripts/g1_expectation.py
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGE = "scint:11.4.2"
TOOL = "build/container-tools/g4data"
OUT = ROOT / "runs_scratch" / "diag" / "g1_expectation.json"
SCRATCH = "runs_scratch/transport_checks_tmp/g1_attenuation.csv"
# The crystal: materials/NaI_Tl.dat (density and exact stoichiometric mass
# fractions), and the 76.2 mm length the photon crosses on the axis.
DENSITY, MASSFRAC, LENGTH_CM = "3.67", "Na:0.153373922,I:0.846626078", 7.62


def main() -> None:
    # The container writes into the repository's git-ignored scratch space,
    # which a fresh checkout does not have yet.
    (ROOT / SCRATCH).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["docker", "run", "--rm", "-v", f"{ROOT}:/work", "-w", "/work", IMAGE, "bash", "-c",
                    f". /opt/geant4/bin/geant4.sh && /work/{TOOL} attenuation --name NaI "
                    f"--density {DENSITY} --massfrac {MASSFRAC} --energies 0.662 --out /work/{SCRATCH}"],
                   check=True, capture_output=True)
    lines = (ROOT / SCRATCH).read_text().splitlines()
    g4_line = next(ln for ln in lines if ln.startswith("# Photon attenuation computed by Geant4"))
    row = next(csv.DictReader(ln for ln in lines if not ln.startswith("#")))
    (ROOT / SCRATCH).unlink()
    mu = float(row["mu_per_cm"])
    terms = {k: float(row[f"mu_{k}_per_cm"]) for k in ("conv", "compt", "phot", "rayl")}
    mu_no_rayl = mu - terms["rayl"]

    verdict = json.loads((ROOT / "runs_scratch" / "gates_final.json").read_text())
    p_sim = float(verdict["header"]["G1_interaction_probability"])
    n_ev = int(verdict["header"]["events"][0])
    se = math.sqrt(p_sim * (1 - p_sim) / n_ev)
    with_r, without_r = 1 - math.exp(-mu * LENGTH_CM), 1 - math.exp(-mu_no_rayl * LENGTH_CM)

    image = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{.Id}}"],
                           capture_output=True, text=True, check=True).stdout.strip()
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, cwd=ROOT).stdout.strip())
    record = {
        "provenance": {"image": f"{IMAGE}@{image}", "tool": TOOL,
                       "tool_sha256": hashlib.sha256((ROOT / TOOL).read_bytes()).hexdigest(),
                       "geant4": g4_line.split("Geant4", 1)[1].strip(),
                       "git_commit": commit, "git_dirty": dirty,
                       "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat()},
        "material": {"density_g_cm3": float(DENSITY), "mass_fractions": MASSFRAC, "length_cm": LENGTH_CM},
        "energy_MeV": 0.662,
        "mu_per_cm": {"total": mu, **terms},
        "rayleigh_fraction_of_total": terms["rayl"] / mu,
        "expected_interaction_probability": {"with_rayleigh": with_r, "without_rayleigh": without_r},
        "simulated": {"interaction_probability": p_sim, "events": n_ev, "binomial_standard_error": se,
                      "source": "runs_scratch/gates_final.json header"},
        "simulated_minus_expected_in_standard_errors": {"with_rayleigh": (p_sim - with_r) / se,
                                                        "without_rayleigh": (p_sim - without_r) / se},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + "\n")
    print(f"mu = {mu:.6f} /cm, Rayleigh {100 * terms['rayl'] / mu:.2f} % of it")
    print(f"expected {with_r:.4f} with Rayleigh, {without_r:.4f} without; simulated {p_sim:.3f} +- {se:.4f}")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
