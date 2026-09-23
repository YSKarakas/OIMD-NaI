#!/usr/bin/env python3
"""Container checks on optical transport, recorded as evidence.

Four checks, each writing one JSON record to runs_scratch/diag/:

  reproduce  Re-run published runs with the archived binary in the container
             image as rebuilt from containers/Dockerfile, and compare each
             per-event output file with the checksum its run recorded.
  escapes    Track every optical photon of a few low-energy events in the
             polished, lossless companion of the closure gate (G5), and record
             how each one ends: detected, absorbed, or leaving the geometry, and
             for those that leave, where and at what incidence they last met
             the wrapped wall.
  lut        Read the LBNL look-up tables Geant4 ships (polished and ground
             Teflon) and record, for each incidence angle, the probability that
             DielectricLUT() draws the edge bin of its polar-angle grid, a
             direction lying in the surface.
  survey     For every configuration in the verdict, track a number of
             662 keV events and record how many photons leave the geometry, and
             an estimate of the light they would have delivered.

All four use the published binary (build/container/scint_optical) and the
published run macros, changing only the event count, the output path, the
tracking verbosity and, for `escapes`, the source energy. They need Docker and
the image scint:11.4.2; its identifier is recorded with every result, together
with the binary checksum and the commit of this repository.

    python3 scripts/transport_checks.py reproduce --runs 616c6dd12c043688 ...
    python3 scripts/transport_checks.py escapes --events 3
    python3 scripts/transport_checks.py lut
    python3 scripts/transport_checks.py survey --events 10 --jobs 2
"""

from __future__ import annotations

import argparse
import bisect
import concurrent.futures as cf
import datetime
import hashlib
import json
import math
import re
import shutil
import subprocess
import zlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGE = "scint:11.4.2"
BINARY = "build/container/scint_optical"
DIAG = ROOT / "runs_scratch" / "diag"
# Scratch space inside the repository, because the container sees the
# repository and nothing else; ignored by git.
WORK = ROOT / "runs_scratch" / "transport_checks_tmp"
VERDICT = ROOT / "runs_scratch" / "gates_final.json"
RADIUS_MM = 38.1          # the 3" x 3" crystal every published run uses
HALF_LENGTH_MM = 38.1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance() -> dict:
    image = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{.Id}}"],
                           capture_output=True, text=True, check=True).stdout.strip()
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                            cwd=ROOT).stdout.strip()
    # Dirty means a tracked file differs from the commit. Untracked files do not
    # count: the records this script writes are themselves untracked until
    # they are committed, and would otherwise mark every record after the first.
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, cwd=ROOT).stdout.strip())
    return {"image": f"{IMAGE}@{image}", "binary": BINARY,
            "binary_sha256": sha256(ROOT / BINARY), "git_commit": commit, "git_dirty": dirty,
            "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat()}


def macro_for(run_id: str, out_stem: str, events: int | None = None,
              verbose: bool = False, energy_keV: float | None = None) -> str:
    """The published macro of a run, redirected; nothing else is changed."""
    mac = (ROOT / "runs" / run_id / "run.mac").read_text()
    mac = mac.replace(f"/scint/output/file /work/runs/{run_id}/output",
                      f"/scint/output/file /work/{out_stem}")
    if energy_keV is not None:
        mac, n = re.subn(r"/scint/source/energy [\d.]+ keV", f"/scint/source/energy {energy_keV} keV", mac)
        assert n == 1, "source energy line not found"
    if events is not None or verbose:
        m = re.search(r"/run/beamOn (\d+)", mac)
        assert m, "beamOn not found"
        n_ev = events if events is not None else int(m.group(1))
        mac = mac.replace(m.group(0), ("/tracking/verbose 1\n" if verbose else "") + f"/run/beamOn {n_ev}")
    return mac


def run_in_container(mac_rel: str, shell_tail: str = "") -> None:
    """Run the published binary on a macro inside the image; `shell_tail` is
    appended to the command line (a pipe into a filter, a redirection)."""
    cmd = f". /opt/geant4/bin/geant4.sh && /work/{BINARY} /work/{mac_rel} {shell_tail}"
    subprocess.run(["docker", "run", "--rm", "-v", f"{ROOT}:/work", "-w", "/work", IMAGE,
                    "bash", "-c", cmd], check=True)


def verdict_runs() -> list[tuple[str, str]]:
    v = json.loads(VERDICT.read_text())
    runs, seen = [], set()
    for r in v["runs"]:
        if r["id"] in seen or r.get("alias_of"):
            continue
        seen.add(r["id"])
        runs.append((r["id"], r["label"]))
    return runs


def label_of(label: str) -> str:
    return next(r for r, lab in verdict_runs() if lab == label)


# --------------------------------------------------------------------------- #
# reproduce
# --------------------------------------------------------------------------- #

def cmd_reproduce(args) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for rid in args.runs:
        status = json.loads((ROOT / "runs" / rid / "status.json").read_text())
        stem = f"runs_scratch/transport_checks_tmp/repro_{rid}"
        (ROOT / f"{stem}.mac").write_text(macro_for(rid, stem))
        run_in_container(f"{stem}.mac", f"> /work/{stem}.log 2>&1")
        out = ROOT / f"{stem}_nt_events.csv"
        got = sha256(out)
        rows.append({"id": rid, "label": status["label"],
                     "published_image": status.get("container_image"),
                     "published_output_sha256": status["output_sha256"],
                     "reproduced_output_sha256": got,
                     "byte_identical": got == status["output_sha256"]})
        print(f"{rid} {status['label']:24s} {'BYTE-IDENTICAL' if rows[-1]['byte_identical'] else 'DIFFERS'}")
        out.unlink(); (ROOT / f"{stem}.log").unlink(); (ROOT / f"{stem}.mac").unlink()
    write(DIAG / "rebuilt_image_reproduction.json", {"provenance": provenance(), "runs": rows})


# --------------------------------------------------------------------------- #
# escapes
# --------------------------------------------------------------------------- #

TRACK = re.compile(r"Particle = (\S+),\s+Track ID = (\d+)")
STEP = re.compile(r"^\s*\d+\s+-?[\d.]")


def tracks_from_log(path: Path):
    """Yield (particle, [step token lists]) for every track of a verbose log."""
    cur = None
    with path.open(errors="replace") as fh:
        for line in fh:
            m = TRACK.search(line)
            if m:
                if cur:
                    yield cur
                cur = (m.group(1), [])
                continue
            if cur is not None and STEP.match(line):
                cur[1].append(line.split())
    if cur:
        yield cur


def incidence_on_wall(prev, cur) -> float:
    """Angle between the step that ends at the lateral wall and the outward normal."""
    x0, y0, z0 = map(float, prev[1:4]); x1, y1, z1 = map(float, cur[1:4])
    d = (x1 - x0, y1 - y0, z1 - z0); norm = math.sqrt(sum(a * a for a in d)) or 1.0
    r = math.hypot(x1, y1) or 1.0
    c = (d[0] * x1 / r + d[1] * y1 / r) / norm
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def on_lateral_wall(tok) -> bool:
    return abs(math.hypot(float(tok[1]), float(tok[2])) - RADIUS_MM) < 0.06


def cmd_escapes(args) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    rid = label_of("G5_trapping_polished")
    stem = "runs_scratch/transport_checks_tmp/escapes"
    (ROOT / f"{stem}.mac").write_text(macro_for(rid, stem, events=args.events, verbose=True,
                                                energy_keV=args.energy))
    run_in_container(f"{stem}.mac", f"> /work/{stem}.log 2>&1")
    fate, where, inc_escape, inc_reflect = Counter(), Counter(), [], []
    photons = 0
    for particle, steps in tracks_from_log(ROOT / f"{stem}.log"):
        if particle != "opticalphoton" or not steps:
            continue
        photons += 1
        nv, proc = steps[-1][-2], steps[-1][-1]
        f = ("detected" if nv == "Photodetector" else "left the geometry" if nv == "OutOfWorld"
             else "absorbed in bulk" if proc == "OpAbsorption" else f"other: {nv}/{proc}")
        fate[f] += 1
        if f == "left the geometry" and len(steps) >= 3:
            last = steps[-2]
            if on_lateral_wall(last):
                where["lateral wall"] += 1
                inc_escape.append(incidence_on_wall(steps[-3], last))
            elif abs(float(last[3]) + HALF_LENGTH_MM) < 0.06:
                where["back face"] += 1
            else:
                where["elsewhere"] += 1
        # ordinary reflections: a lateral-wall hit followed by a zero-length re-entry
        for i in range(1, len(steps) - 1):
            if (steps[i][-2] == "World" and on_lateral_wall(steps[i])
                    and steps[i + 1][6] == "0" and steps[i + 1][-2] == "Crystal"):
                inc_reflect.append(incidence_on_wall(steps[i - 1], steps[i]))

    def summary(a):
        a = sorted(a)
        return {"n": len(a), "median_deg": a[len(a) // 2] if a else None,
                "fraction_above_80_deg": sum(v >= 80 for v in a) / len(a) if a else None}

    rows = [ln.split(",") for ln in (ROOT / f"{stem}_nt_events.csv").read_text().splitlines()
            if ln and not ln.startswith("#")]
    result = {"provenance": provenance(), "configuration": "G5_trapping_polished", "run_id": rid,
              "events": args.events, "energy_keV": args.energy, "photons": photons,
              "ntuple_scintillation": sum(int(r[3]) for r in rows),
              "ntuple_detected": sum(int(r[5]) for r in rows),
              "fate": dict(fate), "left_the_geometry_from": dict(where),
              "incidence_on_the_step_before_leaving": summary(inc_escape),
              "incidence_of_ordinary_reflections": summary(inc_reflect)}
    for p in (f"{stem}.log", f"{stem}.mac", f"{stem}_nt_events.csv"):
        (ROOT / p).unlink()
    write(DIAG / "closure_polished_escapes.json", result)


# --------------------------------------------------------------------------- #
# lut
# --------------------------------------------------------------------------- #

N_INC, N_THETA, N_PHI = 91, 45, 37     # G4OpticalSurface.hh, 11.4.2


def cmd_lut(args) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    rel = "runs_scratch/transport_checks_tmp"
    subprocess.run(["docker", "run", "--rm", "-v", f"{ROOT}:/work", IMAGE, "bash", "-c",
                    # The 11.4 installation compiles the dataset paths in rather than
                    # exporting G4REALSURFACEDATA, so the path is asked of geant4-config.
                    ". /opt/geant4/bin/geant4.sh && d=$(geant4-config --datasets | awk '$1==\"RealSurface\"{print $3}') "
                    f"&& cp \"$d\"/PolishedTeflon.z \"$d\"/GroundTeflon.z /work/{rel}/ && basename \"$d\""],
                   check=True, capture_output=True, text=True)
    tables = {}
    for name in ("PolishedTeflon.z", "GroundTeflon.z"):
        path = ROOT / rel / name
        values = [float(x) for x in zlib.decompress(path.read_bytes()).decode().split()]
        assert len(values) >= N_INC * N_THETA * N_PHI
        rows = {}
        for inc in range(80, 91):
            # DielectricLUT() draws theta index in [0, N_THETA-2] and phi index in
            # [0, N_PHI-2] (G4RandFlat::shootInt(max - 1)); theta index 0 is -90 deg.
            per_theta = [0.0] * (N_THETA - 1)
            for phi in range(N_PHI - 1):
                for th in range(N_THETA - 1):
                    per_theta[th] += values[inc + th * N_INC + phi * N_THETA * N_INC]
            total = sum(per_theta)
            best = max(range(N_THETA - 1), key=lambda t: per_theta[t])
            rows[str(inc)] = {"p_in_surface_bin": per_theta[0] / total if total else None,
                              "most_probable_theta_deg": -90 + 4 * best}
        tables[name] = {"sha256": sha256(path), "by_incidence_deg": rows}
        path.unlink()
    write(DIAG / "lut_grazing.json", {"provenance": provenance(), "tables": tables})


# --------------------------------------------------------------------------- #
# survey
# --------------------------------------------------------------------------- #

AWK = (r"""/Particle = / { if (p != "") print p, last; s = $0; sub(/.*Particle = /, "", s); """
       r"""sub(/,.*/, "", s); p = s; last = "0 ? ?"; next } """
       r"""/^ *[0-9]+ +-?[0-9.]/ { last = $1 " " $(NF-1) " " $NF; next } """
       r"""END { if (p != "") print p, last }""")


def classify(next_volume: str, process: str) -> str:
    if next_volume == "Photodetector":
        return "detected"
    if next_volume == "OutOfWorld":
        return "left the geometry"
    if process == "OpAbsorption":
        return "absorbed in bulk"
    if next_volume == "World":
        return "absorbed at the wrapper"
    return f"other: {next_volume}/{process}"


def survey_one(run, events: int) -> dict:
    rid, label = run
    stem = f"runs_scratch/transport_checks_tmp/survey_{rid}"
    (ROOT / f"{stem}.mac").write_text(macro_for(rid, stem, events=events, verbose=True))
    # Only each track's particle and last step are kept: step number, next
    # volume, process. The full history would be gigabytes.
    run_in_container(f"{stem}.mac",
                     f"2>&1 | awk -f /work/runs_scratch/transport_checks_tmp/last_step.awk > /work/{stem}.tracks")
    fates = []
    for ln in (ROOT / f"{stem}.tracks").read_text().splitlines():
        t = ln.split()
        if len(t) == 4 and t[0] == "opticalphoton":
            fates.append((classify(t[2], t[3]), int(t[1]) if t[1].isdigit() else 0))
    counts = Counter(f for f, _ in fates)
    # Estimate of the detections the leak removed: each photon that left is
    # credited with the probability that a photon which did not leave, and
    # survived at least as many steps, was detected.
    kept = sorted(((f, s) for f, s in fates if f != "left the geometry"), key=lambda x: x[1])
    steps = [s for _, s in kept]
    suffix = [0] * (len(kept) + 1)
    for i in range(len(kept) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + (kept[i][0] == "detected")
    lost = 0.0
    for f, s in fates:
        if f == "left the geometry":
            i = bisect.bisect_left(steps, s)
            m = len(kept) - i
            lost += suffix[i] / m if m else 0.0
    rows = [ln.split(",") for ln in (ROOT / f"{stem}_nt_events.csv").read_text().splitlines()
            if ln and not ln.startswith("#")]
    for p in (f"{stem}.tracks", f"{stem}.mac", f"{stem}_nt_events.csv"):
        (ROOT / p).unlink()
    return {"id": rid, "label": label, "events": events, "photons": len(fates),
            "fate": dict(counts), "lost_detections_estimate": lost,
            "ntuple_scintillation": sum(int(r[3]) for r in rows),
            "ntuple_detected": sum(int(r[5]) for r in rows)}


def cmd_survey(args) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "last_step.awk").write_text(AWK)
    runs = verdict_runs()
    results = []
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for res in ex.map(lambda r: survey_one(r, args.events), runs):
            results.append(res)
            n = res["photons"] or 1
            print(f"{res['label']:28s} photons {n:7d}  left {res['fate'].get('left the geometry', 0) / n:.4%}"
                  f"  lost-detection estimate {res['lost_detections_estimate'] / n:.4%}", flush=True)
    write(DIAG / "leak_survey.json", {"provenance": provenance(), "events_per_configuration": args.events,
                                      "configurations": results})


def write(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1) + "\n")
    print(f"wrote {path.relative_to(ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("reproduce"); p.add_argument("--runs", nargs="+", required=True)
    p = sub.add_parser("escapes"); p.add_argument("--events", type=int, default=3)
    p.add_argument("--energy", type=float, default=60.0, help="source energy [keV]")
    sub.add_parser("lut")
    p = sub.add_parser("survey"); p.add_argument("--events", type=int, default=10)
    p.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    if shutil.which("docker") is None:
        raise SystemExit("docker not found")
    {"reproduce": cmd_reproduce, "escapes": cmd_escapes, "lut": cmd_lut, "survey": cmd_survey}[args.cmd](args)


if __name__ == "__main__":
    main()
