"""Execution layer: a backend abstraction and parameter-sweep expansion.

Only the local backend is implemented. The cluster backends are deliberate
stubs: fixing the interface now -- while there is exactly one implementation to
shape it -- means adding SLURM or HTCondor later is a contained change rather
than a refactor of everything that calls into here.

A browser never launches simulations. Runs are started from the command line or
from Python, write their results into the run registry, and the static report is
generated from that registry afterwards. Keeping execution and presentation
apart is what lets the same pipeline work on a laptop and on a batch farm.
"""

from __future__ import annotations

import itertools
import subprocess
import time
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from scint.registry import Run, create_run, record_status

SWEEP_KEY = "sweep"


@dataclass(frozen=True)
class RunResult:
    """Outcome of executing one run."""

    run: Run
    returncode: int
    wall_seconds: float
    stdout_path: Path
    stderr_path: Path

    @property
    def ok(self) -> bool:
        return self.returncode == 0


class Backend(ABC):
    """Somewhere a run can be executed."""

    name: str = "abstract"

    @abstractmethod
    def execute(self, run: Run, command: Sequence[str], env: dict[str, str] | None = None) -> RunResult:
        """Run ``command`` with the run directory as working directory."""

    def execute_many(
        self,
        jobs: Sequence[tuple[Run, Sequence[str]]],
        env: dict[str, str] | None = None,
    ) -> list[RunResult]:
        """Execute several runs. The default is sequential; backends may override."""
        return [self.execute(run, command, env) for run, command in jobs]


class LocalBackend(Backend):
    """Execute runs as local subprocesses, optionally several at once."""

    name = "local"

    def __init__(self, max_workers: int = 1) -> None:
        if max_workers < 1:
            raise ValueError(f"max_workers must be at least 1, got {max_workers}")
        self.max_workers = max_workers

    def execute(self, run: Run, command: Sequence[str], env: dict[str, str] | None = None) -> RunResult:
        stdout_path = run.path / "stdout.log"
        stderr_path = run.path / "stderr.log"
        record_status(run, state="running", backend=self.name, command=list(command))

        started = time.monotonic()
        with stdout_path.open("w") as out, stderr_path.open("w") as err:
            completed = subprocess.run(
                list(command), cwd=run.path, stdout=out, stderr=err, env=env
            )
        wall = time.monotonic() - started

        record_status(
            run,
            state="completed" if completed.returncode == 0 else "failed",
            returncode=completed.returncode,
            wall_seconds=round(wall, 3),
        )
        return RunResult(
            run=run,
            returncode=completed.returncode,
            wall_seconds=wall,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )

    def execute_many(self, jobs, env=None) -> list[RunResult]:
        if self.max_workers == 1:
            return super().execute_many(jobs, env)
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = [pool.submit(self.execute, run, command, env) for run, command in jobs]
            return [f.result() for f in futures]


class SlurmBackend(Backend):
    """Placeholder for SLURM submission (e.g. TU Dortmund LiDO)."""

    name = "slurm"

    def execute(self, run, command, env=None):
        raise NotImplementedError(
            "SlurmBackend is not implemented yet. Runs currently execute locally; "
            "see scint/runner.py for the interface a batch backend must satisfy."
        )


class CondorBackend(Backend):
    """Placeholder for HTCondor submission (e.g. CERN lxplus)."""

    name = "htcondor"

    def execute(self, run, command, env=None):
        raise NotImplementedError(
            "CondorBackend is not implemented yet. Runs currently execute locally; "
            "see scint/runner.py for the interface a batch backend must satisfy."
        )


BACKENDS: dict[str, type[Backend]] = {
    LocalBackend.name: LocalBackend,
    SlurmBackend.name: SlurmBackend,
    CondorBackend.name: CondorBackend,
}


def expand_sweep(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand a sweep specification into one configuration per point.

    A sweep is a base configuration plus a ``sweep`` mapping of dotted paths to
    value lists; the result is their Cartesian product::

        base:
          geometry: {shape: box, length_mm: 50}
        sweep:
          geometry.surface: [polished, ground]
          geometry.wrapping: [teflon, vm2000]

    yields four configurations. Each carries no trace of the sweep itself, so a
    configuration produced by a sweep and the same configuration written by hand
    hash to the same run identifier.
    """
    base = {k: v for k, v in spec.items() if k != SWEEP_KEY}
    axes = spec.get(SWEEP_KEY) or {}
    if not axes:
        return [base]

    for path, values in axes.items():
        if not isinstance(values, list) or not values:
            raise ValueError(f"sweep axis {path!r} must be a non-empty list, got {values!r}")

    names = list(axes)
    configs: list[dict[str, Any]] = []
    for combination in itertools.product(*(axes[name] for name in names)):
        config = _deepcopy(base)
        for path, value in zip(names, combination):
            _set_path(config, path, value)
        configs.append(config)
    return configs


def _deepcopy(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _deepcopy(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_deepcopy(v) for v in value]
    return value


def _set_path(config: dict[str, Any], dotted: str, value: Any) -> None:
    keys = dotted.split(".")
    node = config
    for key in keys[:-1]:
        nxt = node.get(key)
        if not isinstance(nxt, dict):
            nxt = {}
            node[key] = nxt
        node = nxt
    node[keys[-1]] = value


def prepare_runs(
    configs: Sequence[dict[str, Any]],
    *,
    runs_dir: Path | None = None,
    seed: int | None = None,
    skip_existing: bool = True,
) -> tuple[list[Run], list[str]]:
    """Create run directories, returning (new runs, identifiers already present).

    Skipping existing runs is the default because the identifier is a content
    hash: an existing directory means this exact configuration has already been
    executed, and a sweep that partly overlaps an earlier one should only do the
    new work.
    """
    created: list[Run] = []
    skipped: list[str] = []
    for config in configs:
        try:
            created.append(create_run(config, runs_dir=runs_dir, seed=seed))
        except Exception as exc:  # RunExistsError, kept broad for a clear message
            if skip_existing and type(exc).__name__ == "RunExistsError":
                from scint.registry import run_id

                skipped.append(run_id(config))
                continue
            raise
    return created, skipped
