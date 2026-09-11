"""Provenance and execution-layer behaviour."""

from __future__ import annotations

import json

import pytest

from scint.registry import (
    RunExistsError,
    canonical_yaml,
    checksum,
    create_run,
    iter_runs,
    load_run,
    record_status,
    run_id,
)
from scint.runner import LocalBackend, SlurmBackend, expand_sweep, prepare_runs


# --------------------------------------------------------------------------- #
# Identity
# --------------------------------------------------------------------------- #

def test_run_id_is_independent_of_key_order():
    a = {"b": 2, "a": {"y": 1, "x": [3, 1]}}
    b = {"a": {"x": [3, 1], "y": 1}, "b": 2}
    assert run_id(a) == run_id(b)


def test_run_id_changes_when_any_value_changes():
    base = {"energy_keV": 662, "material": "NaI"}
    assert run_id(base) != run_id({**base, "energy_keV": 663})


def test_run_id_respects_list_order():
    """Lists are ordered data; reordering them is a different configuration."""
    assert run_id({"x": [1, 2]}) != run_id({"x": [2, 1]})


def test_canonical_yaml_is_stable_across_calls():
    config = {"z": 1, "a": {"c": 3, "b": 2}}
    assert canonical_yaml(config) == canonical_yaml(dict(reversed(list(config.items()))))


# --------------------------------------------------------------------------- #
# Run directories
# --------------------------------------------------------------------------- #

def test_create_run_writes_the_full_contract(tmp_path):
    run = create_run({"material": "NaI", "energy_keV": 662}, runs_dir=tmp_path, seed=7)
    assert run.config_path.exists() and run.env_path.exists()
    env = json.loads(run.env_path.read_text())
    assert env["seed"] == 7
    assert env["python"]["version"]
    assert "packages" in env["python"]
    assert env["platform"]["system"]


def test_create_run_refuses_to_clobber_an_existing_run(tmp_path):
    config = {"material": "BGO"}
    create_run(config, runs_dir=tmp_path)
    with pytest.raises(RunExistsError):
        create_run(config, runs_dir=tmp_path)


def test_create_run_can_reuse_when_asked(tmp_path):
    config = {"material": "BGO"}
    first = create_run(config, runs_dir=tmp_path)
    second = create_run(config, runs_dir=tmp_path, exist_ok=True)
    assert first.id == second.id


def test_record_status_merges_rather_than_replaces(tmp_path):
    run = create_run({"m": 1}, runs_dir=tmp_path)
    record_status(run, state="running")
    record_status(run, returncode=0)
    status = run.status()
    assert status["state"] == "running" and status["returncode"] == 0


def test_iter_and_load_runs_round_trip(tmp_path):
    configs = [{"material": m} for m in ("NaI", "BGO", "CsI")]
    for config in configs:
        create_run(config, runs_dir=tmp_path)
    found = list(iter_runs(tmp_path))
    assert len(found) == 3
    for run in found:
        assert load_run(run.id, tmp_path).config == run.config


def test_iter_runs_tolerates_a_missing_directory(tmp_path):
    assert list(iter_runs(tmp_path / "nope")) == []


def test_checksum_detects_a_changed_byte(tmp_path):
    target = tmp_path / "out.csv"
    target.write_text("a")
    before = checksum(target)
    target.write_text("b")
    assert checksum(target) != before


# --------------------------------------------------------------------------- #
# Sweeps
# --------------------------------------------------------------------------- #

def test_sweep_expands_to_the_cartesian_product():
    spec = {
        "geometry": {"shape": "box"},
        "sweep": {"geometry.surface": ["polished", "ground"], "geometry.wrapping": ["teflon", "vm2000"]},
    }
    configs = expand_sweep(spec)
    assert len(configs) == 4
    assert {(c["geometry"]["surface"], c["geometry"]["wrapping"]) for c in configs} == {
        ("polished", "teflon"), ("polished", "vm2000"),
        ("ground", "teflon"), ("ground", "vm2000"),
    }


def test_sweep_output_carries_no_trace_of_the_sweep():
    """A swept configuration must hash identically to the same one written by hand."""
    swept = expand_sweep({"a": 1, "sweep": {"b": [2]}})[0]
    assert "sweep" not in swept
    assert run_id(swept) == run_id({"a": 1, "b": 2})


def test_sweep_without_axes_returns_the_base_configuration():
    assert expand_sweep({"a": 1}) == [{"a": 1}]


def test_sweep_does_not_mutate_shared_nested_structures():
    configs = expand_sweep({"geometry": {"shape": "box"}, "sweep": {"geometry.surface": ["p", "g"]}})
    assert configs[0]["geometry"] is not configs[1]["geometry"]
    assert configs[0]["geometry"]["shape"] == "box"


@pytest.mark.parametrize("bad", [[], "polished", None, {}])
def test_sweep_rejects_a_malformed_axis(bad):
    with pytest.raises(ValueError):
        expand_sweep({"sweep": {"x": bad}})


def test_prepare_runs_skips_configurations_already_executed(tmp_path):
    configs = [{"m": "NaI"}, {"m": "BGO"}]
    created, skipped = prepare_runs(configs, runs_dir=tmp_path)
    assert len(created) == 2 and skipped == []
    created, skipped = prepare_runs(configs + [{"m": "CsI"}], runs_dir=tmp_path)
    assert len(created) == 1 and len(skipped) == 2


# --------------------------------------------------------------------------- #
# Backends
# --------------------------------------------------------------------------- #

def test_local_backend_records_success(tmp_path):
    run = create_run({"m": 1}, runs_dir=tmp_path)
    result = LocalBackend().execute(run, ["python3", "-c", "print('hello')"])
    assert result.ok
    assert result.stdout_path.read_text().strip() == "hello"
    assert run.status()["state"] == "completed"


def test_local_backend_records_failure_without_raising(tmp_path):
    run = create_run({"m": 2}, runs_dir=tmp_path)
    result = LocalBackend().execute(run, ["python3", "-c", "import sys; sys.exit(3)"])
    assert not result.ok and result.returncode == 3
    assert run.status()["state"] == "failed"


def test_local_backend_runs_several_jobs_in_parallel(tmp_path):
    runs = [create_run({"m": i}, runs_dir=tmp_path) for i in range(4)]
    jobs = [(run, ["python3", "-c", f"print({i})"]) for i, run in enumerate(runs)]
    results = LocalBackend(max_workers=4).execute_many(jobs)
    assert all(r.ok for r in results)
    assert {r.stdout_path.read_text().strip() for r in results} == {"0", "1", "2", "3"}


def test_cluster_backends_fail_loudly_rather_than_silently(tmp_path):
    run = create_run({"m": 3}, runs_dir=tmp_path)
    with pytest.raises(NotImplementedError, match="not implemented yet"):
        SlurmBackend().execute(run, ["true"])


def test_local_backend_rejects_zero_workers():
    with pytest.raises(ValueError):
        LocalBackend(max_workers=0)
