"""Smoke tests: confirm the package imports, configs parse, external data exist,
and the ``btmw`` CLI is wired. No raw simulation data is touched.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_package_importable() -> None:
    import btmw

    assert btmw.__version__


def test_paths_module() -> None:
    from btmw import paths

    assert paths.repo_root() == REPO_ROOT
    assert paths.configs_dir().is_dir()
    assert paths.data_dir().is_dir()
    assert paths.cache_dir().is_dir()
    assert paths.external_dir().is_dir()


def test_load_simulations_yields_nine_sims() -> None:
    from btmw.config import load_simulations

    cfg = load_simulations()
    assert cfg.working_root  # non-empty path
    assert len(cfg.simulations) == 9
    # 6 fiducial + 3 high-res
    fiducial = [s for s in cfg.simulations.values() if s.resolution == "fiducial"]
    high     = [s for s in cfg.simulations.values() if s.resolution == "high"]
    assert len(fiducial) == 6
    assert len(high) == 3
    # has_vr only on the 6 fiducial sims
    assert all(s.has_vr for s in fiducial)
    assert not any(s.has_vr for s in high)


def test_simulation_filter() -> None:
    from btmw.config import load_simulations

    cfg = load_simulations()
    m12i_fiducial = cfg.filter(host="m12i", resolution="fiducial")
    assert len(m12i_fiducial) == 3
    assert {s.model for s in m12i_fiducial} == {"PL", "BT_deep", "BT_soft"}


def test_sim_root_env_override(monkeypatch: pytest.MonkeyPatch) -> None:
    from btmw.config import load_simulations

    load_simulations.cache_clear()
    monkeypatch.setenv("BTMW_SIM_ROOT", "/tmp/btmw-data")
    try:
        cfg = load_simulations()
        assert cfg.working_root == Path("/tmp/btmw-data")
    finally:
        load_simulations.cache_clear()


@pytest.mark.parametrize(
    "filename",
    [
        "hellwing2016_smf_r50c.txt",
        "cautun2014_svf_r100c.txt",
        "springel2008_nprofile_r200c.txt",
        "lovell2014_cmf_r200c.txt",
        "grand2021_rmax_r200c.txt",
    ],
)
def test_active_external_reference_data_present(filename: str) -> None:
    from btmw import paths

    p = paths.external_dir() / filename
    assert p.is_file(), f"missing external reference file: {p}"
    # All curves are 2-column whitespace tables.
    rows = [ln for ln in p.read_text().splitlines() if ln.strip()]
    assert len(rows) > 1, f"empty reference file: {p}"
    cols = rows[0].split()
    assert len(cols) >= 2, f"expected >=2 columns, got {cols!r} in {p}"


def test_cli_help_runs() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "btmw.cli", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "btmw" in proc.stdout


def test_cli_list_sims_runs() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "btmw.cli", "list-sims"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "m12i_cdmo" in proc.stdout
    assert "m12f_cdmo" in proc.stdout


def test_figure_subcommands_are_wired() -> None:
    """Figure subcommands parse successfully without touching raw data."""
    from btmw.cli import build_parser

    parser = build_parser()
    commands = [
        ["hmf"],
        ["hmf", "--compare-vr"],
        ["hmf", "--resolution-study"],
        ["hvf"],
        ["hrf", "--bin", "6"],
        ["cmf"],
        ["rvsv", "--host", "m12i"],
        ["mvsv", "--host", "m12i"],
    ]
    for command in commands:
        args = parser.parse_args(command)
        assert callable(args.func)
