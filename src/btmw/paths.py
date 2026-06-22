"""Repo-relative path resolution.

The repo root is detected by walking up from this file until we find a
``pyproject.toml``. This works whether the package is installed editable
(``pip install -e .``) or in-place from the source tree.
"""

from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    """Return the absolute path to the BlueTilt-MWSubhalos repo root."""
    here = Path(__file__).resolve()
    for parent in [here, *here.parents]:
        if (parent / "pyproject.toml").is_file():
            return parent
    raise RuntimeError(
        f"Could not locate repo root (no pyproject.toml above {here})."
    )


def configs_dir() -> Path:
    return repo_root() / "configs"


def data_dir() -> Path:
    return repo_root() / "data"


def cache_dir() -> Path:
    return data_dir() / "cache"


def external_dir() -> Path:
    return data_dir() / "external"


def figures_dir() -> Path:
    return repo_root() / "figures"


def figures_static_dir() -> Path:
    return repo_root() / "figures_static"


def logs_dir() -> Path:
    return repo_root() / "logs"
