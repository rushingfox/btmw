"""Loaders and accessors for ``configs/simulations.yaml``."""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

from .paths import configs_dir


@dataclass(frozen=True)
class Simulation:
    label: str
    host: str            # "m12i" or "m12f"
    model: str           # "PL", "BT_deep", or "BT_soft"
    resolution: str      # "fiducial" or "high"
    snapshot_num: int
    has_vr: bool
    working_root: Path | None
    path_templates: dict[str, str]

    def _format(self, template: str) -> str:
        return template.format(label=self.label, snap=self.snapshot_num)

    def _require_root(self) -> Path:
        if self.working_root is None:
            raise RuntimeError(
                "BTMW_SIM_ROOT is not set. "
                "Export it before running with --refresh: "
                "export BTMW_SIM_ROOT=/path/to/simulation_root"
            )
        return self.working_root

    def snapshot_path(self) -> Path:
        return self._require_root() / self._format(self.path_templates["snapshot"])

    def soap_hbt_path(self) -> Path:
        return self._require_root() / self._format(self.path_templates["soap_hbt"])

    def hbt_subsnap_dir(self) -> Path:
        return self._require_root() / self._format(self.path_templates["hbt_subsnap_dir"])

    def hbt_subsnap_glob(self) -> str:
        return self.path_templates["hbt_subsnap_glob"].format(snap=self.snapshot_num)

    def soap_vr_path(self) -> Path | None:
        if not self.has_vr:
            return None
        return self._require_root() / self._format(self.path_templates["soap_vr"])

    def vr_raw_dir(self) -> Path | None:
        if not self.has_vr:
            return None
        return self._require_root() / self._format(self.path_templates["vr_raw_dir"])

    def vr_raw_basename(self) -> str:
        return self.path_templates["vr_raw_basename"]


@dataclass(frozen=True)
class SimulationsConfig:
    working_root: Path | None
    path_templates: dict[str, str]
    simulations: dict[str, Simulation]
    particle_mass_in_swift_unit: dict[str, float | None]
    defaults: dict[str, Any]

    def __getitem__(self, label: str) -> Simulation:
        return self.simulations[label]

    def filter(
        self,
        *,
        host: str | None = None,
        resolution: str | None = None,
        has_vr: bool | None = None,
    ) -> list[Simulation]:
        out = []
        for sim in self.simulations.values():
            if host is not None and sim.host != host:
                continue
            if resolution is not None and sim.resolution != resolution:
                continue
            if has_vr is not None and sim.has_vr != has_vr:
                continue
            out.append(sim)
        return out


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r") as fh:
        return yaml.safe_load(fh)


@lru_cache(maxsize=1)
def load_simulations() -> SimulationsConfig:
    raw = _load_yaml(configs_dir() / "simulations.yaml")
    sim_root_env = os.environ.get("BTMW_SIM_ROOT")
    working_root = Path(sim_root_env).expanduser() if sim_root_env else None
    templates = raw["path_templates"]
    sims_raw = raw["simulations"]
    sims = {
        label: Simulation(
            label=label,
            host=meta["host"],
            model=meta["model"],
            resolution=meta["resolution"],
            snapshot_num=int(meta["snapshot_num"]),
            has_vr=bool(meta.get("has_vr", False)),
            working_root=working_root,
            path_templates=templates,
        )
        for label, meta in sims_raw.items()
    }
    return SimulationsConfig(
        working_root=working_root,
        path_templates=templates,
        simulations=sims,
        particle_mass_in_swift_unit=raw.get("particle_mass_in_swift_unit", {}),
        defaults=raw.get("defaults", {}),
    )


