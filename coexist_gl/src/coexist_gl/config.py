from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class FamilyConfig:
    path: Path
    family_id: str
    species_tree: Path
    gene_tree: Path
    genes: Path
    gene_tree_rootedness: str
    permutations: int
    clustering_seed: int
    fit_starts: int
    fit_seed: int
    r_values: tuple[float, ...]
    A_values: tuple[float, ...] | None
    species_cut_k: tuple[int, ...]


def _resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def load_config(path: str | Path) -> FamilyConfig:
    config_path = Path(path).resolve()
    raw: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if raw.get("schema_version") != "coexist_gl_family_v2":
        raise ValueError("schema_version must be coexist_gl_family_v2")
    if raw.get("species_tree", {}).get("edge_policy") != "unit":
        raise ValueError("v2 requires species_tree.edge_policy=unit")
    if raw.get("species_tree", {}).get("rooted") is not True:
        raise ValueError("v2 requires species_tree.rooted=true")
    if raw.get("copy_state", {}).get("single") != 1 or raw.get("copy_state", {}).get("coexist_min") != 2:
        raise ValueError("v2 copy states must be single=1 and coexist_min=2")
    if float(raw.get("fit", {}).get("r_lower_bound", 1.0)) != 1.0:
        raise ValueError("v2 requires r_lower_bound=1")

    surface_A = raw.get("surface", {}).get("A_values", "auto_around_mle")
    A_values = None if surface_A == "auto_around_mle" else tuple(float(x) for x in surface_A)
    r_values = tuple(float(x) for x in raw.get("surface", {}).get("r_values", [1, 2, 4, 8, 16, 32, 64]))
    if not r_values or any(x < 1 for x in r_values):
        raise ValueError("surface r_values must all be >=1")
    if A_values is not None and (not A_values or any(x <= 0 for x in A_values)):
        raise ValueError("surface A_values must all be >0")

    base = config_path.parent
    rootedness = str(raw["gene_tree"]["rootedness"])
    if rootedness not in {"rooted", "unrooted", "unknown"}:
        raise ValueError("gene_tree.rootedness must be rooted, unrooted, or unknown")
    return FamilyConfig(
        path=config_path,
        family_id=str(raw["family_id"]),
        species_tree=_resolve(base, raw["species_tree"]["path"]),
        gene_tree=_resolve(base, raw["gene_tree"]["path"]),
        genes=_resolve(base, raw["gene_metadata"]["path"]),
        gene_tree_rootedness=rootedness,
        permutations=int(raw.get("clustering_test", {}).get("permutations", 9999)),
        clustering_seed=int(raw.get("clustering_test", {}).get("seed", 1)),
        fit_starts=int(raw.get("fit", {}).get("starts", 6)),
        fit_seed=int(raw.get("fit", {}).get("seed", 1)),
        r_values=r_values,
        A_values=A_values,
        species_cut_k=tuple(int(x) for x in raw.get("species_cuts", {}).get("topology_k", [])),
    )
