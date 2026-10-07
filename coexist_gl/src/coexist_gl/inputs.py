from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo

from .config import FamilyConfig
from .utils import sha256


@dataclass
class ValidatedInputs:
    config: FamilyConfig
    species_tree: object
    gene_tree: object
    gene_to_host: dict[str, str]
    copy_counts: dict[str, int]
    states: dict[str, int]
    audit: dict[str, object]


def _read_gene_metadata(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "host_id", "include"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("genes.tsv requires gene_id, host_id, include")
        mapping: dict[str, str] = {}
        for row in reader:
            if str(row["include"]).strip() not in {"1", "true", "TRUE", "True"}:
                continue
            gene = str(row["gene_id"]).strip()
            host = str(row["host_id"]).strip()
            if not gene or not host:
                raise ValueError("empty included gene_id or host_id")
            if gene in mapping:
                raise ValueError(f"duplicate included gene_id: {gene}")
            mapping[gene] = host
    if not mapping:
        raise ValueError("no included genes")
    return mapping


def _assert_unique_tip_names(tree, label: str) -> list[str]:
    tips = [tip.name for tip in tree.get_terminals()]
    if any(not tip for tip in tips):
        raise ValueError(f"{label} contains an empty tip name")
    if len(tips) != len(set(tips)):
        raise ValueError(f"{label} contains duplicate tip names")
    return tips


def _assert_binary(tree, unrooted: bool, label: str) -> None:
    for clade in tree.find_clades():
        degree = len(clade.clades)
        if degree == 0:
            continue
        if clade is tree.root and unrooted and degree == 3:
            continue
        if degree != 2:
            raise ValueError(f"{label} is not {'unrooted-' if unrooted else 'rooted-'}binary")


def validate_inputs(config: FamilyConfig) -> ValidatedInputs:
    for path in (config.species_tree, config.gene_tree, config.genes):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    species_tree = Phylo.read(str(config.species_tree), "newick")
    gene_tree = Phylo.read(str(config.gene_tree), "newick")
    gene_to_host = _read_gene_metadata(config.genes)

    species_tips = _assert_unique_tip_names(species_tree, "species tree")
    gene_tips = _assert_unique_tip_names(gene_tree, "gene tree")
    _assert_binary(species_tree, unrooted=False, label="species tree")
    _assert_binary(gene_tree, unrooted=config.gene_tree_rootedness == "unrooted", label="gene tree")

    if set(gene_tips) != set(gene_to_host):
        raise ValueError(
            f"gene tip/mapping mismatch: tree_only={len(set(gene_tips)-set(gene_to_host))} "
            f"mapping_only={len(set(gene_to_host)-set(gene_tips))}"
        )
    mapped_hosts = set(gene_to_host.values())
    if set(species_tips) != mapped_hosts:
        raise ValueError(
            f"species tip/mapped-host mismatch: tree_only={len(set(species_tips)-mapped_hosts)} "
            f"mapping_only={len(mapped_hosts-set(species_tips))}"
        )

    copy_counts = dict(Counter(gene_to_host.values()))
    states = {host: 0 if count == 1 else 1 for host, count in copy_counts.items()}
    audit = {
        "status": "PASS",
        "family_id": config.family_id,
        "species_tips": len(species_tips),
        "gene_tips": len(gene_tips),
        "single_hosts": sum(state == 0 for state in states.values()),
        "coexist_hosts": sum(state == 1 for state in states.values()),
        "species_tree_sha256": sha256(config.species_tree),
        "gene_tree_sha256": sha256(config.gene_tree),
        "genes_sha256": sha256(config.genes),
        "edge_policy": "unit",
    }
    return ValidatedInputs(config, species_tree, gene_tree, gene_to_host, copy_counts, states, audit)


def applicability(validated: ValidatedInputs) -> dict[str, object]:
    n_single = sum(state == 0 for state in validated.states.values())
    n_coexist = sum(state == 1 for state in validated.states.values())
    if n_coexist == 0:
        status = "NOT_APPLICABLE_NO_COEXIST"
    elif n_single == 0:
        status = "NOT_IDENTIFIABLE_NO_STATE_CONTRAST"
    else:
        status = "FIT_ELIGIBLE"
    return {
        "status": status,
        "n_hosts": n_single + n_coexist,
        "n_single": n_single,
        "n_coexist": n_coexist,
        "coexist_fraction": n_coexist / (n_single + n_coexist),
        "coexist_observed": n_coexist > 0,
    }

