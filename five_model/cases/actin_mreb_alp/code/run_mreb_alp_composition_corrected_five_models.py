#!/usr/bin/env python3
"""Fit five composition-corrected MreB/ALP state models on six species trees."""

from __future__ import annotations

import csv
import io
import json
import math
import multiprocessing as mp
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from scipy.linalg import expm
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parents[1]
RETENTION = HERE.parents[1]
ACTIN = RETENTION.parent
STATE_FILE = RETENTION / "mreb_cd_aggregate_retention_lowcost_audit_20260810/data/MAG_MreB_CD_aggregate_states.tsv"
MAP_FILE = RETENTION / "mreb_cd_aggregate_retention_lowcost_audit_20260810/results/gtdb_exact_mapping.tsv"
GTDB_TREE = ACTIN / "actin_origin_species_tree/plot_source_data/gtdb_ar53_expanded_with_user_mags.treefile"
GUIDE_TREE = RETENTION / "actin_retention_exemplar_search_20260727/inputs/species_trees/gtdb_expanded_18261_with_asgard971.tree"
TREE_PATHS = {
    "A1338_addcore": RETENTION / "actin_retention_exemplar_search_20260727/inputs/species_trees/archaea1338/archaea1338_orthofinder_species_tree_rooted.txt",
    "A1338_denovo": RETENTION / "actin_retention_exemplar_search_20260727/inputs/species_trees/archaea1338/archaea1338_denovo_species_tree_unrooted.txt",
    "A1338_ASTRAL64": RETENTION / "archaea1338_full_og_species_marker_screen_20260801/astral_pro_third_species_tree_20260805/archaea1338_64_formal_astral_pro3.tree",
    "Asgard971_addcore": RETENTION / "actin_retention_exemplar_search_20260727/inputs/species_trees/asgard971_orthofinder_species_tree_rooted.txt",
    "Asgard971_ASTRAL89": RETENTION / "asgard971_full_og_species_marker_screen_20260728/formal_gene_trees_richLG_automated1_20260728/astral_pro_third_species_tree_20260730/asgard971_89_formal_astral_pro3.tree",
}
LEVELS = ("M-only", "M+ALP coexist", "ALP-only")
CODE = {state: index for index, state in enumerate(LEVELS)}
SOURCE_TO_FORMAL = {"M-state": "M-only", "coexistence": "M+ALP coexist", "C/D-state": "ALP-only"}
MODELS = [
    ("ARD_unrestricted", "stationary", tuple((a, b) for a in range(3) for b in range(3) if a != b)),
    ("adjacency_coexist_middle", "stationary", ((0, 1), (1, 0), (1, 2), (2, 1))),
    ("adjacency_ALP_middle", "stationary", ((0, 2), (2, 0), (2, 1), (1, 2))),
    ("adjacency_M_middle", "stationary", ((2, 0), (0, 2), (0, 1), (1, 0))),
    ("ancestral_coexistence_strict", 1, ((1, 0), (1, 2))),
]


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def restrict(source, keep):
    def visit(node, root=False):
        if node.is_terminal():
            return Clade(branch_length=node.branch_length, name=node.name) if node.name in keep else None
        children = [value for child in node.clades if (value := visit(child)) is not None]
        if not children:
            return None
        if len(children) == 1 and not root:
            child = children[0]
            child.branch_length = float(child.branch_length or 0) + float(node.branch_length or 0)
            return child
        return Clade(branch_length=node.branch_length, clades=children)

    root = visit(source.root, True)
    if root is None:
        raise RuntimeError("Pruning removed every tip")
    while len(root.clades) == 1 and not root.is_terminal():
        root = root.clades[0]
        root.branch_length = None
    return Tree(root=root, rooted=True)


def leaves(tree_or_clade):
    target = tree_or_clade.root if hasattr(tree_or_clade, "root") else tree_or_clade
    return {tip.name for tip in target.get_terminals()}


def project_operational_root(target, guide):
    universe = leaves(target)
    guide_sides = [leaves(child) & universe for child in guide.root.clades]
    guide_sides = [side for side in guide_sides if side and side != universe]
    guide_side = min(guide_sides, key=len)
    best = None
    for clade in target.find_clades(order="preorder"):
        if clade is target.root:
            continue
        candidate = leaves(clade)
        if not candidate or candidate == universe:
            continue
        score = min(
            len(candidate.symmetric_difference(guide_side)),
            len((universe - candidate).symmetric_difference(guide_side)),
        )
        if best is None or score < best[0]:
            best = (score, clade)
    if best is None:
        raise RuntimeError("No edge available for GTDB root projection")
    target.root_with_outgroup(best[1])
    target.rooted = True
    return best[0]


def to_newick(tree):
    handle = io.StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue().strip()


def prepare_input():
    states = {
        row["mag_id"]: SOURCE_TO_FORMAL[row["state"]]
        for row in read_tsv(STATE_FILE)
        if row["state"] in SOURCE_TO_FORMAL
    }
    guide_source = Phylo.read(str(GUIDE_TREE), "newick")
    guide_tips = leaves(guide_source)
    items = []

    gtdb_source = Phylo.read(str(GTDB_TREE), "newick")
    mapped = {
        row["gtdb_tip"]: states[row["mag_id"]]
        for row in read_tsv(MAP_FILE)
        if row.get("gtdb_tip") and row["mag_id"] in states
    }
    mapped = {tip: state for tip, state in mapped.items() if tip in leaves(gtdb_source)}
    gtdb = restrict(gtdb_source, mapped)
    items.append({"name": "GTDB_r226", "newick": to_newick(gtdb), "states": mapped,
                  "root_projection_mismatch": 0, "source_tree": str(GTDB_TREE)})

    for name, path in TREE_PATHS.items():
        source = Phylo.read(str(path), "newick")
        shared = leaves(source) & set(states) & guide_tips
        mapped = {tip: states[tip] for tip in shared}
        target = restrict(source, mapped)
        guide = restrict(guide_source, mapped)
        mismatch = project_operational_root(target, guide)
        items.append({"name": name, "newick": to_newick(target), "states": mapped,
                      "root_projection_mismatch": mismatch, "source_tree": str(path)})

    for item in items:
        parsed = Phylo.read(io.StringIO(item["newick"]), "newick")
        tips = leaves(parsed)
        if tips != set(item["states"]):
            raise RuntimeError(f"Tip/state mismatch for {item['name']}")
        counts = Counter(item["states"].values())
        item["n_tips"] = len(tips)
        item["state_counts"] = {state: counts[state] for state in LEVELS}
        item["root_semantics"] = "gtdb_projected"
    payload = {"schema_version": "mreb_alp_composition_corrected_five_model_v1",
               "branch_scale": "unit_topology", "trees": items}
    HERE.joinpath("input").mkdir(parents=True, exist_ok=True)
    HERE.joinpath("input/six_tree_states.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def stationary(q):
    values, vectors = np.linalg.eig(q.T)
    vector = np.abs(np.real(vectors[:, np.argmin(np.abs(values))]))
    return vector / vector.sum()


def fit_one(task):
    item, model, root_mode, edges = task
    tree = Phylo.read(io.StringIO(item["newick"]), "newick")
    states = {tip: CODE[state] for tip, state in item["states"].items()}
    nodes = list(tree.find_clades(order="postorder"))
    index = {id(node): i for i, node in enumerate(nodes)}
    children = [[index[id(child)] for child in node.clades] for node in nodes]
    names = [node.name if node.is_terminal() else None for node in nodes]

    def objective(log_rates):
        q = np.zeros((3, 3))
        for value, (source, target) in zip(np.exp(log_rates), edges):
            q[source, target] = value
        np.fill_diagonal(q, -q.sum(axis=1))
        transition = expm(q)
        likelihoods = [np.ones(3) for _ in nodes]
        log_scale = 0.0
        for i, node in enumerate(nodes):
            if node.is_terminal():
                vector = np.zeros(3)
                vector[states[names[i]]] = 1.0
                likelihoods[i] = vector
                continue
            vector = np.ones(3)
            for child_i in children[i]:
                vector *= transition @ likelihoods[child_i]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                return 1e100
            likelihoods[i] = vector / total
            log_scale += math.log(total)
        prior = stationary(q) if root_mode == "stationary" else np.eye(3)[root_mode]
        value = float(prior @ likelihoods[-1])
        return -(math.log(value) + log_scale) if value > 0 and np.isfinite(value) else 1e100

    k = len(edges)
    starts = [np.zeros(k), np.full(k, -1.0), np.full(k, 1.0), np.linspace(-2, 2, k), np.linspace(2, -2, k)]
    fits = [minimize(objective, start, method="L-BFGS-B", bounds=[(-10, 10)] * k,
                     options={"maxiter": 400, "ftol": 1e-11, "gtol": 1e-7}) for start in starts]
    converged = [result for result in fits if result.success and np.isfinite(result.fun)]
    best = min(converged or fits, key=lambda result: result.fun)
    counts = Counter(item["states"].values())
    return {
        "tree": item["name"], "model": model,
        "root_mode": LEVELS[root_mode] if isinstance(root_mode, int) else root_mode,
        "branch_scale": "unit_topology", "k": k,
        "log_likelihood": -float(best.fun), "AIC": 2 * k + 2 * float(best.fun),
        "success": str(bool(best.success)).upper(),
        "gradient_max": float(np.max(np.abs(best.jac))) if best.jac is not None else "",
        "rates": ";".join(f"{LEVELS[a]}->{LEVELS[b]}:{math.exp(value):.10g}" for value, (a, b) in zip(best.x, edges)),
        "mapped_MAGs": len(states), "root_projection_mismatch": item["root_projection_mismatch"],
        "n_M_only": counts["M-only"], "n_coexist": counts["M+ALP coexist"],
        "n_ALP_only": counts["ALP-only"],
    }


def main():
    payload = prepare_input()
    tasks = [(item, *model) for item in payload["trees"] for model in MODELS]
    with mp.Pool(processes=len(tasks)) as pool:
        rows = pool.map(fit_one, tasks)
    for tree in {row["tree"] for row in rows}:
        block = [row for row in rows if row["tree"] == tree]
        minimum = min(row["AIC"] for row in block)
        for row in block:
            row["delta_AIC"] = row["AIC"] - minimum
    out = HERE / "results"
    out.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (out / "MreB_ALP_composition_corrected_five_models_six_trees.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (HERE / "FIT_COMPLETE.marker").write_text("status=complete\nmodels=5\ntrees=6\nbranch_scale=unit_topology\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "rows": len(rows), "trees": len(payload["trees"])}))


if __name__ == "__main__":
    main()
