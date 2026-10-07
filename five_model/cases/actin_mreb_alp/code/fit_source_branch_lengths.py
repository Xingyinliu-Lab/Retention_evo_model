#!/usr/bin/env python3
"""Fit the frozen MreB/ALP five models using source species-tree branch lengths."""

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
from scipy.linalg import expm
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
LEVELS = ("M-only", "M+ALP coexist", "ALP-only")
CODE = {state: index for index, state in enumerate(LEVELS)}
MODELS = [
    ("ARD_unrestricted", "stationary", tuple((a, b) for a in range(3) for b in range(3) if a != b)),
    ("adjacency_coexist_middle", "stationary", ((0, 1), (1, 0), (1, 2), (2, 1))),
    ("adjacency_ALP_middle", "stationary", ((0, 2), (2, 0), (2, 1), (1, 2))),
    ("adjacency_M_middle", "stationary", ((2, 0), (0, 2), (0, 1), (1, 0))),
    ("ancestral_coexistence_strict", 1, ((1, 0), (1, 2))),
]


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
    children = [[(index[id(child)], max(float(child.branch_length or 0.0), 1e-8))
                 for child in node.clades] for node in nodes]
    names = [node.name if node.is_terminal() else None for node in nodes]

    def objective(log_rates):
        q = np.zeros((3, 3))
        for value, (source, target) in zip(np.exp(log_rates), edges):
            q[source, target] = value
        np.fill_diagonal(q, -q.sum(axis=1))
        likelihoods, cache, log_scale = [np.ones(3) for _ in nodes], {}, 0.0
        for i, node in enumerate(nodes):
            if node.is_terminal():
                vector = np.zeros(3); vector[states[names[i]]] = 1.0; likelihoods[i] = vector
                continue
            vector = np.ones(3)
            for child_i, length in children[i]:
                key = round(length, 12)
                if key not in cache:
                    cache[key] = expm(q * length)
                vector *= cache[key] @ likelihoods[child_i]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                return 1e100
            likelihoods[i] = vector / total; log_scale += math.log(total)
        prior = stationary(q) if root_mode == "stationary" else np.eye(3)[root_mode]
        value = float(prior @ likelihoods[-1])
        return -(math.log(value) + log_scale) if value > 0 and np.isfinite(value) else 1e100

    k = len(edges)
    starts = [np.zeros(k), np.full(k, -1.0), np.full(k, 1.0), np.linspace(-2, 2, k), np.linspace(2, -2, k)]
    fits = [minimize(objective, start, method="L-BFGS-B", bounds=[(-10, 10)] * k,
                     options={"maxiter": 400, "ftol": 1e-11, "gtol": 1e-7}) for start in starts]
    best = min(fits, key=lambda result: result.fun)
    counts = Counter(item["states"].values())
    return {
        "tree": item["name"], "model": model,
        "root_mode": LEVELS[root_mode] if isinstance(root_mode, int) else root_mode,
        "branch_scale": "source_branch_lengths", "k": k,
        "log_likelihood": -float(best.fun), "AIC": 2 * k + 2 * float(best.fun),
        "success": str(bool(best.success)).upper(),
        "gradient_max": float(np.max(np.abs(best.jac))) if best.jac is not None else "",
        "rates": ";".join(f"{LEVELS[a]}->{LEVELS[b]}:{math.exp(value):.10g}"
                          for value, (a, b) in zip(best.x, edges)),
        "mapped_MAGs": len(states), "n_M_only": counts["M-only"],
        "n_coexist": counts["M+ALP coexist"], "n_ALP_only": counts["ALP-only"],
    }


def main():
    payload = json.loads((ROOT / "input/six_tree_states.json").read_text(encoding="utf-8"))
    tasks = [(item, *model) for item in payload["trees"] for model in MODELS]
    with mp.Pool(processes=len(tasks)) as pool:
        rows = pool.map(fit_one, tasks)
    for tree in {row["tree"] for row in rows}:
        block = [row for row in rows if row["tree"] == tree]
        minimum = min(row["AIC"] for row in block)
        for row in block:
            row["delta_AIC"] = row["AIC"] - minimum
    rows.sort(key=lambda row: (row["tree"], row["model"]))
    out = ROOT / "results" / "source_branch_lengths_20260820"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "MreB_ALP_five_models_source_branch_lengths.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (out / "COMPLETE.marker").write_text("status=complete\ntrees=6\nmodels=5\nbranch_scale=source_branch_lengths\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "rows": len(rows)}))


if __name__ == "__main__":
    main()
