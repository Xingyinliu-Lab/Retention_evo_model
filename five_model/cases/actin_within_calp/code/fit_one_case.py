#!/usr/bin/env python3
"""Fit the fixed composition-corrected five-state-graph comparison for one case."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

LEVELS = ("A-only", "A+B coexist", "B-only")
CODE = {state: index for index, state in enumerate(LEVELS)}
MODELS = [
    ("ARD_unrestricted", "stationary", tuple((a, b) for a in range(3) for b in range(3) if a != b)),
    ("adjacency_coexist_middle", "stationary", ((0, 1), (1, 0), (1, 2), (2, 1))),
    ("adjacency_B_middle", "stationary", ((0, 2), (2, 0), (2, 1), (1, 2))),
    ("adjacency_A_middle", "stationary", ((2, 0), (0, 2), (0, 1), (1, 0))),
    ("ancestral_coexistence_strict", 1, ((1, 0), (1, 2))),
]


class Node:
    def __init__(self):
        self.name, self.children = "", []
    def terminal(self):
        return not self.children


def parse_newick(text):
    text = text.strip().rstrip(";")
    pos = 0
    def parse():
        nonlocal pos
        node = Node()
        if text[pos:pos + 1] == "(":
            pos += 1
            while True:
                node.children.append(parse())
                if text[pos:pos + 1] == ",":
                    pos += 1
                    continue
                if text[pos:pos + 1] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            pos += 1
        token = text[start:pos].strip()
        node.name = token.rsplit(":", 1)[0] if ":" in token else token
        return node
    return parse()


def stationary(q):
    values, vectors = np.linalg.eig(q.T)
    vector = np.abs(np.real(vectors[:, np.argmin(np.abs(values))]))
    return vector / vector.sum()


def main():
    case = Path(sys.argv[1]).resolve()
    freeze = json.loads((case / "input/INPUT_FREEZE.json").read_text(encoding="utf-8"))
    states = {}
    with (case / "input/host_states.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            states[row["host_id"]] = CODE[row["state"]]
    root = parse_newick((case / "input/species_tree.nwk").read_text(encoding="utf-8"))
    nodes = []
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(root)
    tips = {node.name for node in nodes if node.terminal()}
    if tips != set(states):
        raise RuntimeError("species-tree tips != host-state keys")
    index = {id(node): i for i, node in enumerate(nodes)}
    children = [[index[id(child)] for child in node.children] for node in nodes]
    names = [node.name if node.terminal() else None for node in nodes]
    rows = []

    for model, root_mode, edges in MODELS:
        def objective(log_rates):
            q = np.zeros((3, 3))
            for value, (source, target) in zip(np.exp(log_rates), edges):
                q[source, target] = value
            np.fill_diagonal(q, -q.sum(axis=1))
            transition = expm(q)
            likelihoods = [np.ones(3) for _ in nodes]
            log_scale = 0.0
            for i, node in enumerate(nodes):
                if node.terminal():
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
        converged = [fit for fit in fits if fit.success and np.isfinite(fit.fun)]
        best = min(converged or fits, key=lambda fit: fit.fun)
        rows.append({
            "case_id": freeze["case_id"], "universe": freeze["universe"], "tree": "GTDB_r226",
            "model": model, "root_mode": LEVELS[root_mode] if isinstance(root_mode, int) else root_mode,
            "branch_scale": "unit_topology", "k": k, "log_likelihood": -float(best.fun),
            "AIC": 2 * k + 2 * float(best.fun), "success": str(bool(best.success)).upper(),
            "gradient_max": float(np.max(np.abs(best.jac))) if best.jac is not None else "",
            "rates": ";".join(f"{LEVELS[a]}->{LEVELS[b]}:{math.exp(x):.10g}" for x, (a, b) in zip(best.x, edges)),
            "n_hosts": len(states), "n_A_only": sum(x == 0 for x in states.values()),
            "n_coexist": sum(x == 1 for x in states.values()), "n_B_only": sum(x == 2 for x in states.values()),
        })
    minimum = min(row["AIC"] for row in rows)
    for row in rows:
        row["delta_AIC"] = row["AIC"] - minimum
    results = case / "results"
    results.mkdir(exist_ok=True)
    output = results / "five_model_fit_GTDB.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 5 or any(row["success"] != "TRUE" or not math.isfinite(row["AIC"]) for row in rows):
        raise RuntimeError("five-model completion audit failed")
    (case / "COMPLETE.marker").write_text("status=complete\ntree=GTDB_r226\nmodels=5\nbranch_scale=unit_topology\n", encoding="utf-8")
    print(json.dumps({"case_id": freeze["case_id"], "status": "complete", "best_model": min(rows, key=lambda row: row["AIC"])["model"]}))


if __name__ == "__main__":
    main()
