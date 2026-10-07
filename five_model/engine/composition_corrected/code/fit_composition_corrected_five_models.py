#!/usr/bin/env python3
"""Fit the final composition-corrected five-model CTMC comparison.

This is the portable, supplement-facing implementation. It intentionally
contains only the five final models and accepts only the frozen six-tree input.
"""

import argparse
import csv
import json
import math
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize


STATE = {"dALP-only": 0, "coexist": 1, "cALP-only": 2}
LEVEL = ["dALP-only", "coexist", "cALP-only"]

# The three four-rate models differ only in which state is the middle of the
# reversible adjacency graph. Their root prior is the stationary distribution
# implied by the fitted Q matrix. The strict ancestral-coexistence model fixes
# the root to coexist by definition. ARD allows all six directed transitions.
MODELS = [
    ("adjacency_coexist_middle", "stationary", ((0, 1), (1, 0), (1, 2), (2, 1))),
    ("adjacency_cALP_middle", "stationary", ((0, 2), (2, 0), (2, 1), (1, 2))),
    ("adjacency_dALP_middle", "stationary", ((2, 0), (0, 2), (0, 1), (1, 0))),
    ("ancestral_coexist_divergence", 1, ((1, 0), (1, 2))),
    ("ARD_unrestricted", "stationary", tuple((a, b) for a in range(3) for b in range(3) if a != b)),
]


class Node:
    def __init__(self, name="", branch_length=0.0):
        self.name = name
        self.branch_length = branch_length
        self.children = []

    def terminal(self):
        return not self.children


def parse_newick(text):
    text = text.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if pos < len(text) and text[pos] == ",":
                    pos += 1
                    continue
                if pos >= len(text) or text[pos] != ")":
                    raise ValueError("Malformed Newick: missing closing parenthesis")
                pos += 1
                break
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            pos += 1
        token = text[start:pos].strip()
        if ":" in token:
            node.name, length = token.rsplit(":", 1)
            node.branch_length = float(length or 0.0)
        else:
            node.name = token
        return node

    root = parse_node()
    if pos != len(text):
        raise ValueError(f"Malformed Newick: parser stopped at {pos}/{len(text)}")
    return root


def tree_tip_names(root):
    tips = []

    def walk(node):
        if node.terminal():
            tips.append(node.name)
        else:
            if len(node.children) < 2:
                raise ValueError("Unary internal nodes are not permitted")
            for child in node.children:
                walk(child)

    walk(root)
    if len(tips) != len(set(tips)):
        raise ValueError("Duplicate species-tree tip labels")
    return tips


def validate_input(payload):
    if payload.get("schema_version") != "composition_corrected_five_model_input_v1":
        raise ValueError("Unexpected input schema_version")
    trees = payload.get("trees", [])
    if len(trees) != 6:
        raise ValueError(f"Expected exactly six final trees; observed {len(trees)}")
    names = [item["name"] for item in trees]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate tree names")
    for item in trees:
        if item.get("root_semantics") != "gtdb_projected":
            raise ValueError(f"Non-final root semantics in {item['name']}")
        states = item.get("states", {})
        if set(states.values()) - set(STATE):
            raise ValueError(f"Unexpected state label in {item['name']}")
        tips = tree_tip_names(parse_newick(item["newick"]))
        if set(tips) != set(states):
            missing = sorted(set(tips) - set(states))[:5]
            extra = sorted(set(states) - set(tips))[:5]
            raise ValueError(f"Tip/state mismatch in {item['name']}: missing={missing}, extra={extra}")
    return trees


def stationary(q):
    values, vectors = np.linalg.eig(q.T)
    vector = np.real(vectors[:, np.argmin(np.abs(values))])
    vector = np.abs(vector)
    return vector / vector.sum()


def fit(task):
    item, scale, model, root_mode, transitions = task
    tree = parse_newick(item["newick"])
    tip_states = {name: STATE[state] for name, state in item["states"].items()}
    nodes = []

    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)

    walk(tree)
    index = {id(node): i for i, node in enumerate(nodes)}
    children = [
        [
            (index[id(child)], 1.0 if scale == "unit_topology" else max(child.branch_length, 1e-8))
            for child in node.children
        ]
        for node in nodes
    ]
    names = [node.name if node.terminal() else None for node in nodes]

    def objective(log_rates):
        q = np.zeros((3, 3))
        for rate, (source, target) in zip(np.exp(log_rates), transitions):
            q[source, target] = rate
        np.fill_diagonal(q, -q.sum(axis=1))
        likelihoods = [np.ones(3) for _ in nodes]
        cache = {}
        log_scale = 0.0
        for i, node in enumerate(nodes):
            if node.terminal():
                vector = np.zeros(3)
                vector[tip_states[names[i]]] = 1.0
                likelihoods[i] = vector
                continue
            vector = np.ones(3)
            for child_i, branch_length in children[i]:
                key = round(branch_length, 12)
                if key not in cache:
                    cache[key] = expm(q * branch_length)
                vector *= cache[key] @ likelihoods[child_i]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                return 1e100
            likelihoods[i] = vector / total
            log_scale += math.log(total)
        root_prior = stationary(q) if root_mode == "stationary" else np.eye(3)[root_mode]
        root_like = float(root_prior @ likelihoods[-1])
        if root_like <= 0 or not np.isfinite(root_like):
            return 1e100
        return -(math.log(root_like) + log_scale)

    k = len(transitions)
    starts = [
        np.zeros(k), np.full(k, -1.0), np.full(k, 1.0),
        np.linspace(-2, 2, k), np.linspace(2, -2, k),
    ]
    fits = [
        minimize(
            objective, start, method="L-BFGS-B", bounds=[(-10, 10)] * k,
            options={"maxiter": 400, "ftol": 1e-11, "gtol": 1e-7},
        )
        for start in starts
    ]
    best = min(fits, key=lambda result: result.fun)
    return {
        "tree": item["name"],
        "branch_scale": scale,
        "model": model,
        "root_mode": LEVEL[root_mode] if isinstance(root_mode, int) else root_mode,
        "k": k,
        "log_likelihood": -float(best.fun),
        "AIC": 2 * k + 2 * float(best.fun),
        "success": str(bool(best.success)).upper(),
        "gradient_max": float(np.max(np.abs(best.jac))) if best.jac is not None else "",
        "n_starts": len(starts),
        "rates": ";".join(
            f"{LEVEL[source]}->{LEVEL[target]}:{math.exp(value):.10g}"
            for value, (source, target) in zip(best.x, transitions)
        ),
        "n_tips": len(tip_states),
        "dALP_only": sum(value == 0 for value in tip_states.values()),
        "coexist": sum(value == 1 for value in tip_states.values()),
        "cALP_only": sum(value == 2 for value in tip_states.values()),
    }


def main():
    package_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=package_root / "input/six_tree_tip_states.json")
    parser.add_argument("--output", type=Path, default=package_root / "output/five_model_fits.tsv")
    parser.add_argument("--marker", type=Path, default=package_root / "output/FIVE_MODEL_FIT_COMPLETE.marker")
    parser.add_argument("--processes", type=int, default=24)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    trees = validate_input(payload)
    tasks = [
        (item, scale, *model)
        for item in trees
        for scale in ("unit_topology", "source_branch_lengths")
        for model in MODELS
    ]
    with mp.Pool(processes=min(max(args.processes, 1), len(tasks))) as pool:
        rows = pool.map(fit, tasks)
    for key in {(row["tree"], row["branch_scale"]) for row in rows}:
        block = [row for row in rows if (row["tree"], row["branch_scale"]) == key]
        minimum = min(row["AIC"] for row in block)
        for row in block:
            row["delta_AIC"] = row["AIC"] - minimum

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    args.marker.write_text(
        "status=complete\nroot_prior=stationary_for_adjacency_and_ARD\nmodels=5\n"
        f"trees={len(trees)}\nbranch_scales=2\nrows={len(rows)}\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "complete", "trees": len(trees), "rows": len(rows)}))


if __name__ == "__main__":
    main()
