#!/usr/bin/env python3
"""Fit source-branch-length five models and decompose ARD for fastmode cases."""

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

LEVELS = ("A-only", "A+B coexist", "B-only")
INDEX = {x: i for i, x in enumerate(LEVELS)}
MODELS = [
    ("ARD_unrestricted", "stationary", tuple((a, b) for a in range(3) for b in range(3) if a != b)),
    ("adjacency_coexist_middle", "stationary", ((0, 1), (1, 0), (1, 2), (2, 1))),
    ("adjacency_B_middle", "stationary", ((0, 2), (2, 0), (2, 1), (1, 2))),
    ("adjacency_A_middle", "stationary", ((2, 0), (0, 2), (0, 1), (1, 0))),
    ("ancestral_coexistence_strict", 1, ((1, 0), (1, 2))),
]
GROUPS = {
    "coexist_mediated": ((0, 1), (1, 0), (1, 2), (2, 1)),
    "direct_endpoint": ((0, 2), (2, 0)),
}


class Node:
    def __init__(self):
        self.name, self.length, self.children = "", 0.0, []

    def terminal(self):
        return not self.children


def parse_newick(text: str) -> Node:
    text, pos = text.strip().rstrip(";"), 0

    def parse() -> Node:
        nonlocal pos
        node = Node()
        if text[pos:pos + 1] == "(":
            pos += 1
            while True:
                node.children.append(parse())
                if text[pos:pos + 1] == ",":
                    pos += 1
                    continue
                if text[pos:pos + 1] != ")":
                    raise ValueError("Malformed Newick")
                pos += 1
                break
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            pos += 1
        token = text[start:pos].strip()
        if ":" in token:
            node.name, length = token.rsplit(":", 1)
            node.length = float(length or 0.0)
        else:
            node.name = token
        return node

    root = parse()
    if pos != len(text):
        raise ValueError(f"Newick parser stopped at {pos}/{len(text)}")
    return root


def stationary(q):
    values, vectors = np.linalg.eig(q.T)
    vector = np.abs(np.real(vectors[:, np.argmin(np.abs(values))]))
    return vector / vector.sum()


def load_case(case: Path):
    freeze = json.loads((case / "input/INPUT_FREEZE.json").read_text(encoding="utf-8"))
    states = {}
    with (case / "input/host_states.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            states[row["host_id"]] = INDEX[row["state"]]
    root = parse_newick((case / "input/species_tree.nwk").read_text(encoding="utf-8"))
    return freeze, states, root


class Engine:
    def __init__(self, root, states):
        self.nodes = []

        def walk(node):
            for child in node.children:
                walk(child)
            self.nodes.append(node)

        walk(root)
        idx = {id(node): i for i, node in enumerate(self.nodes)}
        self.children = [
            [(idx[id(child)], max(float(child.length), 1e-8)) for child in node.children]
            for node in self.nodes
        ]
        self.names = [node.name if node.terminal() else None for node in self.nodes]
        self.states = states
        tips = {name for name in self.names if name is not None}
        if tips != set(states):
            raise ValueError("species-tree tips != host-state keys")

    def log_likelihood(self, q, prior):
        likelihoods = [np.ones(3) for _ in self.nodes]
        cache, log_scale = {}, 0.0
        for i, node in enumerate(self.nodes):
            if node.terminal():
                vector = np.zeros(3)
                vector[self.states[self.names[i]]] = 1.0
                likelihoods[i] = vector
                continue
            vector = np.ones(3)
            for child_i, length in self.children[i]:
                key = round(length, 12)
                if key not in cache:
                    cache[key] = expm(q * length)
                vector *= cache[key] @ likelihoods[child_i]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                return -1e100
            likelihoods[i] = vector / total
            log_scale += math.log(total)
        value = float(prior @ likelihoods[-1])
        return math.log(value) + log_scale if value > 0 and np.isfinite(value) else -1e100


def fit_task(task):
    case_text, model, root_mode, edges = task
    case = Path(case_text)
    freeze, states, root = load_case(case)
    engine = Engine(root, states)

    def objective(log_rates):
        q = np.zeros((3, 3))
        for value, (source, target) in zip(np.exp(log_rates), edges):
            q[source, target] = value
        np.fill_diagonal(q, -q.sum(axis=1))
        prior = stationary(q) if root_mode == "stationary" else np.eye(3)[root_mode]
        return -engine.log_likelihood(q, prior)

    k = len(edges)
    starts = [np.zeros(k), np.full(k, -1.0), np.full(k, 1.0), np.linspace(-2, 2, k), np.linspace(2, -2, k)]
    fits = [minimize(objective, start, method="L-BFGS-B", bounds=[(-10, 10)] * k,
                     options={"maxiter": 400, "ftol": 1e-11, "gtol": 1e-7}) for start in starts]
    best = min(fits, key=lambda x: x.fun)
    return {
        "case_id": freeze["case_id"], "universe": freeze["universe"], "tree": "GTDB_r226",
        "model": model, "root_mode": LEVELS[root_mode] if isinstance(root_mode, int) else root_mode,
        "branch_scale": "source_branch_lengths", "k": k, "log_likelihood": -float(best.fun),
        "AIC": 2 * k + 2 * float(best.fun), "success": str(bool(best.success)).upper(),
        "gradient_max": float(np.max(np.abs(best.jac))) if best.jac is not None else "",
        "rates": ";".join(f"{LEVELS[a]}->{LEVELS[b]}:{math.exp(x):.10g}" for x, (a, b) in zip(best.x, edges)),
        "n_hosts": len(states), "n_A_only": sum(x == 0 for x in states.values()),
        "n_coexist": sum(x == 1 for x in states.values()), "n_B_only": sum(x == 2 for x in states.values()),
    }


def parse_q(text):
    q = np.zeros((3, 3))
    for token in text.split(";"):
        edge, value = token.rsplit(":", 1)
        source, target = edge.split("->", 1)
        q[INDEX[source], INDEX[target]] = float(value)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


def tilted(q, edges, theta):
    matrix, factor = q.copy(), math.exp(theta)
    for source, target in edges:
        matrix[source, target] = q[source, target] * factor
    np.fill_diagonal(matrix, np.diag(q))
    return matrix


def ard_task(task):
    case_text, ard = task
    case = Path(case_text)
    _, states, root = load_case(case)
    engine, q = Engine(root, states), parse_q(ard["rates"])
    prior = stationary(q)
    base = engine.log_likelihood(q, prior)
    values = {}
    h = 5e-4
    for name, edges in GROUPS.items():
        delta = {m: engine.log_likelihood(tilted(q, edges, m * h), prior) - base for m in (-2, -1, 1, 2)}
        values[name] = (delta[-2] - 8 * delta[-1] + 8 * delta[1] - delta[2]) / (12 * h)
    total = sum(values.values())
    return {
        "case_id": ard["case_id"], "universe": ard["universe"], "branch_scale": "source_branch_lengths",
        "coexist_mediated_mean": values["coexist_mediated"], "direct_endpoint_mean": values["direct_endpoint"],
        "coexist_mediated_fraction": values["coexist_mediated"] / total,
        "direct_endpoint_fraction": values["direct_endpoint"] / total,
        "expected_total": total, "log_likelihood_recomputed": base,
        "fit_log_likelihood": ard["log_likelihood"],
        "log_likelihood_abs_difference": abs(base - float(ard["log_likelihood"])),
    }


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=16)
    parser.add_argument("--root", type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve() if args.root else Path(__file__).resolve().parents[1]
    cases = sorted(path for path in (root / "cases").iterdir() if path.is_dir())
    if not cases:
        raise ValueError("No prepared cases found")
    tasks = [(str(case), *model) for case in cases for model in MODELS]
    with mp.Pool(processes=min(args.processes, len(tasks))) as pool:
        fits = pool.map(fit_task, tasks)
    for case_id in {row["case_id"] for row in fits}:
        block = [row for row in fits if row["case_id"] == case_id]
        minimum = min(row["AIC"] for row in block)
        for row in block:
            row["delta_AIC"] = row["AIC"] - minimum
    fits.sort(key=lambda row: (row["case_id"], row["model"]))
    out = root / "results" / "source_branch_length_validation"
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "five_model_fit_GTDB_source_branch_lengths.tsv", fits)
    ard_fits = {row["case_id"]: row for row in fits if row["model"] == "ARD_unrestricted"}
    with mp.Pool(processes=min(args.processes, len(cases))) as pool:
        ard_rows = pool.map(ard_task, [(str(case), ard_fits[case.name]) for case in cases])
    ard_rows.sort(key=lambda row: row["case_id"])
    if max(row["log_likelihood_abs_difference"] for row in ard_rows) > 1e-4:
        raise ValueError("ARD likelihood reproduction failed")
    write_tsv(out / "ARD_contribution_GTDB_source_branch_lengths.tsv", ard_rows)

    ard_map = {row["case_id"]: row for row in ard_rows}
    classified = []
    for case_id in sorted(ard_map):
        block = [row for row in fits if row["case_id"] == case_id]
        by_model = {row["model"]: row for row in block}
        seq = by_model["adjacency_coexist_middle"]
        other = max((row for row in block if row["model"] not in {"ARD_unrestricted", "adjacency_coexist_middle"}),
                    key=lambda row: row["log_likelihood"])
        margin = seq["log_likelihood"] - other["log_likelihood"]
        fraction = ard_map[case_id]["coexist_mediated_fraction"]
        compatible = margin >= -1e-6 and fraction >= 0.75
        strong85 = margin >= -1e-6 and fraction >= 0.85
        strong90 = margin >= -1e-6 and fraction >= 0.90
        classified.append({
            "case_id": case_id, "universe": seq["universe"], "sequential_log_likelihood": seq["log_likelihood"],
            "best_other_restricted_model": other["model"], "best_other_restricted_log_likelihood": other["log_likelihood"],
            "sequential_likelihood_margin": margin, "ARD_coexist_mediated_fraction": fraction,
            "retention_compatible_r75": str(compatible).upper(), "retention_strong_r85": str(strong85).upper(),
            "retention_strong_r90": str(strong90).upper(),
        })
    write_tsv(out / "likelihood_ARD_reclassification_source_branch_lengths.tsv", classified)
    summary = {
        "n_cases": len(cases),
        "compatible_or_strong_r75": sum(row["retention_compatible_r75"] == "TRUE" for row in classified),
        "strong_r85": sum(row["retention_strong_r85"] == "TRUE" for row in classified),
        "strong_r90": sum(row["retention_strong_r90"] == "TRUE" for row in classified),
        "all_fits_success": all(row["success"] == "TRUE" for row in fits),
        "max_ARD_reproduction_abs_difference": max(row["log_likelihood_abs_difference"] for row in ard_rows),
    }
    (out / "classification_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (out / "COMPLETE.marker").write_text(f"status=complete\ncases={len(cases)}\nbranch_scale=source_branch_lengths\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
