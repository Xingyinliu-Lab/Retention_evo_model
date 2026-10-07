#!/usr/bin/env python3
"""Fixed-Q ARD transition contributions for one GTDB case and branch scale."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

LEVELS = ("A-only", "A+B coexist", "B-only")
INDEX = {state: index for index, state in enumerate(LEVELS)}
DIRECTIONAL = {
    "A-only->A+B coexist": ((0, 1),),
    "A+B coexist->A-only": ((1, 0),),
    "A+B coexist->B-only": ((1, 2),),
    "B-only->A+B coexist": ((2, 1),),
    "A-only->B-only": ((0, 2),),
    "B-only->A-only": ((2, 0),),
}
GROUPS = {
    "coexist_mediated": ((0, 1), (1, 0), (1, 2), (2, 1)),
    "direct_endpoint": ((0, 2), (2, 0)),
}


class Node:
    def __init__(self):
        self.name, self.children, self.length = "", [], 0.0
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
        if ":" in token:
            node.name, length = token.rsplit(":", 1)
            node.length = float(length or 0.0)
        else:
            node.name = token
        return node
    return parse()


def stationary(q):
    values, vectors = np.linalg.eig(q.T)
    vector = np.abs(np.real(vectors[:, np.argmin(np.abs(values))]))
    return vector / vector.sum()


def parse_q(text):
    q = np.zeros((3, 3))
    for token in text.split(";"):
        edge, value = token.rsplit(":", 1)
        source, target = edge.split("->", 1)
        q[INDEX[source], INDEX[target]] = float(value)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


class Engine:
    def __init__(self, root, states, prior, branch_scale):
        self.nodes = []
        def walk(node):
            for child in node.children:
                walk(child)
            self.nodes.append(node)
        walk(root)
        index = {id(node): i for i, node in enumerate(self.nodes)}
        self.children = [
            [
                (
                    index[id(child)],
                    1.0 if branch_scale == "unit_topology" else max(float(child.length), 1e-8),
                )
                for child in node.children
            ]
            for node in self.nodes
        ]
        self.names = [node.name if node.terminal() else None for node in self.nodes]
        self.states = states
        self.prior = prior
        self.branch_count = sum(map(len, self.children))
        tips = {name for name in self.names if name is not None}
        if tips != set(states):
            raise ValueError("species-tree tips != host-state keys")

    def log_likelihood(self, generator):
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
                    cache[key] = expm(generator * length)
                vector *= cache[key] @ likelihoods[child_i]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                raise FloatingPointError(total)
            likelihoods[i] = vector / total
            log_scale += math.log(total)
        root_value = float(self.prior @ likelihoods[-1])
        return math.log(root_value) + log_scale


def tilted(q, edges, theta):
    matrix = q.copy()
    factor = math.exp(theta)
    for source, target in edges:
        matrix[source, target] = q[source, target] * factor
    np.fill_diagonal(matrix, np.diag(q))
    return matrix


def moments(engine, q, edges, h=5e-4):
    base = engine.log_likelihood(q)
    delta = {}
    for multiplier in (-2, -1, 1, 2):
        delta[multiplier] = engine.log_likelihood(tilted(q, edges, multiplier * h)) - base
    mean = (delta[-2] - 8 * delta[-1] + 8 * delta[1] - delta[2]) / (12 * h)
    variance = max(0.0, (-delta[2] + 16 * delta[1] + 16 * delta[-1] - delta[-2]) / (12 * h * h))
    return mean, variance, base


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument(
        "--branch-scale",
        choices=("unit_topology", "source_branch_lengths"),
        default="unit_topology",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    case = args.case.resolve()
    branch_scale = args.branch_scale
    if branch_scale == "unit_topology":
        fit_path = case / "results/five_model_fit_GTDB.tsv"
        default_results = case / "results"
        suffix = "GTDB"
    else:
        default_results = case / "results/source_branch_length_validation"
        fit_path = default_results / "five_model_fit_GTDB_source_branch_lengths.tsv"
        suffix = "GTDB_source_branch_lengths"
    results = args.output_dir.resolve() if args.output_dir else default_results
    with fit_path.open(encoding="utf-8", newline="") as handle:
        fits = list(csv.DictReader(handle, delimiter="\t"))
    ard = next(row for row in fits if row["model"] == "ARD_unrestricted")
    states = {}
    with (case / "input/host_states.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            states[row["host_id"]] = INDEX[row["state"]]
    root = parse_newick((case / "input/species_tree.nwk").read_text(encoding="utf-8"))
    q = parse_q(ard["rates"])
    engine = Engine(root, states, stationary(q), branch_scale)
    directional_rows, group_rows = [], []
    for name, edges in DIRECTIONAL.items():
        mean, variance, logl = moments(engine, q, edges)
        directional_rows.append({"case_id": ard["case_id"], "branch_scale": branch_scale,
                                 "transition": name, "posterior_mean": mean,
                                 "posterior_variance": variance, "posterior_sd": math.sqrt(variance)})
    values = {}
    for name, edges in GROUPS.items():
        mean, variance, logl = moments(engine, q, edges)
        values[name] = mean
        group_rows.append({"case_id": ard["case_id"], "branch_scale": branch_scale,
                           "transition_group": name, "posterior_mean": mean,
                           "posterior_variance": variance, "posterior_sd": math.sqrt(variance),
                           "log_likelihood_recomputed": logl, "fit_log_likelihood": ard["log_likelihood"],
                           "log_likelihood_abs_difference": abs(logl - float(ard["log_likelihood"]))})
    total = values["coexist_mediated"] + values["direct_endpoint"]
    summary = [{
        "case_id": ard["case_id"], "universe": ard["universe"], "tree": "GTDB_r226",
        "branch_scale": branch_scale,
        "n_hosts": ard["n_hosts"], "n_A_only": ard["n_A_only"], "n_coexist": ard["n_coexist"],
        "n_B_only": ard["n_B_only"], "coexist_mediated_mean": values["coexist_mediated"],
        "direct_endpoint_mean": values["direct_endpoint"],
        "coexist_mediated_fraction": values["coexist_mediated"] / total,
        "direct_endpoint_fraction": values["direct_endpoint"] / total,
        "expected_total": total, "expected_total_per_branch": total / engine.branch_count,
    }]
    if max(float(row["log_likelihood_abs_difference"]) for row in group_rows) > 1e-4:
        raise RuntimeError("ARD likelihood reproduction failed")
    results.mkdir(parents=True, exist_ok=True)
    write_tsv(results / f"ARD_directional_transition_posterior_moments_{suffix}.tsv", directional_rows)
    write_tsv(results / f"ARD_transition_group_posterior_moments_{suffix}.tsv", group_rows)
    contribution_name = (
        "ARD_transition_contribution_GTDB.tsv"
        if branch_scale == "unit_topology"
        else "ARD_contribution_GTDB_source_branch_lengths.tsv"
    )
    write_tsv(results / contribution_name, summary)
    marker = results / f"ARD_CONTRIBUTION_{branch_scale.upper()}_COMPLETE.marker"
    marker.write_text(
        f"status=complete\ntree=GTDB_r226\nbranch_scale={branch_scale}\n"
        "method=fixed_Q_stochastic_history_posterior_moments\n", encoding="utf-8"
    )
    print(json.dumps({"case_id": ard["case_id"], "status": "complete",
                      "branch_scale": branch_scale,
                      "coexist_mediated_fraction": summary[0]["coexist_mediated_fraction"]}))


if __name__ == "__main__":
    main()
