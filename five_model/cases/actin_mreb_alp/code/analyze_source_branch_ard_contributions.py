#!/usr/bin/env python3
"""Paired AIC summaries and ARD route decomposition on source branch lengths."""

import csv, io, itertools, json, math
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/source_branch_lengths_20260820"
FIT = RESULT / "MreB_ALP_five_models_source_branch_lengths.tsv"
INPUT = ROOT / "input/six_tree_states.json"
LEVELS = ("M-only", "M+ALP coexist", "ALP-only")
INDEX = {x: i for i, x in enumerate(LEVELS)}
GROUPS = {
    "coexist_mediated": ((0, 1), (1, 0), (1, 2), (2, 1)),
    "direct_endpoint": ((0, 2), (2, 0)),
}


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


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
    def __init__(self, tree, states, prior):
        self.nodes = list(tree.find_clades(order="postorder"))
        index = {id(node): i for i, node in enumerate(self.nodes)}
        self.children = [[(index[id(child)], max(float(child.branch_length or 0.0), 1e-8))
                          for child in node.clades] for node in self.nodes]
        self.names = [node.name if node.is_terminal() else None for node in self.nodes]
        self.states, self.prior = states, prior
        self.branch_count = sum(map(len, self.children))

    def logl(self, matrix):
        likes, cache, scale = [np.ones(3) for _ in self.nodes], {}, 0.0
        for i, node in enumerate(self.nodes):
            if node.is_terminal():
                vector = np.zeros(3); vector[self.states[self.names[i]]] = 1; likes[i] = vector
                continue
            vector = np.ones(3)
            for child, length in self.children[i]:
                key = round(length, 12)
                if key not in cache:
                    cache[key] = expm(matrix * length)
                vector *= cache[key] @ likes[child]
            total = vector.sum()
            if total <= 0 or not np.isfinite(total):
                raise FloatingPointError(total)
            likes[i] = vector / total; scale += math.log(total)
        root = float(self.prior @ likes[-1])
        return math.log(root) + scale


def tilted(q, edges, theta):
    matrix, factor = q.copy(), math.exp(theta)
    for source, target in edges:
        matrix[source, target] = q[source, target] * factor
    np.fill_diagonal(matrix, np.diag(q))
    return matrix


def moments(engine, q, edges, h=5e-4):
    base = engine.logl(q)
    delta = {k: engine.logl(tilted(q, edges, k * h)) - base for k in (-2, -1, 1, 2)}
    mean = (delta[-2] - 8 * delta[-1] + 8 * delta[1] - delta[2]) / (12 * h)
    variance = max(0.0, (-delta[2] + 16 * delta[1] + 16 * delta[-1] - delta[-2]) / (12 * h * h))
    return mean, variance, base


def exact_signed_rank(values):
    values = [x for x in values if abs(x) > 1e-12]
    n, order = len(values), sorted(range(len(values)), key=lambda i: abs(values[i]))
    ranks, j = [0.0] * n, 0
    while j < n:
        k = j + 1
        while k < n and abs(abs(values[order[k]]) - abs(values[order[j]])) < 1e-12:
            k += 1
        rank = (j + 1 + k) / 2
        for z in order[j:k]: ranks[z] = rank
        j = k
    observed, total = sum(r for r, x in zip(ranks, values) if x > 0), sum(ranks)
    stats = [sum(r for r, sign in zip(ranks, signs) if sign) for signs in itertools.product((0, 1), repeat=n)]
    return sum(abs(x - total / 2) >= abs(observed - total / 2) - 1e-12 for x in stats) / len(stats)


def main():
    rows = read_tsv(FIT)
    trees, benchmark = sorted({row["tree"] for row in rows}), "adjacency_coexist_middle"
    lookup = {(row["tree"], row["model"]): float(row["AIC"]) for row in rows}
    paired = []
    for model in sorted({row["model"] for row in rows} - {benchmark}):
        differences = [lookup[(tree, model)] - lookup[(tree, benchmark)] for tree in trees]
        paired.append({
            "benchmark": benchmark, "comparison_model": model, "n_trees": len(trees),
            "retention_better_trees": sum(x > 0 for x in differences),
            "comparison_better_trees": sum(x < 0 for x in differences),
            "median_AIC_other_minus_retention": float(np.median(differences)),
            "exact_p_two_sided": exact_signed_rank(differences),
        })
    write_tsv(RESULT / "paired_tests_vs_sequential_retention_source_branch_lengths.tsv", paired)

    payload = json.loads(INPUT.read_text(encoding="utf-8")); items = {x["name"]: x for x in payload["trees"]}
    summaries, groups = [], []
    for fit in [row for row in rows if row["model"] == "ARD_unrestricted"]:
        item = items[fit["tree"]]
        tree = Phylo.read(io.StringIO(item["newick"]), "newick")
        states = {key: INDEX[value] for key, value in item["states"].items()}
        q = parse_q(fit["rates"]); engine = Engine(tree, states, stationary(q)); values = {}
        for group, edges in GROUPS.items():
            mean, variance, logl = moments(engine, q, edges); values[group] = (mean, variance)
            groups.append({
                "tree": fit["tree"], "branch_scale": "source_branch_lengths", "transition_group": group,
                "posterior_mean": mean, "posterior_variance": variance, "posterior_sd": math.sqrt(variance),
                "log_likelihood_recomputed": logl, "fit_log_likelihood": fit["log_likelihood"],
                "log_likelihood_abs_difference": abs(logl - float(fit["log_likelihood"])),
            })
        adjacent, direct = values["coexist_mediated"][0], values["direct_endpoint"][0]
        total = adjacent + direct
        summaries.append({
            "tree": fit["tree"], "branch_scale": "source_branch_lengths",
            "coexist_mediated_mean": adjacent, "direct_endpoint_mean": direct,
            "coexist_mediated_fraction": adjacent / total, "direct_endpoint_fraction": direct / total,
            "expected_total": total, "expected_total_per_branch": total / engine.branch_count,
        })
    if max(float(row["log_likelihood_abs_difference"]) for row in groups) > 1e-4:
        raise SystemExit("ARD likelihood reproduction failed")
    write_tsv(RESULT / "MreB_ALP_ARD_transition_group_posterior_moments_source_branch_lengths.tsv", groups)
    write_tsv(RESULT / "MreB_ALP_ARD_transition_contribution_source_branch_lengths.tsv", summaries)
    (RESULT / "ARD_SOURCE_BRANCH_COMPLETE.marker").write_text(
        "status=complete\ntrees=6\nbranch_scale=source_branch_lengths\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "ARD_trees": len(summaries)}))


if __name__ == "__main__":
    main()
