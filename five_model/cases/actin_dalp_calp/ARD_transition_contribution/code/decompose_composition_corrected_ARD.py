#!/usr/bin/env python3
"""Decompose transition routes for the frozen composition-corrected ARD fits.

This script does not fit an evolutionary model. It conditions on the final
six GTDB-projected trees, their frozen tip states, and the fitted unrestricted
ARD Q matrices from the composition-corrected five-model package. Posterior
transition-count moments are obtained from derivatives of the CTMC history
probability-generating function.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
# Packaged layout: the frozen five-model input/output/code directories are the
# parent of ARD_transition_contribution.
FROZEN = ROOT.parent
INPUT_JSON = FROZEN / "input" / "six_tree_tip_states.json"
FIT_TSV = FROZEN / "output" / "five_model_fits.tsv"
FIT_CODE = FROZEN / "code" / "fit_composition_corrected_five_models.py"
OUT = ROOT / "source_data"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


FINAL = load_module(FIT_CODE, "final_five_model_for_ard_decomposition")
STATE_TO_INDEX = FINAL.STATE
LEVEL = FINAL.LEVEL

TRANSITIONS = {
    "dALP-only->coexist": (0, 1),
    "dALP-only->cALP-only": (0, 2),
    "coexist->dALP-only": (1, 0),
    "coexist->cALP-only": (1, 2),
    "cALP-only->dALP-only": (2, 0),
    "cALP-only->coexist": (2, 1),
}
GROUPS = {
    "coexist-mediated_adjacent": ((0, 1), (1, 0), (1, 2), (2, 1)),
    "direct_endpoint": ((0, 2), (2, 0)),
    "all_transitions": tuple(TRANSITIONS.values()),
}


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_q(rate_text: str):
    q = np.zeros((3, 3), dtype=float)
    for token in rate_text.split(";"):
        edge, value = token.rsplit(":", 1)
        source, target = edge.split("->", 1)
        q[STATE_TO_INDEX[source], STATE_TO_INDEX[target]] = float(value)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


class TreeLikelihood:
    def __init__(self, root, tip_states, branch_scale, root_prior):
        self.nodes = []

        def walk(node):
            for child in node.children:
                walk(child)
            self.nodes.append(node)

        walk(root)
        self.index = {id(node): i for i, node in enumerate(self.nodes)}
        self.children = [
            [
                (
                    self.index[id(child)],
                    1.0 if branch_scale == "unit_topology" else max(float(child.branch_length), 1e-8),
                )
                for child in node.children
            ]
            for node in self.nodes
        ]
        self.names = [node.name if node.terminal() else None for node in self.nodes]
        self.tip_states = tip_states
        self.root_prior = root_prior
        self.branch_count = sum(len(x) for x in self.children)

    @staticmethod
    def transition_factory(matrix):
        values, vectors = np.linalg.eig(matrix)
        inverse = np.linalg.inv(vectors)
        cache = {}

        def transition(length):
            key = round(length, 14)
            if key not in cache:
                result = (vectors * np.exp(values * length)) @ inverse
                result = np.real_if_close(result, tol=1000).real
                result[np.abs(result) < 1e-15] = 0.0
                cache[key] = result
            return cache[key]

        return transition

    def log_likelihood(self, matrix):
        transition = self.transition_factory(matrix)
        likelihoods = [np.ones(3, dtype=float) for _ in self.nodes]
        log_scale = 0.0
        for i, node in enumerate(self.nodes):
            if node.terminal():
                vector = np.zeros(3, dtype=float)
                vector[self.tip_states[self.names[i]]] = 1.0
                likelihoods[i] = vector
                continue
            vector = np.ones(3, dtype=float)
            for child_i, length in self.children[i]:
                vector *= transition(length) @ likelihoods[child_i]
            total = float(vector.sum())
            if not math.isfinite(total) or total <= 0:
                raise FloatingPointError(f"Invalid likelihood scale: {total}")
            likelihoods[i] = vector / total
            log_scale += math.log(total)
        root_like = float(np.dot(self.root_prior, likelihoods[-1]))
        if not math.isfinite(root_like) or root_like <= 0:
            raise FloatingPointError(f"Invalid root likelihood: {root_like}")
        return math.log(root_like) + log_scale


def tilted_matrix(q, counted, theta):
    matrix = q.copy()
    multiplier = math.exp(theta)
    for i, j in counted:
        matrix[i, j] = q[i, j] * multiplier
    np.fill_diagonal(matrix, np.diag(q))
    return matrix


def posterior_moments(engine, q, counted, h=5e-4):
    ll0 = engine.log_likelihood(q)
    values = {}
    for multiplier in (-2, -1, 1, 2):
        theta = multiplier * h
        values[multiplier] = engine.log_likelihood(tilted_matrix(q, counted, theta)) - ll0
    mean = (values[-2] - 8 * values[-1] + 8 * values[1] - values[2]) / (12 * h)
    variance = (-values[2] + 16 * values[1] + 16 * values[-1] - values[-2]) / (12 * h * h)
    variance = max(0.0, variance)
    sd = math.sqrt(variance)
    return {
        "posterior_mean": mean,
        "posterior_variance": variance,
        "posterior_sd": sd,
        "posterior_95_low_normal_approx": max(0.0, mean - 1.96 * sd),
        "posterior_95_high_normal_approx": mean + 1.96 * sd,
        "log_likelihood_recomputed": ll0,
        "derivative_step": h,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    payload = json.loads(INPUT_JSON.read_text(encoding="utf-8"))
    items = {x["name"]: x for x in FINAL.validate_input(payload)}
    fits = [x for x in read_tsv(FIT_TSV) if x["model"] == "ARD_unrestricted"]
    if len(fits) != 12:
        raise ValueError(f"Expected 12 final ARD fits; observed {len(fits)}")

    group_rows = []
    transition_rows = []
    summary_rows = []

    for fit in fits:
        item = items[fit["tree"]]
        tree = FINAL.parse_newick(item["newick"])
        coded = {name: STATE_TO_INDEX[state] for name, state in item["states"].items()}
        q = parse_q(fit["rates"])
        engine = TreeLikelihood(tree, coded, fit["branch_scale"], FINAL.stationary(q))

        group_moments = {}
        for group, transitions in GROUPS.items():
            moments = posterior_moments(engine, q, transitions)
            group_moments[group] = moments
            group_rows.append({
                "tree": fit["tree"],
                "branch_scale": fit["branch_scale"],
                "transition_group": group,
                "mapped_tips": len(coded),
                "fitted_log_likelihood": fit["log_likelihood"],
                "log_likelihood_abs_difference": abs(float(fit["log_likelihood"]) - moments["log_likelihood_recomputed"]),
                **moments,
            })

        for label, edge in TRANSITIONS.items():
            moments = posterior_moments(engine, q, (edge,))
            transition_rows.append({
                "tree": fit["tree"],
                "branch_scale": fit["branch_scale"],
                "transition": label,
                "transition_class": "direct_endpoint" if edge in GROUPS["direct_endpoint"] else "coexist-mediated_adjacent",
                "fitted_rate": q[edge],
                **moments,
            })

        adjacent = group_moments["coexist-mediated_adjacent"]
        direct = group_moments["direct_endpoint"]
        total = group_moments["all_transitions"]
        covariance = (total["posterior_variance"] - adjacent["posterior_variance"] - direct["posterior_variance"]) / 2
        denominator = adjacent["posterior_mean"] + direct["posterior_mean"]
        direct_fraction = direct["posterior_mean"] / denominator if denominator > 0 else float("nan")
        grad_d = adjacent["posterior_mean"] / (denominator * denominator)
        grad_a = -direct["posterior_mean"] / (denominator * denominator)
        var_fraction = (
            grad_d * grad_d * direct["posterior_variance"]
            + grad_a * grad_a * adjacent["posterior_variance"]
            + 2 * grad_d * grad_a * covariance
        )
        sd_fraction = math.sqrt(max(0.0, var_fraction))
        max_rate = float(np.max(q - np.diag(np.diag(q))))
        per_branch = total["posterior_mean"] / engine.branch_count
        identifiability = (
            "boundary_saturated"
            if max_rate >= math.exp(10) * 0.99
            else "high_turnover"
            if per_branch > 1.0
            else "stable"
        )
        summary_rows.append({
            "tree": fit["tree"],
            "branch_scale": fit["branch_scale"],
            "mapped_tips": len(coded),
            "coexist_mediated_mean": adjacent["posterior_mean"],
            "direct_endpoint_mean": direct["posterior_mean"],
            "direct_fraction": direct_fraction,
            "direct_fraction_95_low_delta": max(0.0, direct_fraction - 1.96 * sd_fraction),
            "direct_fraction_95_high_delta": min(1.0, direct_fraction + 1.96 * sd_fraction),
            "expected_total_mean": total["posterior_mean"],
            "expected_total_per_branch": per_branch,
            "max_fitted_rate": max_rate,
            "history_identifiability": identifiability,
            "fitted_log_likelihood": float(fit["log_likelihood"]),
            "recomputed_log_likelihood": total["log_likelihood_recomputed"],
            "log_likelihood_abs_difference": abs(float(fit["log_likelihood"]) - total["log_likelihood_recomputed"]),
        })

    write_tsv(OUT / "ARD_transition_group_posterior_moments_composition_corrected.tsv", group_rows)
    write_tsv(OUT / "ARD_individual_transition_posterior_moments_composition_corrected.tsv", transition_rows)
    write_tsv(OUT / "ARD_transition_contribution_summary_composition_corrected.tsv", summary_rows)

    audit = {
        "status": "complete",
        "source_model": "frozen composition-corrected ARD_unrestricted",
        "trees": len(items),
        "tree_branch_sets": len(summary_rows),
        "new_model_fitted": False,
        "root_prior": "stationary(Q)",
        "maximum_log_likelihood_reproduction_error": max(x["log_likelihood_abs_difference"] for x in summary_rows),
        "uncertainty": "fixed-Q stochastic-history posterior moments; normal/delta 95% approximation",
        "parameter_uncertainty_included": False,
        "tree_topology_uncertainty_in_interval": False,
        "identifiability_rule": "boundary_saturated if max rate >= 0.99*exp(10); high_turnover if expected transitions per branch > 1; otherwise stable",
        "frozen_input_json": str(INPUT_JSON),
        "frozen_fit_tsv": str(FIT_TSV),
    }
    (OUT / "ARD_transition_contribution_composition_corrected_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
