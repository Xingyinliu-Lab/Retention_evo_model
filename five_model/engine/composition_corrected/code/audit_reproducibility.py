#!/usr/bin/env python3
"""Audit closure, numerical reproduction, optimizer diagnostics and hashes."""

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
INPUT = PACKAGE / "input/six_tree_tip_states.json"
OBSERVED = PACKAGE / "output/five_model_fits.tsv"
OBSERVED_P = PACKAGE / "output/paired_tests_vs_sequential_retention.tsv"
REFERENCE = PACKAGE / "audit/reference_official_five_model_fits.tsv"
REFERENCE_P = PACKAGE / "audit/reference_official_paired_tests.tsv"
AUDIT_DIR = PACKAGE / "audit"


def read_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def key_fit(row):
    return row["tree"], row["branch_scale"], row["model"]


def key_p(row):
    return row["branch_scale"], row["comparison_model"]


def max_abs_difference(observed, reference, key, field):
    obs = {key(row): float(row[field]) for row in observed}
    ref = {key(row): float(row[field]) for row in reference}
    if set(obs) != set(ref):
        return math.inf
    return max(abs(obs[item] - ref[item]) for item in obs)


def model_orders(rows):
    orders = {}
    blocks = {(row["tree"], row["branch_scale"]) for row in rows}
    for block in blocks:
        selected = [row for row in rows if (row["tree"], row["branch_scale"]) == block]
        orders[block] = tuple(row["model"] for row in sorted(
            selected, key=lambda row: (float(row["AIC"]), row["model"])
        ))
    return orders


def main():
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location(
        "fit_final", PACKAGE / "code/fit_composition_corrected_five_models.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    trees = module.validate_input(payload)

    observed = read_tsv(OBSERVED)
    reference = read_tsv(REFERENCE)
    observed_p = read_tsv(OBSERVED_P)
    reference_p = read_tsv(REFERENCE_P)
    expected_models = {
        "adjacency_coexist_middle", "adjacency_cALP_middle", "adjacency_dALP_middle",
        "ancestral_coexist_divergence", "ARD_unrestricted",
    }
    expected_scales = {"unit_topology", "source_branch_lengths"}
    expected_keys = {
        (item["name"], scale, model)
        for item in trees for scale in expected_scales for model in expected_models
    }
    observed_keys = {key_fit(row) for row in observed}

    max_logl = max_abs_difference(observed, reference, key_fit, "log_likelihood")
    max_aic = max_abs_difference(observed, reference, key_fit, "AIC")
    max_delta = max_abs_difference(observed, reference, key_fit, "delta_AIC")
    max_p = max_abs_difference(observed_p, reference_p, key_p, "exact_p_two_sided")
    observed_orders = model_orders(observed)
    reference_orders = model_orders(reference)
    rank_orders_identical = observed_orders == reference_orders
    best_models_identical = all(
        observed_orders[block][0] == reference_orders[block][0]
        for block in observed_orders
    )
    group_minima_ok = all(
        abs(min(float(row["delta_AIC"]) for row in observed
                if row["tree"] == tree and row["branch_scale"] == scale)) <= 1e-10
        for tree in {row["tree"] for row in observed}
        for scale in expected_scales
    )

    checks = [
        ("input_tree_count", len(trees) == 6, len(trees), "exactly six frozen gtdb_projected trees"),
        ("fit_row_count", len(observed) == 60, len(observed), "6 trees x 2 branch scales x 5 models"),
        ("fit_key_closure", observed_keys == expected_keys, len(observed_keys), "all expected keys and no extras"),
        ("model_set", {row["model"] for row in observed} == expected_models,
         ";".join(sorted({row["model"] for row in observed})), "five final models only"),
        ("branch_scale_set", {row["branch_scale"] for row in observed} == expected_scales,
         ";".join(sorted({row["branch_scale"] for row in observed})), "two prespecified scales"),
        ("optimizer_success", all(row["success"] == "TRUE" for row in observed),
         sum(row["success"] == "TRUE" for row in observed), "all 60 optimizer success flags"),
        ("finite_likelihood_AIC", all(math.isfinite(float(row["log_likelihood"])) and
                                      math.isfinite(float(row["AIC"])) for row in observed),
         len(observed), "all values finite"),
        ("delta_AIC_group_minima", group_minima_ok, group_minima_ok,
         "Delta AIC recomputed within each tree x branch-scale block"),
        ("paired_test_rows", len(observed_p) == 8, len(observed_p), "4 comparisons x 2 scales"),
        ("five_model_rank_orders_identical", rank_orders_identical, rank_orders_identical,
         "identical complete AIC ordering in every tree x branch-scale block"),
        ("best_models_identical", best_models_identical, best_models_identical,
         "identical minimum-AIC model in every tree x branch-scale block"),
        ("max_abs_exact_P_difference", max_p <= 1e-12, f"{max_p:.12g}", "tolerance 1e-12"),
    ]

    with (AUDIT_DIR / "AUDIT_SUMMARY.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status", "observed", "criterion"])
        for check, passed, value, criterion in checks:
            writer.writerow([check, "PASS" if passed else "FAIL", value, criterion])

    lower = math.exp(-10)
    upper = math.exp(10)
    diagnostic_rows = []
    for row in observed:
        rates = [float(token.rsplit(":", 1)[1]) for token in row["rates"].split(";")]
        gradient = float(row["gradient_max"])
        at_lower = sum(value <= lower * (1 + 1e-6) for value in rates)
        at_upper = sum(value >= upper * (1 - 1e-6) for value in rates)
        if gradient > 1e-2 or at_lower or at_upper:
            diagnostic_rows.append([
                row["tree"], row["branch_scale"], row["model"], gradient,
                at_lower, at_upper, "WARN_BOUNDARY_OR_FLAT_SURFACE",
            ])
    with (AUDIT_DIR / "OPTIMIZER_DIAGNOSTICS.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tree", "branch_scale", "model", "gradient_max",
                         "n_rates_at_lower_bound", "n_rates_at_upper_bound", "status"])
        writer.writerows(diagnostic_rows)

    critical_pass = all(passed for _, passed, _, _ in checks)
    audit_json = {
        "status": "PASS" if critical_pass else "FAIL",
        "scope": "composition_corrected_five_models_only",
        "n_trees": len(trees),
        "n_fit_rows": len(observed),
        "n_optimizer_warning_rows": len(diagnostic_rows),
        "max_abs_log_likelihood_difference": max_logl,
        "max_abs_AIC_difference": max_aic,
        "max_abs_delta_AIC_difference": max_delta,
        "max_abs_exact_P_difference": max_p,
        "five_model_rank_orders_identical": rank_orders_identical,
        "best_models_identical": best_models_identical,
        "numerical_difference_policy": (
            "reported descriptively; not used as a pass threshold because deterministic "
            "optimization can terminate at slightly different points across platforms"
        ),
        "input_sha256": sha256(INPUT),
    }
    (AUDIT_DIR / "REPRODUCIBILITY_AUDIT.json").write_text(
        json.dumps(audit_json, indent=2) + "\n", encoding="utf-8"
    )
    if critical_pass:
        (PACKAGE / "REPRODUCIBILITY_PASS.marker").write_text(
            "status=PASS\nscope=composition_corrected_five_models_only\n",
            encoding="utf-8",
        )
    else:
        raise SystemExit("Reproducibility audit failed; inspect audit/AUDIT_SUMMARY.tsv")


if __name__ == "__main__":
    main()
