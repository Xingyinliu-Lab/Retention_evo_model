from __future__ import annotations

import json
from pathlib import Path

from . import CONTRACT_VERSION, __version__
from .abstract_fit import fit_abstract
from .abstract_model import AbstractTreeLikelihood
from .clustering import clustering_test
from .config import FamilyConfig, load_config
from .fit import Point, fit_model, surface_A_values
from .gene_cuts import discover_gene_splits
from .inputs import ValidatedInputs, applicability, validate_inputs
from .model import StateTreeLikelihood
from .species_cuts import discover_topology_cut
from .utils import sha256, write_json, write_tsv


def _point_payload(point: Point) -> dict[str, float]:
    return {"r": point.r, "A": point.A, "normalized_G": point.A, "normalized_L": 1.0, "log_likelihood": point.log_likelihood}


def write_input_outputs(validated: ValidatedInputs, output: Path) -> None:
    write_json(output / "input_audit.json", validated.audit)
    write_tsv(
        output / "validated_gene_to_host.tsv",
        ({"gene_id": gene, "host_id": host} for gene, host in sorted(validated.gene_to_host.items())),
        ["gene_id", "host_id"],
    )
    write_tsv(
        output / "host_copy_counts.tsv",
        (
            {
                "host_id": host,
                "n_genes": validated.copy_counts[host],
                "state": "single" if validated.states[host] == 0 else "coexist_2plus",
            }
            for host in sorted(validated.states)
        ),
        ["host_id", "n_genes", "state"],
    )
    write_tsv(
        output / "gene_tree_tip_audit.tsv",
        ({"gene_id": gene, "mapped_host": validated.gene_to_host[gene], "status": "PASS"} for gene in sorted(validated.gene_to_host)),
        ["gene_id", "mapped_host", "status"],
    )
    write_tsv(
        output / "species_tree_tip_audit.tsv",
        ({"host_id": host, "n_genes": validated.copy_counts[host], "status": "PASS"} for host in sorted(validated.states)),
        ["host_id", "n_genes", "status"],
    )


def write_gene_cut_outputs(validated: ValidatedInputs, output: Path) -> None:
    rows, members = discover_gene_splits(validated.gene_tree)
    write_tsv(
        output / "gene_cut_manifest.tsv",
        rows,
        ["split_id", "n_genes_side", "n_genes_other", "branch_support", "member_sha256", "first_gene"],
    )
    write_tsv(output / "gene_block_members.tsv", members, ["split_id", "gene_id"])
    write_json(
        output / "gene_cut_audit.json",
        {
            "status": "PASS",
            "source": "main_gene_tree_only",
            "n_candidate_splits": len(rows),
            "enters_state_likelihood": False,
        },
    )


def run_analysis(config_path: str | Path, output_dir: str | Path) -> dict[str, object]:
    config: FamilyConfig = load_config(config_path)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    validated = validate_inputs(config)
    write_input_outputs(validated, output)
    app = applicability(validated)
    write_json(output / "applicability.json", app)
    write_gene_cut_outputs(validated, output)
    write_tsv(
        output / "species_cut_manifest.tsv",
        [{"species_cut_id": "full", "method": "full", "n_blocks": len(validated.states), "enters_primary_fit": True}],
        ["species_cut_id", "method", "n_blocks", "enters_primary_fit"],
    )
    write_tsv(
        output / "host_block_members.tsv",
        ({"species_cut_id": "full", "species_block": host, "host_id": host} for host in sorted(validated.states)),
        ["species_cut_id", "species_block", "host_id"],
    )

    summary: dict[str, object] = {
        "contract_version": CONTRACT_VERSION,
        "software_version": __version__,
        "family_id": config.family_id,
        "applicability": app,
    }
    if app["status"] != "FIT_ELIGIBLE":
        summary["status"] = app["status"]
        write_json(output / "summary.json", summary)
        return summary

    cluster = clustering_test(validated.species_tree, validated.states, config.permutations, config.clustering_seed)
    write_json(output / "clustering_test.json", cluster)

    model = StateTreeLikelihood(validated.species_tree, validated.states)
    fit = fit_model(model, config.fit_starts, config.fit_seed)
    null: Point = fit["null"]
    free: Point = fit["free"]
    profile: list[Point] = fit["profile"]
    free_detail = model.evaluate(free.A, free.r, with_expectations=True)
    null_detail = model.evaluate(null.A, 1.0, with_expectations=True)

    fit_payload = {
        "model": "coexistence_aware_collapsed_dt_gain_loss",
        "branch_length_policy": "unit",
        "null": _point_payload(null),
        "free": _point_payload(free),
        "delta_log_likelihood": fit["delta_log_likelihood"],
        "lrt_statistic": fit["lrt_statistic"],
        "p_persistence_one_sided": fit["p_persistence_one_sided"],
        "r_profile_low": fit["r_profile_low"],
        "r_profile_high": fit["r_profile_high"],
        "status": fit["status"],
        "expected_macro_edge_transitions": free_detail.expected,
    }
    write_json(output / "fit.json", fit_payload)
    write_tsv(
        output / "r_profile.tsv",
        (
            {
                "family_id": config.family_id,
                "r": point.r,
                "A_hat": point.A,
                "logL": point.log_likelihood,
                "delta_logL_from_free": point.log_likelihood - free.log_likelihood,
            }
            for point in profile
        ),
        ["family_id", "r", "A_hat", "logL", "delta_logL_from_free"],
    )

    test_row = {
        "family_id": config.family_id,
        "species_cut_id": "full",
        **{k: app[k] for k in ("n_hosts", "n_single", "n_coexist", "coexist_fraction", "coexist_observed")},
        **{k: cluster[k] for k in ("fitch_steps_observed", "permutation_reps", "p_cluster")},
        "r_null": 1.0,
        "r_hat": free.r,
        "A_null": null.A,
        "A_hat": free.A,
        "logL_null": null.log_likelihood,
        "logL_free": free.log_likelihood,
        "delta_logL": fit["delta_log_likelihood"],
        "LRT_stat": fit["lrt_statistic"],
        "p_persistence_one_sided": fit["p_persistence_one_sided"],
        "r_profile_low": fit["r_profile_low"],
        "r_profile_high": fit["r_profile_high"],
        "clustering_status": cluster["status"],
        "persistence_status": fit["status"],
    }
    test_fields = list(test_row)
    write_tsv(output / "coexistence_test.tsv", [test_row], test_fields)

    r_values = sorted(set(config.r_values) | {1.0, free.r})
    A_values = sorted(set(surface_A_values(free.A, config.A_values)) | {null.A, free.A})
    surface_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    for r in r_values:
        for A in A_values:
            point = model.evaluate(A, r, with_expectations=True)
            is_null = abs(r - 1.0) < 1e-12 and abs(A - null.A) / null.A < 1e-8
            is_free = abs(r - free.r) / free.r < 1e-8 and abs(A - free.A) / free.A < 1e-8
            surface_rows.append(
                {
                    "family_id": config.family_id,
                    "species_cut_id": "full",
                    "r": r,
                    "A": A,
                    "normalized_G": A,
                    "normalized_L": 1.0,
                    "kappa_hat": "NA",
                    "logL": point.log_likelihood,
                    "delta_logL_vs_null_mle": point.log_likelihood - null.log_likelihood,
                    "is_null": is_null,
                    "is_free_mle": is_free,
                    "status": "PASS",
                }
            )
            event_rows.append(
                {
                    "family_id": config.family_id,
                    "species_cut_id": "full",
                    "r": r,
                    "A": A,
                    "kappa_hat": "NA",
                    "logL": point.log_likelihood,
                    **point.expected,
                    "status": "PASS",
                }
            )
    write_tsv(output / "likelihood_surface.tsv", surface_rows, list(surface_rows[0]))
    write_tsv(output / "event_expectations_surface.tsv", event_rows, list(event_rows[0]))
    write_tsv(output / "branch_state_posteriors.tsv", free_detail.branch_rows or [], list((free_detail.branch_rows or [{}])[0]))
    write_tsv(output / "branch_state_posteriors_r1.tsv", null_detail.branch_rows or [], list((null_detail.branch_rows or [{}])[0]))
    write_tsv(
        output / "expected_macro_edge_transitions.tsv",
        [{"family_id": config.family_id, "species_cut_id": "full", "r": free.r, "A": free.A, **(free_detail.expected or {})}],
        ["family_id", "species_cut_id", "r", "A", *(free_detail.expected or {}).keys()],
    )

    cut_manifest = [{"species_cut_id": "full", "method": "full", "n_blocks": len(validated.states), "enters_primary_fit": True}]
    cut_members = [
        {"species_cut_id": "full", "species_block": host, "host_id": host}
        for host in sorted(validated.states)
    ]
    cut_rows = [
        {
            "family_id": config.family_id,
            "species_cut_id": "full",
            "n_blocks": len(validated.states),
            "null_r": 1.0,
            "null_A": null.A,
            "null_kappa": "NA",
            "free_A": free.A,
            "free_r": free.r,
            "free_kappa": "NA",
            "null_logL": null.log_likelihood,
            "free_logL": free.log_likelihood,
            "delta_logL": fit["delta_log_likelihood"],
            "LRT_stat": fit["lrt_statistic"],
            "p_persistence_one_sided": fit["p_persistence_one_sided"],
            "interpretation_scope": "primary_full_species_tree",
        }
    ]
    for k in config.species_cut_k:
        cut = discover_topology_cut(validated.species_tree, validated.states, k)
        abstract_fit = fit_abstract(AbstractTreeLikelihood(cut.root), config.fit_starts, config.fit_seed + k)
        cut_manifest.append({"species_cut_id": cut.cut_id, "method": "topology_k", "n_blocks": cut.k, "enters_primary_fit": False})
        cut_members.extend(cut.members)
        cut_rows.append({"family_id": config.family_id, "species_cut_id": cut.cut_id, "n_blocks": cut.k, **abstract_fit})
    write_tsv(
        output / "species_cut_manifest.tsv",
        cut_manifest,
        ["species_cut_id", "method", "n_blocks", "enters_primary_fit"],
    )
    write_tsv(
        output / "host_block_members.tsv",
        cut_members,
        ["species_cut_id", "species_block", "host_id"],
    )
    write_tsv(output / "cut_sensitivity.tsv", cut_rows, list(cut_rows[0]))

    summary.update(
        {
            "status": "COMPLETE",
            "clustering": cluster,
            "fit": fit_payload,
            "outputs": {
                name: {"bytes": (output / name).stat().st_size, "sha256": sha256(output / name)}
                for name in (
                    "coexistence_test.tsv",
                    "likelihood_surface.tsv",
                    "event_expectations_surface.tsv",
                    "expected_macro_edge_transitions.tsv",
                )
            },
        }
    )
    write_json(output / "summary.json", summary)
    write_json(
        output / "SUMMARY_COMPLETE.json",
        {
            "status": "COMPLETE",
            "contract_version": CONTRACT_VERSION,
            "config_sha256": sha256(config.path),
            "species_tree_sha256": validated.audit["species_tree_sha256"],
            "gene_tree_sha256": validated.audit["gene_tree_sha256"],
            "genes_sha256": validated.audit["genes_sha256"],
            "summary_sha256": sha256(output / "summary.json"),
        },
    )
    return summary


def validate_only(config_path: str | Path) -> dict[str, object]:
    return validate_inputs(load_config(config_path)).audit
