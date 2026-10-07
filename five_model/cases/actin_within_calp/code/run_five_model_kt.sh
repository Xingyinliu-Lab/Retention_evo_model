#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE="$ROOT/cases/asgard_cALP_copy_associated_AB_gappyout584_GTDB_r226_full737"
CODE="$ROOT/code"
mkdir -p "$ROOT/logs"
python3 "$CODE/fit_one_case.py" "$CASE" >"$ROOT/logs/unit_five_model.log" 2>&1
python3 "$CODE/compute_ard_contribution_one_case.py" "$CASE" >"$ROOT/logs/unit_ARD.log" 2>&1
python3 "$CODE/run_source_branch_length_validation.py" --root "$ROOT" --processes 24 >"$ROOT/logs/source_five_model_ARD.log" 2>&1
printf 'status=complete\ncase=1\nbranch_scales=2\n' >"$ROOT/COMPLETE.marker"
