#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
case_dir="$here/cases/asgard_cALP_copy_associated_AB_gappyout584_GTDB_r226_full737"
python3 "$here/code/fit_one_case.py" "$case_dir"
python3 "$here/code/compute_ard_contribution_one_case.py" "$case_dir"
python3 "$here/code/run_source_branch_length_validation.py" --root "$here" --processes 1

