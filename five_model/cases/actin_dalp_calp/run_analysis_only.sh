#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
python3 "$here/code/build_six_tree_tip_states.py"
python3 "$here/code/fit_composition_corrected_five_models.py" --processes "${PROCESSES:-24}"
Rscript "$here/code/paired_tests_vs_sequential_retention.R"
python3 "$here/code/audit_reproducibility.py"

