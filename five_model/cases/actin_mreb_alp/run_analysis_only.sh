#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
python3 "$here/code/fit_frozen_input_no_biopython.py"
python3 "$here/code/analyze_pvalues_and_ard_contributions.py"
python3 "$here/code/fit_source_branch_lengths.py"
python3 "$here/code/analyze_source_branch_ard_contributions.py"

