#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
for case_dir in "$here"/cases/archaea1338__OG*; do
  python3 "$here/code/fit_one_case.py" "$case_dir"
  python3 "$here/code/compute_ard_contribution_one_case.py" "$case_dir"
  python3 "$here/code/fit_source_one_case.py" --case "$case_dir" --code-dir "$here/code"
done

