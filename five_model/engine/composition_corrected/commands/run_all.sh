#!/usr/bin/env bash
set -euo pipefail

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_EXE="${PYTHON_EXE:-python3}"
RSCRIPT_EXE="${RSCRIPT_EXE:-Rscript}"
PROCESSES="${PROCESSES:-24}"
mkdir -p "${PACKAGE_ROOT}/logs"

"${PYTHON_EXE}" "${PACKAGE_ROOT}/code/build_six_tree_tip_states.py" \
  2>&1 | tee "${PACKAGE_ROOT}/logs/00_rebuild_input.log"
"${PYTHON_EXE}" "${PACKAGE_ROOT}/code/fit_composition_corrected_five_models.py" \
  --processes "${PROCESSES}" 2>&1 | tee "${PACKAGE_ROOT}/logs/01_fit.log"
"${RSCRIPT_EXE}" "${PACKAGE_ROOT}/code/paired_tests_vs_sequential_retention.R" \
  2>&1 | tee "${PACKAGE_ROOT}/logs/02_paired_tests.log"
"${RSCRIPT_EXE}" "${PACKAGE_ROOT}/code/plot_final_five_model_figures.R" \
  2>&1 | tee "${PACKAGE_ROOT}/logs/03_plot.log"
"${PYTHON_EXE}" "${PACKAGE_ROOT}/code/audit_reproducibility.py" \
  2>&1 | tee "${PACKAGE_ROOT}/logs/04_audit.log"
"${PYTHON_EXE}" "${PACKAGE_ROOT}/code/build_manifest.py" \
  2>&1 | tee "${PACKAGE_ROOT}/logs/05_manifest.log"

printf 'status=complete\nscope=composition_corrected_five_models_only\nmodels=5\ntrees=6\nbranch_scales=2\n' \
  > "${PACKAGE_ROOT}/PACKAGE_COMPLETE.marker"

printf 'status=complete package=%s\n' "${PACKAGE_ROOT}"
