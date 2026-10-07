# Composition-corrected five-model validation

This directory is the frozen, supplement-ready reproduction package for the final five-model comparison of MAG-level cALP/dALP states. It contains only the composition-corrected analysis retained for the manuscript.

## Scope

- States: `dALP-only`, `cALP+dALP coexist`, and `cALP-only`.
- Six GTDB-projected species-tree reconstructions.
- Two prespecified branch treatments: unit topology and source branch lengths.
- Five unique models listed in `input/MODEL_DEFINITIONS.tsv`.
- AIC is calculated within each identical tree, tip-state data set, and branch treatment.
- Sequential retention is the prespecified paired-comparison reference.

This package does not infer which profile is ancestral, date transitions, or identify duplication, transfer, or loss events. cALP/dALP are frozen profile labels, not perfectly monophyletic clade labels.

## Reproduction

Windows PowerShell:

```powershell
.\commands\run_all.ps1 -PythonExe python -RscriptExe "Rscript" -Processes 24
```

Linux/macOS:

```bash
PYTHON_EXE=python3 RSCRIPT_EXE=Rscript PROCESSES=24 bash commands/run_all.sh
```

The pipeline fits all 60 combinations, performs exact paired signed-rank tests, recreates the two final PDFs, and audits the model set, input closure, model ranking, best model and paired-test conclusions against the archived formal run. Last-decimal likelihood differences caused by platform-specific floating-point, matrix-exponential and optimizer termination behavior are recorded descriptively rather than treated as a biological discrepancy. `REPRODUCIBILITY_PASS.marker` is written only after all critical checks pass.

## Directory contract

- `input/`: frozen tip states, Newick trees and model definitions.
- `input/source_trees/`: the six actual GTDB-projected Newick trees used to construct the CTMC input.
- `input/source_tables/`: frozen 6769-sequence metadata and 584-placement assignments used to derive MAG states and the placement-concordant cohort.
- `code/`: the complete final fitting, statistics, plotting and audit code.
- `output/`: machine-readable model fits and paired tests.
- `figures/`: the two final manuscript-support figures.
- `source_data/`: figure source data.
- `audit/`: archived formal references, optimizer diagnostics and reproducibility checks.
- `docs/`: detailed Chinese audit and English supplementary methods.
- `environment/`: tested runtime versions and dependencies.
- `commands/`: one-command reproduction entry points.

See `docs/FINAL_SCOPE_AND_EXCLUSIONS.md` before reusing any result.

`code/build_six_tree_tip_states.py` reconstructs `six_tree_tip_states.json` from the packaged source trees and tables. The pipeline requires the rebuilt JSON to be byte-identical to the frozen input before fitting begins.
