param([string]$PythonExe = 'python', [string]$RscriptExe = 'Rscript', [int]$Processes = 24)
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
& $PythonExe (Join-Path $here 'code/build_six_tree_tip_states.py')
& $PythonExe (Join-Path $here 'code/fit_composition_corrected_five_models.py') --processes $Processes
& $RscriptExe (Join-Path $here 'code/paired_tests_vs_sequential_retention.R')
& $PythonExe (Join-Path $here 'code/audit_reproducibility.py')

