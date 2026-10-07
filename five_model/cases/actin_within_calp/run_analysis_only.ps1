param([string]$PythonExe = 'python')
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$case = Join-Path $here 'cases/asgard_cALP_copy_associated_AB_gappyout584_GTDB_r226_full737'
& $PythonExe (Join-Path $here 'code/fit_one_case.py') $case
& $PythonExe (Join-Path $here 'code/compute_ard_contribution_one_case.py') $case
& $PythonExe (Join-Path $here 'code/run_source_branch_length_validation.py') --root $here --processes 1

