param([string]$PythonExe = 'python')
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
& $PythonExe (Join-Path $here 'code/fit_frozen_input_no_biopython.py')
& $PythonExe (Join-Path $here 'code/analyze_pvalues_and_ard_contributions.py')
& $PythonExe (Join-Path $here 'code/fit_source_branch_lengths.py')
& $PythonExe (Join-Path $here 'code/analyze_source_branch_ard_contributions.py')

