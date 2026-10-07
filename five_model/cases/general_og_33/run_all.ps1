$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$code = Join-Path $here 'code'
$python = 'python'
Get-ChildItem -LiteralPath (Join-Path $here 'cases') -Directory | Sort-Object Name | ForEach-Object {
    & $python (Join-Path $code 'fit_one_case.py') $_.FullName
    & $python (Join-Path $code 'compute_ard_contribution_one_case.py') $_.FullName
    & $python (Join-Path $code 'fit_source_one_case.py') --case $_.FullName --code-dir $code
}

