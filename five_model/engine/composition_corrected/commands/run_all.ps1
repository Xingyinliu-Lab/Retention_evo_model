param(
    [string]$PythonExe = "python",
    [string]$RscriptExe = "Rscript",
    [int]$Processes = 24
)

$ErrorActionPreference = "Stop"
$PackageRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$LogDir = Join-Path $PackageRoot "logs"
New-Item -ItemType Directory -Path $LogDir -Force | Out-Null

& $PythonExe (Join-Path $PackageRoot "code\build_six_tree_tip_states.py") *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "00_rebuild_input.log")
if ($LASTEXITCODE -ne 0) { throw "Six-tree input rebuild failed" }

& $PythonExe (Join-Path $PackageRoot "code\fit_composition_corrected_five_models.py") --processes $Processes *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "01_fit.log")
if ($LASTEXITCODE -ne 0) { throw "Five-model fitting failed" }

& $RscriptExe (Join-Path $PackageRoot "code\paired_tests_vs_sequential_retention.R") *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "02_paired_tests.log")
if ($LASTEXITCODE -ne 0) { throw "Paired tests failed" }

& $RscriptExe (Join-Path $PackageRoot "code\plot_final_five_model_figures.R") *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "03_plot.log")
if ($LASTEXITCODE -ne 0) { throw "Figure generation failed" }

& $PythonExe (Join-Path $PackageRoot "code\audit_reproducibility.py") *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "04_audit.log")
if ($LASTEXITCODE -ne 0) { throw "Reproducibility audit failed" }

& $PythonExe (Join-Path $PackageRoot "code\build_manifest.py") *>&1 |
    Tee-Object -FilePath (Join-Path $LogDir "05_manifest.log")
if ($LASTEXITCODE -ne 0) { throw "Manifest generation failed" }

Set-Content -LiteralPath (Join-Path $PackageRoot "PACKAGE_COMPLETE.marker") -Encoding UTF8 -Value @(
    "status=complete",
    "scope=composition_corrected_five_models_only",
    "models=5",
    "trees=6",
    "branch_scales=2"
)

Write-Output "status=complete package=$PackageRoot"
