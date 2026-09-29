param([switch]$OpenReport)

$ErrorActionPreference = "Stop"
$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Candidates = @()

$SystemPython = Get-Command python -ErrorAction SilentlyContinue
if ($SystemPython) { $Candidates += $SystemPython.Source }

$PythonLauncher = Get-Command py -ErrorAction SilentlyContinue
if ($PythonLauncher) { $Candidates += $PythonLauncher.Source }

if ($Candidates.Count -eq 0) {
    throw "Python was not found. Install Python 3 or update the Python path in run_error_demo.ps1."
}

$Python = $Candidates[0]
Write-Host "Using Python: $Python" -ForegroundColor Cyan
Write-Host "Creating a disposable copy with four controlled errors..." -ForegroundColor Yellow

& $Python (Join-Path $ProjectDir "build_error_demo_data.py")
if ($LASTEXITCODE -ne 0) {
    throw "The controlled error dataset could not be created."
}

Write-Host "`nRunning the same protein-sequence QC workflow..." -ForegroundColor Cyan
& $Python (Join-Path $ProjectDir "sequence_qc.py") `
    --config (Join-Path $ProjectDir "error_demo_config.json") `
    --output "error_demo_outputs"

if ($LASTEXITCODE -ne 0) {
    throw "The error-detection demo stopped with exit code $LASTEXITCODE."
}

$Report = Join-Path $ProjectDir "error_demo_outputs\protein_sequence_qc_report.html"
Write-Host "`nFinished. The red CHECK cards are expected." -ForegroundColor Green
Write-Host "Report: $Report" -ForegroundColor Green

if ($OpenReport) {
    Start-Process -FilePath $Report
}
