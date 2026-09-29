param([switch]$OpenReport)

$ErrorActionPreference = "Stop"
$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Candidates = @()

$SystemPython = Get-Command python -ErrorAction SilentlyContinue
if ($SystemPython) { $Candidates += $SystemPython.Source }

$PythonLauncher = Get-Command py -ErrorAction SilentlyContinue
if ($PythonLauncher) { $Candidates += $PythonLauncher.Source }

if ($Candidates.Count -eq 0) {
    throw "Python was not found. Install Python 3 or update the Python path in run_demo.ps1."
}

$Python = $Candidates[0]
Write-Host "Using Python: $Python" -ForegroundColor Cyan
Write-Host "Running protein-sequence QC demo..." -ForegroundColor Cyan

& $Python (Join-Path $ProjectDir "sequence_qc.py") `
    --config (Join-Path $ProjectDir "config.json") `
    --output "outputs"

if ($LASTEXITCODE -ne 0) {
    throw "The Python demo stopped with exit code $LASTEXITCODE."
}

$Report = Join-Path $ProjectDir "outputs\protein_sequence_qc_report.html"
Write-Host "`nFinished. Report: $Report" -ForegroundColor Green

if ($OpenReport) {
    Start-Process -FilePath $Report
}
