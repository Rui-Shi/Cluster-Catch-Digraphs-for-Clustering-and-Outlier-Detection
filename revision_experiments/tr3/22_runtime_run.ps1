# revision_experiments/tr3/22_runtime_run.ps1
#
# Runs every cell of tab:runtime_comparison with 22_runtime.R, one timing
# process at a time, -Reps passes over all cells (cheapest cell first in each
# pass, a fresh R process per run). Finished (method, n, rep) rows are skipped, so the
# script can be stopped and restarted. It stops at the first failure,
# including "other R processes are running" from the load guard.
#
# From the clone root, in PowerShell:
#   .\revision_experiments\tr3\22_runtime_run.ps1                  # all 12 cells, 3 reps
#   .\revision_experiments\tr3\22_runtime_run.ps1 -Methods UN-CCDs,RK-CCDs
param(
  [string[]] $Methods = @('UN-CCDs', 'RK-CCDs', 'KS-CCDs'),
  [int[]]    $Ns      = @(50, 100, 200, 500),
  [int]      $Reps    = 3,   # passes; -Reps 1 first to project the total time
  [string]   $IterN   = '500',
  [string[]] $Once    = @('KS-CCDs:500'),   # cells run in pass 1 only (author's decision: timed once)
  [string]   $Cores   = ''          # empty: all cores (detectCores(), as the drivers)
)
$Rscript = 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
$env:CCD_QUANTILE_GENERATE = 'FALSE'
if ($Cores -ne '') { $env:TR3_CORES = $Cores } else { Remove-Item Env:TR3_CORES -ErrorAction SilentlyContinue }
# cheapest first (weights from the published minutes: KS about 10x RK about 10x UN; time grows with n^2); KS-CCDs at n = 500 last
$wt = @{ 'UN-CCDs' = 1; 'RK-CCDs' = 10; 'KS-CCDs' = 100 }
$cells = foreach ($m in $Methods) { foreach ($n in $Ns) { [pscustomobject]@{ m = $m; n = $n; w = $wt[$m] * $n * $n } } }
$cells = $cells | Sort-Object w
# Pass r covers every cell once before pass r + 1 starts, so that a stop leaves complete passes.
for ($r = 1; $r -le $Reps; $r++) {
  foreach ($c in $cells) {
    if ($r -gt 1 -and ("$($c.m):$($c.n)" -in $Once)) { continue }
    Write-Host ("{0} {1} n={2} rep={3}" -f (Get-Date -Format 'HH:mm:ss'), $c.m, $c.n, $r)
    # A foreign R process at the start (the load guard) is waited out up to 10 times, 60 s each;
    # any other failure stops the loop.
    for ($try = 1; $try -le 10; $try++) {
      $out = & $Rscript revision_experiments/tr3/22_runtime.R $c.m $c.n $r $IterN 2>&1 | ForEach-Object { "$_" }
      $out | Select-Object -Last 3
      if ($LASTEXITCODE -eq 0) { break }
      if (($out -join "`n") -match 'not timing under load') { Write-Host "load guard, waiting 60 s (try $try)"; Start-Sleep -Seconds 60; continue }
      Write-Host "stopped: exit code $LASTEXITCODE"; exit $LASTEXITCODE
    }
    if ($LASTEXITCODE -ne 0) { Write-Host "stopped: load guard did not clear"; exit 1 }
  }
}
