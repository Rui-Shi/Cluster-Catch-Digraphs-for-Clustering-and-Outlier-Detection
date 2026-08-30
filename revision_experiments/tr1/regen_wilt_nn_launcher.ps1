# regen_wilt_nn_launcher.ps1 -- wilt x {UN-MCCD, SUN-MCCD} (NND-based half of
# the 2x2 design family). Runs concurrently with the RK half; the NN d=5
# quantile table is ~70 KB, so the two processes together stay well inside RAM.
#
# Owned by the main session, not by a subagent: SUN-MCCD on wilt was measured
# at ~948 s, above the 600 s subagent tool cap.

$ErrorActionPreference = "Continue"
$Root = "G:\Submissions\TR1\TR1_Neurocomputing_resubmit\Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection"
Set-Location $Root

$env:WP0_GATE_OUT_CSV     = "$Root\revision_experiments\results\tr1\regen_proposed_wilt_nn.csv"
$env:WP0_GATE_TIMEOUT_SEC = "5400"

$R   = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$Log = "$Root\revision_experiments\results\tr1\regen_wilt_nn.log"

"=== wilt NN half started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8

& $R "revision_experiments/tr1/13_wp0_gate.R" "wilt" "UN-MCCD,SUN-MCCD" "half_contam" *>&1 |
    Out-File -FilePath $Log -Encoding utf8 -Append

"=== full_contam reading for SUN-MCCD $(Get-Date -Format 'HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8 -Append

& $R "revision_experiments/tr1/13_wp0_gate.R" "wilt" "SUN-MCCD" "full_contam" *>&1 |
    Out-File -FilePath $Log -Encoding utf8 -Append

"=== wilt NN half finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8 -Append
"DONE" | Out-File -FilePath "$Root\revision_experiments\results\tr1\regen_wilt_nn.DONE" -Encoding utf8
