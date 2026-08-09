# regen_wilt_rk_launcher.ps1 -- wilt x {U-MCCD, SU-MCCD} (RK-based half of the
# 2x2 design family). Split from the NND half so the 711 MB RK d=5 quantile
# table is loaded once, in one process, instead of twice.
#
# Owned by the main session, not by a subagent: these cells were measured at
# ~1130 s (U-MCCD) and are the only ones that exceed the 600 s subagent tool cap.

$ErrorActionPreference = "Continue"
$Root = "G:\Submissions\TR1\TR1_Neurocomputing_resubmit\Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection"
Set-Location $Root

$env:WP0_GATE_OUT_CSV     = "$Root\revision_experiments\results\regen_proposed_wilt_rk.csv"
$env:WP0_GATE_TIMEOUT_SEC = "5400"

$R   = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$Log = "$Root\revision_experiments\results\regen_wilt_rk.log"

"=== wilt RK half started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8

& $R "revision_experiments/13_wp0_gate.R" "wilt" "U-MCCD,SU-MCCD" "half_contam" *>&1 |
    Out-File -FilePath $Log -Encoding utf8 -Append

"=== full_contam reading for SU-MCCD $(Get-Date -Format 'HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8 -Append

& $R "revision_experiments/13_wp0_gate.R" "wilt" "SU-MCCD" "full_contam" *>&1 |
    Out-File -FilePath $Log -Encoding utf8 -Append

"=== wilt RK half finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" | Out-File -FilePath $Log -Encoding utf8 -Append
"DONE" | Out-File -FilePath "$Root\revision_experiments\results\regen_wilt_rk.DONE" -Encoding utf8
