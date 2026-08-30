# 14_patch_manuscript_tables.ps1
#
# Replaces the four MCCD rows in TR2's two real-data tables with the values
# regenerated for TR1 in this session.
#
# Authorised by the author as a ONE-TIME write outside the TR1 repository.
# CLAUDE.md otherwise forbids writes outside G:\Submissions\TR1\TR1_Neurocomputing_resubmit.
#
# Why this is the right substitution, not a cross-paper contamination:
# the audit established that TR2's four MCCD rows for hepatitis, glass,
# vertebral, stamps, vowels and wilt are verbatim copies of TR1's ORIGINAL
# published table -- every logged discrepancy matches to three decimals. They
# were never independently computed. TR2 also states no S_min and no MC-SRT
# alpha for these rows (its "threshold of 2" and "S_min = 0" apply to the four
# OS methods only, line 1496-1497), so there is no competing configuration to
# conflict with. Replacing them with TR1's corrected values is therefore a
# correction of the same numbers, not an import of foreign ones.
#
# Scope: the four MCCD rows only. OS rows, baseline rows, the description
# table, the prose and tab:highdim are left untouched for the TR2 session.
#
# Line replacement is used rather than string matching because the source rows
# contain literal tab characters. Each target line is verified by its row label
# before being overwritten; the script aborts if any check fails.

$ErrorActionPreference = "Stop"

$TR2 = "C:\Users\shiru\Box\RuiShi_OutlierDetection\Submissions\TR2 Outlyingness Score\Pattern Recognition (Elsevier) - Resubmit\CCDwScores.tex"

if (-not (Test-Path $TR2)) { throw "TR2 manuscript not found: $TR2" }

# Timestamped backup beside the original, so the TR2 session can diff or revert.
$stamp  = Get-Date -Format "yyyyMMdd_HHmmss"
$backup = "$TR2.bak_$stamp"
Copy-Item -LiteralPath $TR2 -Destination $backup
"backup written: $backup"

$lines = [System.IO.File]::ReadAllLines($TR2)
"file has {0} lines" -f $lines.Count

# 1-based line numbers of the eight rows to replace.
$tprBlock = 1505..1508     # Real_Data_Result_OS1, TPR/TNR
$baBlock  = 1527..1530     # Real_Data_Result_OS2, BA/F2
$labels   = @("U-MCCDs", "SU-MCCDs", "UN-MCCDs", "SUN-MCCDs")

# Verify before touching anything.
foreach ($blk in @($tprBlock, $baBlock)) {
    for ($i = 0; $i -lt 4; $i++) {
        $ln  = $blk[$i]
        $txt = $lines[$ln - 1]
        if ($txt.Trim() -notlike ($labels[$i] + "*")) {
            throw "line $ln does not start with $($labels[$i]): '$($txt.Trim())'"
        }
    }
}
"verified: all 8 target lines carry the expected row labels"

# Column order in both TR2 tables:
# hepatitis, lymph, glass, WBC, vertebral, stamps, WDBC, vowels, thyroid, wilt
# Values from revision_experiments/results/.../final_comparison.csv (TR1, S_min = 0.0625).

$newTpr = @(
'  U-MCCDs   & 0.286 & 0.910 & 1.000 & 0.901 & 1.000 & 0.441 & 1.000 & 0.493 & 0.433 & 0.638 & 0.065 & 0.958 & 0.500 & 0.714 & 1.000 & 0.327 & 1.000 & 0.685 & 0.763 & 0.630 \\',
'  SU-MCCDs  & 0.286 & 0.925 & 1.000 & 0.923 & 1.000 & 0.647 & 1.000 & 0.493 & 0.200 & 0.576 & 1.000 & 0.168 & 1.000 & 0.073 & 1.000 & 0.317 & 1.000 & 0.732 & 0.296 & 0.785 \\',
'  UN-MCCDs  & 0.714 & 0.657 & 1.000 & 0.859 & 0.222 & 0.765 & 1.000 & 0.498 & 0.033 & 0.914 & 0.484 & 0.812 & 1.000 & 0.182 & 1.000 & 0.541 & 1.000 & 0.630 & 0.140 & 0.898 \\',
'  SUN-MCCDs & 0.714 & 0.657 & 1.000 & 0.880 & 1.000 & 0.539 & 1.000 & 0.563 & 0.100 & 0.929 & 0.516 & 0.883 & 1.000 & 0.325 & 0.978 & 0.563 & 0.989 & 0.798 & 0.366 & 0.745 \\'
)

$newBa = @(
'  U-MCCDs   & 0.598 & 0.278 & 0.951 & 0.682 & 0.721 & 0.283 & 0.746 & 0.316 & 0.536 & 0.311 & 0.511 & 0.072 & 0.607 & 0.170 & 0.664 & 0.196 & 0.842 & 0.293 & 0.696 & 0.336 \\',
'  SU-MCCDs  & 0.606 & 0.286 & 0.961 & 0.732 & 0.824 & 0.385 & 0.746 & 0.316 & 0.388 & 0.140 & 0.584 & 0.376 & 0.536 & 0.131 & 0.658 & 0.193 & 0.866 & 0.327 & 0.540 & 0.182 \\',
'  UN-MCCDs  & 0.686 & 0.446 & 0.930 & 0.600 & 0.493 & 0.116 & 0.749 & 0.318 & 0.474 & 0.036 & 0.648 & 0.381 & 0.591 & 0.146 & 0.771 & 0.263 & 0.815 & 0.261 & 0.519 & 0.118 \\',
'  SUN-MCCDs & 0.686 & 0.446 & 0.940 & 0.638 & 0.770 & 0.324 & 0.782 & 0.350 & 0.514 & 0.109 & 0.700 & 0.455 & 0.662 & 0.172 & 0.770 & 0.267 & 0.894 & 0.389 & 0.555 & 0.206 \\'
)

for ($i = 0; $i -lt 4; $i++) {
    $lines[$tprBlock[$i] - 1] = $newTpr[$i]
    $lines[$baBlock[$i]  - 1] = $newBa[$i]
}

[System.IO.File]::WriteAllLines($TR2, $lines)
"patched 8 rows"

"--- TPR/TNR block after patch ---"
1505..1508 | ForEach-Object { "{0}: {1}" -f $_, [System.IO.File]::ReadAllLines($TR2)[$_ - 1] }
"--- BA/F2 block after patch ---"
1527..1530 | ForEach-Object { "{0}: {1}" -f $_, [System.IO.File]::ReadAllLines($TR2)[$_ - 1] }
