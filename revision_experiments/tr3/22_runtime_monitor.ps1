# revision_experiments/tr3/22_runtime_monitor.ps1
#
# Load monitor for 22_runtime.R. Started by the driver in the background; it
# samples every -Every seconds until the driver process ends and writes one
# line per sample to -Out:
#   time, foreign_r (number of R processes that do not belong to the driver),
#   foreign_ids, sys_cpu (total CPU use of the machine, percent)
# "Own" processes are the driver (R.exe), its parent (Rscript.exe) and all their
# descendants (the PSOCK workers). Any other R-like process counts as foreign.
param(
  [Parameter(Mandatory = $true)] [int] $DriverPid,
  [Parameter(Mandatory = $true)] [string] $Out,
  [int] $Every = 3
)
$rlike = '^(Rscript|Rterm|R|Rgui|rsession|rsession-utf8)\.exe$'
# Rscript.exe starts the R binary through cmd.exe (bin\Rscript.exe -> cmd.exe -> bind\Rscript.exe)
$launcher = '^(Rscript|Rterm|R|Rgui|cmd)\.exe$'
"time,foreign_r,foreign_ids,sys_cpu" | Set-Content -Path $Out -Encoding ascii
while ($true) {
  if (-not (Get-Process -Id $DriverPid -ErrorAction SilentlyContinue)) { break }
  $all = Get-CimInstance Win32_Process
  $byId = @{}
  foreach ($p in $all) { $byId[[int]$p.ProcessId] = $p }
  # the driver's root: walk up through the launcher chain (Rscript.exe -> cmd.exe -> Rscript.exe)
  $root = $DriverPid
  while ($byId.ContainsKey($root)) {
    $par = [int]$byId[$root].ParentProcessId
    if ($byId.ContainsKey($par) -and $byId[$par].Name -match $launcher) { $root = $par } else { break }
  }
  $own = New-Object 'System.Collections.Generic.HashSet[int]'
  [void]$own.Add($root)
  $grew = $true
  while ($grew) {
    $grew = $false
    foreach ($p in $all) {
      if ($own.Contains([int]$p.ParentProcessId) -and -not $own.Contains([int]$p.ProcessId)) {
        [void]$own.Add([int]$p.ProcessId); $grew = $true
      }
    }
  }
  $foreign = @($all | Where-Object { $_.Name -match $rlike -and -not $own.Contains([int]$_.ProcessId) })
  $cpu = (Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average
  "{0},{1},{2},{3}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $foreign.Count, (($foreign | ForEach-Object { $_.ProcessId }) -join ';'), $cpu |
    Add-Content -Path $Out -Encoding ascii
  Start-Sleep -Seconds $Every
}
