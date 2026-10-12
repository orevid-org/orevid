# Optional, Windows only: helper used during the original runs; not needed to rebuild the project.
# Keeps the computer awake (without changing any setting) while TxGNN trains in WSL
# or an OREVID evaluation runs. Exits by itself when none is left.
#   powershell -ExecutionPolicy Bypass -File evaluation\keep_awake.ps1
Add-Type -Namespace Win32 -Name Power -MemberDefinition @'
[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);
'@
$ES_CONTINUOUS = [uint32]"0x80000000"
$ES_SYSTEM_REQUIRED = [uint32]"0x00000001"
[Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED) | Out-Null
"$(Get-Date -Format T) awake"
try {
    while ($true) {
        Start-Sleep -Seconds 300
        $wsl = (wsl -d Ubuntu -- pgrep -f txgnn_wsl.py) -ne $null
        $ours = Get-CimInstance Win32_Process -Filter "Name='python.exe' OR Name='py.exe' OR Name='powershell.exe'" |
                Where-Object { $_.CommandLine -match "evaluation\.|predictor\.py|web\.build|rebuild_all|txgnn_queue|OREVID-bench" }
        if (-not $wsl -and -not $ours) { break }
        [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED) | Out-Null
    }
} finally {
    [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS) | Out-Null
    "$(Get-Date -Format T) no longer needed: sleep allowed"
}
