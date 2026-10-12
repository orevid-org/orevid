# Optional, Windows + WSL only: helper used during the original runs; not needed to rebuild the project.
# Queue of TxGNN trainings in WSL (see docs/TXGNN_EVALUATION.md). If WSL crashes (e.g. out of
# memory), it restarts it and starts again: finished splits are skipped and already trained
# models are loaded. Up to 6 attempts.
#   powershell -ExecutionPolicy Bypass -File evaluation\txgnn_queue.ps1 > data\eval\txgnn_queue.txt
# (evaluation/txgnn_progress.py reads that file to show the state of the queue.)
$queue = "complex_disease 1 mental_health 1 diabetes 1 neurodigenerative 1 cardiovascular 1 anemia 1 autoimmune 1 adrenal_gland 1 metabolic_disorder 1 cell_proliferation 1 complex_disease 2 complex_disease 3 complex_disease 4 complex_disease 5"
$parts = $queue.Split(" ")
$repo = Split-Path -Parent $PSScriptRoot
$pred = Join-Path $repo "data\eval\txgnn_predictions"
$script = (wsl -d Ubuntu -- wslpath -a ($repo + "\evaluation\train_txgnn.sh")).Trim()  # same file, WSL path
for ($attempt = 1; $attempt -le 6; $attempt++) {
    "$(Get-Date -Format 'dd HH:mm') attempt $attempt"
    wsl -d Ubuntu -- bash -c "sed -i 's/\r$//' $script"
    wsl -d Ubuntu -- bash $script @parts
    $missing = 0
    for ($i = 0; $i -lt $parts.Length; $i += 2) {
        if (-not (Test-Path "$pred\TxGNN_$($parts[$i+1])_$($parts[$i]).parquet")) { $missing++ }
    }
    if ($missing -eq 0) { "$(Get-Date -Format 'dd HH:mm') all done"; break }
    "$(Get-Date -Format 'dd HH:mm') ${missing} missing; restarting WSL"
    wsl --shutdown
    Start-Sleep -Seconds 20
}
