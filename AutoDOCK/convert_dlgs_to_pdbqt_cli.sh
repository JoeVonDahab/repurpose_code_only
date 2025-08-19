#!/bin/bash

# === Default Configuration (overridable via CLI) ===
PYTHON_EXEC="/home/joe/miniconda3/envs/autodock/bin/python"
MGLTOOLS_ROOT="/home/joe/miniconda3/envs/autodock/MGLToolsPckgs"
SOURCE_DIR="docking_results_compiled"
DEST_DIR="docked_converted"
THREADS="${THREADS:-32}"               # <<< how many processes you want
# === End Defaults ===

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -p <path>   Python executable in autodock conda environment
  -m <path>   MGLTools root directory
  -s <dir>    Source directory containing .dlg files
  -d <dir>    Destination directory for .pdbqt files
  -t <num>    Number of parallel threads (default: 32)
  -h          Show this help

Example:
  $(basename "$0") -s results/gpu1_results -d converted_poses -p /path/to/python -t 16
EOF
}

# --- Parse CLI Args ---
while getopts ":p:m:s:d:t:h" opt; do
  case $opt in
    p) PYTHON_EXEC="$OPTARG" ;;
    m) MGLTOOLS_ROOT="$OPTARG" ;;
    s) SOURCE_DIR="$OPTARG" ;;
    d) DEST_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Missing value for -$OPTARG"; usage; exit 1 ;;
    \?) echo "Unknown option -$OPTARG"; usage; exit 1 ;;
  esac
done

CONVERT_SCRIPT="${MGLTOOLS_ROOT}/AutoDockTools/Utilities24/write_lowest_energy_ligand.py"

set -euo pipefail
echo "Starting DLG → PDBQT conversion with $THREADS parallel workers…"
echo "Source DLG directory: $SOURCE_DIR"
echo "Destination PDBQT directory: $DEST_DIR"
echo "Using MGLTools script: $CONVERT_SCRIPT"
echo "-------------------------------------------"

export MGLTOOLS="$MGLTOOLS_ROOT"
export PYTHONPATH="${MGLTOOLS}:${PYTHONPATH:-}"

# --- Validation ---
for p in "$PYTHON_EXEC" "$CONVERT_SCRIPT"; do
  [[ -x $p || -f $p ]] || { echo "ERROR: missing $p"; exit 1; }
done

if [ ! -d "$SOURCE_DIR" ]; then
    echo "ERROR: Source directory '$SOURCE_DIR' not found."
    exit 1
fi

mkdir -p "$DEST_DIR"
if [ ! -d "$DEST_DIR" ]; then
    echo "ERROR: Could not create destination directory '$DEST_DIR'."
    exit 1
fi

# -- build file list ---------------------------------------------------------
mapfile -d '' DLG_FILES < <(find "$SOURCE_DIR" -maxdepth 1 -type f -name "*.dlg" -print0)
[[ ${#DLG_FILES[@]} -gt 0 ]] || { echo "No .dlg files found"; exit 0; }
echo "Found ${#DLG_FILES[@]} .dlg files."

# -- helper used by each worker ---------------------------------------------
convert_one() {
  local dlg="$1"
  local base
  base=$(basename "${dlg%.dlg}")
  local out="$DEST_DIR/${base}_best_pose.pdbqt"

  echo "[PID $$] $(basename "$dlg") → $(basename "$out")"
  "$PYTHON_EXEC" "$CONVERT_SCRIPT" -f "$dlg" -o "$out" \
    && [[ -f $out ]] \
    && echo "  ✔︎ done" \
    || echo "  ✖︎ failed"
}
export -f convert_one
export PYTHON_EXEC CONVERT_SCRIPT DEST_DIR  # needed inside subshells

# -- run in parallel ---------------------------------------------------------
if command -v parallel &>/dev/null; then
  printf '%s\0' "${DLG_FILES[@]}" | \
    parallel -0 -P "$THREADS" convert_one {}
else
  # fallback using xargs if GNU parallel isn't installed
  printf '%s\0' "${DLG_FILES[@]}" | \
    xargs -0 -n1 -P "$THREADS" bash -c 'convert_one "$1"' _
fi

echo "All conversions finished; output in $DEST_DIR"