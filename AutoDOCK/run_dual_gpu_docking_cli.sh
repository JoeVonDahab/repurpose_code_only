#!/bin/bash
# IMPORTANT: Ensure this file uses LF line endings (run: dos2unix AutoDOCK/run_dual_gpu_docking_cli.sh)
set -e
set -u
if (set -o | grep -q pipefail) 2>/dev/null; then
  set -o pipefail 2>/dev/null || true
fi

# === Default Configuration (overridable via CLI) ===
AUTODOCK_GPU_EXEC="/home/joe/projects/AutoDOCK/AutoDock-GPU/bin/autodock_gpu_128wi"
MAP_FILE_FLD_FILENAME="myreceptor_targeted.maps.fld"
LIGAND_DIR="ligands_pdbqt"
BASE_OUTPUT_DIR="docking_results_dual"
NUM_RUNS="10"
GPU_DEVICES_CSV="1"
LOG_DIR="ad_gpu_logs_solA"
# === End Defaults ===

usage() {
  echo "Usage: $(basename "$0") [options]"
  echo "  -e <path>   AutoDock-GPU executable"
  echo "  -m <file>   Receptor map .maps.fld"
  echo "  -l <dir>    Ligand directory (*.pdbqt)"
  echo "  -o <dir>    Base output directory"
  echo "  -n <int>    Runs per ligand"
  echo "  -g <list>   Comma-separated GPU device numbers (1-indexed)"
  echo "  -L <dir>    Log directory"
  echo "  -h          Help"
}

# --- Parse CLI Args ---
while getopts ":e:m:l:o:n:g:L:h" opt; do
  case $opt in
    e) AUTODOCK_GPU_EXEC="$OPTARG" ;;
    m) MAP_FILE_FLD_FILENAME="$OPTARG" ;;
    l) LIGAND_DIR="$OPTARG" ;;
    o) BASE_OUTPUT_DIR="$OPTARG" ;;
    n) NUM_RUNS="$OPTARG" ;;
    g) GPU_DEVICES_CSV="$OPTARG" ;;
    L) LOG_DIR="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Missing value for -$OPTARG"; usage; exit 1 ;;
    \?) echo "Unknown option -$OPTARG"; usage; exit 1 ;;
  esac
done

CURRENT_WORK_DIR=$(pwd)

# --- Sanity Checks (unchanged logic) ---
if [ ! -x "$AUTODOCK_GPU_EXEC" ]; then
  echo "ERROR: AutoDock-GPU binary not executable: $AUTODOCK_GPU_EXEC" >&2
  exit 1
fi

# Test AutoDock-GPU executable briefly
echo "Testing AutoDock-GPU executable..."
if ! "$AUTODOCK_GPU_EXEC" --help >/dev/null 2>&1; then
  echo "WARNING: AutoDock-GPU executable may have issues (--help failed)" >&2
fi

# Check GPU availability if nvidia-smi is available
if command -v nvidia-smi >/dev/null 2>&1; then
  echo "GPU status:"
  nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo "nvidia-smi query failed"
fi

# --- Input File Validation ---
if [[ "$MAP_FILE_FLD_FILENAME" = /* ]]; then MAP_FILE_ABS="$MAP_FILE_FLD_FILENAME"; else MAP_FILE_ABS="$CURRENT_WORK_DIR/$MAP_FILE_FLD_FILENAME"; fi
[ -f "$MAP_FILE_ABS" ] || { echo "ERROR: Map file missing: $MAP_FILE_ABS" >&2; exit 1; }
if [[ "$LIGAND_DIR" = /* ]]; then LIGAND_DIR_ABS="$LIGAND_DIR"; else LIGAND_DIR_ABS="$CURRENT_WORK_DIR/$LIGAND_DIR"; fi
[ -d "$LIGAND_DIR_ABS" ] || { echo "ERROR: Ligand directory missing: $LIGAND_DIR_ABS" >&2; exit 1; }

IFS=',' read -r -a GPU_DEVICES <<< "$GPU_DEVICES_CSV"
[ "${#GPU_DEVICES[@]}" -gt 0 ] || { echo "ERROR: No GPUs specified." >&2; exit 1; }
echo "GPUs: ${GPU_DEVICES[*]}"

echo "Collecting ligands..."
ALL_LIGANDS_FILE_TMP="${CURRENT_WORK_DIR}/all_ligands_abs.tmp"
find "$LIGAND_DIR_ABS" -maxdepth 1 -type f -name "*.pdbqt" -print0 | xargs -0 realpath > "$ALL_LIGANDS_FILE_TMP" || true
TOTAL_LIGANDS=$(wc -l < "$ALL_LIGANDS_FILE_TMP" || echo 0)
if [ "$TOTAL_LIGANDS" -eq 0 ]; then
  echo "ERROR: No .pdbqt ligands found in $LIGAND_DIR_ABS" >&2
  rm -f "$ALL_LIGANDS_FILE_TMP"
  exit 1
fi
echo "Found $TOTAL_LIGANDS ligands."
mapfile -t LIGANDS < "$ALL_LIGANDS_FILE_TMP"
rm -f "$ALL_LIGANDS_FILE_TMP"

# Output / logs
if [[ "$BASE_OUTPUT_DIR" = /* ]]; then BASE_OUTPUT_DIR_ABS="$BASE_OUTPUT_DIR"; else BASE_OUTPUT_DIR_ABS="$CURRENT_WORK_DIR/$BASE_OUTPUT_DIR"; fi
if [[ "$LOG_DIR" = /* ]]; then LOG_DIR_ABS="$LOG_DIR"; else LOG_DIR_ABS="$CURRENT_WORK_DIR/$LOG_DIR"; fi
mkdir -p "$LOG_DIR_ABS"

NUM_GPUS=${#GPU_DEVICES[@]}
BASE_COUNT=$(( TOTAL_LIGANDS / NUM_GPUS ))
REMAINDER=$(( TOTAL_LIGANDS % NUM_GPUS ))
echo "Distributing: base=$BASE_COUNT remainder=$REMAINDER"

declare -a PIDS
declare -A GPU_BATCH_FILES

START_INDEX=0
for idx in "${!GPU_DEVICES[@]}"; do
  DEV="${GPU_DEVICES[$idx]}"
  COUNT=$BASE_COUNT
  if [ "$idx" -lt "$REMAINDER" ]; then
    COUNT=$((COUNT + 1))
  fi
  if [ "$COUNT" -eq 0 ]; then
    echo "Warning: GPU $DEV gets 0 ligands."
    GPU_BATCH_FILES["$DEV"]=""
    continue
  fi
  BATCH_FILE="${CURRENT_WORK_DIR}/ligand_batch_gpu${DEV}.txt"
  
  # Create batch file with map file as first line, then ligands
  echo "$MAP_FILE_ABS" > "$BATCH_FILE"
  printf "%s\n" "${LIGANDS[@]:START_INDEX:COUNT}" >> "$BATCH_FILE"
  
  GPU_BATCH_FILES["$DEV"]="$BATCH_FILE"
  echo "GPU $DEV: $COUNT ligands -> $BATCH_FILE (with map file)"
  START_INDEX=$((START_INDEX + COUNT))
done

echo "-------------------------------------------"
echo "Launching docking jobs..."

for DEV in "${GPU_DEVICES[@]}"; do
  BATCH_FILE="${GPU_BATCH_FILES[$DEV]}"
  [ -n "$BATCH_FILE" ] || continue
  OUT_DIR_GPU="${BASE_OUTPUT_DIR_ABS}/gpu${DEV}_results"
  mkdir -p "$OUT_DIR_GPU"
  RESNAM="${OUT_DIR_GPU}/"
  LOG_FILE="${LOG_DIR_ABS}/autodock_gpu${DEV}.log"
  
  echo "[GPU $DEV] Starting docking..."
  echo "  Batch file: $BATCH_FILE ($(wc -l < "$BATCH_FILE") lines: 1 map + $(($(wc -l < "$BATCH_FILE") - 1)) ligands)"
  echo "  Output: $OUT_DIR_GPU"
  echo "  Log: $LOG_FILE"
  
  # Use -B for batch file (don't specify --ffile separately since map is in batch)
  echo "  Command: $AUTODOCK_GPU_EXEC -B $BATCH_FILE --nrun $NUM_RUNS -resnam $RESNAM --devnum $DEV"
  
  "$AUTODOCK_GPU_EXEC" \
      -B "$BATCH_FILE" \
      --nrun "$NUM_RUNS" \
      -resnam "$RESNAM" \
      --devnum "$DEV" > "$LOG_FILE" 2>&1 &
  PIDS+=("$!:GPU$DEV")
done

echo "-------------------------------------------"
echo "Waiting for ${#PIDS[@]} job(s)..."
OVERALL_OK=1
for entry in "${PIDS[@]}"; do
  PID="${entry%%:*}"
  TAG="${entry##*:}"
  if wait "$PID"; then
    echo "$TAG (PID $PID) OK"
  else
    rc=$?
    echo "$TAG (PID $PID) FAILED (exit $rc)"
    LOG_FILE="${LOG_DIR_ABS}/autodock_${TAG,,}.log"
    if [ -f "$LOG_FILE" ]; then
      echo "  Last 10 lines of $LOG_FILE:"
      tail -10 "$LOG_FILE" | sed 's/^/    /'
    fi
    OVERALL_OK=0
  fi
done

echo "-------------------------------------------"
if [ "$OVERALL_OK" -eq 1 ]; then
  echo "All docking jobs finished successfully."
else
  echo "One or more docking jobs failed. See logs: $LOG_DIR_ABS"
fi
echo "Done."