#!/bin/bash

# === Configuration ===
# 1. FULL PATH to your compiled AutoDock-GPU executable
AUTODOCK_GPU_EXEC="/home/joe/projects/AutoDOCK/AutoDock-GPU/bin/autodock_gpu_128wi"

# 2. Name of your PRE-CALCULATED receptor grid map file (.maps.fld)
MAP_FILE_FLD_FILENAME="myreceptor_targeted.maps.fld" # Make sure this file exists

# 3. Name of the directory containing ALL your ligand PDBQT files
LIGAND_DIR="ligands_pdbqt" # Your main ligand directory

# 4. Name of the BASE directory where docking results will be saved
BASE_OUTPUT_DIR="docking_results_dual" # New output dir name

# 5. Docking run parameters
NUM_RUNS="10"               # Number of docking runs per ligand

# 6. GPU Device Numbers (AutoDock-GPU is 1-indexed)
GPU1_DEVNUM="1"
GPU2_DEVNUM="2"

# 7. Log directory
LOG_DIR="ad_gpu_logs_solA" # New log dir name
# === End of Configuration ===

CURRENT_WORK_DIR=$(pwd) # Script expects to be run from the folder containing the .maps.fld

# --- Sanity Checks ---
if [ ! -x "$AUTODOCK_GPU_EXEC" ]; then
    echo "ERROR: AutoDock-GPU binary not found or not executable at: $AUTODOCK_GPU_EXEC"
    exit 1
fi

MAP_FILE_ABS="$CURRENT_WORK_DIR/$MAP_FILE_FLD_FILENAME" # Absolute path to map file
if [ ! -f "$MAP_FILE_ABS" ]; then
    echo "ERROR: Receptor map file not found: $MAP_FILE_ABS"
    exit 1
fi

LIGAND_DIR_ABS="$CURRENT_WORK_DIR/$LIGAND_DIR"
if [ ! -d "$LIGAND_DIR_ABS" ]; then
    echo "ERROR: Ligand directory not found: $LIGAND_DIR_ABS"
    exit 1
fi

# --- Prepare Ligand Lists for Each GPU (with MAP FILE AS FIRST LINE) ---
echo "Preparing ligand lists for 2 GPUs (map file will be the first line in each batch file)..."
ALL_LIGANDS_FILE_TMP="${CURRENT_WORK_DIR}/all_ligands_abs.tmp"

# Generate list of absolute paths to ligand files
find "$LIGAND_DIR_ABS" -maxdepth 1 -name "*.pdbqt" -type f -print0 | xargs -0 realpath > "$ALL_LIGANDS_FILE_TMP"

TOTAL_LIGANDS=$(wc -l < "$ALL_LIGANDS_FILE_TMP")

if [ "$TOTAL_LIGANDS" -eq 0 ]; then
    echo "ERROR: No PDBQT files found in '$LIGAND_DIR_ABS'."
    rm "$ALL_LIGANDS_FILE_TMP"
    exit 1
fi
echo "Found $TOTAL_LIGANDS total ligands."

LIGANDS_PER_GPU=$(( (TOTAL_LIGANDS + 1) / 2 )) # Ensure all ligands are covered if odd number

BATCH_FILE_GPU1="${CURRENT_WORK_DIR}/ligand_batch_gpu1_map_inc.txt"
BATCH_FILE_GPU2="${CURRENT_WORK_DIR}/ligand_batch_gpu2_map_inc.txt"

# Create batch file for GPU 1 (ONLY ligands, NO map file)
head -n "$LIGANDS_PER_GPU" "$ALL_LIGANDS_FILE_TMP" > "$BATCH_FILE_GPU1"

# Create batch file for GPU 2 (ONLY ligands, NO map file)
tail -n "$((TOTAL_LIGANDS - LIGANDS_PER_GPU))" "$ALL_LIGANDS_FILE_TMP" > "$BATCH_FILE_GPU2"

rm "$ALL_LIGANDS_FILE_TMP" # Clean up temporary file

echo "Batch file for GPU $GPU1_DEVNUM: $BATCH_FILE_GPU1 ($(wc -l < "$BATCH_FILE_GPU1") ligand files)"
echo "Batch file for GPU $GPU2_DEVNUM: $BATCH_FILE_GPU2 ($(wc -l < "$BATCH_FILE_GPU2") ligand files)"
# Optional: Check content of batch files
# echo "--- Content of $BATCH_FILE_GPU1 (first 5 lines) ---"
# head -n 5 "$BATCH_FILE_GPU1"
# echo "--------------------------------------------------"
echo "-------------------------------------------"

# --- Create Output and Log Directories ---
OUTPUT_DIR_GPU1="${CURRENT_WORK_DIR}/${BASE_OUTPUT_DIR}/gpu${GPU1_DEVNUM}_results"
OUTPUT_DIR_GPU2="${CURRENT_WORK_DIR}/${BASE_OUTPUT_DIR}/gpu${GPU2_DEVNUM}_results"
LOG_DIR_ABS="${CURRENT_WORK_DIR}/${LOG_DIR}"

mkdir -p "$OUTPUT_DIR_GPU1"
mkdir -p "$OUTPUT_DIR_GPU2"
mkdir -p "$LOG_DIR_ABS"

echo "Results for GPU $GPU1_DEVNUM will be in: $OUTPUT_DIR_GPU1"
echo "Results for GPU $GPU2_DEVNUM will be in: $OUTPUT_DIR_GPU2"
echo "Logs will be in: $LOG_DIR_ABS"
echo "-------------------------------------------"

# --- Launch Docking Jobs in Parallel ---
# NOTE: Map file must be specified with --ffile, batch files contain ONLY ligands

RESNAM_GPU1="${OUTPUT_DIR_GPU1}/"
RESNAM_GPU2="${OUTPUT_DIR_GPU2}/"

echo "Starting docking on GPU $GPU1_DEVNUM..."
echo "Command: \"$AUTODOCK_GPU_EXEC\" --ffile \"$MAP_FILE_ABS\" -B \"$BATCH_FILE_GPU1\" --nrun \"$NUM_RUNS\" -resnam \"$RESNAM_GPU1\" --devnum \"$GPU1_DEVNUM\""
"$AUTODOCK_GPU_EXEC" \
    --ffile "$MAP_FILE_ABS" \
    -B "$BATCH_FILE_GPU1" \
    --nrun "$NUM_RUNS" \
    -resnam "$RESNAM_GPU1" \
    --devnum "$GPU1_DEVNUM" > "${LOG_DIR_ABS}/autodock_gpu1.log" 2>&1 &
PID_GPU1=$!
echo "AutoDock-GPU process for GPU $GPU1_DEVNUM started with PID $PID_GPU1. Log: ${LOG_DIR_ABS}/autodock_gpu1.log"

echo "Starting docking on GPU $GPU2_DEVNUM..."
echo "Command: \"$AUTODOCK_GPU_EXEC\" --ffile \"$MAP_FILE_ABS\" -B \"$BATCH_FILE_GPU2\" --nrun \"$NUM_RUNS\" -resnam \"$RESNAM_GPU2\" --devnum \"$GPU2_DEVNUM\""
"$AUTODOCK_GPU_EXEC" \
    --ffile "$MAP_FILE_ABS" \
    -B "$BATCH_FILE_GPU2" \
    --nrun "$NUM_RUNS" \
    -resnam "$RESNAM_GPU2" \
    --devnum "$GPU2_DEVNUM" > "${LOG_DIR_ABS}/autodock_gpu2.log" 2>&1 &
PID_GPU2=$!
echo "AutoDock-GPU process for GPU $GPU2_DEVNUM started with PID $PID_GPU2. Log: ${LOG_DIR_ABS}/autodock_gpu2.log"

echo "-------------------------------------------"
echo "Both AutoDock-GPU processes launched. Waiting for completion..."

wait $PID_GPU1
EXIT_STATUS_GPU1=$?
echo "GPU $GPU1_DEVNUM process (PID $PID_GPU1) finished with exit status: $EXIT_STATUS_GPU1"

wait $PID_GPU2
EXIT_STATUS_GPU2=$?
echo "GPU $GPU2_DEVNUM process (PID $PID_GPU2) finished with exit status: $EXIT_STATUS_GPU2"

echo "-------------------------------------------"
if [ "$EXIT_STATUS_GPU1" -eq 0 ] && [ "$EXIT_STATUS_GPU2" -eq 0 ]; then
    echo "All docking jobs seem to have completed. Check logs and output directories."
else
    echo "One or both docking jobs may have failed. Check logs: ${LOG_DIR_ABS}/autodock_gpu1.log and ${LOG_DIR_ABS}/autodock_gpu2.log"
fi

echo "Batch docking complete."
# Keep batch files for inspection if needed
# echo "Cleaning up batch files: $BATCH_FILE_GPU1, $BATCH_FILE_GPU2"
# rm "$BATCH_FILE_GPU1" "$BATCH_FILE_GPU2"