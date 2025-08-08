#!/bin/bash

# === Default Configuration ===
NUM_GPUS=2
JOBS_PER_GPU=19
LIGAND_DIR=""
RECEPTOR_FILE=""
OUTPUT_DIR="scored"
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER="8g"
CPUS_PER_CONTAINER="0.8"

# === Parse Command Line Arguments ===
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    --receptor_file PATH        Path to the prepared receptor file (required)
    --ligand_files PATH         Directory containing ligand .sdf files (required)
    --output_dir PATH           Directory where scored output files will be saved (default: scored)
    --num_gpus NUMBER           Number of GPUs available (default: 2)
    --jobs_per_gpu NUMBER       Concurrent gnina processes per GPU (default: 19)
    --cleanup_interval NUMBER   Cleanup interval to prevent crashes (default: 500)
    --max_container_time NUMBER Timeout per job in seconds (default: 300)
    --gpu_temp_check BOOL       Enable GPU temperature checking (default: true)
    --memory_per_container SIZE Memory limit per container (default: 8g)
    --cpus_per_container NUMBER CPU limit per container (default: 0.8)
    --help                      Show this help message

Example:
    $0 --receptor_file protein.pdb --ligand_files /path/to/ligands --output_dir results --num_gpus 2
EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --receptor_file)
            RECEPTOR_FILE="$2"
            shift 2
            ;;
        --ligand_files)
            LIGAND_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --num_gpus)
            NUM_GPUS="$2"
            shift 2
            ;;
        --jobs_per_gpu)
            JOBS_PER_GPU="$2"
            shift 2
            ;;
        --cleanup_interval)
            CLEANUP_INTERVAL="$2"
            shift 2
            ;;
        --max_container_time)
            MAX_CONTAINER_TIME="$2"
            shift 2
            ;;
        --gpu_temp_check)
            GPU_TEMP_CHECK="$2"
            shift 2
            ;;
        --memory_per_container)
            MEMORY_PER_CONTAINER="$2"
            shift 2
            ;;
        --cpus_per_container)
            CPUS_PER_CONTAINER="$2"
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# === Validate Required Arguments ===
if [ -z "$RECEPTOR_FILE" ]; then
    echo "ERROR: --receptor_file is required"
    show_help
    exit 1
fi

if [ -z "$LIGAND_DIR" ]; then
    echo "ERROR: --ligand_files is required"
    show_help
    exit 1
fi

# --- Sanity Checks & Preparation ---
ACTUAL_GPUS=$(nvidia-smi -L | wc -l)
if [ "$NUM_GPUS" -gt "$ACTUAL_GPUS" ]; then
    echo "ERROR: Your configuration requests NUM_GPUS=$NUM_GPUS, but only $ACTUAL_GPUS GPUs were found."
    echo "Please set NUM_GPUS to $ACTUAL_GPUS or lower."
    exit 1
fi

# Detect host CPUs and allocate cores per container
HOST_CPUS=$(nproc)
CPUS_PER_CONTAINER=$(( HOST_CPUS / NUM_GPUS ))

# fall back to 1 if arithmetic ever gives zero
if [ "$CPUS_PER_CONTAINER" -lt 1 ]; then
    CPUS_PER_CONTAINER=1
fi

echo "Host CPUs: $HOST_CPUS"
echo "Allocating $CPUS_PER_CONTAINER CPU cores per GPU container"

if [ ! -d "$LIGAND_DIR" ]; then
    echo "ERROR: Ligand directory not found: $LIGAND_DIR"
    exit 1
fi

if [ ! -f "$RECEPTOR_FILE" ]; then
    echo "ERROR: Receptor file not found: $RECEPTOR_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Pre-pull Docker image to avoid repeated downloads ---
echo "Ensuring Docker image is available..."
docker pull gnina/gnina:latest
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to pull Docker image gnina/gnina:latest"
    exit 1
fi

# --- Set up crash prevention ---
# Increase file descriptor limits
ulimit -n 65536

# Set up memory monitoring
MEMORY_CHECK_INTERVAL=500
PROCESSED_COUNT=0

# --- Create Master Ligand List & split into GPU chunks ---
MASTER_LIST_FILE="all_ligands_to_process.txt"
find "$LIGAND_DIR" -name "*.sdf" > "$MASTER_LIST_FILE"

TOTAL_LIGANDS=$(wc -l < "$MASTER_LIST_FILE")
if [ "$TOTAL_LIGANDS" -eq 0 ]; then
    echo "ERROR: No .sdf files found in '$LIGAND_DIR'."
    rm "$MASTER_LIST_FILE"
    exit 1
fi

echo "Total ligands: $TOTAL_LIGANDS → splitting into $NUM_GPUS parts"
# split into N roughly-equal parts: ligands_part_aa, ligands_part_ab, ...
split -n l/$NUM_GPUS "$MASTER_LIST_FILE" ligands_part_

TOTAL_JOBS=$((NUM_GPUS * JOBS_PER_GPU))

echo "Receptor:      $RECEPTOR_FILE"
echo "Ligand Dir:    $LIGAND_DIR"
echo "Output Dir:    $OUTPUT_DIR"
echo "Total Ligands: $TOTAL_LIGANDS"
echo "GPUs to use:   $NUM_GPUS"
echo "Jobs per GPU:  $JOBS_PER_GPU"
echo "Total concurrent gnina processes: $TOTAL_JOBS"
echo "Memory per container: $MEMORY_PER_CONTAINER"
echo "CPUs per container: $CPUS_PER_CONTAINER"
echo "----------------------------------------------------"

# --- Launch one worker container per GPU ---
echo "Launching $NUM_GPUS worker containers (one per GPU)..."
i=0
for part in ligands_part_*; do
    GPU_ID=$(( i % NUM_GPUS ))
    echo "  • GPU $GPU_ID → processing $(wc -l < "$part") ligands (file: $part)"

    docker run --rm --ipc=host --gpus "device=$GPU_ID" \
        --memory="$MEMORY_PER_CONTAINER" \
        --cpus="$CPUS_PER_CONTAINER" \
        --log-driver=none \
        -v "$(realpath "$(pwd)")":/work \
        -v "$(realpath "$LIGAND_DIR")":/ligands \
        -v "$(dirname "$(realpath "$RECEPTOR_FILE")")":/receptor \
        -w /work \
        gnina/gnina:latest bash -c '
            set -euo pipefail
            # Allow GNINA to use all cores in this container
            export OMP_NUM_THREADS='"$CPUS_PER_CONTAINER"'
            
            REC="/receptor/'"$(basename "$RECEPTOR_FILE")"'"
            OUTDIR="'"$OUTPUT_DIR"'"
            mkdir -p "$OUTDIR"

            # Create a simple processing function
            process_ligand() {
                ligand_path="$1"
                base=$(basename "$ligand_path" .sdf)
                out="$OUTDIR/${base}_scored.sdf"
                
                if [ -f "$out" ]; then
                    echo "[skip] $base"
                    return 0
                fi
                
                echo "[GPU '"$GPU_ID"'] scoring $base... (threads=$OMP_NUM_THREADS)"
                timeout -s TERM -k 15s '"$MAX_CONTAINER_TIME"'s \
                    gnina --score_only --cnn_scoring all \
                          -r "$REC" \
                          -l "/ligands/$(basename "$ligand_path")" \
                          --autobox_ligand "/ligands/$(basename "$ligand_path")" \
                          -o "$out"
            }
            
            export -f process_ligand
            export REC OUTDIR OMP_NUM_THREADS
            
            # Use xargs to process ligands in parallel
            cat '"$part"' | xargs -n1 -P '"$JOBS_PER_GPU"' -I {} bash -c "process_ligand {}"
        ' &

    ((i++))
done

# wait for all GPU workers to finish
wait

# --- Cleanup temporary split files ---
rm ligands_part_* "$MASTER_LIST_FILE"

echo "----------------------------------------------------"
echo "All jobs completed. Check '$OUTPUT_DIR' for results."