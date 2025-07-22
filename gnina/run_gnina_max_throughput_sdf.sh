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
if ! command -v parallel &> /dev/null; then
    echo "ERROR: GNU Parallel is not installed. Please install it first."
    exit 1
fi

ACTUAL_GPUS=$(nvidia-smi -L | wc -l)
if [ "$NUM_GPUS" -gt "$ACTUAL_GPUS" ]; then
    echo "ERROR: Your configuration requests NUM_GPUS=$NUM_GPUS, but only $ACTUAL_GPUS GPUs were found."
    echo "Please set NUM_GPUS to $ACTUAL_GPUS or lower."
    exit 1
fi

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

# --- Create Master Ligand List ---
MASTER_LIST_FILE="all_ligands_to_process.txt"
find "$LIGAND_DIR" -name "*.sdf" > "$MASTER_LIST_FILE"

TOTAL_LIGANDS=$(wc -l < "$MASTER_LIST_FILE")
if [ "$TOTAL_LIGANDS" -eq 0 ]; then
    echo "ERROR: No .sdf files found in '$LIGAND_DIR'."
    rm "$MASTER_LIST_FILE"
    exit 1
fi

# --- Create processing function ---
process_ligand() {
    ligand_file="$1"
    GPU_ID=$(( (${PARALLEL_SEQ} - 1) % NUM_GPUS ))
    base=$(basename "$ligand_file" .sdf)
    output_file="$OUTPUT_DIR/${base}_scored.sdf"
    
    # Skip if output file already exists
    if [ -f "$output_file" ]; then
        echo "Skipping ligand $base - already scored (output file exists)"
        return 0
    fi
    
    echo "Processing ligand $base on GPU $GPU_ID (job ${PARALLEL_SEQ})"
    
    # Corrected ligand path for use inside the container
    ligand_path_in_container="/ligands/$(basename "$ligand_file")"

    # Add timeout to prevent hanging jobs
    timeout $MAX_CONTAINER_TIME docker run --rm --ipc=host --gpus "device=$GPU_ID" \
      --memory="$MEMORY_PER_CONTAINER" \
      --cpus="$CPUS_PER_CONTAINER" \
      -v "$(realpath $(pwd))":/work \
      -v "$(realpath "$LIGAND_DIR")":"/ligands" \
      -w /work \
      gnina/gnina:latest \
      gnina --score_only --cnn_scoring all \
            -r "$RECEPTOR_FILE" \
            -l "$ligand_path_in_container" \
            --autobox_ligand "$ligand_path_in_container" \
            -o "$output_file"
    
    # Check if the job completed successfully
    if [ $? -ne 0 ]; then
        echo "ERROR: Job for ligand $base failed or timed out"
        return 1
    fi
}

# --- Export variables and function so they are available to the subshells created by parallel ---
export LIGAND_DIR
export RECEPTOR_FILE
export OUTPUT_DIR
export NUM_GPUS
export MAX_CONTAINER_TIME
export MEMORY_PER_CONTAINER
export CPUS_PER_CONTAINER
export -f process_ligand

TOTAL_JOBS=$((NUM_GPUS * JOBS_PER_GPU))

echo "Receptor:      $RECEPTOR_FILE"
echo "Ligand Dir:    $LIGAND_DIR"
echo "Output Dir:    $OUTPUT_DIR"
echo "Total Ligands: $TOTAL_LIGANDS"
echo "GPUs to use:   $NUM_GPUS"
echo "Jobs per GPU:  $JOBS_PER_GPU"
echo "Total concurrent Docker containers: $TOTAL_JOBS"
echo "Memory per container: $MEMORY_PER_CONTAINER"
echo "CPUs per container: $CPUS_PER_CONTAINER"
echo "----------------------------------------------------"
echo "Starting parallel processing... Progress will be shown below."

# --- Run everything using GNU Parallel with periodic cleanup ---
echo "Starting processing with periodic cleanup every $CLEANUP_INTERVAL jobs..."

# Enhanced cleanup function
cleanup_docker() {
    echo "$(date): Running cleanup at job {#}..."
    
    # Check GPU temperatures before continuing
    if [ "$GPU_TEMP_CHECK" = "true" ] && command -v nvidia-smi &> /dev/null; then
        max_temp=$(nvidia-smi --query-gpu=temperature.gpu --format=csv,noheader,nounits | sort -nr | head -1)
        if [ "$max_temp" -gt 75 ]; then
            echo "WARNING: GPU temperature ($max_temp°C) is high. Pausing for 30 seconds..."
            sleep 30
        fi
    fi
    
    # Kill any containers running longer than timeout
    docker ps --format "{{.ID}} {{.RunningFor}}" | grep -E "(hour|[1-9][0-9]+ minute)" | awk '{print $1}' | while read container_id; do
        if [ ! -z "$container_id" ]; then
            echo "Killing long-running container: $container_id"
            docker kill "$container_id" 2>/dev/null || true
        fi
    done
    
    # System cleanup
    docker system prune -f > /dev/null 2>&1 || true
    
    # Brief pause to let system stabilize
    sleep 5
    
    echo "$(date): Cleanup completed"
}

# Export cleanup function
export -f cleanup_docker

# Run with periodic cleanup
cat "$MASTER_LIST_FILE" | parallel -j "$TOTAL_JOBS" --eta --joblog gnina_parallel.log \
  --line-buffer \
  'if (( {#} % '$CLEANUP_INTERVAL' == 0 )); then cleanup_docker; fi; process_ligand {}'

# --- Cleanup ---
rm "$MASTER_LIST_FILE"

echo "----------------------------------------------------"
echo "All jobs completed. Check '$OUTPUT_DIR' for results and 'gnina_parallel.log' for a detailed log."