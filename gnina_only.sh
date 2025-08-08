#!/bin/bash

# add protein & library paths & output paths here
protein_path="jobs/benchmarking/ADRB2/4lde_ADRB2.pdb" # your protein 
output_path="jobs/benchmarking/ADRB2/4lde_ADRB2_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_ADRB2_5P/" # path to the sdf files to be docked
smiles_file="diffdock/ADRB2_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)
x_min=-9.5548
y_min=-23.6649
z_min=-58.0429
x_max=5.4960
y_max=-5.9767
z_max=-44.1756
# =======================================================================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi

# gnina container settings:
NUM_GPUS=2
JOBS_PER_GPU=19
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8

# this means 19 CPU * 0.8 = 15.2 CPU cores per GPU so you need at least 31.2 CPU threads for 2 GPUs
# other settings could be if you have 24 Threads and one gpu:
# NUM_GPUS=1
# JOBS_PER_GPU=28
# which means 28 * 0.8 = 22.4 CPU cores per GPU so you need at least 22.4 CPU threads for 1 GPU

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate diffdock

echo "Debug: Checking paths..."
echo "SDF files dir: $sdf_files"
echo "Output path: $output_path"
echo "Protein path: $protein_path"
echo "SMILES file: $smiles_file"

# python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
# python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"

