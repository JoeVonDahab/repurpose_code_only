#!/bin/bash

# add protein & library paths & output paths here

# KAT2A

protein_path="jobs/benchmarking/KAT2A/5h84_KAT2A.pdb" # your protein
output_path="jobs/benchmarking/KAT2A/5h84_KAT2A_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_KAT2A_5P/" # path to the sdf files to be docked
smiles_file="diffdock/KAT2A_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)
x_min=-11.35946
y_min=-2.50477
z_min=-20.51286
x_max=9.13694
y_max=15.97513
z_max=1.74724
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

python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================



# SECOND PDB 
# MAPK1

#!/bin/bash

# add protein & library paths & output paths here

# FEN1
protein_path="jobs/benchmarking/MAPK1/5v62_MAPK1.pdb" # your protein
output_path="jobs/benchmarking/MAPK1/5v62_MAPK1_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_MAPK1_5P/" # path to the sdf files to be docked
smiles_file="diffdock/MAPK1_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)
x_min=3.56820
y_min=12.89045
z_min=5.48755
x_max=20.50200
y_max=27.00155
z_max=20.60785

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

python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"



# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================



# THIRD PDB 
# MTORC1

#!/bin/bash

# add protein & library paths & output paths here

# MTORC1
protein_path="jobs/benchmarking/MTORC1/4jt5_MTORC1.pdb" # your protein 
output_path="jobs/benchmarking/MTORC1/4jt5_MTORC1_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_MTORC1_5P/" # path to the sdf files to be docked
smiles_file="diffdock/MTORC1_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)
x_min=40.63890
y_min=-10.11256
z_min=-55.66120
x_max=61.45530
y_max=9.06774
z_max=-41.37440

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

python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"


# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================



# FORTH PDB 
# PKM2

#!/bin/bash

# add protein & library paths & output paths here

# PKM2
protein_path="jobs/benchmarking/PKM2/3gr4_PKM2.pdb" # your protein 
output_path="jobs/benchmarking/PKM2/3gr4_PKM2_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_PKM2_5P/" # path to the sdf files to be docked
smiles_file="diffdock/PKM2_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)

x_min=-0.55805
y_min=-6.35544
z_min=5.11660
x_max=22.65825
y_max=10.70176
z_max=20.89540

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

python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================
# ============================ Cleanup ==================================



# FIFTH PDB 
# VDR

#!/bin/bash

# add protein & library paths & output paths here

# VDR
protein_path="jobs/benchmarking/VDR/3a2i_VDR.pdb" # your protein 
output_path="jobs/benchmarking/VDR/3a2i_VDR_output" # output path for the results
sdf_files="diffdock/ligand_sdf_files_VDR_5P/" # path to the sdf files to be docked
smiles_file="diffdock/VDR_5P.smi" # path to the smiles file of the library

# =======================================================================================
best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# ============================ Box Filter Coordinates ===================================
# box filter coordinates (optional)

x_min=1.58200
y_min=-13.61846
z_min=-40.19225
x_max=18.73120
y_max=5.03914
z_max=-24.75795

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

python diffdock/diffdock_using_api.py --input_dir "$sdf_files" --output_dir "$output_path" --receptor_path "$protein_path"

# With box filter - uncomment and modify coordinates as needed
echo "Box coordinates: ${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"
python diffdock/after_diffdock_formatting.py --input_dir "$sdf_files" --output_dir "$output_path" --smiles_file "$smiles_file" --box_filter="${x_min},${y_min},${z_min},${x_max},${y_max},${z_max}"

# Without box filter (original behavior)
# python diffdock/after_diffdock_formatting.py --input_dir $sdf_files --output_dir $output_path --smiles_file $smiles_file

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir $gnina_output_path --csv_path $gnina_csv_path

# Initialize conda and activate environment
eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

