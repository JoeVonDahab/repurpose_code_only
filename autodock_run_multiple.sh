set -euo pipefail

# VDR

job_source_dict="jobs/benchmarking/VDR/"
protein_path="$job_source_dict/autodock_meta/3a2i_VDR.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/3a2i_ligand.mol2"
output_path="$job_source_dict/3a2i_VDR_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"


# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# Second Protein TP53


job_source_dict="jobs/benchmarking/TP53/"
protein_path="$job_source_dict/autodock_meta/5o1i_TP53.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5o1i_ligand.mol2"
output_path="$job_source_dict/5o1i_TP53_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# Third Protein PPARG


job_source_dict="jobs/benchmarking/PPARG/"
protein_path="$job_source_dict/autodock_meta/PPARG_5y2t.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5y2t_ligand.mol2"
output_path="$job_source_dict/5y2t_PPARG_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# Forth Protein PKM2


job_source_dict="jobs/benchmarking/PKM2/"
protein_path="$job_source_dict/autodock_meta/3gr4_PKM2.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/3gr4_ligand.mol2"
output_path="$job_source_dict/3gr4_PKM2_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# Fifth Protein MTORC1


job_source_dict="jobs/benchmarking/MTORC1/"
protein_path="$job_source_dict/autodock_meta/4jt5_MTORC1.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/4jt5_ligand.mol2"
output_path="$job_source_dict/4jt5_MTORC1_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# Sixth Protein MAPK1


job_source_dict="jobs/benchmarking/MAPK1/"
protein_path="$job_source_dict/autodock_meta/5v62_MAPK1.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5v62_ligand.mol2"
output_path="$job_source_dict/5v62_MAPK1_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# 8 Protein IDH1


job_source_dict="jobs/benchmarking/IDH1/"
protein_path="$job_source_dict/autodock_meta/5tqh_IDH1.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5tqh_ligand.mol2"
output_path="$job_source_dict/5tqh_IDH1_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# Seventh Protein KAT2A


job_source_dict="jobs/benchmarking/KAT2A/"
protein_path="$job_source_dict/autodock_meta/5h84_KAT2A.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5h84_ligand.mol2"
output_path="$job_source_dict/5h84_KAT2A_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# ninth Protein GBA


job_source_dict="jobs/benchmarking/GBA/"
protein_path="$job_source_dict/autodock_meta/2v3d_GBA.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5tqh_ligand.mol2"
output_path="$job_source_dict/2v3d_GBA_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# tenth Protein FEN1


job_source_dict="jobs/benchmarking/FEN1/"
protein_path="$job_source_dict/autodock_meta/5fv7_FEN1.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5fv7_ligand.mol2"
output_path="$job_source_dict/5fv7_FEN1_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# 11th Protein ESR1_ant


job_source_dict="jobs/benchmarking/ESR1_ant/"
protein_path="$job_source_dict/autodock_meta/5ufx_ESR1_ant.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5ufx_ligand.mol2"
output_path="$job_source_dict/5ufx_ESR1_ant_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# 12th Protein ESR1_ago


job_source_dict="jobs/benchmarking/ESR1_ago/"
protein_path="$job_source_dict/autodock_meta/2b1v_ESR1_ago.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/2b1v_ligand.mol2"
output_path="$job_source_dict/2b1v_ESR1_ago_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"


# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# 13th Protein ALDH1


job_source_dict="jobs/benchmarking/ALDH1/"
protein_path="$job_source_dict/autodock_meta/5l2m_ALDH1.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/5l2m_ligand.mol2"
output_path="$job_source_dict/5l2m_ALDH1_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"

# ============================
# ============================
# ============================
# ============================
# ============================
# ============================
# ============================

# 14th Protein ADRB2


job_source_dict="jobs/benchmarking/ADRB2/"
protein_path="$job_source_dict/autodock_meta/4lde_ADRB2.pdb" # your protein
ligand_file="$job_source_dict/autodock_meta/4lde_ligand.mol2"
output_path="$job_source_dict/4lde_ADRB2_autodock_output" # output path for the results
smiles_file="$job_source_dict/autodock_meta/smiles.smi"

#=============================

# keep those as they are preferably 

map_receptor_name="${job_source_dict}autodock_meta/myreceptor_targeted"
gpf_file="${map_receptor_name}.gpf"
glg_file="${map_receptor_name}.glg"
param_file="${job_source_dict}autodock_meta/boron-silicon-atom_par.dat"
MAP_FILE_FLD_FILENAME="${map_receptor_name}.maps.fld"
input_pdbqt_ligands="${job_source_dict}autodock_meta/ligands/"
autodock_dlgs_output_base="${output_path}/autodock_dlgs_output/"
autodock_dlgs_output="${output_path}/autodock_dlgs_output/gpu1_results/"
log_dir="$job_source_dict/autodock_meta/autodock_logs/"
pdbqt_output="${output_path}/autodock_pdbqt_output/"
top_100_ligands_autodock="$output_path/top_100_ligands_autodock/"
csv_autodock_rank="$output_path/csv_autodock_rankings.csv"


#==================================================

best_poses_path=$output_path/best_poses/
gnina_output_path=$output_path/gnina_output/
gnina_csv_path=$output_path/gnina_output.csv
nmdn_csv_path=$output_path/nmdn_output.csv

# gnina container settings, customize:
NUM_GPUS=1
JOBS_PER_GPU=28
CLEANUP_INTERVAL=500
MAX_CONTAINER_TIME=300
GPU_TEMP_CHECK=true
MEMORY_PER_CONTAINER=8G
CPUS_PER_CONTAINER=0.8
NUM_CPU=$(($NUM_GPUS * $JOBS_PER_GPU))
padding=4.0

#======================================

if [ ! -d "$best_poses_path" ]; then
    mkdir -p "$best_poses_path"
fi

if [ ! -d "$gnina_output_path" ]; then
    mkdir -p "$gnina_output_path"
fi

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
fi


source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate meeko_py310

echo "[INFO] Creating grid box..."
bash AutoDOCK/creating_box_with_example_cli.sh "$protein_path" "$ligand_file" "$map_receptor_name" "$padding"

echo "[INFO] Cleaning GPF..."
python AutoDOCK/clean_gfp_cli.py "${gpf_file}"

echo "[INFO] Verifying files..."
[ -f "$gpf_file" ] || { echo "[ERROR] Missing GPF: $gpf_file"; exit 1; }
[ -f "$param_file" ] || { echo "[ERROR] Missing parameter file: $param_file"; exit 1; }

echo "[INFO] Running autogrid4 from GPF directory..."
gpf_dir="$(dirname "$gpf_file")"
(
  cd "$gpf_dir"
  autogrid4 -p "$(basename "$gpf_file")" -l "$(basename "$glg_file")"
)

python AutoDOCK/prepare_ligands_parallel_command_line.py "$smiles_file" "$input_pdbqt_ligands" -n "$NUM_CPU"

eval "$(conda shell.bash hook)"
conda activate autodock

bash AutoDOCK/run_dual_gpu_docking_cli.sh -m "$MAP_FILE_FLD_FILENAME" -l "$input_pdbqt_ligands" -o "$autodock_dlgs_output_base" -L "$log_dir"

echo "Docking is done, now converting dlgs to pdbqt..."

eval "$(conda shell.bash hook)"
conda activate meeko_py310

bash AutoDOCK/convert_dlgs_to_pdbqt_cli.sh -s "$autodock_dlgs_output" -d "$pdbqt_output"

echo "converting is done" 

python AutoDOCK/rank_docking_results_no_filter_cli.py --dlg-dir "$autodock_dlgs_output" --pdbqt-dir "$pdbqt_output" --output-dir "$top_100_ligands_autodock" --csv-output "$csv_autodock_rank" --smiles-file "$smiles_file" 

# converting pdbqt to sdfs
python AutoDOCK/convert_pdbqt_to_sdf_cli.py "$pdbqt_output" "$best_poses_path" -n "$NUM_CPU"

echo "[INFO] Done autoDOCK!!!!!!!!!"

bash gnina/run_gnina_max_throughput_sdf_1.sh --receptor_file "$protein_path" --ligand_files "$best_poses_path" --output_dir "$gnina_output_path" --num_gpus "$NUM_GPUS" --jobs_per_gpu "$JOBS_PER_GPU" --cleanup_interval "$CLEANUP_INTERVAL" --max_container_time "$MAX_CONTAINER_TIME" --gpu_temp_check "$GPU_TEMP_CHECK" --memory_per_container "$MEMORY_PER_CONTAINER" --cpus_per_container "$CPUS_PER_CONTAINER" 2>&1 | tee gnina_output.log

python gnina/make_csv.py --input_dir "$gnina_output_path" --csv_path "$gnina_csv_path"

eval "$(conda shell.bash hook)"
conda activate diffdock_nmdn

bash DiffDock-NMDN/run.sh --best_poses_path "$best_poses_path" --protein "$protein_path" --csv_output "$nmdn_csv_path"

echo "Best poses saved to $best_poses_path"
echo "Gnina run completed. Output saved to $gnina_csv_path"
echo "DiffDock-NMDN run completed. Output saved to $nmdn_csv_path"