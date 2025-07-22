#!/usr/bin/env bash

# Initialize variables
input_dir=""
protein_file=""
csv_output=""
best_poses_path=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --protein)
            protein_file="$2"
            shift 2
            ;;
        --csv_output)
            csv_output="$2"
            shift 2
            ;;
        --best_poses_path)
            best_poses_path="$2"
            shift 2
            ;;
        *)
            echo "Unknown option $1"
            echo "Usage: $0 --protein <protein_file> --csv_output <output_csv_file> --best_poses_path <best_poses_directory>"
            echo "Example: $0 --protein receptor_ready_5tbm.pdb --csv_output results_nmdn.csv --best_poses_path /path/to/best_poses"
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$protein_file" ] || [ -z "$csv_output" ] || [ -z "$best_poses_path" ]; then
    echo "Usage: $0 --protein <protein_file> --csv_output <output_csv_file> --best_poses_path <best_poses_directory>"
    echo "Example: $0 --protein receptor_ready_5tbm.pdb --csv_output results_nmdn.csv --best_poses_path /path/to/best_poses"
    exit 1
fi

# Check if protein file exists
if [ ! -f "$protein_file" ]; then
    echo "Error: Protein file '$protein_file' not found"
    exit 1
fi

# Check if best poses directory exists
if [ ! -d "$best_poses_path" ]; then
    echo "Error: Best poses directory '$best_poses_path' not found"
    exit 1
fi

# Use filter_valid_sdfs.py to get valid SDF files and run predict.py only on those
VALID_SDFS=$(python DiffDock-NMDN/filter_valid_sdfs.py --input_dir $best_poses_path)
python DiffDock-NMDN/predict.py --prot "$protein_file" --ligs $VALID_SDFS --save_csv "$csv_output"
