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

# Use filter_valid_sdfs.py to get valid SDF files and save to temporary file
temp_file=$(mktemp)
python DiffDock-NMDN/filter_valid_sdfs.py --input_dir "$best_poses_path" 2>/dev/null > "$temp_file"

# Check if we have any valid files
if [ ! -s "$temp_file" ]; then
    echo "Error: No valid SDF files found in '$best_poses_path'"
    rm "$temp_file"
    exit 1
fi

# Convert space-separated list to newline-separated list
temp_file_lines=$(mktemp)
tr ' ' '\n' < "$temp_file" > "$temp_file_lines"

# Count files for progress
file_count=$(wc -w < "$temp_file")  # Count words (files) instead of lines
echo "Found $file_count valid SDF files"
echo "Processing $file_count SDF files in batches..."

# Remove existing output file to start fresh
[ -f "$csv_output" ] && rm "$csv_output"

# Process files in batches of 1000 to avoid argument list length issues
batch_size=1000
batch_num=1
temp_csv_dir=$(mktemp -d)

# Split files into batches
line_count=0
while IFS= read -r line || [ -n "$line" ]; do
    if [ -n "$line" ]; then  # Skip empty lines
        echo "$line" >> "${temp_csv_dir}/batch_${batch_num}.txt"
        ((line_count++))
        if [ $line_count -eq $batch_size ]; then
            ((batch_num++))
            line_count=0
        fi
    fi
done < "$temp_file_lines"

# Calculate total batches (adjust if last batch is partial)
if [ $line_count -gt 0 ]; then
    total_batches=$batch_num
else
    total_batches=$((batch_num - 1))
fi

echo "Created $total_batches batches of up to $batch_size files each"

for ((i=1; i<=total_batches; i++)); do
    batch_file="${temp_csv_dir}/batch_${i}.txt"
    if [ -f "$batch_file" ]; then
        batch_file_count=$(wc -l < "$batch_file")
        echo "Processing batch $i/$total_batches ($batch_file_count files)..."
        batch_csv="${temp_csv_dir}/batch_${i}_results.csv"
        
        # Read batch files and run predict.py
        batch_files=$(tr '\n' ' ' < "$batch_file")
        eval "python DiffDock-NMDN/predict.py --prot \"$protein_file\" --save_csv \"$batch_csv\" --ligs $batch_files"
        
        # Check if batch processing was successful
        if [ ! -f "$batch_csv" ]; then
            echo "Error: Batch $i failed to produce results"
            continue
        fi
        
        # Append to main CSV (skip header for subsequent batches)
        if [ $i -eq 1 ]; then
            cp "$batch_csv" "$csv_output"
        else
            tail -n +2 "$batch_csv" >> "$csv_output"
        fi
    else
        echo "Warning: Batch file $batch_file not found"
    fi
done

echo "All batches processed. Final results in: $csv_output"

# Clean up temporary files
rm "$temp_file"
rm "$temp_file_lines"
rm -rf "$temp_csv_dir"
