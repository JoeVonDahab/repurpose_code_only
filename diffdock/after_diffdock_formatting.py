import os
import glob
import json
import csv
import re
import argparse
import numpy as np
from rdkit import Chem

def extract_molecule_name_from_sdf(sdf_content):
    """Extract molecule name from SDF content"""
    lines = sdf_content.strip().split('\n')
    if len(lines) > 0:
        # First line usually contains the molecule name
        name = lines[0].strip()
        if name and name != "":
            return name
    return "Unknown"

def extract_smiles_from_sdf(sdf_content):
    """Extract SMILES from SDF content using RDKit"""
    try:
        mol = Chem.MolFromMolBlock(sdf_content)
        if mol:
            return Chem.MolToSmiles(mol)
    except:
        pass
    return "N/A"

def get_geometric_center_from_sdf(sdf_content):
    """Calculate geometric center from SDF content using coordinates."""
    try:
        mol = Chem.MolFromMolBlock(sdf_content)
        if mol is None:
            return None
        
        conf = mol.GetConformer()
        coords = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        
        if not coords:
            return None
        return np.mean(np.array(coords), axis=0)
    except Exception as e:
        print(f"Warning: Could not calculate geometric center from SDF: {e}")
        return None

def is_center_in_box(center_coords, box_min, box_max):
    """Check if the geometric center is within the defined box."""
    if center_coords is None:
        return False
    return (box_min[0] <= center_coords[0] <= box_max[0] and
            box_min[1] <= center_coords[1] <= box_max[1] and
            box_min[2] <= center_coords[2] <= box_max[2])

def parse_box_coordinates(box_string):
    """Parse box coordinates from string format: 'x_min,y_min,z_min,x_max,y_max,z_max'"""
    try:
        coords = [float(x.strip()) for x in box_string.split(',')]
        if len(coords) != 6:
            raise ValueError("Box coordinates must have exactly 6 values")
        box_min = np.array([coords[0], coords[1], coords[2]])
        box_max = np.array([coords[3], coords[4], coords[5]])
        return box_min, box_max
    except Exception as e:
        raise ValueError(f"Invalid box coordinates format: {e}. Expected format: 'x_min,y_min,z_min,x_max,y_max,z_max'")

def load_smiles_mapping(smiles_file_path):
    """Load SMILES mapping from approved_drugs.smi file"""
    smiles_mapping = {}
    try:
        with open(smiles_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # Split on first space to separate SMILES from molecule name
                    parts = line.split(' ', 1)
                    if len(parts) == 2:
                        smiles, molecule_name = parts
                        # Store with lowercase name for case-insensitive matching
                        smiles_mapping[molecule_name.lower()] = smiles
        print(f"Loaded {len(smiles_mapping)} SMILES entries from {smiles_file_path}")
    except Exception as e:
        print(f"Error loading SMILES file {smiles_file_path}: {e}")
    return smiles_mapping

def extract_molecule_name_from_response(response_file_path):
    """Extract molecule name from DiffDock API response file."""
    try:
        with open(response_file_path, 'r') as f:
            content = f.read()
            
        # Parse JSON to get trajectory
        data = json.loads(content)
        if 'trajectory' in data and data['trajectory']:
            trajectory = data['trajectory'][0]
            
            # Extract what's between "COMPND " and "\nHETATM"
            start = trajectory.find("COMPND ") + 7  # +7 to skip "COMPND "
            end = trajectory.find("\nHETATM")
            
            if start > 6 and end > start:  # start > 6 means "COMPND " was found
                molecule_name = trajectory[start:end].strip()
                return molecule_name
                
    except Exception as e:
        print(f"Error reading response file {response_file_path}: {e}")
    return None

def create_ligand_name_mapping(base_path, smiles_mapping=None):
    """Create mapping from ligand indices to molecule information from API responses."""
    mapping = {}
    
    # Find all ligand folders
    ligand_folders = [d for d in os.listdir(base_path) if d.startswith("ligand_") and os.path.isdir(os.path.join(base_path, d))]
    
    for folder in ligand_folders:
        try:
            # Extract ligand index from folder name
            ligand_index = int(folder.split("_")[1])
            
            # Path to response file
            response_file = os.path.join(base_path, folder, "response_text.txt")
            
            if os.path.exists(response_file):
                # Extract molecule name from API response
                mol_name = extract_molecule_name_from_response(response_file)
                
                if mol_name:
                    # Look up SMILES from mapping (case-insensitive)
                    smiles = "N/A"
                    if smiles_mapping:
                        smiles = smiles_mapping.get(mol_name.lower(), "N/A")
                    
                    mapping[ligand_index] = {
                        'name': mol_name,
                        'original_file': f'{mol_name}.sdf',  # Inferred from molecule name
                        'smiles': smiles
                    }
                    print(f"Mapped ligand_{ligand_index} -> {mol_name} (SMILES: {'Found' if smiles != 'N/A' else 'Not Found'})")
                else:
                    print(f"Could not extract molecule name from {response_file}")
                    mapping[ligand_index] = {
                        'name': f'ligand_{ligand_index}',
                        'original_file': f'ligand_{ligand_index}.sdf',
                        'smiles': 'N/A'
                    }
            else:
                print(f"Response file not found: {response_file}")
                mapping[ligand_index] = {
                    'name': f'ligand_{ligand_index}',
                    'original_file': f'ligand_{ligand_index}.sdf',
                    'smiles': 'N/A'
                }
                
        except (ValueError, IndexError) as e:
            print(f"Error processing folder {folder}: {e}")
            continue
    
    return mapping

def main():
    parser = argparse.ArgumentParser(description='Process DiffDock results and format outputs')
    parser.add_argument('--input_dir', required=True, help='Input directory containing SDF files')
    parser.add_argument('--output_dir', required=True, help='Output directory for DiffDock results')
    parser.add_argument('--smiles_file', default='approved_drugs.smi', help='SMILES mapping file (default: approved_drugs.smi)')
    parser.add_argument('--box_filter', default=None, 
                       help='Optional box coordinates for filtering poses: "x_min,y_min,z_min,x_max,y_max,z_max"')
    
    args = parser.parse_args()
    
    # Configuration
    input_dir = args.input_dir
    base_path = args.output_dir
    smiles_file_path = args.smiles_file
    
    # Parse box filter if provided
    box_min = None
    box_max = None
    if args.box_filter:
        try:
            box_min, box_max = parse_box_coordinates(args.box_filter)
            print(f"Box filter enabled: min={box_min}, max={box_max}")
        except ValueError as e:
            print(f"Error parsing box coordinates: {e}")
            return
    else:
        print("No box filter specified - all poses will be processed")
    
    best_poses_dir = os.path.join(base_path, "best_poses")
    os.makedirs(best_poses_dir, exist_ok=True)

    # Load SMILES mapping from approved_drugs.smi
    print("Loading SMILES mapping...")
    smiles_mapping = load_smiles_mapping(smiles_file_path)

    # Create mapping between ligand numbers and molecule names
    print("Creating ligand name mapping...")
    ligand_mapping = create_ligand_name_mapping(base_path, smiles_mapping)

    ligand_confidences = []
    molecule_data = []

    ligand_folders = glob.glob(os.path.join(base_path, "ligand_*"))
    ligand_numbers = []

    for folder in ligand_folders:
        folder_name = os.path.basename(folder)
        if folder_name.startswith("ligand_"):
            try:
                num = int(folder_name.split("_")[1])
                ligand_numbers.append(num)
            except (ValueError, IndexError):
                continue

    ligand_numbers.sort()

    for i in ligand_numbers:
        ligand_folder = os.path.join(base_path, f"ligand_{i}")
        
        input_file = os.path.join(ligand_folder, "response_text.txt")
        output_folder = os.path.join(ligand_folder, "diffdock_actual_outcome")

        if not os.path.exists(input_file):
            print(f"Missing file in {ligand_folder}, skipping.")
            continue

        os.makedirs(output_folder, exist_ok=True)

        try:
            with open(input_file, "r") as f:
                content = f.read().strip()
                if not content:
                    print(f"Empty file in {ligand_folder}, skipping.")
                    continue
                data = json.loads(content)
        except json.JSONDecodeError as e:
            print(f"Invalid JSON in {ligand_folder}: {e}, skipping.")
            continue
        except Exception as e:
            print(f"Error reading file in {ligand_folder}: {e}, skipping.")
            continue

        # Get molecule information
        mol_info = ligand_mapping.get(i, {
            'name': f'ligand_{i}',
            'original_file': f'ligand_{i}.sdf',
            'smiles': 'N/A'
        })
        
        mol_name = mol_info['name']
        mol_smiles = mol_info['smiles']

        # Write PDB files
        for j, pose in enumerate(data.get("trajectory", []), start=1):
            with open(os.path.join(output_folder, f"pose_{j}.pdb"), "w") as pdb_file:
                pdb_file.write(pose)

        # Write SDF files
        sdf_entries = data.get("ligand_positions", [])
        for j, sdf in enumerate(sdf_entries, start=1):
            with open(os.path.join(output_folder, f"ligand_pose_{j}.sdf"), "w") as sdf_file:
                sdf_file.write(sdf)

        # Write confidence scores
        confidences = data.get("position_confidence", [])
        with open(os.path.join(output_folder, "pose_confidences.txt"), "w") as out_file:
            out_file.write("Rank \t Pose Confidence\n\n")
            for j, conf in enumerate(confidences, start=1):
                out_file.write(f"{j} \t {conf}\n")

        # Find best pose and save to best_poses directory
        valid_confidences = [c for c in confidences if c is not None]
        if valid_confidences and sdf_entries:
            # If box filter is specified, filter poses first
            if box_min is not None and box_max is not None:
                # Find poses within the box
                poses_in_box = []
                for idx, sdf in enumerate(sdf_entries):
                    if confidences[idx] is not None:
                        center = get_geometric_center_from_sdf(sdf)
                        if is_center_in_box(center, box_min, box_max):
                            poses_in_box.append((idx, confidences[idx]))
                
                # If poses found in box, use the best one from the box
                if poses_in_box:
                    # Sort by confidence and get the best pose in the box
                    poses_in_box.sort(key=lambda x: x[1], reverse=True)
                    best_pose_idx, highest_conf = poses_in_box[0]
                    print(f"Found {len(poses_in_box)} poses in box for {mol_name}, using best (confidence: {highest_conf:.4f})")
                    
                    # Save best pose with molecule name
                    clean_name = re.sub(r'[<>:"/\\|?*]', '_', mol_name)  # Clean filename
                    best_pose_filename = f"{clean_name}_best_pose.sdf"
                    best_pose_path = os.path.join(best_poses_dir, best_pose_filename)
                    
                    with open(best_pose_path, "w") as f:
                        f.write(sdf_entries[best_pose_idx])
                    
                    ligand_confidences.append((f"ligand_{i}", highest_conf))
                    molecule_data.append({
                        'ligand_id': f"ligand_{i}",
                        'molecule_name': mol_name,
                        'smiles': mol_smiles,
                        'confidence_score': highest_conf,
                        'original_file': mol_info['original_file']
                    })
                    
                    print(f"Processed {mol_name}: confidence {highest_conf:.4f}")
                else:
                    # No poses in box, skip this ligand
                    print(f"No poses within box for {mol_name}, skipping ligand")
                    molecule_data.append({
                        'ligand_id': f"ligand_{i}",
                        'molecule_name': mol_name,
                        'smiles': mol_smiles,
                        'confidence_score': 'N/A - No poses in box',
                        'original_file': mol_info['original_file']
                    })
            else:
                # No box filter, use original logic
                highest_conf = max(valid_confidences)
                best_pose_idx = confidences.index(highest_conf)
                
                # Save best pose with molecule name
                clean_name = re.sub(r'[<>:"/\\|?*]', '_', mol_name)  # Clean filename
                best_pose_filename = f"{clean_name}_best_pose.sdf"
                best_pose_path = os.path.join(best_poses_dir, best_pose_filename)
                
                with open(best_pose_path, "w") as f:
                    f.write(sdf_entries[best_pose_idx])
                
                ligand_confidences.append((f"ligand_{i}", highest_conf))
                molecule_data.append({
                    'ligand_id': f"ligand_{i}",
                    'molecule_name': mol_name,
                    'smiles': mol_smiles,
                    'confidence_score': highest_conf,
                    'original_file': mol_info['original_file']
                })
                
                print(f"Processed {mol_name}: confidence {highest_conf:.4f}")
        else:
            print(f"No valid confidence values in {ligand_folder}, skipping.")
            molecule_data.append({
                'ligand_id': f"ligand_{i}",
                'molecule_name': mol_name,
                'smiles': mol_smiles,
                'confidence_score': 'N/A',
                'original_file': mol_info['original_file']
            })

    # Sort molecules by confidence score (highest first)
    molecule_data_sorted = sorted([m for m in molecule_data if m['confidence_score'] != 'N/A'], 
                                 key=lambda x: x['confidence_score'], reverse=True)

    # Add molecules with N/A confidence at the end
    molecule_data_sorted.extend([m for m in molecule_data if m['confidence_score'] == 'N/A'])

    # Create CSV file with all molecule data
    csv_path = os.path.join(base_path, "molecule_results_ranked.csv")
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['rank', 'molecule_name', 'ligand_id', 'confidence_score', 'smiles', 'original_file']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for rank, mol in enumerate(molecule_data_sorted, 1):
            writer.writerow({
                'rank': rank,
                'molecule_name': mol['molecule_name'],
                'ligand_id': mol['ligand_id'],
                'confidence_score': mol['confidence_score'],
                'smiles': mol['smiles'],
                'original_file': mol['original_file']
            })

    # Sort ligands by highest confidence descending (for backward compatibility)
    ligand_confidences.sort(key=lambda x: x[1], reverse=True)

    # Select top 100 ligands
    top_100_ligands = ligand_confidences[:100]

    if top_100_ligands:
        # Compute average confidence score
        average_confidence = sum(score for _, score in top_100_ligands) / len(top_100_ligands)

        # Output top ligands (backward compatibility)
        output_path = os.path.join(base_path, "top_100_ligands_by_confidence.txt")
        with open(output_path, "w") as out:
            out.write("Top 100 ligands by highest confidence score:\n\n")
            out.write("Ligand ID\tConfidence Score\n")
            for ligand_id, score in top_100_ligands:
                out.write(f"{ligand_id}\t{score:.4f}\n")
            out.write(f"\nAverage confidence score: {average_confidence:.4f}\n")

        # Print summary
        print(f"\n=== SUMMARY ===")
        print(f"Total molecules processed: {len(molecule_data)}")
        print(f"Molecules with valid confidence scores: {len([m for m in molecule_data if m['confidence_score'] != 'N/A'])}")
        print(f"Best poses saved to: {best_poses_dir}")
        print(f"CSV results saved to: {csv_path}")
        print(f"Selected top {len(top_100_ligands)} ligands by confidence.")
        print(f"Average confidence score of top 100 ligands: {average_confidence:.4f}")
        
        if molecule_data_sorted:
            best_mol = molecule_data_sorted[0]
            print(f"Best molecule: {best_mol['molecule_name']} (confidence: {best_mol['confidence_score']:.4f})")
    else:
        print("No valid confidence scores found!")

if __name__ == "__main__":
    main()