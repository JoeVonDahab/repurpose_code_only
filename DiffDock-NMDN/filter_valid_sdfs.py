import os
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

def is_valid_sdf(sdf_file):
    """Check if an SDF file is valid and has conformers"""
    try:
        with SDMolSupplier(sdf_file, removeHs=False, sanitize=False, strictParsing=False) as supp:
            for mol in supp:
                if mol is None:
                    return False
                if mol.GetNumConformers() == 0:
                    return False
                conf = mol.GetConformer()
                if conf is None:
                    return False
                    
                # Additional validation for AutoDock converted files
                try:
                    # Try to sanitize the molecule to catch chemical errors
                    Chem.SanitizeMol(mol)
                    
                    # Check for valid coordinates (not all zeros)
                    coords = conf.GetPositions()
                    if len(coords) == 0:
                        return False
                    
                    # Check if all coordinates are zero (invalid)
                    if all(abs(x) < 1e-6 and abs(y) < 1e-6 and abs(z) < 1e-6 for x, y, z in coords):
                        return False
                        
                    # Try to add hydrogens to catch valence errors
                    mol_with_h = Chem.AddHs(mol, addCoords=True)
                    if mol_with_h is None:
                        return False
                        
                    return True
                except (Chem.AtomValenceException, Chem.KekulizeException, ValueError, AttributeError):
                    return False
        return False
    except Exception:
        return False

def filter_valid_sdfs(input_dir):
    """Filter out corrupted SDF files and return list of valid ones"""
    sdf_files = []
    for file in os.listdir(input_dir):
        if file.endswith('.sdf'):
            full_path = os.path.join(input_dir, file)
            if is_valid_sdf(full_path):
                sdf_files.append(full_path)
            else:
                print(f"Skipping corrupted file: {file}", file=sys.stderr)
    
    return sdf_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter valid SDF files from a directory")
    parser.add_argument("--input_dir", required=True, help="Directory containing SDF files")
    args = parser.parse_args()
    
    input_dir = os.path.expanduser(args.input_dir)
    valid_files = filter_valid_sdfs(input_dir)
    
    print(f"Found {len(valid_files)} valid SDF files out of {len([f for f in os.listdir(input_dir) if f.endswith('.sdf')])}", file=sys.stderr)
    
    # Print valid files separated by spaces
    print(" ".join(valid_files))