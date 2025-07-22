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
                if mol is not None and mol.GetNumConformers() > 0:
                    conf = mol.GetConformer()
                    if conf is not None:
                        return True
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
