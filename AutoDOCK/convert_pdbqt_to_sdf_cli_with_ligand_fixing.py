#!/usr/bin/env python3
import os, subprocess, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import SDMolSupplier, SDWriter

def repair_valence_errors(mol):
    """
    Repair valence errors in a molecule while preserving 3D coordinates.
    
    Strategy:
    1. Find over-valent atoms (N with >3 bonds, C with >4 bonds, etc.)
    2. Remove excess bonds starting with the longest/weakest ones
    3. Preserve the most chemically reasonable connections
    4. Keep all coordinates intact
    
    Returns: (repaired_mol, success_flag, repair_log)
    """
    if mol is None:
        return None, False, "Input molecule is None"
    
    # Make a copy to work with
    mol_copy = Chem.RWMol(mol)
    repair_log = []
    
    # Define normal valence limits
    normal_valences = {
        'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'P': 5, 'S': 6, 'Cl': 1, 'Br': 1, 'I': 1
    }
    
    try:
        # Check each atom for valence violations
        atoms_to_fix = []
        for i, atom in enumerate(mol_copy.GetAtoms()):
            symbol = atom.GetSymbol()
            current_valence = atom.GetTotalValence()
            max_valence = normal_valences.get(symbol, 4)  # Default to 4 if unknown
            
            if current_valence > max_valence:
                atoms_to_fix.append((i, symbol, current_valence, max_valence))
                repair_log.append(f"Over-valent {symbol} at position {i}: {current_valence} bonds (max {max_valence})")
        
        if not atoms_to_fix:
            return mol, True, "No valence errors found"
        
        # Fix each over-valent atom
        for atom_idx, symbol, current_val, max_val in atoms_to_fix:
            atom = mol_copy.GetAtomWithIdx(atom_idx)
            bonds_to_remove = current_val - max_val
            
            # Get all bonds for this atom with their lengths (for prioritizing removal)
            atom_bonds = []
            for bond in atom.GetBonds():
                other_idx = bond.GetOtherAtomIdx(atom_idx)
                # Calculate bond length using 3D coordinates
                conf = mol_copy.GetConformer()
                pos1 = conf.GetAtomPosition(atom_idx)
                pos2 = conf.GetAtomPosition(other_idx)
                bond_length = ((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2)**0.5
                
                atom_bonds.append((bond.GetIdx(), bond_length, other_idx, bond.GetBondType()))
            
            # Sort bonds by length (longest first) - remove longest bonds first
            atom_bonds.sort(key=lambda x: x[1], reverse=True)
            
            # Remove the longest bonds
            bonds_removed = 0
            for bond_idx, bond_length, other_idx, bond_type in atom_bonds:
                if bonds_removed >= bonds_to_remove:
                    break
                
                try:
                    other_atom = mol_copy.GetAtomWithIdx(other_idx)
                    repair_log.append(f"Removing bond {atom_idx}({symbol})-{other_idx}({other_atom.GetSymbol()}) "
                                    f"length={bond_length:.2f}Å type={bond_type}")
                    mol_copy.RemoveBond(atom_idx, other_idx)
                    bonds_removed += 1
                except:
                    continue
        
        # Try to sanitize the repaired molecule
        try:
            Chem.SanitizeMol(mol_copy)
            repair_log.append("Successfully sanitized repaired molecule")
            return mol_copy.GetMol(), True, "; ".join(repair_log)
        except Exception as e:
            repair_log.append(f"Sanitization failed after repair: {str(e)}")
            return None, False, "; ".join(repair_log)
            
    except Exception as e:
        repair_log.append(f"Repair process failed: {str(e)}")
        return None, False, "; ".join(repair_log)

def repair_sdf_file(sdf_file):
    """
    Repair valence errors in an SDF file and save the corrected version.
    Returns: (success, error_message, repair_log)
    """
    try:
        # Read the molecule without sanitization
        mol = Chem.SDMolSupplier(str(sdf_file), sanitize=False, removeHs=False)[0]
        if mol is None:
            return False, "Could not read molecule from SDF", ""
        
        # Attempt repair
        repaired_mol, success, repair_log = repair_valence_errors(mol)
        
        if success and repaired_mol is not None:
            # Write the repaired molecule back to the same file
            writer = SDWriter(str(sdf_file))
            writer.write(repaired_mol)
            writer.close()
            return True, "Successfully repaired", repair_log
        else:
            return False, "Repair failed", repair_log
            
    except Exception as e:
        return False, f"Exception during repair: {str(e)}", ""

def is_valid_sdf_strict(sdf_file):
    """Strict SDF validation - same as filter_valid_sdfs.py with detailed error reporting"""
    try:
        with SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False, strictParsing=False) as supp:
            for mol in supp:
                if mol is None:
                    return False, "Failed to read molecule from SDF"
                if mol.GetNumConformers() == 0:
                    return False, "No conformers found"
                conf = mol.GetConformer()
                if conf is None:
                    return False, "Conformer is None"
                    
                # Additional validation for AutoDock converted files
                try:
                    # Try to sanitize the molecule to catch chemical errors
                    Chem.SanitizeMol(mol)
                    
                    # Check for valid coordinates (not all zeros)
                    coords = conf.GetPositions()
                    if len(coords) == 0:
                        return False, "No coordinates found"
                    
                    # Check if all coordinates are zero (invalid)
                    if all(abs(x) < 1e-6 and abs(y) < 1e-6 and abs(z) < 1e-6 for x, y, z in coords):
                        return False, "All coordinates are zero"
                        
                    # Try to add hydrogens to catch valence errors
                    mol_with_h = Chem.AddHs(mol, addCoords=True)
                    if mol_with_h is None:
                        return False, "Failed to add hydrogens"
                        
                    return True, "Valid"
                except Chem.AtomValenceException as e:
                    return False, f"Valence error: {str(e)}"
                except Chem.KekulizeException as e:
                    return False, f"Kekulize error: {str(e)}"
                except ValueError as e:
                    return False, f"Value error: {str(e)}"
                except AttributeError as e:
                    return False, f"Attribute error: {str(e)}"
        return False, "No molecules found in SDF"
    except Exception as e:
        return False, f"Exception during validation: {str(e)}"

def convert_single_pdbqt(pdbqt_file, out_dir, attempt_repair=True):
    """Convert a single PDBQT file to SDF with optional valence repair"""
    pdbqt_path = Path(pdbqt_file)
    
    if not pdbqt_path.exists():
        return False, "PDBQT file not found"
        
    # Keep the full filename including _best_pose
    sdf_path = Path(out_dir) / f"{pdbqt_path.stem}.sdf"
    
    # Try converting this pose
    obabel_cmd = ["obabel", str(pdbqt_file), "-O", str(sdf_path), "-osdf"]
    try:
        result = subprocess.run(obabel_cmd, check=True, capture_output=True, text=True)
        
        if not sdf_path.exists():
            return False, "SDF file not created"
        
        # Check if repair is needed
        is_valid, error_msg = is_valid_sdf_strict(sdf_path)
        
        if is_valid:
            return True, "Converted successfully"
        elif attempt_repair and ("Valence error" in error_msg):
            # Try to repair valence errors
            repair_success, repair_error, repair_log = repair_sdf_file(sdf_path)
            
            if repair_success:
                # Re-validate after repair
                is_valid_after_repair, _ = is_valid_sdf_strict(sdf_path)
                if is_valid_after_repair:
                    return True, f"Repaired and validated: {repair_log}"
                else:
                    sdf_path.unlink()  # Remove failed repair
                    return False, f"Repair failed validation: {repair_log}"
            else:
                sdf_path.unlink()  # Remove failed file
                return False, f"Repair failed: {repair_error}"
        else:
            sdf_path.unlink()  # Remove invalid file
            return False, error_msg
            
    except subprocess.CalledProcessError as e:
        return False, f"Conversion failed: {e}"

def convert_pdbqt_to_sdf(pdbqt_dir, out_dir, n_workers=None, enable_repair=True):
    """Convert all PDBQT files in a directory to SDF format (best poses only)"""
    
    pdbqt_dir = Path(pdbqt_dir)
    out_dir = Path(out_dir)
    
    if not pdbqt_dir.is_dir():
        raise ValueError(f"PDBQT directory not found: {pdbqt_dir}")
    
    os.makedirs(out_dir, exist_ok=True)
    n_workers = n_workers or os.cpu_count() or 1
    
    # Find all _best_pose.pdbqt files
    best_pose_files = list(pdbqt_dir.glob("*_best_pose.pdbqt"))
    
    if not best_pose_files:
        print(f"No *_best_pose.pdbqt files found in {pdbqt_dir}")
        return
    
    print(f"► {len(best_pose_files)} PDBQT files → {out_dir}  ({n_workers} parallel workers)")
    print("Converting PDBQT files to SDF with repair...")
    
    # Convert all best poses
    conversion_details = {}  # Store conversion details
    successful_conversions = 0
    
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(convert_single_pdbqt, pdbqt_file, out_dir, enable_repair) 
                  for pdbqt_file in best_pose_files]
        
        for future, pdbqt_file in zip(as_completed(futures), best_pose_files):
            file_name = pdbqt_file.stem  # e.g., "87349289_best_pose"
            success, message = future.result()
            conversion_details[file_name] = message
            
            if success:
                print(f"[{file_name}] {message}")
                successful_conversions += 1
            else:
                print(f"[{file_name}] Failed: {message}")
    
    # Print final statistics
    total = len(best_pose_files)
    print(f"\n📊 Conversion Results:")
    print(f"   Successful: {successful_conversions:4d} ({successful_conversions/total*100:.1f}%)")
    print(f"   Failed:     {total - successful_conversions:4d} ({(total - successful_conversions)/total*100:.1f}%)")
    print(f"   Total:      {total:4d}")
    
    # Count repairs
    repaired_count = sum(1 for message in conversion_details.values() 
                        if 'Repaired and validated' in message)
    if repaired_count > 0:
        print(f"🔧 Repaired {repaired_count} molecules with valence errors")
    
    # Print detailed failure analysis
    failed_files = {name: message for name, message in conversion_details.items() 
                   if not message.startswith(('Converted successfully', 'Repaired and validated'))}
    
    if failed_files:
        print(f"\n❌ Detailed Failure Analysis ({len(failed_files)} files):")
        print("="*80)
        
        # Group failures by error type
        error_groups = {}
        for file_name, error_msg in failed_files.items():
            if error_msg not in error_groups:
                error_groups[error_msg] = []
            error_groups[error_msg].append(file_name)
        
        # Print summary by error type
        for error_type, files in error_groups.items():
            print(f"\n📋 {error_type} ({len(files)} files):")
            for file_name in sorted(files):
                print(f"   {file_name}")
        
        print("="*80)

def _build_arg_parser():
    p = argparse.ArgumentParser(
        description="Convert PDBQT files to SDF format with multiple pose fallback.",
        epilog="Converts multiple poses from .pdbqt files to .sdf files, trying alternative poses if validation fails."
    )
    p.add_argument("pdbqt_dir", help="Directory containing PDBQT files with multiple poses")
    p.add_argument("out_dir", nargs="?", help="Output directory (default: <pdbqt_dir>_sdf)")
    p.add_argument("-n", "--workers", type=int, help="Worker processes (default: CPU count)")
    p.add_argument("--no-repair", action="store_true", help="Disable automatic valence repair")
    return p

def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    
    pdbqt_dir = Path(args.pdbqt_dir)
    if not pdbqt_dir.is_dir():
        parser.error(f"PDBQT directory not found: {args.pdbqt_dir}")
    
    out_dir = args.out_dir or f"{pdbqt_dir.name}_sdf"
    enable_repair = not args.no_repair
    convert_pdbqt_to_sdf(args.pdbqt_dir, out_dir, args.workers, enable_repair)

if __name__ == "__main__":
    main()