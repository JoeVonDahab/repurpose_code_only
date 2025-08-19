#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scan *.dlg files, pick the most negative
“Estimated Free Energy of Binding = …” value in each one.
Then, add SMILES strings from a specified file.
it outputs a CSV file with the results,
and copies the top 100 ligands with existing PDBQT files
"""

import re
from pathlib import Path
import csv, shutil, sys

INPUT_DLG_DIR   = Path("docking_results_compiled")
INPUT_PDBQT_DIR = Path("docked_converted")
TOP_OUT_DIR     = Path("top_100_ligands")
CSV_OUT         = Path("docking_results_ranked.csv")
SMILES_INPUT_FILE = Path("smiles.smi") # Added SMILES input file
TOP_N           = 100

line_rx = re.compile(
    r"Estimated\s+Free\s+Energy\s+of\s+Binding\s*=\s*([+-]?\d+\.\d+)",
    re.I,
)

def best_energy(dlg: Path):
    best = None
    with dlg.open("r", errors="ignore") as fh:
        for ln in fh:
            m = line_rx.search(ln)
            if m:
                e = float(m.group(1))
                best = e if best is None else min(best, e)
    return best

def pdbqt_name(stem):
    return f"{stem}_best_pose.pdbqt"


def parse_smile_strings(filename: Path):
    """
    Parse SMILES strings from a file.
    Each line should contain a SMILES string and a name, separated by whitespace.
    Returns a dictionary: {name: smiles}
    """
    smiles_dict = {}
    if not filename.exists():
        print(f"Warning: SMILES file not found: {filename}")
        return smiles_dict
    with filename.open("r") as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                smiles = parts[0]
                name = parts[1]
                smiles_dict[name] = smiles
    return smiles_dict

def create_smiles_for_docking_results(smiles_file: Path, docking_data_list: list):
    """
    Add SMILES strings to a list of docking result dictionaries.
    The dictionaries in docking_data_list should have a 'Compound' key.
    A 'SMILES' key will be added or updated in each dictionary.
    """
    smiles_dict = parse_smile_strings(smiles_file)
    for row_dict in docking_data_list:
        name = row_dict['Compound']
        # Assuming compound name in DLG/PDBQT (dlg.stem) might be base name or need parsing.
        base_name = name.split('_')[0] # Use base name for SMILES lookup
        smile = smiles_dict.get(base_name, "N/A")
        row_dict['SMILES'] = smile
    return docking_data_list


def main():
    if not INPUT_DLG_DIR.is_dir():
        sys.exit(f"DLG directory not found: {INPUT_DLG_DIR.resolve()}")
    if not INPUT_PDBQT_DIR.is_dir():
        sys.exit(f"PDBQT directory not found: {INPUT_PDBQT_DIR.resolve()}")
    TOP_OUT_DIR.mkdir(exist_ok=True)

    raw_rows = []
    for dlg in INPUT_DLG_DIR.glob("*.dlg"):
        e = best_energy(dlg)
        
        # Check for corresponding PDBQT file existence
        pdbqt_file_for_dlg = INPUT_PDBQT_DIR / pdbqt_name(dlg.stem)
        pdbqt_exists = pdbqt_file_for_dlg.exists()

        if e is not None:
            raw_rows.append({
                'Compound': dlg.stem,
                'Affinity_kcal_per_mol': e,
                'DLG_file': dlg.name,
                'Within_Pocket': pdbqt_exists # Store PDBQT existence status
            })
        else:
            print(f"No energy line in {dlg.name}")

    if not raw_rows:
        sys.exit("No energies extracted.")

    # Add SMILES strings
    if SMILES_INPUT_FILE.exists():
        print(f"Adding SMILES from {SMILES_INPUT_FILE}")
        # create_smiles_for_docking_results modifies raw_rows in-place
        # and returns the modified list.
        processed_rows = create_smiles_for_docking_results(SMILES_INPUT_FILE, raw_rows)
    else:
        print(f"SMILES file {SMILES_INPUT_FILE} not found. Skipping SMILES addition.")
        # Add 'SMILES' column with N/A if file not found
        for row in raw_rows:
            row['SMILES'] = "N/A"
        processed_rows = raw_rows # raw_rows has been modified with N/A SMILES

    # Sort all processed rows by affinity (lowest first)
    processed_rows.sort(key=lambda x: x['Affinity_kcal_per_mol'])

    with CSV_OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        # Updated CSV header
        w.writerow(["Rank", "Compound", "SMILES", "Affinity_kcal_per_mol", "DLG_file", "Within_Pocket"])
        for rank, data_dict in enumerate(processed_rows, 1):
            w.writerow([
                rank,
                data_dict['Compound'],
                data_dict.get('SMILES', "N/A"), # Ensure SMILES key exists
                f"{data_dict['Affinity_kcal_per_mol']:.3f}",
                data_dict['DLG_file'],
                data_dict['Within_Pocket'] # Write True/False for PDBQT existence
            ])
    print("CSV written:", CSV_OUT)

    copied = 0
    # No filtering by pocket location - process all top ligands regardless of Within_Pocket status
    # The list 'processed_rows' is already sorted by affinity.
    # Select top N ligands for copying (but only copy those with existing PDBQT files)
    top_ligands_to_copy = processed_rows[:TOP_N]

    for data_dict in top_ligands_to_copy:
        name = data_dict['Compound']
        src = INPUT_PDBQT_DIR / pdbqt_name(name)
        
        # Only copy if PDBQT file exists to avoid errors
        if src.exists():
            shutil.copy2(src, TOP_OUT_DIR / src.name)
            copied += 1
        else:
            print(f"Warning: PDBQT file not found for {name}, skipping copy")
        
    print(f"Copied {copied} PDBQT files from top {TOP_N} ligands to {TOP_OUT_DIR}")

if __name__ == "__main__":
    main()

