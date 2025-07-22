# smiles_to_3d_sdf.py
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

def process_molecule_worker(record, output_dir):
    idx, smiles, name = record
    safe_name = "".join(c if c.isalnum() else "_" for c in name)
    output_path = os.path.join(output_dir, f"{safe_name}.sdf")

    # Skip if the file already exists to allow for resumable runs
    if os.path.exists(output_path):
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return f"Mol '{safe_name}' - INVALID SMILES, skipped."

    mol.SetProp("_Name", name)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = idx 
    if AllChem.EmbedMolecule(mol, params) == -1:
        return f"Mol '{safe_name}' - Conformer embedding FAILED, skipped."

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass 

    try:
        with Chem.SDWriter(output_path) as writer:
            writer.write(mol)
        return None
    except Exception as e:
        return f"Mol '{safe_name}' - File write FAILED: {e}"

def prepare_ligands_from_smiles(smiles_file, out_dir, n_workers=None):
    os.makedirs(out_dir, exist_ok=True)
    n_workers = n_workers or os.cpu_count()
    work_items = []
    
    with open(smiles_file, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                work_items.append((i, parts[0], parts[1]))
            elif len(parts) == 1:
                work_items.append((i, parts[0], f"ligand_{i+1}"))
    
    print(f"► Processing {len(work_items)} molecules into '{out_dir}' using {n_workers} workers.")
    os.environ['OMP_NUM_THREADS'] = '1'
    errors = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(process_molecule_worker, item, out_dir) for item in work_items]
        for future in tqdm(as_completed(futures), total=len(work_items), desc="Creating 3D Ligands"):
            result = future.result()
            if result is not None:
                errors.append(result)

    print("\n✔︎ Preprocessing complete.")
    if errors:
        print(f"\nEncountered {len(errors)} errors. See errors.log for details.")
        with open('errors.log', 'w') as f:
            for error in errors: f.write(f"{error}\n")

if __name__ == "__main__":
    SMILES_INPUT_FILE = "PPARG.smi"
    OUTPUT_DIRECTORY = "ligand_sdf_files_PPARG"
    # Safe number of workers to prevent system crashes. Tune if needed.
    NUMBER_OF_WORKERS = 24
    
    prepare_ligands_from_smiles(SMILES_INPUT_FILE, OUTPUT_DIRECTORY, NUMBER_OF_WORKERS)