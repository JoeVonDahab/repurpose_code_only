# prepare_ligands_parallel.py
import os, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem

# ---------- per-molecule worker ----------
def _smiles_to_pdbqt(record, out_dir):
    """
    `record` is (idx, smiles, name) so we can keep ordering if we want.
    Runs in a *separate* process → do not rely on globals.
    """
    idx, smiles, name = record
    safe_name = "".join(c if c.isalnum() else "_" for c in name) or f"ligand_{idx}"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"[{safe_name}]  INVALID SMILES → skipped"

    mol = Chem.AddHs(mol)

    # 3-D conformer
    params = AllChem.ETKDGv3()
    params.randomSeed = (idx ^ 0xF00D) & 0xFFFFFFFF  # different but reproducible
    if AllChem.EmbedMolecule(mol, params) == -1:
        return f"[{safe_name}]  Embed failed → skipped"

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        # keep going with unoptimised geometry
        pass

    sdf_path = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf").name
    Chem.SDWriter(sdf_path).write(mol)

    pdbqt_path = os.path.join(out_dir, f"{safe_name}.pdbqt")
    obabel_cmd = [
        "obabel", sdf_path, "-O", pdbqt_path,
        "-opdbqt", "-h", "--partialcharge", "gasteiger"
    ]
    try:
        subprocess.run(obabel_cmd, check=True, capture_output=True)
        ok = os.path.exists(pdbqt_path)
    finally:
        os.remove(sdf_path)

    return f"[{safe_name}]  {'OK' if ok else 'Open Babel failed'}"

# ---------- coordinator ----------
def prepare_ligands_from_smiles(smiles_file, out_dir, n_workers=None):
    os.makedirs(out_dir, exist_ok=True)
    n_workers = n_workers or os.cpu_count() or 1

    # read & build work list -------------------------------------------------
    work = []
    with open(smiles_file) as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            parts = line.split()
            work.append((i + 1, parts[0], parts[1] if len(parts) > 1 else f"ligand_{i+1}"))

    print(f"► {len(work)} molecules → {out_dir}  ({n_workers} parallel workers)\n")

    # optional: avoid RDKit/OpenMP over-subscription
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    # parallel map -----------------------------------------------------------
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(_smiles_to_pdbqt, rec, out_dir) for rec in work]
        for f in as_completed(futures):
            print(f.result())

    print("\n✔︎ All done.")

# ---------- CLI ----------
if __name__ == "__main__":
    prepare_ligands_from_smiles(
        smiles_file="smiles.smi",
        out_dir="ligands_pdbqt",
        n_workers=32,          # ← use all 32 logical cores
    )
