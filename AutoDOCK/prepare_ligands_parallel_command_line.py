#!/usr/bin/env python3
import os, subprocess, tempfile, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem

def _smiles_to_pdbqt(record, out_dir):
    idx, smiles, name = record
    safe_name = "".join(c if c.isalnum() else "_" for c in name) or f"ligand_{idx}"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"[{safe_name}]  INVALID SMILES → skipped"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = (idx ^ 0xF00D) & 0xFFFFFFFF
    if AllChem.EmbedMolecule(mol, params) == -1:
        return f"[{safe_name}]  Embed failed → skipped"
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        pass
    sdf_path = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf").name
    Chem.SDWriter(sdf_path).write(mol)
    pdbqt_path = os.path.join(out_dir, f"{safe_name}.pdbqt")
    obabel_cmd = ["obabel", sdf_path, "-O", pdbqt_path, "-opdbqt", "-h", "--partialcharge", "gasteiger"]
    try:
        subprocess.run(obabel_cmd, check=True, capture_output=True)
        ok = os.path.exists(pdbqt_path)
    finally:
        os.remove(sdf_path)
    return f"[{safe_name}]  {'OK' if ok else 'Open Babel failed'}"

def prepare_ligands_from_smiles(smiles_file, out_dir, n_workers=None):
    os.makedirs(out_dir, exist_ok=True)
    n_workers = n_workers or os.cpu_count() or 1
    work = []
    with open(smiles_file) as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            parts = line.split()
            work.append((i + 1, parts[0], parts[1] if len(parts) > 1 else f"ligand_{i+1}"))
    print(f"► {len(work)} molecules → {out_dir}  ({n_workers} parallel workers)\n")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for f in as_completed([pool.submit(_smiles_to_pdbqt, rec, out_dir) for rec in work]):
            print(f.result())
    print("\n✔︎ All done.")

def _build_arg_parser():
    p = argparse.ArgumentParser(
        description="Generate 3D PDBQT ligand files from a SMILES (.smi) list in parallel.",
        epilog="Input lines: '<SMILES> <optional_name>'."
    )
    p.add_argument("smiles_file", help="Path to SMILES (.smi) file")
    p.add_argument("out_dir", nargs="?", help="Output directory (default: <smiles_file_stem>_pdbqt)")
    p.add_argument("-n", "--workers", type=int, help="Worker processes (default: CPU count)")
    return p

def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    if not os.path.isfile(args.smiles_file):
        parser.error(f"SMILES file not found: {args.smiles_file}")
    out_dir = args.out_dir or f"{Path(args.smiles_file).stem}_pdbqt"
    prepare_ligands_from_smiles(args.smiles_file, out_dir, args.workers)

if __name__ == "__main__":
    main()