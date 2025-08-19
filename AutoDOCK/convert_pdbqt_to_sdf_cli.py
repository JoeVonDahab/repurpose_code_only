#!/usr/bin/env python3
import os, subprocess, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def _pdbqt_to_sdf(pdbqt_file, out_dir):
    """Convert a single PDBQT file to SDF format"""
    pdbqt_path = Path(pdbqt_file)
    safe_name = pdbqt_path.stem
    sdf_path = Path(out_dir) / f"{safe_name}.sdf"
    
    obabel_cmd = ["obabel", str(pdbqt_path), "-O", str(sdf_path), "-osdf"]
    try:
        result = subprocess.run(obabel_cmd, check=True, capture_output=True, text=True)
        ok = sdf_path.exists()
        return f"[{safe_name}]  {'OK' if ok else 'Failed to create SDF'}"
    except subprocess.CalledProcessError as e:
        return f"[{safe_name}]  Open Babel failed: {e.stderr.strip()}"
    except Exception as e:
        return f"[{safe_name}]  Error: {str(e)}"

def convert_pdbqt_to_sdf(pdbqt_dir, out_dir, n_workers=None):
    """Convert all PDBQT files in a directory to SDF format"""
    pdbqt_dir = Path(pdbqt_dir)
    out_dir = Path(out_dir)
    
    if not pdbqt_dir.is_dir():
        raise ValueError(f"PDBQT directory not found: {pdbqt_dir}")
    
    os.makedirs(out_dir, exist_ok=True)
    n_workers = n_workers or os.cpu_count() or 1
    
    # Find all PDBQT files
    pdbqt_files = list(pdbqt_dir.glob("*.pdbqt"))
    
    if not pdbqt_files:
        print(f"No PDBQT files found in {pdbqt_dir}")
        return
    
    print(f"► {len(pdbqt_files)} PDBQT files → {out_dir}  ({n_workers} parallel workers)\n")
    
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(_pdbqt_to_sdf, pdbqt_file, out_dir) for pdbqt_file in pdbqt_files]
        for future in as_completed(futures):
            print(future.result())
    
    print("\n✔︎ All done.")

def _build_arg_parser():
    p = argparse.ArgumentParser(
        description="Convert PDBQT files to SDF format in parallel using Open Babel.",
        epilog="Converts all .pdbqt files in the input directory to .sdf files in the output directory."
    )
    p.add_argument("pdbqt_dir", help="Directory containing PDBQT files")
    p.add_argument("out_dir", nargs="?", help="Output directory (default: <pdbqt_dir>_sdf)")
    p.add_argument("-n", "--workers", type=int, help="Worker processes (default: CPU count)")
    return p

def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    
    pdbqt_dir = Path(args.pdbqt_dir)
    if not pdbqt_dir.is_dir():
        parser.error(f"PDBQT directory not found: {args.pdbqt_dir}")
    
    out_dir = args.out_dir or f"{pdbqt_dir.name}_sdf"
    convert_pdbqt_to_sdf(args.pdbqt_dir, out_dir, args.workers)

if __name__ == "__main__":
    main()