#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys

def clean_gpf(path: Path):
    drop_full_lines = {
        "map myreceptor_targeted.Si.map",
        "map myreceptor_targeted.CL.map",
        "map myreceptor_targeted.BR.map",
    }

    with path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        stripped = line.strip()

        if stripped in drop_full_lines:
            continue

        if line.startswith("ligand_types"):
            line = line.replace(" CL", "").replace(" BR", "").replace(" Si", "")

        new_lines.append(line)

    with path.open("w", encoding="utf-8") as f:
        f.writelines(new_lines)

def main():
    p = argparse.ArgumentParser(description="Clean a GPF file.")
    p.add_argument("gpf_file", help="Path to the .gpf file to clean")
    args = p.parse_args()

    fp = Path(args.gpf_file)
    if not fp.is_file():
        print(f"Error: {fp} not found", file=sys.stderr)
        sys.exit(1)

    clean_gpf(fp)

if __name__ == "__main__":
    main()
