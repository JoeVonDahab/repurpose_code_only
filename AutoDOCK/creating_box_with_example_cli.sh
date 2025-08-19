#!/bin/bash

# Usage: ./creating_box_with_example_cli.sh <receptor_pdb> <ligand_mol2> <output_prefix> [padding]
# Example: ./creating_box_with_example_cli.sh 4lde_ADRB2.pdb 4lde_ligand.mol2 myreceptor_targeted 4.0

RECEPTOR_PDB="$1"
LIGAND_MOL2="$2"
OUTPUT_PREFIX="$3"
PADDING="${4:-4.0}" # Default padding is 4.0 if not provided

mk_prepare_receptor.py \
  --read_pdb "$RECEPTOR_PDB" \
  -o "$OUTPUT_PREFIX" -p -g -v \
  --box_enveloping "$LIGAND_MOL2" \
  --padding "$PADDING" -a
