how to setup envionment 

```bash
conda env create -f environment.yml
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv --no-cache-dir -f https://data.pyg.org/whl/torch-2.2.0+cu121.html
pip install tensorboard

```
* Open run_gnina_max_throughput.sh and edit the part:
```bash

RECEPTOR_FILE="2GQG_one_chain_docking.pdb"
LIGAND_DIR="all_sdf" # adjust also to the recptors list
```
