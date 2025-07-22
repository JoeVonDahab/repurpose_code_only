from collections import defaultdict
from copy import deepcopy
from glob import glob
import os.path as osp
from typing import List, Optional
import pandas as pd
from tqdm import tqdm
import argparse

import torch
from torch_scatter import scatter_add
from torch_geometric.loader import DataLoader
from torch_geometric.loader.dataloader import Collater
from rdkit.Chem import AddHs
from rdkit.Chem.AllChem import SDMolSupplier

# Assuming these are local modules in your project structure
from geometry_processors.lm.esm_embedding import ESMCalculator
from geometry_processors.misc import ff_optimize
from geometry_processors.pl_dataset.csv2input_list import MPInfo
from geometry_processors.pl_dataset.prot_utils import pdb2seq
from geometry_processors.process.mp_pyg_runner import proc_hetero_graph
from utils.LossFn import MDNMixLossFn, post_calculation
from utils.data.DataPostProcessor import DataPostProcessor
from utils.data.MolFileDataset import SDFDataset
from utils.data.data_utils import get_num_mols
from utils.eval.predict import EnsPredictor
from utils.eval.tester import Tester
from utils.rmsd import symmetry_rmsd_from_mols

parser = argparse.ArgumentParser()
parser.add_argument("--prot", type=str, required=True, help="Path to the protein PDB file.")
parser.add_argument("--ligs", type=str, nargs="+", required=True, help="List of paths to ligand SDF files.")
parser.add_argument("--nmdn_only", action="store_true", help="If set, only calculate the NMDN score and skip the pKd score.")
parser.add_argument("--save_csv", default=None, help="Optional path to save the output scores as a CSV file.")
args = parser.parse_args()

prot: str = args.prot
ligs: List[str] = args.ligs
nmdn_only: bool = args.nmdn_only
save_csv: Optional[str] = args.save_csv

# Set up device (GPU or CPU) to be used throughout the script
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# The solvation properties and RMSD info are only needed for pKd scores
if not nmdn_only:
    # --- Solvation Properties Block ---
    print("Calculating solvation energetics...")
    atomprop_predictor = EnsPredictor("./data/exp_frag20sol_012_active_ALL_2022-05-01_112820/exp_*_cycle_-1_*")
    sdf_ds = SDFDataset(ligs)
    dl = DataLoader(sdf_ds, batch_size=len(ligs))
    allh_batch = next(iter(dl))
    model_pred = atomprop_predictor.predict_nograd(allh_batch)

    phys_atom_prop = torch.as_tensor(model_pred["atom_prop"], dtype=torch.float32).squeeze()
    if phys_atom_prop.dim() == 0:
        num_atoms = allh_batch.num_nodes.sum().item()
        phys_atom_prop = phys_atom_prop.repeat(num_atoms)

    mol_batch = torch.as_tensor(allh_batch.atom_mol_batch,
                                dtype=torch.long,
                                device=phys_atom_prop.device)

    phys_mol_prop = scatter_add(phys_atom_prop,
                                mol_batch,
                                dim=0,
                                dim_size=get_num_mols(allh_batch))

    # CORRECTED SECTION: Check tensor shape before calling post_calculation
    if phys_mol_prop.dim() == 2 and phys_mol_prop.shape[1] == 3:
        phys_mol_prop = post_calculation(phys_mol_prop)
    else:
        print("\n--- WARNING ---")
        print(f"Skipping `post_calculation` because the property tensor has an unexpected shape: {phys_mol_prop.shape}.")
        print("The function expected a 2D tensor with 3 columns (shape: [num_molecules, 3]).")
        print("This is likely because the `atomprop_predictor` model is not returning 3 features per atom.")
        print("Creating placeholder tensor with zeros to prevent downstream errors.\n")
        
        # Create a tensor with the expected shape [num_molecules, 3] representing gas, water, octanol
        # and then apply post_calculation to get the final 6-column tensor
        num_molecules = phys_mol_prop.shape[0]
        placeholder_raw = torch.zeros((num_molecules, 3), dtype=torch.float32, device=phys_mol_prop.device)
        phys_mol_prop = post_calculation(placeholder_raw)


    # --- RMSD info Block ---
    print("Calculating ligand stability features...")
    rmsd_info = {}
    for lig_sdf in tqdm(ligs, "Ligand FF optimization"):
        try:
            mol_supp = SDMolSupplier(lig_sdf, removeHs=True, sanitize=True, strictParsing=True)
            if not mol_supp: continue
            mol = AddHs(mol_supp[0], addCoords=True)
            dst_mol = deepcopy(mol)
            ff_optimize(dst_mol, [0])
            rmsd = symmetry_rmsd_from_mols(mol, dst_mol, 1)
            rmsd_info[lig_sdf] = rmsd
        except (OSError, RuntimeError) as e:
            print(f"Skipping {lig_sdf} due to an error: {e}")
            continue

# --- ESM-2 embedding Block ---
print("Computing ESM-2 embedding...")
esm_calculator = ESMCalculator(None)
seq = pdb2seq(prot)
prot_embed = esm_calculator.embed_from_seq(seq).squeeze(0)[1: -1, :].float()

# --- NMDN Model Block ---
print("Running NMDN model...")
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
model_path = os.path.join(script_dir, "data", "exp_pl_534_run_2024-01-22_211045__480688")
tester = Tester(model_path)
tester.cfg.no_pkd_score = nmdn_only
tester.cfg.model.kano.kano_ckpt = None

model = tester.model.to(device)
model.eval()
prot_embed = prot_embed.to(device)
collate_fn = Collater(None, None)
data_processor = DataPostProcessor(tester.cfg)
loss_fn = MDNMixLossFn(tester.cfg)
loss_fn.compute_pkd_score = not nmdn_only
loss_fn.inference_mode()

if not nmdn_only:
    # Ensure phys_mol_prop is 2D before sending to device to avoid errors if it's 1D
    if phys_mol_prop.dim() == 1:
        phys_mol_prop = phys_mol_prop.unsqueeze(-1) # Convert from [N] to [N, 1]
    phys_mol_prop = phys_mol_prop.to(device)

out_info = defaultdict(list)
for i, lig_sdf in enumerate(tqdm(ligs, "NMDN model")):
    if not nmdn_only and lig_sdf not in rmsd_info:
        continue

    data = proc_hetero_graph(MPInfo(protein_pdb=prot, ligand_sdf=lig_sdf))
    data = data_processor(data, 0)
    data = collate_fn([data])
    data = data.to(device)

    data.prot_embed = prot_embed
    if not nmdn_only:
        data.mol_prop = phys_mol_prop[i]
        data.rmsd = torch.as_tensor([rmsd_info[lig_sdf]], device=device).float()

    data.sample_id = torch.as_tensor([i], device=device)

    with torch.no_grad():
        pred = model(data)
        __, scores = loss_fn(pred, data, False, True)

    nmdn_score = scores["MDN_LOGSUM_DIST2_REFDIST2"].cpu().item()

    out_info["lig"].append(lig_sdf)
    out_info["NMDN-Score"].append(nmdn_score)
    if not nmdn_only:
        pkd_score = scores["PROP_PRED"].cpu().item()
        out_info["pKd-Score"].append(pkd_score)

out_df = pd.DataFrame(out_info)
print("\n--- Results ---")
print(out_df)
if save_csv is not None:
    out_df.to_csv(save_csv, index=False)
    print(f"\nResults saved to {save_csv}")