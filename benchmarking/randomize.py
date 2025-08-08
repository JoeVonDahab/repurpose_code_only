import random
import os

# this code reads through folders of drugs smiles picks inactives 10% random and all actives and merge them into one

def read_folder(folder_path):
    actives = f'{folder_path}/actives.smi'
    inactives = f'{folder_path}/inactives.smi'

    with open(actives, 'r') as f:
        active_smiles = f.readlines()

    with open(inactives, 'r') as f:
        inactive_smiles = f.readlines()
    ten_percent = max(1, int(0.05 * len(inactive_smiles)))
    inactive_smiles = inactive_smiles[:ten_percent]

    merged = active_smiles + inactive_smiles
    random.shuffle(merged)
    return merged

def read_all_folders():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    folders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]
    
    for folder in folders:
            folder_path = os.path.join(base_dir, folder)
            merged_smiles = read_folder(folder_path)
            with open(f'{folder}/merged.smi', 'w') as f:
                f.writelines(merged_smiles)

if __name__ == "__main__":
    read_all_folders()
    
