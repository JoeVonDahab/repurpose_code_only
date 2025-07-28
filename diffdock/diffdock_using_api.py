import os
import requests
import time
import random
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run DiffDock API for multiple ligands')
    parser.add_argument('--input_dir', required=True, help='Directory containing SDF files')
    parser.add_argument('--output_dir', required=True, help='Output directory for results')
    parser.add_argument('--receptor_path', required=True, help='Path to receptor PDB file')
    parser.add_argument('--force_reprocess', action='store_true', 
                       help='Force reprocessing of all ligands, even if already completed')
    return parser.parse_args()

# ---- CONFIG ----
args = parse_arguments()
input_dir = args.input_dir
output_dir = args.output_dir
receptor_path = args.receptor_path

url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
header_auth = "Bearer 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"

# ---- ASSET UPLOAD FUNCTION ----
def _upload_asset(input_data):
    assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"
    headers = {
        "Authorization": header_auth,
        "Content-Type": "application/json",
        "accept": "application/json",
    }
    s3_headers = {
        "x-amz-meta-nvcf-asset-description": "diffdock-file",
        "content-type": "text/plain",
    }
    payload = {
        "contentType": "text/plain",
        "description": "diffdock-file"
    }

    for attempt in range(5):  # retry up to 5 times
        try:
            response = requests.post(assets_url, headers=headers, json=payload, timeout=30)
            response.raise_for_status()
            asset_url = response.json()["uploadUrl"]
            asset_id = response.json()["assetId"]

            response = requests.put(asset_url, data=input_data, headers=s3_headers, timeout=300)
            response.raise_for_status()

            return asset_id
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                wait = 2 ** attempt + random.uniform(0, 1)
                print(f"[WARN] Rate limited. Retrying after {wait:.2f}s...")
                time.sleep(wait)
            else:
                raise
    raise RuntimeError("Failed to upload asset after multiple attempts")

# ---- UPLOAD PROTEIN ONCE ----
with open(receptor_path, "rb") as f:
    protein_id = _upload_asset(f.read())
print(f"Protein uploaded: {protein_id}")

# ---- PROCESS ONE LIGAND ----
def process_ligand(idx, sdf_file):
    try:
        ligand_path = os.path.join(input_dir, sdf_file)
        out_folder = os.path.join(output_dir, f"ligand_{idx}")
        
        # Check if ligand has already been processed successfully (unless force_reprocess is True)
        if not args.force_reprocess:
            response_file = os.path.join(out_folder, "response_text.txt")
            if os.path.exists(response_file):
                # Check if the response file contains valid data (not empty and has expected content)
                try:
                    with open(response_file, 'r') as f:
                        content = f.read().strip()
                        if content and '"trajectory"' in content:  # Basic check for valid response
                            print(f"Skipping ligand_{idx} ({sdf_file}) - already processed")
                            return
                except Exception:
                    pass  # If we can't read the file, proceed with processing
        
        os.makedirs(out_folder, exist_ok=True)

        with open(ligand_path, "rb") as f:
            ligand_id = _upload_asset(f.read())

        print(f"Ligand {sdf_file} uploaded: {ligand_id}")

        headers = {
            "Content-Type": "application/json",
            "NVCF-INPUT-ASSET-REFERENCES": f"{protein_id},{ligand_id}",
            "Authorization": header_auth
        }

        payload = {
            "ligand": ligand_id,
            "ligand_file_type": "sdf",
            "protein": protein_id,
            "num_poses": 20,
            "time_divisions": 20,
            "steps": 18,
            "save_trajectory": True,
            "is_staged": True
        }

        for attempt in range(5):  # Retry logic for rate-limited inference
            response = requests.post(url, headers=headers, json=payload)
            if response.status_code != 429:
                break
            wait = 2 ** attempt + random.uniform(0, 1)
            print(f"[WARN] Inference rate-limited. Retrying after {wait:.2f}s...")
            time.sleep(wait)

        with open(os.path.join(out_folder, "response_status.txt"), "w") as f:
            f.write(str(response))

        with open(os.path.join(out_folder, "request_url.txt"), "w") as f:
            f.write(url)

        with open(os.path.join(out_folder, "response_text.txt"), "w") as f:
            f.write(response.text)

        print(f"Completed ligand_{idx}: {response.status_code}")
    except Exception as e:
        print(f"[ERROR] Ligand {idx} failed: {e}")

# ---- MULTITHREADING EXECUTION ----
os.makedirs(output_dir, exist_ok=True)
sdf_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".sdf")])

print(f"Found {len(sdf_files)} SDF files to process")

# Check how many are already processed (unless force_reprocess is True)
if not args.force_reprocess:
    already_processed = 0
    for idx, sdf_file in enumerate(sdf_files):
        out_folder = os.path.join(output_dir, f"ligand_{idx}")
        response_file = os.path.join(out_folder, "response_text.txt")
        if os.path.exists(response_file):
            try:
                with open(response_file, 'r') as f:
                    content = f.read().strip()
                    if content and '"trajectory"' in content:
                        already_processed += 1
            except Exception:
                pass

    print(f"Already processed: {already_processed}, To process: {len(sdf_files) - already_processed}")
else:
    print("Force reprocessing enabled - will process all ligands")

# Reduce concurrency to avoid 429 errors
max_workers = min(3, len(sdf_files))  # Try 2–3 threads instead of 8

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(process_ligand, idx, sdf_file) for idx, sdf_file in enumerate(sdf_files)]
    for future in as_completed(futures):
        pass  # Ensures we wait for all tasks
