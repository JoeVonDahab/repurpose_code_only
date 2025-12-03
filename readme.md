# 💊 Multi-Modal Drug Repurposing for Parkinson's Disease

A fully reproducible computational pipeline developed at the Zhao Lab, UCSF. This repository integrates four distinct data modalities—Molecular Docking, Transcriptomics, Real-World Evidence (RWE), and Knowledge Graphs—to identify and validate drug candidates for Parkinson's Disease (PD).

## 🎯 Project Objective

To provide a robust framework for drug repurposing by converging evidence from multiple sources:

- **Virtual Screening**: Benchmarking classical vs. ML docking methods on PD targets (DYRK1A, GCase, LRRK2, USP30).
- **Transcriptomics**: Identifying drugs that reverse PD-specific gene expression signatures (scRNA-seq).
- **Real-World Evidence (RWE)**: Analyzing EHR data (25,548 patients) for survival signals.
- **Pharmacology Graph**: Predicting novel links using a Transformer-based Knowledge Graph.

## 📋 Repository Structure

| Folder / File | Description |
|---------------|-------------|
| `AutoDOCK/` | Scripts and logs for classical AutoDock Vina pipelines. |
| `DiffDock-NMDN/` | Environment for Diffusion Model docking + Neural Mixture Density Network rescoring. |
| `benchmarking/` | Scripts for the LIT-PCBA dataset benchmark (15 protein targets). |
| `diffdock/` | Core execution scripts and inputs for DiffDock runs. |
| `gnina/` | Scripts for GNINA, used for rescoring AutoDock poses (Best performing method). |
| `autodock_run_multiple.sh` | Batch execution script for AutoDock jobs. |
| `gnina_only.sh` | Script to run standalone GNINA rescoring on existing poses. |
| `run_pipeline.sh` | Master script for launching the complete screening workflow. |
| `setup.sh` | Environment and dependency setup script (resets without large files). |

## 🌐 Related Repositories & Data

**Pharmacology Graph Transformer**: [JoeVonDahab/pharmacology-graph](https://github.com/JoeVonDahab/pharmacology-graph)

## 🧪 Docking Strategy & Insights

Our benchmarking on the LIT-PCBA dataset revealed critical insights for pipeline design:

- **Classical + ML Rescoring is Superior**: The most reliable method was GNINA rescoring of AutoDock poses (median EF1% = 2.14).
- **ML Docking Limitations**: Pure ML methods (like DiffDock) were highly target-dependent and failed completely on some targets (e.g., OPRK1, GBA).
- **Consensus Scoring**: Averaging multiple scores (Global Consensus) yielded no meaningful improvement over the best single method.

## 🔧 Installation & Setup

### 1. Dependencies

This pipeline requires the following external tools and repositories:

#### Public Dependencies

- **DiffDock-NMDN**: [SongXia-NYU/DiffDock-NMDN](https://github.com/SongXia-NYU/DiffDock-NMDN)
- **AutoDock-GPU**: [ccsb-scripps/AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU)
- **GNINA**: [Docker Image](https://hub.docker.com/r/gnina/gnina)
- **DiffDock**: [gcorso/DiffDock](https://github.com/gcorso/DiffDock)

#### Private Dependencies (Request Access Required)

- **Single Cell Seq**: [JoeVonDahab/single_cell_seq](https://github.com/JoeVonDahab/single_cell_seq) *(private)*
- **Atorvastatin Project**: [JoeVonDahab/atovastatin_project](https://github.com/JoeVonDahab/atovastatin_project) *(private)*

> **Note**: For access to private repositories, please email the author.

### 2. Environment Setup

```bash
# Clone the repository
git clone https://github.com/JoeVonDahab/repurpose_code_only.git
cd repurpose_code_only

# Run the setup script to install dependencies
./setup.sh
```

### 3. Data Download (Large Files via rclone & gdown)

Large model weights (e.g., for DiffDock) are hosted externally. Use the following commands to manage them.

#### ⬇️ Retrieve Large File (Download)

To download the required model weights (diffdock.tar.gz) from Google Drive:

```bash
set -e

echo "📦 Installing gdown..."
pip install -q gdown

echo "🔽 Downloading diffdock.tar.gz from Google Drive…"
gdown --id 14_Lce88Vb1hL4vuL4KHYlnNlcYA9hzf0 -O diffdock.tar.gz
```

#### ⬆️ Share Large File (Upload)

If you need to upload new weights or datasets to the shared drive:

```bash
# 1. Install rclone
sudo apt update && sudo apt install unzip -y
curl https://rclone.org/install.sh | sudo bash
rclone config

# 2. Copy file to Drive and generate link
rclone copy diffdock.tar.gz drive:/SharedFiles/ --progress
rclone link drive:/SharedFiles/diffdock.tar.gz
```

## 🚀 Quick Start

To run the full virtual screening pipeline (Ligand Prep -> AutoDock -> GNINA Rescoring):

```bash
# Activate the environment
conda activate pd_repurposing

# Execute the master pipeline script
./run_pipeline.sh
```

## ✍️ Contact

**Youssef Abo-Dahab, Pharm.D.** — Student Intern, Zhao Lab, UCSF

GitHub: [@JoeVonDahab](https://github.com/JoeVonDahab)
