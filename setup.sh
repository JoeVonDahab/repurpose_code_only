#!/bin/bash
set -e

echo "📦 Installing gdown..."
pip install -q gdown

echo "🔽 Downloading diffdock.tar.gz from Google Drive…"
gdown --id 14_Lce88Vb1hL4vuL4KHYlnNlcYA9hzf0 -O diffdock.tar.gz

echo "📂 Extracting archive…"
tar -xzf diffdock.tar.gz

echo "📦 Creating conda environment…"
conda env create -f diffdock/environment.yml -n diffdock

echo "✅ Done. Activate with: conda activate diffdock"
