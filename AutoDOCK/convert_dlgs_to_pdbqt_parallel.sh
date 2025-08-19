#!/usr/bin/env bash
###############################################################################
#  DLG → PDBQT batch converter, parallel edition
#  - Uses GNU parallel if available (best)
#  - Falls back to xargs -P if not
#  - THREADS defaults to 32 but honours $THREADS env variable
###############################################################################

# === Configuration (unchanged) =============================================
PYTHON_EXEC="/home/joe/miniconda3/envs/autodock/bin/python"
MGLTOOLS_ROOT="/home/joe/miniconda3/envs/autodock/MGLToolsPckgs"
CONVERT_SCRIPT="${MGLTOOLS_ROOT}/AutoDockTools/Utilities24/write_lowest_energy_ligand.py"
SOURCE_DIR="docking_results_compiled"  # <<< where your .dlg files are
DEST_DIR="docked_converted"
THREADS="${THREADS:-32}"               # <<< how many processes you want
# ============================================================================

set -euo pipefail
echo "Starting DLG → PDBQT conversion with $THREADS parallel workers…"

# -- sanity checks (unchanged) ----------------------------------------------
export MGLTOOLS="$MGLTOOLS_ROOT"
export PYTHONPATH="${MGLTOOLS}:${PYTHONPATH:-}"
for p in "$PYTHON_EXEC" "$CONVERT_SCRIPT"; do
  [[ -x $p || -f $p ]] || { echo "ERROR: missing $p"; exit 1; }
done
[[ -d $SOURCE_DIR ]] || { echo "ERROR: $SOURCE_DIR not found"; exit 1; }
mkdir -p "$DEST_DIR"

# -- build file list ---------------------------------------------------------
mapfile -d '' DLG_FILES < <(find "$SOURCE_DIR" -maxdepth 1 -type f -name "*.dlg" -print0)
[[ ${#DLG_FILES[@]} -gt 0 ]] || { echo "No .dlg files found"; exit 0; }
echo "Found ${#DLG_FILES[@]} .dlg files."

# -- helper used by each worker ---------------------------------------------
convert_one() {
  local dlg="$1"
  local base
  base=$(basename "${dlg%.dlg}")
  local out="$DEST_DIR/${base}_best_pose.pdbqt"

  echo "[PID $$] $(basename "$dlg") → $(basename "$out")"
  "$PYTHON_EXEC" "$CONVERT_SCRIPT" -f "$dlg" -o "$out" \
    && [[ -f $out ]] \
    && echo "  ✔︎ done" \
    || echo "  ✖︎ failed"
}
export -f convert_one
export PYTHON_EXEC CONVERT_SCRIPT DEST_DIR  # needed inside subshells

# -- run in parallel ---------------------------------------------------------
if command -v parallel &>/dev/null; then
  printf '%s\0' "${DLG_FILES[@]}" | \
    parallel -0 -P "$THREADS" convert_one {}
else
  # fallback using xargs if GNU parallel isn’t installed
  printf '%s\0' "${DLG_FILES[@]}" | \
    xargs -0 -n1 -P "$THREADS" bash -c 'convert_one "$1"' _
fi

echo "All conversions finished; output in $DEST_DIR"
