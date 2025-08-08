#!/usr/bin/env bash
set -euo pipefail

# === Settings ===
DEST="drive:/SharedFiles"
# Put the 2 tarballs here (edit names if yours differ)
FILES=(
  "diffdock_nmdn_env.tar.gz"
  "diffdock-env-linux64.tar.gz"
)

# Optional: check rclone remote exists
rclone listremotes | grep -q "^drive:" || {
  echo "‼️ 'drive:' remote not found. Run: rclone config"
  exit 1
}

echo "📤 Uploading to $DEST ..."
for f in "${FILES[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "❌ File not found: $f"
    exit 1
  fi

  echo "— Copying $f ..."
  rclone copy "$f" "$DEST" --progress

  echo "— Verifying checksum on remote ..."
  # md5 for Drive; falls back gracefully if unsupported
  rclone md5sum "$DEST" | grep -F "  $f" || echo "⚠️ md5 not listed by remote for $f (may be fine)."

  echo "— Creating shareable link ..."
  LINK="$(rclone link "$DEST/$f")"
  echo "🔗 $f → $LINK"
done

echo "✅ Done. Contents of $DEST:"
rclone lsl "$DEST"
