#!/usr/bin/env bash
# =============================================================================
# Create symbolic links from the generic sample IDs used by this repository
# to your real raw-data folders/files, so no real sample names appear in code.
#
#   cp config/raw_data_links.example.tsv config/raw_data_links.tsv   # git-ignored
#   # edit: one line per link  ->  <link path>  <TAB>  <real path>
#   bash scripts/link_raw_data.sh
#
# Link paths are relative to the repository root (or absolute).
# =============================================================================
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MAP="${1:-$REPO_DIR/config/raw_data_links.tsv}"
[[ -f "$MAP" ]] || { echo "ERROR: $MAP not found (copy config/raw_data_links.example.tsv)" >&2; exit 1; }

cd "$REPO_DIR"
n=0
while IFS=$'\t' read -r link target; do
  [[ -z "${link// }" || "$link" == \#* ]] && continue
  [[ -e "$target" ]] || echo "WARNING: target does not exist: $target" >&2
  mkdir -p "$(dirname "$link")"
  ln -sfn "$target" "$link"
  echo "  $link -> $target"
  n=$((n + 1))
done < "$MAP"
echo "Created $n links."
