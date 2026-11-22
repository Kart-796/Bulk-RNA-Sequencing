#!/usr/bin/env bash
# =============================================================================
# Build the STAR genome index from the reference genome + GENCODE GTF.
# --sjdbOverhang should be (read length - 1); set star: sjdb_overhang in the config.
#
# Usage:  bash alignment/01_build_star_index.sh
# =============================================================================
set -euo pipefail

# ---- Configuration: edit config/config.yaml (run from the repository root) ----
REPO_DIR="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cfg() { python3 "$REPO_DIR/utils/get_config.py" "$@"; }
REF_DIR="$(cfg path alignment_dir)/references"
GENCODE_RELEASE="$(cfg value star.gencode_release)"
THREADS="$(cfg value star.threads)"
OVERHANG="$(cfg value star.sjdb_overhang)"
INDEX="$REF_DIR/star_index"
mkdir -p "$INDEX"

STAR --runMode genomeGenerate \
  --genomeDir "$INDEX" \
  --genomeFastaFiles "$REF_DIR/Homo_sapiens_assembly38.fasta" \
  --sjdbGTFfile "$REF_DIR/gencode.v${GENCODE_RELEASE}.annotation.gtf" \
  --sjdbOverhang "$OVERHANG" \
  --runThreadN "$THREADS"

echo "STAR index written to $INDEX"
