#!/usr/bin/env bash
# =============================================================================
# Download the hg38 genome (UCSC) and GENCODE gene annotation into
# <alignment_dir>/references.
#
# Usage:  bash alignment/00_download_references.sh
# =============================================================================
set -euo pipefail

# ---- Configuration: edit config/config.yaml (run from the repository root) ----
REPO_DIR="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cfg() { python3 "$REPO_DIR/utils/get_config.py" "$@"; }
REF_DIR="$(cfg path alignment_dir)/references"
GENCODE_RELEASE="$(cfg value star.gencode_release)"
mkdir -p "$REF_DIR" && cd "$REF_DIR"

# hg38 reference genome (UCSC)
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -f hg38.fa.gz
mv hg38.fa Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta

# GTF annotation (GENCODE, hosted by EBI)
wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
gunzip -f "gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"

echo "References ready in $REF_DIR"
