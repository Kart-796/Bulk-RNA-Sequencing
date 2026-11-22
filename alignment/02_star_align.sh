#!/usr/bin/env bash
# =============================================================================
# STAR alignment of paired-end RNA-seq reads, one sample at a time
# (same STAR settings as the first alignment pass of the GATK RNA-seq pipeline).
#
# Input : <fastq_dir>/<sample_id>_R1.fastq.gz and _R2.fastq.gz
#         sample IDs from config/samples_bulk.csv
# Output: <alignment_dir>/aligned/<sample_id>_sorted.bam (+ .bai) and STAR logs
#
# Usage:  bash alignment/02_star_align.sh
# =============================================================================
set -euo pipefail

# ---- Configuration: edit config/config.yaml (run from the repository root) ----
REPO_DIR="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cfg() { python3 "$REPO_DIR/utils/get_config.py" "$@"; }
FASTQ_DIR="$(cfg path fastq_dir)"
BASE="$(cfg path alignment_dir)"
GENCODE_RELEASE="$(cfg value star.gencode_release)"
THREADS="$(cfg value star.threads)"
INDEX="$BASE/references/star_index"
GTF="$BASE/references/gencode.v${GENCODE_RELEASE}.annotation.gtf"
OUT="$BASE/aligned"
mkdir -p "$OUT" "$BASE/logs"
read -r -a SAMPLES <<< "$(cfg samples bulk)"

for SAMPLE in "${SAMPLES[@]}"; do
  if [[ -s "$OUT/${SAMPLE}_sorted.bam" ]]; then
    echo "[$SAMPLE] already aligned, skipping"; continue
  fi
  echo "[$SAMPLE] aligning ..."
  STAR --runThreadN "$THREADS" \
    --genomeDir "$INDEX" \
    --readFilesIn "$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" "$FASTQ_DIR/${SAMPLE}_R2.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$OUT/${SAMPLE}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS NM MD \
    --outSAMattrRGline ID:${SAMPLE} LB:lib1 PL:ILLUMINA SM:${SAMPLE} PU:unit1 \
    --sjdbGTFfile "$GTF" \
    --twopassMode Basic \
    --limitBAMsortRAM 31000000000 \
    > "$BASE/logs/${SAMPLE}_star.log" 2>&1
  mv "$OUT/${SAMPLE}_Aligned.sortedByCoord.out.bam" "$OUT/${SAMPLE}_sorted.bam"
  samtools index "$OUT/${SAMPLE}_sorted.bam"
  grep -E "Uniquely mapped reads %|% of reads mapped to multiple loci" "$OUT/${SAMPLE}_Log.final.out" | sed "s/^/  [$SAMPLE] /"
done
