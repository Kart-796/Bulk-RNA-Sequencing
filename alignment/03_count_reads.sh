#!/usr/bin/env bash
# =============================================================================
# Gene-level read counts with featureCounts, written in the two formats used by
# the analysis scripts (Ensembl IDs without version, one column per sample_id):
#   <bulk_rna_workdir>/counts.txt       tab-separated, column Geneid
#                                       (analysis/01_limma_voom_analysis.R)
#   <bulk_rna_workdir>/gene_counts.csv  columns gene_id, gene_name, samples
#                                       (02_deseq2, 03_wgcna, 04 notebook)
# The raw featureCounts output and summary are kept next to them.
#
# Usage:  bash alignment/03_count_reads.sh
# =============================================================================
set -euo pipefail

# ---- Configuration: edit config/config.yaml (run from the repository root) ----
REPO_DIR="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cfg() { python3 "$REPO_DIR/utils/get_config.py" "$@"; }
BASE="$(cfg path alignment_dir)"
WORKDIR="$(cfg path bulk_rna_workdir)"
GENCODE_RELEASE="$(cfg value star.gencode_release)"
THREADS="$(cfg value star.threads)"
STRAND="$(cfg value featurecounts.strandedness)"   # 0 unstranded, 1 stranded, 2 reverse
GTF="$BASE/references/gencode.v${GENCODE_RELEASE}.annotation.gtf"
mkdir -p "$WORKDIR"
read -r -a SAMPLES <<< "$(cfg samples bulk)"

BAMS=()
for SAMPLE in "${SAMPLES[@]}"; do BAMS+=("$BASE/aligned/${SAMPLE}_sorted.bam"); done

featureCounts -T "$THREADS" -p --countReadPairs -s "$STRAND" \
  -a "$GTF" -t exon -g gene_id \
  -o "$WORKDIR/featureCounts_raw.txt" "${BAMS[@]}"

# Tidy: drop Chr/Start/End/Strand/Length, strip Ensembl version, name columns by sample_id,
# and add gene names from the GTF for gene_counts.csv
python3 - "$WORKDIR/featureCounts_raw.txt" "$GTF" "$WORKDIR" "${SAMPLES[@]}" <<'PY'
import csv, re, sys
src, gtf, workdir, samples = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]

names = {}
with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if len(f) == 9 and f[2] == "gene":
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            gname = re.search(r'gene_name "([^"]+)"', f[8])
            if gid:
                names.setdefault(gid.group(1).split(".")[0], gname.group(1) if gname else "")

with open(src) as fin, open(f"{workdir}/counts.txt", "w") as ftxt, \
     open(f"{workdir}/gene_counts.csv", "w", newline="") as fcsv:
    w = csv.writer(fcsv)
    ftxt.write("\t".join(["Geneid"] + samples) + "\n")
    w.writerow(["gene_id", "gene_name"] + samples)
    for line in fin:
        if line.startswith("#") or line.startswith("Geneid"):
            continue
        f = line.rstrip("\n").split("\t")
        gid = f[0].split(".")[0]
        ftxt.write("\t".join([gid] + f[6:]) + "\n")
        w.writerow([gid, names.get(gid) or gid] + f[6:])
PY
echo "Counts written to $WORKDIR/counts.txt and $WORKDIR/gene_counts.csv"
