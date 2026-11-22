# Bulk RNA-seq analysis: STAR alignment, gene counts and downstream analysis

From paired-end FASTQ files to normalised expression, PCA, SOM metagenes and pathway enrichment,
for an experiment with **three genotypes × two conditions**.

| Folder | Language | What it does |
|---|---|---|
| `alignment/` | bash | hg38 + GENCODE references, STAR index, STAR alignment, featureCounts → `counts.txt` |
| `analysis/` | R | Normalisation, QC plots, PCA, SOM metagenes, enrichment (clusterProfiler, ReactomePA, enrichR, fgsea) |

## Repository structure

```
├── config/
│   ├── config.example.yaml        # template: copy to config.yaml (git-ignored) and edit
│   ├── samples_bulk.csv           # one row per sequenced library
│   └── raw_data_links.example.tsv # template for scripts/link_raw_data.sh
├── alignment/                     # 00-03: references -> counts.txt
├── analysis/bulk_rna_analysis.R
├── utils/                         # shared config loaders (R, Python, bash)
├── scripts/link_raw_data.sh       # links generic sample IDs to your FASTQ files
├── envs/                          # conda environment + R package installer
└── data/, results/                # created locally, git-ignored
```

## Sample naming

| Factor | Labels used in the code |
|---|---|
| Genotype | `Control`, `MutA`, `MutB` |
| Condition | `Vehicle`, `Treated` |
| Sample ID | `OrgB_<genotype>_<condition>_<replicate>`, e.g. `OrgB_MutA_Treated_n2` |

Edit `config/samples_bulk.csv` to match your libraries.

## Setup

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>

cp config/config.example.yaml config/config.yaml        # then edit config/config.yaml

# optional: link your FASTQ files to the generic sample IDs
cp config/raw_data_links.example.tsv config/raw_data_links.tsv   # then edit
bash scripts/link_raw_data.sh

conda env create -f envs/bulk_rnaseq.yml
Rscript envs/install_r_packages.R
```

**Run every script from the repository root.** Scripts read `config/config.yaml`; to use another
file, set `PROJECT_CONFIG=/path/to/config.yaml`.

## Running the analysis

### 1. Alignment and counting

Input: `<fastq_dir>/<sample_id>_R1.fastq.gz` and `_R2.fastq.gz` for every row of `config/samples_bulk.csv`.

```bash
conda activate bulk_rnaseq
bash alignment/00_download_references.sh   # hg38 (UCSC) + GENCODE GTF
bash alignment/01_build_star_index.sh      # sjdbOverhang = read length - 1 (config: star.sjdb_overhang)
bash alignment/02_star_align.sh            # sorted, indexed BAM per sample; skips finished samples
bash alignment/03_count_reads.sh           # featureCounts -> <bulk_rna_workdir>/counts.txt
```

`03_count_reads.sh` counts read pairs per gene (GENCODE `gene_id`), removes the Ensembl version
suffix and names the columns by `sample_id`. Set `featurecounts.strandedness` in the config to match
your library kit (0 unstranded, 1 stranded, 2 reverse-stranded).

### 2. Downstream analysis

Place `Metadata.csv` (columns `Sample`, `Genotype`, `Treatment`, `Batch`; `Sample` must match the
column names of `counts.txt`) in `<bulk_rna_workdir>`, then:
```bash
Rscript analysis/bulk_rna_analysis.R
```
