# Bulk RNA-seq analysis: STAR alignment, gene counts, DESeq2, limma-voom and WGCNA

From paired-end FASTQ files to differential expression, co-expression modules and pathway enrichment,
for an experiment with **three genotypes × two conditions × three replicates** (18 libraries).

| Folder | Language | What it does |
|---|---|---|
| `alignment/` | bash | hg38 + GENCODE references, STAR index, STAR alignment, featureCounts → count tables |
| `analysis/` | R, Python | limma-voom, DESeq2 (incl. interaction contrasts), WGCNA, PyDESeq2 with PCA-weighted re-analysis; volcano/MA plots, Venn diagrams, Enrichr and GSEA |

## Repository structure

```
├── config/
│   ├── config.example.yaml        # template: copy to config.yaml (git-ignored) and edit
│   ├── samples_bulk.csv           # one row per library (replicate)
│   └── raw_data_links.example.tsv # template for scripts/link_raw_data.sh
├── alignment/                     # 00-03: references -> count tables
├── analysis/                      # 01-04: downstream analyses
├── utils/                         # shared config loaders (R, Python, bash)
├── scripts/link_raw_data.sh       # links generic sample IDs to your FASTQ files
├── envs/                          # conda environments + R package installer
└── data/, results/                # created locally, git-ignored
```

## Sample naming

| Factor | Labels used in the code |
|---|---|
| Genotype | `Control`, `MutA`, `MutB` |
| Treatment | `Vehicle`, `Treated` |
| Condition | `<genotype>_<treatment>`, e.g. `MutA_Treated` |
| Sample ID | `OrgB_<genotype>_<treatment>_<replicate>`, e.g. `OrgB_MutA_Treated_n2` |

Edit `config/samples_bulk.csv` to match your libraries.

## Setup

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>

cp config/config.example.yaml config/config.yaml        # then edit config/config.yaml

# optional: link your FASTQ files to the generic sample IDs
cp config/raw_data_links.example.tsv config/raw_data_links.tsv   # then edit
bash scripts/link_raw_data.sh

conda env create -f envs/bulk_rnaseq.yml     # alignment + counting
conda env create -f envs/bulk_python.yml     # PyDESeq2 notebook
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
bash alignment/02_star_align.sh            # sorted, indexed BAM per library; skips finished libraries
bash alignment/03_count_reads.sh           # featureCounts -> count tables in <bulk_rna_workdir>
```

`03_count_reads.sh` counts read pairs per gene (GENCODE `gene_id`), removes the Ensembl version suffix
and names the columns by `sample_id`. It writes:
- `counts.txt`: tab-separated, column `Geneid`, used by the limma-voom analysis.
- `gene_counts.csv`: columns `gene_id`, `gene_name` and the samples, used by DESeq2, WGCNA and PyDESeq2.

Set `featurecounts.strandedness` in the config to match your library kit (0 unstranded, 1 stranded,
2 reverse-stranded).

### 2. Downstream analysis

Place `Metadata.csv` in `<bulk_rna_workdir>`, with one row per sample. `Sample` must match the
count-table column names. Columns used:

| Column | Example | Used by |
|---|---|---|
| `Sample` | `OrgB_MutA_Treated_n2` | all |
| `Condition` | `MutA_Treated` | DESeq2, WGCNA |
| `Genotype`, `Treatment`, `Batch` | `MutA`, `Treated`, `batch1` | limma-voom (Genotype also WGCNA) |
| `Age`, `Sex`, `Type` | | WGCNA traits |

| Step | Script | What it does | Output folder |
|---|---|---|---|
| 1 | `analysis/01_limma_voom_analysis.R` | TMM + voom normalisation, QC, PCA, limma DE (`~ Genotype + Treatment + Batch`), volcano plots, Venn diagrams, Enrichr; SOM metagenes, sample MSTs, Reactome ORA | `<bulk_rna_workdir>` |
| 2 | `analysis/02_deseq2_analysis.R` | DESeq2 (`~ Condition`), VST QC, PCA + loadings, pairwise and interaction contrasts, volcano plots, Enrichr | `<deseq2_workdir>` |
| 3 | `analysis/03_wgcna.R` | WGCNA signed network on VST data, module merging, module-trait correlations, gene significance per condition | `<wgcna_workdir>` |
| 4 | `analysis/04_pydeseq2_pca_weighted_analysis.ipynb` | PyDESeq2, dispersion and PCA diagnostics, DE per contrast, Venn and MA plots, Enrichr/GSEA (gseapy), DE re-analysis on PC2/PC3-weighted counts | `<pydeseq2_workdir>` |

```bash
Rscript analysis/01_limma_voom_analysis.R
Rscript analysis/02_deseq2_analysis.R
Rscript analysis/03_wgcna.R
# then open analysis/04_pydeseq2_pca_weighted_analysis.ipynb (environment: bulk_python)
```

**Significance.** DE genes are called at **BH-adjusted p < 0.05** in all four analyses.

**DESeq2 contrasts.**
- Pairwise contrasts:
  - Treated vs Vehicle within each genotype.
  - MutA and MutB vs Control (Vehicle).
  - MutA vs MutB (Vehicle).
- Interaction contrasts: whether the treatment response differs between genotypes.

Both lists are defined in the "Define Contrasts" section of `02_deseq2_analysis.R`.

**WGCNA.** The soft-thresholding power (`soft_power`) and the module colours examined per condition are
set in the script. Choose them from the scale-free-fit plot and the module-trait heatmap of your data.

Enrichment needs an internet connection (Enrichr).

## Reproducibility

- Record exact package versions before publishing: `renv::snapshot()` (R) and
  `conda env export --no-builds -n <env> > envs/<env>.lock.yml`.
- Keep notebooks free of outputs: `pip install nbstripout && nbstripout --install`.
- Raw data, alignments and results are not part of this repository (see `.gitignore`).

## Related repositories

Other analyses of the same experiment:
- [multiomics-analysis](https://github.com/Kart-796/multiomics-analysis): 10x multiome (snRNA + snATAC), MultiVelo and SCENIC+
- [rnaseq-variant-calling](https://github.com/Kart-796/rnaseq-variant-calling): GATK4 RNA-seq variant calling
- [proteomics-analysis](https://github.com/Kart-796/proteomics-analysis): DIA-NN proteomics

## Citation

If you use this code, please cite: *<paper reference / DOI>*.

## License

MIT, see [LICENSE](LICENSE).
