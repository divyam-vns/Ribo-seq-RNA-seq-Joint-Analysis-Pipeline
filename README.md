# Ribo-seq + RNA-seq Joint Analysis Pipeline

**Purpose:** Reproducible pipeline to analyze matched Ribo-seq and RNA-seq samples, compute Translation Efficiency (TE), and test differential translation using a DESeq2 interaction model. Includes QC steps (FastQC, read length, periodicity), rRNA removal, alignment, counting, P-site assignment (RiboWaltz), and statistical testing.

**Repo structure**
- `environment.yml` — conda environment with required tools
- `config.yaml` — user-editable parameters and sample table path
- `01_trim_and_rRNA_removal.sh` — trim adapters & remove rRNA
- `02_align_star.sh` — align reads with STAR
- `03_counts_featureCounts.sh` — count reads (RNA: gene/exon; Ribo: CDS/P-sites)
- `04_ribo_rscript.R` — RiboWaltz QC and P-site processing (generates P-site counts)
- `05_deseq2_interaction.R` — builds combined matrix, runs DESeq2 interaction model, and produces plots
- `LICENSE` — MIT

---

## Quickstart

1. Clone repo and edit `config.yaml`:
```bash
git clone <this-repo>
cd <this-repo>
# edit config.yaml with your paths, sample table, genome indexes, GTF, etc.

_ _ _

2. Create conda environment:

conda env create -f environment.yml
conda activate ribo_rna_env


3. Run pipeline (examples):

bash 01_trim_and_rRNA_removal.sh config.yaml
bash 02_align_star.sh config.yaml
bash 03_counts_featureCounts.sh config.yaml
Rscript 04_ribo_rscript.R config.yaml
Rscript 05_deseq2_interaction.R config.yaml

## Notes & tips

##Ensure sample table lists sample IDs and maps to Ribo/RNA FASTQs (see config.yaml example).
##Ribo and RNA samples must be matched to biological replicates.
##Adjust STAR parameters (read length, --sjdbOverhang) according to read length.
##For small footprint reads, P-site offsets depend on read length; 04_ribo_rscript.R uses RiboWaltz to compute offsets.
