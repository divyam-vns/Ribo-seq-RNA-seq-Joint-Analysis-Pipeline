#!/usr/bin/env Rscript
# Usage: Rscript 05_deseq2_interaction.R config.yaml
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1) stop("Usage: Rscript 05_deseq2_interaction.R config.yaml")
library(yaml); library(DESeq2); library(tidyverse); library(data.table); library(pheatmap)

cfg <- yaml.load_file(args[1])
proj <- cfg$project_dir
outdir <- file.path(proj,"step05_deseq2")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Read counts from featureCounts outputs
rna_counts_file <- file.path(proj, "step03_counts", "rna_gene_counts.txt")
ribo_counts_file <- file.path(proj, "step03_counts", "ribo_cds_counts.txt")

# Parse featureCounts outputs (skip headers)
read_counts <- function(f) {
  df <- fread(f, skip=1, header=TRUE)
  # featureCounts has first 6 columns as meta; geneID in column 1
  colnames(df)[1] <- "GeneID"
  return(df)
}
rna_df <- read_counts(rna_counts_file)
ribo_df <- read_counts(ribo_counts_file)

# Extract counts matrix
get_counts_mat <- function(df) {
  # columns after 6th are sample counts
  mat <- as.matrix(df[,7:ncol(df)])
  rownames(mat) <- df$GeneID
  return(mat)
}
rna_mat <- get_counts_mat(rna_df)
ribo_mat <- get_counts_mat(ribo_df)

# Ensure column ordering matches sample table
samples <- read_tsv(cfg$sample_table, col_types = cols())
# Build sample metadata: one column per sample per assay
# We'll create combined matrices where columns are labeled sampleID_assay
rna_ids <- paste0(samples$sample_id, "_RNA")
ribo_ids <- paste0(samples$sample_id, "_Ribo")
# Rename count columns to these IDs (assumes featureCounts order matches sample order)
colnames(rna_mat) <- rna_ids
colnames(ribo_mat) <- ribo_ids

# Combine counts: counts_combined = cbind(ribo_mat, rna_mat)
# Ensure genes consistent
common_genes <- intersect(rownames(rna_mat), rownames(ribo_mat))
rna_mat <- rna_mat[common_genes, , drop=FALSE]
ribo_mat <- ribo_mat[common_genes, , drop=FALSE]
counts_combined <- cbind(ribo_mat, rna_mat)

# Build colData
n <- nrow(samples)
cond <- samples$condition
assay <- rep(c("Ribo","RNA"), each=n)
sample_colnames <- c(ribo_ids, rna_ids)
colData <- data.frame(
  sample = sample_colnames,
  condition = rep(cond, 2),
  assay = assay
)
rownames(colData) <- colData$sample

# DESeq2 design with interaction
dds <- DESeqDataSetFromMatrix(countData = counts_combined,
                              colData = colData,
                              design = ~ condition + assay + condition:assay)

# Prefilter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)

# The interaction term name depends on factor levels; get resultsNames
print(resultsNames(dds))
# Find the name that corresponds to the interaction, e.g., "condition_trt.assay_Ribo"
rnames <- resultsNames(dds)
int_name <- rnames[grep(":", rnames)]
if(length(int_name)==0) {
  # fallback: construct
  int_name <- paste0("condition", unique(colData$condition)[2], ".assayRibo")
}
res_int <- results(dds, name=int_name)
res_table <- as.data.frame(res_int) %>% rownames_to_column("GeneID")
write_tsv(res_table, file.path(outdir,"differential_translation_interaction_results.tsv"))

# MA plot for interaction
pdf(file.path(outdir,"interaction_MA.pdf"))
plotMA(res_int, main="Interaction term (differential translation)", ylim=c(-2,2))
dev.off()

# Top genes heatmap (top 50 by padj)
top <- res_table %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(50) %>% pull(GeneID)
mat <- assay(rlog(dds))[top, ]
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, filename=file.path(outdir,"top50_heatmap.png"))

cat("DESeq2 interaction analysis complete. Results in", outdir, "\n")
