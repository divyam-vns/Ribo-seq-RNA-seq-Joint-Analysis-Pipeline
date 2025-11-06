#!/usr/bin/env Rscript
# Usage: Rscript 04_ribo_rscript.R config.yaml
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1) stop("Usage: Rscript 04_ribo_rscript.R config.yaml")
library(yaml); library(riboWaltz); library(data.table); library(tidyverse)

cfg <- yaml.load_file(args[1])
proj <- cfg$project_dir
gtf <- cfg$gtf
fastqdir <- cfg$fastq_dir
outdir <- file.path(proj, "step04_ribowaltz")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

samples <- read_tsv(cfg$sample_table, col_types = cols())

# For each Ribo BAM, produce P-site offsets and per-gene P-site counts
# Note: riboWaltz expects BAM -> reads, ideally paired-end or single-end; here we assume single-end footprints.

cat("Processing Ribo BAMs for P-site assignment and periodicity...\n")
ribo_bams <- list()
for(i in 1:nrow(samples)) {
  ribo_fastq <- samples$ribo_fastq[i]
  base <- sub(".fastq.gz$","", basename(ribo_fastq))
  bam <- file.path(proj, "step02_star", base, "Aligned.sortedByCoord.out.bam")
  ribo_bams[[samples$sample_id[i]]] <- bam
}

# Build reads list for riboWaltz (this step needs BAM -> BED-like read info with read length, start,end,strand)
# We'll use a simple approach: generate a reads frame using readBam from GenomicAlignments (if available).
library(GenomicAlignments); library(Rsamtools)
reads_list <- list()
for(s in names(ribo_bams)) {
  bamfile <- ribo_bams[[s]]
  if(!file.exists(bamfile)) stop(paste("Missing:", bamfile))
  gal <- readGAlignments(bamfile)
  # convert to data.frame for riboWaltz
  df <- data.frame(seqnames=as.character(seqnames(gal)),
                   start=start(gal),
                   end=end(gal),
                   strand=as.character(strand(gal)))
  # keep fragment length
  df$length <- width(gal)
  reads_list[[s]] <- df
  cat("Loaded", s, "num reads:", nrow(df), "\n")
}

# Load CDS annotation from GTF (riboWaltz expects transcripts list)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(gtf, format="gtf")
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

# compute P-site offsets per read length
psite_offsets <- riboWaltz::psite(reads_list, cds)
saveRDS(psite_offsets, file=file.path(outdir,"psite_offsets.rds"))

# Periodicity plots and read-length distributions
pdf(file.path(outdir,"ribo_qc_plots.pdf"))
riboWaltz::length_dist(reads_list)
riboWaltz::frame_psite(psite_offsets)
dev.off()

# Generate P-site counts per gene (using psite() outputs)
# riboWaltz has functions to get counts per transcript region; here we extract psite_df
psite_df <- riboWaltz::psite_df(psite_offsets)
write_tsv(psite_df, file.path(outdir,"psite_table.tsv"))
cat("Wrote psite_table.tsv\n")

# Optionally: aggregate psite per gene (map transcripts -> gene via GTF attributes)
# Create simple gene mapping
gff <- rtracklayer::import(gtf)
tx2gene <- as.data.frame(gff[gff$type=="transcript", c("transcript_id","gene_id")])
# If psite has transcript names, perform mapping
if("transcript" %in% colnames(psite_df)) {
  psite_gene <- psite_df %>% left_join(tx2gene, by=c("transcript"="transcript_id"))
  write_tsv(psite_gene, file.path(outdir,"psite_table_gene_mapped.tsv"))
}
cat("04_ribo_rscript.R done\n")
