#!/usr/bin/env bash
# Usage: bash 03_counts_featureCounts.sh config.yaml
CONFIG=$1
if [ -z "$CONFIG" ]; then
  echo "Usage: $0 config.yaml"; exit 1
fi

project_dir=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['project_dir'])
PY)
gtf=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['gtf'])
PY)
threads=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['featureCounts_threads'])
PY)

outdir=${project_dir}/step03_counts
mkdir -p ${outdir}

# collect BAMs
python - <<PY
import yaml, pandas as pd, os
cfg=yaml.safe_load(open("$CONFIG"))
df=pd.read_csv(cfg['sample_table'], sep='\\t')
ribo_bams=[]
rna_bams=[]
for idx,row in df.iterrows():
    for col,arr in [('ribo_fastq','ribo_bams'),('rna_fastq','rna_bams')]:
        base=os.path.basename(row[col]).replace('.fastq.gz','')
        bam=os.path.join(cfg['project_dir'],'step02_star', base, 'Aligned.sortedByCoord.out.bam')
        if not os.path.exists(bam):
            print("ERROR: missing bam:", bam)
            sys.exit(1)
        if col=='ribo_fastq': ribo_bams.append(bam)
        else: rna_bams.append(bam)
# write file lists
open(os.path.join(cfg['project_dir'],'step03_counts','ribo_bams.txt'),'w').write("\\n".join(ribo_bams))
open(os.path.join(cfg['project_dir'],'step03_counts','rna_bams.txt'),'w').write("\\n".join(rna_bams))
print("Wrote lists to step03_counts/")
PY

# featureCounts for RNA (gene-level)
featureCounts -T ${threads} -a ${gtf} -o ${outdir}/rna_gene_counts.txt $(cat ${project_dir}/step03_counts/rna_bams.txt)

# featureCounts for Ribo (CDS exons only) -- requires GTF filtered to CDS or use -t CDS -g gene_id
featureCounts -T ${threads} -t CDS -g gene_id -a ${gtf} -o ${outdir}/ribo_cds_counts.txt $(cat ${project_dir}/step03_counts/ribo_bams.txt)

echo "Step 3 complete. Counts in ${outdir}"
