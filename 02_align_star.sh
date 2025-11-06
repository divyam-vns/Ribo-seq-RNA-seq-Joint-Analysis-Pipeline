#!/usr/bin/env bash
# Usage: bash 02_align_star.sh config.yaml
CONFIG=$1
if [ -z "$CONFIG" ]; then
  echo "Usage: $0 config.yaml"; exit 1
fi

project_dir=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['project_dir'])
PY)
fastq_dir=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['fastq_dir'])
PY)
outdir=${project_dir}/step02_star
mkdir -p ${outdir}

star_genome=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['star_genome_dir'])
PY)
threads=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['threads'])
PY)

# iterate over sample table
python - <<PY
import yaml, pandas as pd, os, subprocess, sys
cfg=yaml.safe_load(open("$CONFIG"))
df=pd.read_csv(cfg['sample_table'], sep='\\t')
for idx,row in df.iterrows():
    for col in ['rna_fastq','ribo_fastq']:
        fq=os.path.join(cfg['fastq_dir'], row[col])
        base=os.path.basename(row[col]).replace('.fastq.gz','')
        outdir=os.path.join(cfg['project_dir'],'step02_star', base)
        os.makedirs(outdir, exist_ok=True)
        # STAR options (adjust for short Ribo reads if necessary)
        cmd = [
            "STAR","--runThreadN",str(cfg['threads']),
            "--genomeDir", cfg['star_genome_dir'],
            "--readFilesIn", fq,
            "--readFilesCommand","zcat",
            "--outFileNamePrefix", outdir + "/",
            "--outSAMtype","BAM","SortedByCoordinate",
            "--outReadsUnmapped","Fastx"
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd)
PY

echo "Step 2 complete. Aligned BAMs in ${project_dir}/step02_star"
