#!/usr/bin/env bash
# Usage: bash 01_trim_and_rRNA_removal.sh config.yaml

CONFIG=$1
if [ -z "$CONFIG" ]; then
  echo "Usage: $0 config.yaml"; exit 1
fi

# parse YAML (simple)
project_dir=$(python - <<PY
import sys,yaml
print(yaml.safe_load(open("$CONFIG"))['project_dir'])
PY)
fastq_dir=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['fastq_dir'])
PY)
rrna_index=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['rrna_index'])
PY)
adapter_seq=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['adapter_seq'])
PY)
threads=$(python - <<PY
import yaml
print(yaml.safe_load(open("$CONFIG"))['threads'])
PY)

outdir=${project_dir}/step01_trim
mkdir -p ${outdir}

# sample table
samples=$(python - <<PY
import yaml, pandas as pd
cfg=yaml.safe_load(open("$CONFIG"))
df=pd.read_csv(cfg['sample_table'], sep='\\t')
# collect all fastq columns
files = set(df['rna_fastq'].tolist() + df['ribo_fastq'].tolist())
print('\\n'.join(files))
PY)

echo "Trimming and rRNA removal; output -> $outdir"

for f in $samples; do
  fq=${fastq_dir}/${f}
  base=$(basename ${f} .fastq.gz)
  trimmed=${outdir}/${base}.trimmed.fastq.gz
  no_rRNA=${outdir}/${base}.norrna.fastq.gz

  echo "Processing $fq -> $trimmed"
  cutadapt -a ${adapter_seq} -m 15 -o ${trimmed} ${fq}

  echo "Removing rRNA via bowtie2 -> $no_rRNA"
  bowtie2 -p ${threads} -x ${rrna_index} -U ${trimmed} --un ${no_rRNA} -S /dev/null
done

echo "Step 1 complete."
