#!/bin/bash
#SBATCH -J cfr10kb_mef
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 168:00:00
#SBATCH -o cfr10kb_mef.%j.out
#SBATCH -e cfr10kb_mef.%j.err

set -euo pipefail

ASMDIR="/public/home/220322096/FYM/cfr-SRR3/2499/cfr10kb/split_short_fasta"
OUTDIR="/public/home/220322096/FYM/cfr-SRR3/2499/cfr10kb/MGE_out"
THREADS=16

mkdir -p "$OUTDIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mefinder

shopt -s nullglob

total=$(find "$ASMDIR" -maxdepth 1 -name "*.fasta" | wc -l)
i=0

for fa in "$ASMDIR"/*.fasta; do
  i=$((i+1))
  sample=$(basename "${fa%.fasta}")

  # mefinder 会生成多个文件，这里用 gff 前缀作为是否完成的标记
  if [ -s "$OUTDIR/$sample.gff" ]; then
    echo "[$(date)] [$i/$total] skip existing: $sample"
    continue
  fi

  echo "[$(date)] [$i/$total] RUN $sample"
  mefinder find --contig "$fa" -t $THREADS --gff "$OUTDIR/$sample"
done

echo "[$(date)] all done"
