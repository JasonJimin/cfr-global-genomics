#!/bin/bash
#SBATCH -J is10kb
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -t 168:00:00
#SBATCH -o is10kb.%j.out
#SBATCH -e is10kb.%j.err

set -euo pipefail

INDIR="/public/home/220322096/FYM/cfr-SRR3/2499/cfr10kb/split_short_fasta"
OUTDIR="/public/home/220322096/FYM/cfr-SRR3/2499/cfr10kb/isescan_10kb_single"
THREADS=1

mkdir -p "$OUTDIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate amrfinder

shopt -s nullglob

total=$(find "$INDIR" -maxdepth 1 -name "*.fasta" | wc -l)
i=0

: > "$OUTDIR/hit_samples.txt"
: > "$OUTDIR/nohit_samples.txt"
: > "$OUTDIR/error_samples.txt"

for fa in "$INDIR"/*.fasta; do
  i=$((i+1))
  base="${fa##*/}"
  sample="${base%.fasta}"
  workdir="$OUTDIR/$sample"
  log="$workdir/${sample}.isescan.log"

  mkdir -p "$workdir"

  # 已经判定过就跳过
  if [ -f "$workdir/status.done" ]; then
    echo "[$(date)] [$i/$total] skip existing: $sample"
    continue
  fi

  echo "[$(date)] [$i/$total] RUN $sample"

  cp "$fa" "$workdir/"
  cd "$workdir"

  # 不加 --removeShortIS，先统计 complete + partial
  if isescan.py \
      --seqfile "${sample}.fasta" \
      --output out \
      --nthread $THREADS \
      > "$log" 2>&1; then
    :
  else
    # 有些情况下 isescan.py 仍会非 0，但日志能反映真实状态
    :
  fi

  # 判定状态
  if grep -q "No IS element was found" "$log"; then
    echo "$sample" >> "$OUTDIR/nohit_samples.txt"
    echo "no_hit" > "$workdir/status.done"
  elif find "$workdir/out" -maxdepth 1 -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.sum" -o -name "*.raw" \) | grep -q .; then
    echo "$sample" >> "$OUTDIR/hit_samples.txt"
    echo "hit" > "$workdir/status.done"
  elif grep -q "ISEScan ends" "$log"; then
    # 正常结束但没标准输出文件，先记为 no_hit_like
    echo "$sample" >> "$OUTDIR/nohit_samples.txt"
    echo "no_hit_like" > "$workdir/status.done"
  else
    echo "$sample" >> "$OUTDIR/error_samples.txt"
    echo "error" > "$workdir/status.done"
  fi

  cd "$INDIR"
done

hit=$(sort -u "$OUTDIR/hit_samples.txt" | wc -l)
nohit=$(sort -u "$OUTDIR/nohit_samples.txt" | wc -l)
err=$(sort -u "$OUTDIR/error_samples.txt" | wc -l)

{
  echo "Total fasta: $total"
  echo "ISEScan hit (complete or partial): $hit"
  echo "ISEScan no hit: $nohit"
  echo "ISEScan error: $err"
} > "$OUTDIR/isescan_10kb_summary.txt"

sort -u "$OUTDIR/hit_samples.txt" -o "$OUTDIR/hit_samples.txt"
sort -u "$OUTDIR/nohit_samples.txt" -o "$OUTDIR/nohit_samples.txt"
sort -u "$OUTDIR/error_samples.txt" -o "$OUTDIR/error_samples.txt"

echo "[$(date)] all done"
