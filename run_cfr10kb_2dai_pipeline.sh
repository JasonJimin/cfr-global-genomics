#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/public/home/220322096/FYM/cfr-SRR3/2499/cfr10kb"
cd "$WORKDIR"

echo "[$(date)] Step 1: rename headers"
rm -f cfr2d_header_map.tsv
awk '
BEGIN{n=0}
/^>/{
  n++
  old=substr($0,2)
  print "cfr2d_" n "\t" old >> "cfr2d_header_map.tsv"
  print ">cfr2d_" n
  next
}
{print}
' cfr_flanks_10kb.oriented.fasta > cfr_flanks_10kb.oriented.short.fa

echo "[$(date)] Step 2: prokka"
prokka \
  --outdir prokka_10kb \
  --prefix cfr2d \
  --locustag CFR2D \
  --cpus 8 \
  --force \
  cfr_flanks_10kb.oriented.short.fa

echo "[$(date)] Step 3: integron_finder"
integron_finder \
  --local-max \
  --func-annot \
  --union-integrases \
  --cpu 8 \
  --outdir integron_finder_out \
  cfr_flanks_10kb.oriented.short.fa

echo "[$(date)] Step 4: ISEScan"
isescan.py \
  --removeShortIS \
  --seqfile cfr_flanks_10kb.oriented.short.fa \
  --output isescan_out \
  --nthread 8

echo "[$(date)] Done"
