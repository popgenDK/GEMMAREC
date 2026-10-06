#!/bin/bash
# -minor: allele-swapped chr1 vs R reference (flip to minor allele, 0/1 coding as BIMBAM)
# and vs the original chr1 .bed
set -e
G=${GEMMA:-../../bin/gemma}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8}
for m in rec dom; do
  f=""; [ $m = dom ] && f="-dom"
  for b in swap chr1; do
    $G -bfile $b $f -minor -gk 1 -o ${b}_${m}_kin >/dev/null
    $G -bfile $b $f -minor -lm 4 -o ${b}_${m}_lm >/dev/null
    $G -bfile $b $f -minor -k output/swap_${m}_kin.cXX.txt -lmm 4 -o ${b}_${m}_lmm >/dev/null
  done
  $G -g ${m}.geno.txt -p pheno.txt -gk 1 -o bb_${m}_kin >/dev/null
  $G -g ${m}.geno.txt -p pheno.txt -lm 4 -o bb_${m}_lm >/dev/null
  $G -g ${m}.geno.txt -p pheno.txt -k output/swap_${m}_kin.cXX.txt -lmm 4 -o bb_${m}_lmm >/dev/null
done
$G -bfile swap -minor -lm 4 -maf 0 -hom-maf 0.05 -o swap_rec_hommaf >/dev/null
