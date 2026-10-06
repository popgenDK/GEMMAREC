#!/bin/bash
# chr1 subset: PLINK input with recoding vs BIMBAM 0/1 input (no recoding = plain GEMMA)
set -e
G=${GEMMA:-../../bin/gemma}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8}
$G -bfile chr1 -gk 1 -o k_rec >/dev/null
$G -bfile chr1 -gk 1 -dom -o k_dom >/dev/null
for m in rec dom; do
  f=""; [ $m = dom ] && f="-dom"
  $G -bfile chr1 $f -lm 4 -o bed_${m}_lm >/dev/null
  $G -bfile chr1 $f -k output/k_rec.cXX.txt -lmm 4 -o bed_${m}_lmm >/dev/null
  $G -bfile chr1 $f -n 1 6 -k output/k_rec.cXX.txt -lmm 1 -o bed_${m}_mvlmm >/dev/null
  $G -g ${m}.geno.txt -p pheno.txt -lm 4 -o bb_${m}_lm >/dev/null
  $G -g ${m}.geno.txt -p pheno.txt -k output/k_rec.cXX.txt -lmm 4 -o bb_${m}_lmm >/dev/null
  $G -g ${m}.geno.txt -p pheno.txt -n 1 6 -k output/k_rec.cXX.txt -lmm 1 -o bb_${m}_mvlmm >/dev/null
  $G -bfile chr1 $f -maf 0 -hom-maf 0.05 -lm 4 -o bed_${m}_hommaf >/dev/null
done
