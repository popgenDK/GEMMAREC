#!/bin/bash
# gemmaPopgen runs for test/gemmapopgen/run_all.sh (run inside the work dir)
set -e
P=${GEMMAPOPGEN:?}
U=${GEMMA_PLAIN:-$P}   # plain GEMMA; without it, gemmaPopgen is used (its BIMBAM path is plain GEMMA)
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8}
q() { "$@" > /dev/null; }

# 1. GRM: plain GEMMA vs gemmaPopgen with each model
q $U -bfile rnd -gk 1 -o u_gk1
q $U -bfile rnd -gk 2 -o u_gk2
for m in add rec dom; do
  q $P -bfile rnd -model $m -gk 1 -o p_gk1_$m
  q $P -bfile rnd -model $m -minor -gk 2 -o p_gk2_$m
done
K=output/u_gk1.cXX.txt

# 2. additive model: plain GEMMA vs gemmaPopgen, without and with covariates
for n in nocov cov; do
  c=""; [ $n = cov ] && c="-c cov.txt"
  q $U -bfile rnd $c -lm 4 -o u_lm_$n
  q $U -bfile rnd $c -k $K -lmm 4 -o u_lmm_$n
  q $P -bfile rnd $c -lm 4 -o p_lm_$n
  q $P -bfile rnd $c -model add -k $K -lmm 4 -o p_lmm_$n
done

# 3. rec/dom: gemmaPopgen on PLINK vs R-coded 0/1 genotypes as BIMBAM (-maf 0; SNP set checked in R)
for m in rec dom; do
  for pol in a1 minor; do
    f=""; [ $pol = minor ] && f="-minor"
    for a in lm lmm; do
      k=""; [ $a = lmm ] && k="-k $K"
      q $P -bfile rnd -model $m $f $k -$a 4 -o p_${m}_${pol}_$a
      q $U -g ${m}_${pol}.geno.txt -p pheno.txt -maf 0 $k -$a 4 -o r_${m}_${pol}_$a
    done
  done
  q $P -bfile rnd -model $m -minor -c cov.txt -k $K -lmm 4 -o p_${m}_minor_cov
  q $U -g ${m}_minor.geno.txt -p pheno.txt -c cov.txt -maf 0 -k $K -lmm 4 -o r_${m}_minor_cov
  q $P -bfile rnd -model $m -minor -hom-maf 0.1 -lm 4 -o p_${m}_hommaf
done
