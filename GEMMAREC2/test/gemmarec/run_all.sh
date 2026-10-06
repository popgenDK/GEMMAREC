#!/bin/bash
# GEMMAREC check on chromosome 1 of example/mouse_hs1940 (runs in ~1 min).
# PLINK input (recoded inside GEMMA) is compared to BIMBAM input holding the
# same 0/1 recessive/dominant coding made in R (BIMBAM is never recoded).
# Also tests -minor on a copy of chr1 with A1/A2 swapped.
# Needs plink 1.9 and R. Run from the repository root: test/gemmarec/run_all.sh
set -e
ROOT=$(pwd)
W=$(mktemp -d)
cd $W
plink --bfile $ROOT/example/mouse_hs1940 --chr 1 --keep-allele-order --make-bed --out chr1 >/dev/null
plink --bfile chr1 --keep-allele-order --recode A --out chr1 >/dev/null
cp $ROOT/example/mouse_hs1940.fam chr1.fam   # keep all phenotype columns
Rscript --vanilla $ROOT/test/gemmarec/mkbimbam.R
GEMMA=$ROOT/bin/gemma bash $ROOT/test/gemmarec/run_small.sh
Rscript --vanilla $ROOT/test/gemmarec/cmp_small.R

# -minor: swap A1/A2 of every SNP and check that -minor gives the R reference
# and the same results as the original file
mkdir -p minor && cd minor
cp ../chr1.bed ../chr1.bim ../chr1.fam ../pheno.txt .
awk '{print $2, $6}' chr1.bim > a2.txt
plink --bfile chr1 --a1-allele a2.txt 2 1 --make-bed --out swap >/dev/null
cp chr1.fam swap.fam
plink --bfile swap --keep-allele-order --recode A --out swap >/dev/null
Rscript --vanilla $ROOT/test/gemmarec/minor_mkref.R
GEMMA=$ROOT/bin/gemma bash $ROOT/test/gemmarec/minor_run.sh
Rscript --vanilla $ROOT/test/gemmarec/minor_cmp.R
cd ..
echo "work dir: $W"
