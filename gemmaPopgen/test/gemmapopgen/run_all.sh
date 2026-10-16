#!/bin/bash
# gemmaPopgen check on 250 SNPs of chromosome 1 of example/mouse_hs1940 (about 1 min).
# About half of the SNPs get A1 = major allele, so -minor is tested too.
# Needs plink 1.9 and R. Run from the gemmaPopgen directory:
#   test/gemmapopgen/run_all.sh
# Set GEMMA_PLAIN=/path/to/unmodified/gemma to also check that the GRM and the
# additive results are identical to plain GEMMA.
set -e
ROOT=$(pwd)
export GEMMAPOPGEN=$ROOT/bin/gemmaPopgen
W=$(mktemp -d)
cd $W
E=$ROOT/example/mouse_hs1940
awk '$1==1' $E.bim | head -250 | cut -f2 > snps.txt
plink --bfile $E --extract snps.txt --keep-allele-order --make-bed --out base >/dev/null
Rscript --vanilla $ROOT/test/gemmapopgen/prep.R
plink --bfile base --a1-allele a1_rnd.txt 2 1 --make-bed --out rnd >/dev/null
cp $E.fam rnd.fam   # keep the phenotypes
plink --bfile rnd --keep-allele-order --recode A --out rnd >/dev/null
cp $E.fam rnd.fam
awk '{print 1, $5}' rnd.fam > cov.txt
Rscript --vanilla $ROOT/test/gemmapopgen/mkref.R
bash $ROOT/test/gemmapopgen/run.sh
Rscript --vanilla $ROOT/test/gemmapopgen/cmp.R
[ -z "$GEMMA_PLAIN" ] && echo "(GEMMA_PLAIN not set: parts 1-2 compared gemmaPopgen with itself)"
echo "work dir: $W"
