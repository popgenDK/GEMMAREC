#!/bin/bash
# Small examples of when GEMMAREC (gemmarec1) and gemmaPopgen (-model rec/dom)
# give identical results and when they differ.
# Uses 250 SNPs of chromosome 1 of gemmaPopgen/example/mouse_hs1940 (a few seconds).
# Needs plink 1.9. Run from the repository root after building the binaries:
#   examples/compare_versions.sh
# Binary locations can be changed with GEMMAREC1 and GEMMAPOPGEN.
set -e
ROOT=$(pwd)
R1=${GEMMAREC1:-$ROOT/bin/gemmarec1}
P=${GEMMAPOPGEN:-$ROOT/gemmaPopgen/bin/gemmaPopgen}
B=$ROOT/gemmaPopgen/example/mouse_hs1940
W=$(mktemp -d)
cd $W
q() { "$@" > /dev/null 2>&1; }

# 250 SNPs, and two pre-polarised copies: A1 = minor allele, and A1 = major
# allele, in the analysed individuals (those with phenotype 1)
awk '$1==1 && $4>0' $B.bim | head -250 | cut -f2 > snps.txt   # plink skips SNPs without a position
awk '$6!="NA" && $6!="-9" {print $1, $2}' $B.fam > analysed.txt
plink --bfile $B --extract snps.txt --keep analysed.txt --nonfounders --freq --out frq > /dev/null
awk 'NR>1 {print $2, $3}' frq.frq > minor.txt   # plink's A1 is the minor allele
awk 'NR>1 {print $2, $4}' frq.frq > major.txt
for a in minor major; do
  plink --bfile $B --extract snps.txt --a1-allele $a.txt 2 1 --make-bed --out $a > /dev/null
  cp $B.fam $a.fam   # keep all phenotype columns
done
cp $B.bed orig.bed; cp $B.bim orig.bim; cp $B.fam orig.fam
O="-bfile orig -snps snps.txt"   # original file (A1 as in the .bim)
echo "SNPs: $(wc -l < minor.txt); A1 is the major allele in the original file for" \
     "$(awk 'NR==FNR{m[$1]=$2; next} ($2 in m) && m[$2]!=$5' minor.txt orig.bim | wc -l) of them"

# GRMs: GEMMA's (additive, made by gemmaPopgen) and GEMMAREC's (recessive-coded)
q $P  $O -gk 1 -o K_gemma
q $R1 $O -gk 1 -o K_rec
KG=output/K_gemma.cXX.txt
KR=output/K_rec.cXX.txt

# compare two association files by rs: SNPs in each, shared SNPs, and on how
# many shared SNPs the tested allele and all statistics (beta ... p_score) are
# identical. With "flip", the second file must have the other allele as A1 and
# the opposite sign of beta (dominant model for the minor allele = 1 - recessive
# model for the major allele).
cmp() {
  awk -v lab="$1" -v flip="$4" 'BEGIN{FS="\t"}
    FNR==1 { next }
    FILENAME==ARGV[1] { na++; r[$2]=(flip ? $6 : $5); b=$8; if (flip) b=(b ~ /^-/) ? substr(b,2) : "-" b
                        r[$2]=r[$2] FS b; for (i=9;i<=15;i++) r[$2]=r[$2] FS $i; next }
    { nb++; if ($2 in r) { ns++; s=$5; for (i=8;i<=15;i++) s=s FS $i; if (s==r[$2]) nid++ } }
    END { printf "%-66s SNPs %3d vs %3d, identical %3d of %3d shared\n", lab, na, nb, nid, ns }' \
    output/$2.assoc.txt output/$3.assoc.txt
}

q $R1 -bfile minor -k $KG -maf 0 -lmm 4 -o rec1_minor
q $R1 -bfile major -k $KG -maf 0 -lmm 4 -o rec1_major
q $R1 $O -k $KG -maf 0 -lmm 4 -o rec1
q $P  $O -k $KG -maf 0 -lmm 4 -model rec -minor -o pop_rec_minor
q $P  $O -k $KG -maf 0 -lmm 4 -model dom -minor -o pop_dom_minor
q $P  $O -k $KG -maf 0 -lmm 4 -model rec -o pop_rec
echo
echo "== IDENTICAL results (same GRM, -maf 0)"
cmp "1. rec: gemmarec1 on file with A1 = minor vs gemmaPopgen -minor" rec1_minor pop_rec_minor
cmp "2. dom: gemmarec1 on file with A1 = major vs gemmaPopgen -minor" rec1_major pop_dom_minor flip
cmp "3. rec: gemmarec1 vs gemmaPopgen, both on the original file" rec1 pop_rec

echo
echo "== DIFFERENT results"
cmp "4. rec: gemmarec1 on original file vs gemmaPopgen -minor" rec1 pop_rec_minor
q $R1 -bfile minor -k $KG -lmm 4 -o rec1_minor_maf
q $P  $O -k $KG -lmm 4 -model rec -minor -o pop_rec_minor_maf
cmp "5. as 1, but default -maf 0.01" rec1_minor_maf pop_rec_minor_maf
q $R1 -bfile minor -k $KR -maf 0 -lmm 4 -o rec1_minor_ownK
cmp "6. as 1, but gemmarec1 uses its own GRM (-gk)" rec1_minor_ownK pop_rec_minor

echo
echo "Frequency columns for the first SNP (gemmarec1 af = freq_HO/2; gemmaPopgen af = A1 freq, freq_nonref = freq_HO):"
head -2 output/rec1_minor.assoc.txt | cut -f2,5,6,7 | sed 's/^/  gemmarec1    /'
head -2 output/pop_rec_minor.assoc.txt | awk 'BEGIN{FS=OFS="\t"} {print $2,$5,$6,$7,$NF}' | sed 's/^/  gemmaPopgen  /'
cd $ROOT && rm -rf $W
