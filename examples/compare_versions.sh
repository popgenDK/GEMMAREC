#!/bin/bash
# Small examples of when GEMMAREC (gemmarec1), GEMMAREC2 (gemmarec2) and
# gemmaPopgen give the same results and when they differ.
# Uses 250 SNPs of chromosome 1 of gemmaPopgen/example/mouse_hs1940 (a few seconds).
# Run from the repository root after building the binaries:
#   examples/compare_versions.sh
# Binary locations can be changed with GEMMAREC1, GEMMAREC2 and GEMMAPOPGEN.
set -e
ROOT=$(pwd)
R1=${GEMMAREC1:-$ROOT/bin/gemmarec1}
R2=${GEMMAREC2:-$ROOT/GEMMAREC2/bin/gemmarec2}
P=${GEMMAPOPGEN:-$ROOT/gemmaPopgen/bin/gemmaPopgen}
B=$ROOT/gemmaPopgen/example/mouse_hs1940
W=$(mktemp -d)
cd $W
q() { "$@" > /dev/null 2>&1; }

awk '$1==1' $B.bim | head -250 | cut -f2 > snps.txt
S="-bfile $B -snps snps.txt"

# GRMs: GEMMA's (additive) and GEMMAREC's (recessive-coded)
q $P  $S -gk 1 -o K_gemma
q $R1 $S -gk 1 -o K_rec
KG=output/K_gemma.cXX.txt
KR=output/K_rec.cXX.txt

# compare two association files by rs: number of SNPs in each, shared SNPs,
# and the largest difference in beta and p_wald on shared SNPs
cmp() {  # cmp label fileA fileB [sign of beta in B]
  awk -v lab="$1" -v sgn="${4:-1}" 'BEGIN{FS="\t"}
    FNR==1 { for (i=1;i<=NF;i++) col[FILENAME,$i]=i; next }
    FILENAME==ARGV[1] { na++; b[$2]=$(col[FILENAME,"beta"]); p[$2]=$(col[FILENAME,"p_wald"]); next }
    { nb++; if ($2 in b) { ns++
        db=b[$2]-sgn*$(col[FILENAME,"beta"]); if (db<0) db=-db; if (db>mb) mb=db
        dp=p[$2]-$(col[FILENAME,"p_wald"]); if (dp<0) dp=-dp; if (dp>mp) mp=dp } }
    END { printf "%-52s SNPs %3d vs %3d, shared %3d, max|diff| beta %.2g, p_wald %.2g\n", lab, na, nb, ns, mb, mp }' \
    output/$2.assoc.txt output/$3.assoc.txt
}

echo "== SAME results"
q $R1 $S -k $KG -maf 0 -lmm 4 -o rec1
q $P  $S -k $KG -maf 0 -lmm 4 -model rec -o pop_rec
cmp "1. rec, GEMMA's GRM, -maf 0: gemmarec1 vs gemmaPopgen" rec1 pop_rec
q $R2 $S -k $KG -maf 0 -lmm 4 -dom -o rec2_dom
q $P  $S -k $KG -maf 0 -lmm 4 -model dom -o pop_dom
cmp "2. dom, GEMMA's GRM, -maf 0: gemmarec2 vs gemmaPopgen" rec2_dom pop_dom

echo
echo "== DIFFERENT results"
q $R1 $S -k $KG -lmm 4 -o rec1_maf
q $P  $S -k $KG -lmm 4 -model rec -o pop_rec_maf
cmp "3. default -maf 0.01 (different SNP filter)" rec1_maf pop_rec_maf
q $R1 $S -k $KR -maf 0 -lmm 4 -o rec1_ownK
cmp "4. each with its own GRM (-gk), -maf 0" rec1_ownK pop_rec
q $P  $S -k $KG -maf 0 -lmm 4 -model rec -minor -o pop_rec_minor
cmp "5. gemmaPopgen with -minor vs gemmarec1" rec1 pop_rec_minor
echo
echo "Details:"
awk 'BEGIN{FS="\t"} FNR==1{next} FILENAME==ARGV[1]{a[$2]=1; next} !($2 in a){n++; if (n<=2) ex=ex" "$2} END{printf "  3: %d SNPs kept only by gemmaPopgen (real MAF >= 0.01 but freq_HO/2 < 0.01), e.g.%s\n", n, ex}' output/rec1_maf.assoc.txt output/pop_rec_maf.assoc.txt
awk 'BEGIN{FS="\t"} FNR==1{next} FILENAME==ARGV[1]{a[$2]=$5; next} ($2 in a) && a[$2]!=$5 {n++; if (n<=2) ex=ex" "$2" ("a[$2]"->"$5")"} END{printf "  5: -minor swapped A1 for %d SNPs where A1 was the major allele, e.g.%s\n", n, ex}' output/rec1.assoc.txt output/pop_rec_minor.assoc.txt
echo "  frequency columns for the first SNP (gemmarec1 af = freq_HO/2; gemmaPopgen af = A1 freq, freq_nonref = freq_HO):"
head -2 output/rec1.assoc.txt | cut -f2,5,6,7 | sed 's/^/    gemmarec1    /'
head -2 output/pop_rec.assoc.txt | awk 'BEGIN{FS=OFS="\t"} {print $2,$5,$6,$7,$NF}' | sed 's/^/    gemmaPopgen  /'
rm -rf $W
