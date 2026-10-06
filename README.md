# GEMMAREC: GEMMA with recessive (and dominant) association models

This repository contains two versions:

* **GEMMAREC2** (folder [`GEMMAREC2/`](GEMMAREC2)) is the recommended version. It is based on GEMMA 0.98.6 and supports recessive and dominant models, reports correct genotype frequencies, and can test the minor allele.
* **GEMMAREC** (this top folder) is the original version, based on GEMMA 0.98.5. It supports only the recessive model.

Both versions give identical association results for single-trait GWAS (`-gk`, `-lm`, univariate `-lmm`). The models are only tested for single-trait GWAS.

## GEMMAREC2: quick start

```sh
cd GEMMAREC2
make                     # needs GSL and OpenBLAS; binary is bin/gemmarec2

# recessive model (default), testing the minor allele
bin/gemmarec2 -bfile data -minor -gk 1 -o kin
bin/gemmarec2 -bfile data -minor -k output/kin.cXX.txt -lmm 4 -o rec

# dominant model: add -dom (also when making the kinship matrix)
bin/gemmarec2 -bfile data -minor -dom -gk 1 -o kin_dom
bin/gemmarec2 -bfile data -minor -dom -k output/kin_dom.cXX.txt -lmm 4 -o dom
```

* Input must be PLINK (`-bfile`). Genotypes are recoded as A1/A1 → 1, others → 0 (recessive), or A1/A1 and A1/A2 → 1 (dominant, `-dom`). BIMBAM input is not recoded.
* `-minor` tests the minor allele. Without it, A1 from the `.bim` file is tested, as in GEMMAREC. `allele1` in the output is the tested allele.
* The output column `freq_HO` (recessive) or `freq_carrier` (dominant) replaces `af`.
* `-hom-maf X` removes SNPs with `freq_HO` (or `freq_carrier`) below X or above 1-X. `-maf` works as in GEMMAREC: it filters on `freq_HO/2`. Use `-maf 0 -hom-maf X` to filter only on the genotype frequency.

See [`GEMMAREC2/README.md`](GEMMAREC2/README.md) for details and tests.

## GEMMAREC (original)

This is a modified version of [GEMMA v. 0.98.5](https://github.com/genetics-statistics/GEMMA/releases/tag/v0.98.5). The software is modified to perform association test with a recessive model instead of the default additive. See changesToRecessive.txt to see what has been changed. Note that the allele frequencies in the output should not be used! (In the output, `af` is half the frequency of A1/A1 homozygotes.) The binary is `bin/gemmarec1` (run `make` in this folder to build it).
