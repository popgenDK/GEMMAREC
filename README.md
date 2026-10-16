# GEMMAREC: GEMMA with recessive and dominant association models

This repository contains three modified versions of
[GEMMA](https://github.com/genetics-statistics/GEMMA) for genome-wide
association with non-additive models.

| | **gemmaPopgen** (recommended) | **GEMMAREC** (original) |
|---|---|---|
| folder | [`gemmaPopgen/`](gemmaPopgen) | this top folder |
| binary | `gemmaPopgen/bin/gemmaPopgen` | `bin/gemmarec1` |
| based on | GEMMA 0.98.6 | GEMMA 0.98.5 |
| models | additive, recessive, dominant (`-model add\|rec\|dom`) | recessive |
| GRM (`-gk`) | always as in GEMMA (additive) | recoded (recessive) |
| `-maf` filters on | minor allele frequency, as in GEMMA | `freq_HO/2` |
| frequency column | `af` / `freq_HO` / `freq_carrier` | `af` = `freq_HO/2` (do not use) |
| test the minor allele | `-minor` | no (A1 from the `.bim`) |

A third version, **GEMMAREC2**, is in [`GEMMAREC2/`](GEMMAREC2). It reproduces
the GEMMAREC results with a corrected frequency column and adds a dominant
model. See [`GEMMAREC2/README.md`](GEMMAREC2/README.md).

All versions recode genotypes only for PLINK input (`-bfile`). A1 is the tested
allele (`allele1` in the output). **All models are only tested for single-trait
GWAS** (`-gk`, `-lm`, univariate `-lmm`).

## Which version to use

* **gemmaPopgen** for new analyses. The model changes only the association
  test. The GRM and the SNP filters (`-maf`, `-miss`, `-hwe`, `-r2`) are
  exactly as in GEMMA, so the same GRM is used for every model, and
  `-model add` gives GEMMA's results.
* **GEMMAREC2** to reproduce old GEMMAREC analyses; see
  [`GEMMAREC2/README.md`](GEMMAREC2/README.md).
* **GEMMAREC** is kept unchanged for reference.

## gemmaPopgen: quick start

```sh
cd gemmaPopgen
make                     # needs GSL and OpenBLAS; binary is bin/gemmaPopgen

bin/gemmaPopgen -bfile data -gk 1 -o grm                                            # GRM as in GEMMA
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -o add                     # additive (= GEMMA)
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -model rec -minor -o rec   # recessive
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -model dom -minor -o dom   # dominant
```

* `-minor` tests the minor allele. Without it, A1 from the `.bim` file is tested.
* `-hom-maf X` also removes SNPs with `freq_HO` (rec) or `freq_carrier` (dom)
  below X or above 1−X.
* `-gk` ignores `-model`, `-minor` and `-hom-maf`.

See [`gemmaPopgen/README.md`](gemmaPopgen/README.md) for details and tests.

## GEMMAREC (original)

This is a modified version of [GEMMA v. 0.98.5](https://github.com/genetics-statistics/GEMMA/releases/tag/v0.98.5).
The software is modified to perform association tests with a recessive model
instead of the default additive. See `changesToRecessive.txt` for what has been
changed. Note that the allele frequencies in the output should not be used: `af`
is half the frequency of A1/A1 homozygotes.

The pre-built binary `bin/gemmarec1` was linked against GSL 2.5
(`libgsl.so.23`) and does not run on systems with a newer GSL. In that case,
run `make` in this folder to build it again.
