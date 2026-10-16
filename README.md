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
| output columns | GEMMA's columns (with GEMMA's `af`) + last column `freq_nonref` | `af` = `freq_HO/2` (do not use) |
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
* The output has GEMMA's columns plus a last column `freq_nonref`: the
  frequency of individuals coded 1 (A1/A1 for rec, A1 carriers for dom and add).
* `-hom-maf X` also removes SNPs with `freq_nonref` below X or above 1−X (rec/dom).
* `-gk` ignores `-model`, `-minor` and `-hom-maf`.

See [`gemmaPopgen/README.md`](gemmaPopgen/README.md) for details and tests.

## Examples: when GEMMAREC and gemmaPopgen give the same results

[`examples/compare_versions.sh`](examples/compare_versions.sh) runs all three
versions on 250 SNPs of the mouse example data and compares the outputs. It
takes a few seconds. These are its results:

| | comparison | SNPs | beta, p_wald |
|---|---|---|---|
| **same** | 1. recessive, same GRM (`-k`), `-maf 0`: gemmarec1 vs `gemmaPopgen -model rec` | 212 = 212 | identical |
| **same** | 2. dominant, same GRM, `-maf 0`: `gemmarec2 -dom` vs `gemmaPopgen -model dom` | 234 = 234 | identical |
| differ | 3. default `-maf 0.01` | 190 vs 212 | identical for the 190 shared SNPs |
| differ | 4. each program uses its own GRM from `-gk` | 212 = 212 | different (p_wald up to 0.96 apart) |
| differ | 5. `gemmaPopgen -minor` vs gemmarec1 | 212 = 212 | different for 4 SNPs |

So the association test is the same in all versions. The results differ only
because of:

* **SNP filter (3).** GEMMAREC's `-maf` filters on `freq_HO/2`, while
  gemmaPopgen's `-maf` uses the real minor allele frequency, as GEMMA does.
  Here gemmaPopgen keeps 22 extra SNPs, which have MAF ≥ 0.01 but
  `freq_HO/2` < 0.01.
* **GRM (4).** `gemmarec1 -gk` builds the GRM from recessive-coded genotypes,
  while gemmaPopgen always builds GEMMA's additive GRM. Here the GRM uses only
  250 SNPs, which makes the difference large. Pass the same GRM with `-k` to get
  identical results.
* **Tested allele (5).** GEMMAREC always tests A1 from the `.bim`. With
  `-minor`, gemmaPopgen tests the minor allele. In these data A1 is the major
  allele for 4 SNPs (e.g. rs6316488), so those SNPs test different homozygotes.
  Without `-minor` the results are the same.

The frequency columns also differ. For the first SNP (A1 = A), gemmarec1 prints
`af` = 0.093, which is `freq_HO/2`. gemmaPopgen prints GEMMA's `af` = 0.443 (the
A1 allele frequency) and `freq_nonref` = 0.185 (the A1/A1 frequency).

The script uses `bin/gemmarec1`, `GEMMAREC2/bin/gemmarec2` and
`gemmaPopgen/bin/gemmaPopgen`. Other locations can be given with the
`GEMMAREC1`, `GEMMAREC2` and `GEMMAPOPGEN` environment variables. If the
pre-built `bin/gemmarec1` does not run (see below), build it in a copy of this
folder and point `GEMMAREC1` to it.

## GEMMAREC (original)

This is a modified version of [GEMMA v. 0.98.5](https://github.com/genetics-statistics/GEMMA/releases/tag/v0.98.5).
The software is modified to perform association tests with a recessive model
instead of the default additive. See `changesToRecessive.txt` for what has been
changed. Note that the allele frequencies in the output should not be used: `af`
is half the frequency of A1/A1 homozygotes.

The pre-built binary `bin/gemmarec1` was linked against GSL 2.5
(`libgsl.so.23`) and does not run on systems with a newer GSL. In that case,
run `make` in this folder to build it again.
