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

| Same? | comparison | plink allele order |
|---|---|---|
| **yes** | 1. recessive: gemmarec1 vs `gemmaPopgen -model rec -minor` | minor |
| **yes** | 2. dominant: gemmarec1 vs `gemmaPopgen -model dom -minor`, with beta of opposite sign and alleles swapped | major |
| **yes** | 3. recessive: gemmarec1 vs `gemmaPopgen -model rec` (no `-minor`) | mixed |
| no | 4. recessive: gemmarec1 vs `gemmaPopgen -model rec -minor`. The 3 SNPs with A1 = major allele differ. | mixed |
| no | 5. as 1, with the default `-maf 0.01`. gemmaPopgen keeps 16 more SNPs; the 199 shared SNPs are identical. | minor |
| no | 6. as 1, but gemmarec1 uses its own GRM (`gemmarec1 -gk`). All SNPs differ. | minor |

"Same" means the same SNPs, the same tested allele and identical association
statistics.

[`examples/compare_versions.sh`](examples/compare_versions.sh) runs gemmarec1
(GEMMAREC) and gemmaPopgen on 250 SNPs of the mouse example data and compares
the outputs. It takes a few seconds and needs plink 1.9.

**To get the same results from gemmarec1 as from gemmaPopgen, make the GRM with
original GEMMA** (`gemma -gk 1`, or `gemmaPopgen -gk 1`, which gives the
identical GRM) and pass it to gemmarec1 with `-k`. A GRM made with
`gemmarec1 -gk` is built from recessive-coded genotypes and gives different
results. All examples use GEMMA's GRM and `-maf 0` unless stated otherwise.

"plink allele order" is the A1 allele of the PLINK file given to gemmarec1:
the minor allele, the major allele, or mixed (the original file, where A1 is
the major allele for 3 of the 250 SNPs). gemmaPopgen always gets the original
(mixed) file.

* **Dominant (2).** gemmarec1 has no dominant model, but carrying the minor
  allele is the opposite of being homozygous for the major allele. So
  gemmarec1 on a file with A1 = major gives the dominant test for the minor
  allele, with the opposite sign of beta.
* **Allele order (3, 4).** gemmarec1 always tests A1 from the `.bim`. To match
  `gemmaPopgen -minor`, the file must be pre-polarised so that A1 is the minor
  allele, for example with `plink --a1-allele`. Without `-minor`, gemmaPopgen
  also tests A1 from the `.bim`.
* **SNP filter (5).** GEMMAREC's `-maf` filters on `freq_HO/2`, while
  gemmaPopgen's `-maf` uses the real minor allele frequency, as GEMMA does.
* **GRM (6).** Here the GRM uses only 250 SNPs, which makes the difference
  large.

The frequency columns also differ. For the first SNP (A1 = A), gemmarec1 prints
`af` = 0.093, which is `freq_HO/2`. gemmaPopgen prints GEMMA's `af` = 0.443 (the
A1 allele frequency) and `freq_nonref` = 0.185 (the A1/A1 frequency).

The script uses `bin/gemmarec1` and `gemmaPopgen/bin/gemmaPopgen`. Other
locations can be given with the `GEMMAREC1` and `GEMMAPOPGEN` environment
variables. If the pre-built `bin/gemmarec1` does not run (see below), build it
in a copy of this folder and point `GEMMAREC1` to it.

## GEMMAREC (original)

This is a modified version of [GEMMA v. 0.98.5](https://github.com/genetics-statistics/GEMMA/releases/tag/v0.98.5).
The software is modified to perform association tests with a recessive model
instead of the default additive. See `changesToRecessive.txt` for what has been
changed. Note that the allele frequencies in the output should not be used: `af`
is half the frequency of A1/A1 homozygotes.

The pre-built binary `bin/gemmarec1` was linked against GSL 2.5
(`libgsl.so.23`) and does not run on systems with a newer GSL. In that case,
run `make` in this folder to build it again.
