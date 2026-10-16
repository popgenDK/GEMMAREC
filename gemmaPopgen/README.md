# gemmaPopgen: GEMMA with additive, recessive and dominant models

gemmaPopgen is [GEMMA](https://github.com/genetics-statistics/GEMMA) 0.98.6
(master) with a choice of genetic model for single-trait association tests on
PLINK input. **The model only changes the association test.** The kinship
matrix (GRM) and the SNP filters are always computed exactly as in GEMMA, so
the GRM is the same for every model.

| `-model` | A1/A1 | A1/A2 | A2/A2 | last column `freq_nonref` |
|---|---|---|---|---|
| `add` (default) | 2 | 1 | 0 | A1 carrier frequency (A1/A1 + A1/A2) |
| `rec` | 1 | 0 | 0 | A1/A1 frequency |
| `dom` | 1 | 1 | 0 | A1 carrier frequency (A1/A1 + A1/A2) |

A1 is the tested allele (`allele1` in the output). By default it is the first
allele in the `.bim`. With `-minor` it is always the minor allele.

**Output.** The association output has exactly GEMMA's columns, including `af`
(the A1 allele frequency), plus one extra last column, `freq_nonref`. It is the
frequency of the non-reference individuals, i.e. those coded 1 in the model
(see the table). Frequencies are computed over the analysed individuals with a
non-missing genotype. `freq_nonref` is only added for PLINK input.

> **Only tested for single-trait GWAS** (`-lm`, univariate `-lmm`).

## Usage

```sh
make            # needs GSL and OpenBLAS; binary is bin/gemmaPopgen

# GRM, identical to GEMMA's (-model is not needed and is ignored here)
bin/gemmaPopgen -bfile data -gk 1 -o grm

# the same GRM for every model
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -o add
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -model rec -minor -o rec
bin/gemmaPopgen -bfile data -k output/grm.cXX.txt -lmm 4 -model dom -minor -o dom

# linear model without GRM
bin/gemmaPopgen -bfile data -lm 4 -model rec -minor -o rec_lm
```

| option | meaning |
|---|---|
| `-model add\|rec\|dom` | genetic model for `-lm` and `-lmm` (default `add`) |
| `-minor` | test the minor allele: A1/A2 are swapped for SNPs where A1 has frequency > 0.5 in the analysed individuals |
| `-hom-maf X` | with `rec`/`dom`: also remove SNPs with `freq_nonref` below X or above 1−X |

### SNP filtering

* `-maf`, `-miss`, `-hwe` and `-r2` work exactly as in GEMMA, on the additive
  genotypes, for all models. So `-maf` is a filter on the real minor allele
  frequency.
* With `rec`/`dom`, a SNP is also removed if the coded genotype does not vary
  (for example no A1/A1 individuals for `rec`). `-hom-maf` can remove rare
  coded genotypes. With covariates, the `-r2` check is also applied to the
  coded genotype.

### What `-model` and `-minor` do not change

* `-gk` ignores `-model`, `-minor` and `-hom-maf` and prints a note, so the GRM
  is always GEMMA's.
* With `-model add` and without `-minor`, the output is identical to GEMMA's,
  apart from the extra last column `freq_nonref`.
* `-model rec/dom` and `-minor` are only allowed with `-lm` and `-lmm` on
  PLINK input (`-bfile`). Other analyses (`-bslmm`, BIMBAM input, ...) give an error.
  The multivariate LMM also uses the model, but it is not tested.

## Differences from GEMMAREC / GEMMAREC2

GEMMAREC (`..`) and GEMMAREC2 (`../GEMMAREC2`) recode the
genotypes everywhere, including the GRM and `-maf`. Their `-maf` filters on
`freq_HO/2`, not on the minor allele frequency. gemmaPopgen recodes only for the
association test. With the same GRM (`-k`) and the same SNPs, the association
results of `-model rec` are the same as GEMMAREC2's. The SNP set can differ,
because `-maf` filters differently.

## Tests

`test/gemmapopgen/run_all.sh` takes about a minute and needs `plink` 1.9 and R.
It uses 250 SNPs of `example/mouse_hs1940`, about half with A1 = major allele.
With `GEMMA_PLAIN=/path/to/gemma`, it checks that:

1. the GRM (`-gk 1`, `-gk 2`) is identical to plain GEMMA's for every model
   (also with `-minor`)
2. with `-model add`, `-lm`/`-lmm` output files are identical to plain GEMMA's
   apart from the last column, with and without covariates
3. `rec` and `dom`, with and without `-minor`, give exactly the same statistics,
   alleles and SNP set as plain GEMMA run on 0/1 genotypes coded in R, and
   `af` and `freq_nonref` agree with R
4. `-hom-maf` keeps exactly the expected SNPs

The code changes are listed in `changesToPopgen.txt`. The original GEMMA README
is in `README_originalGEMMAREADME.md`.
