# GEMMAREC2: GEMMA 0.98.6 with recessive and dominant association models

This is a modified version of [GEMMA](https://github.com/genetics-statistics/GEMMA)
(master, version 0.98.6). It updates the original GEMMAREC (in the parent
directory of this repository), which was based on GEMMA 0.98.5.

For PLINK input (`-bfile`), GEMMA recodes the genotypes as it reads them, so
association tests and the kinship matrix use the chosen genetic model instead
of the additive one. By default the tested allele is A1, the first allele in the
`.bim` file (`allele1` in the output).

| model | option | A1/A1 | A1/A2 | A2/A2 |
|---|---|---|---|---|
| recessive | (default) | 1 | 0 | 0 |
| dominant | `-dom` | 1 | 1 | 0 |

BIMBAM input (`-g`) is never recoded.

> **Only tested for single-trait GWAS.** The recessive and dominant models
> have only been validated for single-trait association: `-lm`, univariate
> `-lmm`, and the kinship matrix from `-gk`. Other analyses also use the
> recoded genotypes, including the multivariate LMM, `-bslmm` and `-calccor`,
> but they have not been tested and should be used with caution.

## Usage

Build with `make`. GSL and OpenBLAS are required; see `INSTALL.md`. The
binary is `bin/gemmarec2`. Use it like GEMMA, with three extra options:

| option | meaning |
|---|---|
| `-dom` | dominant model (default is recessive) |
| `-minor` | test the minor allele instead of A1 (see below) |
| `-hom-maf [num]` | remove SNPs where `freq_HO` (or `freq_carrier` with `-dom`) is below `num` or above `1-num` |

```sh
# kinship matrix and single-trait LMM, recessive model, minor allele tested
bin/gemmarec2 -bfile data -minor -gk 1 -o kin
bin/gemmarec2 -bfile data -minor -k output/kin.cXX.txt -lmm 4 -maf 0 -hom-maf 0.01 -o rec

# the same with the dominant model
bin/gemmarec2 -bfile data -minor -dom -gk 1 -o kin_dom
bin/gemmarec2 -bfile data -minor -dom -k output/kin_dom.cXX.txt -lmm 4 -o dom

# linear model without kinship
bin/gemmarec2 -bfile data -minor -lm 4 -o rec_lm
```

Use the same `-dom`/`-minor` options for `-gk` as for the association run,
so that the kinship matrix uses the same genotype coding.

The output has the usual GEMMA columns, but `af` is replaced by
* `freq_HO`: frequency of A1/A1 homozygotes (recessive model)
* `freq_carrier`: frequency of A1/A1 plus A1/A2 (dominant model)

Frequencies are computed over the analysed individuals with non-missing
genotypes. The model and tested allele are written to the `.log.txt` file.

### Which allele is tested (`-minor`)

GEMMA does not polarise alleles. It always tests A1 from the `.bim`. PLINK 1.9
usually makes A1 the minor allele, but files made with `--keep-allele-order`,
or with PLINK 2, may have the major allele as A1. In that case the recessive
model would test homozygotes for the major allele.

With `-minor`, A1 and A2 are swapped for every SNP where A1 has a frequency
above 0.5 in the analysed individuals. `allele1` in the output is then always
the tested minor allele. Without `-minor`, A1 is tested as before.

## Changes compared with the old GEMMAREC

* **Frequency column.** The old `af` column showed half the frequency of
  A1/A1 homozygotes and should not be used. It is now `freq_HO` or
  `freq_carrier` (see above).
* **New options:** `-dom`, `-minor` and `-hom-maf`.
* **`-maf` is unchanged.** It still filters on the internal value, which is
  `freq_HO/2` (or `freq_carrier/2`), exactly as in the old GEMMAREC. With the
  default `-maf 0.01`, SNPs with `freq_HO < 0.02` are removed. Use
  `-maf 0 -hom-maf X` to filter only on the genotype frequency.
* **Multivariate LMM (`-lmm` with `-n` and several phenotypes) and
  `-calccor` now use the recoded genotypes too.** The old GEMMAREC did not
  recode genotypes in these code paths. Its multivariate LMM results on PLINK
  input used the additive model, even though the frequency column was
  computed from the recessive coding.

## Reproducing old GEMMAREC results

Without `-dom` and `-minor`, the new version gives results identical to the
old GEMMAREC for `-gk`, `-lm` and univariate `-lmm`. That includes the same
SNPs, `beta`, `se`, likelihoods and p-values (checked on
`example/mouse_hs1940`). The only difference is the frequency column, where
`freq_HO = 2 * af_old`.

## Tests

`test/gemmarec/run_all.sh` takes about a minute and needs `plink` 1.9 and R.
Run it from this directory. It uses chromosome 1 of `example/mouse_hs1940`.

* **Recessive and dominant models:** the script runs `-gk`, `-lm`, `-lmm` and
  the multivariate LMM with both models. For each analysis it compares two
  inputs:
  * the PLINK file, which GEMMA recodes when it reads it
  * BIMBAM files with the same 0/1 coding made in R, which GEMMA does not recode
* **`-minor`:** the script swaps A1 and A2 for every SNP, runs `-minor` on the
  swapped file, and compares it with an R reference (flipped to the minor allele
  and then recoded) and with the original file.

All association statistics, tested alleles and kinship matrices match exactly.
The script also checks the frequencies and `-hom-maf` against R. For the
multivariate LMM, this only shows that GEMMA reads the genotypes with the right
coding. It does not validate multivariate results.

The code changes are listed in `changesToRecessive.txt`. The original GEMMA
README is in `README_originalGEMMAREADME.md`.
