rd <- function(f) read.table(file.path("output", f), header=TRUE, stringsAsFactors=FALSE, colClasses=c(allele1="character", allele0="character"))
stats <- function(A, B) {
  cols <- setdiff(intersect(names(A), names(B)), c("chr","rs","ps","allele1","allele0","n_miss","n_mis","n_obs","af","freq_HO","freq_carrier"))
  max(sapply(cols, function(cn) max(abs(A[[cn]] - B[[cn]]))))
}
mn <- read.table("minor_allele.txt", header=TRUE, stringsAsFactors=FALSE)
for (m in c("rec","dom")) for (a in c("lm","lmm")) {
  S <- rd(sprintf("swap_%s_%s.assoc.txt", m, a)); O <- rd(sprintf("chr1_%s_%s.assoc.txt", m, a)); B <- rd(sprintf("bb_%s_%s.assoc.txt", m, a))
  fc <- if (m == "rec") "freq_HO" else "freq_carrier"
  cat(sprintf("%s %-3s | swapped vs R ref: SNPs %d/%d same=%s maxdiff=%.3g alleles=%s freq maxdiff=%.4f | swapped vs original bed: identical=%s | allele1 is minor: %s\n",
    m, a, nrow(S), nrow(B), identical(S$rs, B$rs), stats(S, B),
    identical(S$allele1, B$allele1) && identical(S$allele0, B$allele0), max(abs(S[[fc]] - 2*B$af)),
    isTRUE(all.equal(S, O)), all(S$allele1 == mn$minor[match(S$rs, mn$rs)])))
}
for (m in c("rec","dom")) {
  k <- function(f) as.matrix(read.table(file.path("output", f)))
  cat(m, "kinship: swapped vs R ref", max(abs(k(sprintf("swap_%s_kin.cXX.txt",m)) - k(sprintf("bb_%s_kin.cXX.txt",m)))),
      "| swapped vs original", max(abs(k(sprintf("swap_%s_kin.cXX.txt",m)) - k(sprintf("chr1_%s_kin.cXX.txt",m)))), "\n")
}
h <- rd("swap_rec_hommaf.assoc.txt"); cat(sprintf("-minor + -hom-maf 0.05: %d SNPs, freq_HO range %.3f-%.3f\n", nrow(h), min(h$freq_HO), max(h$freq_HO)))
