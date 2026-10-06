rd <- function(f) read.table(file.path("output", f), header=TRUE, stringsAsFactors=FALSE)
for (m in c("rec","dom")) for (a in c("lm","lmm","mvlmm")) {
  A <- rd(sprintf("bed_%s_%s.assoc.txt", m, a)); B <- rd(sprintf("bb_%s_%s.assoc.txt", m, a))
  fc <- intersect(names(A), c("freq_HO","freq_carrier"))
  num <- setdiff(intersect(names(A), names(B)), c("chr","rs","ps","allele1","allele0","n_miss","n_mis","n_obs"))
  rel <- sapply(num, function(cn) max(abs(A[[cn]]-B[[cn]]) / pmax(abs(B[[cn]]), 1e-12)))
  cat(sprintf("%s %-5s: SNPs bed=%d bimbam=%d same=%s | freq col=%s, max|freq-2*af_bimbam|=%.4f | max rel diff stats=%.2g (%s)\n",
      m, a, nrow(A), nrow(B), identical(A$rs, B$rs), fc, max(abs(A[[fc]] - 2*B$af)),
      max(rel), names(which.max(rel))))
}
# independent frequencies on individuals used by -lm (phenotype 1 non-missing)
raw <- read.table("chr1.raw", header=TRUE, check.names=FALSE); fam <- read.table("chr1.fam")
G <- as.matrix(raw[!is.na(fam$V6), -(1:6)])
ho <- colMeans(G == 2, na.rm=TRUE); ca <- colMeans(G >= 1, na.rm=TRUE)
r <- rd("bed_rec_lm.assoc.txt"); d <- rd("bed_dom_lm.assoc.txt")
idx <- match(r$rs, sub("_.*", "", colnames(G))); idd <- match(d$rs, sub("_.*", "", colnames(G)))
cat(sprintf("freq_HO vs R: max diff %.4f; freq_carrier vs R: max diff %.4f\n",
    max(abs(r$freq_HO - ho[idx])), max(abs(d$freq_carrier - ca[idd]))))
# -hom-maf 0.05 with -maf 0
miss <- colMeans(is.na(G))
for (m in c("rec","dom")) {
  f <- if (m == "rec") ho else ca
  exp <- sub("_.*", "", colnames(G))[f >= 0.05 & f <= 0.95 & miss <= 0.05]
  got <- rd(sprintf("bed_%s_hommaf.assoc.txt", m))$rs
  fcol <- rd(sprintf("bed_%s_hommaf.assoc.txt", m))[[if (m=="rec") "freq_HO" else "freq_carrier"]]
  cat(sprintf("hom-maf %s: expected %d, got %d, identical=%s, freq range %.3f-%.3f\n",
      m, length(exp), length(got), identical(sort(exp), sort(got)), min(fcol), max(fcol)))
}
