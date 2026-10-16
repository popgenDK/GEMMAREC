# Compare gemmaPopgen with plain GEMMA (GRM, additive) and with R-coded BIMBAM references (rec/dom)
rd <- function(f) read.table(file.path("output", f), header=TRUE, stringsAsFactors=FALSE,
                             colClasses=c(allele1="character", allele0="character"))
K <- function(f) as.matrix(read.table(file.path("output", f)))
cat("1. GRM vs plain GEMMA (max|diff|):\n")
for (m in c("add","rec","dom")) cat(sprintf("   -model %s: -gk 1 %g, -gk 2 (with -minor) %g\n", m,
    max(abs(K(sprintf("p_gk1_%s.cXX.txt", m)) - K("u_gk1.cXX.txt"))), max(abs(K(sprintf("p_gk2_%s.sXX.txt", m)) - K("u_gk2.sXX.txt")))))
cat("\n2. additive vs plain GEMMA (assoc files byte-identical):\n")
for (f in c("lm_nocov","lmm_nocov","lm_cov","lmm_cov")) cat(sprintf("   %-10s %s\n", f,
    identical(readLines(sprintf("output/u_%s.assoc.txt", f)), readLines(sprintf("output/p_%s.assoc.txt", f)))))
cat("\n3. rec/dom vs R-coded BIMBAM in plain GEMMA:\n")
raw <- read.table("rnd.raw", header=TRUE, check.names=FALSE); fam <- read.table("rnd.fam")
G <- as.matrix(raw[!is.na(fam$V6), -(1:6)]); rs <- sub("_.*", "", colnames(G))
p <- colMeans(G, na.rm=TRUE)/2; flip <- p > 0.5
gemma_set <- function(cov) rd(sprintf("u_lm_%s.assoc.txt", cov))$rs   # SNPs passing GEMMA's own filters
for (m in c("rec","dom")) for (pol in c("a1","minor")) for (a in c("lm","lmm","cov")) {
  if (a == "cov" && pol == "a1") next
  pf <- if (a == "cov") sprintf("p_%s_minor_cov", m) else sprintf("p_%s_%s_%s", m, pol, a)
  rf <- if (a == "cov") sprintf("r_%s_minor_cov", m) else sprintf("r_%s_%s_%s", m, pol, a)
  P <- rd(paste0(pf, ".assoc.txt")); R <- rd(paste0(rf, ".assoc.txt"))
  X <- G; if (pol == "minor") X[, flip] <- 2 - X[, flip]
  C <- if (m == "rec") X == 2 else X >= 1
  fC <- colMeans(C, na.rm=TRUE)
  exp <- intersect(gemma_set(if (a == "cov") "cov" else "nocov"), rs[fC > 0 & fC < 1])
  R <- R[match(P$rs, R$rs), ]
  cols <- setdiff(intersect(names(P), names(R)), c("chr","rs","ps","allele1","allele0","n_miss","n_mis","n_obs","af","freq_HO","freq_carrier"))
  fcol <- if (m == "rec") "freq_HO" else "freq_carrier"
  cat(sprintf("   %s %-5s %-3s: SNPs %d, expected set %s | stats max|diff| %.2g | alleles %s | %s vs R %.4f\n",
      m, pol, a, nrow(P), identical(sort(P$rs), sort(exp)), max(sapply(cols, function(cn) max(abs(P[[cn]] - R[[cn]])))),
      identical(P$allele1, R$allele1), fcol, max(abs(P[[fcol]] - fC[match(P$rs, rs)]))))
}
cat("\n4. -hom-maf 0.1 (with -minor):\n")
for (m in c("rec","dom")) {
  P <- rd(sprintf("p_%s_hommaf.assoc.txt", m)); X <- G; X[, flip] <- 2 - X[, flip]
  fC <- colMeans(if (m == "rec") X == 2 else X >= 1, na.rm=TRUE)
  exp <- intersect(gemma_set("nocov"), rs[fC >= 0.1 & fC <= 0.9])
  cat(sprintf("   %s: %d SNPs, expected set %s\n", m, nrow(P), identical(sort(P$rs), sort(exp))))
}
