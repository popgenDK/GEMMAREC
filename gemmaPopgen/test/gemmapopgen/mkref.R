# R reference: 0/1 rec/dom coding (A1 as in .bim, and flipped to the minor allele) as BIMBAM
raw <- read.table("rnd.raw", header=TRUE, check.names=FALSE)
bim <- read.table("rnd.bim", stringsAsFactors=FALSE); fam <- read.table("rnd.fam")
G <- as.matrix(raw[, -(1:6)]); keep <- !is.na(fam$V6)
p <- colMeans(G[keep, ], na.rm=TRUE) / 2
flip <- p > 0.5
for (pol in c("a1", "minor")) {
  X <- G; a1 <- bim$V5; a0 <- bim$V6
  if (pol == "minor") { X[, flip] <- 2 - X[, flip]; a1 <- ifelse(flip, bim$V6, bim$V5); a0 <- ifelse(flip, bim$V5, bim$V6) }
  for (m in c("rec", "dom")) {
    C <- if (m == "rec") (X == 2) * 1 else (X >= 1) * 1
    write.table(data.frame(bim$V2, a1, a0, t(C)), sprintf("%s_%s.geno.txt", m, pol), sep=",",
                quote=FALSE, row.names=FALSE, col.names=FALSE, na="NA")
  }
}
write.table(fam$V6, "pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, na="NA")
