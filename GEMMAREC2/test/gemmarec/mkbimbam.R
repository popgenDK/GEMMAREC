# Write 0/1 recessive and dominant codings of chr1 as BIMBAM mean-genotype files
raw <- read.table("chr1.raw", header=TRUE, check.names=FALSE)
bim <- read.table("chr1.bim", stringsAsFactors=FALSE)
fam <- read.table("chr1.fam", stringsAsFactors=FALSE)
G <- as.matrix(raw[, -(1:6)])            # count of A1 (bim col 5)
stopifnot(ncol(G) == nrow(bim))
wr <- function(X, f) {
  X[is.na(X)] <- NA
  d <- data.frame(bim$V2, bim$V5, bim$V6, t(X))
  write.table(d, f, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, na="NA")
}
wr((G == 2) * 1, "rec.geno.txt")
wr((G >= 1) * 1, "dom.geno.txt")
write.table(fam[, 6:11], "pheno.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, na="NA")
# independent HO / carrier frequencies
fr <- data.frame(rs=bim$V2, ho=colMeans(G == 2, na.rm=TRUE), carrier=colMeans(G >= 1, na.rm=TRUE))
write.table(fr, "freq_R.txt", quote=FALSE, row.names=FALSE)
