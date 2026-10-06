# R reference for -minor: flip to minor allele (among individuals with phenotype 1), then 0/1 coding
raw <- read.table("swap.raw", header=TRUE, check.names=FALSE)
bim <- read.table("swap.bim", stringsAsFactors=FALSE); fam <- read.table("chr1.fam")
G <- as.matrix(raw[, -(1:6)])
f <- colMeans(G[!is.na(fam$V6), ], na.rm=TRUE) / 2
flip <- f > 0.5
G[, flip] <- 2 - G[, flip]
a1 <- ifelse(flip, bim$V6, bim$V5); a0 <- ifelse(flip, bim$V5, bim$V6)
wr <- function(X, fn) write.table(data.frame(bim$V2, a1, a0, t(X)), fn, sep=",",
  quote=FALSE, row.names=FALSE, col.names=FALSE, na="NA")
wr((G == 2) * 1, "rec.geno.txt"); wr((G >= 1) * 1, "dom.geno.txt")
write.table(data.frame(rs=bim$V2, minor=a1), "minor_allele.txt", quote=FALSE, row.names=FALSE)
