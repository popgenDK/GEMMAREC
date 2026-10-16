# Give about half of the SNPs a random allele orientation (A1 = major allele)
bim <- read.table("base.bim", stringsAsFactors=FALSE)
set.seed(1); sw <- runif(nrow(bim)) < 0.5
write.table(data.frame(bim$V2, ifelse(sw, bim$V6, bim$V5)), "a1_rnd.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
