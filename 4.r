dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$BETA <- log(dat$OR)
write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
q()
