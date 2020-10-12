example <- read.csv("example.csv")
snp.mat <- example[,1:20]
phes <- example[,21:22]
source("Multi-CMDR.R")
MRMDR(phes, snp.mat, K=2, cv=10, nperm=0, sele.type='cvc', test.type="mvr", covrt=NULL, trim=T)
