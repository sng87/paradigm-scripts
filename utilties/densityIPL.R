#!/usr/bin/env Rscript
args = commandArgs(TRUE)
paradigmFile = args[1]

pData = read.table(paradigmFile, sep = "\t", header = T)
X = c()
snRatio = c()
for (i in 1:(dim(pData)[1])) {
    X = c(X, as.numeric(pData[i, 2:(dim(pData)[2])]))
    #snRatio = c(snRatio, sd(as.numeric(pData[i, 2:(dim(pData)[2])]))/abs(mean(as.numeric(pData[i, 2:(dim(pData)[2])]))))
}

pdf("IPL-density.pdf")
plot(density(abs(X)), xlim = c(0,1))
dev.off()

#plot(density(snRatio[!is.na(snRatio)]))