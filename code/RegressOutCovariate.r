#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
output <- args[1]

library(glmnet)

data <- read.table("../data/CPMMatrix.txt")
clinic <- read.delim("../data/Clinic_dummy.txt",header=T,sep="\t")
clinic <- as.matrix(clinic[,1:8])

mrna <- data[,grep("case",colnames(data))]

for ( i in 1:dim(mrna)[[1]]){
	re <- 0
	y <- t(mrna[i,])
	fit <- lm(y[,1] ~ clinic)
	rowres <- c(colnames(y),residuals(fit))
        write.table(rbind(rowres),output,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
}
