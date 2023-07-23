#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmnet)
library(nortest)
library(matrixStats)

#data <- read.table("CPMMatrixForDNAMethylationFiltered_id.txt")
data <- read.table("./result-confounder/r1/res.txt",header=F)
rownames(data) <- data[,1]
data <- data[,-1]
data <- as.matrix(data)
data <- (data - rowSums(data))/rowSds(data)
data <- as.data.frame(data)

res <- read.table("./result-TF-2/r1/res.txt",header=F)
rownames(res) <- res[,1]
res <- res[,-1]
#res <- as.matrix(res)
#res <- (res - rowSums(res))/rowSds(res)
#res <- as.data.frame(res)

#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
#lncrna_id <- read.table("lncRNA_id.txt",header=FALSE)[,1]
#pcrna_id <- read.table("ProteinCodingGeneID.txt",header=FALSE)[,1]
#clinic = read.table("Clinic.txt",,header=F)

#lncrna <- t(data[(rownames(data) %in% lncrna_id),])
#n <- dim(lncrna)[1]

#mrna <- data[(rownames(data) %in% pcrna_id),grep("case",colnames(data))]
#mrna <- res[(rownames(res) %in% pcrna_id),]

#lncrna <- as.matrix(cbind(lncrna,clinic))
#colnames(lncrna) <- c(colnames(lncrna)[1:1079],"age","gender","american indian or alaska native","asian","black or african american","white")


#var_list <- read.table("./data/in.50k.txt",header=F)
var_list <- read.table("./data/full.model.txt",header=F)
mrna <- res[(rownames(res) %in% var_list[,1]),]


result <- c("id","cor_1","cor_2","cor_3","cor_4","cor_5","cor_ave","cor_all")
ols <- c("id","R2")
for ( i in args[1]:args[2]){
	re <- 0
	y <- t(mrna[i,])
	
	lncrna_id <- as.character(var_list[var_list[,1]==as.character(colnames(y)),2])	
	lncrna <- t(data[(rownames(data) %in% lncrna_id),,drop=FALSE])
	nvar <- dim(lncrna)[2]
	
	rows <- which(!is.na(y))
        #pred <- cbind(colnames(y),length(rows)-floor(0.8*length(rows)))
        row <- c(colnames(y))
        xx <- args[4]

	if (FALSE){
		all.fit <- cv.glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
		c<-coef(all.fit,s='lambda.min',exact=TRUE)
        	inds<-which(c!=0)
        	variables<-row.names(c)[inds]
		lncrna <- t(data[(rownames(data) %in% variables),])
	}
	


	#lncrna_ols <- t(data[(rownames(data) %in% variables),grep("case",colnames(data))])
        #lncrna_ols <- cbind(rep(1,nrow(lncrna)),lncrna)
	fit <- lm(y[,1] ~ lncrna)
        r2 <- summary(fit)$r.squared
	rss <- sum(resid(fit)^2)
	y_hat <- as.numeric(data[colnames(y),])
	tss2 <- sum((y_hat-mean(y_hat))^2)
	tss1 <- sum((y-mean(y))^2)

	rowols <- c(colnames(y),r2,(tss1-rss)/tss2,nvar)
	#rowols <- cbind(rowols,r2)
	write.table(rbind(rowols),args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	#ols <- rbind(ols,rowols)

	#x <- summary(fit)$coefficients[-1,]
	var <- gsub("lncrna_ols","",rownames(summary(fit)$coefficients)[-1])
	p <- summary(fit)$coefficients[-1,4]	

	for (n in 1:length(var)){
		tmp <- c(colnames(y),var[n],p[n])
		write.table(rbind(tmp),args[5],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	}	

	rowres <- c(colnames(y),residuals(fit))
	write.table(rbind(rowres),args[7],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

}

#write.table(result,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#write.table(ols,args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
