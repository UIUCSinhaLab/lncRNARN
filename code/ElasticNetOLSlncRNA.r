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
#res <- read.table("./result-TF-2/in1m/res.txt",header=F)
rownames(res) <- res[,1]
res <- res[,-1]
#res <- as.matrix(res)
#res <- (res - rowSums(res))/rowSds(res)
#res <- as.data.frame(res)

#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
lncrna_id <- read.table("lncRNA_id.txt",header=FALSE)[,1]
#lncrna_id <- read.table("RBP.1.txt",header=FALSE)[,1]
pcrna_id <- read.table("ProteinCodingGeneID.txt",header=FALSE)[,1]
#clinic = read.table("Clinic.txt",,header=F)

lncrna <- t(data[(rownames(data) %in% lncrna_id),])
#n <- dim(lncrna)[1]

#mrna <- data[(rownames(data) %in% pcrna_id),grep("case",colnames(data))]
mrna <- res[(rownames(res) %in% pcrna_id),]
mrna <- mrna[!(rownames(res) %in% lncrna_id),]

#lncrna <- as.matrix(cbind(lncrna,clinic))
#colnames(lncrna) <- c(colnames(lncrna)[1:1079],"age","gender","american indian or alaska native","asian","black or african american","white")

result <- c("id","cor_1","cor_2","cor_3","cor_4","cor_5","cor_ave","cor_all")
ols <- c("id","R2")
for ( i in args[1]:args[2]){
	re <- 0
	y <- t(mrna[i,])
	
	rows <- which(!is.na(y))
	pred <- cbind(colnames(y),length(rows)-floor(0.8*length(rows)))
	row <- c(colnames(y))
	if (FALSE){
	for ( j in 1:5){
		train_ids <- sample(1:length(rows),floor(0.8*length(rows)))
		train_rows <- rows[train_ids]
		test_rows <- rows[-train_ids]
		x.train <- lncrna[train_rows,!(colnames(lncrna) %in% colnames(y))]
		x.text <- lncrna[test_rows,!(colnames(lncrna) %in% colnames(y))]
		y.train <- y[train_rows,]
		y.text <- y[test_rows,]
		#l.sequence <- glmnet(x.train,y.train,alpha=0.5,family="gaussian",dfmax=5)$lambda
		alpha0.5.fit <- cv.glmnet(x.train,y.train,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
		alpha0.5.predicted <- predict(alpha0.5.fit,s=alpha0.5.fit$lambda.min,newx=x.text)
		pred <- c(pred,alpha0.5.predicted)
		cc <- cor(y.text, alpha0.5.predicted)[1]
		print(j)
		print(cc)
		row <- cbind(row,cc)
		re <- re + cc
	}
	}
	xx <- args[4]
	#write(pred,args[4],ncolumns=2+5*length(test_rows),append=TRUE)
	print("all")
	#print(re)
	#row <- cbind(row,re/5)

	#all.l.sequence <- glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,alpha=0.5,family="gaussian",dfmax=5)$lambda
	all.fit <- cv.glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
	all.pred <- predict(all.fit,s=all.fit$lambda.min,newx=lncrna[,!(colnames(lncrna) %in% colnames(y))])
	all.cc <- cor(y,all.pred)
	row <- cbind(row,all.cc)
	write.table(rbind(row),args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	#result <- rbind(result,row)
	
	c<-coef(all.fit,s='lambda.min',exact=TRUE)
	inds<-which(c!=0)
	variables<-row.names(c)[inds]

	#lncrna_ols <- t(data[(rownames(data) %in% variables),grep("case",colnames(data))])
	lncrna_ols <- lncrna[,(colnames(lncrna) %in% variables)]
        lncrna_ols <- cbind(rep(1,nrow(lncrna_ols)),lncrna_ols)
	fit <- lm(y[,1] ~ lncrna_ols)
        r2 <- summary(fit)$r.squared
	rss <- sum(resid(fit)^2)
	y_hat <- as.numeric(data[colnames(y),])
	tss2 <- sum((y_hat-mean(y_hat))^2)
	tss1 <- sum((y-mean(y))^2)
	rowols <- c(colnames(y),r2,(tss1-rss)/tss2)
	#rowols <- cbind(rowols,r2)
	write.table(rbind(rowols),args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	#ols <- rbind(ols,rowols)

	#x <- summary(fit)$coefficients[-1,]
	#var <- gsub("lncrna_ols","",rownames(x))
	p <- summary(fit)$coefficients[-1,4]	
	var <- names(p)
        var <- gsub("lncrna_ols","",var)

	if(length(var)>0){
	for (n in 1:length(var)){
		tmp <- c(colnames(y),var[n],p[n])
		write.table(rbind(tmp),args[5],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	}
	}

	rowres <- c(colnames(y),residuals(fit))
	write.table(rbind(rowres),args[7],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

}

#write.table(result,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#write.table(ols,args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
