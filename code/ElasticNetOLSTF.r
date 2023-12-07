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

#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
#lncrna_id <- read.table("lncRNA_id.txt",header=FALSE)[,1]
lncrna_id <- read.table("HumanTFList.txt",header=FALSE)[,1]
pcrna_id <- read.table("ProteinCodingGeneID.txt",header=FALSE)[,1]

lncrna <- t(data[(rownames(data) %in% lncrna_id),])

mrna <- data[(rownames(data) %in% pcrna_id),]

lines <- c()
for ( i in 1:13963 ){
	p <- lillie.test(as.numeric(mrna[i,]))$p
	if ( p > 0.0000000001 ){
		lines <- c(lines,i)
	}
}
mrna <- mrna[lines,]

#lncrna <- as.matrix(cbind(lncrna,clinic))
#colnames(lncrna) <- c(colnames(lncrna)[1:1255],"age","gender","american indian or alaska native","asian","black or african american","white")
#mrna <- t(cbind(t(mrna),clinic))


result <- c("id","cor_1","cor_2","cor_3","cor_4","cor_5","cor_ave","cor_all")
ols <- c("id","R2")
for ( i in args[1]:args[2]){
	re <- 0
	y <- t(mrna[i,])

	
	rows <- which(!is.na(y[,1]))
	pred <- cbind(colnames(y),length(rows)-floor(0.8*length(rows)))
	row <- c(colnames(y))

	for ( j in 1:5){
		train_ids <- sample(1:length(rows),floor(0.8*length(rows)))
		train_rows <- rows[train_ids]
		test_rows <- rows[-train_ids]
		x.train <- lncrna[train_rows,!(colnames(lncrna) %in% colnames(y))]
		x.text <- lncrna[test_rows,!(colnames(lncrna) %in% colnames(y))]
		y.train <- y[train_rows,]
		y.text <- y[test_rows,]
		#l.sequence <- glmnet(x.train,y.train,alpha=0.5,family="gaussian",dfmax=5)$lambda
		#alpha0.5.fit <- glmnet(x.train,y.train,alpha=0.5,family="gaussian",dfmax=5)
		alpha0.5.fit <- cv.glmnet(x.train,y.train,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
		alpha0.5.predicted <- predict(alpha0.5.fit,s=alpha0.5.fit$lambda.min,newx=x.text)
		pred <- c(pred,alpha0.5.predicted)
		cc <- cor(y.text, alpha0.5.predicted)[1]
		print(j)
		print(cc)
		row <- cbind(row,cc)
		re <- re + cc
	}
	#write(pred,args[4],ncolumns=2+5*length(test_rows),append=TRUE)
	print("all")
	print(re)
	row <- cbind(row,re/5)

	#all.l.sequence <- glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,alpha=0.5,family="gaussian",dfmax=5)$lambda
	#print(all.l.sequence)
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
	tss <- sum((y-mean(y))^2)

	rowols <- c(colnames(y),r2,all.fit$lambda.min,(1-rss/tss))
	#rowols <- cbind(rowols,r2)
	write.table(rbind(rowols),args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	#ols <- rbind(ols,rowols)

	rowres <- c(colnames(y),residuals(fit))
	write.table(rbind(rowres),args[7],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

	#c <- c[inds]
	#var <- c(colnames(y))
	#for (n in 2:length(c)){
	#	var <- c(var,c[n],variables[n])
	#}

	p <- summary(fit)$coefficients[-1,4]
	var <- names(p)
	var <- gsub("lncrna_ols","",var)
	for (n in 1:length(var)){
                tmp <- c(colnames(y),var[n],as.character(p[n]))
                write.table(rbind(tmp),args[5],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
        }	

	#write(var,args[5],ncolumns=length(var),append=TRUE)
}

#write.table(result,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#write.table(ols,args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
