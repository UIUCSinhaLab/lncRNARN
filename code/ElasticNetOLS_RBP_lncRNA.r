#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmnet)
library(MASS)
library(matrixStats)
library(nortest)


#data <- read.table("CPMMatrixForDNAMethylationFiltered_id.txt")
data <- read.table("./result-confounder/r1/res.txt",header=F)
rownames(data) <- data[,1]
data <- data[,-1]
data <- as.matrix(data)
data <- (data - rowSums(data))/rowSds(data)
data <- as.data.frame(data)


#lncrna_id <- read.table("lncRNA.txt",header=FALSE)[,1]
lncrna_id <- read.table("lncRNA_id.txt",header=FALSE)[,1]
rbp_id <- read.table("RBP.1.txt",header=FALSE)[,1]
pcrna_id <- read.table("ProteinCodingGeneID.txt",header=FALSE)[,1]

lncrna <- t(data[(rownames(data) %in% lncrna_id),])
#n <- dim(lncrna)[1]

rbp <- t(data[(rownames(data) %in% rbp_id),])

res <- read.table("./result-TF-2/r1/res.txt",header=F)
rownames(res) <- res[,1]
res <- res[,-1]
#mrna <- data[(rownames(data) %in% pcrna_id),]
mrna <- res[(rownames(res) %in% pcrna_id),]

lncvar <- read.table("./result-TF-2/r3/var.txt",header=F)
rbpvar <- read.table("./result-TF-2/rbp/var.txt",header=F)
#lines <- c()
#for ( i in 1:dim(mrna)[1] ){
#        p <- lillie.test(as.numeric(mrna[i,]))$p
#        if ( p > 0.0000000001 ){
#                lines <- c(lines,i)
#        }
#}
#mrna <- mrna[lines,]

for ( i in args[1]:args[2]){
	re <- 0
	print(paste0("i=",i))
	y <- t(mrna[i,])
	
	rows <- which(!is.na(y[,1]))

	#all.l.sequence <- glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,alpha=0.5,family="gaussian",dfmax=5)$lambda
	##all.fit <- cv.glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
	##c<-coef(all.fit,s='lambda.min',exact=TRUE)
        ##inds<-which(c!=0)
        #lncrnas<-row.names(c)[inds]
	lncrnas <- as.character(lncvar[lncvar[,1]==as.character(colnames(y)),2])
	##lncrnas <- row.names(c)[inds][grep("ENSG",row.names(c)[inds])]
	
	#rbps <- colnames(rbp)[!(colnames(rbp)==colnames(y))]
	rbps <- as.character(rbpvar[rbpvar[,1]==as.character(colnames(y)),2])

	variables <- c(lncrnas,rbps)

	if( length(variables) > 0){
	lncrna_ols <- t(data[(rownames(data) %in% lncrnas),])
        lncrna_ols <- cbind(rep(1,nrow(lncrna_ols)),lncrna_ols)
	lnc.fit <- lm(y[,1] ~ lncrna_ols)
        lnc.r2 <- summary(lnc.fit)$r.squared

	rbp_ols <- t(data[(rownames(data) %in% rbps),])
        rbp_ols <- cbind(rep(1,nrow(rbp_ols)),rbp_ols)
        rbp.fit <- lm(y[,1] ~ rbp_ols)
        rbp.r2 <- summary(rbp.fit)$r.squared

	#data_tmp <- cbind(y,lncrna_ols,rbp)	

	both_ols <- t(data[(rownames(data) %in% variables),])
        both_ols <- cbind(rep(1,nrow(both_ols)),both_ols)
        both.fit <- lm(y[,1] ~ both_ols)
        both.r2 <- summary(both.fit)$r.squared


	data_tmp <- cbind(y,both_ols)

	f_0 = paste(colnames(y),"~")
	#for (e_lncrna in lncrnas){
	#	f_0 <- paste(f_0,"+",e_lncrna)
	#}
	for (e_var in variables){
                f_0 <- paste(f_0,"+",e_var)
        }
	f_1 <- f_0
	f_0 <- as.formula(f_0)
	null <- lm(f_0,data = as.data.frame(data_tmp))
	
	for (e_lncrna in lncrnas){
                for (e_rbp in rbps){
			e_int <- paste0(e_lncrna,":",e_rbp)
			f_1 <- paste(f_1,"+",e_int)
		}
	}
	f_1 <- as.formula(f_1)	
	#glm(f_1,data = as.data.frame(data_tmp))
	#fit <- stepAIC(null,direction="forward",scope=list(lower=f_0, upper=f_1),trace=F,steps=5)
	fit <- stepAIC(null,direction="forward",scope=list(lower=f_0, upper=f_1),trace=F)
	max.r2 <- summary(fit)$r.square

	#row <- c(colnames(y),lnc.r2,max.r2)
        row <- c(colnames(y),lnc.r2,rbp.r2,both.r2,max.r2)
	write.table(rbind(row),args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

	
	if (FALSE){
	x <- summary(fit)$coefficients[-1,4] < 0.05
	relevant.x <- names(x)[x == TRUE]

	var <- relevant.x[grep(":",relevant.x)]
	p <- summary(fit)$coefficients[-1,4][x==TRUE][grep(":",relevant.x)]

	for (n in 1:length(var)){
		tmp <- c(colnames(y),var[n],as.character(p[n]))
		write.table(rbind(tmp),args[4],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	}
	}

	p <- summary(fit)$coefficients[-1,4]
	var <- names(p)
	#var <- relevant.x[grep("ENSG",relevant.x)]
	for (n in 1:length(var)){
                #tmp <- c(colnames(y),var[n],as.character(p[n]))
		tmp <- c(colnames(y),var[n],as.character(p[n]),length(lncrnas),length(rbps),as.numeric(p[n])*length(lncrnas)*length(rbps))
                write.table(rbind(tmp),args[4],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
        }	
	}
}

#write.table(result,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#write.table(ols,args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
