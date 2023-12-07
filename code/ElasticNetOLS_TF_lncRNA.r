#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmnet)
library(MASS)
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
lncrna_id <- read.table("lncRNA_id.txt",header=FALSE)[,1]
tf_id <- read.table("HumanTFList.txt",header=FALSE)[,1]
pcrna_id <- read.table("ProteinCodingGeneID.txt",header=FALSE)[,1]
#clinic = read.table("Clinic.txt",header=F)

lncrna <- t(data[(rownames(data) %in% lncrna_id),])
#lncrna <- as.matrix(cbind(lncrna,clinic))
#colnames(lncrna) <- c(colnames(lncrna)[1:1079],"age","gender","american indian or alaska native","asian","black or african american","white")
#n <- dim(lncrna)[1]
tf <- t(data[(rownames(data) %in% tf_id),])
#tf <- as.matrix(cbind(tf,clinic))
#colnames(tf) <- c(colnames(tf)[1:1255],"age","gender","american indian or alaska native","asian","black or african american","white")
#n <- dim(lncrna)[1]

mrna <- data[(rownames(data) %in% pcrna_id),]
tfvar <- read.table("./result-TF-2/r1/var.txt",header=F)
lncvar <- read.table("./result-TF-2/r3/var.txt",header=F)

lines <- c()
for ( i in 1:dim(mrna)[1] ){
        p <- lillie.test(as.numeric(mrna[i,]))$p
        if ( p > 0.0000000001 ){
                lines <- c(lines,i)
        }
}
mrna <- mrna[lines,]




#lncrna <- cbind(lncrna,clinic)

#tf <- cbind(tf,clinic)

for ( i in args[1]:args[2]){
	re <- 0
	print(paste0("i=",i))
	y <- t(mrna[i,])
	
	rows <- which(!is.na(y[,1]))

	##print("1\n")
	#all.l.sequence <- glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,alpha=0.5,family="gaussian",dfmax=5)$lambda
	##print("2\n")
	##all.fit <- cv.glmnet(lncrna[,!(colnames(lncrna) %in% colnames(y))],y,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
	##print("3\n")
	##all.pred <- predict(all.fit,s=all.fit$lambda.min,newx=lncrna[,!(colnames(lncrna) %in% colnames(y))])
	##print("4\n")
	##all.lncrna.cc <- cor(y,all.pred)
	##c<-coef(all.fit,s='lambda.min',exact=TRUE)
        ##inds<-which(c!=0)
        #lncrnas<-row.names(c)[inds]
	##lncrnas <- row.names(c)[inds][grep("ENSG",row.names(c)[inds])]
	lncrnas <- as.character(lncvar[lncvar[,1]==as.character(colnames(y)),2])
	
	#all.l.sequence <- glmnet(tf[,!(colnames(tf) %in% colnames(y))],y,alpha=0.5,family="gaussian",dfmax=5)$lambda
        ##all.fit <- cv.glmnet(tf[,!(colnames(tf) %in% colnames(y))],y,dfmax=10,type.measure="mse",alpha=0.5,family="gaussian")
        ##all.pred <- predict(all.fit,s=all.fit$lambda.min,newx=tf[,!(colnames(tf) %in% colnames(y))])
        ##all.tf.cc <- cor(y,all.pred)
	##c<-coef(all.fit,s='lambda.min',exact=TRUE)
        ##inds<-which(c!=0)
        #tfs<-row.names(c)[inds]
	##tfs <- row.names(c)[inds][grep("ENSG",row.names(c)[inds])]
	tfs <- as.character(tfvar[tfvar[,1]==as.character(colnames(y)),2])

	variables <- c(lncrnas,tfs)

	#row <- c(colnames(y),all.lncrna.cc,all.tf.cc)
	#write.table(rbind(row),args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	#result <- rbind(result,row)


	# the varianble name is misleading. lncrna_ols is the matrix for all selected variables.
	lncrna_ols <- t(data[(rownames(data) %in% lncrnas),])
        lncrna_ols <- cbind(rep(1,nrow(lncrna_ols)),lncrna_ols)
	lnc.fit <- lm(y[,1] ~ lncrna_ols)
        lnc.r2 <- summary(lnc.fit)$r.squared

	tf_ols <- t(data[(rownames(data) %in% tfs),])
        tf_ols <- cbind(rep(1,nrow(tf_ols)),tf_ols)
        tf.fit <- lm(y[,1] ~ tf_ols)
        tf.r2 <- summary(tf.fit)$r.squared

	both_ols <- t(data[(rownames(data) %in% variables),])
        both_ols <- cbind(rep(1,nrow(both_ols)),both_ols)
        both.fit <- lm(y[,1] ~ both_ols)
        both.r2 <- summary(both.fit)$r.squared

	#row <- c(colnames(y),lnc.r2,tf.r2,both.r2)
	#write.table(rbind(row),args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

	if (F) {
	for (e_lncrna in lncrnas){
		for (e_tf in tfs){
			if (length(both_ols[,colnames(both_ols) %in% e_lncrna]) > 0 && length(both_ols[,colnames(both_ols) %in% e_tf]) > 0){
				int_ols <- both_ols
				x <- both_ols[,colnames(both_ols) %in% e_lncrna] * both_ols[,colnames(both_ols) %in% e_tf]
				#name <- colnames(int_ols)
				int_ols <- cbind(int_ols,x)
				#colnames(int_ols) <- c(name,paste0(lncrna,"_",tf))
				int.fit <- lm(y[,1] ~ int_ols)	
				if (summary(int.fit)$coefficients["int_olsx",4] <0.05){
					tmp <- c(colnames(y),e_lncrna,e_tf,summary(int.fit)$coefficients["int_olsx",4],summary(int.fit)$r.squared)
					write.table(rbind(tmp),args[4],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
				}
			}
		}
	}
	}

	data_tmp <- cbind(y,both_ols)	

	f_0 = paste(colnames(y),"~")
	for (e_var in variables){
		f_0 <- paste(f_0,"+",e_var)
	}
	f_1 <- f_0
	f_0 <- as.formula(f_0)
	null <- lm(f_0,data = as.data.frame(data_tmp))
	
	for (e_lncrna in lncrnas){
                for (e_tf in tfs){
			e_int <- paste0(e_lncrna,":",e_tf)
			f_1 <- paste(f_1,"+",e_int)
		}
	}
	f_1 <- as.formula(f_1)	
	#glm(f_1,data = as.data.frame(data_tmp))
	#fit <- stepAIC(null,direction="forward",scope=list(lower=f_0, upper=f_1),trace=F)
	fit <- stepAIC(null,direction="forward",scope=list(lower=f_0, upper=f_1),trace=F,steps=5)
	max.r2 <- summary(fit)$r.square

	row <- c(colnames(y),lnc.r2,tf.r2,both.r2,max.r2)
        write.table(rbind(row),args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)

	
	if (FALSE){
	x <- summary(fit)$coefficients[-1,4] < 0.05
	relevant.x <- names(x)[x == TRUE]

	var <- relevant.x[grep(":",relevant.x)]
	p <- summary(fit)$coefficients[-1,4][x==TRUE][grep(":",relevant.x)]

	for (n in 1:length(var)){
		tmp <- c(colnames(y),var[n],as.character(p[n]),length(lncrnas),length(tfs),as.numeric(p[n])*length(lncrnas)*length(tfs))
		write.table(rbind(tmp),args[4],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
	}
	}

	p <- summary(fit)$coefficients[-1,4]
	var <- names(p)
	#var <- relevant.x[grep("ENSG",relevant.x)]
	for (n in 1:length(var)){
                tmp <- c(colnames(y),var[n],as.character(p[n]),length(lncrnas),length(tfs),as.numeric(p[n])*length(lncrnas)*length(tfs))
                write.table(rbind(tmp),args[4],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=T)
        }	
}

#write.table(result,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
#write.table(ols,args[6],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
