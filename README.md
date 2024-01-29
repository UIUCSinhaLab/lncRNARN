# lncRNARN

Xiaoman Xie [xiaoman2@illinois.edu] and Saurabh Sinha [saurabh.sinha@bme.gatech.edu]

University of Illinois Urbana-Champaign
Georgia Institute of Technology

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Tutorial](#tutorial)


## Introduction

This repository provides script and data to repeat analyses performed in the paper "Quantitative estimates of the regulatory influence of long non-coding RNAs on global gene expression variation using TCGA breast cancer transcriptomic data".



![Method Overview](images/Figure1A.png)

In this study, we analyzed transcriptome data of 13,963 protein-coding genes and 1,079 lncRNAs from 1,217 breast cancer samples from TCGA. We used linear regression models to determine how lncRNAs influence gene expression, considering potential confounding effects of transcription factors and covariates. We further refined our analysis by focusing on lncRNAs associated with target genes through specific mechanisms such as lncRNAs that overlap the gene, or lncRNAs that are located in the same Topologically Associating Domain (TAD) as the target gene, or lncRNAs that may compete with the target gene for the same microRNA, i.e., the ceRNA mechanism. To manage model complexity, we applied the Elastic Net approach and validated our models on unseen data. We also explored interactions between lncRNAs and RNA-binding proteins/transcription factors. Combining these analyses, we created a lncRNA-mRNA regulatory network, identifying key lncRNA-mRNA interactions supported by our models and additional regulatory evidence.

[Return to TOC](#table-of-contents)

## Installation
Please first clone this repository from Github: 
```
git clone https://github.com/UIUCSinhaLab/lncRNARN.git
```
The expression data is too large for inclusion in this repository. Please download the expression matrices from: [ExpMatrix.zip](https://drive.google.com/file/d/1jV-kezgQVlZndelWc0N6gYV3i0kMjM4R/view?usp=sharing) After downloading, unzip the file. Place all files from the unzipped directory into the 'data' folder.

## Tutorial
This section of the README is meant to walk a user through the process of recreating the lncRNA regulatory network reported in the paper.

### Download and preprocess TCGA RNA-Seq data
Please note that the TCGA GDC portal no longer supports HTSeq counts data. The file ID used in this study is listed in [file](data/TCGA_metadata.txt). We have included the merged expression matrix in [ExpMatrix.zip](https://drive.google.com/file/d/1jV-kezgQVlZndelWc0N6gYV3i0kMjM4R/view?usp=sharing). The preprocessing script used to create this expression matrix can be found [here](code/).

### Step1: Regress out covariates
Regress out impact of covariates (age, sex and race) on gene expression. This step will create a expression matrix same as CPMMatrix.noCov.txt in [ExpMatrix.zip](https://drive.google.com/file/d/1jV-kezgQVlZndelWc0N6gYV3i0kMjM4R/view?usp=sharing).
```
Rscript RegressOutCovariate.r CPMMatrix.noCov.txt
```

### Step2: Indentify and regress out significant TF regulators
ElasticNet regression was used to select the top 10 most significant TF regulators for each protein. This step will create a expression matrix (res.txt) same as CPMMatrix.noCov.noTF.txt in [ExpMatrix.zip](https://drive.google.com/file/d/1jV-kezgQVlZndelWc0N6gYV3i0kMjM4R/view?usp=sharing), a list of significant TF regulators for each gene as in [TFRegulators.txt](data/TFRegulators.txt), and training and test performance with 5-fold cross-validation as in [CVPerformanceTF.txt](data/CVPerformanceTF.txt)

```
Rscript ElasticNetOLSTF.r
```

### Step3: Indentify significant lncRNA regulators
ElasticNet regression was used to select the top 10 most significant lncRNA regulators for each protein. This step will generate a list of significant lncRNA regulators  for each gene as in [lncRNARegulators.txt](data/lncRNARegulators.txt), and training and test performance with 5-fold cross-validation as in [CVPerformance_lncRNA.txt](data/CVPerformance_lncRNA.txt)

```
Rscript ElasticNetOLSlncRNA.r
```

### Step4: Indentify significant RBP regulators
ElasticNet regression was used to select the top 10 most significant RBP regulators for each protein. This step will generate a list of significant RBP regulators  for each gene as in [RBPRegulators.txt](data/RBPRegulators.txt), and training and test performance with 5-fold cross-validation as in [CVPerformanceRBP.txt](data/CVPerformanceRBP.txt)

```
Rscript ElasticNetOLSRBP.r
```

### Step5: Using StepAIC to indentify significant TF:lncRNA and RBP:lncRNA interaction terms
Considering significant lncRNA, TF and RBP regulators detected from the previous steps, AIC was used to identify significant lncRNA-TF and lncRNA-RBP interaction terms. This step will generate a list of significant interaction terms as in [TFInteractions.txt](data/TFInteractions.txt) and [RBPInteractions.txt](data/RBPInteractions.txt).

```
Rscript ElasticNetOLS_TF_lncRNA.r
Rscript ElasticNetOLS_RBP_lncRNA.r
```

### Step6: OLS regression detect significant lncRNA regulators
This step used the significant lncRNA and interactive regulators from previous steps to generate the lncRNA-mRNA regulatory edges as in [FinalNetwork.txt](data/FinalNetwork.txt) and in Supplementary Table S18.

```
Rscript OLS_res_Condition.r
```




