# lncRNARN

Xiaoman Xie [xiaoman2@illinois.edu] and Saurabh Sinha [saurabh.sinha@bme.gatech.edu]

University of Illinois Urbana-Champaign
Georgia Institute of Technology

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Tutorial](#tutorial)
4. [Acknowledgments and credits](#acknowledgments-and-credits)


## Introduction

This repository provides script and data to repeat analyses performed in the paper "Quantitative estimates of the regulatory influence of long non-coding RNAs on global gene expression variation using TCGA breast cancer transcriptomic data".



![Method Overview](images/Figure1A.png)

In this study, we analyzed transcriptome data of 13,963 protein-coding genes and 1,079 lncRNAs from 1,217 breast cancer samples from TCGA. We used linear regression models to determine how lncRNAs influence gene expression, considering potential confounding effects of transcription factors and covariates. We further refined our analysis by focusing on lncRNAs associated with target genes through specific mechanisms such as lncRNAs that overlap the gene, or lncRNAs that are located in the same Topologically Associating Domain (TAD) as the target gene, or lncRNAs that may compete with the target gene for the same microRNA, i.e., the ceRNA mechanism. To manage model complexity, we applied the Elastic Net approach and validated our models on unseen data. We also explored interactions between lncRNAs and RNA-binding proteins/transcription factors. Combining these analyses, we created a lncRNA-mRNA regulatory network, identifying key lncRNA-mRNA interactions supported by our models and additional regulatory evidence.

### Installation

### Data
Expression profile downloaded from TCGA portal is in file ./data/TCGA_matadata.txt

List of lncRNA Ensembl IDs used in the study: lncRNA.txt

List of protein coding gene Ensembl IDs used in the study: ProteinCodingGene.txt

List of human transcription factors (TF): TF.txt

List of RNA binding proteins (RBP): RBP.txt


### Code
ElasticNetOLS_res.r: Script used to training and select the most significant lncRNA predictor.

OLS_res_Condition.r: Script used to training the full model.
