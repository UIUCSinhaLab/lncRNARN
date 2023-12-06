# lncRNARN

Xiaoman Xie [xiaoman2@illinois.edu] and Saurabh Sinha

University of Illinois Urbana-Champaign

This repository provides script and data to repeat analyses performed in the paper "Quantitative estimates of the regulatory influence of long non-coding RNAs on global gene expression variation using TCGA breast cancer transcriptomic data".

![Method Overview](images/Figure1A.png)

### Data
Expression profile downloaded from TCGA portal is in file ./data/TCGA_matadata.txt

List of lncRNA Ensembl IDs used in the study: lncRNA.txt

List of protein coding gene Ensembl IDs used in the study: ProteinCodingGene.txt

List of human transcription factors (TF): TF.txt

List of RNA binding proteins (RBP): RBP.txt


### Code
ElasticNetOLS_res.r: Script used to training and select the most significant lncRNA predictor.

OLS_res_Condition.r: Script used to training the full model.
