# Project: GBMpanel

# Description:
There are the scripts for scientific article titled "Integrative analysis of gene expression and DNA methylation based on the one-class logistic regression machine learning approach identifies stemness features in medulloblastoma". The scripts are encoded using R. They can generate the result of figures and tables from the article.

# Requirement:
R(https://www.r-project.org/, version >= 3.4) and R package as below:
* survival
* survminer
* forestplot
* survivalROC
* survMisc
* ggplot2
* GGally
* parallel
* glmnet
* caret
* dplyr

# Usage:
Please put the supplementary_append.Rdata and the \*.R scripts in the same directory and change the workdir of R in the same directory.
## Run the scripts: 
Rscript \<the R script\>
## Or open the R software then enter on the command line:
load("\<the R script\>")
