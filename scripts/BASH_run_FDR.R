#!/usr/bin/env Rscript


source("~/SurveyPaper/scripts/Sample_FDR_WYassociations.R")

perms<-read.delim("~/SurveyPaper/data/Substrate_Specific_Associations_10k_permutations.tsv", header=TRUE, stringsAsFactors=FALSE)

permsFDR<-Sampling_FDR(Permutation_df=perms, permutations=10000, pvals=c(0.0001, 0.001, 0.01), perms_to_sample=100)

write.table(permsFDR, "~/SurveyPaper/data/Substrate_Specific_Associations_FDR.tsv", sep="\t", quote=FALSE, row.names=FALSE)