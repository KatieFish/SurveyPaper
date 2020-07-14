Analyze_asin_Opulenteetal2019<-function(raw_permutation_df){
 require(Hmisc)
 
  if("pval" %in% colnames(raw_permutation_df)){
    Permdf<-raw_permutation_df[-which(colnames(raw_permutation_df)=="pval")]
  } else {
    Permdf<-raw_permutation_df
  }
  
  Obs_Exp_df<-Permdf[c(1,2,3)]
  Obs_Exp_df$Expected<-NA
  for(i in 1:nrow(Permdf)){
    permvec<-as.numeric(Permdf[i, 4:ncol(Permdf)])
    Obs_Exp_df$Expected[i]<-mean(permvec)
  }
  Obs_Exp_df$Difference<-Obs_Exp_df$Observed-Obs_Exp_df$Expected
 
 #binomial confidence intervals using hmisc
 #In opulente et al binconf was used to generate a p value 
 # in a manner identical to our generation of a p value 
 # based on permutations. I'm not redoing this because 
 # it is the exact same calculation.
 
results<-merge(raw_permutation_df[c(1,2,4)], Obs_Exp_df, 
                           by=(colnames(raw_permutation_df)[1:2]))

#in instances where observed ==1 and permutations never 
#exceed 1, the pval and the expected count will be the same
#i.e. both = # of perms ==1 / 10000

#This is the correction Opulente et al applied. I'm not 
#sure it is ok since these are not independent pvals. 
results$BH_Padj<-p.adjust(results$pval, method = "BH")  
results<-results[order(results$pval),]
return(results)
}