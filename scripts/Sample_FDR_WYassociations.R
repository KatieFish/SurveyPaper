
Sampling_FDR<-function(Permutation_df, permutations, pvals, perms_to_sample){
  #going to change this from critical-value based approach to determination of p values.
  
  
  
  if("pval" %in% colnames(Permutation_df)){
  Permdf<-Permutation_df[-which(colnames(Permutation_df)=="pval")]
  } else {
    Permdf<-Permutation_df
  }
  
  FDR_df<-data.frame(pvals)
  FDR_df$observed_significant<-NA
  FDR_df[3:(2+perms_to_sample)]<-NA

  #randomly choose permutation columns to sample
  set.seed(43)
  cols<-c(3,sample(4:10003, perms_to_sample, replace=FALSE))
    
  for (i in 1:nrow(FDR_df)){
    for(j in 2:(perms_to_sample+2)){
      sigs<-0
      for(k in 1:nrow(Permdf)){
        sims<-as.numeric(Permdf[k, 4:ncol(Permdf)])
        simpval<-(length(which(sims>=Permdf[k, cols[j-1]])))/permutations
      if(simpval<=FDR_df$pvals[i]){
      sigs<-sigs+1
      }
      }
    FDR_df[i,j]<-sigs
    }
  }
return(FDR_df)
}
