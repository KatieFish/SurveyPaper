
Sampling_FDR<-function(Permutation_df, permutations, pvals, perms_to_sample){
  Permdf<-Permutation_df[c(1:(permutations+3))]
  
  FDR_df<-data.frame(pvals)
  FDR_df$observed_significant<-NA
  FDR_df[3:(2+perms_to_sample)]<-NA

  critvals<-data.frame(Permdf[c(1,2)]) 
  critvals[3:(length(pvals)+2)]<-NA
  colnames(critvals)[3:(length(pvals)+2)]<-pvals
  for (i in 1:nrow(critvals)){
    simvals<-sort(as.numeric(Permdf[i, 4:ncol(Permdf)]),
                  decreasing=TRUE)
    for (j in 3:ncol(critvals)){
    critvals[i,j]<-simvals[permutations*as.numeric(colnames(critvals)[j])]
  }
  }
  
  for (i in 1:nrow(FDR_df)){
    for(j in 2:(perms_to_sample+2)){
      critval_col<-which(colnames(critvals)==FDR_df$pvals[i])
      FDR_df[i,j]<-length(which(Permdf[,(j+1)]>=critvals[,critval_col]))
      z<-Permdf[which(Permdf[,(j+1)]>=critvals[,critval_col]), c(1,2,3)]
    }
}
return(FDR_df)
}
