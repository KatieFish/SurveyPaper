
Run_WY_association_permutations<-function(Permutation_df, permutations){

#each iteration of the loop will set up a new dataframe by permuting the yeast species
#then count each combination, then merge those counts into the Permutation_df
for (i in 1:permutations){
  #the first column (plant genus) will stay the same. So we'll set that 
  #up in a new temporary dataframe. 
  temporary_df<-Permutation_df[1]
  #We're going to shuffle up yeast using the sample function and use that 
  #to create a new column in the dataframe.
  temporary_df$randomized_colum2<-sample(Permutation_df[,2],
                                        replace=FALSE, prob=NULL)
  permuted_df<-data.frame(table(temporary_df[,1],
                                temporary_df[,2]))
  #use the table function to count up our permuted counts. 
  #The new column of counts will be named i, which will correspond 
  # to the permutaiton. 
  colnames(permuted_df)<-c(colnames(Permutation_df[1:2]), i)
  #now we tack it onto our Permutation_df using the merge function
  Permutation_df<-merge(Permutation_df, permuted_df, by=colnames(Permutation_df[1:2]))
}

observed_pairs<-Permutation_df[which(Permutation_df$Observed>0),]
observed_pairs$pval<-NA
for(i in 1:nrow(observed_pairs)){
  permuted_counts<-as.numeric(observed_pairs[i, 4:(ncol(observed_pairs)-1)])
  observed_pairs$pval[i]<-length(which(permuted_counts>=observed_pairs$Observed[i]))/
    (permutations)
}

return(observed_pairs)
}
