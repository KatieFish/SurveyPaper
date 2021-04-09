##Last annotated 04-09-21 KJF WJS


#This function was used to completely permute combinations of OTUs and substrate, 
#substrate type, or isolation temperature without replacement. 

#Arguments

#all_observerations_dataframe - raw data (table S1)

#permutations - number of desired permutations (10,000)

#colnames_to_permute - list of names of 2 columns for specific analysis (e.g. OTU and specific substrate)



Run_WY_association_permutations<-function(all_observations_dataframe, permutations, colnames_to_permute){


col1<-which(colnames(all_observations_dataframe)==colnames_to_permute[1])
col2<-which(colnames(all_observations_dataframe)==colnames_to_permute[2])
#
Permutation_df<-data.frame(table(all_observations_dataframe[,col1],
                                   all_observations_dataframe[,col2]))
colnames(Permutation_df)<-c(colnames_to_permute[1], colnames_to_permute[2], "Observed")
#set up column to tally totals in

to_permute_df<-all_observations_dataframe[c(col1, col2)]
#raw data to permute
  
#each iteration of the loop will set up a new dataframe by permuting the yeast species
#then count each combination, then merge those counts into the Permutation_df
for (i in 1:permutations){
  #the first column will stay the same. So we'll set that 
  #up in a new temporary dataframe. 
  temporary_df<-all_observations_dataframe[col1]
  #We're going to shuffle up yeast using the sample function and use that 
  #to create a new column in the dataframe.
  temporary_df$randomized_colum2<-sample(all_observations_dataframe[,col2],
                                        replace=FALSE, prob=NULL)
  permuted_df<-data.frame(table(temporary_df[,1],
                                temporary_df[,2]))
  #use the table function to count up our permuted counts. 
  #The new column of counts will be named i, which will correspond 
  # to the permutation. 
  colnames(permuted_df)<-c(colnames(Permutation_df[1:2]), i)
  #now we tack it onto our Permutation_df using the merge function
  Permutation_df<-merge(Permutation_df, permuted_df, by=colnames(Permutation_df[1:2]))
}

#retain only combinations that were observed in the data
observed_pairs<-Permutation_df[which(Permutation_df$Observed>0),]
#calculate pval as the # of permuations >= than the observed values / total permutations
observed_pairs$pval<-NA
for(i in 1:nrow(observed_pairs)){
  permuted_counts<-as.numeric(observed_pairs[i, 4:(ncol(observed_pairs)-1)])
  observed_pairs$pval[i]<-length(which(permuted_counts>=observed_pairs$Observed[i]))/
    (permutations)
}

observed_pairs<-observed_pairs[c(1:3, ncol(observed_pairs),4:(ncol(observed_pairs)-1))]

return(observed_pairs)
#returns raw pvalues which were later corrected using BH adjustment
}
