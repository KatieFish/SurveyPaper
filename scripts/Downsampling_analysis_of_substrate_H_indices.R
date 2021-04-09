
require(vegan)

#retain unique combinations of setID-Associated substrate-species
Diversity_df<-unique(raw_WY_dataframe[c(8, 22, 29)])
downsample_H_indices<-data.frame(unique(Diversity_df$Associated))
colnames(downsample_H_indices)[1]<-"Associated"


for(k in 1:100){
####downsample Diversity_df
downsampled_df<-data.frame(Associated=character(),
                           SetID=character(), 
                           Species=character(), 
                           stringsAsFactors=FALSE) 

subs_above_20<-c("Plant","Tree","Soil","Fungus") 

for(i in 1:length(subs_above_20)){
  ndf<-Diversity_df[which(Diversity_df$Associated==subs_above_20[i]),]
  rows<-sample(1:nrow(ndf), size = 20, replace=FALSE)
  downsampled_df<-rbind(downsampled_df, ndf[rows,])
  }


##downsampling done

#eliminate NA Associated categories
downsampled_df<-downsampled_df[which(!is.na(downsampled_df$Associated)), ]

downsampled_df<-data.frame(table(downsampled_df$Species, downsampled_df$Associated))
colnames(downsampled_df)<-c("Species", "Associated", "number")

Diversity_matrix<-matrix(nrow=length(unique(downsampled_df$Species)), ncol=length(unique(downsampled_df$Associated)))
colnames(Diversity_matrix)<-unique(downsampled_df$Associated)
row.names(Diversity_matrix)<-unique(downsampled_df$Species)

#populate matrix
for (i in 1:nrow(Diversity_matrix)){
  ndf<-downsampled_df[which(downsampled_df$Species==row.names(Diversity_matrix)[i]), ]
  for (j in 1:ncol(Diversity_matrix)) {
    ndf1<-ndf[which(ndf$Associated==colnames(Diversity_matrix)[j]),]
    Diversity_matrix[i,j]<-ndf1$number[1]
  }
}
rm(ndf, ndf1, i, j)
#flip matrix so rows are substrates and columns are species
#values are #of indendent times each species is found on each substrate
Diversity_matrix<-t(Diversity_matrix)

Shannon_Wiener_H<-data.frame(diversity(Diversity_matrix, index="shannon"))
colnames(Shannon_Wiener_H)[1]<-sprintf("%d ds",k)
Shannon_Wiener_H$Associated<-row.names(Shannon_Wiener_H)
downsample_H_indices<-merge(downsample_H_indices, Shannon_Wiener_H, by="Associated", all=TRUE)
}

library(reshape2)
downsample_H_indices<-melt(downsample_H_indices)
boxplot(downsample_H_indices$value ~ downsample_H_indices$Associated)

