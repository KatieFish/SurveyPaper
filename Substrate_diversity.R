raw_WY_dataframe<-read.csv("~/SurveyPaper/data/WY_df_2018-02-08.csv",
                           header=TRUE, stringsAsFactors=FALSE,
                           na.strings=c("","NA"),
                           strip.white=TRUE)

Diversity_df<-unique(raw_WY_dataframe[c(10, 22, 29)])
Diversity_df<-data.frame(table(Diversity_df$Species, Diversity_df$Specific))

colnames(Diversity_df)<-c("Species", "Specific", "number")
#Set up matrix
#rows are (OTUs)
#columns are substratespecies 
Diversity_matrix<-matrix(nrow=length(unique(Diversity_df$Species)), ncol=length(unique(Diversity_df$Specific)))
colnames(Diversity_matrix)<-unique(Diversity_df$Specific)
row.names(Diversity_matrix)<-unique(Diversity_df$Species)

#populate matrix
for (i in 1:nrow(Diversity_matrix)){
  ndf<-Diversity_df[which(Diversity_df$Species==row.names(Diversity_matrix)[i]), ]
  for (j in 1:ncol(Diversity_matrix)) {
    ndf1<-ndf[which(ndf$Specific==colnames(Diversity_matrix)[j]),]
    Diversity_matrix[i,j]<-ndf1$number[1]
  }
}
rm(ndf, ndf1, i, j)
Diversity_matrix<-t(Diversity_matrix)

###Species Richness using Menhinicks_D
OTUs<-apply(Diversity_matrix>0,1,sum)
Total_isolates<-apply(Diversity_matrix,1,sum)
Menhinicks_D<-data.frame(OTUs/(sqrt(Total_isolates)))
par(mar=c(8, 4, 4, 4))
plot(Menhinicks_D[,1], xaxt= 'n', ylab = "Menhinick's Index (D)", 
     main = "Species Richness by Substrate", xlab=NA, cex=1, pch=16)
axis(side=1, at=c(1:nrow(Menhinicks_D)), labels = row.names(Menhinicks_D), las=2)


###Species Diversity using Shannon-Weiner###
require(vegan)
Shannon_Wiener_H<-data.frame(diversity(Diversity_matrix, index="shannon"))
plot(Shannon_Wiener_H[,1], xaxt='n', ylab="Shannon-Wiener (H')", 
     main="Species Diversity by Substrate", xlab=NA, cex=.5, pch=16)
axis(side=1, at=c(1:nrow(Shannon_Wiener_H)),labels=row.names(Shannon_Wiener_H), las=2)
