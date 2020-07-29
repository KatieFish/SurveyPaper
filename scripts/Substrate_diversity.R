raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
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
Shannon_Wiener_H$Specific<-row.names(Shannon_Wiener_H)
Shannon_Wiener_H<-Shannon_Wiener_H[order(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..), ]
par(mar=c(10,5,2,5))
plot(Shannon_Wiener_H[,2], xaxt='n', ylab="Shannon-Wiener (H')", 
     main=NA, xlab=NA, cex=1, pch=16)
axis(side=1, at=c(1:nrow(Shannon_Wiener_H)),labels=Shannon_Wiener_H[,1], cex.axis=.65, las=2)
par(new=T)
plot(Shannon_Wiener_H[,3], xaxt='n', ylab=NA, xlab=NA, pch=1, cex=1,
     yaxt='n')
axis(side=4, at=seq(50, 500, 50), labels=seq(50, 500, 50))
mtext("Sampling density", side=4,cex.lab=1, line=2.5)
legend("topleft", legend=c("H'", "sampling density"), pch=c(16, 1))
quartz.save("~/SurveyPaper/Figures/Shannon_weiner_by_substrate.pdf", type="pdf")

###Checking to make sure SW index does not correllate with sampling
substrate_sampling<-data.frame(table(Diversity_df$Specific))
colnames(substrate_sampling)[1]<-"Specific"
Shannon_Wiener_H<-merge(Shannon_Wiener_H, substrate_sampling, by="Specific")

summary(lm(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon.. ~
           Shannon_Wiener_H$Freq))


ggplot(Shannon_Wiener_H, aes(x=Freq, y=diversity.Diversity_matrix..index....shannon..)) +
  geom_point(shape=1) +   
  geom_smooth(method=lm) 