raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                           header=TRUE, stringsAsFactors=FALSE,
                           na.strings=c("","NA"),
                           strip.white=TRUE)

Diversity_df<-unique(raw_WY_dataframe[c(10, 22, 29)])
Diversity_df<-Diversity_df[which(!is.na(Diversity_df$Specific)), ]
store<-Diversity_df
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


###Species Diversity using Shannon-Weiner###
require(vegan)
Shannon_Wiener_H<-data.frame(diversity(Diversity_matrix, index="shannon"))
Shannon_Wiener_H$Specific<-row.names(Shannon_Wiener_H)
Shannon_Wiener_H<-Shannon_Wiener_H[order(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..), ]

###Checking to make sure SW index does not correllate with sampling
substrate_sampling<-data.frame(table(unique(store[c(1,2)])[1]))
colnames(substrate_sampling)[1]<-"Specific"
Shannon_Wiener_H<-merge(Shannon_Wiener_H, substrate_sampling, by="Specific")

summary(lm(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon.. ~
           Shannon_Wiener_H$Freq))
#new adjusted R2 is .2862 (less than before)


write.table(Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_substr_table.tsv", sep="\t", quote=FALSE, row.name=FALSE)
###same exact code was replicated with $IsoTemp instead of $Substrate
write.table(Shannon_Wiener_H_2, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp.tsv", sep="\t", quote=FALSE, row.name=FALSE)









#### for main text figure ###
par(mar=c(6,5,1,5))
#par(mfrow=c(2,1))

layout(matrix(c(1:2, 0, 0), nrow=1, ncol=2, byrow=TRUE), widths = c(2.5,1.5))
layout.show(n=2)  # to inspect layout                   

plot(Shannon_Wiener_H[,2], xaxt='n', ylab="Shannon-Wiener (H')", 
     main=NA, xlab=NA, cex=1, pch=16)
axis(side=1, at=c(1:nrow(Shannon_Wiener_H)),labels=Shannon_Wiener_H[,1], cex.axis=.65, las=2)
par(new=T)
plot(Shannon_Wiener_H[,3], xaxt='n', ylab=NA, xlab=NA, pch=1, cex=1,
     yaxt='n')
axis(side=4, at=seq(50, 500, 50), labels=seq(50, 500, 50))
mtext("Sampling density", side=4,cex.lab=1, line=2.5)
legend("topleft", legend=c("H'", "sampling density"), pch=c(16, 1))

plot(x=Shannon_Wiener_H_2$IsoTemp, y=Shannon_Wiener_H_2[,2], ylab="Shannon-Wiener (H')", 
     main=NA,cex=1, pch=16, xlab= "Isolation temperature")

quartz.save("~/SurveyPaper/Figures/Shannon_weiner_by_substrate_and_temperature.pdf", type="pdf")

###### RAREFACTION ANALYSES - 2-2-21 KJF

#Setting up a matrix appropriate for specaccum curves (setID as sites)
unique_samples<-unique(raw_WY_dataframe[c(29, 22)])

spec_accum_mat<-matrix(nrow=length(unique(unique_samples$SetID)),
                                   ncol=length(unique(unique_samples$Species)))
#Not sure if the data should be binary or cumulative. 
#Will start with binary (present/absent in sample)
colnames(spec_accum_mat)<-unique(unique_samples$Species)
row.names(spec_accum_mat)<-unique(unique_samples$SetID)
for(i in 1:nrow(unique_samples)){
  col<-which(colnames(spec_accum_mat)==unique_samples$Species[i])
  row<-which(row.names(spec_accum_mat)==unique_samples$SetID[i])
  spec_accum_mat[row, col]<-1
}
spec_accum_mat[is.na(spec_accum_mat)] <- 0


#I wanted to see if the cumulative matrix gives a different curve
unique_samples<-data.frame(table(raw_WY_dataframe$Species, raw_WY_dataframe$SetID))
colnames(unique_samples)[1:2]<-c("Species", "SetID")

cum_spec_accum_mat<-matrix(nrow=length(unique(unique_samples$SetID)),
                       ncol=length(unique(unique_samples$Species)))

colnames(cum_spec_accum_mat)<-unique(unique_samples$Species)
row.names(cum_spec_accum_mat)<-unique(unique_samples$SetID)
for(i in 1:nrow(unique_samples)){
  col<-which(colnames(cum_spec_accum_mat)==unique_samples$Species[i])
  row<-which(row.names(cum_spec_accum_mat)==unique_samples$SetID[i])
  cum_spec_accum_mat[row, col]<-unique_samples[i,3]
}

#with each additional sample (setID)
rarefied<-specaccum(cum_spec_accum_mat, method="rarefaction")
#Ok - so exact is sample-based and rarefaction is individual based
#and the difference seems to be the CI. 
plot(rarefied, ci.type = "polygon", ci.lty = 0, ci.col = "grey", ylab="Yeast OTUs", xlab="Samples")
#I think the correct curve to use is the rarefied individual one. 
quartz.save("~/SurveyPaper/Figures/rarefaction_all_data_by_individual.pdf", type="pdf")



#setID rarefaction curves. These are not really 
#relevant since only one of each morphotype was selected
#we would not expect these to saturate. 
rarecurve(cum_spec_accum_mat, label = FALSE, col = "grey45", ylab = "No. yeast OTUs", xlab="No. isolates in sample")
quartz.save("~/SurveyPaper/Figures/rarefaction_by_individual_by_setID.pdf", type="pdf")


#What I should also do is rarefy by subphylum (higher taxanomic group)
unique_samples<-data.frame(table(raw_WY_dataframe$Species, raw_WY_dataframe$Subphylum))
colnames(unique_samples)[1:2]<-c("Species", "Subphylum")

phy_cum_spec_accum_mat<-matrix(nrow=length(unique(unique_samples$Subphylum)),
                           ncol=length(unique(unique_samples$Species)))

colnames(phy_cum_spec_accum_mat)<-unique(unique_samples$Species)
row.names(phy_cum_spec_accum_mat)<-unique(unique_samples$Subphylum)
for(i in 1:nrow(unique_samples)){
  col<-which(colnames(phy_cum_spec_accum_mat)==unique_samples$Species[i])
  row<-which(row.names(phy_cum_spec_accum_mat)==unique_samples$Subphylum[i])
  phy_cum_spec_accum_mat[row, col]<-unique_samples[i,3]
}

#by_phy<-
rarecurve(phy_cum_spec_accum_mat,label = FALSE, col =color_pal, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5 )
#take away labels and fix by color, lty. 

##Ok - now by substrate type. 
by_substrate<-
rarecurve(Diversity_matrix,label = FALSE, ylab = "No. yeast OTUs", xlab="No. isolates" )
#break plot up into bins for different sampling densities?
deep<-Shannon_Wiener_H$Specific[which(Shannon_Wiener_H$Freq>50)]
deep_DM<-Diversity_matrix[which(row.names(Diversity_matrix) %in% deep),]
shallow<-Shannon_Wiener_H$Specific[which(Shannon_Wiener_H$Freq<50)]
shallow_DM<-Diversity_matrix[which(row.names(Diversity_matrix) %in% shallow),]
##
rarecurve(deep_DM,label = FALSE, ylab = "No. yeast OTUs", xlab="No. isolates", lty=)
rarecurve(shallow_DM,label = FALSE, ylab = "No. yeast OTUs", xlab="No. isolates", lty=1)



color_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


