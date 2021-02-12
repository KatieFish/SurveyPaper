require(vegan)
require(ggplot2)

#color pal for plotting downstream
color_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#read in raw data
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)


###################################################
#SHANNON WEINER ANALYSES
###################################################

####################################################
#Diversity by substrate
#####################################################

#retain unique combinations of setID-specific substrate-species
Diversity_df<-unique(raw_WY_dataframe[c(10, 22, 29)])
#eliminate NA specific categories
Diversity_df<-Diversity_df[which(!is.na(Diversity_df$Specific)), ]
#keep this data for later
store<-Diversity_df
# quantify how many times each species-specific substrate is present
Diversity_df<-data.frame(table(Diversity_df$Species, Diversity_df$Specific))
colnames(Diversity_df)<-c("Species", "Specific", "number")


#Set up matrix for vegan analyses
#rows are (OTUs)
#columns are substrates 
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
#flip matrix so rows are substrates and columns are species
#values are #of indendent times each species is found on each substrate
Diversity_matrix<-t(Diversity_matrix)

#find unique times each substrate was sampled (is greater than number of unique samples
# or setIDs by 41. There are 41 setIDs that are assigned to more than one specific substrate)
substrate_sampling<-data.frame(table(unique(store[c(1,2)])[1]))
colnames(substrate_sampling)[1]<-"Specific"
# removing the substrates that are sampled less than 5 times. 
# leaves the 15 top-sampled substrates
over5<-substrate_sampling$Specific[which(substrate_sampling$Freq>5)]
Diversity_matrix<-Diversity_matrix[which(
  row.names(Diversity_matrix) %in% over5),]
#store as seperate dataframe for downstream
Diversity_matrix->substrate_spp_accum_mat

###Species Diversity using Shannon-Weiner###
# find Shannon-Wiener index using Vegan
Shannon_Wiener_H<-data.frame(diversity(Diversity_matrix, index="shannon"))
Shannon_Wiener_H$Specific<-row.names(Shannon_Wiener_H)
Shannon_Wiener_H<-Shannon_Wiener_H[order(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..), ]

###Add sampling depth to dataframe 
Shannon_Wiener_H<-merge(Shannon_Wiener_H, substrate_sampling, by="Specific")

#Linear regression to see how sampling density affects H'
plot(x=Shannon_Wiener_H$Freq, y=Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..)
summary(lm(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon.. ~
             Shannon_Wiener_H$Freq))

write.table(Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_substr_table.tsv", sep="\t", quote=FALSE, row.name=FALSE)


##########################################################
#Diversity by IsoTemp
##########################################################

### same exact code was replicated with $IsoTemp instead of $Substrate
#retain unique combinations of setID-IsoTemp-species
Diversity_df<-unique(raw_WY_dataframe[c(23, 22, 29)])
#eliminate unkown IsoTemp (only 1 isolate)
Diversity_df<-Diversity_df[which(!Diversity_df$IsoTemp=="Unknown"), ]
#keep this data for later
store<-Diversity_df
# quantify how many times each species-IsoTemp substrate is present
Diversity_df<-data.frame(table(Diversity_df$Species, Diversity_df$IsoTemp))
colnames(Diversity_df)<-c("Species", "IsoTemp", "number")


#Set up matrix for vegan analyses
#rows are (OTUs)
#columns are IsoTemps 
Diversity_matrix<-matrix(nrow=length(unique(Diversity_df$Species)),
                         ncol=length(unique(Diversity_df$IsoTemp)))
colnames(Diversity_matrix)<-unique(Diversity_df$IsoTemp)
row.names(Diversity_matrix)<-unique(Diversity_df$Species)

#populate matrix
for (i in 1:nrow(Diversity_matrix)){
  ndf<-Diversity_df[which(
    Diversity_df$Species==row.names(Diversity_matrix)[i]), ]
  for (j in 1:ncol(Diversity_matrix)) {
    ndf1<-ndf[which(
      ndf$IsoTemp==colnames(Diversity_matrix)[j]),]
    Diversity_matrix[i,j]<-ndf1$number[1]
  }
}
rm(ndf, ndf1, i, j)
#flip matrix so rows are IsoTemps and columns are species
#values are #of indendent times each species is found with each IsoTemp
Diversity_matrix<-t(Diversity_matrix)
#store as seperate dataframe for downstream
Diversity_matrix->isoTemp_spp_accum_mat

#find unique times each IsoTemp was sampled. 
#There are 421 setIDs (samples) that are assigned to more than one IsoTemp)
substrate_sampling<-data.frame(table(unique(store[c(1,2)])[1]))
colnames(substrate_sampling)[1]<-"IsoTemp"

###Species Diversity using Shannon-Weiner###
# find Shannon-Wiener index using Vegan
Shannon_Wiener_H_2<-data.frame(diversity(Diversity_matrix, index="shannon"))
Shannon_Wiener_H_2$IsoTemp<-row.names(Shannon_Wiener_H_2)
Shannon_Wiener_H_2<-Shannon_Wiener_H_2[order(Shannon_Wiener_H_2$IsoTemp), ]

###Add sampling depth to dataframe 
Shannon_Wiener_H_2<-merge(Shannon_Wiener_H_2, substrate_sampling, by="IsoTemp")

#Not as concerned about sampling density here as all 
#temps were sampled extensively.
write.table(Shannon_Wiener_H_2, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp.tsv", sep="\t", quote=FALSE, row.name=FALSE)


#####################################################
# SW by temp and phylum
#####################################################
#First pull out unique temp-species-setID combinations again
Diversity_df<-unique(raw_WY_dataframe[c(23, 22, 29)])
Diversity_df<-Diversity_df[which(!Diversity_df$IsoTemp=="Unknown"),]
#generate lists of spp. belonging to each group
Basidios<-unique(raw_WY_dataframe$Species[which(raw_WY_dataframe$Subphylum=="Basidiomycota")])
Pezizos<-unique(raw_WY_dataframe$Species[which(raw_WY_dataframe$Subphylum=="Pezizomycotina")])
Saccharos<-unique(raw_WY_dataframe$Species[which(raw_WY_dataframe$Subphylum=="Saccharomycotina")])

store<-Diversity_df
# quantify how many times each species-IsoTemp substrate is present
Diversity_df<-data.frame(table(Diversity_df$Species, Diversity_df$IsoTemp))
colnames(Diversity_df)<-c("Species", "IsoTemp", "number")


#Set up matrix for vegan analyses
#rows are (OTUs)
#columns are IsoTemps 
Diversity_matrix<-matrix(nrow=length(unique(Diversity_df$Species)),
                         ncol=length(unique(Diversity_df$IsoTemp)))
colnames(Diversity_matrix)<-unique(Diversity_df$IsoTemp)
row.names(Diversity_matrix)<-unique(Diversity_df$Species)

#populate matrix
for (i in 1:nrow(Diversity_matrix)){
  ndf<-Diversity_df[which(
    Diversity_df$Species==row.names(Diversity_matrix)[i]), ]
  for (j in 1:ncol(Diversity_matrix)) {
    ndf1<-ndf[which(
      ndf$IsoTemp==colnames(Diversity_matrix)[j]),]
    Diversity_matrix[i,j]<-ndf1$number[1]
  }
}
rm(ndf, ndf1, i, j)
#flip matrix so rows are IsoTemps and columns are species
#values are #of indendent times each species is found with each IsoTemp
Diversity_matrix<-t(Diversity_matrix)

###Basidio specific analysis
Basidio_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Basidios)]
Basidio_Shannon_Wiener_H<-data.frame(diversity(Basidio_div_matrix, index="shannon"))
Basidio_Shannon_Wiener_H$IsoTemp<-row.names(Basidio_Shannon_Wiener_H)
Basidio_Shannon_Wiener_H<-Basidio_Shannon_Wiener_H[order(Basidio_Shannon_Wiener_H$IsoTemp), ]

write.table(Basidio_Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp_Basidios_only.tsv", sep="\t", quote=FALSE, row.name=FALSE)

###Pezizo specific analysis
Pezizo_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Pezizos)]
Pezizo_Shannon_Wiener_H<-data.frame(diversity(Pezizo_div_matrix, index="shannon"))
Pezizo_Shannon_Wiener_H$IsoTemp<-row.names(Pezizo_Shannon_Wiener_H)
Pezizo_Shannon_Wiener_H<-Pezizo_Shannon_Wiener_H[order(Pezizo_Shannon_Wiener_H$IsoTemp), ]

write.table(Pezizo_Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp_Pezizos_only.tsv", sep="\t", quote=FALSE, row.name=FALSE)

###Saccharo specific analysis
Saccharo_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Saccharos)]
Saccharo_Shannon_Wiener_H<-data.frame(diversity(Saccharo_div_matrix, index="shannon"))
Saccharo_Shannon_Wiener_H$IsoTemp<-row.names(Saccharo_Shannon_Wiener_H)
Saccharo_Shannon_Wiener_H<-Saccharo_Shannon_Wiener_H[order(Saccharo_Shannon_Wiener_H$IsoTemp), ]

write.table(Saccharo_Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp_Saccharos_only.tsv", sep="\t", quote=FALSE, row.name=FALSE)


#########################################################
#RAREFACTION CURVES
#########################################################

#############################
#Total dataset rarefaction 
#############################

#quantify sp-setID combinations
unique_samples<-data.frame(table(raw_WY_dataframe$Species, raw_WY_dataframe$SetID))
colnames(unique_samples)[1:2]<-c("Species", "SetID")

#Matrix with SetIDs (sampleID) as rows and species as columns
cum_spec_accum_mat<-matrix(nrow=length(unique(unique_samples$SetID)),
                           ncol=length(unique(unique_samples$Species)))
colnames(cum_spec_accum_mat)<-unique(unique_samples$Species)
row.names(cum_spec_accum_mat)<-unique(unique_samples$SetID)

#populate matrix with species found in each sample
for(i in 1:nrow(unique_samples)){
  col<-which(colnames(cum_spec_accum_mat)==unique_samples$Species[i])
  row<-which(row.names(cum_spec_accum_mat)==unique_samples$SetID[i])
  cum_spec_accum_mat[row, col]<-unique_samples[i,3]
}


#Individual-based rarefaction curve across all data
#According to vegan docs, 
#"rarefaction" finds the mean when accumulating individuals instead of sites
rarefied<-specaccum(cum_spec_accum_mat, method="rarefaction") 
plot(rarefied, ci.type = "polygon", ci.lty = 0, ci.col = "grey", ylab="No. yeast OTUs", xlab="Samples")
#confidence intervals come from standard deviation 
#sd is found analytically within the function
quartz.save("~/SurveyPaper/Figures/rarefaction_all_data_by_individual.pdf", type="pdf")

#######################
#Rarefaction by sample 
#######################

#setID rarefaction curves. These are not really 
#relevant since only one of each morphotype was selected
#we would not expect these to saturate. 
rarecurve(cum_spec_accum_mat, label = FALSE, col = "grey45", ylab = "No. yeast OTUs", xlab="No. isolates in sample")
quartz.save("~/SurveyPaper/Figures/rarefaction_by_individual_by_setID.pdf", type="pdf")


#########################
#Rarefaction by subphylum 
#########################
#Table of species recovered by subphylum
unique_samples<-data.frame(table(raw_WY_dataframe$Species, raw_WY_dataframe$Subphylum))
colnames(unique_samples)[1:2]<-c("Species", "Subphylum")

#matrix with subphyla as rows and species as columns
phy_cum_spec_accum_mat<-matrix(nrow=length(unique(unique_samples$Subphylum)),
                               ncol=length(unique(unique_samples$Species)))
colnames(phy_cum_spec_accum_mat)<-unique(unique_samples$Species)
row.names(phy_cum_spec_accum_mat)<-unique(unique_samples$Subphylum)
#populate matrix with totals
for(i in 1:nrow(unique_samples)){
  col<-which(colnames(phy_cum_spec_accum_mat)==unique_samples$Species[i])
  row<-which(row.names(phy_cum_spec_accum_mat)==unique_samples$Subphylum[i])
  phy_cum_spec_accum_mat[row, col]<-unique_samples[i,3]
}
rarecurve(phy_cum_spec_accum_mat,label = FALSE, col =color_pal, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5 )
quartz.save("~/SurveyPaper/Figures/rarefaction_by_subphylum.pdf", type="pdf")


#############################################
#Rarefaction to accompany diversity estimates
#############################################

#########################
#Rarefaction by substrate
#########################
rarecurve(substrate_spp_accum_mat,label = FALSE, col="grey45", ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5 )
quartz.save("~/SurveyPaper/Figures/rarefaction_by_substrate.pdf", type="pdf")


#######################
#Rarefaction by IsoTemp
#######################
rarecurve(isoTemp_spp_accum_mat, label=TRUE, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5)
quartz.save("~/SurveyPaper/Figures/rarefaction_by_IsoTemp.pdf", type="pdf")


##############
#Figures
##############

#rarefactionfigure1:
par(mar=c(6,5,1,5))
#par(mfrow=c(2,1))
layout(matrix(c(1:2, 0, 0), nrow=1, ncol=2, byrow=TRUE), widths = c(2.5,2.5))
layout.show(n=2)  # to inspect layout                   
plot(rarefied, ci.type = "polygon", ci.lty = 0, ci.col = "grey", ylab="No. yeast OTUs", xlab="Samples")
rarecurve(phy_cum_spec_accum_mat,label = FALSE, col =color_pal, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5 )

quartz.save("~/SurveyPaper/Figures/Sup_rarefaction_Figure_1.pdf", type="pdf")


#H' plots with rarefaction

par(mar=c(5,5,1,3))
layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE), widths = c(2.5, 2.5))
layout.show(n=4)  # to inspect layout                   

#1
plot(Shannon_Wiener_H[,2], xaxt='n', ylab="Shannon-Wiener (H')", 
     main=NA, xlab=NA, cex=1, pch=16)
axis(side=1, at=c(1:nrow(Shannon_Wiener_H)),labels=Shannon_Wiener_H[,1], cex.axis=.65, las=2)
par(new=T)
plot(Shannon_Wiener_H[,3], xaxt='n', ylab=NA, xlab=NA, pch=1, cex=1,
     yaxt='n')
axis(side=4, at=seq(50, 500, 50), labels=seq(50, 500, 50))
mtext("Sampling density", side=4, line=2.5, cex = .8)
legend("topleft", legend=c("H'", "sampling density"), pch=c(16, 1))
#2
rarecurve(substrate_spp_accum_mat,label = FALSE, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=.5 )
#3
plot(x=as.character(Shannon_Wiener_H_2$IsoTemp), y=Shannon_Wiener_H_2[,2], ylab="Shannon-Wiener (H')", 
     main=NA,cex=1, pch=16, xlab= "Isolation temperature")
#4
rarecurve(isoTemp_spp_accum_mat, label=TRUE, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=1)

quartz.save("~/SurveyPaper/Figures/updatedfig3_2-9-21.pdf", type="pdf")



