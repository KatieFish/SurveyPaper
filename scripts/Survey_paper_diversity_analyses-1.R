require(vegan)
require(ggplot2)
library(pheatmap)
library(ggplotify)



#color pal for plotting downstream
color_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
colfunc <- colorRampPalette(c("khaki", "brown4"))
colfunc(15)->col_grad
plot(rep(1,8),col=color_pal,pch=19,cex=3)
library(RColorBrewer)
col_grad_2<-brewer.pal(name="YlOrBr", n=9)
col_pal_2<-c(brewer.pal(name="BrBG", n=11),
             brewer.pal(name="RdYlBu", n=4))

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
Diversity_df<-unique(raw_WY_dataframe[c(8, 22, 29)])
#eliminate NA specific categories
Diversity_df<-Diversity_df[which(!is.na(Diversity_df$Associated)), ]
#keep this data for later
store<-Diversity_df
# quantify how many times each species-specific substrate is present
Diversity_df<-data.frame(table(Diversity_df$Species, Diversity_df$Associated))
colnames(Diversity_df)<-c("Species", "Associated", "number")


#Set up matrix for vegan analyses
#rows are (OTUs)
#columns are substrates 
Diversity_matrix<-matrix(nrow=length(unique(Diversity_df$Species)), ncol=length(unique(Diversity_df$Associated)))
colnames(Diversity_matrix)<-unique(Diversity_df$Associated)
row.names(Diversity_matrix)<-unique(Diversity_df$Species)

#populate matrix
for (i in 1:nrow(Diversity_matrix)){
  ndf<-Diversity_df[which(Diversity_df$Species==row.names(Diversity_matrix)[i]), ]
  for (j in 1:ncol(Diversity_matrix)) {
    ndf1<-ndf[which(ndf$Associated==colnames(Diversity_matrix)[j]),]
    Diversity_matrix[i,j]<-ndf1$number[1]
  }
}
rm(ndf, ndf1, i, j)
#flip matrix so rows are substrates and columns are species
#values are #of indendent times each species is found on each substrate
Diversity_matrix<-t(Diversity_matrix)

#find unique times each substrate was sampled (is greater than number of unique samples
# or setIDs by 41. There are 41 setIDs that are assigned to more than one Associated substrate)
substrate_sampling<-data.frame(table(unique(store[c(1,2)])[1]))
colnames(substrate_sampling)[1]<-"Associated"
# removing the substrates that are sampled less than 5 times. 
# leaves the 15 top-sampled substrates
over20<-substrate_sampling$Associated[which(substrate_sampling$Freq>20)]
Diversity_matrix<-Diversity_matrix[which(
  row.names(Diversity_matrix) %in% over20),]
#store as seperate dataframe for downstream
Diversity_matrix->substrate_spp_accum_mat

###Species Diversity using Shannon-Weiner###
# find Shannon-Wiener index using Vegan
Shannon_Wiener_H<-data.frame(diversity(Diversity_matrix, index="shannon"))
Shannon_Wiener_H$Associated<-row.names(Shannon_Wiener_H)
Shannon_Wiener_H<-Shannon_Wiener_H[order(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..), ]

###Add sampling depth to dataframe 
Shannon_Wiener_H<-merge(Shannon_Wiener_H, substrate_sampling, by="Associated")


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


rm(substrate_sampling, store, over20, Diversity_df, Diversity_matrix)
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

###Basidio Associated analysis
Basidio_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Basidios)]
Basidio_Shannon_Wiener_H<-data.frame(diversity(Basidio_div_matrix, index="shannon"))
Basidio_Shannon_Wiener_H$IsoTemp<-row.names(Basidio_Shannon_Wiener_H)
Basidio_Shannon_Wiener_H<-Basidio_Shannon_Wiener_H[order(Basidio_Shannon_Wiener_H$IsoTemp), ]

write.table(Basidio_Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp_Basidios_only.tsv", sep="\t", quote=FALSE, row.name=FALSE)

###Pezizo Associated analysis
Pezizo_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Pezizos)]
Pezizo_Shannon_Wiener_H<-data.frame(diversity(Pezizo_div_matrix, index="shannon"))
Pezizo_Shannon_Wiener_H$IsoTemp<-row.names(Pezizo_Shannon_Wiener_H)
Pezizo_Shannon_Wiener_H<-Pezizo_Shannon_Wiener_H[order(Pezizo_Shannon_Wiener_H$IsoTemp), ]

write.table(Pezizo_Shannon_Wiener_H, "~/SurveyPaper/data/Shannon_Weiner_IsoTemp_Pezizos_only.tsv", sep="\t", quote=FALSE, row.name=FALSE)

###Saccharo Associated analysis
Saccharo_div_matrix<-Diversity_matrix[,which(colnames(Diversity_matrix) %in% Saccharos)]
Saccharo_Shannon_Wiener_H<-data.frame(diversity(Saccharo_div_matrix, index="shannon"))
Saccharo_Shannon_Wiener_H$IsoTemp<-row.names(Saccharo_Shannon_Wiener_H)
Saccharo_Shannon_Wiener_H<-Saccharo_Shannon_Wiener_H[order(Saccharo_Shannon_Wiener_H$IsoTemp), ]

rm(Basidio_div_matrix, Basidios, Saccharo_div_matrix, Saccharos, Pezizo_div_matrix, Pezizos, store)

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
par(mar=c(6,4,1,1))
#par(mfrow=c(2,1))
layout(matrix(c(1:2, 0, 0), nrow=1, ncol=2, byrow=TRUE), widths = c(2.5,2.5))
layout.show(n=2)  # to inspect layout                   
plot(rarefied, ci.type = "polygon", ci.lty = 0, ci.col = "grey", ylab="No. yeast OTUs", xlab="Samples")
rarecurve(phy_cum_spec_accum_mat,label = FALSE, col =color_pal, ylab = "No. yeast OTUs", xlab="No. isolates", lwd=2.5 )

quartz.save("~/SurveyPaper/Figures/Sup_rarefaction_Figure_1.pdf", type="pdf")


#H' plots with rarefaction
Shannon_Wiener_H<-Shannon_Wiener_H[order(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..), ]
SW_order<-c(13,8,15,7,14,9,11,6,4,2,10,3,5,12,1)
substrate_spp_accum_mat<-substrate_spp_accum_mat[SW_order, ]
s_s_a_m_rowsums<-rowSums(substrate_spp_accum_mat)
zoom_in<-substrate_spp_accum_mat[1:11, ]



#####
#Jaccard distance estimates for beta diversity
#####
#starting with substrate_spp_accum_mat which is not binary
substrate_distance<-as.matrix(
  vegdist(substrate_spp_accum_mat, method = "jaccard", binary = T))
#checked to make sure a binary matrix would give the same result. it did.
substrate_distance[lower.tri(substrate_distance)] <- NA
library(pheatmap)
pheatmap(substrate_distance, cluster_rows = FALSE,
         cluster_cols = FALSE, na_col="white", color=col_grad_2)

##Now temperature
isoTemp_spp_accum_mat_1<-isoTemp_spp_accum_mat#[c(4,1,2,3),]
temp_distance<-as.matrix(
  vegdist(isoTemp_spp_accum_mat_1, method="jaccard", binary=T))
write.table(temp_distance,"~/SurveyPaper/data/isoTemp_jaccard_distance.tsv", sep="t", quote=FALSE, row.names=FALSE)
temp_distance[lower.tri(temp_distance)]<-NA
pheatmap(temp_distance, cluster_rows = FALSE,
         cluster_cols = FALSE, na_col="white", color=col_grad_2)
quartz.save("jaccard_distance_heatmap_for_fig4.pdf", type="pdf")

####Figure 3 - substrate H' and rarefaction###

par(mar=c(4,4,1,1))
layout(matrix(c(1, 2, 1, 3), nrow=2, ncol=2, byrow=TRUE), widths = c(3, 2))
layout.show(n=3)  # to inspect layout                   

#1
plot(y=log(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..),
     x=log(Shannon_Wiener_H$Freq), col=col_grad[1:4], pch=19, cex=2, 
     ylab="log(Shannon-Weiner H')", xlab="log(Sampling density)")
text(y=log(Shannon_Wiener_H$diversity.Diversity_matrix..index....shannon..),
     x=log(Shannon_Wiener_H$Freq), labels=Shannon_Wiener_H$Associated,
     adj=c(.5,1.9), cex=.7)
#2
rarecurve(substrate_spp_accum_mat,label = FALSE, 
          ylab = "No. yeast OTUs", xlab="No. isolates", lwd=.5, col=col_grad)
#3
rarecurve(zoom_in, label=FALSE,
          ylab = "No. yeast OTUs", xlab="No. isolates", lwd=.5, col=col_grad[1:11])
quartz.save("~/SurveyPaper/Figures/updatedfig3_2-9-21.pdf", type="pdf")
######

######################
#Figure 4 - Isotemp H' and rarefaction
#####################

Shannon_Wiener_H_2<-Shannon_Wiener_H_2[order(as.numeric(Shannon_Wiener_H_2$IsoTemp)), ]
sw_isotemp_order<-c(4,1,2,3)
isoTemp_spp_accum_mat<-isoTemp_spp_accum_mat[sw_isotemp_order, ]

par(mar=c(4,4,1,1))
layout(matrix(c(1, 1, 0, 0,
                2, 3, 4, 5), nrow=2, byrow=TRUE), heights=c(2,1.25))
layout.show(n=5)  # to inspect layout                   

#1
rarecurve(isoTemp_spp_accum_mat, label=FALSE, 
          ylab = "No. yeast OTUs", xlab="No. isolates", lwd=1, col=col_grad[c(1,5,10,15)])
#2
plot(x=Shannon_Wiener_H_2$IsoTemp, y=Shannon_Wiener_H_2[,2], ylab="Shannon-Wiener (H')", 
     main=NA,cex=2, pch=16,ylim=c(1,4.5), xlim=c(0,32), xlab= "Isolation temperature", col=col_grad[c(1,5,10,15)])
text(x=4,y=4.5, labels="All taxa", cex=.65)
#3#Basdio
plot(x=Basidio_Shannon_Wiener_H$IsoTemp, y=Basidio_Shannon_Wiener_H[,1], ylab="Shannon-Wiener (H')", 
     main=NA,cex=2, ylim=c(1,4.5), xlim=c(0,32), pch=16, xlab= "Isolation temperature", col=color_pal[1])
text(x=7.75,y=4.5, labels="Basidiomycota", cex=.65)
#4#Pezizo
plot(x=Pezizo_Shannon_Wiener_H$IsoTemp, y=Pezizo_Shannon_Wiener_H[,1], ylab="Shannon-Wiener (H')", 
     main=NA,cex=2, pch=16,ylim=c(1,4.5),xlim=c(0,32), xlab= "Isolation temperature", col=color_pal[2])
text(x=8.25,y=4.5, labels="Pezizomycotina", cex=.65)
#5#Saccharo
plot(x=Saccharo_Shannon_Wiener_H$IsoTemp, y=Saccharo_Shannon_Wiener_H[,1], ylab="Shannon-Wiener (H')", 
     main=NA,cex=2, pch=16, ylim=c(1,4.5),xlim=c(0,32),xlab= "Isolation temperature", col=color_pal[3])
text(x=9.75,y=4.5, labels="Saccharomycotina", cex=.65)


quartz.save("~/SurveyPaper/Figures/updatedfig4_2-12-21.pdf", type="pdf")




