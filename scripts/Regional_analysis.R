

Species_by_location_DF<-raw_WY_dataframe[c(1,25,29)]
Species_by_location_DF<-Species_by_location_DF[which(!is.na(Species_by_location_DF$State)),]
#1958 isolates
locations_by_sp<-data.frame(table(Species_by_location_DF$Species,
                                      Species_by_location_DF$State))
locations_by_sp<-locations_by_sp[which(locations_by_sp$Freq>0),]
colnames(locations_by_sp)<-c("Species", "State", "Freq")
#621 unique state x species combinations
#360 spp. isolated from more than one state
climate_regions<-read.delim("~/SurveyPaper/data/NOAA_US_Climate_Regions.txt",
                            stringsAsFactors=FALSE)
climate_regions<-climate_regions[c(2,3)]
regions_by_sp<-merge(Species_by_location_DF, climate_regions, by="State")
tbl<-data.frame(table(regions_by_sp$Species, regions_by_sp$Region))
tbl<-tbl[which(tbl$Freq>0),]
spp_accross_multiple_regions<-data.frame(table(tbl$Var1))
colnames(spp_accross_multiple_regions)<-c("Species", "No.Regions")
require(ggplot2)
ggplot(spp_accross_multiple_regions, aes(No.Regions))+
  geom_histogram(bins = 7, fill="lightpink", col="lightgrey")+
  theme_bw()+
  geom_text(stat='count', aes(label=..count..),
            position = position_stack(vjust = 0.5),size=4)+
  xlab("No. of regions species are isolated from")+
  ggtitle("Number of regions yeast spp. are isolated from")
quartz.save("~/SurveyPaper/data/Geographical_analysis/Histogram_of_spp_by_number_of_regions.pdf", type="pdf")

#Need at least 2 isolates in at least 2 different regions for analysis
tbl<-tbl[which(tbl$Freq>1),]
spp_accross_multiple_regions_2plus<-data.frame(table(tbl$Var1))
colnames(spp_accross_multiple_regions_2plus)<-c("Species", "No.Regions")
spp_accross_multiple_regions_2plus<-spp_accross_multiple_regions_2plus[
  which(spp_accross_multiple_regions_2plus$No.Regions>0),]

ggplot(spp_accross_multiple_regions_2plus, aes(No.Regions))+
  geom_histogram(bins=7, fill="lightpink", col="lightgrey")+
  theme_bw()+
  geom_text(stat='count', aes(label=..count..),
            position = position_stack(vjust = 0.5),size=4)+
  xlab("No. of regions species are isolated from")+
  ggtitle("Number of regions yeast spp. isolated\n at least twice in a region are isolated from")
quartz.save("~/SurveyPaper/data/Geographical_analysis/Histogram_of_spp_w/greaterthan2isolatesperregion_by_number_of_regions.pdf", type="pdf")

#52 spp. with potential analysis

spp_for_regional_analysis<-spp_accross_multiple_regions_2plus[which(spp_accross_multiple_regions_2plus$No.Regions>1),]

### Lachancea
Lachancea_spp_for_analysis<-spp_for_regional_analysis[which(grepl("Lachancea",spp_for_regional_analysis$Species)),]
Lachancea_spp_for_analysis<-merge(Lachancea_spp_for_analysis,
                                  Species_by_location_DF, by="Species")
Lachancea_spp_for_analysis<-merge(Lachancea_spp_for_analysis, climate_regions,
                                  by="State")
write.table(Lachancea_spp_for_analysis, "~/SurveyPaper/data/Geographical_analysis/Lachancea_spp.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
### Hanseniaspora
Hanseniaspora_spp_for_analysis<-spp_for_regional_analysis[which(grepl("Hanseniaspora",spp_for_regional_analysis$Species)),]
Hanseniaspora_spp_for_analysis<-merge(Hanseniaspora_spp_for_analysis,
                                  Species_by_location_DF, by="Species")
Hanseniaspora_spp_for_analysis<-merge(Hanseniaspora_spp_for_analysis, climate_regions,
                                  by="State")
write.table(Hanseniaspora_spp_for_analysis, "~/SurveyPaper/data/Geographical_analysis/Hanseniaspora_spp.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
### Kluyveromyces
Kluyveromyces_spp_for_analysis<-spp_for_regional_analysis[which(grepl("Kluyveromyces",spp_for_regional_analysis$Species)),]
Kluyveromyces_spp_for_analysis<-merge(Kluyveromyces_spp_for_analysis,
                                      Species_by_location_DF, by="Species")
Kluyveromyces_spp_for_analysis<-merge(Kluyveromyces_spp_for_analysis, climate_regions,
                                      by="State")
write.table(Kluyveromyces_spp_for_analysis, "~/SurveyPaper/data/Geographical_analysis/Kluyveromyces_spp.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

####tree beginnings
#color pallette
library(RColorBrewer)
n <- length(unique(climate_regions$Region))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(2423)
col=sample(col_vector, n)
pie(rep(1,n), col=col)
regional_color_key<-data.frame(unique(climate_regions$Region), col)
colnames(regional_color_key)<-c("Region", "col")
write.table(regional_color_key, "~/SurveyPaper/data/Geographical_analysis/regional_color_key.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
###

require(tidyr)
require(ape)
tree<-read.tree("~/Downloads/phylo.io.nwk")
tips<-data.frame(tree$tip.label)
tips<-tips %>% separate(tree.tip.label, c("Gen", "Sp", "StrainID"))
tips$order<-c(1:nrow(tips))
tips2<-merge(tips, regions_by_sp[c(2,4)], by="StrainID")
tips2<-merge(tips2, regional_color_key, by="Region")
tips2$col<-as.character(tips2$col)
tips2<-tips2[order(tips2$order),]

plot(tree, edge.color = tips2$col)
