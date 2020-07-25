raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                           header=TRUE, stringsAsFactors=FALSE,
                           strip.white=TRUE, na.strings=c("","NA"))

Species_by_location_DF<-raw_WY_dataframe[c(1,25,29)]
Species_by_location_DF<-Species_by_location_DF[which(!is.na(Species_by_location_DF$State)),]
#1958 isolates
locations_by_sp<-data.frame(table(Species_by_location_DF$Species,
                                      Species_by_location_DF$State))
locations_by_sp<-locations_by_sp[which(locations_by_sp$Freq>0),]
colnames(locations_by_sp)<-c("Species", "State", "Freq")
#621 unique state x species combinations
#360 spp. isolated from more than one state
climate_regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
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
require(plot3D)
n <- 10
col_vector = brewer.pal(n, "Paired")
regional_color_key<-data.frame(unique(climate_regions$Region), col_vector)
colnames(regional_color_key)<-c("Region", "col")
write.table(regional_color_key, "~/SurveyPaper/data/Geographical_analysis/regional_color_key.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
par(mar=c(1,1,1,1))
pie(rep(1, length(col_vector)), 
    labels=regional_color_key$Region, col=as.character(regional_color_key$col))
#quartz.save("~/SurveyPaper/data/Geographical_analysis/Regional_color_key.pdf", type="pdf")

###
#Figured out how to pull all ITS out of labDB: 

all_ITS<-read.delim("~/SurveyPaper/data/tables_for_scripts/All_ITS_seq.txt",header=TRUE, stringsAsFactors=FALSE)

StrainIDs<-raw_WY_dataframe[c(1,29)]
spp_for_regional_analysis<-merge(spp_for_regional_analysis, StrainIDs, by="Species")
spp_for_regional_analysis<-merge(spp_for_regional_analysis, all_ITS, by="StrainID")
spp_for_regional_analysis<-spp_for_regional_analysis[c(2,3,4,1,5)]

spp_for_regional_analysis<-spp_for_regional_analysis[order(
  spp_for_regional_analysis$No.Regions, decreasing = TRUE, spp_for_regional_analysis$Species),]

write.table(spp_for_regional_analysis, "~/SurveyPaper/data/tables_for_scripts/Species_and_ITS_seq_for_regional_analyses.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

