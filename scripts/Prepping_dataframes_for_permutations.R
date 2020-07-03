#Contains code used to prep datasets for each permutation analysis
#Contains descriptive figures for each dataset. 


require(ggplot2)

#Read in the dataframe and convert blank cells to NA
raw_WY_dataframe<-read.csv("~/WY_df_2018-02-08.csv",
                  header=TRUE, stringsAsFactors=FALSE,
                  na.strings=c("","NA"),
                  strip.white=TRUE)


#General descriptive figures
individual_isolates<-data.frame(table(raw_WY_dataframe$Species))
nrow(individual_isolates)
sum(individual_isolates$Freq)
length(which(individual_isolates$Freq>50))
a<-ggplot(individual_isolates, aes(Freq))+
  geom_histogram()+
  theme_bw()+
  ggtitle("262 unique species accross 1964 isolates", 
          subtitle = "104 singletons and 8 spp. isolated >50x")+
  xlab("No. times isolated")


set_isolates<-data.frame(table(unique(raw_WY_dataframe[c(22,29)])[,2]))
nrow(set_isolates)
sum(set_isolates$Freq)
length(which(set_isolates$Freq==1))
length(which(set_isolates$Freq>50))
b<-ggplot(set_isolates, aes(Freq))+
  geom_histogram()+
  theme_bw()+
  ggtitle("262 unique species accross 1520 unique isolations", 
          subtitle = "116 singletons and 6 spp. isolated >50x")+
  xlab("No. times isolated from independent isolations")

ggarrange(a, b, nrow=2)
quartz.save("~/SurveyPaper/Isolate_isolations_hist.pdf", type="pdf")

#phylum descriptive figures
phytbl<-unique(raw_WY_dataframe[c(22,29,30)])
phyls<-data.frame(table(phytbl$Subphylum))
a<-ggplot(phyls, aes(x=Var1, y=Freq))+
  geom_col()+
  theme_bw()+
  ylim(0, 1200)+
  ylab("No.")+
  xlab("")+
  geom_text(aes(label = Freq), vjust = -0.5)+
  ggtitle("subphylum representation accross\n1520 unique isolations")

phyls_div<-unique(phytbl[c(2,3)])
phyls_div1<-data.frame(table(phyls_div$Subphylum))
b<-ggplot(phyls_div1, aes(x=Var1, y=Freq))+
  geom_col()+
  theme_bw()+
  ylim(0, 200)+
  ylab("No.")+
  xlab("")+
  geom_text(aes(label = Freq), vjust = -0.5)+
  ggtitle("unique species isolated\nin each subphylum", )

phyls_temps<-unique(raw_WY_data[c(22,29,30,23)])
phyls_temps<-phyls_temps[which(!phyls_temps$IsoTemp=="Unknown"),]
phyls_temps$IsoTemp<-as.factor(phyls_temps$IsoTemp)
phyls_temps$IsoTemp = factor(phyls_temps$IsoTemp,levels(phyls_temps$IsoTemp)[c(4,1,2,3)])
c<-ggplot(phyls_temps, aes(x=IsoTemp, fill=Subphylum))+
  geom_bar()+
  ylab("No.")+
  theme_bw()

ggarrange(a,b,c, nrow=2, ncol = 2)

quartz.save("~/SurveyPaper/Subphylum_descriptive_figure.pdf", type="pdf")
####### Association test 1 - plant genus - yeast spp. associations
#adding a step to retain only unique SetIDs
plant_genus_WY_dataframe<-unique(raw_WY_dataframe[c(15,22,29)])


#seperate out strains that have a plant genus associated with them
#I'm using the which function and a boolean operator.  
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!is.na(plant_genus_WY_dataframe$Plant.Genus)), ]
#
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Unknown"), ]
#returns just part of dataframe where Plant.Genus does NOT equal Uknown
#descriptive figure
plants<-length(unique(plant_genus_WY_dataframe$Plant.Genus))
yeast<-length(unique(plant_genus_WY_dataframe$Species))
pairs<-data.frame(table(plant_genus_WY_dataframe$Plant.Genus, plant_genus_WY_dataframe$Species))
comb<-length(which(pairs$Freq>0))
for_descrip_fig<-data.frame(c("Yeast sp.", "Plant genera", "Combinations observed"),
                            c(yeast, plants, comb))
colnames(for_descrip_fig)<-c("V1", "V2")

a<-ggplot(for_descrip_fig, aes(x=V1, y=V2))+
  geom_col(fill="lightgreen")+
  theme_bw()+
  xlab("")+
  ylim(0, 800)+
  ylab("No.")+
  geom_text(aes(label = V2), vjust = -0.5)+
  ggtitle("Raw data for plant genus association analysis")



rm(plants, yeast, pairs, comb)



#plant_genus_WY_dataframe gets fed into permutation script: 
#plant_genus_associations<-Run_WY_association_permutations(all_observations_dataframe=plant_genus_WY_dataframe,
 #                       permutations=10000, colnames_to_permute=c("Plant.Genus", "Species"))



####### Association test 1 - isolation temperature - yeast spp. associations
#adding a step to retain only unique SetIDs
iso_temp_WY_dataframe<-unique(raw_WY_dataframe[c(22,23,29)])

#
iso_temp_WY_dataframe<-iso_temp_WY_dataframe[which(!iso_temp_WY_dataframe$IsoTemp=="Unknown"), ]
#returns just part of dataframe where IsoTemp does NOT equal Uknown

#iso_temp_WY_dataframe gets fed into permutation script: 
#iso_temp_associations<-Run_WY_association_permutations(all_observations_dataframe=iso_temp_WY_dataframe,
 #                                                         permutations=10000, colnames_to_permute=c("IsoTemp", "Species"))

#descriptive figure
isotemp<-length(unique(iso_temp_WY_dataframe$IsoTemp))
yeast<-length(unique(iso_temp_WY_dataframe$Species))
pairs<-data.frame(table(iso_temp_WY_dataframe$IsoTemp, iso_temp_WY_dataframe$Species))
comb<-length(which(pairs$Freq>0))
for_descrip_fig<-data.frame(c("Yeast sp.", "Isolation temps", "Combinations observed"),
                            c(yeast, isotemp, comb))
colnames(for_descrip_fig)<-c("V1", "V2")

b<-ggplot(for_descrip_fig, aes(x=V1, y=V2))+
  geom_col(fill="skyblue")+
  theme_bw()+
  xlab("")+
  ylab("No.")+
  ylim(0, 800)+
  geom_text(aes(label = V2), vjust = -0.5)+
  ggtitle("Raw data for isolation temperature association analysis")

rm(isotemp, yeast, pairs, comb)

####### Association test 3 - substrate - yeast spp. associations
#adding a step to retain only unique SetIDs
Substrate_specific_WY_dataframe<-unique(raw_WY_dataframe[c(22,10,29)])

#
Substrate_specific_WY_dataframe<-Substrate_specific_WY_dataframe[which(!is.na(Substrate_specific_WY_dataframe$Specific)),]


#descriptive figure
substr<-length(unique(Substrate_specific_WY_dataframe$Specific))
yeast<-length(unique(Substrate_specific_WY_dataframe$Species))
pairs<-data.frame(table(Substrate_specific_WY_dataframe$Specific, Substrate_specific_WY_dataframe$Species))
comb<-length(which(pairs$Freq>0))
for_descrip_fig<-data.frame(c("Yeast sp.", "Specific substrates", "Combinations observed"),
                            c(yeast, substr, comb))
colnames(for_descrip_fig)<-c("V1", "V2")

c<-ggplot(for_descrip_fig, aes(x=V1, y=V2))+
  geom_col(fill="brown")+
  theme_bw()+
  ylim(0, 800)+
  xlab("")+
  ylab("No.")+
  geom_text(aes(label = V2), vjust = -0.5)+
  ggtitle("Raw data for substrate association analysis")

rm(substr, yeast, pairs, comb)
# 
#Substrate_Specific_Associations<-Run_WY_association_permutations(all_observations_dataframe=Substrate_specific_WY_dataframe,
 #                                                      permutations=10000, colnames_to_permute=c("Specific", "Species"))
#write.table(Substrate_Specific_Associations, "~/SurveyPaper/data/Substrate_Specific_Associations_10k_permutations.tsv", sep="\t",
 #           quote=FALSE, row.names=FALSE)


ggarrange(a, b, c, ncol=1)
quartz.save("~/SurveyPaper/Descriptive_figure_for_raw_data_fed_into_permutations.pdf", type="pdf")

