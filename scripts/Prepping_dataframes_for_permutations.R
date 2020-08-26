#Contains code used to prep datasets for each permutation analysis
#Contains descriptive figures for each dataset. 


require(ggplot2)
require(ggpubr)
require(RColorBrewer)


#Read in the dataframe and convert blank cells to NA
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
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
  ggtitle("262 unique OTUs accross 1962 isolates")+
  xlab("No. times isolated")


set_isolates<-data.frame(table(unique(raw_WY_dataframe[c(22,29)])[,2]))
nrow(set_isolates)
sum(set_isolates$Freq)
length(which(set_isolates$Freq==1))
length(which(set_isolates$Freq>50))
b<-ggplot(set_isolates, aes(Freq))+
  geom_histogram()+
  theme_bw()+
  ggtitle("262 unique OTUs accross 1518 unique isolations", 
          subtitle = "116 singletons and 6 OTUs isolated >50x")+
  xlab("No. times isolated from independent isolations")

ggarrange(a, b, nrow=2)
quartz.save("~/SurveyPaper/Figures/Isolate_isolations_hist.pdf", type="pdf")


####Singleton descriptive figure
singletons_df<-read.delim("~/SurveyPaper/data/Singletons.tsv", header=TRUE,
                       stringsAsFactors=FALSE)
toplot<-unique(singletons_df[c(29,30,31)])
toplot$significant<-NA
toplot$significant[3]<-"*"
a<-ggplot(toplot, aes(x=Subphylum))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  ylim(0, 70)+
  geom_text(stat='count', aes(label=..count..), vjust=-.5, size=2)+
  geom_text(aes(y=65, label=significant), size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Singletons\nby subphylum")

toplot<-unique(singletons_df[c(22,29,10)])
b<-ggplot(toplot, aes(x=Specific))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  ylim(0, 40)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_text(stat='count', aes(label=..count..), vjust=-.5, size=2)+
  ggtitle("Singletons by substrate")

#run mapping script lines 25:31
#need to adjust long and lat values to be plotted with us_map package using us_transform
mappable<-singletons_df[which(!is.na(singletons_df$Long)),]
temp<-usmap_transform(mappable[c(7,6)])
isolates_to_plot<-merge(mappable, temp, by=c("Long", "Lat"))
isolates_to_plot<-merge(isolates_to_plot, regions[c(2,3)], by="State")

p0 <- ggplot(data = us_states,
             mapping = aes(x =long, y = lat,
                           group = group, fill = Region,
                           alpha=.9))
p1<- p0 + geom_polygon(color = "gray90", size = 0.1)+
  scale_fill_manual(values = color_key$col)+
  theme_bw()+
  geom_jitter(data=isolates_to_plot, aes(x=Long.1, y=Lat.1), 
              inherit.aes = FALSE, width=.5, shape=1)+ 
  theme(legend.position = "none")+
  coord_equal()

ggarrange(b, a, p1, widths=c(4,1))
quartz.save("~/SurveyPaper/Figures/Singleton_descriptive_figure.pdf", type="pdf")

#cosmopolitan descriptive figure/table #Updated 08-24-20 to reflect new set of cosmo spp.
# definition - those species that were isolated when they were expected 2 or more times and 
#not isolated when expected 1 or less times. 

regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                    header=TRUE, stringsAsFactors=FALSE)
cosmodf<-read.delim("~/SurveyPaper/data/Cosmopolitans_foundexp>1_notfoundexp<=1.tsv",
                    header=TRUE, stringsAsFactors=FALSE)

toplot<-unique(cosmodf[c(29,30)])
a<-ggplot(toplot, aes(x=Subphylum))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  ylim(0, 10)+
  geom_text(stat='count', aes(label=..count..), vjust=-.5, size=2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Cosmopolitan\nby subphylum")

toplot<-unique(cosmodf[c(22,29,10)])
toplot$pvalue<-NA
toplot$pvalue[which(toplot$Specific=="Soil")]<-"***"
b<-ggplot(toplot, aes(x=Specific))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  geom_text(aes(x=Specific, label= pvalue, y=155), size=2)+
  ylim(0, 175)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_text(stat='count', aes(label=..count..), vjust=-.5, size=2)+
  ggtitle("Cosmopolitan by substrate")

#run mapping script lines 25:32
mappable<-cosmodf[which(!is.na(cosmodf$Long)),]
temp<-usmap_transform(mappable[c(7,6)])
isolates_to_plot<-merge(mappable, temp, by=c("Long", "Lat"))
isolates_to_plot<-merge(isolates_to_plot, regions[c(2,3)], by="State")


p0 <- ggplot(data = us_states,
             mapping = aes(x = long, y = lat,
                           group = group, fill = Region,
                           alpha=.9), show.legend = FALSE)
p1<-  p0 + geom_polygon(color = "gray90", size = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = color_key$col)+
  scale_color_manual(values = brewer.pal(name = "Paired", n=11))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_jitter(data=isolates_to_plot, aes(x=Long.1, y=Lat.1, col=Species), 
              inherit.aes = FALSE, width=.5, shape=1)+ 
  ggtitle("Cosmopolitan isolation locations")+
coord_equal()

p1
quartz.save("~/SurveyPaper/Figures/Cosmopolitan_map.pdf", type="pdf")


ggarrange(b, a, p1, widths=c(4,1))
quartz.save("~/SurveyPaper/Figures/Cosmopolitan_descriptive_figure.pdf", type="pdf")


#regionally restricted descriptive figures
regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                    header=TRUE, stringsAsFactors=FALSE)
regional<-read.delim("~/SurveyPaper/data/Regionally_restricted_foundexp=1_notfoundexp>2.tsv",
                    header=TRUE, stringsAsFactors=FALSE)
#run mapping script lines 25:32
mappable<-regional[which(!is.na(regional$Long)),]
temp<-usmap_transform(mappable[c(7,6)])
isolates_to_plot<-merge(mappable, temp, by=c("Long", "Lat"))
isolates_to_plot<-merge(isolates_to_plot, regions[c(2,3)], by="State")

toplot<-unique(regional[c(29,30)])
a<-ggplot(toplot, aes(x=Subphylum))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  ylim(0, 15)+
  geom_text(stat='count', aes(label=..count..), vjust=-.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Regionally restricted\nby subphylum")

toplot<-unique(regional[c(22,29,10)])
b<-ggplot(toplot, aes(x=Specific))+
  geom_bar()+
  ylab("No.")+
  xlab("")+
  #geom_text(aes(x=Specific, label= pvalue, y=155), size=2)+
  ylim(0, 20)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_text(stat='count', aes(label=..count..), vjust=-.5, size=2)+
  ggtitle("Regionally restricted by substrate")

#run mapping script lines 7:28
p0 <- ggplot(data = us_states,
             mapping = aes(x = long, y = lat,
                           group = group, fill = Region,
                           alpha=.9), show.legend = FALSE)
p1<-  p0 + geom_polygon(color = "gray90", size = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = color_key$col)+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(16))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_jitter(data=isolates_to_plot, aes(x=Long.1, y=Lat.1, col=Species), 
              inherit.aes = FALSE, width=.5, shape=1)+ 
  ggtitle("Regionally restricted isolation locations")+
  coord_equal()

p1
quartz.save("~/SurveyPaper/Figures/Regional_map.pdf", type="pdf")


ggarrange(b, a, p1, widths=c(4,1), heights = c(1.5,2))
quartz.save("~/SurveyPaper/Figures/Regional_descriptive_figure.pdf", type="pdf")





#phylum descriptive figures
phytbl<-unique(raw_WY_dataframe[c(22,29,30)])
phyls<-data.frame(table(phytbl$Subphylum))
colnames(phyls)[2]<-"Unique.isolates"
phyls_div<-unique(phytbl[c(2,3)])
phyls_div1<-data.frame(table(phyls_div$Subphylum))
colnames(phyls_div1)[2]<-"Unique.OTUs"
phyls_div1<-merge(phyls_div1, phyls, by="Var1")
b<-ggplot(phyls_div1, aes(x=Var1, y=Unique.isolates))+
  geom_col(fill="grey90")+
  geom_text(aes(label = Unique.isolates ), vjust = -2, size=2.5)+
  geom_col(aes(x=Var1, y=Unique.OTUs), fill="grey43")+
  geom_text(aes(x=Var1, y=Unique.OTUs, label=Unique.OTUs), vjust=-0.5, size=2.5)+
  theme_bw()+
  ylim(0, 1200)+
  ylab("No.")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

quartz.save("~/SurveyPaper/Figures/Subphylum_descriptive_figure.pdf", type="pdf")


phyls_temps<-unique(raw_WY_dataframe[c(22,29,30,23)])
phyls_temps<-phyls_temps[which(!phyls_temps$IsoTemp=="Unknown"),]
phyls_temps$IsoTemp<-as.factor(phyls_temps$IsoTemp)
phyls_temps$IsoTemp = factor(phyls_temps$IsoTemp,levels(phyls_temps$IsoTemp)[c(4,1,2,3)])
c<-ggplot(phyls_temps, aes(x=IsoTemp, fill=Subphylum))+
  geom_bar()+
  ylab("No.")+
  theme_bw()

quartz.save("~/SurveyPaper/Figures/Temperature_by_subphylum.pdf", type="pdf")


####### Association test 1 - plant genus - yeast spp. associations
#adding a step to retain only unique SetIDs
plant_genus_WY_dataframe<-unique(raw_WY_dataframe[c(15,22,29)])


#seperate out strains that have a plant genus associated with them
#I'm using the which function and a boolean operator.  
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!is.na(plant_genus_WY_dataframe$Plant.Genus)), ]
#
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Unknown"), ]
#returns just part of dataframe where Plant.Genus does NOT equal Uknown
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Uknown"), ]

####
#descriptive figure of plant genus distribution
plant_genus_WY_dataframe<-unique(raw_WY_dataframe[c(15,22,29,10)])
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!is.na(plant_genus_WY_dataframe$Plant.Genus)), ]
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Unknown"), ]
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Uknown"), ]
toplot<-data.frame(table(plant_genus_WY_dataframe$Plant.Genus))
toplot<-toplot[order(toplot$Freq, decreasing = TRUE),]
toplot$order<-c(1:nrow(toplot))
colnames(toplot)<-c("Plant.Genus", "count", "order")
plant_genus_WY_dataframe<-merge(toplot, plant_genus_WY_dataframe, by="Plant.Genus")
plant_genus_WY_dataframe$Plant.Genus<-reorder(plant_genus_WY_dataframe$Plant.Genus, plant_genus_WY_dataframe$order)
colnames(plant_genus_WY_dataframe)[6]<-"Substrate.category"
S_subsB<-ggplot(plant_genus_WY_dataframe, aes(x=Plant.Genus, fill=Substrate.category))+
  geom_bar()+
  theme_bw()+
  ylim(0, 200)+
  xlab("")+
  ylab("No.")+
  geom_text(data=toplot, aes(x=Plant.Genus, y=count, label = count), inherit.aes = FALSE, vjust = -0.5, size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(legend.key.size = unit(0.075, "cm"), legend.key.width = unit(.10, "cm"),
        legend.text=element_text(size=5))

ggarrange(nrow = 2, S_subsA, S_subsB)
quartz.save("~/SurveyPaper/Figures/Descriptive_fig_of_substr_plantgenus_categories.pdf", type="pdf")

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
iso_temp_WY_dataframe<-unique(raw_WY_dataframe[c(22,23,29,30)])

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

#descriptive figure about substrate variable distribution
toplot<-data.frame(table(Substrate_specific_WY_dataframe$Specific))
toplot<-toplot[order(toplot$Freq, decreasing = TRUE),]
toplot$order<-c(1:nrow(toplot))
toplot$Var1<-reorder(toplot$Var1, toplot$order)
S_subsA<-ggplot(toplot, aes(x=Var1, y=Freq))+
  geom_col()+
  theme_bw()+
  ylim(0, 550)+
  xlab("")+
  ylab("No.")+
  geom_text(aes(label = Freq), vjust = -0.5, size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
quartz.save("~/SurveyPaper/Figures/Substrate_specific_category_descriptive_fig.pdf", type="pdf")  



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


ggarrange(c, a, b, ncol=1)
quartz.save("~/SurveyPaper/Figures/Descriptive_figure_for_raw_data_fed_into_permutations.pdf", type="pdf")



##########
#####
#IDing possible "Endemic spp" as those isolated z or more independent times from the same region
z=3
ind_iso_sp_by_state<-unique(raw_WY_dataframe[c(22, 29, 25)])

climate_regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                            stringsAsFactors=FALSE)
climate_regions<-climate_regions[c(2,3)]
ind_iso_sp_by_state<-merge(ind_iso_sp_by_state, climate_regions, by="State")
tbl<-data.frame(table(ind_iso_sp_by_state$Species, ind_iso_sp_by_state$Region))
tbl<-tbl[which(tbl$Freq>0),]
colnames(tbl)<-c("Species", "Region", "Iso_in_region")
tbl2<-data.frame(table(tbl$Species))
colnames(tbl2)<-c("Species", "Number_of_regions")
tbl3<-merge(tbl, tbl2, by="Species")
Endemics<-tbl3[which(tbl3$Iso_in_region>=z & tbl3$Number_of_regions==1),]


probsRegions<-data.frame(table(ind_iso_sp_by_state$Region))
colnames(probsRegions)<-c("Region", "No. independent isolations")
independent_isolations_by_region<-merge

probsRegions$Prob<- probsRegions$Freq/sum(probsRegions$Freq)

