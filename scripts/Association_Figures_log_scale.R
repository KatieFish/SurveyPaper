color_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


phyls_temps<-unique(raw_WY_dataframe[c(22,29,30,23)])
phyls_temps<-phyls_temps[which(!phyls_temps$Subphylum=="Taphrinomycotina"),]
phyls_temps<-phyls_temps[which(!phyls_temps$IsoTemp=="Unknown"),]
phyls_temps$IsoTemp<-as.factor(phyls_temps$IsoTemp)
phyls_temps$IsoTemp = factor(phyls_temps$IsoTemp,levels(phyls_temps$IsoTemp)[c(4,1,2,3)])
c<-ggplot(phyls_temps, aes(x=IsoTemp, fill=Subphylum))+
  geom_bar()+
  scale_fill_manual(values=color_pal)+
  ylab("No.")+
  ggtitle("")+
  theme_bw()

require(reshape2)

Plant_associations<-read.delim("~/SurveyPaper/data/Plant_Genus_Association_results.tsv", header=TRUE)

significant<-Plant_associations[which(Plant_associations$BH_Padj<=0.05),c(1,2,4,7,5)]

raw_WY_data<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                        header=TRUE, stringsAsFactors=FALSE,
                      strip.white=TRUE)

plant_genus_dataframe<-unique(raw_WY_data[c(15,22,29)])
yeast_freq<-data.frame(table(plant_genus_dataframe$Species))
colnames(yeast_freq)[1]<-"Yeast_sp"
significant<-merge(significant, yeast_freq, by="Yeast_sp")
colnames(significant)[6]<-"OTU frequency"

plant_freq<-data.frame(table(plant_genus_dataframe$Plant.Genus))
colnames(plant_freq)[1]<-"Plant_genus"
significant<-merge(significant, plant_freq, by="Plant_genus")
colnames(significant)[7]<-"Substrate genus\nfrequency"

#write.table(significant, "~/SurveyPaper/data/Significant_BHadj_0_05_Plant_genus_Yeast_sp_associations.tsv",
            #sep="\t", quote=FALSE, row.names=FALSE)

significant$pair<-paste(significant$Yeast_sp, significant$Plant_genus, sep="-")
significant<-significant[order(significant$Observed, decreasing=TRUE), ]
significant$order<-c(1:nrow(significant))
significant$pair<-reorder(significant$pair, significant$order)
toplot<-significant[c(8, 3, 6, 5, 7)]
xP<-melt(toplot)

xP$variable <- factor(xP$variable, levels = c("Expected", "Observed", "OTU frequency", "Substrate genus\nfrequency"))

genus<-ggplot(xP, aes(x=pair, y=value, col=variable, shape=variable))+
  geom_jitter(size=3, stroke=1, width=.10)+
  xlab("")+
  ylab("Frequency observed in dataset (log10)")+
  scale_color_manual(values=color_pal[4:7])+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,face = "italic"))+
  scale_y_continuous(trans= 'log10')
  

quartz.save("~/SurveyPaper/Figures/BHadj_log_Significant_observed_frequencies_plant_associations.pdf", type="pdf")

###########Temperature

temperature_associations<-read.delim("~/SurveyPaper/data/Temperature_Associations_results.tsv",
                                     header=TRUE)

significant<-temperature_associations[which(temperature_associations$BH_Padj<=0.05),c(1,2,4,7,5)]

temp_sp_dataframe<-unique(raw_WY_data[c(22,23,29)])

yeast_freq<-data.frame(table(temp_sp_dataframe$Species))
colnames(yeast_freq)<-c("Yeast_sp", "OTU frequency")
significant<-merge(significant, yeast_freq, by="Yeast_sp")
significant$pair<-paste(significant$Yeast_sp, significant$IsoTemp, sep="-")
significant<-significant[order(significant$IsoTemp, significant$Observed, decreasing = TRUE), ]
significant$order<-c(1:nrow(significant))
significant$pair<-reorder(significant$pair, significant$order)

toplot<-significant[c(7,3,6,5)]
require(reshape2)
xT<-melt(toplot)
merge(xT, significant[c(1, 2, 7)], by="pair")->xT
xT$variable <- factor(xT$variable, levels = c("Expected", "Observed", "OTU frequency"))

temperature<-ggplot(xT, aes(x=pair, y=value, col=variable, shape=variable))+
  geom_jitter(size=3, stroke=1, width=.1)+
  xlab("")+
  ylab("Frequency observed in dataset (log10)")+
  scale_color_manual(values=color_pal[4:7])+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(
  labels=c(expression(paste(italic("Metschnikowia pulcherrima")," sp. complex-30")), 
  expression(paste(italic("Lachancea kluyveri"),"-30")), 
  expression(paste(italic("Saccharomyces cerevisiae"),"-30")),
  expression(paste(italic("Kluyveromyces marxianus"),"-30")), 
  expression(paste(italic("Candida pseudolambica"),"-30")), 
  expression(paste(italic("Lecythophora")," sp.-30")), 
  expression(paste(italic("Mrakia")," sp.-10")), 
  expression(paste(italic("Cystofilobasidium capitatum"),"-10")), 
  expression(paste(italic("Mrakia gelida"),"-10")), 
  expression(paste(italic("Candida sake"),"-10")), 
  expression(paste(italic("Rhodotorula fujisanensis"),"-10")), 
  expression(paste(italic("Mrakia")," sp.-4")), 
  expression(paste(italic("Mrakia gelida"),"-4")), 
  expression(paste(italic("Cystofilobasidium capitatum"),"-4")), 
  expression(paste(italic("Mrakia blollopis"),"-4")), 
  expression(paste(italic("Mrakiella cryoconiti"),"-4"))))+
  scale_y_continuous(trans= 'log10')

quartz.save("~/SurveyPaper/Figures/BHadj_log_Significant_observed_frequencies_isotemp_associations.pdf", type="pdf")


#write.table(significant, "~/SurveyPaper/data/Significant_BHadj_0_1_Isotemp_Yeast_sp_associations.tsv", 
           # sep="\t", quote=FALSE, row.name=FALSE)
  
#########
require(gridExtra)
require(ggpubr)

ggarrange(a,b, nrow = 1)
########

#Specific Substrate

Specific_associations<-read.delim("~/SurveyPaper/data/Substrate_Specific_Associations_results.tsv", header=TRUE)

significant<-Specific_associations[which(Specific_associations$BH_Padj<=.05), c(1,2,4,5,7)]

specific_sp_dataframe<-unique(raw_WY_data[c(22,10,29)])

yeast_freq<-data.frame(table(specific_sp_dataframe$Species))
colnames(yeast_freq)<-c("Species", "OTU frequency")
significant<-merge(significant, yeast_freq, by="Species")
significant$pair<-paste(significant$Species, significant$Specific, sep="-")
significant<-significant[order(significant$Observed, decreasing=TRUE), ]
significant$order<-c(1:nrow(significant))
significant$pair<-reorder(significant$pair, significant$order)


Specific_freq<-data.frame(table(specific_sp_dataframe$Specific))
colnames(Specific_freq)<-c("Specific", "Specific substrate\nfrequency")
significant<-merge(significant, Specific_freq, by="Specific")


toplot<-significant[c(7,9,3,4,6)]
require(reshape2)
xT<-melt(toplot)
xT$variable = factor(xT$variable,levels(xT$variable)[c(4,1,3,2)])

xT$variable <- factor(xT$variable, levels = c("Expected", "Observed", "OTU frequency", "Specific substrate\nfrequency"))


substrate<-ggplot(xT, aes(x=pair, y=value, col=variable, shape=variable))+
  geom_jitter(size=3,stroke=1, width=.1)+
  xlab("")+
  ylab("Frequency observed in dataset (log10)")+
  scale_color_manual(values=color_pal[4:7])+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels=c(expression(paste(italic("Torulaspora delbrueckii"), "-Soil")), 
                            expression(paste(italic("Mrakia"), " sp.-Soil")), 
                            expression(paste(italic("Saccharomyces paradoxus"), "-Soil")), 
                            expression(paste(italic("Cyberlindnera saturnus"), "-Soil")), 
                            expression(paste(italic("Hanseniaspora uvarum"), "-Fungus")), 
                            expression(paste(italic("Mrakia gelida"), "-Leaves")), 
                            expression(paste(italic("Lachancea kluyveri"), "-Bark")), 
                            expression(paste(italic("Trichosporon porosum"), "-Soil")), 
                            expression(paste(italic("Scheffersomyces ergatensis"), "-Bark")), 
                            expression(paste(italic("Curvibasidium cygneicollum"), "-Fruit")), 
                            expression(paste(italic("Rhodotorula nothofagi"), "-Leaves")), 
                            expression(paste(italic("Suhomyces bolitotheri"), "-Fungus")), 
                            expression(paste(italic("Teunomyces cretensis/kruisii"), "complex-Fungus")), 
                            expression(paste(italic("Kazachstania serrabonitensis"), "-Sand")), 
                            expression(paste(italic("Pichia kudriavzevii"), "-Fruit")), 
                            expression(paste(italic("Candida mycetangii"), "-Plant Matter")), 
                            expression(paste(italic("Pichia scaptomyzae"), "-Sand")), 
                            expression(paste(italic("Scleroconidioma sphagnicola"), "-Lichen")), 
                            expression(paste(italic("Sydowia polyspora"), "-Needles")), 
                            expression(paste(italic("Zygowilliopsis californica"), "-Flower")), 
                            expression(paste(italic("Kwoniella newhampshirensis"), "-Insect")), 
                            expression(paste(italic("Peterozyma toletana"), "-Feather"))))+
  scale_y_continuous(trans= 'log10')

  
  
  quartz.save("~/SurveyPaper/Figures/BH-adj_log_Significant_observed_frequencies_Specific_associations.pdf", type="pdf")


addSmallLegend(substrate, pointSize = 1, textSize = 8, spaceLegend = 0.1)->substrate
addSmallLegend(genus, pointSize = 1, textSize = 8, spaceLegend = 0.1)->genus
ggarrange(substrate, genus, nrow=1, labels = c("A", "B"))
quartz.save("~/SurveyPaper/Manuscript_main_figs/ToBeFig1.pdf", type="pdf")
addSmallLegend(temperature, pointSize = 1, textSize = 8, spaceLegend = 0.1)->temperature
addSmallLegend(c, textSize = 9, spaceLegend = 0.3)->c
ggarrange(c, temperature, nrow=1, labels=c("A","B"))
quartz.save("~/SurveyPaper/Manuscript_main_figs/tobefig2.pdf", type="pdf")
