color_pal<-(values = c("#00AFBB", "#E7B800", "#FC4E07"))

Plant_associations<-read.delim("~/SurveyPaper/data/Plant_Genus_Associations_10k_permutations.tsv", header=TRUE)

significant<-Plant_associations[which(Plant_associations$pval<=0.001),1:4]

raw_WY_data<-read.csv("~/SurveyPaper/data/WY_df_2018-02-08.csv",
                        header=TRUE, stringsAsFactors=FALSE,
                      strip.white=TRUE)

plant_genus_dataframe<-unique(raw_WY_data[c(15,22,29)])
yeast_freq<-data.frame(table(plant_genus_dataframe$Species))
colnames(yeast_freq)[1]<-"Yeast_sp"
significant<-merge(significant, yeast_freq, by="Yeast_sp")
colnames(significant)[5]<-"Yeast sp. frequency"

plant_freq<-data.frame(table(plant_genus_dataframe$Plant.Genus))
colnames(plant_freq)[1]<-"Plant_genus"
significant<-merge(significant, plant_freq, by="Plant_genus")
colnames(significant)[6]<-"Plant genus frequency"

#write.table(significant, "~/SurveyPaper/data/Significant_0_001_Plant_genus_Yeast_sp_associations.tsv",
            #sep="\t", quote=FALSE, row.names=FALSE)

significant$pair<-paste(significant$Yeast_sp, significant$Plant_genus, sep="-")
significant<-significant[order(significant$Observed, decreasing=TRUE), ]
significant$order<-c(1:nrow(significant))
significant$pair<-reorder(significant$pair, significant$order)
toplot<-significant[c(7, 3, 5, 6)]
xP<-melt(toplot)


ggplot(xP, aes(x=pair, y=value, col=variable,))+
  geom_jitter(size=3, shape=1, stroke=1, width=.15)+
  xlab("")+
  ylab("Frequency observed in dataset")+
  scale_color_manual(values=color_pal)+
  ggtitle("Significant plant genus associations")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

quartz.save("~/SurveyPaper/Figures/Significant_observed_frequencies_plant_associations.pdf", type="pdf")

###########Temperature

temperature_associations<-read.delim("~/SurveyPaper/data/Temperature_Associations_10k_permutations_results.txt",
                                     header=TRUE)

significant<-temperature_associations[which(temperature_associations$pval<=0.01),]

temp_sp_dataframe<-unique(raw_WY_data[c(22,23,29)])

yeast_freq<-data.frame(table(temp_sp_dataframe$Species))
colnames(yeast_freq)<-c("Yeast_sp", "Yeast sp frequency")
significant<-merge(significant, yeast_freq, by="Yeast_sp")
significant$pair<-paste(significant$Yeast_sp, significant$IsoTemp, sep="-")
significant<-significant[order(significant$IsoTemp), ]
significant$order<-c(1:nrow(significant))
significant$pair<-reorder(significant$pair, significant$order)

toplot<-significant[c(6,3,5)]
require(reshape2)
xT<-melt(toplot)

ggplot(xT, aes(x=pair, y=value, col=variable,))+
  geom_jitter(size=3, shape=1, stroke=1, width=.15)+
  xlab("")+
  ylab("Frequency observed in dataset")+
  scale_color_manual(values=color_pal)+
  ggtitle("Significant isolation temperature associations")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
quartz.save("~/SurveyPaper/Significant_observed_frequencies_isotemp_associations.pdf", type="pdf")


#write.table(significant, "~/SurveyPaper/data/Significant_0_01_Isotemp_Yeast_sp_associations.tsv", 
           # sep="\t", quote=FALSE, row.name=FALSE)
  
#########
require(gridExtra)
require(ggpubr)

ggarrange(a,b, nrow = 1)

