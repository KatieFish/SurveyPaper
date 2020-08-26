#Examining regionally restricted taxa

regional<-read.delim("~/SurveyPaper/data/Regionally_restricted_foundexp=1_notfoundexp>2.tsv",
                     header=TRUE, stringsAsFactors=FALSE)
cosmodf<-read.delim("~/SurveyPaper/data/Cosmopolitans_foundexp>1_notfoundexp<=1.tsv",
                    header=TRUE, stringsAsFactors=FALSE)

temp<-read.delim("~/SurveyPaper/data/Temperature_Associations_results.tsv")
temsig<-as.character(temp$Yeast_sp[which(temp$BH_Padj<=0.05)])
genus<-read.delim("~/SurveyPaper/data/Plant_Genus_Association_results.tsv")
genussig<-as.character(genus$Yeast_sp[which(genus$BH_Padj<=0.05)])
subst<-read.delim("~/SurveyPaper/data/Substrate_Specific_Associations_results.tsv")
substsig<-as.character(subst$Species[which(subst$BH_Padj<=0.05)])
regionals<-as.character(unique(regional$Species))
cosmos<-as.character(unique(cosmodf$Species))

unique_sp<-append(temsig, genussig)
unique_sp<-append(unique_sp, substsig)
unique_sp<-append(unique_sp, cosmos)
unique_sp<-unique(append(unique_sp, regionals))

regionalupsetdf<-data.frame(unique_sp)
regionalupsetdf$"Regionally restricted"<-0
regionalupsetdf$"Cosmopolitan"<-0
regionalupsetdf$"Temperature-associated"<-0
regionalupsetdf$"Substrate-associated"<-0
regionalupsetdf$"Substrate genus-associated"<-0



for (i in 1:nrow(regionalupsetdf)){
  if (regionalupsetdf$unique_sp[i] %in% regionals){
    regionalupsetdf[i, 2]<-1
  }
  if (regionalupsetdf$unique_sp[i] %in% temsig){
    regionalupsetdf[i, 4]<-1
  }
  if (regionalupsetdf$unique_sp[i] %in% substsig){
    regionalupsetdf[i, 5]<-1
  }
  if (regionalupsetdf$unique_sp[i] %in% genussig){
    regionalupsetdf[i, 6]<-1
  }
  if (regionalupsetdf$unique_sp[i] %in% cosmos){
    regionalupsetdf[i, 3]<-1
  }
}

require(UpSetR)

upset(regionalupsetdf, sets = colnames(regionalupsetdf)[2:6],
      keep.order=TRUE, 
      queries=list(list(query=intersects, params=list("Cosmopolitan", "Substrate-associated"), color="dodgerblue", active=T),
                   list(query=intersects, params=list("Cosmopolitan", "Substrate genus-associated"), color="dodgerblue", active=T),
                   list(query=intersects, params=list("Cosmopolitan", "Temperature-associated"), color="dodgerblue", active=T),
                   list(query=intersects, params=list("Cosmopolitan"), color="dodgerblue", active=T)
      ))

quartz.save("~/SurveyPaper/Figures/Cosmopolitan_associated_upset.pdf", type="pdf")

upset(regionalupsetdf, sets = colnames(regionalupsetdf)[2:6],
      keep.order=TRUE, 
      queries=list(list(query=intersects, params=list("Regionally restricted", "Substrate-associated"), color="seagreen4", active=T), 
                   list(query=intersects, params=list("Regionally restricted", "Substrate genus-associated"), color="seagreen4", active=T),
                   list(query=intersects, params=list("Regionally restricted", "Substrate-associated"), color="seagreen4", active=T),
                   list(query=intersects, params=list("Regionally restricted", "Substrate-associated", "Substrate genus-associated"), color="seagreen4", active=T),
                   list(query=intersects, params=list("Regionally restricted"), color="seagreen4", active=T)
      ))

quartz.save("~/SurveyPaper/Figures/Regional_associated_upset.pdf", type="pdf")

#looking at the intersections
cosmos[which(cosmos %in% substsig)]->x
subst[which(subst$Species %in% x),]->y
#intersection of 4 - 2 soils and 1 bark
cosmos[which(cosmos %in% temsig)]->x
temp[which(temp$Yeast_sp %in% x),]->y
# 2 - both 10 degree (rhodotorula sp and cadida sake)
cosmos[which(cosmos %in% genussig)]->x
genus[which(genus$Yeast_sp %in% x),]->y
# Debaryomyces hansenii - Thuja - yew trees - do appear to be Northern spp.

regionals[which(regionals %in% substsig)]->x
subst[which(subst$Species %in% x),]->y
# Fungus  Suhomyces bolitotheri	
# Lichen	Scleroconidioma sphagnicola
# Sand	Kazachstania serrabonitensis
# Sand	Pichia scaptomyzae
# Needles	Sydowia polyspora
regionals[which(regionals %in% temsig)]->x
temp[which(temp$Species %in% x),]->y
# 0
regionals[which(regionals %in% genussig)]->x
genus[which(genus$Yeast_sp %in% x),]->y
# Candida coipomoensis	Taxus	- yew- maybe restricted to PacNW
# Sydowia polyspora	Picea - spruce - restricted to north?
# Wickerhamomyces onychis	Festuca - cosmo plant genus
# Kwoniella betulae	Betula - birches - resticted to north?
regionals[which(regionals %in% genussig & regionals %in% substsig)]->x
genus[which(genus$Yeast_sp %in% x),]->y
# Sydowia polyspora-Picea (spruce) & needles

#quick map of Sydowia polyspora
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)
#run mapping lines 25:32
mappable<-raw_WY_dataframe[which(!is.na(raw_WY_dataframe$Long)),]
mappable<-mappable[which(mappable$Species=="Sydowia polyspora"),]
temp<-usmap_transform(mappable[c(7,6)])
isolates_to_plot<-merge(mappable, temp, by=c("Long", "Lat"))
isolates_to_plot<-merge(isolates_to_plot, regions[c(2,3)], by="State")


p0 <- ggplot(data = us_states,
             mapping = aes(x = long, y = lat,
                           group = group, fill = Region,
                           alpha=.9), show.legend = FALSE)
p1<-  p0 + geom_polygon(color = "gray90", size = 0.1, show.legend = FALSE)+
  scale_fill_manual(values = color_key$col)+
 # scale_color_manual(values = brewer.pal(name = "Paired", n=11))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_jitter(data=isolates_to_plot, aes(x=Long.1, y=Lat.1), 
              inherit.aes = FALSE, width=.5, shape=1)+ 
  ggtitle("Sydowia polyspora isolation locations")+
  coord_equal()

quartz.save("~/SurveyPaper/Figures/Sydowia_polyspora_iso_locations.pdf", 
            type="pdf")


which(substsig %in% genussig & substsig %in% temsig)->x
substsig[x]
#Lachancea kluyveri associated with 30 degr, bark, & Cercis
