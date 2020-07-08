
###################################################

#STEP 1 - set the species (spell it right)
species<-"Kluyveromyces marxianus"

require(tidyr)
require(ape)
require(usmap)
require(ggplot2)
require(maps)
require(ggthemes)
require(mapproj)
require(cowplot)
require("seqinr")
require("ips")


### for map
regional_color_key<-read.delim("~/SurveyPaper/data/Geographical_analysis/regional_color_key.tsv",
                      stringsAsFactors=FALSE, header=TRUE)
color_key<-regional_color_key
colnames(color_key)[1]="Clim.Region"
color_key<-color_key[which(!color_key$Clim.Region=="Arctic"),]
color_key<-color_key[order(color_key$Clim.Region),]

us_states<-read.delim("~/SurveyPaper/data/tables_for_scripts/state_coord_for_mapping.tsv",
                      header=TRUE, stringsAsFactors=FALSE)
raw_WY_dataframe<-read.csv("~/SurveyPaper/data/WY_df_2018-02-08.csv",
                           header=TRUE, stringsAsFactors=FALSE,
                           strip.white=TRUE)
regions_by_sp<-read.delim("~/SurveyPaper/data/tables_for_scripts/regions_by_sp.tsv",
                          stringsAsFactors=FALSE, header=TRUE)

spp_of_interest<-raw_WY_dataframe[which(raw_WY_dataframe$Species_Old==species),]

p0 <- ggplot(data = us_states,
             mapping = aes(x = long, y = lat,
                           group = group, fill = Clim.Region))
p1<- p0 + geom_polygon(color = "gray90", size = 0.1)+
  scale_fill_manual(values = color_key$col)+
  theme_bw()+
  geom_jitter(data=spp_of_interest, aes(x=Long, y=Lat), 
              inherit.aes = FALSE, width=.5, shape=1)+ 
  coord_equal()

####STEP 2 - SET THE TREE
tree<-read.tree("/Users/katiefisher/SurveyPaper/data/Geographical_analysis/Intra_species/FASTA_files/trees/Kluyveromyces_marxianus.tree")
tips<-data.frame(tree$tip.label)
tips<-tips %>% separate(tree.tip.label, c("Gen", "Sp", "StrainID"))
tips$order<-c(1:nrow(tips))
tips2<-merge(tips, regions_by_sp[c(2,4)], by="StrainID")
tips2<-merge(tips2, regional_color_key, by="Region")
tips2$col<-as.character(tips2$col)
tips2<-tips2[order(tips2$order),]
par(mar=c(.1,.1,.1,.1))
plot(tree, tip.color = tips2$col,
     use.edge.length = FALSE, edge.width=2)
p2<-recordPlot()

#STEP 3 - CHANGE THE NAME OF THE FIGURE
pdf("~/SurveyPaper/data/Geographical_analysis/Intra_species/Kluyveromyces_marxianus.pdf", height=8, width=6)
plot_grid(p2, p1,
          nrow = 2, rel_heights = c(2, 1))
dev.off()
