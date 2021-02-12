
###################################################

require(tidyr)
require(ape)
require(usmap)
require(ggplot2)
require(maps)
require(ggthemes)
require(mapproj)

### DO NOT NEED TO RUN LINES 13:22 - JUST READ IN THE US STATES DF
us_states<-usmap::us_map()
us_states<-us_states[which(!us_states$full == "Hawaii"),]
us_states$full[which(us_states$full=="District of Columbia")]<-"Maryland"
regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
          header=TRUE, stringsAsFactors=FALSE)
colnames(us_states)[c(9,1,2)]<-c("State", "long", "lat")
us_states<-merge(us_states, regions[c(2,3)], by="State", all=TRUE)
us_states<-us_states[order(us_states$order),]
write.table(us_states,"~/SurveyPaper/data/tables_for_scripts/state_coord_for_mapping.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
################

regional_color_key<-read.delim("~/SurveyPaper/data/Geographical_analysis/regional_color_key.tsv",
                               stringsAsFactors=FALSE, header=TRUE)
color_key<-regional_color_key
colnames(color_key)[1]="Clim.Region"
us_states<-read.delim("~/SurveyPaper/data/tables_for_scripts/state_coord_for_mapping.tsv",
                      header=TRUE, stringsAsFactors=FALSE)
regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                    header=TRUE, stringsAsFactors=FALSE)


#######script below plots the number of all independent isolations on the map
TotalsByState<-read.delim("~/SurveyPaper/data/tables_for_scripts/Number_ind_isolations_by_region.tsv")
#adjust lat and long 
temp<-usmap_transform(TotalsByState[c(6,5)])
TotalsByState<-merge(TotalsByState, temp, by=c("LatCenter", "LongCenter"))

p0 <- ggplot(data = us_states,
             mapping = aes(x = long, y = lat,
                           group = group, fill = Region))
 p0 + geom_polygon(color = "gray90", size = 0.1)+
  scale_fill_manual(values = color_key$col)+
  theme_bw()+
  geom_text(data=TotalsByState, aes(x=LongCenter.1, y=LatCenter.1, label=No..independent.isolations, fontface=2), 
            inherit.aes = FALSE, size=3)+ 
  geom_text(data=TotalsByState, aes(x=LongCenter.1, y=LatCenter.1, label=No..unique.species, fontface=3),
            nudge_y=-100500, inherit.aes=FALSE, size=2)+
  coord_equal()

quartz.save("~/SurveyPaper/Figures/Map_of_all_isolates_USA.pdf", type="pdf")



mappable<-raw_WY_dataframe[which(raw_WY_dataframe$Species=="Candida sake"),]
mappable<-mappable[which(!is.na(mappable$Long)),]
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
  geom_text(data=TotalsByState, aes(x=LongCenter.1, y=LatCenter.1, label=No..independent.isolations, fontface=2), 
            inherit.aes = FALSE, size=3)+
  theme(legend.position = "none")+
  coord_equal()

quartz.save("~/SurveyPaper/Figures/Regional_detection_example.pdf", type="pdf")
