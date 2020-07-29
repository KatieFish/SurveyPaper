


#IDing possible "Endemic spp" as those isolated z or more independent times from the same region
z<-5


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

#Are 22 endemics real? All are in more heavily sampled areas. 

