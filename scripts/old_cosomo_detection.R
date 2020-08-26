#cosomo def - current - species isolated z times or more from 2 different regions

#read in data
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)

regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                    header=TRUE, stringsAsFactors=FALSE)

raw_WY_dataframe<-merge(raw_WY_dataframe, regions, by="State", all=TRUE)
raw_WY_dataframe<-raw_WY_dataframe[which(!is.na(raw_WY_dataframe$StrainID)),]

#table of unique isolations per species (unique setID)
set_isolates<-data.frame(table(unique(raw_WY_dataframe[c(23,29)])[,2]))
#Set up table of unique region X species combos
set_isolate_regions<-data.frame(table(unique(raw_WY_dataframe[c(22,29,33)])[,c(2,3)]))
set_isolate_regions<-set_isolate_regions[which(set_isolate_regions$Freq>0),]
set_isolate_regions$Species<-as.character(set_isolate_regions$Species)
set_iso_regions<-data.frame(table(set_isolate_regions$Species))


#isolations and regions vectors are parameter space to test for each variable
isos<-seq(2,100,1)
rgns<-seq(1, 10, 1)

#Creating data matrix to populate with results
cosmo_result_mat<-matrix(nrow=length(isos), ncol=length(rgns))
colnames(cosmo_result_mat)<-rgns
row.names(cosmo_result_mat)<-isos

for(i in 1:nrow(cosmo_result_mat)){
  n<-as.numeric(row.names(cosmo_result_mat)[i])
  n_or_more<-as.character(set_isolates$Var1[
    which(set_isolates$Freq>=n)])
  n_or_more_df<-raw_WY_dataframe[which(raw_WY_dataframe$Species %in% n_or_more),]
  for(j in 1:ncol(cosmo_result_mat)){
    z<-as.numeric(colnames(cosmo_result_mat)[j])
    z_or_more<-as.character(set_iso_regions$Var1[
      which(set_iso_regions$Freq>=z)])
    z_or_more_df<-n_or_more_df[which(n_or_more_df$Species %in% z_or_more), ]
    
    #remaining observations are the cosmo yeasts
    cosmo_result_mat[i,j]<-length(unique(z_or_more_df$Species))
  }
}

require(reshape2)
require(ggplot2)

cosmo_parameters_allyeasts<-melt(cosmo_result_mat)
colnames(cosmo_parameters_allyeasts)<-c("Min.Isolations", "Min.Regions", "No.Cosmo")

a<-ggplot(cosmo_parameters_allyeasts,
       aes(x=Min.Isolations, y=No.Cosmo, group=Min.Regions, col=as.factor(Min.Regions)))+
  geom_line()+
  xlab("No. unique isolations")+
  ylab("No. cosmopolitan yeasts")+
  ggtitle("Parameters to identify cosmopolitan yeasts - all subphyla")+
  geom_rect(inherit.aes=FALSE, aes(xmin=0,xmax=25, ymin=0,
                                   ymax=50), col="black", fill=NA)+
  labs(col = "Min. Regions")+
  theme_bw()+
  theme(plot.title = element_text(size=10))
  

  

zoom<-subset(cosmo_parameters_allyeasts, cosmo_parameters_allyeasts$Min.Isolations <= 25 &
               cosmo_parameters_allyeasts$Min.Regions<=50)

b<-ggplot(zoom, aes(x=Min.Isolations, y=No.Cosmo, group=Min.Regions, col=as.factor(Min.Regions)))+
  geom_line()+
  xlab("No. unique isolations")+
  ylab("No. cosmopolitan yeasts")+
  labs(col = "Min. Regions")+
  theme_bw()



####REPEAT - just saccharomycotina
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)

regions<-read.delim("~/SurveyPaper/data/tables_for_scripts/NOAA_US_Climate_Regions.txt",
                    header=TRUE, stringsAsFactors=FALSE)

raw_WY_dataframe<-merge(raw_WY_dataframe, regions, by="State", all=TRUE)
raw_WY_dataframe<-raw_WY_dataframe[which(!is.na(raw_WY_dataframe$StrainID)),]
raw_WY_dataframe<-raw_WY_dataframe[which(raw_WY_dataframe$Subphylum=="Saccharomycotina"),]

#table of unique isolations per species (unique setID)
set_isolates<-data.frame(table(unique(raw_WY_dataframe[c(23,29)])[,2]))
#Set up table of unique region X species combos
set_isolate_regions<-data.frame(table(unique(raw_WY_dataframe[c(22,29,33)])[,c(2,3)]))
set_isolate_regions<-set_isolate_regions[which(set_isolate_regions$Freq>0),]
set_isolate_regions$Species<-as.character(set_isolate_regions$Species)
set_iso_regions<-data.frame(table(set_isolate_regions$Species))


#isolations and regions vectors are parameter space to test for each variable
isos<-seq(2,100,1)
rgns<-seq(1, 10, 1)

#Creating data matrix to populate with results
cosmo_result_mat<-matrix(nrow=length(isos), ncol=length(rgns))
colnames(cosmo_result_mat)<-rgns
row.names(cosmo_result_mat)<-isos

for(i in 1:nrow(cosmo_result_mat)){
  n<-as.numeric(row.names(cosmo_result_mat)[i])
  n_or_more<-as.character(set_isolates$Var1[
    which(set_isolates$Freq>=n)])
  n_or_more_df<-raw_WY_dataframe[which(raw_WY_dataframe$Species %in% n_or_more),]
  for(j in 1:ncol(cosmo_result_mat)){
    z<-as.numeric(colnames(cosmo_result_mat)[j])
    z_or_more<-as.character(set_iso_regions$Var1[
      which(set_iso_regions$Freq>=z)])
    z_or_more_df<-n_or_more_df[which(n_or_more_df$Species %in% z_or_more), ]
    
    #remaining observations are the cosmo yeasts
    cosmo_result_mat[i,j]<-length(unique(z_or_more_df$Species))
  }
}

require(reshape2)
require(ggplot2)

cosmo_parameters_saccharomycotina<-melt(cosmo_result_mat)
colnames(cosmo_parameters_saccharomycotina)<-c("Min.Isolations", "Min.Regions", "No.Cosmo")

c<-ggplot(cosmo_parameters_saccharomycotina,
          aes(x=Min.Isolations, y=No.Cosmo, group=Min.Regions, col=as.factor(Min.Regions)))+
  geom_line()+
  xlab("No. unique isolations")+
  ylab("No. cosmopolitan yeasts")+
  ggtitle("Parameters to identify cosmopolitan yeasts - Saccharomycotina")+
  geom_rect(inherit.aes=FALSE, aes(xmin=0,xmax=25, ymin=0,
                                   ymax=50), col="black", fill=NA)+
  labs(col = "Min. Regions")+
  theme_bw()+
  theme(plot.title = element_text(size=10))

zoom<-subset(cosmo_parameters_saccharomycotina, cosmo_parameters_saccharomycotina$Min.Isolations <= 25 &
               cosmo_parameters_saccharomycotina$Min.Regions<=50)

d<-ggplot(zoom, aes(x=Min.Isolations, y=No.Cosmo, group=Min.Regions, col=as.factor(Min.Regions)))+
  geom_line()+
  xlab("No. unique isolations")+
  ylab("No. cosmopolitan yeasts")+
  labs(col = "Min. Regions")+
  theme_bw()

####
require(ggpubr)
ggarrange(nrow=2, ncol=2, a, c, b, d)
####

quartz.save("Cosmopolitan_by_isolation_and_region_parameter_figure.pdf", type="pdf")

