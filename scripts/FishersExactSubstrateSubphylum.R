
# 1 tailed Fishers Exact 2x2 tables tested for
#enrichment amongst singletons or cosmopolitans for substrate
#enrichment amongst singletons or cosmopolitans for subphylum


raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)

singletons<-read.delim("~/SurveyPaper/data/Singletons.tsv",
                       header=TRUE, stringsAsFactors=FALSE)

cosmopolitans<-read.delim("~/SurveyPaper/data/Cosmopolitan_z=3.tsv",
                          header=TRUE, stringsAsFactors=FALSE)
#create df of all NON-SINGLETON SPP. 
nonsingletons<-raw_WY_dataframe[which(!raw_WY_dataframe$Species %in%
                                        singletons$Species), ]

#create df of all NON-COSMOPOLITAN SPP.
noncosmos<-raw_WY_dataframe[which(!raw_WY_dataframe$Species %in%
                                        cosmopositons$Species), ]

###SPECIFIC SUBSTRATE FISHERS EXACTS
tbl<-data.frame(table(unique(nonsingletons[c(22,29,10)])$Specific))
colnames(tbl)<-c("Specific", "NonsingletonObs")
tbl1<-data.frame(table(unique(singletons[c(22,29,10)])$Specific))
colnames(tbl1)<-c("Specific", "SingletonObs")
tbl2<-data.frame(table(unique(noncosmos[c(22,29,10)])$Specific))
colnames(tbl2)<-c("Specific", "NoncosmopolitanObs")
tbl3<-data.frame(table(unique(cosmopolitans[c(23,29,11)])$Specific))
colnames(tbl3)<-c("Specific", "CosmopolitanObs")
substrate_Fisher_df<-merge(tbl, tbl1, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl2, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl3, by="Specific", all=TRUE)
substrate_Fisher_df[is.na(substrate_Fisher_df)]<-0

##FISHER LOOP for SUBSTRATE
#for each substrate, 2 matrices are indpendently tested:
# ########### substrate # Not substrate
#singleton              #
########################################
#Nonsingleton           #

# ########### substrate # Not substrate
#cosmopolitan           #
########################################
#Noncosmopolitan        #

substrate_Fisher_df$SingletonEnrichment<-NA
substrate_Fisher_df$CosmopolitanEnrichment<-NA
for(i in 1:nrow(substrate_Fisher_df)){
  for (j in c(2,4)){
    substr<-c(substrate_Fisher_df[i, (j+1)], substrate_Fisher_df[i, j])
    NOTsubstr<-c(sum(substrate_Fisher_df[j+1]), sum(substrate_Fisher_df[,j]))-substr
    mat<-as.matrix(data.frame(substr, NOTsubstr))
    if(j==2){
    substrate_Fisher_df$SingletonEnrichment[i]<-as.numeric(fisher.test(mat, alternative="greater")$p.value)
    } else if (j==4){
      substrate_Fisher_df$CosmopolitanEnrichment[i]<-as.numeric(fisher.test(mat, alternative="greater")$p.value)
    }
  } 
}
#P values for 1 tailed fishers exact adjusted using BH
substrate_Fisher_df$SingletonEnrichmentBHadj<-p.adjust(substrate_Fisher_df$SingletonEnrichment, method = "BH")
substrate_Fisher_df$CosmopolitanEnrichmentBHadju<-p.adjust(substrate_Fisher_df$CosmopolitanEnrichment, method="BH")

length(which(substrate_Fisher_df$SingletonEnrichmentBHadj<0.05))
# 0 significant enrichment for singletons in substrates
length(which(substrate_Fisher_df$CosmopolitanEnrichmentBHadj<0.05))
# 1 significant enrichment for cosmopolitan spp.
# cosmopolitan spp. are enriched for soil associations (padj=0.01)

write.table(substrate_Fisher_df, "~/SurveyPaper/data/Singleton_Cosmo_substrate_enrichment_pvals.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)



######SUBPHYUM FISHERS EXACTS
tbl<-data.frame(table(unique(nonsingletons[c(22,29,30)])$Subphylum))
colnames(tbl)<-c("Subphylum", "NonsingletonObs")
tbl1<-data.frame(table(unique(singletons[c(22,29,30)])$Subphylum))
colnames(tbl1)<-c("Subphylum", "SingletonObs")
tbl2<-data.frame(table(unique(noncosmos[c(22,29,30)])$Subphylum))
colnames(tbl2)<-c("Subphylum", "NoncosmopolitanObs")
tbl3<-data.frame(table(unique(cosmopolitans[c(23,29,30)])$Subphylum))
colnames(tbl3)<-c("Subphylum", "CosmopolitanObs")
subphylum_Fisher_df<-merge(tbl, tbl1, by="Subphylum", all=TRUE)
subphylum_Fisher_df<-merge(subphylum_Fisher_df, tbl2, by="Subphylum", all=TRUE)
subphylum_Fisher_df<-merge(subphylum_Fisher_df, tbl3, by="Subphylum", all=TRUE)
subphylum_Fisher_df[is.na(subphylum_Fisher_df)]<-0

##FISHER LOOP for SUBPHYLUM
#for each subphylum, 2 matrices are indpendently tested:
# ########### subphylum # Not subphylum
#singleton              #
########################################
#Nonsingleton           #

# ########### subphylum # Not subphylum
#cosmopolitan           #
########################################
#Noncosmopolitan        #

subphylum_Fisher_df$SingletonEnrichment<-NA
subphylum_Fisher_df$CosmopolitanEnrichment<-NA
for(i in 1:nrow(subphylum_Fisher_df)){
  for (j in c(2,4)){
    substr<-c(subphylum_Fisher_df[i, (j+1)], subphylum_Fisher_df[i, j])
    NOTsubstr<-c(sum(subphylum_Fisher_df[j+1]), sum(subphylum_Fisher_df[,j]))-substr
    mat<-as.matrix(data.frame(substr, NOTsubstr))
    if(j==2){
      subphylum_Fisher_df$SingletonEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }else if(j==4){
      subphylum_Fisher_df$CosmopolitanEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }
  } 
}

#P values for 1 tailed fishers exact adjusted using BH
subphylum_Fisher_df$SingletonEnrichmentBHadj<-p.adjust(subphylum_Fisher_df$SingletonEnrichment, method = "BH")
subphylum_Fisher_df$CosmopolitanEnrichmentBHadju<-p.adjust(subphylum_Fisher_df$CosmopolitanEnrichment, method="BH")

length(which(subphylum_Fisher_df$SingletonEnrichmentBHadj<0.05))
# 2 significant enrichment for singletons in subphylums
# Singletons are enriched for Basidios (p=0.003)
# Singletons enriched for Pezizos (p=0.0005)
length(which(subphylum_Fisher_df$CosmopolitanEnrichmentBHadj<0.05))
# Cosmopolitans are enriched for Saccharomycotina (p=0.0003)

write.table(subphylum_Fisher_df, "~/SurveyPaper/data/Singleton_Cosmo_subphylum_enrichment_pvals.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)



##FISHER LOOP for SUBSTRATE SUBPHYLUM
#for each substrate, 4 matrices are indpendently tested:
# ########### substrate # Not substrate
#Subphylum             #
########################################
#Nonsubphlum           #

raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)

Pezizos<-raw_WY_dataframe[which(raw_WY_dataframe$Subphylum=="Pezizomycotina"), ]
NONPezizos<-raw_WY_dataframe[which(!raw_WY_dataframe$Subphylum=="Pezizomycotina"), ]
Basidios<-raw_WY_dataframe[which(raw_WY_dataframe$Subphylum=="Basidiomycota"), ]
NONBasidios<-raw_WY_dataframe[which(!raw_WY_dataframe$Subphylum=="Basidiomycota"), ]
Saccharos<-raw_WY_dataframe[which(raw_WY_dataframe$Subphylum=="Saccharomycotina"), ]
NONSaccharos<-raw_WY_dataframe[which(!raw_WY_dataframe$Subphylum=="Saccharomycotina"), ]
Taprinos<-raw_WY_dataframe[which(raw_WY_dataframe$Subphylum=="Taphrinomycotina"), ]
NONTaphrinos<-raw_WY_dataframe[which(!raw_WY_dataframe$Subphylum=="Taphrinomycotina"), ]

tbl<-data.frame(table(unique(NONBasidios[c(22,29,10)])$Specific))
colnames(tbl)<-c("Specific", "NonBasidioObs")
tbl1<-data.frame(table(unique(Basidios[c(22,29,10)])$Specific))
colnames(tbl1)<-c("Specific", "BasidioObs")
tbl2<-data.frame(table(unique(NONPezizos[c(22,29,10)])$Specific))
colnames(tbl2)<-c("Specific", "NonPezizoObs")
tbl3<-data.frame(table(unique(Pezizos[c(22,29,10)])$Specific))
colnames(tbl3)<-c("Specific", "PezizoObs")
tbl4<-data.frame(table(unique(NONSaccharos[c(22,29,10)])$Specific))
colnames(tbl4)<-c("Specific", "NonSaccharosObs")
tbl5<-data.frame(table(unique(Saccharos[c(22,29,10)])$Specific))
colnames(tbl5)<-c("Specific", "SaccharoObs")
tbl6<-data.frame(table(unique(NONTaphrinos[c(22,29,10)])$Specific))
colnames(tbl6)<-c("Specific", "NonTaphrinoObs")
tbl7<-data.frame(table(unique(Taprinos[c(22,29,10)])$Specific))
colnames(tbl7)<-c("Specific", "TaphrinoObs")

substrate_Fisher_df<-merge(tbl, tbl1, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl2, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl3, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl4, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl5, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl6, by="Specific", all=TRUE)
substrate_Fisher_df<-merge(substrate_Fisher_df, tbl7, by="Specific", all=TRUE)
substrate_Fisher_df[is.na(substrate_Fisher_df)]<-0

substrate_Fisher_df$BasidioEnrichment<-NA
substrate_Fisher_df$PezizoEnrichment<-NA
substrate_Fisher_df$SaccharoEnrichment<-NA
substrate_Fisher_df$TaphrinoEnrichment<-NA

for(i in 1:nrow(substrate_Fisher_df)){
  for (j in c(2,4,6,8)){
    substr<-c(substrate_Fisher_df[i, (j+1)], substrate_Fisher_df[i, j])
    NOTsubstr<-c(sum(substrate_Fisher_df[j+1]), sum(substrate_Fisher_df[,j]))-substr
    mat<-as.matrix(data.frame(substr, NOTsubstr))
    if(j==2){
      substrate_Fisher_df$BasidioEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }else if(j==4){
      substrate_Fisher_df$PezizoEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }else if(j==6){
      substrate_Fisher_df$SaccharoEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }else if(j==8){
      substrate_Fisher_df$TaphrinoEnrichment[i]<-as.numeric(fisher.test(mat, alternative = "greater")$p.value)
    }
  } 
}

substrate_Fisher_df$BH_BasidioEnrichment<-p.adjust(substrate_Fisher_df$BasidioEnrichment, method = "BH")
substrate_Fisher_df$BH_PezizoEnrichment<-p.adjust(substrate_Fisher_df$PezizoEnrichment, method = "BH")
substrate_Fisher_df$BH_SaccharoEnrichment<-p.adjust(substrate_Fisher_df$SaccharoEnrichment, method = "BH")
substrate_Fisher_df$BH_TaphrinoEnrichment<-p.adjust(substrate_Fisher_df$TaphrinoEnrichment, method = "BH")


write.table(substrate_Fisher_df, "~/SurveyPaper/data/Subphylum_x_substrate_FishersExact.tsv", sep="\t", quote=FALSE, row.names=FALSE)
