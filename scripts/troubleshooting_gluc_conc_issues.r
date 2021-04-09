
#Downloaded data from strain collection. Cleaned up by hand first. 
gluc_conc_comp_datatable <- read.delim(
  "~/Desktop/gluc_conc_comp_datatable.txt", header=FALSE, na.strings=c("","NA"),
  strip.white=TRUE)
x<-which(gluc_conc_comp_datatable$V2=="unknown")
y<-which(gluc_conc_comp_datatable$V2=="Unknown")
z<-which(gluc_conc_comp_datatable$V2=="UNKNOWN")
a<-which(gluc_conc_comp_datatable$V2=="redo")
b<-which(gluc_conc_comp_datatable$V2=="REDO")
c<-which(grepl("bacteria", gluc_conc_comp_datatable$V2))
d<-which(grepl("Bacteria", gluc_conc_comp_datatable$V2))
e<-which(is.na(gluc_conc_comp_datatable$V3))
torm<-c(x,y,z,a,b,c,d,e)
gluc_conc_comp_datatable<-gluc_conc_comp_datatable[-torm, ]

colnames(gluc_conc_comp_datatable)<-c("StrainID","Species","gluc_percent")

length(which(gluc_conc_comp_datatable$gluc_percent==8))
#1097 at 8%
length(which(gluc_conc_comp_datatable$gluc_percent==0.8))
#2238 at 0.8%
### Looking for assns in ALL data - not just the dataset we're publishing. 
gluc_conc_perms<-Run_WY_association_permutations(
  gluc_conc_comp_datatable, permutations = 10000,
  colnames_to_permute = c("Species", "gluc_percent"))

gluc_conc_sigs<-gluc_conc_perms[c(1:4)]
gluc_conc_sigs$FDR<-p.adjust(gluc_conc_sigs$pval, method="BH")
gluc_conc_sigs<-gluc_conc_sigs[order(gluc_conc_sigs$FDR), ]

## TAXA ARE NOT UNIFORMLY DISTRIBUTED AMONGST GLUCOSE CONCENTRATIONS ##

######

raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                                               header=TRUE, stringsAsFactors=FALSE,
                                               na.strings=c("","NA"),
                                               strip.white=TRUE)

#merged into raw wy dataframe
gluc_comp<-merge(raw_WY_dataframe,
                 gluc_conc_comp_datatable[c(1,3)], by="StrainID", all=TRUE)

#496 total isolations (about 25%) of isolations are accompanied by gluc percent data
length(which(gluc_comp$gluc_percent==8))
#only 2 instances of 8% in the entire dataset
length(which(gluc_comp$gluc_percent==0.8))
#494 instances of .8%

# I can't really do anything with these data. 

#Can I assume that isolates without a sugar indication are 8%?
x<-gluc_conc_comp_datatable[which(gluc_conc_comp_datatable$StrainID %in% raw_WY_dataframe$StrainID), ]
gluc_comp<-merge(raw_WY_dataframe,
                x[c(1,3)], by="StrainID", all=TRUE)
gluc_comp$gluc_percent[which(is.na(gluc_comp$gluc_percent))]<-8
length(which(gluc_comp$gluc_percent==8))
#1468 isolates 
length(which(gluc_comp$gluc_percent==0.8))
#still 494 instances of .8%

#assuming 8% for no data isolates - performing permutations for spp-sugar
species_gluc_conc_assns<-Run_WY_association_permutations(
  all_observations_dataframe = unique(gluc_comp[c(22, 29, 32)]), 
  permutations = 10000,colnames_to_permute = c("Species", "gluc_percent"))

species_gluc_conc_assns<-species_gluc_conc_assns[c(1:4)]
species_gluc_conc_assns$FDR<-p.adjust(species_gluc_conc_assns$pval, method="BH")
species_gluc_conc_assns<-species_gluc_conc_assns[order(species_gluc_conc_assns$FDR), ]

# In this analysis only Torulaspora delbrueckii is significantly ass 
# w/ sugar percent (8%)

#Sacch cerevisiae no enrichment in this dataset. 

# what about quercus? equally sampled in both? 
plntgns_gluc_conc_assns<-Run_WY_association_permutations(
  all_observations_dataframe = unique(gluc_comp[c(22, 15, 32)]), 
  permutations = 10000,colnames_to_permute = c("Plant.Genus", "gluc_percent"))

plntgns_gluc_conc_assns<-plntgns_gluc_conc_assns[c(1:4)]
plntgns_gluc_conc_assns$FDR<-p.adjust(plntgns_gluc_conc_assns$pval, method="BH")
plntgns_gluc_conc_assns<-plntgns_gluc_conc_assns[order(plntgns_gluc_conc_assns$FDR), ]

#quercus is equally sampled in both. there are 3 plant genera that are not
#Fagus  8
#	Acer	8
# Malus	0.8

#none of these genera are in the significant assn's we found. 

spcsub_gluc_conc_assns<-Run_WY_association_permutations(
  all_observations_dataframe = unique(gluc_comp[c(22, 10, 32)]), 
  permutations = 10000,colnames_to_permute = c("Specific", "gluc_percent"))

spcsub_gluc_conc_assns<-spcsub_gluc_conc_assns[c(1:4)]
spcsub_gluc_conc_assns$FDR<-p.adjust(spcsub_gluc_conc_assns$pval, method="BH")
spcsub_gluc_conc_assns<-spcsub_gluc_conc_assns[order(spcsub_gluc_conc_assns$FDR), ]


# Bark enriched for 8% samples
# Bark does show up 2x in specific substrate associations