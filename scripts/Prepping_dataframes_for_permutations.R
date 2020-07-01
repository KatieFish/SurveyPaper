#We discussed 2 ways of permuting this data. The first is brute force randomizing
#the list of combinations. The second is matrix-based. Here is a loop that will 
#do the first. 

#Read in the dataframe and convert blank cells to NA
raw_WY_dataframe<-read.csv("~/WY_df_2018-02-08.csv",
                  header=TRUE, stringsAsFactors=FALSE,
                  na.strings=c("","NA"),
                  strip.white=TRUE)


####### Association test 1 - plant genus - yeast spp. associations
#adding a step to retain only unique SetIDs
plant_genus_WY_dataframe<-unique(raw_WY_dataframe[c(15,22,29)])

#seperate out strains that have a plant genus associated with them
#I'm using the which function and a boolean operator.  
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!is.na(plant_genus_WY_dataframe$Plant.Genus)), ]
#
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Unknown"), ]
#returns just part of dataframe where Plant.Genus does NOT equal Uknown

#plant_genus_WY_dataframe gets fed into permutation script: 
plant_genus_associations<-Run_WY_association_permutations(all_observations_dataframe=plant_genus_WY_dataframe,
                        permutations=10000, colnames_to_permute=c("Plant.Genus", "Species"))

plant_genus_FDR<-Sampling_FDR(Permutation_df=plant_genus_associations,permutations = 10000,
                              perms_to_sample = 100, pvals = c(.0001, .0005, .001, .005, .01, .05))


####### Association test 1 - isolation temperature - yeast spp. associations
#adding a step to retain only unique SetIDs
iso_temp_WY_dataframe<-unique(raw_WY_dataframe[c(22,23,29)])

#
iso_temp_WY_dataframe<-iso_temp_WY_dataframe[which(!iso_temp_WY_dataframe$IsoTemp=="Unknown"), ]
#returns just part of dataframe where IsoTemp does NOT equal Uknown

#iso_temp_WY_dataframe gets fed into permutation script: 
iso_temp_associations<-Run_WY_association_permutations(all_observations_dataframe=iso_temp_WY_dataframe,
                                                          permutations=10000, colnames_to_permute=c("IsoTemp", "Species"))

#iso_temp_FDR<-Sampling_FDR(Permutation_df=iso_temp_associations,permutations = 10000,
                              perms_to_sample = 100, pvals = c(.0001, .0005, .001, .005, .01, .05))
