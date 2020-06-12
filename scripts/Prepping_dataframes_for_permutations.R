#We discussed 2 ways of permuting this data. The first is brute force randomizing
#the list of combinations. The second is matrix-based. Here is a loop that will 
#do the first. 

#Read in the dataframe and convert blank cells to NA
raw_WY_dataframe<-read.csv("~/WY_df_2018-02-08.csv",
                  header=TRUE, stringsAsFactors=FALSE,
                  na.strings=c("","NA"),
                  strip.white=TRUE)

#adding a step to retain only unique SetIDs
plant_genus_WY_dataframe<-unique(raw_WY_dataframe[c(15,22,29)])

#seperate out strains that have a plant genus associated with them
#I'm using the which function and a boolean operator.  
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!is.na(plant_genus_WY_dataframe$Plant.Genus)), ]
#
plant_genus_WY_dataframe<-plant_genus_WY_dataframe[which(!plant_genus_WY_dataframe$Plant.Genus=="Unknown"), ]
#returns just part of dataframe where Plant.Genus does NOT equal Uknown


#set up a dataframe with the table function
Permutation_df<-data.frame(table(plant_genus_WY_dataframe$Plant.Genus,
                           plant_genus_WY_dataframe$Species))
#rename the columns
colnames(Permutation_df)<-c("Plant_genus", "Yeast_sp", "Observed")


