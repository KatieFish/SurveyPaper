raw_WY_data<-read.delim("~/SurveyPaper/data/WY_dataframe_w_regions.tsv",
                        header=TRUE, stringsAsFactors=FALSE,
                        strip.white=TRUE)

Regions<- unique(raw_WY_data[c(32)])

#Arctic Unique
ArcticRegion<- raw_WY_data[which(grepl("Arctic", raw_WY_data$Region)), ]
ArcticRegion<- unique(ArcticRegion[c(29,23,32)])
#Unique Frequency
arcticspfreq<- data.frame(table(ArcticRegion$Species))
colnames(arcticspfreq)<- c("Species", "UniqueFreq")
#Total Unique isolations
ArcticUnique<-length(unique(ArcticRegion$SetID))
ArcticRegion$TotalUnique<- ArcticUnique
#Discovery rate
ArcticRegion<-merge(ArcticRegion, arcticspfreq, by="Species")
ArcticRegion$DiscRate<- ArcticRegion$TotalUnique/ArcticRegion$UniqueFreq


#Southeast
SoutheastRegion<- raw_WY_data[which(grepl("Southeast", raw_WY_data$Region)), ]
SoutheastRegion<- unique(SoutheastRegion[c(29, 23, 32)])
southeastspfreq<- data.frame(table(SoutheastRegion$Species))
colnames(southeastspfreq)<-c("Species", "UniqueFreq")

SoutheastUnique<-length(unique(SoutheastRegion$SetID))
SoutheastRegion$TotalUnique<- SoutheastUnique
#Discovery rate
SoutheastRegion<-merge(SoutheastRegion, southeastspfreq, by="Species")
SoutheastRegion$DiscRate<- SoutheastRegion$TotalUnique/SoutheastRegion$UniqueFreq


#OhioValley
OhioValley<- raw_WY_data[which(grepl("Ohio Valley", raw_WY_data$Region)), ]
OhioValley<- unique(OhioValley[c(29, 23, 32)])
ohiovalspfreq<- data.frame(table(OhioValley$Species))
colnames(ohiovalspfreq)<- c("Species", "UniqueFreq")
OhioValleyUnique<-length(unique(OhioValley$SetID))
OhioValley$TotalUnique<- OhioValleyUnique
#Discovery rate
OhioValley<-merge(OhioValley, ohiovalspfreq, by="Species")
OhioValley$DiscRate<- OhioValley$TotalUnique/OhioValley$UniqueFreq


#Upper Midwest
UpperMidwest<- raw_WY_data[which(grepl("Upper Midwest", raw_WY_data$Region)), ]
UpperMidwest<- unique(UpperMidwest[c(29, 23,32)])
uppermidspfrq<- data.frame(table(UpperMidwest$Species))
colnames(uppermidspfrq)<- c("Species", "UniqueFreq")

UpperMidwestUnique<-length(unique(UpperMidwest$SetID))
UpperMidwest$TotalUnique<- UpperMidwestUnique
#Discovery rate
UpperMidwest<-merge(UpperMidwest, uppermidspfrq, by="Species")
UpperMidwest$DiscRate<- UpperMidwest$TotalUnique/UpperMidwest$UniqueFreq

#Northern Rockies and Plains
NorthernRock<- raw_WY_data[which(grepl("Northern Rockies and Plains", raw_WY_data$Region)), ]
NorthernRock<- unique(NorthernRock[c(29,23,32)])
northernrockspfreq<- data.frame(table(NorthernRock$Species))
colnames(northernrockspfreq)<- c("Species", "UniqueFreq") 
NorthernRockUnique<-length(unique(NorthernRock$SetID))
NorthernRock$TotalUnique<- NorthernRockUnique
#Discovery rate
NorthernRock<-merge(NorthernRock, northernrockspfreq, by="Species")
NorthernRock$DiscRate<- NorthernRock$TotalUnique/NorthernRock$UniqueFreq


#Northeast
NortheastRegion<- raw_WY_data[which(grepl("Northeast", raw_WY_data$Region)), ]
NortheastRegion<- unique(NortheastRegion[c(29,23,32)])
northeastspfreq<- data.frame(table(NortheastRegion$Species))
colnames(northeastspfreq)<- c("Species", "UniqueFreq") 
NortheastRegionUnique<-length(unique(NortheastRegion$SetID))
NortheastRegion$TotalUnique<- NortheastRegionUnique
#Discovery rate
NortheastRegion<-merge(NortheastRegion, northeastspfreq, by="Species")
NortheastRegion$DiscRate<- NortheastRegion$TotalUnique/NortheastRegion$UniqueFreq

#Northwest
NorthwestRegion<- raw_WY_data[which(grepl("Northwest", raw_WY_data$Region)), ]
NorthwestRegion<- unique(NorthwestRegion[c(29,23,32)])
northwestspfreq<- data.frame(table(NorthwestRegion$Species))
colnames(northwestspfreq)<- c("Species", "UniqueFreq")
NorthwestRegionUnique<-length(unique(NorthwestRegion$SetID))
NorthwestRegion$TotalUnique<- NorthwestRegionUnique
#Discovery rate
NorthwestRegion<-merge(NorthwestRegion, northwestspfreq, by="Species")
NorthwestRegion$DiscRate<- NorthwestRegion$TotalUnique/NorthwestRegion$UniqueFreq


#West
WestRegion<- raw_WY_data[which(grepl("West", raw_WY_data$Region)), ]
WestRegion<- unique(WestRegion[c(29,23,32)])
westspfreq<- data.frame(table(WestRegion$Species))
colnames(westspfreq)<-c("Species", "UniqueFreq")
WestRegionUnique<-length(unique(WestRegion$SetID))
WestRegion$TotalUnique<- WestRegionUnique
#Discovery rate
WestRegion<-merge(WestRegion, westspfreq, by="Species")
WestRegion$DiscRate<- WestRegion$TotalUnique/WestRegion$UniqueFreq


#South 
SouthRegion<- raw_WY_data[which(grepl("Arkansas", raw_WY_data$State)), ]
Louisiana<- raw_WY_data[which(grepl("Louisiana", raw_WY_data$State)), ]
Mississippi<- raw_WY_data[which(grepl("Mississippi", raw_WY_data$State)), ]
Texas<- raw_WY_data[which(grepl("Texas", raw_WY_data$State)), ]
SouthRegion<-rbind(SouthRegion, Louisiana, Mississippi, Texas)
SouthRegion<- unique(SouthRegion[c(29,23,32)])
southspfreq<- data.frame(table(SouthRegion$Species))
colnames(southspfreq)<-c("Species", "UniqueFreq")
SouthRegionUnique<-length(unique(SouthRegion$SetID))
SouthRegion$TotalUnique<- SouthRegionUnique
#Discovery rate
SouthRegion<-merge(SouthRegion, southspfreq, by="Species")
SouthRegion$DiscRate<- SouthRegion$TotalUnique/SouthRegion$UniqueFreq


##Result Table
RegionalDiscoveryRate<- rbind(ArcticRegion, SoutheastRegion, SouthRegion, UpperMidwest, WestRegion, OhioValley, NorthernRock, NortheastRegion, NorthwestRegion)

##Remove singletons
Singletons<- read.delim("~/SurveyPaper/data/Singletons.tsv", header=TRUE, stringsAsFactors = FALSE)
RegionalDiscoveryRate<- RegionalDiscoveryRate[which(!RegionalDiscoveryRate$Species %in% Singletons$Species), ]
RegionalDiscoveryRate<-unique(RegionalDiscoveryRate[c(1,3,4,5,6)])


##Pipeline
y<-NULL
## loop to find most sampled region each species was found in and determine the discovery rate used for next step 
for (i in 1:length(unique(RegionalDiscoveryRate$Species))) {
  uniquesp<-unique(RegionalDiscoveryRate$Species)
  spmatch<- RegionalDiscoveryRate[which(RegionalDiscoveryRate$Species %in% uniquesp[i]), ] 
  for (j in 1:length(matched)) {
    maxsampled= max(spmatch$TotalUnique)
    spmatch$DiscRateTrue<- spmatch[which(spmatch$TotalUnique>=maxsampled), 5]
  } 
  y<- rbind(y, spmatch)
  
}

x<- unique(y[c(1,6)])
x$ArcticExpected<-NA
x$ArcticFound<-NA
x$SouthExpected<-NA
x$SouthFound<-NA
x$WestExpected<-NA
x$WestFound<-NA
x$SoutheastExpected<-NA
x$SoutheastFound<-NA
x$OhioValleyExpected<-NA
x$OhioValleyFound<-NA
x$UpperMidExpected<-NA
x$UpperMidFound<- NA
x$NorthernRockExpected<-NA
x$NorthernRockFound<-NA
x$NortheastExpected<-NA
x$NortheastFound<-NA
x$NorthwestExpected<-NA
x$NorthwestFound<-NA

#loop to fill in whether or not the species is expected and found for each region
for (i in 1:nrow(x)){
  #Arctic
  if (x$DiscRateTrue[i]<=ArcticUnique) {
  x$ArcticExpected[i]=TRUE
} else {
  x$ArcticExpected[i]=FALSE
}
  if (length(which(x$Species[i] %in% ArcticRegion$Species))>0) {
    x$ArcticFound[i]=TRUE
  } else {
    x$ArcticFound[i]=FALSE
  }
  #South
  if (x$DiscRateTrue[i]<=SouthRegionUnique) {
    x$SouthExpected[i]=TRUE
  } else {
    x$SouthExpected[i]=FALSE
  }
  if (length(which(x$Species[i] %in% SouthRegion$Species))>0) {
    x$SouthFound[i]=TRUE
  } else {
    x$SouthFound[i]=FALSE
  }
  #West
  if (x$DiscRateTrue[i]<=WestRegionUnique) {
    x$WestExpected[i]=TRUE
  } else {
    x$WestExpected[i]=FALSE
  }
  if (length(which(x$Species[i] %in% WestRegion$Species))>0) {
    x$WestFound[i]=TRUE
  } else {
    x$WestFound[i]=FALSE
  }
  #Southeast
  if (x$DiscRateTrue[i]<=SoutheastUnique) {
    x$SoutheastExpected[i]=TRUE
  } else {
    x$SoutheastExpected[i]=FALSE
  }
  if (length(which(x$Species[i] %in% SoutheastRegion$Species))>0) {
    x$SoutheastFound[i]=TRUE
  } else {
    x$SoutheastFound[i]=FALSE
  }
  #Ohio Valley
  if (x$DiscRateTrue[i]<=OhioValleyUnique) {
    x$OhioValleyExpected[i]=TRUE
  } else {
    x$OhioValleyExpected[i]=FALSE
  } 
  if (length(which(x$Species[i] %in% OhioValley$Species))>0) {
    x$OhioValleyFound[i]=TRUE
  } else {
    x$OhioValleyFound[i]=FALSE
  }
  #Upper Mid
  if (x$DiscRateTrue[i]<=UpperMidwestUnique) {
    x$UpperMidExpected[i]=TRUE
  } else {
    x$UpperMidExpected[i]=FALSE
  } 
  if (length(which(x$Species[i] %in% UpperMidwest$Species))>0) {
    x$UpperMidFound[i]=TRUE
  } else {
    x$UpperMidFound[i]=FALSE
  }
  #Northern Rockies and Plains
  if (x$DiscRateTrue[i]<=NorthernRockUnique) {
    x$NorthernRockExpected[i]=TRUE
  } else {
    x$NorthernRockExpected[i]=FALSE
  } 
  if (length(which(x$Species[i] %in% NorthernRock$Species))>0) {
    x$NorthernRockFound[i]=TRUE
  } else {
    x$NorthernRockFound[i]=FALSE
  }
  #Northeast
  if (x$DiscRateTrue[i]<=NortheastRegionUnique) {
    x$NortheastExpected[i]=TRUE
  } else {
    x$NortheastExpected[i]=FALSE
  } 
  if (length(which(x$Species[i] %in% NortheastRegion$Species))>0) {
    x$NortheastFound[i]=TRUE
  } else {
    x$NortheastFound[i]=FALSE
  }
  #Northwest
  if (x$DiscRateTrue[i]<=NorthwestRegionUnique) {
    x$NorthwestExpected[i]=TRUE
  } else {
    x$NorthwestExpected[i]=FALSE
  } 
  if (length(which(x$Species[i] %in% NorthwestRegion$Species))>0) {
    x$NorthwestFound[i]=TRUE
  } else {
    x$NorthwestFound[i]=FALSE
  }
}



Rowstoremove<-c()
#Remove rows where only expected in one region
for (i in 1:nrow(x)) {
  booleanvec<-unlist(x[i,(seq(3,19,2))])
  if (length(which(booleanvec))<2) {
    Rowstoremove=append(Rowstoremove, i)
  }
}

x<-x[-Rowstoremove, ]
Cosmoquery<-data.frame(x[1])
Cosmoquery$FoundExpect<- NA
Cosmoquery$NotFoundExpect<-NA

for (i in 1:nrow(x)) {
  foundwhenexpect<-0
  notfoundwhenexpect<-0
  for (j in seq(3,19,2)) {
  arcticvec<-unlist(x[i,c(j,j+1)])
  if (arcticvec[1]) {
   if (arcticvec[2]) {
     foundwhenexpect<- foundwhenexpect+1
   } else {
     notfoundwhenexpect<- notfoundwhenexpect+1
   }
  }
  }
  Cosmoquery$FoundExpect[i]<-foundwhenexpect
  Cosmoquery$NotFoundExpect[i]<-notfoundwhenexpect
}


length(which(Cosmoquery$FoundExpect>0 & Cosmoquery$NotFoundExpect==0))
#5
length(which(Cosmoquery$FoundExpect>Cosmoquery$NotFoundExpect))
#21
length(which(Cosmoquery$FoundExpect>0 & Cosmoquery$NotFoundExpect<=1))
#22

hist(Cosmoquery$FoundExpect)
hist(Cosmoquery$NotFoundExpect)

Restricted<-Cosmoquery[which(Cosmoquery$NotFoundExpect>5), ]

Cosmos<-Cosmoquery[which(Cosmoquery$FoundExpect>Cosmoquery$NotFoundExpect), ]
Cosmos<-Cosmoquery[which(Cosmoquery$FoundExpect>1 & Cosmoquery$NotFoundExpect<=1), ]
Restricted<-Cosmoquery[which(Cosmoquery$FoundExpect==1 & Cosmoquery$NotFoundExpect>2), ]

write.table(Cosmoquery, "~/SurveyPaper/data/Cosmoquery.tsv", sep = "\t")



##Hist of discovery rates
require(ggplot2)
toplot<-x[c(1,2)]
qplot(toplot$DiscRateTrue,
      geom = "histogram", 
      main = "Histogram of Discovery Rates",
      binwidth=10,
      xlab= "Discovery Rate",
      ylab= "Count",
      col=I("gray"),
      fill=I("#00AFBB"),
      ylim=c(0,15))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black"))+
  scale_x_continuous(breaks= seq(0,100,10))

quartz.save("~/SurveyPaper/Figures/HistDiscoveryRates.pdf", type = "pdf")


###Graph of Restricted unique isos
RestrictedFreq<-raw_WY_data[which(raw_WY_data$Species %in% Restricted$Species), ]
RestrictedFreq<-unique(RestrictedFreq[c(29,23)])
RestrictedFreq<-data.frame(table(RestrictedFreq$Species))
colnames(RestrictedFreq)<-c("Species", "IsolationCount")


ggplot(data=RestrictedFreq, aes(x=Species, y=IsolationCount))+
  geom_bar(stat = "identity", fill= "#00AFBB")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Restricted Isolation Count")+
  scale_y_continuous(breaks = seq(0,10,1))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black"))

quartz.save("~/SurveyPaper/Figures/Restricted_iso_freq.pdf", type = "pdf")

##Cosmopolitan unique iso
CosmoFreq<-raw_WY_data[which(raw_WY_data$Species %in% Cosmos2$Species), ]
CosmoFreq<-unique(CosmoFreq[c(29,23)])
CosmoFreq<-data.frame(table(CosmoFreq$Species))
colnames(CosmoFreq)<-c("Species", "IsolationCount")

ggplot(data=CosmoFreq, aes(x=Species, y=IsolationCount))+
  geom_bar(stat = "identity", fill= "#00AFBB")+
  ggtitle("Cosmopolitan Isolation Count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_continuous(breaks = seq(0,110,10))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black"))

quartz.save("~/SurveyPaper/Figures/Cosmopolitan_Iso_Freq.pdf", type = "pdf")

##Expected Histogram
x$ExpectedHist<-NA
for (i in 1:nrow(x)) {
  freqvec<- unlist(x[i, seq(3,19,2)])
  x$ExpectedHist[i] <- length(which(freqvec))
}

qplot(x$ExpectedHist, geom="histogram", 
      main = "Histogram of Expected Regions",
      binwidth=1,
      xlab= "Expected Regions",
      ylab= "Count",
      col=I("gray"),
      fill=I("#00AFBB"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black"))

quartz.save("~/SurveyPaper/Figures/Expected_Regions_Hist.pdf", type = "pdf")





