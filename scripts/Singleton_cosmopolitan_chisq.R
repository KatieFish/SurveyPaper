
raw_WY_dataframe<-read.delim("~/SurveyPaper/data/WY_df_2018-02-08.tsv",
                             header=TRUE, stringsAsFactors=FALSE,
                             na.strings=c("","NA"),
                             strip.white=TRUE)

#expected substrate frequencies
tbl<-data.frame(table(unique(raw_WY_dataframe[c(22,29,10)])$Specific))
tbl$Freq<-tbl$Freq/(sum(tbl$Freq))
colnames(tbl)<-c("Specific", "Expected_frequency")
tbl->expected_substrate_frequencies

#expected subphylum frequencies
tbl<-data.frame(table(unique(raw_WY_dataframe[c(22,29,30)])$Subphylum))
tbl$Freq<-tbl$Freq/(sum(tbl$Freq))
colnames(tbl)<-c("Subphylum", "Expected_frequency")
tbl->expected_subphylum_frequencies

#singletons
singletons<-read.delim("~/SurveyPaper/data/Singletons.tsv",
                       header=TRUE, stringsAsFactors=FALSE)





#cosmopolitans
cosmopolitans<-read.delim("~/SurveyPaper/data/Cosmopolitan_z=3.tsv",
                          header=TRUE, stringsAsFactors=FALSE)




