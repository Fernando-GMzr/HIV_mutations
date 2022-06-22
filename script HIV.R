##script HIV
library(tidyverse)
setwd("/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv")
ls
list.files()

ResultHIV = read.csv("patterns_5Jun_.csv")
names(ResultHIV)
names(resultOLD)


resultOLD = read.csv("/home/fernando/Documentos/Epitopes_pipi/bajado_drive/HIV_patterns/RESULTS/DATA SETS /11-04-21_Clean Data Set CTL analysis.csv", sep= ";")
head(resultOLD)
nw =unique(resultOLD$Donor.id)
od=unique(ResultHIV$Donor.id)
nrow(ResultHIV)
nrow(resultOLD)
extra = left_join(resultOLD, ResultHIV,  by= c("Donor.id","pattern","sequences"))
extra = left_join( ResultHIV,resultOLD,  by= c("Donor.id","pattern","sequences"))
nrow(extra)
setdiff(y = od, x = nw)
head(extra)
view(extra )

filtrado <-  filter(extra, is.na(Epitope_recover))
nrow(filtrado)
view(filtrado)
filtradold <-  filter(resultOLD, pattern_recover!= "no")
nrow(filtradold )
filtradonw<-  filter(ResultHIV, pattern_recover!= "no")
filtradonw
nrow(filtradonw)
exB <- left_join(filtradonw,filtradold,by= c("Donor.id","pattern","sequences"))
exB
nrow(distinct(exB))
na= filter(exB, is.na(Epitope_recover))
na
view(na)


#####fusion datasets

