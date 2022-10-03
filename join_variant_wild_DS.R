##script HIV

#setwd("/home/fernando/Documentos/Epitopes_pipi/redesign_inmuno_hiv")
ls
list.files()

#####
library(dplyr)
library(tidyverse)


##Open donnor dataset, containing variant epitope information
donnor_variant= unique(read.csv("donnor_hla_variant_epitope.csv"))
donnor_wild= unique(read.csv("donnor_hla_wild_epitope.csv" ))
names(donnor_variant)
names(donnor_variant)
'Next, the datasets will be merged according to the columns "Patient_ID", "allele", "HLA", "Protein" in order to gather the information of the two variants in a single dataset.
As a consequence of combining the variants and the information of each patient, duplicate epitope patterns will be generated,
but which will then be filtered in the python script'

fusion_w_d = left_join(donnor_variant, donnor_wild,  by= c("Patient_ID","allele",'HLA', "Protein"))####fusion dataset by donnor 
fusion_w_d= (unique(fusion_w_d))
names(example)
names(fusion_w_d)[4]<- "Variant_Epitope"
names(fusion_w_d)[6] <-"Epitope_WT"
nrow(unique(fusion_w_d))
nrow(unique(example))

write.csv(fusion_w_d,"fusion_w_d.csv", row.names=FALSE)
nrow(extra)
names(example)
setdiff(  unique(example$Epitope_WT  ), unique(fusion_w_d$Epitope_WT ))
head(extra)
view(extra )

