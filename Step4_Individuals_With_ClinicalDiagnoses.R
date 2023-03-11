

### Generate file for individuals with clinical diagnoses, based on ICD9 and ICD10 


#------------------------------------------------------------------------------------
#---Set environments-----------------------------------------------------------------
#------------------------------------------------------------------------------------

setwd("/home/aoxliu/LOX/pheWASR10/2023Feb22")
# install.packages("data.table")
library(data.table)
library(dplyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))





####################################################################################################
#                   SCAs with/without clinical diagnoses (ICD9 758 and ICD10 Q97-99)               #
####################################################################################################

## Read in file contains icd codes ------------
phe_icd_all <- data.frame(get(load("20230222_phe_icd_WithIntensity.Rdata")))
dim(phe_icd_all)    # 41,287,599        4


## SCAs diagnosed by ICD9 ------------
phe_sca_icd9 <- phe_icd_all %>% filter(ICDVER==9 & substr(CODE,1,3)==758)
dim(phe_sca_icd9)    # 156   4


## SCAs diagnosed by ICD10 ------------
phe_sca_icd10 <- phe_icd_all %>% filter(ICDVER==10 & substr(CODE,1,3) %in% c("Q97","Q98","Q99"))
dim(phe_sca_icd10)    # 2820    4


## All SCAs diagnosed by either ICD9 or ICD10 ------------
phe_sca_all <- rbind(phe_sca_icd9, phe_sca_icd10)
dim(phe_sca_all)    # 2976    4
write.table(phe_sca_all, "20220310_SCA_FinnGenR10.phe_sca_all.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)


## How many individuals? ------------
phe_sca_all %>% group_by(FINNGENID) %>% count() %>% arrange(desc(n))  # some with many diagnoses
phe_sca_all %>% group_by(FINNGENID) %>% count() %>% arrange(desc(n)) %>% dim()   # 309 




