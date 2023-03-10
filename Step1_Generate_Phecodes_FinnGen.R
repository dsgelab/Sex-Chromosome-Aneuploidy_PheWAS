
#------------------------------------------------------------------------------------
#---Set environments-----------------------------------------------------------------
#------------------------------------------------------------------------------------

setwd("/home/aoxliu/LOX/pheWASR10/2023Feb22")
library(data.table)
# install.packages("R.utils")
library(R.utils)
library(dplyr)

'%!in%' <- function(x,y)!('%in%'(x,y))


## format ICD codes (remove "." in the middle) ------------
phecode_map <- get(load("/home/aoxliu/PheWAS-master/data/phecode_map.rda")) %>% mutate(code=gsub("\\.", "", code))
phecode_rollup_map <- get(load("/home/aoxliu/PheWAS-master/data/phecode_rollup_map.rda"))
pheinfo.rda <- get(load("/home/aoxliu/PheWAS-master/data/pheinfo.rda"))


## map to phecodes function ------------
mapCodesToPhecodes <-
  function(input, 
           vocabulary.map=phecode_map,
           rollup.map=phecode_rollup_map,
           make.distinct=TRUE) {
    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[["code"]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
    
    if(!is.null(vocabulary.map)){
      #Perform the direct map
      withCallingHandlers(output <- inner_join(input,vocabulary.map,by=c("vocabulary_id","code")), 
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
      #Remove old columns
      output = output %>% select(-code,-vocabulary_id) %>% rename(code=phecode) 
    } else {
      #Warn if the vocabulary IDs are not phecodes
      if(sum(input$vocabulary_id!="phecode")!=0) {warning("Phecode mapping was not requested, but the vocabulary_id of all codes is not 'phecode'")}
      #Prepare for just the phecode expansion
      output=input %>% filter(vocabulary_id=="phecode") %>% select(-vocabulary_id)
    }
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    #Perform the rollup
    if(!is.null(rollup.map)) {
      withCallingHandlers(output <- inner_join(output ,rollup.map,by="code"),
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}}) 
      output = output %>% select(-code) %>% rename(phecode=phecode_unrolled)
      #Make distinct
      if(make.distinct) {output = distinct(output)}
    } else {
      #Rename output column to phecode
      output = output %>% rename(phecode=code)
    }
    
    #Return the output
    output 
  }




#########################################################################################
#                         Extract ICD codes from longitudinal data                      #
#########################################################################################

## read in longitudinal ICD codes -----------
# phe_icd <- fread("gsutil cat gs://finngen-production-library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal_1.0.txt.gz|zcat", header=T, sep="\t") 
phe_icd <- fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_detailed_longitudinal_internal_1.0.txt.gz|zcat", header=T, sep="\t") 
dim(phe_icd)   # 259,233,268        11


## Kela drug reimbursement register (only CODE2 is ICD) -----------
phe_icd_kela <- phe_icd %>% filter(SOURCE=="REIMB" & ICDVER %in% c(9,10))
dim(phe_icd_kela)   # 727,310     11

phe_icd_kela_2 <- phe_icd_kela %>% filter(!is.na(CODE2)) %>% rename(CODE="CODE2") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_kela_2)   # 661,106      4
phe_icd_kela_2 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # only 3,4,5


## Cause of death register (only CODE1 is ICD) -----------
phe_icd_death <- phe_icd %>% filter(SOURCE=="DEATH" & ICDVER %in% c(9,10))
dim(phe_icd_death)   # 134,257     11

phe_icd_death_1 <- phe_icd_death %>% filter(!is.na(CODE1)) %>% rename(CODE="CODE1") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_death_1)   # 134,182      4
phe_icd_death_1 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # only 3,4,5


## Inpatient Hilmo register (only CODE1 and CODE2 are ICD) -----------
phe_icd_inpat <- phe_icd %>% filter(SOURCE=="INPAT" & ICDVER %in% c(9,10))
dim(phe_icd_inpat)   # 6,289,696      11

phe_icd_inpat_1 <- phe_icd_inpat %>% filter(!is.na(CODE1)) %>% rename(CODE="CODE1") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_inpat_1)   # 6,289,696       4
phe_icd_inpat_1 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # also have 1,2,6
phe_icd_inpat_1 <- phe_icd_inpat_1 %>% mutate(code_len=nchar(CODE)) %>% filter(code_len>=3 & code_len<=5) %>% select(-code_len)
dim(phe_icd_inpat_1)   # 6,288,555       4


phe_icd_inpat_2 <- phe_icd_inpat %>% filter(!is.na(CODE2)) %>% rename(CODE="CODE2") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_inpat_2)   # 56,916     4
phe_icd_inpat_2 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # also have 1,2,6
phe_icd_inpat_2 <- phe_icd_inpat_2 %>% mutate(code_len=nchar(CODE)) %>% filter(code_len>=3 & code_len<=5) %>% select(-code_len)
dim(phe_icd_inpat_2)   # 56,915     4


## Inpatient Hilmo register - (only Operations code, so ignore this register) -----------
phe_icd_operin <- phe_icd %>% filter(SOURCE=="OPER_IN" & ICDVER %in% c(9,10))
dim(phe_icd_operin)   # 4,806,874      11
rm(phe_icd_operin)


## Outpatient Hilmo register - (only Operations code, so ignore this register) -----------
phe_icd_operout <- phe_icd %>% filter(SOURCE=="OPER_OUT" & ICDVER %in% c(9,10)) 
dim(phe_icd_operout)   # 10,421,024      11
rm(phe_icd_operout)


## Specialist outpatient Hilmo register -----------
phe_icd_outpat <- phe_icd %>% filter(SOURCE=="OUTPAT" & ICDVER %in% c(9,10))
dim(phe_icd_outpat)   # 33,441,914       11

phe_icd_outpat_1 <- phe_icd_outpat %>% filter(!is.na(CODE1)) %>% rename(CODE="CODE1") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_outpat_1)   # 33,441,914        4
phe_icd_outpat_1 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # also have 1,2,6
phe_icd_outpat_1 <- phe_icd_outpat_1 %>% mutate(code_len=nchar(CODE)) %>% filter(code_len>=3 & code_len<=5) %>% select(-code_len)
dim(phe_icd_outpat_1)   # 33,424,653        4

phe_icd_outpat_2 <- phe_icd_outpat %>% filter(!is.na(CODE2)) %>% rename(CODE="CODE2") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_outpat_2)   # 1,563,190       4
phe_icd_outpat_2 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # also have 1,2,6
phe_icd_outpat_2 <- phe_icd_outpat_2 %>% mutate(code_len=nchar(CODE)) %>% filter(code_len>=3 & code_len<=5) %>% select(-code_len)
dim(phe_icd_outpat_2)   # 1,546,044       4


## Primary health care outpatient visits register -----------
phe_icd_primout <- phe_icd %>% filter(SOURCE=="PRIM_OUT" & ICDVER %in% c(9,10))
dim(phe_icd_primout)   # 82,400,850       11

phe_icd_primout <- phe_icd_primout %>% filter(substr(CATEGORY,1,3)=="ICD")
dim(phe_icd_primout)   # 13,618,414       11

phe_icd_primout_1 <- phe_icd_primout %>% filter(!is.na(CODE1)) %>% rename(CODE="CODE1") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_primout_1)   # 13,596,417        4
phe_icd_primout_1 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   # also have 1,2,6
phe_icd_primout_1 <- phe_icd_primout_1 %>% mutate(code_len=nchar(CODE)) %>% filter(code_len>=3 & code_len<=5) %>% select(-code_len)
dim(phe_icd_primout_1)   # 13,580,262        4


phe_icd_primout_2 <- phe_icd_primout %>% filter(!is.na(CODE2)) %>% rename(CODE="CODE2") %>% select(FINNGENID, EVENT_AGE, CODE, ICDVER)
dim(phe_icd_primout_2)   # 59,241     4
phe_icd_primout_2 %>% mutate(code_len=nchar(CODE)) %>% group_by(code_len) %>% count()   #  only 3,4,5


## Combine all ICD codes -----------
phe_icd_all <- rbind(phe_icd_kela_2, 
                     phe_icd_death_1, 
                     phe_icd_inpat_1, phe_icd_inpat_2, 
                     phe_icd_outpat_1, phe_icd_outpat_2,
                     phe_icd_primout_1, phe_icd_primout_2)
dim(phe_icd_all)  # 55,750,958        4




#########################################################################################
#             Only keep ICD codes for individuals with intensity data                   #
#########################################################################################

sca_dat <- data.frame(fread("gsutil cat gs://dsge-aoxing/FinnGenReleaseSandbox/finngen_R10/sca_1.0/data/finngen_R10_sca.tsv", sep="\t", header=T))
dim(sca_dat)   # 345,923      2


PhenoBase <- data.frame(fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_baseline_internal_1.0.txt.gz|zcat", sep="\t", header=T))
dim(PhenoBase)   # 462,166      5


phe_icd_all <- phe_icd_all %>% filter(FINNGENID %in% sca_dat$FINNGENID) %>% filter(FINNGENID %in% PhenoBase$FINNGENID)
dim(phe_icd_all)   # 41,287,599        4


phe_icd_all %>% group_by(ICDVER) %>% count()
#   ICDVER        n
# 1 10     40632804
# 2 9        654795


save(phe_icd_all, file="20230222_phe_icd_WithIntensity.Rdata")
# phe_icd_all <- data.frame(get(load("20230222_phe_icd_WithIntensity.Rdata")))
dim(phe_icd_all)   # 41,287,599        4




#########################################################################################
#          Convert from original ICD codes to phecode using phewas package              #
#########################################################################################

## add variables required by the phewas package -----------
phe_icd_all <- phe_icd_all %>% mutate(vocabulary_id=paste0("ICD",ICDVER,"CM"), code=as.character(CODE) ) %>% 
     # select(FINNGENID, EVENT_AGE, BL_AGE, vocabulary_id, code)
      select(FINNGENID, EVENT_AGE, vocabulary_id, code)
dim(phe_icd_all)    # 41,287,599        4


## convert from ICD to phecodes -----------
phecode_all <- mapCodesToPhecodes(phe_icd_all)
dim(phecode_all)   # 39,561,624        3 


## first event for each phecode -----------
phecode_all_uniq <- phecode_all %>% select(FINNGENID, EVENT_AGE, phecode) %>% 
                  # arrange(desc(-EVENT_AGE)) %>% 
                    group_by(FINNGENID, phecode) %>% slice_min(EVENT_AGE, n=1)
dim(phecode_all_uniq)   # 9,325,959       3
save(phecode_all_uniq, file="20230222_phecode_all_uniq.Rdata")


## Add BL_AGE -----------
phecode_all_uniq <- phecode_all_uniq %>% inner_join(PhenoBase[,c("FINNGENID","BL_AGE")], by="FINNGENID")
dim(phecode_all_uniq)   # 4,774,069       4
save(phecode_all_uniq, file="20230222_phecode_all_uniq.Rdata")
# phecode_all_uniq <- data.frame(get(load("20230222_phecode_all_uniq.Rdata")))

