
### Given the pheWAS will be a matched case-control studies (instead of a population-based analysis), this script is to match 5 controls for each trisomy case 
# Variable used to match, including sex, birth year, and region of birth (better control for SES)
# For QC, we removed individuals who ever moved abroad 



#------------------------------------------------------------------------------------
#---Set environments-----------------------------------------------------------------
#------------------------------------------------------------------------------------

setwd("/home/aoxliu/SCA/SCA_FinnGenR10")
# install.packages("data.table")
library(data.table)
library(dplyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))





####################################################################################
#                    Read in related files and basic QC                            #
####################################################################################

## File with SCA type (1st column is ID, 2nd column is SCA type (XX, XY, Undetermined, XXY, XO, XXX, XYY, XXYY) ) --------------
sca_dat <- data.frame(fread("gsutil cat gs://dsge-aoxing/FinnGenReleaseSandbox/finngen_R10/sca_1.0/data/finngen_R10_sca.tsv", sep="\t", header=T))
dim(sca_dat)   # 345,923      2


## Samples with baseline information --------------
PhenoBase <- data.frame(fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_baseline_internal_1.0.txt.gz|zcat", sep="\t", header=T)) %>% 
      filter(FINNGENID %in% sca_dat$FINNGENID)
dim(PhenoBase)   # 345,898      5


## Only keep individuals with both SCA and baseline information available --------------
sca_dat <- sca_dat %>% filter(FINNGENID %in% PhenoBase$FINNGENID)
dim(sca_dat)   # 345,898      2


## Only keep individuals with SCA as XX, XY, and trisomy (XXY/XXX/XYY) --------------
sca_dat <- sca_dat %>% filter(SCA %in% c("XX","XY","XXY","XXX","XYY"))
dim(sca_dat)   # 343,454      2


## Read in variables used to match (baseline year, baseline age, and region of birth) --------------
cov_dat <- fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_minimum_internal_1.0.txt.gz|zcat", sep="\t", header=T) %>% 
      filter(FINNGENID %in% sca_dat$FINNGENID) %>% 
      mutate(birth_year=BL_YEAR-round(BL_AGE)) %>%    # generate birth year based on the difference between baseline year and baseline age 
      filter(is.na(movedabroad)) %>% 
      select(-HEIGHT, -HEIGHT_AGE, -WEIGHT, -WEIGHT_AGE, -movedabroad, -NUMBER_OF_OFFSPRING)  
dim(cov_dat)   # 342,492     11





####################################################################################
#               The matching process for female sex chromosome trisomy             #
####################################################################################

## Loop for females ---------------------
for (i in c("XXX")) {
	print(i)
	
	
	# List of IDs for individuals with sex chromosome trisomy -----------
	sca_i <- sca_dat %>% filter(SCA==i)
	dim(sca_i)   # 158   2
	
	
	# Extract covariate information ("birth_year" and "regionofbirth") for individuals with sex chromosome trisomy -----------
	cov_sca <- cov_dat %>% filter(FINNGENID %in% sca_i$FINNGENID)
	dim(cov_sca)   # 158  10
	
	
	# Extract XX females (controls) and their covariate information -----------
	cov_xx <- cov_dat %>% filter(FINNGENID %!in% sca_i$FINNGENID) %>% filter(SEX=="female") 
	dim(cov_xx)    # 194,121     11
	
	
	# Generate files as the start step -----------
	d_cov_xx <- cov_xx
	dat_match <- NULL
	
	
	# Matching is performed for each case -----------
	for (j in 1:nrow(cov_sca)) {
		print(j)
		
		# Extract matched variable for each SCA case -----------
		d <- cov_sca[j,]
		d_birth_year <- unlist(d[,"birth_year"])
		d_regionofbirth <- unlist(d[,"regionofbirth"])
		d_BL_YEAR <- unlist(d[,"BL_YEAR"])
		
		# match, when birth year and region of birth matched, try to select samples had closest birth year -----------
		d_match <- d_cov_xx %>% filter(birth_year==d_birth_year & regionofbirth==d_regionofbirth) %>% 
		      mutate(diff_BL_YEAR=abs(BL_YEAR-d_BL_YEAR)) %>% arrange(diff_BL_YEAR)
		d_match_5 <- head(d_match, n=5)
		
		# generate column for the matched group (I used ID of the SCA case), later can be used as strata for the conditional logistic regression -----------
		d_match_5[,"grp"] <- unlist(d[,"FINNGENID"])
		dat_match <- rbind(dat_match, d_match_5)
		d_cov_xx <- d_cov_xx %>% filter(FINNGENID %!in% d_match_5$FINNGENID)
		print(nrow(dat_match))
	}
	
	
	# double check the sex column ---
	dat_match %>% group_by(SEX) %>% count()   # 776 females 
	
	
	# add control to sca phenotype column ---
	dat_match <- dat_match %>% mutate(sca_i=0) %>% 
	      select(FINNGENID, BL_YEAR, BL_AGE, regionofbirth, birth_year, grp, sca_i)
	
	dat_sca <- cov_sca %>% mutate(sca_i=1, grp=FINNGENID) %>% 
	      select(FINNGENID, BL_YEAR, BL_AGE, regionofbirth, birth_year, grp, sca_i)
	

	dat_match1to5 <- rbind(dat_match, dat_sca) 
	dim(dat_match1to5)    # 934   7
	dat_match1to5 %>% group_by(sca_i) %>% count()  # 158 cases and 776 controls 
	
	
	# write out the match list ---
	write.table(dat_match1to5, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	
	
	# full phecodes file ---------------------
	phecode_match1to5 <- data.frame(get(load("20230222_phecode_all_uniq.Rdata"))) %>% filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(phecode_match1to5)    # 28,460     4
	save(phecode_match1to5, file=paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.phecode.Rdata"))
	
	
	# full finngen endpoint file ---------------------
	ep_match1to5 <- fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_endpoint_internal_1.0.txt.gz|zcat", sep="\t", header=T) %>% 
	      filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(ep_match1to5)    # 934 19989
	save(ep_match1to5, file=paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.endpoint.Rdata"))
	
	
	# list of phecodes for running phewas ---------------------
	phecodes <- phecode_match1to5 %>% mutate(phecode=as.character(phecode)) %>% distinct(phecode) %>% unlist()
	length(phecodes)   # 1,236
	write.table(phecodes, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), append=F, quote=F, sep="\t", row.names=F, col.names=F)
	
	
	# list of endpoints for running phewas ---------------------
	endpoints <- colnames(ep_match1to5)[colnames(ep_match1to5) %!in% c("BL_AGE","FU_END_AGE") & substr(colnames(ep_match1to5),nchar(colnames(ep_match1to5))-3,nchar(colnames(ep_match1to5))) %in% "_AGE"]
	endpoints <- gsub("_AGE", "", endpoints)
	length(endpoints)   # 4,996
	write.table(endpoints, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".endpoints.lst"), append=F, quote=F, sep="\t", row.names=F, col.names=F)
}





####################################################################################
#               The matching process for male sex chromosome trisomy               #
####################################################################################

## Loop for males ---------------------
for (i in c("XXY", "XYY")) {
	print(i)
	
	
	# List of IDs for individuals with sex chromosome trisomy -----------
	sca_i <- sca_dat %>% filter(SCA==i)
	dim(sca_i)   # 233   2
	
	
	# Extract covariate information ("birth_year" and "regionofbirth") for individuals with sex chromosome trisomy -----------
	cov_sca <- cov_dat %>% filter(FINNGENID %in% sca_i$FINNGENID)
	dim(cov_sca)   # 158  10
	
	
	# Extract XX females (controls) and their covariate information -----------
	cov_xy <- cov_dat %>% filter(FINNGENID %!in% sca_i$FINNGENID) %>% filter(SEX=="male") 
	dim(cov_xy)    # 147,980     11
	
	
	# Generate files as the start step -----------
	d_cov_xy <- cov_xy
	dat_match <- NULL
	
	
	# Matching is performed for each case -----------
	for (j in 1:nrow(cov_sca)) {
		print(j)
		
		# Extract matched variable for each SCA case -----------
		d <- cov_sca[j,]
		d_birth_year <- unlist(d[,"birth_year"])
		d_regionofbirth <- unlist(d[,"regionofbirth"])
		d_BL_YEAR <- unlist(d[,"BL_YEAR"])
		
		# match, when birth year and region of birth matched, try to select samples had closest birth year -----------
		d_match <- d_cov_xy %>% filter(birth_year==d_birth_year & regionofbirth==d_regionofbirth) %>% 
		      mutate(diff_BL_YEAR=abs(BL_YEAR-d_BL_YEAR)) %>% arrange(diff_BL_YEAR)
		d_match_5 <- head(d_match, n=5)
		
		# generate column for the matched group (I used ID of the SCA case), later can be used as strata for the conditional logistic regression -----------
		d_match_5[,"grp"] <- unlist(d[,"FINNGENID"])
		dat_match <- rbind(dat_match, d_match_5)
		d_cov_xy <- d_cov_xy %>% filter(FINNGENID %!in% d_match_5$FINNGENID)
		print(nrow(dat_match))
	}
	
	
	# double check the sex column ---
	dat_match %>% group_by(SEX) %>% count()   # 1,158 males 
	
	
	# add control to sca phenotype column ---
	dat_match <- dat_match %>% mutate(sca_i=0) %>% 
	      select(FINNGENID, BL_YEAR, BL_AGE, regionofbirth, birth_year, grp, sca_i)
	
	dat_sca <- cov_sca %>% mutate(sca_i=1, grp=FINNGENID) %>% 
	      select(FINNGENID, BL_YEAR, BL_AGE, regionofbirth, birth_year, grp, sca_i)
	

	dat_match1to5 <- rbind(dat_match, dat_sca) 
	dim(dat_match1to5)    # 1391   7
	dat_match1to5 %>% group_by(sca_i) %>% count()  # 233 cases and 1158 controls 
	
	
	# write out the match list ---
	write.table(dat_match1to5, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	
	
	# full phecodes file ---------------------
	phecode_match1to5 <- data.frame(get(load("20230222_phecode_all_uniq.Rdata"))) %>% filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(phecode_match1to5)    # 35,311     4
	save(phecode_match1to5, file=paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.phecode.Rdata"))
	
	
	# full finngen endpoint file ---------------------
	ep_match1to5 <- fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_endpoint_internal_1.0.txt.gz|zcat", sep="\t", header=T) %>% 
	      filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(ep_match1to5)    # 1,391 19989
	save(ep_match1to5, file=paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.endpoint.Rdata"))
	
	
	# list of phecodes for running phewas ---------------------
	phecodes <- phecode_match1to5 %>% mutate(phecode=as.character(phecode)) %>% distinct(phecode) %>% unlist()
	length(phecodes)   # 1,230
	write.table(phecodes, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), append=F, quote=F, sep="\t", row.names=F, col.names=F)
	
	
	# list of endpoints for running phewas ---------------------
	endpoints <- colnames(ep_match1to5)[colnames(ep_match1to5) %!in% c("BL_AGE","FU_END_AGE") & substr(colnames(ep_match1to5),nchar(colnames(ep_match1to5))-3,nchar(colnames(ep_match1to5))) %in% "_AGE"]
	endpoints <- gsub("_AGE", "", endpoints)
	length(endpoints)   # 4,996
	write.table(endpoints, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".endpoints.lst"), append=F, quote=F, sep="\t", row.names=F, col.names=F)
}



