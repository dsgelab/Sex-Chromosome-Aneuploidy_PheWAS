
## PheWAS Kela drug purchase register (ATC codes) -----------


#------------------------------------------------------------------------------------
#---Set environments-----------------------------------------------------------------
#------------------------------------------------------------------------------------

setwd("/home/aoxliu/LOX/pheWASR10/2023Feb22")
library(data.table)
# install.packages("survival")
library("survival")
# install.packages("R.utils")
library(R.utils)
library(dplyr)
# library(lubridate)

'%!in%' <- function(x,y)!('%in%'(x,y))




#########################################################################################
#                         Extract Kela drug purchase register                           #
#########################################################################################

## Read in longitudinal ICD codes -----------
phe_icd <- fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_detailed_longitudinal_internal_1.0.txt.gz|zcat", header=T, sep="\t") 
dim(phe_icd)   # 259,233,268        11


## Extract kela drug purchase register -----------
phe_drug_purchase <- phe_icd %>% filter(SOURCE=="PURCH")
dim(phe_drug_purchase)   # 117,211,255        11


## Remove phe_icd to save memory -----------
rm(phe_icd)



##############################################################################################################
#                 First Kela drug purchase event for SCA matched cases and controls                          #
##############################################################################################################

### Since there are too many drug purchase records, therefore only keept those for matched SCA cases and controls 
## Individuals that are either matched cases or controls -----------
dat_match1to5_XXX <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.XXX.dat_match1to5.tsv"), header=T)
dim(dat_match1to5_XXX)   # 934   7

dat_match1to5_XXY <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.XXY.dat_match1to5.tsv"), header=T)
dim(dat_match1to5_XXY)   # 1391    7

dat_match1to5_XYY <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.XYY.dat_match1to5.tsv"), header=T)
dim(dat_match1to5_XYY)   # 1157    7

dat_match1to5_all <- rbind(dat_match1to5_XXX, dat_match1to5_XXY, dat_match1to5_XYY)
dim(dat_match1to5_all)   # 3482    7


## Drug purchase records for matched cases and controls -----------
drug_purchase_matched <- phe_drug_purchase %>% filter(FINNGENID %in% dat_match1to5_all$FINNGENID)
dim(drug_purchase_matched)   # 802,925     11
# rm(phe_drug_purchase)

save(drug_purchase_matched, file="20230222_drug_purchase.dat_match1to5.Rdata")
# drug_purchase_matched <- data.frame(get(load("20230222_drug_purchase.dat_match1to5.Rdata")))


## First drug purchase record for matched cases and controls -----------
drug_purchase_uniq <- drug_purchase_matched %>% select(FINNGENID, EVENT_AGE, CODE1) %>% 
                    group_by(FINNGENID, CODE1) %>% slice_min(EVENT_AGE, n=1)
dim(drug_purchase_uniq)   # 116,410      3


## Add BL_AGE -----------
PhenoBase <- data.frame(fread("gsutil cat gs://dsge-aoxing/FromFinnGenTeam/PhenotypeBeforeGenotypeQC/data/finngen_R10_baseline_internal_1.0.txt.gz|zcat", sep="\t", header=T))
dim(PhenoBase)   # 462,166      5

drug_purchase_uniq <- drug_purchase_uniq %>% inner_join(PhenoBase[,c("FINNGENID","BL_AGE")], by="FINNGENID")
dim(drug_purchase_uniq)   # 116,410       4
save(drug_purchase_uniq, file="20230222_drug_purchase_uniq.Rdata")
# drug_purchase_uniq <- data.frame(get(load("20230222_drug_purchase_uniq.Rdata")))


# list of phecodes for running phewas ---------------------
drugcodes_uniq <- drug_purchase_uniq %>% mutate(CODE1=as.character(CODE1)) %>% group_by(CODE1) %>% count() %>% filter(n>=5)
drugcodes <- drugcodes_uniq$CODE1
length(drugcodes)   # 1,236
write.table(drugcodes, paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), append=F, quote=F, sep="\t", row.names=F, col.names=F)





####################################################################################################
#                    PheWAS for all SCA cases identified from SNP array data                       #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20230222_drug_purchase_uniq.Rdata")))) %>% filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(ep_match1to5)   # 38,057     4
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.drugcodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.drugcodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(CODE1==ep)
		
		## generate input dat and run the analysis -------------------
		if (nrow(d_case)>=5) {
			dat <- dat_match1to5 %>% mutate(disease=ifelse(FINNGENID %in% d_case$FINNGENID, 1, 0))
			dat$grp <- factor(dat$grp)
			dat$sca_i <- factor(dat$sca_i)
			mod <- clogit(disease ~ sca_i + strata(grp), data=dat) 
			
			df <- as.data.frame(t(as.data.frame(summary(mod)$coeff["sca_i1",])))
			rownames(df) <- "Endpoint"
			df <- df[,c("coef", "exp(coef)", "se(coef)", "z","Pr(>|z|)")]
			colnames(df) <- c("beta", "HR", "se","Z", "P_val")
			df$OR <- summary(mod)$conf.int["sca_i1","exp(coef)"]
			df$OR_025 <- summary(mod)$conf.int["sca_i1","lower .95"]
			df$OR_975 <- summary(mod)$conf.int["sca_i1","upper .95"]
			
			df[1,"n_total"] <- length(unique(dat[,c("FINNGENID")]))   
			df[1,"n_cases"] <- sum(dat[,"disease"]==1)   
			df[1,"n_controls"] <- sum(dat[,"disease"]==0)   
			df[1,"n_strata"] <- length(unique(dat[,c("grp")]))  
			df[1,"phenotype"] <- ep
			df[1,"sca"] <- i
			df[1,"sca_i"] <- "all"
			
			if (df[1,"OR_975"]!="Inf" & !is.na(df[1,"OR_975"])){
				print(df)
			
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.drugcodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.drugcodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.drugcodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			}
			
		}
	}
}

system("gsutil cp  20220310_SCA_FinnGenR10.*All.match1to5.drugcodes.tsv   gs://dsge-aoxing/mocha/FinnGenXXY/")





####################################################################################################
#                PheWAS for SCA cases with clinical diagnoses (ICD9 758; ICD10 Q97-99)             #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## only keep those with clinical diagnoses -------------------
	phe_sca_all <- read.table("20220310_SCA_FinnGenR10.phe_sca_all.tsv", header=T)
	dim(phe_sca_all)   # 2976    4
	
	dat_match1to5 <- dat_match1to5 %>% filter(grp %in% phe_sca_all$FINNGENID)
	dim(dat_match1to5)  # 419   7
	
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20230222_drug_purchase_uniq.Rdata")))) %>% filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(ep_match1to5)   # 38,057     4
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.drugcodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.drugcodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(CODE1==ep)
		
		## generate input dat and run the analysis -------------------
		if (nrow(d_case)>=5) {
			dat <- dat_match1to5 %>% mutate(disease=ifelse(FINNGENID %in% d_case$FINNGENID, 1, 0))
			dat$grp <- factor(dat$grp)
			dat$sca_i <- factor(dat$sca_i)
			mod <- clogit(disease ~ sca_i + strata(grp), data=dat) 
			
			df <- as.data.frame(t(as.data.frame(summary(mod)$coeff["sca_i1",])))
			rownames(df) <- "Endpoint"
			df <- df[,c("coef", "exp(coef)", "se(coef)", "z","Pr(>|z|)")]
			colnames(df) <- c("beta", "HR", "se","Z", "P_val")
			df$OR <- summary(mod)$conf.int["sca_i1","exp(coef)"]
			df$OR_025 <- summary(mod)$conf.int["sca_i1","lower .95"]
			df$OR_975 <- summary(mod)$conf.int["sca_i1","upper .95"]
			
			df[1,"n_total"] <- length(unique(dat[,c("FINNGENID")]))   
			df[1,"n_cases"] <- sum(dat[,"disease"]==1)   
			df[1,"n_controls"] <- sum(dat[,"disease"]==0)   
			df[1,"n_strata"] <- length(unique(dat[,c("grp")]))  
			df[1,"phenotype"] <- ep
			df[1,"sca"] <- i
			df[1,"sca_i"] <- "Diagnosed"
			
			if (df[1,"OR_975"]!="Inf" & !is.na(df[1,"OR_975"])){
				print(df)
			
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.drugcodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.drugcodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.drugcodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			}
			
		}
	}
}

system("gsutil cp  20220310_SCA_FinnGenR10.*Diagnosed.match1to5.drugcodes.tsv   gs://dsge-aoxing/mocha/FinnGenXXY/")





####################################################################################################
#                PheWAS for SCA cases without clinical diagnoses (ICD9 758; ICD10 Q97-99)             #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## only keep those without clinical diagnoses -------------------
	phe_sca_all <- read.table("20220310_SCA_FinnGenR10.phe_sca_all.tsv", header=T)
	dim(phe_sca_all)   # 2976    4
	
	dat_match1to5 <- dat_match1to5 %>% filter(grp %!in% phe_sca_all$FINNGENID)
	dim(dat_match1to5)  # 419   7
	
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20230222_drug_purchase_uniq.Rdata")))) %>% filter(FINNGENID %in% dat_match1to5$FINNGENID)
	dim(ep_match1to5)   # 38,057     4
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.drugcodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.drugcodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.drugcodes.lst"), header=F) %>% mutate(V1=as.character(V1))
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(CODE1==ep)
		
		## generate input dat and run the analysis -------------------
		if (nrow(d_case)>=5) {
			dat <- dat_match1to5 %>% mutate(disease=ifelse(FINNGENID %in% d_case$FINNGENID, 1, 0))
			dat$grp <- factor(dat$grp)
			dat$sca_i <- factor(dat$sca_i)
			mod <- clogit(disease ~ sca_i + strata(grp), data=dat) 
			
			df <- as.data.frame(t(as.data.frame(summary(mod)$coeff["sca_i1",])))
			rownames(df) <- "Endpoint"
			df <- df[,c("coef", "exp(coef)", "se(coef)", "z","Pr(>|z|)")]
			colnames(df) <- c("beta", "HR", "se","Z", "P_val")
			df$OR <- summary(mod)$conf.int["sca_i1","exp(coef)"]
			df$OR_025 <- summary(mod)$conf.int["sca_i1","lower .95"]
			df$OR_975 <- summary(mod)$conf.int["sca_i1","upper .95"]
			
			df[1,"n_total"] <- length(unique(dat[,c("FINNGENID")]))   
			df[1,"n_cases"] <- sum(dat[,"disease"]==1)   
			df[1,"n_controls"] <- sum(dat[,"disease"]==0)   
			df[1,"n_strata"] <- length(unique(dat[,c("grp")]))  
			df[1,"phenotype"] <- ep
			df[1,"sca"] <- i
			df[1,"sca_i"] <- "Undiagnosed"
			
			if (df[1,"OR_975"]!="Inf" & !is.na(df[1,"OR_975"])){
				print(df)
			
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.drugcodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.drugcodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.drugcodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			}
			
		}
	}
}

system("gsutil cp  20220310_SCA_FinnGenR10.*Undiagnosed.match1to5.drugcodes.tsv   gs://dsge-aoxing/mocha/FinnGenXXY/")




