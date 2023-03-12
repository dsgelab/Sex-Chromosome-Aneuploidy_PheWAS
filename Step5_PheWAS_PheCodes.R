

### pheWAS using phecodes (ICD9/ICD10)


#------------------------------------------------------------------------------------
#---Set environments-----------------------------------------------------------------
#------------------------------------------------------------------------------------

setwd("/home/aoxliu/LOX/pheWASR10/2023Feb22")
# install.packages("data.table")
library(data.table)
# install.packages("survival")
library("survival")
library(dplyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))




####################################################################################################
#                    PheWAS for all SCA cases identified from SNP array data                       #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.phecode.Rdata"))))
	dim(ep_match1to5)   # 28,460     4
	
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.phecodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.phecodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(phecode==ep)
		
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
			
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.phecodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.phecodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"All.match1to5.phecodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			}
			
		}
	}
}





####################################################################################################
#                PheWAS for SCA cases with clinical diagnoses (ICD9 758; ICD10 Q97-99)             #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.phecode.Rdata"))))
	dim(ep_match1to5)   # 28,460     4
	
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## only keep those with clinical diagnoses -------------------
	phe_sca_all <- read.table("20220310_SCA_FinnGenR10.phe_sca_all.tsv", header=T)
	dim(phe_sca_all)   # 2976    4
	
	dat_match1to5 <- dat_match1to5 %>% filter(grp %in% phe_sca_all$FINNGENID)
	dim(dat_match1to5)  # 419   7
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.phecodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.phecodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		print(ep)
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(phecode==ep)
		
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
				
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.phecodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.phecodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Diagnosed.match1to5.phecodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			} 
			
		}
	}
}





####################################################################################################
#                PheWAS for SCA cases without clinical diagnoses (ICD9 758; ICD10 Q97-99)          #
####################################################################################################

for (i in c("XXX", "XYY", "XXY")) {
	print(i)
	
	## phenotype file -------------------
	ep_match1to5 <- data.frame(get(load(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.phecode.Rdata"))))
	dim(ep_match1to5)   # 28,460     4
	
	
	## matched case-control file -------------------
	dat_match1to5 <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".dat_match1to5.tsv"), header=T)
	dim(dat_match1to5)   # 934   7
	
	
	## only keep those with clinical diagnoses -------------------
	phe_sca_all <- read.table("20220310_SCA_FinnGenR10.phe_sca_all.tsv", header=T)
	dim(phe_sca_all)   # 2976    4
	
	dat_match1to5 <- dat_match1to5 %>% filter(grp %!in% phe_sca_all$FINNGENID)
	dim(dat_match1to5)  # 419   7
	
	
	## phecodes used for analysis -------------------
	codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
	codes <- codes_all$V1
	length(codes)   # 1236
	
	if (file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.phecodes.tsv"))){
		res_done <- read.table(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.phecodes.tsv"), header=T) %>% mutate(phenotype=as.character(phenotype))
		codes_all <- read.table(paste0("20220310_SCA_FinnGenR10.phe_FinnGenR10.",i,".phecodes.lst"), header=F) 
		codes_all <- codes_all %>% filter(V1 %!in% res_done$phenotype) %>% mutate(V1=as.character(V1))
		codes <- codes_all$V1
	}
	
	for (ep in codes){
		print(ep)
		## cases -------------------
		d_case <- ep_match1to5 %>% filter(phecode==ep)
		
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
				
				if (!file.exists(paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.phecodes.tsv"))){
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.phecodes.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
				} else {
					write.table(df, paste0("20220310_SCA_FinnGenR10.",i,"Undiagnosed.match1to5.phecodes.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
				} 
			}
			
		}
	}
}




