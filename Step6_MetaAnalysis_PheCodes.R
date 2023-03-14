
### Meta-analysis across biobanks and report P from heterogeneity test 

## Given now ukb results are missing, mvp only having results from considering all cases identified from SNP array data
## The following codes are for meta-analysing results from main-analyses (all cases identified from SNP array data) and for mvp and finngen using phecodes 


setwd("/Users/aoxliu/Documents/Project5_SexChromosomeAneuploidy/Script/FinnGen/FinnGenR10_Affymetrix/MetaAnalysisData/")
library(meta)
# install.packages("data.table")
library(data.table)
# install.packages("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))






#########################################################################################
#               Read in results from different biobanks                                 #
#########################################################################################

## Skip this part
#    ## Read in MVP phecode pheWAS results ------------------
#    mvp_xxy_main <- read_excel("RES_MVP_xxy and xyy phewas outcomes 1sep2022.xlsx", sheet="XXY", range=cell_cols("A:O"), col_type="text")
#    colnames(mvp_xxy_main) <- mvp_xxy_main[4,]
#    mvp_xxy_main <- mvp_xxy_main[-(1:4), ]
#    dim(mvp_xxy_main)   # 1,875   15
#    mvp_xxy_main <- read_excel("RES_MVP_xxy and xyy phewas outcomes 1sep2022.xlsx", sheet="XXY", range="A5:P1880", col_type="text")
#    dim(mvp_xxy_main)


## Copy to tsv file and directly read tsv file ------------------
mvp_xxy_main <- fread("RES_MVP_xxy_main_1sep2022.tsv", sep="\t", colClasses=c("character"), header=T)
dim(mvp_xxy_main)   # 1,875   16


mvp_xyy_main <- fread("RES_MVP_xyy_main_1sep2022.tsv", sep="\t", colClasses=c("character"), header=T)  # the sheet missing name of the column L
dim(mvp_xyy_main)   # 1,875   14


mvp_xxx_main <- fread("RES_MVP_xxx_v2_out_with_desc_with_counts.tsv", sep="\t", colClasses=c("character"), header=T)
dim(mvp_xxx_main)   # 1,875   14



## Read in FinnGen phecode pheWAS results ------------------
finngen_xxy_main <- fread("20220310_SCA_FinnGenR10.XXYAll.match1to5.phecodes.tsv", sep="\t", colClasses=c("character"), header=T)
dim(finngen_xxy_main)  # 711  15


finngen_xyy_main <- fread("20220310_SCA_FinnGenR10.XYYAll.match1to5.phecodes.tsv", sep="\t", colClasses=c("character"), header=T)
dim(finngen_xyy_main)  # 662  15


finngen_xxx_main <- fread("20220310_SCA_FinnGenR10.XXXAll.match1to5.phecodes.tsv", sep="\t", colClasses=c("character"), header=T)
dim(finngen_xxx_main)  # 663  15



## Read in UKB phecode pheWAS results [Results are not available yet] ------------------




#########################################################################################
#                    Format the summary stats for MVP                                   #
#########################################################################################

## XXY in MVP ------------------
mvp_xxy_f <- mvp_xxy_main %>% 
      mutate(Disease=phecode,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(SE)), 
             OR=as.numeric(as.character(OR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(p)), 
             N_Cases=as.numeric(as.character(casesWithPhen))+as.numeric(as.character(controlsWithPhen)), 
             N_Controls=as.numeric(as.character(casesWithoutPhen))+as.numeric(as.character(controlsWithoutPhen)), 
             study="MVP",
             SCA="XXY") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(mvp_xxy_f)    # 1087   11


## XYY in MVP ------------------
mvp_xyy_f <- mvp_xyy_main %>% 
      mutate(Disease=phecode,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(SE)), 
             OR=as.numeric(as.character(OR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(p)), 
             N_Cases=as.numeric(as.character(casesWithPhen))+as.numeric(as.character(controlsWithPhen)), 
             N_Controls=as.numeric(as.character(casesWithoutPhen))+as.numeric(as.character(controlsWithoutPhen)), 
             study="MVP",
             SCA="XYY") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(mvp_xyy_f)    # 1,019   11


## XXX in MVP ------------------
mvp_xxx_f <- mvp_xxx_main %>% 
      mutate(Disease=Phecode,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(SE)), 
             OR=as.numeric(as.character(OR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(p)), 
             N_Cases=as.numeric(as.character(casesWithPhen))+as.numeric(as.character(controlsWithPhen)), 
             N_Controls=as.numeric(as.character(casesWithoutPhen))+as.numeric(as.character(controlsWithoutPhen)), 
             study="MVP",
             SCA="XXX") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(mvp_xxx_f)    # 264  11




#########################################################################################
#                    Format the summary stats for FinnGen                               #
#########################################################################################

## XXY in FinnGen ------------------
finngen_xxy_f <- finngen_xxy_main %>% 
      mutate(Disease=phenotype,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(se)), 
             OR=as.numeric(as.character(HR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(P_val)), 
             N_Cases=as.numeric(as.character(n_cases)), 
             N_Controls=as.numeric(as.character(n_controls)), 
             study="FinnGen",
             SCA="XXY") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(finngen_xxy_f)    # 711  11


## XYY in FinnGen ------------------
finngen_xyy_f <- finngen_xyy_main %>% 
      mutate(Disease=phenotype,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(se)), 
             OR=as.numeric(as.character(HR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(P_val)), 
             N_Cases=as.numeric(as.character(n_cases)), 
             N_Controls=as.numeric(as.character(n_controls)), 
             study="FinnGen",
             SCA="XYY") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(finngen_xyy_f)    # 662  11


## XXX in FinnGen ------------------
finngen_xxx_f <- finngen_xxx_main %>% 
      mutate(Disease=phenotype,
             Estimate=as.numeric(as.character(beta)), 
             SE=as.numeric(as.character(se)), 
             OR=as.numeric(as.character(HR)), 
             OR_025=exp(Estimate-1.96*SE), 
             OR_975=exp(Estimate+1.96*SE), 
             P_val=as.numeric(as.character(P_val)), 
             N_Cases=as.numeric(as.character(n_cases)), 
             N_Controls=as.numeric(as.character(n_controls)), 
             study="FinnGen",
             SCA="XXX") %>% 
      filter(N_Cases>=5 & !is.na(Estimate)) %>% 
      select(Disease, Estimate, SE, OR, OR_025, OR_975, P_val, N_Cases, N_Controls, study, SCA)
dim(finngen_xxx_f)    # 663  11






#########################################################################################
#                                     Meta-analysis                                     #
#########################################################################################

## meta-analysis for XXY -------
i <- "XXY"
d <- rbind(mvp_xxy_f, finngen_xxy_f)
dim(d)   # 1798   11

eps <- unique(d$Disease)
for (ep_n in 1:length(eps)){
	print(paste0(ep_n, ": ", eps[ep_n]))
	dat <- d %>% filter(Disease==eps[ep_n])
	m <- metagen(Estimate, SE, data=dat, studlab=paste(study), fixed=TRUE, random=FALSE, prediction=TRUE, sm="SMD")
	
	res <- matrix(c(unlist(dat[1,] %>% select(Disease)), 
		            m$TE.fixed, m$seTE.fixed, m$pval.fixed, exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed), m$pval.Q, sum(dat$N_Cases, na.rm=T), sum(dat$N_Controls, na.rm=T), length(dat$study), 
		           rep(NA, 8*2)), nrow=1)
	
	colnames(res) <- c("Disease",
			           "Estimate.meta", "SE.meta", "P_val.meta", "OR.meta", "OR_025.meta", "OR_975.meta", "Heterogeneity.meta", "N_Cases.meta", " N_Controls.meta", "N_Biobank.meta", 
			           paste0(rep(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"),2), rep(c(".MVP", ".FinnGen"),each=8)))
	
	if (sum(dat$study %in% "MVP")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".MVP")] <- unlist(dat %>% filter(study=="MVP") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (sum(dat$study %in% "FinnGen")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".FinnGen")] <- unlist(dat %>% filter(study=="FinnGen") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (ep_n==1){
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
	}
}



## meta-analysis for XYY -------
i <- "XYY"
d <- rbind(mvp_xyy_f, finngen_xyy_f)
dim(d)   # 1798   11

eps <- unique(d$Disease)
for (ep_n in 1:length(eps)){
	print(paste0(ep_n, ": ", eps[ep_n]))
	dat <- d %>% filter(Disease==eps[ep_n])
	m <- metagen(Estimate, SE, data=dat, studlab=paste(study), fixed=TRUE, random=FALSE, prediction=TRUE, sm="SMD")
	
	res <- matrix(c(unlist(dat[1,] %>% select(Disease)), 
		            m$TE.fixed, m$seTE.fixed, m$pval.fixed, exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed), m$pval.Q, sum(dat$N_Cases, na.rm=T), sum(dat$N_Controls, na.rm=T), length(dat$study), 
		           rep(NA, 8*2)), nrow=1)
	
	colnames(res) <- c("Disease",
			           "Estimate.meta", "SE.meta", "P_val.meta", "OR.meta", "OR_025.meta", "OR_975.meta", "Heterogeneity.meta", "N_Cases.meta", " N_Controls.meta", "N_Biobank.meta", 
			           paste0(rep(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"),2), rep(c(".MVP", ".FinnGen"),each=8)))
	
	if (sum(dat$study %in% "MVP")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".MVP")] <- unlist(dat %>% filter(study=="MVP") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (sum(dat$study %in% "FinnGen")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".FinnGen")] <- unlist(dat %>% filter(study=="FinnGen") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (ep_n==1){
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
	}
}



## meta-analysis for XXX -------
i <- "XXX"
d <- rbind(mvp_xxx_f, finngen_xxx_f)
dim(d)   # 1798   11

eps <- unique(d$Disease)
for (ep_n in 1:length(eps)){
	print(paste0(ep_n, ": ", eps[ep_n]))
	dat <- d %>% filter(Disease==eps[ep_n])
	m <- metagen(Estimate, SE, data=dat, studlab=paste(study), fixed=TRUE, random=FALSE, prediction=TRUE, sm="SMD")
	
	res <- matrix(c(unlist(dat[1,] %>% select(Disease)), 
		            m$TE.fixed, m$seTE.fixed, m$pval.fixed, exp(m$TE.fixed), exp(m$lower.fixed), exp(m$upper.fixed), m$pval.Q, sum(dat$N_Cases, na.rm=T), sum(dat$N_Controls, na.rm=T), length(dat$study), 
		           rep(NA, 8*2)), nrow=1)
	
	colnames(res) <- c("Disease",
			           "Estimate.meta", "SE.meta", "P_val.meta", "OR.meta", "OR_025.meta", "OR_975.meta", "Heterogeneity.meta", "N_Cases.meta", " N_Controls.meta", "N_Biobank.meta", 
			           paste0(rep(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"),2), rep(c(".MVP", ".FinnGen"),each=8)))
	
	if (sum(dat$study %in% "MVP")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".MVP")] <- unlist(dat %>% filter(study=="MVP") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (sum(dat$study %in% "FinnGen")==1){
		res[1,paste0(c("Estimate","SE","P_val","OR","OR_025","OR_975","N_Cases","N_Controls"), ".FinnGen")] <- unlist(dat %>% filter(study=="FinnGen") %>% select(Estimate, SE, P_val, OR, OR_025, OR_975, N_Cases, N_Controls))
	}
	
	if (ep_n==1){
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table(res, paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
	}
}




#########################################################################################
#                              Format meta-analysis results                             #
#########################################################################################


## Read in PheCode, Phenotype, and Category information -----------------
phecoe_def <- read.table("/Users/aoxliu/Documents/Project3_Finngen_mCA/Analysis_PheWAS/phecode_definitions1.2.csv", sep=",", header=T) %>% 
      select(phecode, phenotype, category)
dim(phecoe_def)    # 1867    3


## Add Phenotype and Category information to meta-analysis results -----------------
for (i in c("XXY","XYY","XXX")) {
	print(i)
	
	res <- fread(paste0("META_pheWAS_",i,"_MVP_FinnGen.tsv"), header=T)
	dim(res)   # 1,184   27
	
	res_out <- res %>% inner_join(phecoe_def, by=c("Disease"="phecode"))
	dim(res_out)  # 1,177   29
	
	res %>% filter(Disease %!in% res_out$Disease)   # All phecodes with missing name are 10**, which are for "other test" or null, so ignore them 
	
	res_out <- res_out %>% arrange(desc(-P_val.meta)) %>% 
	     mutate(Estimate.meta=round(Estimate.meta, 3), Estimate.MVP=round(Estimate.MVP, 3), Estimate.FinnGen=round(Estimate.FinnGen, 3), 
	            SE.meta=round(SE.meta, 3), SE.MVP=round(SE.MVP, 3), SE.FinnGen=round(SE.FinnGen, 3), 
	            OR.meta=round(OR.meta, 3), OR.MVP=round(OR.MVP, 3), OR.FinnGen=round(OR.FinnGen, 3), 
	            OR_025.meta=round(OR_025.meta, 3), OR_025.MVP=round(OR_025.MVP, 3), OR_025.FinnGen=round(OR_025.FinnGen, 3), 
	            OR_975.meta=round(OR_975.meta, 3), OR_975.MVP=round(OR_975.MVP, 3), OR_975.FinnGen=round(OR_975.FinnGen, 3), 
	            P_val.meta=ifelse(round(P_val.meta,2)<0.01, formatC(P_val.meta, format="e", digits=2), round(P_val.meta,2)),
	            Heterogeneity.meta=ifelse(round(Heterogeneity.meta,2)<0.01, formatC(Heterogeneity.meta, format="e", digits=2), round(Heterogeneity.meta,2)),
	            P_val.MVP=ifelse(round(P_val.MVP,2)<0.01, formatC(P_val.MVP, format="e", digits=2), round(P_val.MVP, 2)),
	            P_val.FinnGen=ifelse(round(P_val.FinnGen,2)<0.01, formatC(P_val.FinnGen, format="e", digits=2), round(P_val.FinnGen, 2)) )
	write.table(res_out, paste0("META_pheWAS_",i,"_MVP_FinnGen.formatted.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
}



