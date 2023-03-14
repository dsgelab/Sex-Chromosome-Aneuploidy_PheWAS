setwd("/Users/aoxliu/Documents/Project5_SexChromosomeAneuploidy/Script/Plots/figures")
library(data.table)
library(meta)
library(tidyverse)
library(dplyr)
library(htmlwidgets)
library(plotly)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(egg)
'%!in%' <- function(x,y)!('%in%'(x,y))





######################################################
#                  Set  Variables                    #
######################################################

## One corlor for one disease category
nb.cols <- 17 
mycolors <- data.frame(matrix(c(colorRampPalette(brewer.pal(8, "Set1"))(nb.cols), 
      c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other")), byrow=F, ncol=2))
colnames(mycolors) <- c("mycol","grp_text")

mycolors <- mycolors %>% mutate(grp_text=factor(grp_text, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other"))) %>% arrange(-desc(grp_text))



###################################################################################################
###                    Read in all phewas results first                                           #
###################################################################################################

## Read in all phewas results ------------------
res_dir <- "/Users/aoxliu/Documents/Project5_SexChromosomeAneuploidy/Result/FinnGen/FinnGenR11_Affymetrix/"
res9_dir <- "/Users/aoxliu/Documents/Project5_SexChromosomeAneuploidy/Result/FinnGen/FinnGenR9_Affymetrix/"


# xxx ------
res_xxx <- fread(paste0(res_dir, "PheWAS_SCA_FinnGenR11.phe_FinnGenR10.XXXAll.dat_match1to5.disease.tsv"), sep="\t", header=T) %>% 
      filter(!is.na(P_val) & OR_975!="Inf" & OR_025!="Inf" & OR!="Inf" ) %>% distinct()
dim(res_xxx)   # 1,554   15


# xo ------
res_xo <- fread(paste0(res_dir, "PheWAS_SCA_FinnGenR11.phe_FinnGenR10.XOAll.dat_match1to5.disease.tsv"), sep="\t", header=T) %>% 
      filter(!is.na(P_val) & OR_975!="Inf" & OR_025!="Inf" & OR!="Inf" ) %>% distinct()
dim(res_xo)   # 1,921   15


# xxy ------
res_xxy <- fread(paste0(res9_dir, "PheWAS_FinnGenR9.klinefelterAll.dat_match1to5.disease.tsv"), sep="\t", header=T) %>% 
      filter(!is.na(P_val) & OR_975!="Inf" & OR_025!="Inf" & OR!="Inf" ) %>% distinct()
dim(res_xxy)   # 1,642   14


# xyy ------
res_xyy <- fread(paste0(res9_dir, "PheWAS_FinnGenR9.XYYAll.dat_match1to5.disease.tsv"), sep="\t", header=T) %>% 
      filter(!is.na(P_val) & OR_975!="Inf" & OR_025!="Inf" & OR!="Inf" ) %>% distinct()
dim(res_xyy)   # 1,542   14


## endpoints with phewas results ------------------
ep_phewas <- rbind(res_xxx %>% select(phenotype), 
                   res_xo %>% select(phenotype),
                   res_xxy %>% select(phenotype), 
                   res_xyy %>% select(phenotype)) %>% distinct()
dim(ep_phewas)   # 2,265    1




###################################################################################################
###                   Format endpoints, add ICD group information                                 #
###################################################################################################

## read in endpoint file ------------------
ep_df10 <- fread(paste0(res_dir, "FinnGen_Endpoints_DF10_Final_2022_05_16.tsv"), sep="\t", header=T) %>% 
      select(TAGS, NAME, LONGNAME, SEX, HD_ICD_10, HD_ICD_9, HD_ICD_8) %>% 
      filter(NAME %in% ep_phewas$phenotype) # some with "_EXMORE" were not included in endpoint file
dim(ep_df10)   # 2,097    7


## change covid to J10 ------------------
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(NAME,1,11)=="U22_COVID19", "J10", TAGS))
dim(ep_df10)   # 2,097    7


## change alcohol-related to F5 or K11 ------------------
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(NAME,1,4)=="ALCO" & substr(HD_ICD_10,1,1)=="F", "F5", TAGS))
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(NAME,1,4)=="ALCO" & substr(HD_ICD_10,1,1)=="K", "K11", TAGS))
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(NAME=="ALCO_PRESCRIPTION", "RX", TAGS))
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(NAME=="ALCOHOL_RELATED", "Other", TAGS))
ep_rm <- ep_df10 %>% filter(TAGS=="ALCO" & substr(NAME,1,4)!="ALCO")
ep_df10 <- ep_df10 %>% filter(NAME %!in% ep_rm$NAME)


## change covid to Other ------------------
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(TAGS,1,4)=="MISC", "Other", TAGS))
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(TAGS,1,6)=="SHARED", "Other", TAGS))
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(substr(TAGS,1,5)=="OTHER", "Other", TAGS))

## chaneg to I9 ------------------
ep_df10 <- ep_df10 %>% mutate(TAGS=ifelse(NAME=="CARDIAC_ARRHYTM", "I9", TAGS))


## filter out endpoints that are due to external cause ------------------
ep_df10 <- ep_df10 %>% filter(TAGS %!in% c("ST19", "R18", "Z21", "K11,DENTAL,AVO", "VWXY20", "K11,DENTAL", "DENTAL"))
ep_df10 <- ep_df10 %>% filter(TAGS %!in% c("NEURO_CM", "", "ST19", "Z21,AVO", "ST19,GASTRO_CM,NEURO", "R18,ASTHMA_CM,COPD_CM", "U22", "OPER"))
dim(ep_df10)   # 1,714    7


## filter out endpoints that are difficult to put into a specific group, also there are already many related endpoints included ------------------
ep_df10 <- ep_df10 %>% filter(NAME %!in% c("ST19_TOXIC_EFFECT_ALCOHOL", "ST19_TOXIC_EFFECT_ETHAN", "VWXY20_ACCIDENTAL_POISO_EXPOS_ALCOHOL","Z21_EXTRACO_DIALY", "Z21_PREPAR_CARE_DIALY", "ALCOHOLACUTE9"))
dim(ep_df10)   # 1,708    7


## exclude endpoints with "EXCLUSIONS" in their names ------------------
ep_df10 <- ep_df10 %>% mutate(name_len=nchar(NAME), 
                              name_check=substr(NAME, name_len-9, name_len)) %>% 
                       filter(name_check!="EXCLUSIONS") %>% 
                       select(-name_len, -name_check)
dim(ep_df10)   # 1,706    7


## drop D3_ANAEMIANAS since all cases are totally overlapped with "D3_OTHERANAEMIA" ------------------
ep_df10 <- ep_df10 %>% filter(NAME %!in% c("D3_ANAEMIANAS"))
dim(ep_df10)   # 1,705    7


## drop O15_DELIVERY_PROXY since there is already "O15_DELIVERY" ------------------
ep_df10 <- ep_df10 %>% filter(NAME %!in% c("O15_DELIVERY_PROXY"))
dim(ep_df10)   # 1,704    7


## filter out some endpoints that are not core endpoints (or too many endpoints) ------------------


## endpoint groups (perhaps one additional plot for drug purchase) ------------------
ep_grp <- ep_df10 %>% mutate(TAGS_prefix=gsub(",.*$", "", TAGS)) %>% 
     # group_by(TAGS_prefix) %>% count() %>% arrange(desc(n)) %>% data.frame() %>% 
      mutate(icd_grp=TAGS_prefix) %>% 
      mutate(icd_grp=ifelse(icd_grp=="ICDMAIN", "ICD-10 main chapters", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp=="AB1", "Infectious & parasitic", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("CD2", "C3", "C3,NEURO_CM"), "Neoplasms", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("D3"), "Blood & immune mechanism", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("E4", "DIAB", "DIAB_CM"), "Endocrine, nutritional & metabolic", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("F5", "F5,NEURO_CM", "KRA_PSY"), "Mental & behavioural", icd_grp)) %>% 
      mutate(icd_grp=ifelse((icd_grp %in% c("G6", "G6,GASTRO_CM", "NEURO")) | (NAME %in% c("ALCONEURODEGEN")), "Nervous system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("H7", "H7,OPHTAL"), "Eye & adnexa", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp=="H8", "Ear & mastoid process", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("I9", "I9,CARDIO,NEURO_CM", "CARDIO"), "Circulatory system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("J10", "ASTHMA", "ASTHMA_CM", "COPD", "COPD_CM", "ILD", "ILD_CM"), "Respiratory system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("K11", "K11,GASTRO", "GASTRO", "GASTRO_CM"), "Digestive system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp=="L12", "Skin & subcutaneous tissue", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("M13", "M13,GASTRO_CM", "RHEUMA", "RHEUMA_CM"), "Musculoskeletal system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp=="N14", "Genitourinary system", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("O15","P16"), "Pregnancy, childbirth & puerperium", icd_grp)) %>% 
      mutate(icd_grp=ifelse((icd_grp %in% c("Q17")) | (NAME %in% c("K11_CLEFT_SYNDR")), "Congenital anomalies", icd_grp)) %>% 
      mutate(icd_grp=ifelse(icd_grp %in% c("RX", "RX,I9"), "Drug purchase endpoints", icd_grp)) %>% 
      select(NAME, LONGNAME, SEX, icd_grp)



## update ICD chapter grp names 
icd_grpname <- data.frame(matrix(c("AB1_INFECT_PARASIT", "Infectious & parasitic", 
                               "CD2_NEOPLASM", "Neoplasms", 
                               "D3_BLOOD_IMMUN", "Blood & immune mechanism", 
                               "E4_ENDONUTRMET", "Endocrine, nutritional & metabolic", 
                               "F5_PSYCH", "Mental & behavioural", 
                               "G6_NEURO", "Nervous system", 
                               "H7_EYE_ADNEXA", "Eye & adnexa", 
                               "H8_EAR_MASTOID", "Ear & mastoid process", 
                               "I9_CVD", "Circulatory system",
                               "J10_RESPIRATORY", "Respiratory system", 
                               "K11_GIDISEASES", "Digestive system", 
                               "L12_SKIN_SUBCUTAN", "Skin & subcutaneous tissue", 
                               "M13_MUSCULOSKELETAL", "Musculoskeletal system", 
                               "N14_GENITOURINARY", "Genitourinary system", 
                               "O15_PREG_BIRTH_PUERP", "Pregnancy, childbirth & puerperium",
                               "Q17_CONGEN_MALFO_DEFORMAT_CHROMOSOMAL_ABNORMALITI", "Congenital anomalies"), ncol=2, byrow=T ))
dim(icd_grpname)   # 22   19
colnames(icd_grpname) <- c("phenotype","grpname")




#######################################################################################################################
###           N of diseases being significant for each phenotype, with endpoint info available                        #
#######################################################################################################################

# xxx ------
res_xxx_ep <- res_xxx %>% inner_join(ep_grp, by=c("phenotype"="NAME"))  # 1169   19
res_xxx_ep <- res_xxx_ep %>% mutate(al=ifelse(P_val<=0.05/length(P_val), 1, 0.8))  
res_xxx_ep %>% group_by(al) %>% count()  # 1,077 not significant, and 88 significant 
res_xxx %>% filter(phenotype %!in% ep_grp$NAME)   # check which endpoints are dropped due to missing endpoint info 


# xo ------
res_xo_ep <- res_xo %>% inner_join(ep_grp, by=c("phenotype"="NAME"))  # 1169   19
res_xo_ep <- res_xo_ep %>% mutate(al=ifelse(P_val<=0.05/length(P_val), 1, 0.8))  
res_xo_ep %>% group_by(al) %>% count()  # 1,440 not significant, and 28 significant 
res_xo %>% filter(phenotype %!in% ep_grp$NAME)   # check which endpoints are dropped due to missing endpoint info 


# xxy ------
res_xxy_ep <- res_xxy %>% inner_join(ep_grp, by=c("phenotype"="NAME"))  # 1169   19
res_xxy_ep <- res_xxy_ep %>% mutate(al=ifelse(P_val<=0.05/length(P_val), 1, 0.8))  
res_xxy_ep %>% group_by(al) %>% count()  # 1,151 not significant, and 63 significant 
res_xxy %>% filter(phenotype %!in% ep_grp$NAME)   # check which endpoints are dropped due to missing endpoint info 


# xyy ------
res_xyy_ep <- res_xyy %>% inner_join(ep_grp, by=c("phenotype"="NAME"))  # 1169   19
res_xyy_ep <- res_xyy_ep %>% mutate(al=ifelse(P_val<=0.05/length(P_val), 1, 0.8))  
res_xyy_ep %>% group_by(al) %>% count()  # 1,102 not significant, and 31 significant 
res_xyy %>% filter(phenotype %!in% ep_grp$NAME)   # check which endpoints are dropped due to missing endpoint info 



###################################################################################################
###               Figure 1A: Manhattan plot for Trisomy X syndrome                                #
###################################################################################################

## first panel exclude "ICD-10 main chapters", "Drug purchase endpoints" ------------------
res_xxx_p1 <- res_xxx_ep %>% filter(icd_grp %!in% c("ICD-10 main chapters","Drug purchase endpoints"))
dim(res_xxx_ep)   # 1,169   19
dim(res_xxx_p1)   # 1,080   19


## order by icd group ------------------
res_xxx_p1$icd_grp <- factor(res_xxx_p1$icd_grp, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other"))


## Add variables later used as x-axis and xlabs_text ------------------
ep_seq <- res_xxx_p1 %>% distinct(icd_grp, phenotype) %>% arrange(desc(icd_grp, phenotype)) %>% 
       group_by(icd_grp) %>% mutate(Disease_grp = -row_number()) %>% ungroup() %>% select(phenotype,Disease_grp) %>% mutate(Disease=-row_number())

res_xxx_p1 <- res_xxx_p1 %>% inner_join(ep_seq, by="phenotype") 

hlines <- aggregate(Disease ~ icd_grp, res_xxx_p1, max)[,"Disease"] + 0.5
xlabs <- 0.5 * aggregate(Disease ~ icd_grp, res_xxx_p1, min)[,"Disease"] + 0.5 * aggregate(Disease ~ icd_grp, res_xxx_p1, max)[,"Disease"]
xlabs_text <- aggregate(Disease ~ icd_grp, res_xxx_p1, min)[,"icd_grp"]


## disease to highlight -----------
res_xxx_p1_label <- res_xxx_p1 %>% filter(al==1) %>% 
      filter(phenotype %in% c("AB1_ERYSIPELAS",   # 
                              "G6_EPIPAROX", "G6_XTRAPYROTH", 
                              "H7_LENS",   # "H7_CATARACTSENILE" was dropped due to too high correlation
                              "AUD", "KRA_PSY_ANYMENTAL",   # need to update longname of AUD
                              "I9_VTE", "I9_VARICVE", "FG_CVD", "CARDIAC_ARRHYTM",  # need to update longname of FG_CVD
                              "E4_DM2",   # all dropped are different definitions of T2D 
                              "K11_CHOLELITH",   # all dropped are related
                              "M13_SCIATICA", "M13_KNEEDERANGEMENTS",  # soft tissue; need to update longname of M13_SCIATICA
                              "J10_LOWCHRON", "J10_ASTHMA", "J10_PNEUMONIA",
                              "ALCOHOL_RELATED"))


res_xxx_p1_label %>% filter(phenotype %in% c("AUD", "ALCOHOL_RELATED", "FG_CVD"))
res_xxx_p1_label <- res_xxx_p1_label %>% mutate(LONGNAME=ifelse(phenotype=="AUD", "Alcohol use disorder", LONGNAME), 
                                                LONGNAME=ifelse(phenotype=="ALCOHOL_RELATED", "Alcohol related diseases and deaths", LONGNAME), 
                                                LONGNAME=ifelse(phenotype=="FG_CVD", "Cardiovascular diseases", LONGNAME), 
                                                LONGNAME=ifelse(phenotype=="M13_SCIATICA", "Sciatica with lumbago", LONGNAME) )


# tit <- paste0("Relationship of disease diagnoses with Trisomy X syndrome in women")
tit <- paste0("Relationship of disease diagnoses with Trisomy X syndrome")
pp_xxx <- ggplot(data=res_xxx_p1, aes(x=Disease, y=OR, label=LONGNAME)) + 
             geom_point(shape=25, aes(color=icd_grp, fill=icd_grp, alpha=al), size=2) + 
             geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.6) +
             geom_vline(xintercept=hlines, linetype="dotted", color="grey", size=0.6) + 
             labs(y="Odds ratio",x="",title=tit) + 
             scale_x_continuous(breaks=xlabs, label=xlabs_text) +
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
             scale_y_log10(limits=c(0.5,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
             scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_shape_manual(values=c(25))+ 
             theme_classic() + 
            # theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
             theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) + 
             theme(legend.position="none") + 
             theme(strip.background=element_blank(), strip.text.x=element_blank()) + 
             geom_text_repel(data=res_xxx_p1_label, size=3.5, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black", min.segment.length=0)
pp_xxx





###################################################################################################
###                Figure 1B: Manhattan plot for Turner syndrome                                  #
###################################################################################################

## first panel exclude "ICD-10 main chapters", "Drug purchase endpoints" ------------------
res_xo_p1 <- res_xo_ep %>% filter(icd_grp %!in% c("ICD-10 main chapters","Drug purchase endpoints"))
dim(res_xo_ep)   # 1,169   19
dim(res_xo_p1)   # 1,080   19


## order by icd group ------------------
res_xo_p1$icd_grp <- factor(res_xo_p1$icd_grp, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other"))


## Add variables later used as x-axis and xlabs_text ------------------
ep_seq <- res_xo_p1 %>% distinct(icd_grp, phenotype) %>% arrange(desc(icd_grp, phenotype)) %>% 
       group_by(icd_grp) %>% mutate(Disease_grp = -row_number()) %>% ungroup() %>% select(phenotype,Disease_grp) %>% mutate(Disease=-row_number())

res_xo_p1 <- res_xo_p1 %>% inner_join(ep_seq, by="phenotype") 

hlines <- aggregate(Disease ~ icd_grp, res_xo_p1, max)[,"Disease"] + 0.5
xlabs <- 0.5 * aggregate(Disease ~ icd_grp, res_xo_p1, min)[,"Disease"] + 0.5 * aggregate(Disease ~ icd_grp, res_xo_p1, max)[,"Disease"]
xlabs_text <- aggregate(Disease ~ icd_grp, res_xo_p1, min)[,"icd_grp"]


## disease to highlight -----------
res_xo_p1_label <- res_xo_p1 %>% filter(al==1) %>% 
      filter(phenotype %in% c("DEATH", 
                              "GOUT",   # need to update longname of GOUT 
                              "D3_OTHERANAEMIA",   # drop D3_ANAEMIANAS since all cases are overlapped 
                              "O15_DELIVERY",    # drop "O15_DELIVERY_PROXY" and "O15_DELIV_SPONT"
                              "H8_CONSENHEARINGLOSS", "H8_MIDDLEMASTOID",   # "H8_HL_SEN_NAS" is included in "H8_CONSENHEARINGLOSS" 
                              "Q17_CHROMOSOMAL_ABNORMALITI_NOT_ELSEW_CLASSIFIED", "Q17_LVOTO_NARROW", "Q17_LVOTO_BROAD", "Q17_CONGEN_MALFO_AORTIC_MITRAL_VALVES",  # drop "Q17_CONGEN_MALFO_CIRCULATO_SYSTEM" and "Q17_VACTREL" since there are already simliar endpoints highlight
                              "I9_CHD", "I9_VALVES", "I9_AORTANEUR", "I9_DOAAC", "CONGEN_HEART_ARTER"))  # drop "I9_VHD" since there are already simliar endpoints highlight
dim(res_xo_p1_label)   # 18 21
res_xo_p1_label <- res_xo_p1_label %>% mutate(LONGNAME=ifelse(phenotype=="GOUT", "Gout", LONGNAME) )


# tit <- paste0("Relationship of disease diagnoses with Turner syndrome in women")
tit <- paste0("Relationship of disease diagnoses with Turner syndrome")
pp_xo <- ggplot(data=res_xo_p1, aes(x=Disease, y=OR, label=LONGNAME)) + 
             geom_point(shape=25, aes(color=icd_grp, fill=icd_grp, alpha=al), size=2) + 
             geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.6) +
             geom_vline(xintercept=hlines, linetype="dotted", color="grey", size=0.6) + 
             labs(y="Odds ratio",x="",title=tit) + 
             scale_x_continuous(breaks=xlabs, label=xlabs_text) +
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
             scale_y_log10(limits=c(0.5,111),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100)) + 
             scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_shape_manual(values=c(25))+ 
             theme_classic() + 
            # theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
             theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) + 
             theme(legend.position="none") + 
             theme(strip.background=element_blank(), strip.text.x=element_blank()) + 
             geom_text_repel(data=res_xo_p1_label, size=3.5, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black", min.segment.length=0)
pp_xo




###################################################################################################
###           Figure 1C: Manhattan plot for Trisomy X syndrome by disease category                #
###################################################################################################

## second panel for "ICD-10 main chapters" ------------------
res_xxx_p2 <- res_xxx_ep %>% filter(icd_grp %in% c("ICD-10 main chapters")) %>% 
      inner_join(icd_grpname, by="phenotype") 
dim(res_xxx_p2)   # 16   19

## order by icd group ------------------
res_xxx_p2$grpname <- factor(res_xxx_p2$grpname, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies"))

pp_xxx_category <- ggplot(data=res_xxx_p2, aes(x=grpname, y=OR, ymin=OR_025, ymax=OR_975 )) + 
		geom_pointrange(size=0.5, shape=22, aes(color=grpname, fill=grpname)) +
		geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) + 
		theme(legend.position="none") + 
		labs(y="Odds ratio", x="", title=paste0("")) +
		# scale_y_log10(limits=c(0.3,25), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		# scale_y_log10(limits=c(0.2,10), breaks=c(0.5,1,2,3,4,6,8,10), labels=c(0.5,1,2,3,4,6,8,10)) + 
		scale_y_log10(limits=c(0.2,30), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() +
		# theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) +
		theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) +
		theme(legend.position="none")
pp_xxx_category




###################################################################################################
###           Figure 1D: Manhattan plot for Turner syndrome by disease category                   #
###################################################################################################

## second panel for "ICD-10 main chapters" ------------------
res_xo_p2 <- res_xo_ep %>% filter(icd_grp %in% c("ICD-10 main chapters")) %>% 
      inner_join(icd_grpname, by="phenotype") 
dim(res_xo_p2)   # 16   20

## order by icd group ------------------
res_xo_p2$grpname <- factor(res_xo_p2$grpname, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies"))

pp_xo_category <- ggplot(data=res_xo_p2, aes(x=grpname, y=OR, ymin=OR_025, ymax=OR_975 )) + 
		geom_pointrange(size=0.5, shape=22, aes(color=grpname, fill=grpname)) +
		geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) + 
		theme(legend.position="none") + 
		labs(y="Odds ratio", x="", title=paste0("")) +
		# scale_y_log10(limits=c(0.3,25), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_y_log10(limits=c(0.2,10), breaks=c(0.5,1,2,3,4,6,8,10), labels=c(0.5,1,2,3,4,6,8,10)) + 
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() +
		# theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) +
		theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) +
		theme(legend.position="none")
pp_xo_category



###################################################################################################
###                                   Align ABCD for Figure 1                                     #
###################################################################################################

## align plots of Trisomy X syndrome and Turner syndrome
p_xxx <- grid.arrange(pp_xxx, pp_xxx_category, ncol=2, widths=c(2.5,1), top=" ", bottom="", left=" ", right=" ")
p_xo <- grid.arrange(pp_xo, pp_xo_category, ncol=2, widths=c(2.5,1), top=" ", bottom="", left=" ", right=" ")

p <- grid.arrange(p_xxx, p_xo, nrow=2)
# ggsave(paste0("Figure1.Relationship_WithDisease_XXX_XO.FinnGen.png"), p, width=9.6, height=12.5)
ggsave(paste0("Figure1.Relationship_WithDisease_XXX_XO.FinnGen.png"), p, width=13.5, height=16)





###################################################################################################
###                   Figure 2A: Manhattan plot for Klinefelter syndrome                          #
###################################################################################################

## first panel exclude "ICD-10 main chapters", "Drug purchase endpoints" ------------------
res_xxy_p1 <- res_xxy_ep %>% filter(icd_grp %!in% c("ICD-10 main chapters","Drug purchase endpoints","Pregnancy, childbirth & puerperium"))
dim(res_xxy_ep)   # 1,169   19
dim(res_xxy_p1)   # 1,080   19


## order by icd group ------------------
res_xxy_p1$icd_grp <- factor(res_xxy_p1$icd_grp, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other"))


## Add variables later used as x-axis and xlabs_text ------------------
ep_seq <- res_xxy_p1 %>% distinct(icd_grp, phenotype) %>% arrange(desc(icd_grp, phenotype)) %>% 
       group_by(icd_grp) %>% mutate(Disease_grp = -row_number()) %>% ungroup() %>% select(phenotype,Disease_grp) %>% mutate(Disease=-row_number())

res_xxy_p1 <- res_xxy_p1 %>% inner_join(ep_seq, by="phenotype") 

hlines <- aggregate(Disease ~ icd_grp, res_xxy_p1, max)[,"Disease"] + 0.5
xlabs <- 0.5 * aggregate(Disease ~ icd_grp, res_xxy_p1, min)[,"Disease"] + 0.5 * aggregate(Disease ~ icd_grp, res_xxy_p1, max)[,"Disease"]
xlabs_text <- aggregate(Disease ~ icd_grp, res_xxy_p1, min)[,"icd_grp"]


## disease to highlight -----------
res_xxy_p1_label <- res_xxy_p1 %>% filter(al==1) %>% 
      filter(phenotype %in% c("AB1_ERYSIPELAS", "AB1_INTESTINAL_INFECTIONS", 
                              "I9_VARICVE", "I9_HEARTFAIL_NS", 
                              "E4_DM2",   # all dropped are different definitions of T2D 
                              "M13_KNEEDERANGEMENTS",  # soft tissue; need to update longname of M13_SCIATICA
                              "J10_LOWCHRON",
                              "E4_TESTICHYPO", "N14_MALEINFERT", "E4_OBESITY", "ASTHMA_PNEUMONIA", "PAIN", "ANY_REIMB", "K11_GULC", "F5_PSYCHDEV"))
dim(res_xxy_p1_label)   # 18 21
res_xxy_p1_label <- res_xxy_p1_label %>% mutate(LONGNAME=ifelse(phenotype=="PAIN", "Pain", LONGNAME), 
                                                LONGNAME=ifelse(phenotype=="ANY_REIMB", "Any reimbursement code in KELA", LONGNAME))



# tit <- paste0("Relationship of disease diagnoses with Klinefelter syndrome in men")
tit <- paste0("Relationship of disease diagnoses with Klinefelter syndrome")
pp_xxy <- ggplot(data=res_xxy_p1, aes(x=Disease, y=OR, label=LONGNAME)) + 
             geom_point(shape=25, aes(color=icd_grp, fill=icd_grp, alpha=al), size=2) + 
             geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.6) +
             geom_vline(xintercept=hlines, linetype="dotted", color="grey", size=0.6) + 
             labs(y="Odds ratio",x="",title=tit) + 
             scale_x_continuous(breaks=xlabs, label=xlabs_text) +
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
             scale_y_log10(limits=c(0.5,111),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100)) + 
             scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_shape_manual(values=c(25))+ 
             theme_classic() + 
            # theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
             theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) + 
             theme(legend.position="none") + 
             theme(strip.background=element_blank(), strip.text.x=element_blank()) + 
             geom_text_repel(data=res_xxy_p1_label, size=3.5, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black", min.segment.length=0)
pp_xxy





###################################################################################################
###                   Figure 2B: Manhattan plot for Jacobs syndrome                               #
###################################################################################################

## first panel exclude "ICD-10 main chapters", "Drug purchase endpoints" ------------------
res_xyy_p1 <- res_xyy_ep %>% filter(icd_grp %!in% c("ICD-10 main chapters","Drug purchase endpoints","Pregnancy, childbirth & puerperium"))
dim(res_xyy_ep)   # 1,169   19
dim(res_xyy_p1)   # 1,080   19


## order by icd group ------------------
res_xyy_p1$icd_grp <- factor(res_xyy_p1$icd_grp, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies", "Other"))


## Add variables later used as x-axis and xlabs_text ------------------
ep_seq <- res_xyy_p1 %>% distinct(icd_grp, phenotype) %>% arrange(desc(icd_grp, phenotype)) %>% 
       group_by(icd_grp) %>% mutate(Disease_grp = -row_number()) %>% ungroup() %>% select(phenotype,Disease_grp) %>% mutate(Disease=-row_number())

res_xyy_p1 <- res_xyy_p1 %>% inner_join(ep_seq, by="phenotype") 

hlines <- aggregate(Disease ~ icd_grp, res_xyy_p1, max)[,"Disease"] + 0.5
xlabs <- 0.5 * aggregate(Disease ~ icd_grp, res_xyy_p1, min)[,"Disease"] + 0.5 * aggregate(Disease ~ icd_grp, res_xyy_p1, max)[,"Disease"]
xlabs_text <- aggregate(Disease ~ icd_grp, res_xyy_p1, min)[,"icd_grp"]


## disease to highlight -----------
res_xyy_p1_label <- res_xyy_p1 %>% filter(al==1) %>% 
      filter(phenotype %in% c("Q17_CHROMOSOMAL_ABNORMALITI_NOT_ELSEW_CLASSIFIED", 
                              "J10_LOWCHRON", 
                              "ASTHMA_INFECTIONS", 
                              "E4_OBESITY", "E4_DM2", 
                              "AB1_ERYSIPELAS", "AB1_INTESTINAL_INFECTIONS", 
                               "I9_VARICVE", "I9_HEARTFAIL_NS", 
                               "E4_DM2",   # all dropped are different definitions of T2D 
                               "M13_KNEEDERANGEMENTS",  # soft tissue; need to update longname of M13_SCIATICA
                               "J10_LOWCHRON",
                               "E4_TESTICHYPO", "N14_MALEINFERT", "E4_OBESITY", "ASTHMA_PNEUMONIA", "PAIN", "ANY_REIMB", "K11_GULC", "F5_PSYCHDEV"))
dim(res_xyy_p1_label)   # 18 21
res_xyy_p1_label <- res_xyy_p1_label %>% mutate(LONGNAME=ifelse(phenotype=="PAIN", "Pain", LONGNAME), 
                                                LONGNAME=ifelse(phenotype=="ANY_REIMB", "Any reimbursement code in KELA", LONGNAME))


# tit <- paste0("Relationship of disease diagnoses with Jacobs syndrome in men")
tit <- paste0("Relationship of disease diagnoses with Jacobs syndrome")
pp_xyy <- ggplot(data=res_xyy_p1, aes(x=Disease, y=OR, label=LONGNAME)) + 
             geom_point(shape=25, aes(color=icd_grp, fill=icd_grp, alpha=al), size=2) + 
             geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.6) +
             geom_vline(xintercept=hlines, linetype="dotted", color="grey", size=0.6) + 
             labs(y="Odds ratio",x="",title=tit) + 
             scale_x_continuous(breaks=xlabs, label=xlabs_text) +
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
            # scale_y_log10(limits=c(0.1,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
             scale_y_log10(limits=c(0.5,111),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25,50,100)) + 
             scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
             scale_shape_manual(values=c(25))+ 
             theme_classic() + 
            # theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
             theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) + 
             theme(legend.position="none") + 
             theme(strip.background=element_blank(), strip.text.x=element_blank()) + 
             geom_text_repel(data=res_xyy_p1_label, size=3.5, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black", min.segment.length=0)
pp_xyy




###################################################################################################
###           Figure 2C: Manhattan plot for Klinefelter syndrome by disease category              #
###################################################################################################

## second panel for "ICD-10 main chapters" ------------------
res_xxy_p2 <- res_xxy_ep %>% filter(icd_grp %in% c("ICD-10 main chapters")) %>% 
      inner_join(icd_grpname, by="phenotype") 
dim(res_xxy_p2)   # 16   20

## order by icd group ------------------
res_xxy_p2$grpname <- factor(res_xxy_p2$grpname, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies"))

pp_xxy_category <- ggplot(data=res_xxy_p2, aes(x=grpname, y=OR, ymin=OR_025, ymax=OR_975 )) + 
		geom_pointrange(size=0.5, shape=22, aes(color=grpname, fill=grpname)) +
		geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) + 
		theme(legend.position="none") + 
		labs(y="Odds ratio", x="", title=paste0("")) +
		# scale_y_log10(limits=c(0.3,25), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_y_log10(limits=c(0.2,30), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() +
		# theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) +
		theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) +
		theme(legend.position="none")
pp_xxy_category




###################################################################################################
###           Figure 2D: Manhattan plot for Jacobs syndrome by disease category                   #
###################################################################################################

## second panel for "ICD-10 main chapters" ------------------
res_xyy_p2 <- res_xyy_ep %>% filter(icd_grp %in% c("ICD-10 main chapters")) %>% 
      inner_join(icd_grpname, by="phenotype") 
dim(res_xyy_p2)   # 16   20

## order by icd group ------------------
res_xyy_p2$grpname <- factor(res_xyy_p2$grpname, levels=c("Infectious & parasitic","Neoplasms", "Blood & immune mechanism",
                                                          "Endocrine, nutritional & metabolic", "Mental & behavioural", 
                                                          "Nervous system", "Eye & adnexa", "Ear & mastoid process", 
                                                          "Circulatory system", "Respiratory system", "Digestive system", 
                                                          "Skin & subcutaneous tissue", "Musculoskeletal system", 
                                                          "Genitourinary system", "Pregnancy, childbirth & puerperium", 
                                                          "Congenital anomalies"))

pp_xyy_category <- ggplot(data=res_xyy_p2, aes(x=grpname, y=OR, ymin=OR_025, ymax=OR_975 )) + 
		geom_pointrange(size=0.5, shape=22, aes(color=grpname, fill=grpname)) +
		geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) + 
		theme(legend.position="none") + 
		labs(y="Odds ratio", x="", title=paste0("")) +
		# scale_y_log10(limits=c(0.3,25), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_y_log10(limits=c(0.2,30), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() +
		# theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-70, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) +
		theme(axis.text.x=element_text(hjust=0, vjust=0, size=9, angle=-80, face="bold"), axis.title.y=element_text(size=12, face="bold"), plot.title=element_text(size=14, face="bold")) +
		theme(legend.position="none")
pp_xyy_category



###################################################################################################
###                                   Align XXX,XXY for Figure 2                                     #
###################################################################################################

## align plots of Trisomy X syndrome and Turner syndrome
p_xxx <- grid.arrange(pp_xxx, pp_xxx_category, ncol=2, widths=c(2.5,1), top=" ", bottom="", left=" ", right=" ")
p_xxy <- grid.arrange(pp_xxy, pp_xxy_category, ncol=2, widths=c(2.5,1), top=" ", bottom="", left=" ", right=" ")
p_xyy <- grid.arrange(pp_xyy, pp_xyy_category, ncol=2, widths=c(2.5,1), top=" ", bottom="", left=" ", right=" ")

ggsave(paste0("Figure1.Relationship_WithDisease_XXX.FinnGen.png"), p_xxx, width=13.5, height=8)
ggsave(paste0("Figure1.Relationship_WithDisease_XXY.FinnGen.png"), p_xxy, width=13.5, height=8)
ggsave(paste0("Figure1.Relationship_WithDisease_XYY.FinnGen.png"), p_xyy, width=13.5, height=8)



############################################################################################################################
###                 Figure 3B. XXX versus XXY (female with an extra versus male with an extra X)                           #
############################################################################################################################

## FDR for significance test --------
res_xxx_xxy <- inner_join(res_xxx_ep, res_xxy_ep %>% select(-icd_grp, -LONGNAME), by=c("phenotype")) %>% filter(al.x==1|al.y==1) %>% 
	mutate(zp=(beta.x-beta.y)/sqrt(se.x^2+se.y^2), 
	       P_val=2*pnorm(abs(zp),lower.tail=F), sig=ifelse(P_val<0.05/length(P_val),"T","F"),
	       FDR=p.adjust(P_val, method="fdr", n=length(P_val)),
	       sig.FDR=ifelse(FDR<0.05,"T","F"))
dim(res_xxx_xxy)   # 109  42

res_xxx_xxy_label <- res_xxx_xxy %>% filter(sig=="T")
# res_xxx_xxy_label <- res_xxx_xxy %>% filter(P_val<0.05)

## (women as x-axis and men as y-axis)
p_xxx_xxy <- ggplot(res_xxx_xxy, aes(x=OR.x, y=OR.y, label=LONGNAME)) + 
     # geom_point(shape=25, size=2, aes(alpha=al, color=grp_text_f, fill=grp_text_f)) + 
     geom_point(shape=25, size=2, aes(alpha=sig, color=icd_grp, fill=icd_grp)) + 
     geom_errorbar(aes(xmin=OR_025.x, xmax=OR_975.x, colour=icd_grp, alpha=sig), size=0.15, linetype=1) + 
     geom_errorbar(aes(ymin=OR_025.y, ymax=OR_975.y, colour=icd_grp, alpha=sig), size=0.15, linetype=1) + 
     geom_abline(intercept=0, slope=1, linetype="dotted", alpha=1, size=0.5, colour="grey") + 
     geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.35) +
     geom_vline(xintercept=1, linetype="dotted", color="grey", size=0.35) +
     labs(y="Odds ratio of disease diagnoses in men with an extra X", x="Odds ratio of disease diagnoses in women with an extra X", title="") +
     # xlim(-3.8,0.61) + ylim(-3.8,0.61) +
     scale_x_log10(limits=c(0.5,60), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,40), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,40))+
     scale_y_log10(limits=c(0.5,60), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,40), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,40))+
     scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
     scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
     theme_classic() + 
   #  theme(axis.text=element_text(size=8, hjust=0, vjust=0, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
     theme(axis.text=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
     theme(legend.position="none") + 
     # geom_text_repel(data=res_sex_label, size=2.8, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black")
     # geom_text_repel(data=subset(res_both, P_val<0.05), color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1)
  ##   geom_text_repel(data=res_sex_label, color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1, max.overlaps=Inf)
     geom_text_repel(data=res_xxx_xxy_label, color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1, max.overlaps=Inf, min.segment.length=0)
   # geom_text_repel(data=res,           color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, min.segment.length=0)
p_xxx_xxy

ggsave(paste0("Figure1.Relationship_WithDisease_Compare_XXX_XXY.FinnGen.png"), p_xxx_xxy, width=5, height=5)




############################################################################################################################
###                 Figure 3C. XYY versus XXY (male with an extra Y versus male with an extra X)                           #
############################################################################################################################

## FDR for significance test --------
res_xyy_xxy <- inner_join(res_xyy_ep, res_xxy_ep %>% select(-icd_grp, -LONGNAME), by=c("phenotype")) %>% filter(al.x==1|al.y==1) %>% 
	mutate(zp=(beta.x-beta.y)/sqrt(se.x^2+se.y^2), 
	       P_val=2*pnorm(abs(zp),lower.tail=F), sig=ifelse(P_val<0.05/length(P_val),"T","F"),
	       FDR=p.adjust(P_val, method="fdr", n=length(P_val)),
	       sig.FDR=ifelse(FDR<0.05,"T","F"))
dim(res_xyy_xxy)   # 109  42

res_xyy_xxy_label <- res_xyy_xxy %>% filter(sig=="T")
# res_xyy_xxy_label <- res_xyy_xxy %>% filter(P_val<0.05)

## (women as x-axis and men as y-axis)
p_xyy_xxy <- ggplot(res_xyy_xxy, aes(x=OR.x, y=OR.y, label=LONGNAME)) + 
     # geom_point(shape=25, size=2, aes(alpha=al, color=grp_text_f, fill=grp_text_f)) + 
     geom_point(shape=25, size=2, aes(alpha=sig, color=icd_grp, fill=icd_grp)) + 
     geom_errorbar(aes(xmin=OR_025.x, xmax=OR_975.x, colour=icd_grp, alpha=sig), size=0.15, linetype=1) + 
     geom_errorbar(aes(ymin=OR_025.y, ymax=OR_975.y, colour=icd_grp, alpha=sig), size=0.15, linetype=1) + 
     geom_abline(intercept=0, slope=1, linetype="dotted", alpha=1, size=0.5, colour="grey") + 
     geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.35) +
     geom_vline(xintercept=1, linetype="dotted", color="grey", size=0.35) +
     labs(y="Odds ratio of disease diagnoses in men with an extra X", x="Odds ratio of disease diagnoses in men with an extra Y", title="") +
     # xlim(-3.8,0.61) + ylim(-3.8,0.61) +
     scale_x_log10(limits=c(0.5,60), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,40), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,40))+
     scale_y_log10(limits=c(0.5,60), breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,40), labels=c(0.5,1,2,3,4,6,8,10,12,15,20,40))+
     scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
     scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
     theme_classic() + 
   #  theme(axis.text=element_text(size=8, hjust=0, vjust=0, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
     theme(axis.text=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold")) + 
     theme(legend.position="none") + 
     # geom_text_repel(data=res_sex_label, size=2.8, box.padding=0.65, point.padding=0.85, segment.size=0.36, max.overlaps=Inf, segment.color="grey", color="black")
     # geom_text_repel(data=subset(res_both, P_val<0.05), color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1)
  ##   geom_text_repel(data=res_sex_label, color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1, max.overlaps=Inf)
     geom_text_repel(data=res_xyy_xxy_label, color="black", size=2, box.padding=0.6, point.padding=NA, segment.size=0.26, force=1, max.overlaps=Inf, min.segment.length=0)
   # geom_text_repel(data=res,           color="black", size=2.8, box.padding=0.6, point.padding=NA, segment.size=0.26, min.segment.length=0)
p_xyy_xxy

ggsave(paste0("Figure1.Relationship_WithDisease_Compare_XYY_XXY.FinnGen.png"), p_xyy_xxy, width=5, height=5)







