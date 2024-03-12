library(tidyverse)

Kinase_Substrate_Dataset <- read.delim("~/Documents/Git/AGS_synergies_multiomics/Phosphoproteomics/Data/External/Kinase_Substrate_Dataset", row.names=NULL)
DPPs_long_fdr005 <- read.csv("~/Documents/Git/AGS_synergies_multiomics/Phosphoproteomics/Data/DPPs_long_fdr005.csv")

Kinase <- "CAMK2A"

Int_flt <- Kinase_Substrate_Dataset %>% 
  filter(row.names == Kinase) %>% 
  unite(PTM, SUB_GENE,SUB_MOD_RSD, sep = "_")

sites <- DPPs_long_fdr005 %>% 
  separate(PTM_collapse_key, c("Prot","Site","M"), sep = "_") %>% 
  unite(PTM, Prot,Site,remove = F)

intersect(Int_flt$PTM, sites$PTM)


Reg_sites <- Regulatory_sites %>% 
  filter(ORGANISM == "human" & grepl("-p",MOD_RSD)) %>% 
  mutate(MOD_RSD = gsub("-p","",MOD_RSD)) %>% 
  unite(PTM, GENE, MOD_RSD ) %>% 
  select(-PROTEIN,-ACC_ID,-GENE_ID,-HU_CHR_LOC,-ORGANISM,-SITE_GRP_ID,-SITE_...7_AA) 

PIPD_DPPs <- vapply(strsplit(PIPD_specific_DPPs$x, '_'), function(x)
  paste(x[seq.int(2)], collapse='_'), character(1L))

All_ratios_Limma$PTM <- vapply(strsplit(All_ratios_Limma$PTM_collapse_key, '_'), function(x)
  paste(x[seq.int(2)], collapse='_'), character(1L)) 

PIPD_reg_sites <- All_ratios_Limma %>% 
  right_join(Reg_sites) %>% 
  filter(PTM %in%PIPD_DPPs) %>% 
  select(colnames(Reg_sites),contains("PIPD")) %>% 
  mutate(logFC_PIPD_30m = ifelse(adj.P.Val_PIPD_30m > 0.05,NA,logFC_PIPD_30m)) %>% 
  mutate(logFC_PIPD_2h = ifelse(adj.P.Val_PIPD_2h > 0.05,NA,logFC_PIPD_2h)) %>% 
  mutate(logFC_PIPD_8h = ifelse(adj.P.Val_PIPD_8h > 0.05,NA,logFC_PIPD_8h)) %>% 
  select(-contains("adj"))

