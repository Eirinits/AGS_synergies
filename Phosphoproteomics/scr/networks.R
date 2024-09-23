library(tidyverse)

All_ratios_Limma <- read.csv("Phosphoproteomics/Data/All_ratios_Limma.txt")
Regulatory_sites <- read.delim("Phosphoproteomics/Data/External/Regulatory_sites")
PIPD_specific <- read.csv("Phosphoproteomics/Output/PIPD_specific_DPPs.csv")
PI5Z_specific <- read.csv("Phosphoproteomics/Output/PI5Z_specific_DPPs.csv")

filter_nan <- function(df){
  n_nas <- rowSums(is.na(df[,1:3]))
  valid_rows <- which(n_nas <= 2)
  return(df[valid_rows,])
}

Reg_sites <- Regulatory_sites %>% 
  filter(ORGANISM == "human" & grepl("-p",MOD_RSD)) %>%
  mutate(MOD_RSD = gsub("-p","",MOD_RSD)) %>% 
  unite(PTM, GENE, MOD_RSD) %>% 
  select(PTM, contains("ON")) %>% 
  mutate(Main_effect =ifelse(grepl("activity, induced",ON_FUNCTION) & !grepl("activity, inhibited",ON_FUNCTION), "Activates",
                             ifelse(grepl("activity, inhibited",ON_FUNCTION) & !grepl("activity, induced",ON_FUNCTION), "Inhibits",
                                    ifelse(grepl("activity, induced",ON_FUNCTION) & grepl("activity, inhibited",ON_FUNCTION), "Both",NA)))) %>% 
  mutate(Main_effect = ifelse(is.na(Main_effect) & grepl("molecular", ON_FUNCTION),"Other",
                              ifelse(is.na(Main_effect) & grepl("localization", ON_FUNCTION),"Other",Main_effect)))

nodes <- c("AKT1", "AKT1S1", "ARAF", "Apoptosis", "Autophagy", "BRAF", "CASP2", "CASP4", "CASP7", 
           "CSKN2A1", "CTNNA1", "DUSP16", "EIF2S2", "EIF4E", "EIF4EBP1", "FOXO1", "FOXO3", "GAB1", 
           "GAB2", "GSK3A", "GSK3B", "GYS1", "Glycogen", "IKBKB", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", 
           "MDM2", "MKNK1", "MTOR", "NFKB1", "PDPK1", "PI3K", "PI3KCA", "PIK3R1", "Proliferation", 
           "Protein", "RAF1", "RICTOR", "RPS6KA1", "RPS6KB1", "RPTOR", "TNFAIP1", "USP7", "YAP1", 
           "mTORC1", "mTORC2","MAP3K7","MAP2K3","MAP2K4","MAPK14","RPS6KA4")

PI5Z <- All_ratios_Limma %>% 
  select(contains("PI5Z") & contains("og"),PTM_collapse_key) %>%
  filter_nan(.) %>% 
  mutate(Specific = ifelse(PTM_collapse_key %in% PI5Z_specific$PTM_collapse_key,1,0)) %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"), sep = "_", remove = F) %>% 
  unite(Site_M , Site,M, remove = F) %>% 
  unite(PTM, Protein, Site, remove = F) %>% 
  filter(Protein %in% nodes) %>% 
  left_join(Reg_sites)  %>% 
  select(-contains("ON"))

write.table(PI5Z, file = "Phosphoproteomics/Output/Response_networks/PI5Z_psites_char.csv", sep = ",",row.names = F)

PIPD <- All_ratios_Limma %>% 
  select(contains("PIPD") & contains("og"),PTM_collapse_key) %>%
  filter_nan(.)%>% 
  mutate(Specific = ifelse(PTM_collapse_key %in% PIPD_specific$PTM_collapse_key,1,0))%>% 
  separate(PTM_collapse_key, c("Protein","Site","M"), sep = "_", remove = F) %>% 
  unite(Site_M , Site,M, remove = F) %>% 
  unite(PTM, Protein, Site, remove = F) %>% 
  filter(Protein %in% nodes) %>% 
  left_join(Reg_sites)%>% 
  select(-contains("ON"))

write.table(PIPD, file = "Phosphoproteomics/Output/Response_networks/PIPD_psites_char.csv", sep = ",",row.names = F)


