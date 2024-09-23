library(tidyverse)

counts<-read.table("Transcriptomics/Data/gene_counts_flt_norm_wCTRL.csv", sep = ",")
gene_info <- read.table("Transcriptomics/Data/gene_info.tsv", sep = "\t") %>% 
  rownames_to_column("Gene_name")

DPPs_long_fdr005 <- read.csv("Phosphoproteomics/Data/DPPs_long_fdr005.csv")

Regulatory_sites_formatted <- read.csv("Phosphoproteomics/Data/External/Regulatory_sites_formatted.csv")%>% 
  separate(PTM, c("Protein","Site"),sep = "_", remove = F)

TF_sites <- Regulatory_sites_formatted %>% 
  filter(grepl("transcription",ON_PROCESS)|grepl("DNA",ON_OTHER_INTERACT)) %>% 
  separate(PTM, c("Protein","Site"),sep = "_", remove = F)

drug <- "PI5Z"

DPPs_long_fdr005_PIPD <- DPPs_long_fdr005 %>% 
  separate(PTM_collapse_key, c("Prot","Site","M"), sep = "_", remove = F) %>% 
  unite(PTM, Prot, Site, sep = "_",remove = F) %>% 
  filter(grepl(drug, Condition)) %>% 
  select(-p_val) %>% 
  spread(., Condition, LogFC)

gene_counts_log_PIPD <- log2(counts+1) %>% 
  as.data.frame(.) %>% 
  dplyr::select(contains(drug), -contains("DMSO"),-contains("untreated")) %>% 
  rownames_to_column("Gene_name")  %>% 
#  filter(Gene_name %in% rownames(prot_coding)) %>% 
#  filter(Gene_name %in% PIPD_degs$Gene_name) %>% 
  gather(., Condition, counts, 2:ncol(.)) %>% 
  separate(Condition, c("Treatment","Time","Rep"), sep = "_") %>% 
  group_by(Gene_name,Treatment,Time) %>% 
  summarise(mean_counts = mean(counts)) %>% 
  ungroup() %>% 
  unite(Condition, Treatment, Time, sep = "_") %>% 
  spread(Condition, mean_counts) %>% 
  column_to_rownames("Gene_name") %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame(.) %>% 
  rownames_to_column("Gene_name") %>% 
  left_join(gene_info)

unique(intersect(TF_sites$PTM, DPPs_long_fdr005_PIPD$PTM))

load("~/Documents/Git/AGS_synergies_multiomics/Transcriptomics/Output/res_decoupler_all_DEGs.RData")
res_decoupler_flt <-  res_decoupler %>% discard(is.null)
res_decoupler_flt <- lapply(res_decoupler_flt, function(x) x %>% dplyr::filter(statistic == "consensus" & p_value < 0.05 & source %in% c(gene_info$gene_name,"NFKB","AP1") & condition =="V1"))
res_decoupler_flt <- res_decoupler_flt[sapply(res_decoupler_flt, nrow)>0]

res_decoupler_flt <- Map(cbind, res_decoupler_flt, Condition = names(res_decoupler_flt)) %>%
  purrr::reduce(rbind) %>% 
  mutate(Condition = gsub("Diff_","",Condition)) %>% 
  separate(Condition,c("Drug","Time"), sep = "_", remove = F) %>% 
  rename(pval = "p_value") %>% 
  rename(Prot = "source") %>% 
  #  mutate(score = round(score,2))%>% 
  select(Prot,Condition,score,pval) %>% 
  filter(grepl(drug,Condition)) %>% 
  select(-pval) %>%
  spread(., Condition, score)


annotated_DPPs <- left_join(DPPs_long_fdr005_PIPD, Regulatory_sites_formatted, by = "PTM", relationship ="many-to-many")
DPPs_and_activities <- inner_join(annotated_DPPs, res_decoupler_flt, by ="Prot")





AKT1 <- DPPs_long_fdr005 %>% 
  filter(grepl("AKT1_", PTM_collapse_key)) %>% 
  filter(grepl("PI5Z", Condition)) %>% 
  as.data.frame() %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"),remove = F) %>% 
  unite(PTM, Protein,Site, remove = F)

specific <- intersect(AKT1$PTM_collapse_key, PI5Z_specific_DPPs$PTM_collapse_key)
specific
annot <- intersect(AKT1$PTM, Regulatory_sites_formatted$PTM)
annot

annotations <- Regulatory_sites_formatted %>% 
  filter(PTM %in% annot)

dpps <- All_ratios_Limma %>% 
  filter(PTM_collapse_key %in% AKT1$PTM_collapse_key) %>% 
  select(contains("PI5Z"), PTM_collapse_key) %>%
  mutate(logFC_PI5Z_30m = ifelse(adj.P.Val_PI5Z_30m < 0.05,logFC_PI5Z_30m, NA)) %>% 
  mutate(logFC_PI5Z_2h = ifelse(adj.P.Val_PI5Z_2h < 0.05,logFC_PI5Z_2h, NA)) %>% 
  mutate(logFC_PI5Z_8h = ifelse(adj.P.Val_PI5Z_8h < 0.05,logFC_PI5Z_8h, NA)) %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"),remove = F) %>% 
  unite(PTM, Protein,Site, remove = F) %>% 
  filter(PTM %in% annot)

AKT1_df <- data.frame("PTM" = annot, "Target" = "AKT1") 
net <- left_join(AKT1_df, dpps) %>% 
  filter(Protein == "AKT1") 

write.table(net,"AKT1_net.csv",quote = F, row.names = F)








syn_spec <- PIPD_specific_DPPs %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"),remove = F, sep = "_") %>% 
  unite(PTM, Protein,Site, remove = F)

specific <- intersect(syn_spec$PTM, Regulatory_sites_formatted$PTM)
specific

dpps <- DPPs_long_fdr005 %>% 
  #  filter(grepl("FOXO3_", PTM_collapse_key)) %>% 
  filter(grepl("PIPD", Condition)) %>% 
  as.data.frame() %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"),remove = F,sep = "_") %>% 
  unite(PTM, Protein,Site, remove = F)

specific <- intersect(dpps$PTM_collapse_key, PIPD_specific_DPPs$PTM_collapse_key)
specific

annot <- intersect(specific, Regulatory_sites_formatted$PTM)
annot

annotations <- Regulatory_sites_formatted %>% 
  filter(PTM %in% syne)



dpps <- All_ratios_Limma %>% 
  select(contains("PI5Z"), PTM_collapse_key) %>%
  mutate(logFC_PI5Z_30m = ifelse(adj.P.Val_PI5Z_30m < 0.05,logFC_PI5Z_30m, NA)) %>% 
  mutate(logFC_PI5Z_2h = ifelse(adj.P.Val_PI5Z_2h < 0.05,logFC_PI5Z_2h, NA)) %>% 
  mutate(logFC_PI5Z_8h = ifelse(adj.P.Val_PI5Z_8h < 0.05,logFC_PI5Z_8h, NA)) %>% 
  select(-contains("adj")) %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"),remove = F,sep = "_") %>% 
  unite(PTM, Protein,Site, remove = F) %>% 
  unite(Site_M, Site, M, remove = F)

dpps <-  dpps[rowSums(is.na(dpps)) != 3, ]

FOXO3_df <- data.frame("PTM" = annot, "Target" = "FOXO3") 
net <- left_join(FOXO3_df, dpps) %>% 
  filter(Protein == "FOXO3") 

write.table(net,"FOXO3_net.csv",quote = F, row.names = F)

write.table(dpps, file = "PI5Z_psites_for_viz.csv", row.names = F, quote = F)

