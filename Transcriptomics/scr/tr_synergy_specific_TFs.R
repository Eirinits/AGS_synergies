library(tidyverse)
library(magrittr)
library(decoupleR)
library(clusterProfiler)
library(org.Mm.eg.db)

TCF2_categories <- read.delim("Transcriptomics/Data/External/TCF2_categories.csv")
psite_metadata <- read.delim("Phosphoproteomics/Data/Raw_data/psite_metadata.csv", header=T, sep = " ") 
CollecTRI <- read.csv("Transcriptomics/Data/External/signed_CollecTRI.csv")
PIPD_specific_DEGs <- read.csv("Transcriptomics/Output/PIPD_specific_DEGs.csv") 
PI5Z_specific_DEGs <- read.csv("Transcriptomics/Output/PI5Z_specific_DEGs.csv")

PIPD_specific_DPPs <- read.csv("Phosphoproteomics/Output/PIPD_specific_DPPs.csv") %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"), sep = "_", remove = F) %>% 
  left_join(psite_metadata) 

ROKAI_PIPD <- select(PIPD_specific_DPPs, Protein,Site,logFC)
colnames(ROKAI_PIPD) <- c("Protein","Position","Quantification")

#write.table(ROKAI_PIPD, file = "Phosphoproteomics/Output/ROKAI_inp_PIPD.csv",sep = ",", row.names = F)

PI5Z_specific_DPPs <- read.csv("Phosphoproteomics/Output/PI5Z_specific_DPPs.csv") %>% 
  separate(PTM_collapse_key, c("Protein","Site","M"), sep = "_", remove = F) %>% 
  left_join(psite_metadata)

ROKAI_PI5Z <- select(PI5Z_specific_DPPs, PTM_collapse_key,Protein,logFC, Prot_description)

CollecTRI <- read.csv("Transcriptomics/Data/External/signed_CollecTRI.csv")

intersect(TCF2_categories$Associated.Gene.Name,PIPD_specific_DEGs$gene_name)
intersect(TCF2_categories$Associated.Gene.Name,PI5Z_specific_DEGs$gene_name)
intersect(TCF2_categories$Associated.Gene.Name,PIPD_specific_DPPs$Protein)
intersect(TCF2_categories$Associated.Gene.Name,PI5Z_specific_DPPs$Protein)


run_decoupler <- function(mat,network){
  tryCatch(run_ulm(as.matrix(mat), 
                   network = network, 
                   .source='source', 
                   .target='target',
                   minsize = 5), error=function(e) NULL)
  
}

PI5Z_DEGs_split <- PI5Z_specific_DEGs %>% 
  split(~ Time)

PI5Z_DEGs_split <-lapply(PI5Z_DEGs_split, function(x) x[!duplicated(x$gene_name),])
PI5Z_DEGs_split <- lapply(PI5Z_DEGs_split, function(x) as.data.frame(x) %>% set_rownames(.$gene_name) %>% dplyr::select(logFC))

res_decoupler_PI5Z <- lapply(PI5Z_DEGs_split, function(x) run_decoupler(x,CollecTRI))

res_decoupler_PI5Z <-  res_decoupler_PI5Z %>% discard(is.null)
res_decoupler_PI5Z <- lapply(res_decoupler_PI5Z, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.1))
res_decoupler_PI5Z <- res_decoupler_PI5Z[sapply(res_decoupler_PI5Z, nrow)>0]

res_decoupler_PI5Z_flt <- Map(cbind, res_decoupler_PI5Z, Condition = names(res_decoupler_PI5Z))

res_decoupler_PI5Z_flt <- res_decoupler_PI5Z_flt  %>%
  purrr::reduce(rbind) 

write.table(res_decoupler_PI5Z_flt, file = "Transcriptomics/Output/PI5Z_specific_TFs.csv", sep = ",")

PIPD_DEGs_split <- PIPD_specific_DEGs %>% 
  split(~ Time)

PIPD_DEGs_split <-lapply(PIPD_DEGs_split, function(x) x[!duplicated(x$gene_name),])
PIPD_DEGs_split <- lapply(PIPD_DEGs_split, function(x) as.data.frame(x) %>% set_rownames(.$gene_name) %>% dplyr::select(logFC))

res_decoupler_PIPD <- lapply(PIPD_DEGs_split, function(x) run_decoupler(x,CollecTRI))

res_decoupler_PIPD_flt <-  res_decoupler_PIPD %>% discard(is.null)
res_decoupler_PIPD_flt <- lapply(res_decoupler_PIPD_flt, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.1))
res_decoupler_PIPD_flt <- res_decoupler_PIPD_flt[sapply(res_decoupler_PIPD_flt, nrow)>0]

res_decoupler_PIPD_flt <- Map(cbind, res_decoupler_PIPD_flt, Condition = names(res_decoupler_PIPD_flt))

res_decoupler_PIPD_flt <- res_decoupler_PIPD_flt  %>%
  purrr::reduce(rbind) 

write.table(res_decoupler_PIPD_flt, file = "Transcriptomics/Output/PIPD_specific_TFs.csv", sep = ",")

intersect(res_decoupler_PIPD_flt$source, PIPD_specific_DPPs$Protein)
