library(tidyverse)
library(magrittr)
library(clusterProfiler)

library(org.Hs.eg.db)

PIPD_specific_DEGs <- read.csv("Transcriptomics/Output/PIPD_specific_DEGs.csv") 
PI5Z_specific_DEGs <- read.csv("Transcriptomics/Output/PI5Z_specific_DEGs.csv")
background_DEGs <- read.csv("Transcriptomics/Output/DEGs_long.csv") %>% 
  dplyr::select(gene_name) %>% distinct(gene_name)

GO_clusters <- function(genes,background_genes){
  tryCatch(enrichGO(gene      = genes,
                    universe      = background_genes,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE), error=function(e) NULL)
}

PIPD <- PIPD_specific_DEGs %>% 
  split(~Condition)

PIPD <- lapply(PIPD, function(x) x %>% distinct(gene_name))
PIPD_GO <- lapply(PIPD, function(x) GO_clusters(x$gene_name,background_DEGs$gene_name))

PIPD_8 <- PIPD_GO[["PIPD_8h"]]@result %>% 
  filter(p.adjust < 0.05)
PIPD_24 <- PIPD_GO[["PIPD_24h"]]@result %>% 
  filter(p.adjust < 0.05)

PI5Z <- PI5Z_specific_DEGs %>% 
  filter(grepl("PI5Z",Condition)) %>% 
  split(~Condition)

PI5Z <- lapply(PI5Z, function(x) x %>% distinct(gene_name))
PI5Z_GO <- lapply(PI5Z, function(x) GO_clusters(x$gene_name,background_DEGs$gene_name))

PI5Z_8 <- PI5Z_GO[["PI5Z_8h"]]@result %>% 
  filter(p.adjust < 0.05)

PI5Z_24 <- PI5Z_GO[["PI5Z_24h"]]@result %>% 
  filter(p.adjust < 0.05)

gene_list <- unique(unlist(strsplit(c(PI5Z_8$geneID, PIPD_8$geneID), "/")))
gene_list

RNA_PIPD <- PIPD_specific_DEGs %>% filter(gene_name %in% gene_list)
RNA_PI5Z <- PI5Z_specific_DEGs %>% filter(gene_name %in% gene_list)
