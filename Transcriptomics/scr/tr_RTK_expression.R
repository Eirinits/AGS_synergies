library(ComplexHeatmap)
library(tidyverse)

RTKs <- read.table("Transcriptomics/Data/External/RTKs", quote="\"", comment.char="")
pi3k <- c("AKT","AKT1S1", "DEPTOR", "MLST8", "RPTOR", "MTOR","EEF1A1", "DEPTOR", "MLST8", "PRR5", "RICTOR", 
          "AKT1", "AKT2", "AKT3", "AMPK", "BRAF", "CREB1", "EIF4E", "EIF4EBP1", "ERK1/2", "FOXO", "GAB1", "GF", "GRB2", "GSK3B", "GYS1", "HRAS", "INS", "INSR", "IRS1", "LAMTOR", "MEK1/2", "MKNK1", "MTOR", "mTORC1", "mTORC2", "PDPK1", "PI3K", "PIK3CA", "PIK3R1", "PIP3", "PPARGC1A", "PPP2CA", "PTEN", "RAGAC", "RHEB", "RPS6K", "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KB1", "RTKs", "SOS1", "SREBF1", "TSC", "ULK1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG2", "PRKAG3", "PRKAA1", "PRKAG1", "FOXO1", "FOXO3", "FOXO4", "FOXO6", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")

RTK_logFC <- DEGs_wide %>% 
  filter(gene_name %in% RTKs$V1) %>% 
  column_to_rownames("gene_name") 

RTK_logFC[is.na(RTK_logFC)] <- 0

Heatmap(as.matrix(RTK_logFC[,1:25]), cluster_columns = F)


RTK_logFC <- DEGs_wide %>% 
  filter(gene_name %in% pi3k) %>% 
  column_to_rownames("gene_name") 

RTK_logFC[is.na(RTK_logFC)] <- 0

Heatmap(as.matrix(RTK_logFC[,1:25]), cluster_columns = F)
