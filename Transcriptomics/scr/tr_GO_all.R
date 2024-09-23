library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)

PIPD_specific_DEGs <- read.csv("Transcriptomics/Output/PIPD_specific_DEGs.csv") 
PI5Z_specific_DEGs <- read.csv("Transcriptomics/Output/PI5Z_specific_DEGs.csv")
all_DEGs <- read.csv("Transcriptomics/Output/DEGs_long.csv")
all_genes <- read.csv("Transcriptomics/Data/gene_info.tsv", sep = "\t")

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

DEGs_per_condition <- all_DEGs %>% 
  split(~Contrast)

all_GO <- lapply(DEGs_per_condition, function(x) GO_clusters(x$gene_name,unique(all_genes$gene_name)))
all_GO_flt <- all_GO[sapply(all_GO, nrow)>0]
all_GO_flt <- lapply(all_GO_flt, function(x) x@result %>% 
                       filter(p.adjust < 0.01))

all_GO_flt_df <- do.call("rbind", all_GO_flt) %>% 
  rownames_to_column("rowID") %>% 
  mutate(Condition = sub("\\..*$", "", rowID))

all_GO_simplified <- lapply(all_GO, function(x) simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))
all_GO_sim_flt <- all_GO[sapply(all_GO_simplified, nrow)>0]
all_GO_sim_flt <- lapply(all_GO_sim_flt, function(x) x@result %>% 
                       filter(p.adjust < 0.01))

all_GO_sim_flt_df <- do.call("rbind", all_GO_sim_flt) %>% 
  rownames_to_column("rowID") %>% 
  mutate(Condition = sub("\\..*$", "", rowID))


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


compare_input <- lapply(DEGs_per_condition, function(x) x %>% dplyr::select(gene_name) %>% unlist()) 

ck <- compareCluster(geneCluster = compare_input, fun = enrichGO,  OrgDb= org.Hs.eg.db,  
                     keyType       = "SYMBOL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dotplot(ck)

ck_sim <- simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)

ck_sim@compareClusterResult$Cluster <- factor(ck_sim@compareClusterResult$Cluster, levels = c("Diff_PI_1h","Diff_PD_1h","Diff_PIPD_1h",
                                                                                              "Diff_OXO_2h", "Diff_PI_2h","Diff_PD_2h","Diff_PIPD_2h","Diff_PIOXO_2h",
                                                                                              "Diff_OXO_4h", "Diff_PI_4h","Diff_PD_4h","Diff_PIPD_4h","Diff_PIOXO_4h",
                                                                                              "Diff_OXO_8h","Diff_PD_8h","Diff_PIPD_8h","Diff_PIOXO_8h",
                                                                                              "Diff_OXO_24h", "Diff_PI_24h","Diff_PD_24h","Diff_PIPD_24h","Diff_PIOXO_24h"),
                                              labels = c("PI_1h","PD_1h","PIPD_1h",
                                                         "OXO_2h", "PI_2h","PD_2h","PIPD_2h","PIOXO_2h",
                                                         "OXO_4h", "PI_4h","PD_4h","PIPD_4h","PIOXO_4h",
                                                         "OXO_8h","PD_8h","PIPD_8h","PIOXO_8h",
                                                         "OXO_24h", "PI_24h","PD_24h","PIPD_24h","PIOXO_24h"))

ck_sim@compareClusterResult <- ck_sim@compareClusterResult %>% 
  filter(!grepl("development",Description)) %>% 
  filter(!grepl("organ",Description)) %>% 
  filter(!grepl("differentiation",Description)) %>% 
  filter(!grepl("skeletal",Description)) %>% 
  filter(!grepl("muscle",Description)) %>% 
  filter(!grepl("blood",Description)) %>% 
  filter(!grepl("lymphocyte",Description)) %>% 
  filter(!grepl("fibroblast",Description))  
  
  
pdf("Transcriptomics/Output/all_GO_terms_comparison.pdf", width = 20, height = 20)
dotplot(ck_sim)
dev.off()
