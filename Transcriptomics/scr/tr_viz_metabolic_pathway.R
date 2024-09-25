library(tidyverse)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(circlize)

DEGs_wide <- read.csv("Transcriptomics/Output/DEGs_wide.csv")

Human.GEM <- read_csv("Transcriptomics/Data/External/Human-GEM.csv")
Translation_initiation <- read.delim("Transcriptomics/Data/External/GO_Translation_initiation.tsv")

gene_associations <- Human.GEM %>% 
  separate_rows("GENE ASSOCIATION", sep = "and|or") %>% 
  mutate(across(where(is.character), str_trim))

gene_associations$Symbol <- unname(mapIds(org.Hs.eg.db, keys=gene_associations$`GENE ASSOCIATION`, column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
gene_associations <- dplyr::select(gene_associations, Symbol, `GENE ASSOCIATION`, SUBSYSTEM) %>% 
drop_na() %>%   as.data.frame()

## Plot all pathways
pdf("Transcriptomics/Output/Metabolism.pdf")
for (i in unique(gene_associations$SUBSYSTEM)) {
  path <- DEGs_wide %>% 
    filter(gene_name %in% gene_associations$Symbol[gene_associations$SUBSYSTEM== i]) %>% 
    #  dplyr::select(contains("PI5Z"), gene_name) %>% 
    dplyr::select(contains("logFC"), gene_name) %>% 
    dplyr::select(-contains("1h"),-contains("_2h"),-contains("_4h")) %>% 
    column_to_rownames("gene_name") 
  
  colnames(path) <- gsub("logFC_","",colnames(path))
h <-  Heatmap(path, cluster_rows = F,cluster_columns = F, name = i)
print(h)
}
dev.off()

sum(unique(gene_associations$Symbol[gene_associations$SUBSYSTEM== "Nucleotide metabolism"]) %in% DEGs_wide$gene_name)
sum(unique(gene_associations$Symbol[gene_associations$SUBSYSTEM== "Aminoacyl-tRNA biosynthesis"]) %in% DEGs_wide$gene_name)

for(nr in c(10, 20)) {
  png(paste0("test_heatmap_nr_", nr, ".png"), width = 5, height = 0.1969*nr + 1.2744, 
      units = "in", res = 100)
  draw(Heatmap(random_mat(nr), height = unit(5, "mm")*nr, 
               column_title = "foo", # column title can be any one-line string
               top_annotation = HeatmapAnnotation(bar = 1:10)))
  dev.off()
}

plot_metabolic_heatmap <- function(path_name,timepoints_to_remove) {
  nr <-sum(unique(gene_associations$Symbol[gene_associations$SUBSYSTEM== path_name]) %in% DEGs_wide$gene_name)

  met <- DEGs_wide %>% 
    filter(gene_name %in% gene_associations$Symbol[gene_associations$SUBSYSTEM== path_name]) %>% 
    dplyr::select(contains("logFC"), gene_name) %>% 
    dplyr::select(-contains((timepoints_to_remove))) %>% 
    column_to_rownames("gene_name") 
  
  met[is.na(met)] <- 0
  colnames(met) <- gsub("logFC_","", colnames(met))
  
  values <-colnames(met)
  values_sorted <- values[order(sub("_.*", "", values), sub(".*_", "", sub("_.*", "", values)))]
  
  met <- met[values_sorted]
  sample_annot <- sub("_.*", "",  colnames(met))
  names(sample_annot) <- colnames(met)
  
  ha = HeatmapAnnotation(Drug= sample_annot,
                         col = list(Drug = c("PD" = "#E69F00", "PI" = "#009E73", "5Z" ="#CC79A7", "PIPD" = "#0072B2", "PI5Z" = "#D55E00")))
  
  col_fun = colorRamp2(c(-2, 0, 2), c("navyblue", "white", "red"))
  labels =sub(".*_", "",colnames(met))
  colnames(met) <- labels
  png(paste0(path_name, ".png"), res = 300, height = 10,width = 10, units = "in")
  
  h <- Heatmap(met, cluster_rows = T, show_row_dend = F, top_annotation = ha, column_split = sample_annot,
               cluster_columns = F,col = col_fun, name = "logFC", rect_gp = gpar(col= "black",lwd = 0.5),
       width = ncol(met)*unit(5, "mm"), 
       height = nrow(met)*unit(5, "mm"))
  print(h)
  dev.off()
}


plot_metabolic_heatmap("Nucleotide metabolism",c("1h","2h"))
plot_metabolic_heatmap("Aminoacyl-tRNA biosynthesis",c("1h","2h","_4h"))

prot_init <- DEGs_wide %>% 
  filter(gene_name %in% Translation_initiation$SYMBOL) %>% 
  dplyr::select(contains("logFC"), gene_name) %>% 
  dplyr::select(-contains("1h"),-contains("_2h"),-contains("_4h"),-contains("_5Z_"),-contains("_PI_"),-contains("_PD_")) %>% 
  column_to_rownames("gene_name") 

prot_init[is.na(prot_init)] <- 0
colnames(prot_init) <- gsub("logFC_","", colnames(prot_init))

values <-colnames(prot_init)
values_sorted <- values[order(sub("_.*", "", values), sub(".*_", "", sub("_.*", "", values)))]

prot_init <- prot_init[values_sorted]
sample_annot <- sub("_.*", "",  colnames(prot_init))
names(sample_annot) <- colnames(prot_init)
ha = HeatmapAnnotation(Drug= sample_annot,
                       col = list(Drug = c("PD" = "#E69F00", "PI" = "#009E73", "5Z" ="#CC79A7", "PIPD" = "#0072B2", "PI5Z" = "#D55E00")))

col_fun = colorRamp2(c(-2, 0, 2), c("navyblue", "white", "red"))
labels =sub(".*_", "",colnames(prot_init))
colnames(prot_init) <- labels

png("translation_initiation.png", res = 300, height = 10,width = 10, units = "in")

h <- Heatmap(prot_init, cluster_rows = T, show_row_dend = F, top_annotation = ha, column_split = sample_annot,
             cluster_columns = F,col = col_fun, name = "logFC", rect_gp = gpar(col= "black",lwd = 0.5),
             width = ncol(prot_init)*unit(5, "mm"), 
             height = nrow(prot_init)*unit(5, "mm"))
print(h)
dev.off()
