library(tidyverse)
library(decoupleR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(magrittr)
library(biomaRt)
library(gridExtra)

cluster_names <- read.delim("Transcriptomics/Output/Clustering/tr_clusters_names.csv",sep = " ")
gene_info <- read.delim("Transcriptomics/Data/gene_info.tsv")
CollecTRI <- read.csv("Transcriptomics/Data/External/signed_CollecTRI.csv")
DEGs_long <- read.csv("Transcriptomics/Output/DEGs_long.csv")

temp = list.files(pattern="__clusters",path = "Transcriptomics/Output/Clustering", full.names = T)
list2env(
  lapply(setNames(temp, make.names(gsub("Transcriptomics/Output/Clustering|.csv", "", temp))), 
         read.csv), envir = .GlobalEnv)


split_clusters_per_timepoint <- function(cluster_df, drug_name){
  
  cluster_df <- cluster_df %>% 
    dplyr::select(Row.names, hclust_analysis) 
  
  DEGs <- DEGs_long %>% 
    dplyr::filter(grepl(paste0("_",drug_name,"_"), Contrast)) %>% 
    mutate(Contrast = gsub("Diff_","",Contrast)) %>% 
    filter(gene_name != "") %>% 
    left_join(cluster_df, by = c("Gene_name" = "Row.names")) %>% 
    unite(Split_fct, Contrast, hclust_analysis) %>% 
    dplyr::select(gene_name,t,Split_fct)  %>% 
    split(~Split_fct) 
  
  var_name <- paste0(drug_name,"_clusters")
  assign(var_name, DEGs, env=.GlobalEnv)
  }

split_clusters_per_timepoint(X.OXO__clusters,"OXO")
split_clusters_per_timepoint(X.PD__clusters,"PD")
split_clusters_per_timepoint(X.PI__clusters,"PI")
split_clusters_per_timepoint(X.PIPD__clusters,"PIPD")
split_clusters_per_timepoint(X.PIOXO__clusters,"PIOXO")

all_clusters <- c(PI_clusters,PD_clusters,OXO_clusters,PIOXO_clusters,PIPD_clusters)


# decoupleR for TF activities ---------------------------------------------

all_clusters_flt <- lapply(all_clusters, function(x) as.data.frame(x) %>% set_rownames(.$gene_name) %>% dplyr::select(t))
all_clusters_flt <- all_clusters_flt[sapply(all_clusters_flt, nrow)>5]

run_decoupler <- function(mat,network){
  tryCatch(decouple(as.matrix(mat),network = network,.source='source', 
                    .target='target', args = list(wsum = list(times = 1000)),
                    minsize = 5), error=function(e) NULL)
  
}


#res_decoupler <- lapply(all_clusters_flt, function(x) run_decoupler(x,CollecTRI))
load("Transcriptomics/Output/Clustering/res_decoupler_clusters.RData")
res_decoupler_flt <-  res_decoupler %>% discard(is.null)
res_decoupler_flt <- lapply(res_decoupler_flt, function(x) x %>% dplyr::filter(statistic == "consensus" & p_value < 0.05 & condition == "t"))
res_decoupler_flt <- res_decoupler_flt[sapply(res_decoupler_flt, nrow)>0]
#save(res_decoupler_flt, file = "Transcriptomics/Output/Clustering/res_decoupler_flt.RData")
res_decoupler_list <- Map(cbind, res_decoupler_flt, Condition = names(res_decoupler_flt))

res_decoupler_df <- res_decoupler_list  %>%
  purrr::reduce(rbind) %>% 
  separate(Condition, into = c("Drug","Time","Cluster"), remove = F) %>% 
  unite(Drug_cl, Drug, Cluster, remove = F)

pdf("Transcriptomics/Output/Clustering/TF_activity_clusters_consensus.pdf", height = 21)
ggplot(data=res_decoupler_df,aes( y = source, x=Time, color = score, size = -log10(p_value))) +
  geom_point() +
  theme_minimal() + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                         high="red", space ="Lab" ) +
  ylab("-log10(q-value)") + 
  facet_wrap(~Drug_cl, scales = "free", nrow = 20)
  
dev.off() 

res_decoupler_per_cluster <- res_decoupler_df %>% 
  split(~ Drug_cl)

save(res_decoupler, file = "Transcriptomics/Output/Clustering/res_decoupler_clusters.RData")
# TF targets --------------------------------------------------------------
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

split_genes_per_cluster <- function(cluster_df,drug_name){
  
  genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                 values = cluster_df$Row.names, mart= mart)

  cluster_df <- cluster_df %>%
    dplyr::select(Row.names, hclust_analysis) %>%
    full_join(genes, by =c("Row.names" = "ensembl_gene_id")) %>%
    filter(hgnc_symbol != "") %>%
    split(~hclust_analysis) %>%
    setNames(., paste0(drug_name, "_",names(.)))
  
  var_name <- paste0(drug_name,"_clusters_per_drug")
  assign(var_name, cluster_df, env=.GlobalEnv)

}

split_genes_per_cluster(X.OXO__clusters,"OXO")
split_genes_per_cluster(X.PI__clusters,"PI")
split_genes_per_cluster(X.PD__clusters,"PD")
split_genes_per_cluster(X.PIOXO__clusters,"PIOXO")
split_genes_per_cluster(X.PIPD__clusters,"PIPD")

all_clusters_splitted <- c(PI_clusters_per_drug,
                           PD_clusters_per_drug,
                           OXO_clusters_per_drug,
                           PIOXO_clusters_per_drug,
                           PIPD_clusters_per_drug)

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

gene_info <- read.delim("Transcriptomics/Data/gene_info.tsv")

#enriched_GO_cl <- lapply(all_clusters_splitted, function(x) GO_clusters(x$hgnc_symbol,gene_info$gene_name))
#save(enriched_GO_cl, file = "Transcriptomics/Output/Clustering/enriched_GO_clusters.RData")
load("Transcriptomics/Output/Clustering/enriched_GO_clusters.RData")

enriched_GO_cl_flt <-  enriched_GO_cl %>% discard(is.null)
enriched_GO_cl_flt <- lapply(enriched_GO_cl_flt, function(x) x@result %>% dplyr::filter(p.adjust < 0.01))
enriched_GO_cl_flt <- enriched_GO_cl_flt[sapply(enriched_GO_cl_flt, nrow)>0]

enriched_GO_cl_flt <- map2(enriched_GO_cl_flt, names(enriched_GO_cl_flt), ~ mutate(.x, Condition = .y))
#save(enriched_GO_cl_flt, file = "Transcriptomics/Output/Clustering/enriched_GO_clusters_flt.RData")
load("Transcriptomics/Output/Clustering/enriched_GO_clusters_flt.RData")

enriched_GO_df <- do.call("rbind", enriched_GO_cl_flt)

#write.table(enriched_GO_df, file = "Transcriptomics/Output/enrich_GO_clusters.csv", row.names = F)

enriched_GO_df$GeneRatio <- eval(parse(text = enriched_GO_df$GeneRatio))

ggplot(data=enriched_GO_df[1:200,],aes( x=reorder(Description,qvalue), y = -log10(qvalue))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  facet_wrap(~Condition, scales = "free", nrow = 20)

simplifyGO_clusters <- function(res){
  tryCatch(simplify(res), error=function(e) NULL)
}
simplified_GOs <- lapply(enriched_GO_cl, function(x) simplifyGO_clusters(x))

simplified_GOs_flt <-  simplified_GOs %>% discard(is.null)
simplified_GOs_flt <- lapply(simplified_GOs_flt, function(x) x@result %>% dplyr::filter(p.adjust < 0.01))
simplified_GOs_flt <- simplified_GOs_flt[sapply(simplified_GOs_flt, nrow)>0]

simplified_GOs_flt <- map2(simplified_GOs_flt, names(simplified_GOs_flt), ~ mutate(.x, Condition = .y))

save(simplified_GOs_flt, file = "Transcriptomics/Output/Clustering/simplified_GOs_clusters_flt.RData")

simplified_GO_df <- do.call("rbind", simplified_GOs_flt) %>% 
  left_join(cluster_names, by = c("Condition"="Cluster_ID")) %>% 
  dplyr::select(Description,GeneRatio,Full_name,Condition)

enriched_GO_df$GeneRatio <- eval(parse(text = enriched_GO_df$GeneRatio))

ggplot(data=enriched_GO_df[1:10,],aes( x=reorder(Description,qvalue), y = -log10(qvalue))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  facet_wrap(~Condition, scales = "free", nrow = 20)

# TF target overrepresentation --------------------------------------------

run_enrich_clusters <- function(genes,network,background_genes){
  tryCatch(enricher(genes,
                    TERM2GENE = network,
                    minGSSize = 10,
                    universe = background_genes), error=function(e) NULL)
}

enrich_TFs_cl <- lapply(all_clusters_splitted, function(x) run_enrich_clusters(x$hgnc_symbol,CollecTRI, gene_info$gene_name))
#save(enrich_TFs_cl, file = "Transcriptomics/Output/Clustering/enriced_targets_clusters.RData")
load("Transcriptomics/Output/Clustering/enriced_targets_clusters.RData")

enrich_TFs_cl_flt <-  enrich_TFs_cl %>% discard(is.null)
enrich_TFs_cl_flt <- lapply(enrich_TFs_cl_flt, function(x) x@result %>% dplyr::filter(p.adjust < 0.05))
enrich_TFs_cl_flt <- enrich_TFs_cl_flt[sapply(enrich_TFs_cl_flt, nrow)>0]
enrich_TFs_cl_flt <- map2(enrich_TFs_cl_flt, names(enrich_TFs_cl_flt), ~ mutate(.x, Condition = .y))

enriched_TFs_df <- do.call("rbind", enrich_TFs_cl_flt) %>% 
  left_join(cluster_names, by = c("Condition"="Cluster_ID")) %>% 
  dplyr::select(ID, Condition,Direction,GeneRatio, Full_name, Diff_than_DMSO, geneID) %>% 
  filter(!grepl("OXO",Condition))

  

ggplot(data=enriched_TFs_df, aes(x=reorder(ID, qvalue, decreasing = T), y=-log10(qvalue))) +
  geom_bar(stat="identity") +
  coord_flip() + theme_minimal() + 
  ylab("-log10(q-value)") + 
  xlab("TF") + 
  facet_wrap(~Condition, scales = "free")
# Visualize all cluster characterization ----------------------------------

temp = list.files(path ="Transcriptomics/Output/Clustering/", pattern="*.csv", full.names = T)
cluster_files = lapply(temp, read.csv)
cluster_files = setNames(cluster_files, temp)
names(cluster_files) <- gsub("Transcriptomics/Output/Clustering//|__clusters.csv","", names(cluster_files))

for (i in 1:length(cluster_files)) {
  centroid <- all_centroids[[i]]
  cluster <- cluster_files[[i]]
  colnames(centroid) <- colnames(cluster)[3:14]
  for (j in 1:12) {
    centroid[,j]<-tapply(as.matrix(cluster[,j+2]), as.factor(cluster$hclust_analysis), function(x) median(x, na.rm=T))
    var_name <- paste0("centroids_",names(all_centroids[i]))
    assign(var_name, centroid, env=.GlobalEnv)
  }
}

group.colors <- c("PD" = "#E69F00", "PI" = "#009E73", "OXO" ="#CC79A7", "PIPD" = "#0072B2", "PIOXO" = "#D55E00", "DMSO" = "#999999")

plot_centroids_only <- function(centroid_df,drug){
    matplot(x = (1:6),t(centroid)[1:6], type="l", lwd=4, col=group.colors[drug], ylim = c(-3,3) ,xaxt="n",ylab = "Z-score",xlab = "Time")
    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
    lines(t(centroid)[7:12], type="l", lwd=1, col="black", lty = 2)
    # Adding points
    points(x = (1:6),t(centroid)[1:6],      # Coordinates
           pch = 21,      # Symbol
           cex = 1,       # Size of the symbol
           bg = group.colors[drug],  # Background color of the symbol
           col = group.colors[drug],  # Border color of the symbol
           lwd = 3)       # Border width of the symbol
  }


pdf("clusters_test.pdf")
for (cluster in names(res_decoupler_per_cluster)) {
  
  drug <- str_extract(cluster, "[^_]+")
  cluster_no <- str_extract(cluster, "(?<=_)[^_]*")
  
  cluster_trend_name <- ls(pattern = paste0("^X.",drug,"__cl"), envir = .GlobalEnv)
  cluster_trend<-do.call("data.frame",mget(cluster_trend_name))
  colnames(cluster_trend) <- gsub(paste0(cluster_trend_name,"."),"", colnames(cluster_trend))
  
  centroid_name <- ls(pattern =  paste0("^centroids_",drug,"$"), envir = .GlobalEnv)
  centroid<-do.call("data.frame",mget(centroid_name))
  centroid <- centroid[cluster_no,]
  centroid <- data.frame("Z" = t(centroid), Time = c(0,1,2,4,8,24))
  colnames(centroid) <- c("Z","Time")
  TF <- res_decoupler_per_cluster[[cluster]]
  
  par(mfrow=c(1,2))
  par(pty = "s")

  p1 <- ggplot(centroid[1:6,], aes(x= Time, y = Z)) +
    geom_line() + theme(aspect.ratio=1) +
    ggtitle(cluster)
  
  p2 <- ggplot(data=TF,aes( y = source, x=Time, color = score, size = -log10(p_value))) +
    geom_point() +
    theme_minimal() + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" ) +
    ylab("-log10(q-value)") + theme(aspect.ratio=1)

  
  grid.arrange(p1, p2, nrow = 1)
}
dev.off()
``





# Load results --------------------------------------

load("Transcriptomics/Output/Clustering/enriched_GO_clusters_flt.RData")
load("Transcriptomics/Output/Clustering/res_decoupler_flt.RData")

