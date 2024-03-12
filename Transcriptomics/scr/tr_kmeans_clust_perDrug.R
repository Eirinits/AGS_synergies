library(dplyr)
library(multiClust)

drugs=c("OXO_", "PIPD_", "PIOXO_", "PI_")

DEGs_perContrast <- read.csv("Transcriptomics/Output/DEGs_perContrast_all.csv")
counts<-read.table("Transcriptomics/Data/gene_counts_flt_norm_wCTRL.csv", row.names = 1, sep = ",")

col_palettes <- list(PIOXO_ = colorRampPalette(c("gray", 'yellow', "#ffc599", "#D55E00"), alpha=T),
                     OXO_ =  colorRampPalette(c("gray", 'yellow', "#80ffdd", "#CC79A7"), alpha=T),
                     PI_ = colorRampPalette(c("gray", 'yellow', "#eac8da", "#009E73"), alpha=T),
                     PIPD_ = colorRampPalette(c("gray", 'yellow', "#beebf3", "#0072B2"), alpha=T),
                     PD_ = colorRampPalette(c("gray", 'yellow', "#ffdf99", "#E69F00"), alpha=T))

group.colors <- c("PD" = "#E69F00", "PI" = "#009E73", "OXO" ="#CC79A7", "PIPD" = "#0072B2", "PIOXO" = "#D55E00", "DMSO" = "#999999")

cl_num <- 20

for(j in 1:length(drugs)){
  
  drug<-drugs[j]
  treatment<-dplyr::select(counts,starts_with(c(drug, "DMSO", "untreated")))
  treatment_sig <- subset(DEGs_perContrast, grepl(paste0("Diff_",drug),Contrast) & adj.P.Val < 0.05 &  abs(logFC) > 0.5)
  
  treatment <- treatment[which(rownames(treatment) %in% treatment_sig$Gene_name),]
  treatment_only<-dplyr::select(treatment,contains(c(drug, "untreated")))
  
  groups<-c(rep("1h",3),rep("2h",3),rep("24h",3),rep("4h",3),rep("8h",3), rep("C", 5))
  treatment_median<-t(apply(treatment_only, 1, function(x) tapply(x, groups, median)))
  
  DMSO<-dplyr::select(treatment,contains(c("DMSO", "untreated")))
  groups_dmso<-c(rep("DMSO_1h",3),rep("DMSO_2h",3),rep("DMSO_24h",3),rep("DMSO_4h",3),rep("DMSO_8h",3),rep("C",5))
  DMSO_median<-t(apply(DMSO, 1, function(x) tapply(x, groups_dmso, median)))
  
  treatment_median<-treatment_median[,c("C", "1h", "2h", "4h","8h","24h")]
  treatment_median<-t(scale(t(treatment_median)))
  treatment_median<-as.data.frame(treatment_median)
  
  DMSO_median<-DMSO_median[,c("C", "DMSO_1h", "DMSO_2h","DMSO_4h", "DMSO_8h","DMSO_24h")]
  DMSO_median<-t(scale(t(DMSO_median)))
  DMSO_median<-as.data.frame(DMSO_median)
  
  hclust_analysis <- cluster_analysis(sel.exp=t(treatment_median),
                                      cluster_type="Kmeans",
                                      distance=NULL,
                                      linkage_type=NULL, 
                                      gene_distance="pearson",
                                      num_clusters=cl_num, data_name=paste(drug, "clusters",cl_num), seed = 1234,
                                      cluster_num_selection="Fixed_Clust_Num")
  hclust_analysis<-as.data.frame(hclust_analysis)
  # clusters_transl <- hclust_analysis %>%
  #   rownames_to_column("ENS") %>% 
  #   left_join(features)

  treatment_cluster<-merge(hclust_analysis, treatment_median, by=0)
  row.names(treatment_cluster)<-treatment_cluster$Row.names
  treatment_cluster<-within(treatment_cluster, rm(Row.names))
  treatment_cluster<-merge(treatment_cluster, DMSO_median, by=0)
  
  treatment_cluster[,3:8]<-treatment_cluster[,3:8]-treatment_cluster[,3]
  treatment_cluster[,9:14]<-treatment_cluster[,9:14]-treatment_cluster[,9]
  
  
  all_centroids<-data.frame(matrix(0, nrow=cl_num, ncol=12)) #rows=clusters col=samples
  colnames(all_centroids)<-colnames(treatment_cluster)[3:14]
  
  for(i in 1:12){
    all_centroids[,i]<-tapply(as.matrix(treatment_cluster[,i+2]), as.factor(treatment_cluster$hclust_analysis), function(x) median(x, na.rm=T))
  }
  
  for(i in 1:nrow(treatment_cluster)){
    for(j in 1:cl_num){
      if(treatment_cluster[i,"hclust_analysis"]==j){
        treatment_cluster[i,15]<-cor(as.numeric(treatment_cluster[i,3:14]),as.numeric(all_centroids[j,1:12]), method="pearson")}
    }
  }  
  colnames(treatment_cluster)[15]<-"CORR"
  treatment_cluster<-treatment_cluster[order(treatment_cluster$hclust_analysis, treatment_cluster$CORR),]

  write.table(treatment_cluster, file = paste0("Transcriptomics/Output/Clustering/",drug,"_clusters.csv"), sep = ",", row.names = F)
  pdf(paste0("Transcriptomics/Output/Clustering/",drug, "clusters.pdf"), width=10, height=7)
  par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
  par(pty = "s")

    for(i in 1:cl_num){
      
      centroid<-as.data.frame(all_centroids[i,1:12])
      df<-treatment_cluster[which(treatment_cluster$hclust_analysis==i),3:14]
      colors<-treatment_cluster[which(treatment_cluster$hclust_analysis==i), "CORR"]
      palette_gr <- col_palettes[[drug]]
      colors_gr <- palette_gr(12)[as.numeric(cut(colors,breaks = 10))]
      matplot(t(df[,1:6]), type="l",lwd = 2, lty=1, col=alpha(colors_gr,0.3), ylim=c(-3,3),ylab="", xaxt='n')
      axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
      lines(t(centroid)[1:6], type="l", lwd=4, col="white")
      lines(t(centroid)[1:6], type="l", lwd=2, col="black")
  #    matplot(t(df[,7:12]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-4,4),ylab="", xaxt='n')
  #    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
      lines(t(centroid)[7:12], type="l", lwd=1, col="red", lty = 2)
    }
  dev.off()
}

drug<-"PD_"
treatment<-dplyr::select(counts,starts_with(c(drug, "DMSO", "untreated")))
treatment_sig <- subset(DEGs_perContrast, grepl(paste0("Diff_",drug),Contrast) & adj.P.Val < 0.05 &  abs(logFC) > 0.5)

treatment <- treatment[which(rownames(treatment) %in% treatment_sig$Gene_name),]
treatment_only<-dplyr::select(treatment,contains(c(drug, "untreated")))

groups<-c(rep("1h",2),rep("2h",3),rep("24h",3),rep("4h",3),rep("8h",3), rep("C", 5))
treatment_median<-t(apply(treatment_only, 1, function(x) tapply(x, groups, median)))

DMSO<-dplyr::select(treatment,contains(c("DMSO", "untreated")))
groups_dmso<-c(rep("DMSO_1h",3),rep("DMSO_2h",3),rep("DMSO_24h",3),rep("DMSO_4h",3),rep("DMSO_8h",3),rep("C",5))
DMSO_median<-t(apply(DMSO, 1, function(x) tapply(x, groups_dmso, median)))

treatment_median<-treatment_median[,c("C", "1h", "2h", "4h","8h","24h")]
treatment_median<-t(scale(t(treatment_median)))
treatment_median<-as.data.frame(treatment_median)

DMSO_median<-DMSO_median[,c("C", "DMSO_1h", "DMSO_2h","DMSO_4h", "DMSO_8h","DMSO_24h")]
DMSO_median<-t(scale(t(DMSO_median)))
DMSO_median<-as.data.frame(DMSO_median)

hclust_analysis <- cluster_analysis(sel.exp=t(treatment_median),
                                    cluster_type="Kmeans",
                                    distance=NULL,
                                    linkage_type=NULL, 
                                    gene_distance="pearson",
                                    num_clusters=cl_num, data_name=paste(drug, "clusters"), seed = 1234,
                                    cluster_num_selection="Fixed_Clust_Num")
hclust_analysis<-as.data.frame(hclust_analysis)

treatment_cluster<-merge(hclust_analysis, treatment_median, by=0)
row.names(treatment_cluster)<-treatment_cluster$Row.names
treatment_cluster<-within(treatment_cluster, rm(Row.names))
treatment_cluster<-merge(treatment_cluster, DMSO_median, by=0)

treatment_cluster[,3:8]<-treatment_cluster[,3:8]-treatment_cluster[,3]
treatment_cluster[,9:14]<-treatment_cluster[,9:14]-treatment_cluster[,9]

all_centroids<-data.frame(matrix(0, nrow=cl_num, ncol=12)) #rows=clusters col=samples
colnames(all_centroids)<-colnames(treatment_cluster)[3:14]

for(i in 1:12){
  all_centroids[,i]<-tapply(as.matrix(treatment_cluster[,i+2]), as.factor(treatment_cluster$hclust_analysis), function(x) median(x, na.rm=T))
}

for(i in 1:nrow(treatment_cluster)){
  for(j in 1:cl_num){
    if(treatment_cluster[i,"hclust_analysis"]==j){
      treatment_cluster[i,15]<-cor(as.numeric(treatment_cluster[i,3:14]),as.numeric(all_centroids[j,1:12]), method="pearson")}
  }
}  
colnames(treatment_cluster)[15]<-"CORR"
treatment_cluster<-treatment_cluster[order(treatment_cluster$hclust_analysis, treatment_cluster$CORR),]
write.table(treatment_cluster, file = paste0("Transcriptomics/Output/Clustering/",drug,"_clusters.csv"), sep = ",", row.names = F)

#write.table(centroid, paste(drug, "centroid_clusters.txt"), sep="\t")
pdf(paste("Transcriptomics/Output/Clustering/",drug, "clusters.pdf", sep=""), width=10, height=7)
par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
par(pty = "s")

for(i in 1:cl_num){
  
  centroid<-as.data.frame(all_centroids[i,1:12])
  df<-treatment_cluster[which(treatment_cluster$hclust_analysis==i),3:14]
  colors<-treatment_cluster[which(treatment_cluster$hclust_analysis==i), "CORR"]
  palette_gr <- col_palettes[[drug]]
  colors_gr <- palette_gr(12)[as.numeric(cut(colors,breaks = 10))]
  matplot(t(df[,1:6]), type="l",lwd = 2, lty=1, col=alpha(colors_gr,0.3), ylim=c(-4,4),ylab="", xaxt='n')
  axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
  lines(t(centroid)[1:6], type="l", lwd=4, col="white")
  lines(t(centroid)[1:6], type="l", lwd=2, col="black")
  #    matplot(t(df[,7:12]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-4,4),ylab="", xaxt='n')
  #    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
  lines(t(centroid)[7:12], type="l", lwd=1, col="red", lty = 2)
}
dev.off()


# Visualize clusters without rerunning  -----------------------------------

temp = list.files(path ="Transcriptomics/Output/Clustering/", pattern="*_clusters.csv", full.names = T)
cluster_files = lapply(temp, read.csv)
cluster_files = setNames(cluster_files, temp)
names(cluster_files) <- gsub("Transcriptomics/Output/Clustering//|__clusters.csv","", names(cluster_files))

plot_clusters_in_separate_files <- function(drug)
  {
  col_palettes <- list(PIOXO = colorRampPalette(c("gray", '#e6b793', "#D55E00"), alpha=T),
                       OXO =  colorRampPalette(c("gray","#d1b6c5", "#CC79A7"), alpha=T),
                       PI = colorRampPalette(c("gray", "#e0fff7", "#009E73"), alpha=T),
                       PIPD = colorRampPalette(c("gray", '#95bacf', "#0072B2"), alpha=T),
                       PD = colorRampPalette(c("gray", '#e3c174', "#E69F00"), alpha=T))
  
  treatment_cluster <- cluster_files[[drug]]
  all_centroids<-data.frame(matrix(0, nrow=20, ncol=12)) #rows=clusters col=samples
  colnames(all_centroids)<-colnames(treatment_cluster)[3:14]
  
  for(i in 1:12){
    all_centroids[,i]<-tapply(as.matrix(treatment_cluster[,i+2]), as.factor(treatment_cluster$hclust_analysis), function(x) median(x, na.rm=T))
  }
  
  for(i in 1:nrow(treatment_cluster)){
    for(j in 1:20){
      if(treatment_cluster[i,"hclust_analysis"]==j){
        treatment_cluster[i,15]<-cor(as.numeric(treatment_cluster[i,3:14]),as.numeric(all_centroids[j,1:12]), method="pearson")}
    }
  }  
  colnames(treatment_cluster)[15]<-"CORR"
  treatment_cluster<-treatment_cluster[order(treatment_cluster$hclust_analysis, treatment_cluster$CORR,decreasing = F),]
  
#  pdf(paste0("Transcriptomics/Output/Clustering/",drug, "clusters_PI_test.pdf"), width=10, height=7)
#  par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
#  par(pty = "s")
  
  for(i in 1:20){
    pdf(paste0("Transcriptomics/Output/Clustering/Individual_cluster_figures/",drug,"_",i, ".pdf"), width=10, height=7)
    par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
    par(pty = "s")
    centroid<-as.data.frame(all_centroids[i,1:12])
    df<-treatment_cluster[which(treatment_cluster$hclust_analysis==i),3:14]
    colors<-treatment_cluster[which(treatment_cluster$hclust_analysis==i), "CORR"]
    palette_gr <- col_palettes[[drug]]
    colors_gr <- palette_gr(12)[as.numeric(cut(colors,breaks = 10))]
    matplot(t(df[,1:6]), type="l",lwd = 2, lty=1, col=alpha(colors_gr,0.3), ylim=c(-3,3),ylab="", xaxt='n')
    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
    lines(t(centroid)[1:6], type="l", lwd=4, col="white")
    lines(t(centroid)[1:6], type="l", lwd=2, col="black")
    #    matplot(t(df[,7:12]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-4,4),ylab="", xaxt='n')
    #    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
    lines(t(centroid)[7:12], type="l", lwd=1, col="red", lty = 2)
    
    dev.off()
  }
}
plot_clusters_in_separate_files("PD")
plot_clusters_in_separate_files("OXO")
plot_clusters_in_separate_files("PIPD")
plot_clusters_in_separate_files("PIOXO")
plot_clusters_in_separate_files("PI")

col_palettes <- list(PIOXO = colorRampPalette(c("gray", 'yellow', "#ffc599", "#D55E00"), alpha=T),
                     OXO =  colorRampPalette(c("gray","#d1b6c5", "#CC79A7"), alpha=T),
                     PI = colorRampPalette(c("gray", "#e0fff7", "#009E73"), alpha=T),
                     PIPD = colorRampPalette(c("gray", 'yellow', "#beebf3", "#0072B2"), alpha=T),
                     PD = colorRampPalette(c("gray", 'yellow', "#ffdf99", "#E69F00"), alpha=T))

treatment_cluster <- cluster_files[["PI"]]
all_centroids<-data.frame(matrix(0, nrow=20, ncol=12)) #rows=clusters col=samples
colnames(all_centroids)<-colnames(treatment_cluster)[3:14]

for(i in 1:12){
  all_centroids[,i]<-tapply(as.matrix(treatment_cluster[,i+2]), as.factor(treatment_cluster$hclust_analysis), function(x) median(x, na.rm=T))
}

for(i in 1:nrow(treatment_cluster)){
  for(j in 1:20){
    if(treatment_cluster[i,"hclust_analysis"]==j){
      treatment_cluster[i,15]<-cor(as.numeric(treatment_cluster[i,3:14]),as.numeric(all_centroids[j,1:12]), method="pearson")}
  }
}  
colnames(treatment_cluster)[15]<-"CORR"
treatment_cluster<-treatment_cluster[order(treatment_cluster$hclust_analysis, treatment_cluster$CORR,decreasing = F),]

drug = "PI"
pdf(paste0("Transcriptomics/Output/Clustering/",drug, "clusters_PI_test.pdf"), width=10, height=7)
par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
par(pty = "s")

for(i in 1:20){
  pdf(paste0("Transcriptomics/Output/Clustering/Individual_cluster_figures/",drug,"_",i, ".pdf"), width=10, height=7)
  par(mfrow=c(3,5),mai = c(0.5, 0.3, 0.5, 0.3))
  par(pty = "s")
  centroid<-as.data.frame(all_centroids[i,1:12])
  df<-treatment_cluster[which(treatment_cluster$hclust_analysis==i),3:14]
  colors<-treatment_cluster[which(treatment_cluster$hclust_analysis==i), "CORR"]
  palette_gr <- col_palettes[[drug]]
  colors_gr <- palette_gr(12)[as.numeric(cut(colors,breaks = 10))]
  matplot(t(df[,1:6]), type="l",lwd = 2, lty=1, col=alpha(colors_gr,0.3), ylim=c(-3,3),ylab="", xaxt='n')
  axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
  lines(t(centroid)[1:6], type="l", lwd=4, col="white")
  lines(t(centroid)[1:6], type="l", lwd=2, col="black")
  #    matplot(t(df[,7:12]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-4,4),ylab="", xaxt='n')
  #    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
  lines(t(centroid)[7:12], type="l", lwd=1, col="red", lty = 2)
  
  dev.off()
}


plot_centroids_only <- function(centroid_df,drug){
  pdf(paste0("Transcriptomics/Output/Clustering/",drug,"_centroid_plots.pdf"))
  
  for(i in 1:20){
    centroid<-as.data.frame(centroid_df[i,1:12])
  #  palette_gr <- col_palettes[[drug]]
  #  colors_gr <- palette_gr(12)[as.numeric(cut(colors,breaks = 10))]
    matplot(x = (1:6),t(centroid)[1:6], type="l", lwd=4, col=group.colors[drug], ylim = c(-3,3) ,xaxt="n",ylab = "Z-score",xlab = "Time")
    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
    #lines(t(centroid)[1:6], type="l", lwd=2, col="black")
    #    matplot(t(df[,7:12]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-4,4),ylab="", xaxt='n')
    #    axis(side=1,at=1:6,labels=c("0","1","2","4","8","24"))
    lines(t(centroid)[7:12], type="l", lwd=1, col="black", lty = 2)
    
    # Adding points
    points(x = (1:6),t(centroid)[1:6],      # Coordinates
           pch = 21,      # Symbol
           cex = 1,       # Size of the symbol
           bg = group.colors[drug],  # Background color of the symbol
           col = group.colors[drug],  # Border color of the symbol
           lwd = 3)       # Border width of the symbol
  }
  dev.off()
}

plot_centroids_only(centroids_OXO,"OXO")
plot_centroids_only(centroids_PIOXO,"PIOXO")
plot_centroids_only(centroids_PI,"PI")
plot_centroids_only(centroids_PD,"PD")
plot_centroids_only(centroids_PIPD,"PIPD")

