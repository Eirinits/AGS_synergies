library(dplyr)
library(multiClust)

# Function to filter by SIG in at least one treatment and measured in three time conditions

filter_sig<-function(treatment_sig, treatment_fc){
  treatment_sig$SIG <- dplyr::case_when((treatment_sig[,1]<=0.05 | treatment_sig[,2]<=0.05 | treatment_sig[,3]<=0.05  ) ~ TRUE, TRUE ~ FALSE)
  treatment_fc$vv<-apply(treatment_fc[,1:3], 1, function(x) sum(!is.na(x)))
  treatment<-merge(treatment_sig, treatment_fc, by="PTM_collapse_key")
  treatment<-treatment[treatment$SIG==TRUE,]
  treatment<-treatment[treatment$vv==3,]
  row.names(treatment)<-treatment$PTM_collapse_key
  treatment_fc<-dplyr::select(treatment,contains(c("logFC", "PTM")))
  treatment_fc$logFC_zero<-0 #add zero column for reference at time 0
  treatment_fc<-treatment_fc[,c(4,5,1,2,3)]
  treatment_fc[,2:5]<-t(scale(t(treatment_fc[,2:5])))
  return(treatment_fc)
}


# Clusters based on the drug combination, side by side with each individual drug.
# Read the tables (limma results) and separate into q-value data (e.g.: PI5Z_sig) and logFC data (e.g.: PI5Z_fc)
# Using "filter_sig" function, filter the tables to keep sites with at least one q-val <0.05 and vv>=3 in all time points.

PI5Z<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/XPI5Z_LimmaRobust.txt", head=T, sep="\t")
PI<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/XPI_LimmaRobust.txt", head=T, sep="\t")
X5Z<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/X5Z_LimmaRobust.txt", head=T, sep="\t")

PI5Z_sig<-dplyr::select(PI5Z,contains(c("adj.P.Val", "PTM")))
PI5Z_fc<-dplyr::select(PI5Z,contains(c("logFC", "PTM")))

PI_sig<-dplyr::select(PI,contains(c("adj.P.Val", "PTM")))
PI_fc<-dplyr::select(PI,contains(c("logFC", "PTM")))

X5Z_sig<-dplyr::select(X5Z,contains(c("adj.P.Val", "PTM")))
X5Z_fc<-dplyr::select(X5Z,contains(c("logFC", "PTM")))

X5Z<-filter_sig(X5Z_sig, X5Z_fc)
colnames(X5Z)[2]="logFC_5Z_0h"
PI<-filter_sig(PI_sig, PI_fc)
colnames(PI)[2]="logFC_PI_0h"
PI5Z<-filter_sig(PI5Z_sig, PI5Z_fc)
colnames(PI5Z)[2]="logFC_PI5Z_0h"

# Calculate clusters based on the PI+5Z drug combination

hclust_analysis <- cluster_analysis(sel.exp=t(PI5Z[,2:5]),
                                    cluster_type="Kmeans",
                                    distance=NULL,
                                    linkage_type=NULL, 
                                    gene_distance="pearson",
                                    num_clusters=10, data_name="PI5Z", 
                                    cluster_num_selection="Fixed_Clust_Num")
hclust_analysis<-as.data.frame(hclust_analysis)

#for re-doing the plots, stat here
hclust_analysis<-read.csv("PI5Z Kmeans SD_Rank Fixed_Probe_Num Fixed_Clust_Num Samples.Clusters.csv", head=T, row.names = 1)
PI5Z<-merge(hclust_analysis, PI5Z, by=0)
row.names(PI5Z)<-PI5Z$Row.names
PI5Z<-PI5Z[,c(2,4:7)]

# Merge the data from the PI5Z clusters with each independent treatment (keep only sites quantified in PI5Z).

synergy<-merge(PI5Z, PI, by=0, all.x=T)
synergy<-merge(synergy, X5Z, by.x="Row.names", by.y=0, all.x=T)
synergy<-dplyr::select(synergy,contains(c("logFC", "Row.names", "hclust_analysis")))
row.names(synergy)<-synergy$Row.names
synergy<-synergy[,c(14,1:12)]

# Calculate centroids of the clusters

all_centroids<-data.frame(matrix(0, nrow=10, ncol=12)) #rows=clusters col=samples
colnames(all_centroids)<-colnames(synergy)[2:13]

for(i in 1:12){
  all_centroids[,i]<-tapply(as.matrix(synergy[,i+1]), as.factor(synergy$hclust_analysis), function(x) median(x, na.rm=T))
}



# Calculate correlation of each site to the centroid of its cluster.
# Order them by correlation within in cluster. 
# That way they will be plot in order: from back (less correlation) to front (higher correlation).

for(i in 1:nrow(synergy)){
  for(j in 1:10){
    if(synergy[i,"hclust_analysis"]==j){
      synergy[i,14]<-cor(as.numeric(synergy[i,2:13]),as.numeric(all_centroids[j,1:12]), method="pearson")}
  }
}  
colnames(synergy)[14]<-"CORR"
synergy<-synergy[order(synergy$hclust_analysis, synergy$CORR),]

# Function to calculate the color gradient, which will be mapped to the CORR column.

ybPal <- colorRampPalette(c("gray", 'yellow', "#ffc599", "#D55E00"), alpha=T)
ybPal2 <- colorRampPalette(c("gray", 'yellow', "#80ffdd", "#009E73"), alpha=T)
ybPal3 <- colorRampPalette(c("gray", 'yellow', "#eac8da", "#CC79A7"), alpha=T)

# Plotting: for each cluster, first the PI5Z will be plot, then the 5Z and then the PI.

pdf(paste("PI5Z_PI_X5Z_FC_new_colors.pdf"), width=10, height=7)
par(mfrow=c(5,6),mai = c(0.5, 0.3, 0.1, 0.1))

for(i in 1:10){
  centroid<-all_centroids[i,1:12]
  df<-synergy[which(synergy$hclust_analysis==i),2:13]
  colors<-synergy[which(synergy$hclust_analysis==i), "CORR"]
  colors_gr <- ybPal(12)[as.numeric(cut(colors,breaks = 10))]
  colors_gr_2 <- ybPal2(12)[as.numeric(cut(colors,breaks = 10))]
  colors_gr_3 <- ybPal3(12)[as.numeric(cut(colors,breaks = 10))]
  matplot(t(df[,1:4]), type="l", lty=1, col=alpha(colors_gr,0.3), ylim=c(-2,2),ylab="")
  lines(t(centroid[1:4]), type="l", lwd=4, col="white")
  lines(t(centroid[1:4]), type="l", lwd=2, col="black")
  matplot(t(df[,5:8]), type="l", lty=1, col=alpha(colors_gr_2,0.3), ylim=c(-2,2),ylab="")
  lines(t(centroid[5:8]), type="l", lwd=4, col="white")
  lines(t(centroid[5:8]), type="l", lwd=2, col="black")
  matplot(t(df[,9:12]), type="l", lty=1, col=alpha(colors_gr_3,0.3), ylim=c(-2,2),ylab="")
  lines(t(centroid[9:12]), type="l", lwd=4, col="white")
  lines(t(centroid[9:12]), type="l", lwd=2, col="black")
}

dev.off()


