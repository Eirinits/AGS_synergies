library(dplyr)
library(multiClust)

setwd("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/fc_clustering")

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

PIPD<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/XPIPD_LimmaRobust.txt", head=T, sep="\t")
PI<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/XPI_LimmaRobust.txt", head=T, sep="\t")
PD<-read.table("N:/SUN-CPR-proteomics_jvo/AMV/00-COLLABS/05-Norway/02-Drug Screening Time Course/PHOSPHO/SN v16/limma_results/merged_tables/XPD_LimmaRobust.txt", head=T, sep="\t")

PIPD_sig<-dplyr::select(PIPD,contains(c("adj.P.Val", "PTM")))
PIPD_fc<-dplyr::select(PIPD,contains(c("logFC", "PTM")))

PI_sig<-dplyr::select(PI,contains(c("adj.P.Val", "PTM")))
PI_fc<-dplyr::select(PI,contains(c("logFC", "PTM")))

PD_sig<-dplyr::select(PD,contains(c("adj.P.Val", "PTM")))
PD_fc<-dplyr::select(PD,contains(c("logFC", "PTM")))

PD<-filter_sig(PD_sig, PD_fc)
colnames(PD)[2]="logFC_PD_0h"
PI<-filter_sig(PI_sig, PI_fc)
colnames(PI)[2]="logFC_PI_0h"
PIPD<-filter_sig(PIPD_sig, PIPD_fc)
colnames(PIPD)[2]="logFC_PIPD_0h"

# Calculate clusters based on the PI+5Z drug combination

hclust_analysis <- cluster_analysis(sel.exp=t(PIPD[,2:5]),
                                    cluster_type="Kmeans",
                                    distance=NULL,
                                    linkage_type=NULL, 
                                    gene_distance="pearson",
                                    num_clusters=10, data_name="PIPD", 
                                    cluster_num_selection="Fixed_Clust_Num")
hclust_analysis<-as.data.frame(hclust_analysis)

#for re-doing the plots, stat here
hclust_analysis<-read.csv("PIPD Kmeans SD_Rank Fixed_Probe_Num Fixed_Clust_Num Samples.Clusters.csv", head=T, row.names = 1)

PIPD<-merge(hclust_analysis, PIPD, by=0)
row.names(PIPD)<-PIPD$Row.names
PIPD<-PIPD[,c(2,4:7)]

# Merge the data from the PI5Z clusters with each independent treatment (keep only sites quantified in PI5Z).

synergy<-merge(PIPD, PI, by=0, all.x=T)
synergy<-merge(synergy, PD, by.x="Row.names", by.y=0, all.x=T)
synergy<-dplyr::select(synergy,contains(c("logFC", "Row.names", "hclust")))
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

ybPal <- colorRampPalette(c("gray", 'yellow', "#beebf3", "#0072B2"), alpha=T)
ybPal2 <- colorRampPalette(c("gray", 'yellow', "#80ffdd", "#009E73"), alpha=T)
ybPal3 <- colorRampPalette(c("gray", 'yellow', "#ffdf99", "#E69F00"), alpha=T)

# Plotting: for each cluster, first the PI5Z will be plot, then the 5Z and then the PI.

pdf(paste("PIPD_PI_PD_FC.pdf"), width=7, height=7)
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
