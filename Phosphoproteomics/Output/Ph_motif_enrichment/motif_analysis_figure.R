setwd("~/Documents/Git/AGS_synergies_multi-omics/Phospho_newest/output/Clustering_Ana/")
library(tidyverse)
library(RColorBrewer)
X5Z<-read.table("X5Z_motifs.txt", head=T, sep="\t")
X5Z$Selection.column<-"X5Z"
XPI<-read.table("XPI_motifs.txt", head=T, sep="\t")
XPI$Selection.column<-"PI"
XPD<-read.table("XPD_motifs.txt", head=T, sep="\t")
XPD$Selection.column<-"PD"
XPIPD<-read.table("XPIPD_motifs.txt", head=T, sep="\t")
XPIPD$Selection.column<-"PIPD"
XPI5Z<-read.table("XPI5Z_motifs.txt", head=T, sep="\t")
XPI5Z$Selection.column<-"PI5Z"

all<-rbind( XPI, XPD, XPIPD)
all<-all[all$Category.column!="NA",]
all<-droplevels(all)

motif_plot<-function(input_table){
input_table$Category.value<-gsub(" substrate motif", "", input_table$Category.value)
input_table$Category.value<-gsub(" binding motif", "", input_table$Category.value)
input_table$Category.value<-gsub(" domain", "", input_table$Category.value)
input_table$Category.value<-gsub(" kinase", "", input_table$Category.value)
input_table$Category.column<-gsub("30m", "0.5h", input_table$Category.column)

input_table$Category.value<-factor(input_table$Category.value, levels=sort(unique(input_table$Category.value)))
input_table$Category.column<-factor(input_table$Category.column, levels=rev(c("EARLY UP 0.5h",
                                                                          "EARLY UP 2h",
                                                                          "EARLY UP TRANSITORY",
                                                                          "LATE UP",
                                                                          "EARLY DOWN 0.5h",
                                                                          "EARLY DOWN 2h",
                                                                          "EARLY DOWN TRANSITORY",
                                                                          "LATE DOWN", "TRANSITORY")))

motif_barplot<-ggplot(input_table, aes(y=Category.column, x=Category.value))+
  geom_point(aes(size=Enrichment.factor, col=-log10(Benj..Hoch..FDR)))+
  #geom_point(aes(y=(-log10(Benj..Hoch..FDR))/5, x=Category.value))+
  #geom_line(aes(y=(-log10(Benj..Hoch..FDR))/5, x=Category.value, group=Selection.value))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+
  scale_color_gradientn(name="-log10 q-val",
                        colors = brewer.pal(5,"YlGnBu"),
                       limit = c(0,5),oob=scales::squish)+
  scale_size_continuous(name = "Enr.Factor",
                        limits = c(1, 20),
                        range = c(1, 7))+
  facet_grid(~Selection.column) 
print(motif_barplot)

}

motif_plot(all)       
