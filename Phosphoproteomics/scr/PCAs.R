library(ggfortify)
library(stringr)
setwd("D:/01-PROJECTs/00-CPR/04-SYNERGIES")

input=read.table("Phosphoproteomics/Data/PTM_collapsed_log2_norm_imp_18517.txt", head=T, sep="\t", quote="")
input<-t(input)
colnames(input)<-input[1,]
input<-input[-1,]
input<-as.data.frame(input)
input$condition<-row.names(input)
input$condition<-gsub("0_5", "0.5", input$condition)
input$condition<-gsub("XCTRL", "XCONTROL", input$condition)
input$condition<-gsub("X", "", input$condition)
input[,c("drug", "time", "replicate")]<-str_split_fixed(input$condition, "_",3)
input$drug <- gsub("CONTROL","Untreated", input$drug)
input$drug <- gsub("PD","MEKi", input$drug)
input$drug <- gsub("PI","PI3Ki", input$drug)
input$drug <- gsub("5Z","TAKi", input$drug)
input$drug <- gsub("PI3KiMEKi","PI3Ki MEKi", input$drug)
input$drug <- gsub("PI3KiTAKi","PI3Ki TAKi", input$drug)
input$drug<-factor(input$drug, levels=c("Untreated", "DMSO", "MEKi", "PI3Ki", "TAKi", "PI3Ki MEKi", "PI3Ki TAKi"))
input$time<-factor(input$time, levels=c("0h","0.5h","2h","8h"))

input[,1:18517]<-apply(input[,1:18517], 2, as.numeric)
pca<-prcomp(input[,1:18517])

group.colors <- c("MEKi" = "#E69F00", "PI3Ki" = "#009E73", "TAKi" ="#CC79A7", "PI3Ki MEKi" = "#0072B2", "PI3Ki TAKi" = "#D55E00", "DMSO" = "#999999", "Untreated" = "black")

autoplot(pca, data=input, colour="drug", shape="time")+facet_grid(~drug)+theme_bw() +  scale_color_manual(values=group.colors) + coord_flip()

