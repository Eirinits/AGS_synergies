# Load libraries & data ---------------------------------------------------

library(tidyverse)
library(Biobase)
library(cluster)
library(stringi)
library(PCAtools)
library(plotly)
library(stats)

gene_counts <- read.delim("Transcriptomics/Data/Raw_data/gene_counts_renamed.tsv")
norm_gene_counts <- read.delim("Transcriptomics/Data/gene_counts_flt_norm_wCTRL.csv", sep = ",")
samples <- read_delim("Transcriptomics/Data/sample_info_to_use.tsv", delim = "\t", escape_double = FALSE)

samples <- samples %>% 
  mutate(treatment = factor(treatment)) %>% 
  mutate(treatment = relevel(treatment, ref="DMSO")) %>% 
  mutate(Timepoint = factor(paste0(Timepoint,"h"), levels=c("0h", "1h", "2h", "4h", "8h", "24h"))) %>% 
  filter(ID %in% colnames(norm_gene_counts)) %>% 
  column_to_rownames("ID") 

# PCA ---------------------------------------------------------------------


p <- pca(norm_gene_counts, metadata = samples, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

pairsplot(p,
          components = getComponents(p, c(1:3)),
          colby = "treatment",shape = "Timepoint")

eigencorplot(p,
             metavars = c('Timepoint','treatment','Condition'))



pca<-prcomp(t(log10(norm_gene_counts)))

components <- pca[["x"]]

components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

components <- components %>% 
  rownames_to_column("Condition") %>% 
  separate(Condition, into = c("Treatment","Timepoint"), remove = F) %>% 
  mutate(Timepoint = as.numeric(Timepoint))

grp.colors <- c("PD" = "#E69F00", "PI" = "#009E73", "OXO" ="#CC79A7", "PIPD" = "#0072B2", "PIOXO" = "#D55E00", "DMSO" = "#999999", "untreated" = "black")

p1 <- plot_ly(components, x = ~PC1, y = ~PC2,
               color = ~components$Treatment, colors = grp.colors,
               size = ~components$Time) %>% 
  layout(showlegend = F) %>% 
  layout(xaxis = list(title = 'PC1 - 37.4%'),
         yaxis = list(title = 'PC2 - 16.7%')) %>% 
  layout( xaxis = list(titlefont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22)))
p1
save_image(p1, file = "Figures/tr_PC1_PC2.svg")

p2 <- plot_ly(components, x = ~PC3, y = ~PC2,
        color = ~components$Treatment, colors = grp.colors,
        size = ~components$Time) %>% 
        layout(yaxis = list(title = 'PC2 - 16.7%'),
               xaxis = list(title = 'PC3 - 8.95%'))%>% 
  layout( yaxis = list(titlefont = list(size = 22)),
          xaxis = list(titlefont = list(size = 22))) %>% 
  layout(showlegend = F) 
  
p2
save_image(p2, file = "Figures/tr_PC2_PC3.svg")


