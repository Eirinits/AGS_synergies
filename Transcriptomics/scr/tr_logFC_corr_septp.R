library(tidyverse)
library(gridExtra)
library(grid)

DEGs_perContrast <- read.csv("Transcriptomics/Output/DEGs_perContrast_all.csv")

DEGs_perContrast_flt <- DEGs_perContrast %>% 
  mutate(Contrast = gsub("Diff_","",Contrast)) %>% 
  separate(Contrast, into = c("Treatment","Timepoint"), remove = F, sep = "_")  %>% 
  dplyr::filter(adj.P.Val < 0.05)

PIPD_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "PIPD") %>% 
  unite(Gene_tp, Gene_name, Timepoint )

PI_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "PI") %>% 
  unite(Gene_tp, Gene_name, Timepoint ) %>% 
  dplyr::filter(Gene_tp %in% PIPD_DEGS$Gene_tp) 

PD_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "PD")%>% 
  unite(Gene_tp, Gene_name, Timepoint ) %>% 
  dplyr::filter(Gene_tp %in% PIPD_DEGS$Gene_tp) 

PI_and_PIPD <- left_join(PI_DEGS, PIPD_DEGS, by = "Gene_tp",suffix = c(".PI", ".PIPD")) %>% 
  separate(Gene_tp, c("Gene","Timepoint"), sep = "_") %>% 
  mutate(Timepoint = factor(Timepoint, levels = c("1h","2h","4h","8h","24h")))

PD_and_PIPD <- left_join(PD_DEGS, PIPD_DEGS, by = "Gene_tp",suffix = c(".PD", ".PIPD"))%>% 
  separate(Gene_tp, c("Gene","Timepoint"), sep = "_") %>% 
  mutate(Timepoint = factor(Timepoint, levels = c("1h","2h","4h","8h","24h")))

plot1 <- ggplot(PI_and_PIPD, aes(x = logFC.PI, y = logFC.PIPD)) + 
  geom_hex(bins= 25) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500)), ) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1) +
  theme_minimal() +
  ylim(-10,10)+
  xlim(-10,10)+
 # facet_grid(~Timepoint, scales = "free") +
  labs( x ="LogFC PI3Ki", y = "LogFC MEKi+PI3Ki")

plot1 
plot2 <-ggplot(PD_and_PIPD, aes(x = logFC.PD, y = logFC.PIPD)) + 
  geom_hex(bins= 25) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500)), ) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1) +
  theme_minimal() +
  ylim(-10,10)+
  xlim(-10,10)+
 # facet_grid(~Timepoint, scales = "free") +
  labs( x ="LogFC MEKi", y = "LogFC MEKi+PI3Ki")

plot2

PIOXO_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "PIOXO") %>% 
  unite(Gene_tp, Gene_name, Timepoint )

PI_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "PI") %>% 
  unite(Gene_tp, Gene_name, Timepoint ) %>% 
  dplyr::filter(Gene_tp %in% PIOXO_DEGS$Gene_tp) 

OXO_DEGS <- DEGs_perContrast_flt %>% 
  dplyr::filter(Treatment == "OXO")%>% 
  unite(Gene_tp, Gene_name, Timepoint ) %>% 
  dplyr::filter(Gene_tp %in% PIOXO_DEGS$Gene_tp) 

PI_and_PIOXO <- left_join(PI_DEGS, PIOXO_DEGS, by = "Gene_tp",suffix = c(".PI", ".PIOXO"))%>% 
  separate(Gene_tp, c("Gene","Timepoint"), sep = "_") %>% 
  mutate(Timepoint = factor(Timepoint, levels = c("1h","2h","4h","8h","24h")))

OXO_and_PIOXO <- left_join(OXO_DEGS, PIOXO_DEGS, by = "Gene_tp",suffix = c(".OXO", ".PIOXO")) %>% 
  separate(Gene_tp, c("Gene","Timepoint"), sep = "_") %>% 
  mutate(Timepoint = factor(Timepoint, levels = c("1h","2h","4h","8h","24h")))

plot3 <-ggplot(PI_and_PIOXO, aes(x = logFC.PI, y = logFC.PIOXO)) + 
  geom_hex(bins= 25) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500,2000)), ) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1) +
  theme_minimal() +
  ylim(-10,10)+
  xlim(-10,10)+
#  facet_grid(~Timepoint, scales = "free") +
  labs( x ="LogFC PI3Ki", y = "LogFC TAKi+PI3Ki")

plot3
plot4 <-ggplot(OXO_and_PIOXO, aes(x = logFC.OXO, y = logFC.PIOXO)) + 
  geom_hex(bins= 25) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500,2000)), ) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1) +
  ylim(-10,10)+
  xlim(-10,10)+
  theme_minimal() +
  #facet_grid(~Timepoint, scales = "free") +
  labs( x ="LogFC TAKi", y = "LogFC TAKi+PI3Ki")

plot4
grid.arrange(plot1, plot2, plot3, plot4, ncol=4,
             top = textGrob("Transcriptomics",gp=gpar(fontsize=15)))

g <- arrangeGrob(plot1, plot2, plot3, plot4, ncol=4, nrow =2,
                 top = textGrob("Transcriptomics",gp=gpar(fontsize=20)))
ggsave("log2_correlation_transcriptomics.png",g,
       path = paste0(wd$figures),
       width = 15,
       height = 15,
       units = "cm",
       dpi = 300)

cor(PI_and_PIPD$logFC.PI, PI_and_PIPD$logFC.PIPD, method=c("pearson"))
cor(PD_and_PIPD$logFC.PD, PD_and_PIPD$logFC.PIPD, method=c("pearson"))
cor(PI_and_PIOXO$logFC.PI, PI_and_PIOXO$logFC.PIOXO, method=c("pearson"))
cor(OXO_and_PIOXO$logFC.OXO, OXO_and_PIOXO$logFC.PIOXO, method=c("pearson"))

cor(PI_and_PIPD$logFC.PI[PI_and_PIPD$Timepoint == "1h"], PI_and_PIPD$logFC.PIPD[PI_and_PIPD$Timepoint == "1h"], method=c("pearson"))
cor(PI_and_PIPD$logFC.PI[PI_and_PIPD$Timepoint == "2h"], PI_and_PIPD$logFC.PIPD[PI_and_PIPD$Timepoint == "2h"], method=c("pearson"))
cor(PI_and_PIPD$logFC.PI[PI_and_PIPD$Timepoint == "4h"], PI_and_PIPD$logFC.PIPD[PI_and_PIPD$Timepoint == "4h"], method=c("pearson"))
cor(PI_and_PIPD$logFC.PI[PI_and_PIPD$Timepoint == "8h"], PI_and_PIPD$logFC.PIPD[PI_and_PIPD$Timepoint == "8h"], method=c("pearson"))
cor(PI_and_PIPD$logFC.PI[PI_and_PIPD$Timepoint == "24h"], PI_and_PIPD$logFC.PIPD[PI_and_PIPD$Timepoint == "24h"], method=c("pearson"))

cor(PD_and_PIPD$logFC.PD[PD_and_PIPD$Timepoint == "1h"], PD_and_PIPD$logFC.PIPD[PD_and_PIPD$Timepoint == "1h"], method=c("pearson"))
cor(PD_and_PIPD$logFC.PD[PD_and_PIPD$Timepoint == "2h"], PD_and_PIPD$logFC.PIPD[PD_and_PIPD$Timepoint == "2h"], method=c("pearson"))
cor(PD_and_PIPD$logFC.PD[PD_and_PIPD$Timepoint == "4h"], PD_and_PIPD$logFC.PIPD[PD_and_PIPD$Timepoint == "4h"], method=c("pearson"))
cor(PD_and_PIPD$logFC.PD[PD_and_PIPD$Timepoint == "8h"], PD_and_PIPD$logFC.PIPD[PD_and_PIPD$Timepoint == "8h"], method=c("pearson"))
cor(PD_and_PIPD$logFC.PD[PD_and_PIPD$Timepoint == "24h"], PD_and_PIPD$logFC.PIPD[PD_and_PIPD$Timepoint == "24h"], method=c("pearson"))

