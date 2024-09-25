library(tidyverse)
require("corrplot")

# LogFC correlation -------------------------------------------------------

DPPs_vs_DMSO_long <- read.csv("Phosphoproteomics/Data/DPPs_long.csv")

DPPs_vs_DMSO_long <- DPPs_vs_DMSO_long %>% 
  separate(Condition, c("Treatment","Timepoint"), sep = "_", remove = F) %>%
#  filter(p_val < 0.05) %>% 
  rename(PTM_id = "PTM_collapse_key")

PIPD_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "PIPD") %>% 
  unite(PTM_tp, PTM_id, Timepoint, remove = F )

PI_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "PI") %>% 
  unite(PTM_tp, PTM_id, Timepoint, remove = F  ) %>% 
  dplyr::filter(PTM_tp %in% PIPD_DPPS$PTM_tp) 

PD_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "PD")%>% 
  unite(PTM_tp, PTM_id, Timepoint , remove = F ) %>% 
  dplyr::filter(PTM_tp %in% PIPD_DPPS$PTM_tp) 

PD_and_PIPD <- left_join(PD_DPPS, PIPD_DPPS, by = c("PTM_tp","Timepoint"),  suffix = c(".PD", ".PIPD"))

ggplot(PD_and_PIPD, aes(x = LogFC.PD, y = LogFC.PIPD)) + 
  geom_hex(bins= 70) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500))) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1) +
  theme_minimal() +
  facet_grid(~Timepoint) +
  ylim(-10,10)+
  xlim(-10,10)+
  labs( x ="Log2FC MEKi", y = "Log2FC PI3Ki-MEKi") +
  theme(text = element_text(size=10))
        

ggsave("Figures/ph_PD_PIPD_corr.png",
       width = 10,
       height = 5,
       units = "cm")

cor(PD_and_PIPD$LogFC.PD, PD_and_PIPD$LogFC.PIPD, method=c("pearson"))
cor(PD_and_PIPD$LogFC.PD[PD_and_PIPD$Timepoint == "0.5h"], PD_and_PIPD$LogFC.PIPD[PD_and_PIPD$Timepoint == "0.5h"], method=c("pearson"))
cor(PD_and_PIPD$LogFC.PD[PD_and_PIPD$Timepoint == "2h"], PD_and_PIPD$LogFC.PIPD[PD_and_PIPD$Timepoint == "2h"], method=c("pearson"))
cor(PD_and_PIPD$LogFC.PD[PD_and_PIPD$Timepoint == "8h"], PD_and_PIPD$LogFC.PIPD[PD_and_PIPD$Timepoint == "8h"], method=c("pearson"))

PI_and_PIPD <- left_join(PI_DPPS, PIPD_DPPS, by = c("PTM_tp","Timepoint"),  suffix = c(".PI", ".PIPD")) 

cor(PI_and_PIPD$LogFC.PI, PI_and_PIPD$LogFC.PIPD, method=c("pearson"))
cor(PI_and_PIPD$LogFC.PI[PI_and_PIPD$Timepoint == "0.5h"], PI_and_PIPD$LogFC.PIPD[PI_and_PIPD$Timepoint == "0.5h"], method=c("pearson"))
cor(PI_and_PIPD$LogFC.PI[PI_and_PIPD$Timepoint == "2h"], PI_and_PIPD$LogFC.PIPD[PI_and_PIPD$Timepoint == "2h"], method=c("pearson"))
cor(PI_and_PIPD$LogFC.PI[PI_and_PIPD$Timepoint == "8h"], PI_and_PIPD$LogFC.PIPD[PI_and_PIPD$Timepoint == "8h"], method=c("pearson"))

ggplot(PI_and_PIPD, aes(x = LogFC.PI, y = LogFC.PIPD)) + 
  geom_hex(bins= 70) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500)) ) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1)+
  theme_minimal() +
  ylim(-10,10)+
  xlim(-10,10)+
  facet_grid(~Timepoint) +
  labs( x ="Log2FC PI3Ki", y = "Log2FC PI3Ki-MEKi") +
  theme(text = element_text(size=10))

ggsave("Figures/ph_PI_PIPD_corr.png",
       width = 10,
       height = 5,
       units = "cm")

# PI5Z --------------------------------------------------------------------

PIOXO_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "PI5Z") %>% 
  unite(PTM_tp, PTM_id, Timepoint ,remove = F)

PI_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "PI") %>% 
  unite(PTM_tp, PTM_id, Timepoint, remove = F) %>% 
  dplyr::filter(PTM_tp %in% PIOXO_DPPS$PTM_tp) 

OXO_DPPS <- DPPs_vs_DMSO_long %>% 
  dplyr::filter(Treatment == "5Z")%>% 
  unite(PTM_tp, PTM_id, Timepoint, remove = F ) %>% 
  dplyr::filter(PTM_tp %in% PIOXO_DPPS$PTM_tp) 

PI_and_PIOXO <- left_join(PI_DPPS, PIOXO_DPPS, by = c("PTM_tp","Timepoint"),  suffix = c(".PI", ".PIOXO")) 
OXO_and_PIOXO <- left_join(OXO_DPPS, PIOXO_DPPS, by = c("PTM_tp","Timepoint"),  suffix = c(".OXO", ".PIOXO")) 

ggplot(PI_and_PIOXO, aes(x = LogFC.PI, y = LogFC.PIOXO)) + 
  geom_hex(bins= 70) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500))) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1)+
  theme_minimal() +
  facet_grid(~Timepoint) +
  ylim(-10,10)+
  xlim(-10,10)+
  labs( x ="Log2FC PI3Ki", y = "Log2FC PI3Ki-TAKi") +theme(text = element_text(size=10))

ggsave("Figures/ph_PI_PI5Z_corr.png",
       width = 10,
       height = 5,
       units = "cm")

ggplot(OXO_and_PIOXO, aes(x = LogFC.OXO, y = LogFC.PIOXO)) + 
  geom_hex(bins = 70) +
  scale_fill_gradientn(colours = c("beige", "red", "darkblue"),values = scales::rescale(c(5, 100, 500))) +
  geom_smooth(method = "lm", colour = "red", lwd = 0.3) +
  geom_abline(intercept =0 , slope = 1)+
  theme_minimal(base_size = 11) +
  facet_grid(~Timepoint) +
  ylim(-10,10)+
  xlim(-10,10)+
  labs( x ="Log2FC TAKi", y = "Log2FC PI3Ki-TAKi" ) + theme(text = element_text(size=10))

ggsave("Figures/ph_5Z_PI5Z_corr.png",
       width = 10,
       height = 5,
       units = "cm")
cor(PI_and_PIOXO$LogFC.PI, PI_and_PIOXO$LogFC.PIOXO, method=c("pearson"))

cor(PI_and_PIOXO$LogFC.PI[PI_and_PIOXO$Timepoint == "0.5h"], PI_and_PIOXO$LogFC.PIOXO[PI_and_PIOXO$Timepoint == "0.5h"], method=c("pearson"))
cor(PI_and_PIOXO$LogFC.PI[PI_and_PIOXO$Timepoint == "2h"], PI_and_PIOXO$LogFC.PIOXO[PI_and_PIOXO$Timepoint == "2h"], method=c("pearson"))
cor(PI_and_PIOXO$LogFC.PI[PI_and_PIOXO$Timepoint == "8h"], PI_and_PIOXO$LogFC.PIOXO[PI_and_PIOXO$Timepoint == "8h"], method=c("pearson"))

cor(OXO_and_PIOXO$LogFC.OXO, OXO_and_PIOXO$LogFC.PIOXO, method=c("pearson"))

cor(OXO_and_PIOXO$LogFC.OXO[OXO_and_PIOXO$Timepoint == "0.5h"], OXO_and_PIOXO$LogFC.PIOXO[OXO_and_PIOXO$Timepoint == "0.5h"], method=c("pearson"))
cor(OXO_and_PIOXO$LogFC.OXO[OXO_and_PIOXO$Timepoint == "2h"], OXO_and_PIOXO$LogFC.PIOXO[OXO_and_PIOXO$Timepoint == "2h"], method=c("pearson"))
cor(OXO_and_PIOXO$LogFC.OXO[OXO_and_PIOXO$Timepoint == "8h"], OXO_and_PIOXO$LogFC.PIOXO[OXO_and_PIOXO$Timepoint == "8h"], method=c("pearson"))


# Corr plots per timepoint ------------------------------------------------

All_ratios_Limma <- read.csv("~/Documents/Git/AGS_synergies_multiomics/Phosphoproteomics/Data/All_ratios_Limma.txt")

PIPD_05 <- All_ratios_Limma %>% 
  select(logFC_PD_8h, logFC_PI_8h,logFC_PIPD_8h,PTM_collapse_key) %>% 
  column_to_rownames("PTM_collapse_key") %>% 
  as.data.frame()
PIPD_05[is.na(PIPD_05)] <- 0

res <- cor(PIPD_05)
corrplot(res,order = "original",type = "lower", addCoef.col = "black",cl.pos = "r")


PI5Z_05 <- All_ratios_Limma %>% 
  select(logFC_5Z_8h, logFC_PI_8h,logFC_PI5Z_8h,PTM_collapse_key) %>% 
  column_to_rownames("PTM_collapse_key") %>% 
  as.data.frame()

PI5Z_05[is.na(PI5Z_05)] <- 0

res_PI5Z <- cor(PI5Z_05)
corrplot(res_PI5Z,order = "original",type = "lower", diag = F, addCoef.col = "black", cl.pos = "r", cl.length = 5)
