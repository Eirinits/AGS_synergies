library(decoupleR)
library(OmnipathR)
library(tidyverse)
library(magrittr)
library(data.table)
library(ggh4x)
library(dorothea)

DEGs <- read.csv("Transcriptomics/Output/DEGs_long.csv")
gene_info <- read.delim("Transcriptomics/Data/gene_info.tsv")

long_to_wide  <- function(df) {
  dcast(setDT(df),
        source ~ Condition,
        value.var = c('score', 'pval'),
        fun.aggregate = mean)
}

run_decoupler <- function(mat,network){
  tryCatch(decouple(as.matrix(mat),network = network,.source='source', 
                    .target='target', args = list(wsum = list(times = 1000)),
                    minsize = 5), error=function(e) NULL)
  
}

CollecTRI <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)

DEGs_long <- DEGs %>% 
  split(~ Contrast)

DEGs_long <-lapply(DEGs_long, function(x) x[!duplicated(x$gene_name),])
DEGs_long <- lapply(DEGs_long, function(x) as.data.frame(x) %>% set_rownames(.$gene_name) %>% dplyr::select(t))
DEGs_long <- DEGs_long[sapply(DEGs_long, nrow)>5]

res_decoupler <- lapply(DEGs_long, function(x) run_ulm(x,CollecTRI))
save(res_decoupler, file="Transcriptomics/Output/res_decoupler_all_DEGs.RData")
load("Transcriptomics/Output/20240105_res_ulm_all_DEGs.RData")

res_decoupler_flt <-  res_decoupler %>% discard(is.null)
res_decoupler_flt <- lapply(res_decoupler_flt, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.05 & source %in% c(gene_info$gene_name,"NFKB","AP1") & condition == "t"))
res_decoupler_flt <- res_decoupler_flt[sapply(res_decoupler_flt, nrow)>0]

res_decoupler_flt <- Map(cbind, res_decoupler_flt, Condition = names(res_decoupler_flt)) %>%
  purrr::reduce(rbind) %>% 
  mutate(Condition = gsub("Diff_","",Condition)) %>% 
  separate(Condition,c("Drug","Time"), sep = "_", remove = F) %>% 
  rename(pval = "p_value") %>% 
#  mutate(score = round(score,2))%>% 
  select(source,Condition,score,pval)  %>% 
  filter(abs(score) > 2)

tf_n <- res_decoupler_flt %>% 
  group_by(Condition) %>% 
  summarise(count = n()) %>% 
  separate(Condition, c("Drug","Time"), remove = F)

tf_n <- tf_n %>% 
  mutate(Condition = factor(Condition, levels = c("OXO_2h","OXO_4h","OXO_8h","OXO_24h",
                                                  "PD_4h","PD_8h","PD_24h",
                                                  "PI_2h","PI_4h","PI_8h","PI_24h",
                                                  "PIOXO_2h","PIOXO_4h","PIOXO_8h","PIOXO_24h",
                                                  "PIPD_2h","PIPD_4h","PIPD_8h","PIPD_24h"),
                            labels = c("2h","4h","8h","24h",
                                       "4h","8h","24h",
                                       "2h","4h","8h","24h",
                                       "2h","4h","8h","24h",
                                       "2h","4h","8h","24h")))
  
group.colors <- c("OXO" = "#CC79A7",
                  "PD" = "#E69F00",
                  "PI" = "#009E73",
                  "PIOXO" = "#D55E00",
                  "PIPD" = "#0072B2")
drug_labels <- c("TAKi", "MEKi","PI3Ki","PI3Ki-TAKi","PI3Ki-MEKi")
names(drug_labels) <- c("OXO","PD","PI","PIOXO","PIPD")

ggplot(tf_n, aes(x=Condition, y=count), fill = Drug) + 
  geom_bar(stat = "identity") +
#  geom_boxplot()+
  ylab("No. of TFs") +
  facet_grid(~ Drug, scales = "free_x")+
 # facet_grid2(~Drug, scales = "free",  labeller=as_labeller(drug_labels)) +
  theme(strip.background = element_rect(colour = "black", fill = group.colors)) +
  theme_minimal()

ggsave("Figures/No_of_TFs.png",
       width = 15,
       height = 5,
       units = "cm",
       dpi = 300)

ggplot(res_decoupler_flt, aes(x=Condition, y=score)) + 
#  geom_bar(stat = "identity") +
  geom_boxplot()+
#  facet_grid(~ Drug, scales = "free_x")+
#  facet_grid2(~Drug, scales = "free", strip = strip, labeller=as_labeller(drug_labels)) +
  theme(strip.background = element_rect(colour = "black", fill = group.colors))

plot_TFs<- function(drugs_vec, res_df){
  res_decoupler_wide <- res_df  %>%
    separate(Condition,c("Drug","Time"), sep = "_", remove = F) %>% 
    filter(Drug %in% drugs_vec) %>% 
    select(-Drug,-Time) %>% 
    long_to_wide(.) 
  
  res_decoupler_wide$Type <- 0
  
  columns_to_check <- grep("^score_(PI|PD|OXO)_", names(res_decoupler_wide), value = TRUE)
  
  res_decoupler_wide$Type <- apply(res_decoupler_wide, 1, function(row) {
    if(all(is.na(row[columns_to_check]))) {
      return(1)
    } else {
      return(0)
    }
  })
  
  res_decoupler_long <- gather(res_decoupler_wide, v, value, -c(source,Type)) %>% 
    separate(v, c("var", "Drug","Time"))  %>% 
    spread(var, value) %>% 
    filter(!is.na(pval) & !is.na(score)) %>% 
    unite(Condition, Drug, Time, remove = F) %>% 
    mutate(Time = gsub("h","",Time)) %>% 
    mutate(Time = factor(Time, levels = c("2","4","8","24")))
  
  if ("OXO" %in% drugs_vec && "PD" %in% drugs_vec) {
    group.colors <- c("OXO" = "#CC79A7",
                      "PD" = "#E69F00",
                      "PI" = "#009E73",
                      "PIOXO" = "#D55E00",
                      "PIPD" = "#0072B2")
  } else if ("PD" %in% drugs_vec) {
    group.colors <- c("PD" = "#E69F00",
                      "PI" = "#009E73",
                      "PIPD" = "#0072B2")
  }else {
    group.colors <- c("PI" = "#009E73",
                      "OXO" = "#CC79A7",
                      "PIOXO" = "#D55E00")
  }

  drug_labels <- c("TAKi", "MEKi","PI3Ki","PI3Ki-TAKi","PI3Ki-MEKi")
  names(drug_labels) <- c("OXO","PD","PI","PIOXO","PIPD")
  
  strip <- strip_themed(background_x = elem_list_rect(fill = group.colors))
  
  ggplot(res_decoupler_long,aes(x=Time, y = reorder(source,Type, decreasing = T),color = score, size = -log10(pval))) +
  #  geom_rect(color = NA, aes(fill = Drug),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) +
    scale_fill_manual(values = group.colors, guide ="none")+
    scale_colour_gradient2(low = "navyblue",mid = "white",
                           high = "red", 
                           midpoint = 0, limits = c(-7,5)) +
    geom_point() + 
    geom_point(shape = 1,colour = "black")+
    facet_grid(~Drug, scales = "free") +
   # facet_grid2(~Drug, scales = "free", strip = strip, labeller=as_labeller(drug_labels)) +
    ylab("Transcription Factors") + 
    theme(
   #   panel.background = element_rect(fill = NA),
  #    panel.ontop = F,
      legend.direction="horizontal",
      legend.position = "bottom"
    )+
    theme(strip.background = element_rect(colour = "black", fill = group.colors))
  }
p <-plot_TFs(c("PI","PD","PIPD","OXO","PIOXO"),res_decoupler_flt)
p

plot_TFs(c("PIPD","OXO","PIOXO"),res_decoupler_flt)



ggsave("Figures/20230921_TFs_all_DEGs_noPanels.png",
       width = 15,
       height = 45,
       units = "cm",
       dpi = 300)

top <- res_decoupler_flt %>% 
  filter(grepl("PIOXO|PIPD", Condition)) %>% 
  top_n(n = 10,wt = score) %>% 
  distinct(source)

bottom <- res_decoupler_flt %>% 
  filter(grepl("PIOXO|PIPD", Condition)) %>% 
  top_n(n = -10,wt = score) %>% 
  distinct(source)

sel_TFs <- res_decoupler_flt %>% 
  filter(source %in% top$source | source %in% bottom$source)

plot_TFs(c("PI","PD","PIPD","OXO","PIOXO"),sel_TFs)

ggsave("Figures/TF_activities_top_TFs.png",
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300)

plot_TFs(c("PI","OXO","PIOXO"),res_decoupler_flt)

ggsave("Figures/20230921_TFs_Pi3KTAK.png",
       width = 10,
       height = 33,
       units = "cm",
       dpi = 300)

plot_TFs(c("PI","PD","PIPD"),res_decoupler_flt)


ggsave("Figures/20230921_TFs_Pi3KMEK.png",
       width = 10,
       height = 33,
       units = "cm",
       dpi = 300)
