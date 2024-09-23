library(decoupleR)
library(Omnipa  )
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
load("Transcriptomics/Output/res_decoupler_all_DEGs.RData")

res_decoupler_flt <-  res_decoupler %>% discard(is.null)
res_decoupler_flt <- lapply(res_decoupler_flt, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.05 & source %in% c(gene_info$gene_name,"NFKB","AP1") & condition == "t"))
res_decoupler_flt <- res_decoupler_flt[sapply(res_decoupler_flt, nrow)>0]

res_decoupler_flt <- Map(cbind, res_decoupler_flt, Condition = names(res_decoupler_flt)) %>%
  purrr::reduce(rbind) %>% 
  mutate(Condition = gsub("Diff_","",Condition)) %>% 
  separate(Condition,c("Drug","Time"), sep = "_", remove = F) %>% 
  rename(pval = "p_value") %>% 
#  mutate(score = round(score,2))%>% 
  select(source,Condition,score,pval)  

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
    facet_grid2(~Drug, scales = "free", strip = strip, labeller=as_labeller(drug_labels)) +
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

ggsave("Figures/20230921_TFs_all_DEGs_noPanels.png",
       width = 15,
       height = 45,
       units = "cm",
       dpi = 300)

sel_TFs <- res_decoupler_flt %>% 
  filter(source %in% c("MYC", "FOXO3","FOXO1","FOXO4", "FOSB","E2F1","E2F2","E2F3","E2F4", "CREB1", "ATF4", "STAT3","AP1","RUNX3","JUN", "NFKB","REL"))
plot_TFs(c("PI","PD","PIPD","OXO","PIOXO"),sel_TFs)

ggsave("Figures/TF_activities_sel_TFs.png",
       width = 13,
       height = 17,
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
