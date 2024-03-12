library(tidyverse)
library(decoupleR)
library(magrittr)

results_decoupler <- list()
run_decoupler_in_list <- function(DEGs_df, network, condition_name){
  df_list <- DEGs_df %>% 
    dplyr::select(gene_name, logFC,Time)  %>% 
    split(~Time)
  
  df_list <-lapply(df_list, function(x) x[!duplicated(x$gene_name),])
  df_list <- lapply(df_list, function(x) as.data.frame(x) %>% set_rownames(.$gene_name) %>% dplyr::select(-gene_name)%>% dplyr::select(logFC))
  
  res_decoupler <- lapply(df_list, function(x) decoupler(x,network))
  
  res_decoupler <-  res_decoupler %>% discard(is.null)
  res_decoupler <- lapply(res_decoupler, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.05))
  res_decoupler <- res_decoupler[sapply(res_decoupler, nrow)>0]
  
  res_decoupler <- Map(cbind, res_decoupler, Condition = names(res_decoupler))
  
  res_decoupler_df <- res_decoupler  %>%
    purrr::reduce(rbind) 
  
  results_decoupler[[condition_name]] <<- res_decoupler_df
  
  res_decoupler_df %>%
    mutate(Condition = factor(Condition,levels = c("2h","4h","8h","24h"))) %>% 
    ggplot(aes(Condition, source)) +
    geom_point(aes(colour = score, size = -log10(p_value))) +
    scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                          high="red")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = "Estimated TF Activities", subtitle = "FDR < 0.05 \nmin. 5 targets") + 
    ylab("Transcription Factor")
}

decoupler <- function(mat,network){
  tryCatch(run_ulm(as.matrix(mat), 
                   network = network,
                   .source = "tf",
                   minsize = 5), error=function(e) NULL)
  
}

# working directory
PI5Z_specific_DEGs <- read.csv("Transcriptomics/Output/PI5Z_specific_DEGs.csv")
PIPD_specific_DEGs <- read.csv("Transcriptomics/Output/PIPD_specific_DEGs.csv")
CollecTRI_20221020 <- read.delim("Transcriptomics/Data/External/CollecTRI_20221020.tsv")

PI5Z_specific_DEGs$Time <- factor(PI5Z_specific_DEGs$Time, levels = c("1h","2h","4h","8h","24h"))
PIPD_specific_DEGs$Time <- factor(PIPD_specific_DEGs$Time, levels = c("1h","2h","4h","8h","24h"))

run_decoupler_in_list(PIPD_specific_DEGs, CollecTRI_20221020, "PIPD")
run_decoupler_in_list(PI5Z_specific_DEGs, CollecTRI_20221020, "PI5Z")
