library(tidyverse)

PI5Z_all <- read.csv("Transcriptomics/Output/PI5Z_specific_barplots.csv")
PIPD_all <- read.csv("Transcriptomics/Output/PIPD_specific_barplots.csv")

PI5Z_specific_DEGs <- read.csv("Transcriptomics/Output/PI5Z_specific_DEGs.csv")
PIPD_specific_DEGs <- read.csv("Transcriptomics/Output/PIPD_specific_DEGs.csv")

PIPD_all$Time <- factor(PIPD_all$Time, levels = c("1h","2h","4h","8h","24h"))
PI5Z_all$Time <- factor(PI5Z_all$Time, levels = c("1h","2h","4h","8h","24h"))


cls <- c("5Z_driven" = "#CC79A7", "PI_driven" = "#009E73","PI5Z_specific" = "#D55E00" )

ggplot(PI5Z_all, aes(x=Time, y=count, fill = Type)) + 
  scale_fill_manual(values= cls,
                    labels = c("TAKi-driven","PI3Ki-driven","PI3Ki-TAKi specific")) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Figures/tr_change_type_PI5Z.png",width = 15,
       height = 15,
       units = "cm",
       dpi = 300)

cls <- c("PD_driven" = "#E69F00", "PI_driven" = "#009E73","PIPD_specific" = "#0072B2" )
ggplot(PIPD_all, aes(x=Time, y=count, fill = Type)) + 
  scale_fill_manual(values= cls,
                    labels = c("MEKi-driven","PI3Ki-driven","PI3Ki-MEKi specific")) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal()+
  theme(legend.position = "bottom")

ggsave("Figures/tr_change_type_PIPD.png",width = 15,
       height = 15,
       units = "cm",
       dpi = 300)
intersect_by_time <- function(df1, df2) {
  # Split both dataframes by time
  df1_split <- split(df1, df1$Time)
  df2_split <- split(df2, df2$Time)
  # Create a empty list to store the results
  result <- list()
  # loop over the elements in both split dataframes and find the common gene 
  for(i in 1:length(df1_split)){
    name <- names(df1_split)[i]
    result[[name]] <-intersect(df1_split[[i]]$gene_name, df2_split[[i]]$gene_name)
  }
  # return the result
  return(result)
}

synergy_specific_overlap <- intersect_by_time(PI5Z_specific_DEGs, PIPD_specific_DEGs)
