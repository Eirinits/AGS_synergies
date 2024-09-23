library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

counts_norm <- read.csv("Transcriptomics/Data/gene_counts_flt_norm_wCTRL.csv") 
counts_long <-  merge(counts_norm,gene_info,by="row.names",all.x=TRUE) %>% 
  filter(gene_name != "") %>% 
  distinct(gene_name,.keep_all = T) %>% 
  select(contains("Rep"),"gene_name") %>% 
  pivot_longer(cols = contains("Rep"), names_to = "ID", values_to = "counts") %>% 
  mutate(Condition = vapply(strsplit(ID, '_'), function(x) paste(x[seq.int(2)], collapse='_'), character(1L))) %>% 
  split(~Condition)

# Collapse replicates and scale
gene_counts_scaled <- map(counts_long, ~ .x %>%
                            group_by(gene_name) %>%
                            mutate(counts_scaled = (counts - mean(counts))/sd(counts)) %>%
                            group_by(gene_name, Condition) %>%
                            summarise(mean_counts = mean(counts_scaled)) %>%
                            ungroup() %>%
                            filter(!is.na(mean_counts)))

gene_counts_sc_mat <- bind_rows(gene_counts_scaled)

gene_counts_sc_mat <- gene_counts_sc_mat %>%
                      spread(., gene_name, mean_counts, drop = TRUE) %>%
                      t() %>%
                      row_to_names(row_number = 1) %>% 
  as.data.frame()
gene_counts_sc_mat <- mutate_all(gene_counts_sc_mat, function(x) as.numeric(as.character(x))) %>% 
  as.matrix()

design <- data.frame(Sample = colnames(gene_counts_sc_mat), Condition = colnames(gene_counts_sc_mat))  
design$Condition <- vapply(strsplit(design$Sample, '_'), function(x)
    paste(x[seq.int(2)], collapse='_'), character(1L))
net <- get_progeny(organism = 'human', top = 100)

# Transform to wide matrix


sample_acts <- run_wmean(mat=gene_counts_sc_mat, net=net, .source='source', .target='target',
                         .mor='weight', times = 100, minsize = 5)
sample_acts
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  separate(condition, c("drug","time"),remove = F) %>% 
  mutate(time = factor(time, levels = c("0","1","2","4","8","24"))) %>% 
  arrange(time) %>% 
  column_to_rownames('condition') %>%
  select(-drug,-time) %>% 
  as.matrix()
# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

PI5Z <- sample_acts_mat %>% 
  as.data.frame() %>% 
  filter(grepl("PI_|OXO", rownames(.))) %>% 
  as.matrix()%>%
  t()
pheatmap(PI5Z, border_color = NA, color=my_color, breaks = my_breaks, cluster_cols =  F) 


PIPD <- sample_acts_mat %>% 
  as.data.frame() %>% 
  filter(grepl("PI_|PD", rownames(.))) %>% 
  as.matrix()%>%
  t()
pheatmap(PIPD, border_color = NA, color=my_color, breaks = my_breaks, cluster_cols =  F) 
