# Load libraries & data ---------------------------------------------------

library(tidyverse)
library(limma)
library(edgeR)
library(stringi)
library(UpSetR)
library(grid)
library(vsn)


gene_counts <- read.table("Transcriptomics/Data/Raw_data/gene_counts.tsv", row.names = 1, sep = "\t")
features <- read.delim("Transcriptomics/Data/gene_info.tsv", sep="\t")
features <- features[,c('gene_name','gene_biotype')]

samples <- read.table("Transcriptomics/Data/Raw_data/sample_info.tsv", sep = "\t", header = TRUE)

samples <- samples %>% 
  select(Sample_ID, Timepoint, treatment) %>% 
  mutate(Sample_ID = paste0("X", Sample_ID)) %>% 
  unite(Condition, treatment, Timepoint, sep = "_", remove = F) %>% 
  mutate(Condition = gsub("5Z","OXO",Condition)) %>% 
  mutate(Condition = gsub("OXO_PI","PIOXO",Condition)) %>% 
  mutate(Condition = gsub("PD_PI","PIPD",Condition)) 

samples <- samples[samples$Sample_ID != "X101",]
samples <- samples[samples$Sample_ID != "X63",]

count_df <- gene_counts[,samples$Sample_ID]
lev <- unique(samples$Condition)
fac <- factor(samples$Condition, levels=lev)
design <- model.matrix(~ 0 + fac)
colnames(design) <- gsub("fac","",colnames(design))

table(rowSums(count_df==0)==94)

x <- DGEList(counts=count_df)
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs, keep.lib.sizes=FALSE]
dim(x)

## Normalized counts for visualization
counts_flt <- x$counts
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

fit <- vsnMatrix(as.matrix(x$counts)) #train vsn parameters

meanSdPlot(fit)

gene_counts_vsn <- as.data.frame(vsn::predict(fit,as.matrix(x$counts)))
#write.table(gene_counts_vsn, "Transcriptomics/Data/gene_counts_flt_norm_wCTRL.csv", quote = F)

## Differential expression with limma-voom
v <- voom(x, design, plot=TRUE)

contrasts_list <- makeContrasts(
  Diff_OXO_1h =  (OXO_1) - (DMSO_1),
  Diff_PI_1h = (PI_1) - (DMSO_1),
  Diff_PD_1h = (PD_1) - (DMSO_1),
  Diff_PIPD_1h = (PIPD_1) - (DMSO_1),
  Diff_PIOXO_1h = (PIOXO_1)  - (DMSO_1),
  Diff_OXO_2h =  (OXO_2) - (DMSO_2),
  Diff_PI_2h =  (PI_2) - (DMSO_2),
  Diff_PD_2h =  (PD_2) - (DMSO_2),
  Diff_PIPD_2h =  (PIPD_2) - (DMSO_2),
  Diff_PIOXO_2h =  (PIOXO_2) - (DMSO_2),
  Diff_OXO_4h =  (OXO_4) - (DMSO_4),
  Diff_PI_4h =  (PI_4) - (DMSO_4),
  Diff_PD_4h = (PD_4) - (DMSO_4),
  Diff_PIPD_4h = (PIPD_4) - (DMSO_4),
  Diff_PIOXO_4h =  (PIOXO_4) - (DMSO_4),
  Diff_OXO_8h =  (OXO_8) - (DMSO_8),
  Diff_PI_8h =  (PI_8) - (DMSO_8),
  Diff_PD_8h =  (PD_8) - (DMSO_8),
  Diff_PIPD_8h =  (PIPD_8) - (DMSO_8),
  Diff_PIOXO_8h =  (PIOXO_8) - (DMSO_8),
  Diff_OXO_24h = (OXO_24) - (DMSO_24),
  Diff_PI_24h = (PI_24) - (DMSO_24),
  Diff_PD_24h = (PD_24) - (DMSO_24),
  Diff_PIPD_24h = (PIPD_24) - (DMSO_24),
  Diff_PIOXO_24h = (PIOXO_24) - (DMSO_24),
  levels=design)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrasts_list)
efit <- eBayes(vfit)

extract_DEGs <- function(x1){
  contr_name <- colnames(efit$coefficients)[[x1]]
  tt <- topTreat(efit, x1, number = Inf) %>%
    rownames_to_column("Gene_name") %>%
    mutate(Contrast = contr_name) %>%
    as_tibble() 
  tt
}

Coeffs <- dimnames(efit$coefficients)[["Contrasts"]]

DEGs_perContrast <- map_df(seq_along(Coeffs) , extract_DEGs) %>% 
  filter(adj.P.Val < 0.05)

DEGs_perContrast_all <- map_df(seq_along(Coeffs) , extract_DEGs) 

DEGs_n <- DEGs_perContrast %>%
  dplyr::count(Contrast) %>% 
  separate(Contrast,into = c("Diff","Drug","Time"),"_", remove = F)

DEGs_n$Drug <- gsub("OXO","5Z", DEGs_n$Drug)
x_order <- factor(DEGs_n$Time,levels = c("1h","2h","4h","8h","24h"))

ggplot(DEGs_n,aes(fill=Time,  y=n, x=x_order)) +
  geom_bar(stat="identity") +
  facet_wrap(Drug ~ ., scales = "free_x", ncol = 5) + 
  ylab("No. of differentially expressed genes") +
  xlab("Timepoint") +
  scale_fill_discrete(limits = c("1h","2h","4h","8h","24h")) +
  ggtitle("No. of differentially expressed genes", subtitle = "Limma results - FDR < 0.05 and abs(logFC) > 1")

write.csv(DEGs_perContrast_all, file = "Transcriptomics/Output/20230904_DEGs_perContrast_all.csv", row.names = F, quote = F)
write.csv(DEGs_perContrast, file = "Transcriptomics/Output/20230904_DEGs_perContrast_FDR005_FC1.csv", row.names = F, quote = F)

# Overlap calculations ----------------------------------------------------
DEGs_perContrast$Up_Down <- ifelse(DEGs_perContrast$logFC < 0, "DOWN", "UP")
DEGs_perContrast_UP <- subset(DEGs_perContrast, Up_Down == "UP")
DEGs_perContrast_DW <- subset(DEGs_perContrast, Up_Down == "DOWN")

DEGs_list_UP <- with(DEGs_perContrast_UP, split(Gene_name, Contrast))
#DEGs_list_UP <- setNames(DEGs_list_UP, paste0(names(DEGs_list_UP),"_UP"))
DEGs_list_UP <- setNames(DEGs_list_UP, gsub("Diff_","",names(DEGs_list_UP)))

DEGs_list_DOWN <- with(DEGs_perContrast_DW, split(Gene_name, Contrast))
#DEGs_list_DOWN <- setNames(DEGs_list_DOWN, paste0(names(DEGs_list_DOWN),"_DOWN"))
DEGs_list_DOWN <- setNames(DEGs_list_DOWN, gsub("Diff_","",names(DEGs_list_DOWN)))

# upset(fromList(DEGs_list_UP), 
#       sets = c("OXO_1h","PI_1h","PD_1h","PIOXO_1h","PIPD_1h",
#                "OXO_2h","PI_2h","PD_2h","PIOXO_2h","PIPD_2h",
#                "OXO_4h","PI_4h","PD_4h","PIOXO_4h","PIPD_4h",
#                "OXO_8h","PI_8h","PD_8h","PIOXO_8h","PIPD_8h",
#                "OXO_24h","PI_24h","PD_24h","PIOXO_24h","PIPD_24h"))

upset(fromList(DEGs_list_UP), order.by = "freq", keep.order = T, nintersects = 20, sets = c("OXO_4h","PI_4h","PD_4h","PIOXO_4h","PIPD_4h",
                                                                                           "OXO_8h","PI_8h","PD_8h","PIOXO_8h","PIPD_8h",
                                                                                           "OXO_24h","PI_24h","PD_24h","PIOXO_24h","PIPD_24h"),
      text.scale = 2,
      main.bar.color = "darkblue")
grid.text("Overlap between upregulated genes",x = 0.65, y=0.95, gp = gpar(
  fontsize = 15))

upset(fromList(DEGs_list_DOWN), order.by = "freq", keep.order = T, nintersects = 20, sets = c("OXO_4h","PI_4h","PD_4h","PIOXO_4h","PIPD_4h",
                                                                                            "OXO_8h","PI_8h","PD_8h","PIOXO_8h","PIPD_8h",
                                                                                            "OXO_24h","PI_24h","PD_24h","PIOXO_24h","PIPD_24h"),
      text.scale = 2,
      main.bar.color = "darkblue")
grid.text("Overlap between downregulated genes",x = 0.65, y=0.95, gp = gpar(
  fontsize = 15))

library(dplyr)
DEGs_perContrast_FDR005_FC1_HUGOSYMBOLS %>% 
  dplyr::filter(Contrast == "Diff_PIOXO_24h") %>% 
  arrange(abs(logFC))
