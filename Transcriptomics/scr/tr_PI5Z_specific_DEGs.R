library(tidyverse)
library(gridExtra)
library(data.table)
library(purrr)

long_to_wide  <- function(df) {
  dcast(setDT(df),
        gene_name ~ Condition,
        value.var = c('logFC', 'adj.P.Val'))
}

delete.na <- function(DF, n = 0) {
  DF[rowSums(is.na(DF)) <= n, ]
}

All_ratios_Limma <- read.csv("Transcriptomics/Output/DEGs_wide.csv")

PI5Z_long = All_ratios_Limma %>% 
  pivot_longer(-gene_name) %>% 
  separate(name, into = c("name", "Condition"), sep = "_", extra = "merge")  %>% 
  filter(!is.nan(value)) %>% 
  pivot_wider(names_from = "name", values_from = "value") %>% 
  filter(adj.P.Val < 0.05) %>% 
  separate(Condition, into = c("Drug","Time"), remove = F,sep = "_") %>% 
  filter(!grepl("PD",Drug)) %>% 
  split(~gene_name) 

# Remove p-sites that are regulated not regulated in the synergy
PI5Z_regulated = lapply(PI5Z_long, function(d) {
  if(!("PI5Z" %in% d$Drug)) {
    d = NULL
  }
  return(d)
}) 

PI5Z_regulated <- PI5Z_regulated %>%
  discard(is.null)

## Select p-sites that are only regulated in the synergies in at least one time point.
PI5Z_specific = lapply(PI5Z_regulated, function(d) {
  if(length(unique(d$Drug)) > 1) {
    d = NULL
  }
  return(d)
}) 

PI5Z_specific <- PI5Z_specific %>%
  discard(is.null) %>% 
  do.call("rbind", .)

# Select p-sites that are regulated in more than 2 conditions (synergy and at least one more drug).

mult_cond = lapply(PI5Z_regulated, function(d) {
  if(length(unique(d$Drug)) == 1) {
    d = NULL
  }
  return(d)
}) 

mult_cond <- mult_cond %>%
  discard(is.null) 

# Select p-sites that are regulated in a single time point in the synergy and at least one more drug.

single_tp = lapply(mult_cond, function(d) {
  if(length(unique(d$Time)) > 1) {
    d = NULL
  }
  return(d)
}) 

single_tp <- single_tp %>%
  discard(is.null) %>% 
  do.call("rbind", .) %>% 
  split(., ~ Time)

PI5Z_wide <- lapply(single_tp, function(x) long_to_wide(x) %>%
                      dplyr::select(-contains("adj")) %>% 
                      delete.na(.,2)  %>% 
                      pivot_longer(-gene_name) %>% 
                      mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
                      select(-value) %>% 
                      spread(., name, Sign) %>% 
                      rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
                                  .cols = starts_with("logFC")))


PI5Z_wide <- lapply(PI5Z_wide, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

PI5Z_wide <- map2(PI5Z_wide, names(PI5Z_wide), ~ mutate(.x, Timepoint = .y))

PI5Z_wide <- do.call("rbind", PI5Z_wide)

synergy_specific <- c("Up_Unaffected_Down","Unaffected_Up_Down",
                      "Down_Unaffected_Up","Unaffected_Down_Up",
                      "Up_Up_Down","Down_Down_Up")


PI5Z_specific_opposite <- subset(PI5Z_wide, `5Z_PI_PI5Z` %in% synergy_specific)

coordinated <- c("Down_Down_Down","Up_Up_Up")
coordinated_single_tp <- subset(PI5Z_wide, `5Z_PI_PI5Z` %in% coordinated)

# Qualitative -------------------------------------------------------------
# 
# higher_than_both_PI5Z <- lapply(single_tp, function(x){
#   x %>%
#     group_by(gene_name) %>% 
#     filter(n() >= 3) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "PI"]) & abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "5Z"])) 
# })
# 
# higher_than_PI <- lapply(single_tp, function(x){
#   x %>% 
#     group_by(gene_name) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PI5Z", "PI")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "PI"])) 
# })
# 
# higher_than_5Z <- lapply(single_tp, function(x){
#   x %>%
#     group_by(gene_name) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PI5Z", "5Z")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "5Z"])) 
# })


# Multiple timepoints  ----------------------------------------------------

mult_cond_split <- do.call("rbind", mult_cond)  
mult_cond_split <- mult_cond_split %>% 
  unite(Split_fct,gene_name, Time, sep = "_", remove = F) %>% 
  split(~Split_fct)

mult_cond_split <- lapply(mult_cond_split, function(x) long_to_wide(x) %>%
                            dplyr::select(-contains("adj")) %>% 
                            pivot_longer(-gene_name) %>% 
                            mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
                            select(-value) %>% 
                            spread(., name, Sign) %>% 
                            rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
                                        .cols = starts_with("logFC")))

mult_cond_named <- do.call("bind_rows",mult_cond_split)

mult_cond_named_1h <- mult_cond_named %>% 
  select(gene_name, contains("1h")) %>% 
  delete.na(2) %>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_1h <- mult_cond_named_1h[,c(1,4,3,2)]

mult_cond_named_2h <- mult_cond_named %>% 
  select(gene_name, contains("2h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_2h <- mult_cond_named_2h[,c(1,4,2,3)]

mult_cond_named_4h <- mult_cond_named %>% 
  select(gene_name, contains("_4h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_4h <- mult_cond_named_4h[,c(1,4,3,2)]

mult_cond_named_8h <- mult_cond_named %>% 
  select(gene_name, contains("8h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_8h <- mult_cond_named_8h[,c(1,2,4,3)]

mult_cond_named_24h <- mult_cond_named %>% 
  select(gene_name, contains("24h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_24h <- mult_cond_named_24h[,c(1,2,4,3)]

mult_cond_named <- list("1h" = mult_cond_named_1h, "2h" = mult_cond_named_2h, "4h" = mult_cond_named_4h,"8h" = mult_cond_named_8h, "24h" = mult_cond_named_24h)

mult_cond_united <- lapply(mult_cond_named, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

mult_cond_united <- map2(mult_cond_united, names(mult_cond_united), ~ mutate(.x, Timepoint = .y))

mult_cond_united <- do.call("rbind", mult_cond_united)

coordinated_multi_tp <- subset(mult_cond_united, `5Z_PI_PI5Z` %in% coordinated)

changes_PI5Z <- read.csv("Phosphoproteomics/Output/changes_PI5Z.csv")
mult_cond_united <- left_join(mult_cond_united, changes_PI5Z, by = c("5Z_PI_PI5Z" = "Change_type"))

mult_cond_united_named <-  mult_cond_united %>% 
  select(-`5Z_PI_PI5Z`) %>% 
  spread(., Timepoint, Grouped_changes_I)  %>% 
  split(~gene)

mult_cond_united_named <- lapply(mult_cond_united_named, function(x) {
  x <- x[,c("gene","1h","2h","4h","8h","24h")]
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:6], collapse = '_')] <- do.call(paste, c(x[2:6], sep = '_'))
  x <- x[,c(1,7)]
  x
})

mult_cond_all_tps <- do.call("rbind",mult_cond_united_named) 

PI5Z_specific_mult_tps <- mult_cond_all_tps %>% 
  filter(!grepl("Single|Same", `1h_2h_4h_8h_24h`)) 

mult_cond_all_tps_sep <- mult_cond_all_tps %>% 
  separate(`1h_2h_4h_8h_24h`, into = c("X1h","X2h","X4h","X8h","X24h"), sep = "_") %>% 
  filter(!gene %in% PI5Z_specific_mult_tps$gene)

# Single drug responses ---------------------------------------------------

# Single drug-driven in a single time point

PI_driven_names <- c("Unaffected_Up_Up","Unaffected_Down_Down","Down_Up_Up","Up_Down_Down")
PI_specific_names <- c("Unaffected_Up_Unaffected","Unaffected_Down_Unaffected")
PI_driven <- subset(PI5Z_wide, `5Z_PI_PI5Z` %in% PI_driven_names)

X5Z_driven_names <- c("Up_Unaffected_Up","Down_Unaffected_Down","Up_Down_Up","Down_Up_Down")
X5Z_specific_names <- c("Up_Unaffected_Unaffected","Down_Unaffected_Unaffected")
X5Z_driven <- subset(PI5Z_wide, `5Z_PI_PI5Z` %in% X5Z_driven_names)

# Single drug-driven in multiple time points

SD_driven_mult <- mult_cond_united
SD_driven_mult <- SD_driven_mult %>% 
  mutate(Grouped_changes_I = ifelse( `5Z_PI_PI5Z` %in% PI_driven_names, "PI_driven", Grouped_changes_I)) %>% 
  mutate(Grouped_changes_I = ifelse( `5Z_PI_PI5Z` %in% X5Z_driven_names, "5Z_driven", Grouped_changes_I)) %>% 
  mutate(Grouped_changes_I = ifelse( `5Z_PI_PI5Z` %in% PI_specific_names, "PI_specific", Grouped_changes_I)) %>% 
  mutate(Grouped_changes_I = ifelse( `5Z_PI_PI5Z` %in% X5Z_specific_names, "5Z_specific", Grouped_changes_I)) 

SD_driven_mult <-  SD_driven_mult %>% 
  select(- `5Z_PI_PI5Z`) %>% 
  spread(., Timepoint, Grouped_changes_I)  %>% 
  split(~gene)

SD_driven_mult <- lapply(SD_driven_mult, function(x) {
  x <- x[,c("gene","1h","2h","4h","8h","24h")]
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:6], collapse = '_')] <- do.call(paste, c(x[2:6], sep = '_'))
  x <- x[,c(1,7)]
  x
})

SD_all_tps <- do.call("rbind",SD_driven_mult)   

PI_driven_mult <- SD_all_tps %>% 
  filter(!grepl("Synergy|5Z|Same|PI_specific|Level|Single",`1h_2h_4h_8h_24h`))

X5Z_driven_mult <- SD_all_tps%>% 
  filter(!grepl("Synergy|PI|Same|5Z_specific|Level|Single",`1h_2h_4h_8h_24h`))

# Coordinated changes

Same_in_all <- SD_all_tps %>% 
  filter(!grepl("PI_|5Z_|Single|Synergy|driven|Level", `1h_2h_4h_8h_24h`))

PI5Z_specific_all <- do.call("rbind", PI5Z_long) %>% 
  filter(gene_name %in% PI5Z_specific$gene_name | gene_name %in% PI5Z_specific_mult_tps$gene | gene_name %in% PI5Z_specific_opposite$gene) %>% 
  filter(Drug == "PI5Z") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PI5Z_specific")

PI_driven_all <- do.call("rbind", PI5Z_long) %>% 
  filter(gene_name %in% PI_driven$gene | gene_name %in% PI_driven_mult$gene) %>% 
  filter(Drug == "PI") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Type = "PI_driven")

X5Z_driven_all <- do.call("rbind", PI5Z_long) %>% 
  filter(gene_name %in% X5Z_driven$gene | gene_name %in% X5Z_driven_mult$gene) %>% 
  filter(Drug == "5Z") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "5Z_driven")

Same_all <- do.call("rbind", PI5Z_long) %>% 
  filter(gene_name %in% Same_in_all$gene) %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "Coordinated")
coordinated_multi_tp <- coordinated_multi_tp[,1:3]
Coordinated_all <- rbind(coordinated_multi_tp, coordinated_single_tp) %>% 
  group_by(Timepoint) %>% 
  summarize(count = n()) %>% 
  rename(Time = "Timepoint") %>% 
  mutate(Type = "Concerted")

PI5Z_all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(PI5Z_specific_all, PI_driven_all, X5Z_driven_all,Coordinated_all))
PI5Z_all$Time <- factor(PI5Z_all$Time, levels = c("1h","2h","4h","8h","24h"))
cls <- c("5Z_driven" = "#CC79A7", "PI_driven" = "#009E73","PI5Z_specific" = "#D55E00", Concerted = "darkgray")
ggplot(PI5Z_all, aes(x=Time, y=count, fill = Type)) + 
  scale_fill_manual(values= cls) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal()

ggsave("PI5Z_specific_DEGs_concerted.png", dpi = 300, height = 5, width = 5)
