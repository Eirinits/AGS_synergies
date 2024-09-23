library(tidyverse)
library(gridExtra)
library(data.table)
library(purrr)

select <- dplyr::select
long_to_wide  <- function(df) {
  dcast(setDT(df),
        PTM_collapse_key ~ Condition,
        value.var = c('logFC', 'adj.P.Val'))
}

delete.na <- function(DF, n = 0) {
  DF[rowSums(is.na(DF)) <= n, ]
}

All_ratios_Limma <- read.csv("Phosphoproteomics/Data/All_ratios_Limma.txt")

PI5Z_long = All_ratios_Limma %>% 
  pivot_longer(-PTM_collapse_key) %>% 
  separate(name, into = c("name", "Condition"), sep = "_", extra = "merge")  %>% 
  filter(!is.nan(value)) %>% 
  pivot_wider(names_from = "name", values_from = "value") %>% 
  filter(adj.P.Val < 0.05) %>% 
  separate(Condition, into = c("Drug","Time"), remove = F,sep = "_") %>% 
  filter(!grepl("PD",Drug)) %>% 
  split(~PTM_collapse_key) 

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
                      delete.na(.,2) %>% 
                      pivot_longer(-PTM_collapse_key) %>% 
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

# Qualitative -------------------------------------------------------------
# 
# higher_than_both_PI5Z <- lapply(single_tp, function(x){
#   x %>%
#     group_by(PTM_collapse_key) %>% 
#     filter(n() >= 3) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "PI"]) & abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "5Z"])) 
# })
# 
# higher_than_PI <- lapply(single_tp, function(x){
#   x %>% 
#     group_by(PTM_collapse_key) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PI5Z", "PI")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "PI"])) 
# })
# 
# higher_than_5Z <- lapply(single_tp, function(x){
#   x %>%
#     group_by(PTM_collapse_key) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PI5Z", "5Z")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PI5Z"]) > 2 * abs(logFC[Drug == "5Z"])) 
# })


# Multiple timepoints & Accelerated/delayed responses ----------------------------------------------------

# The following code characterizes each phosphosite as unaffected, down, up for each time point
# Then it uses a predefined file with hyper-categories (synergy-specific, single drug-driven etc) to map these changes
# Finally, it filters those phosphosites that are regulated uniquely in the synergy, in multiple timepoints

mult_cond_split <- do.call("rbind", mult_cond)  
mult_cond_split <- mult_cond_split %>% 
  unite(Split_fct,PTM_collapse_key, Time, sep = "_", remove = F) %>% 
  split(~Split_fct)

mult_cond_split <- lapply(mult_cond_split, function(x) long_to_wide(x) %>%
                            dplyr::select(-contains("adj")) %>% 
                            pivot_longer(-PTM_collapse_key) %>% 
                            mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
                            select(-value) %>% 
                            spread(., name, Sign) %>% 
                            rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
                                        .cols = starts_with("logFC")))

mult_cond_named <- do.call("bind_rows",mult_cond_split)

mult_cond_named_05h <- mult_cond_named %>% 
  select(PTM_collapse_key, contains("30m")) %>% 
  delete.na(2) %>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_05h <- mult_cond_named_05h[,c(1,2,4,3)]

mult_cond_named_2h <- mult_cond_named %>% 
  select(PTM_collapse_key, contains("2h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_2h <- mult_cond_named_2h[,c(1,3,4,2)]

mult_cond_named_8h <- mult_cond_named %>% 
  select(PTM_collapse_key, contains("8h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named <- list("05h" = mult_cond_named_05h, "2h" = mult_cond_named_2h, "8h" = mult_cond_named_8h)

mult_cond_united <- lapply(mult_cond_named, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

mult_cond_united <- map2(mult_cond_united, names(mult_cond_united), ~ mutate(.x, Timepoint = .y))

mult_cond_united <- do.call("rbind", mult_cond_united)

changes_PI5Z <- read.csv("Phosphoproteomics/Output/changes_PI5Z.csv")
mult_cond_united <- left_join(mult_cond_united, changes_PI5Z, by = c("5Z_PI_PI5Z" = "Change_type"))

mult_cond_united_named <-  mult_cond_united %>% 
  select(-`5Z_PI_PI5Z`) %>% 
  spread(., Timepoint, Grouped_changes_I)  %>% 
  split(~PTM)

mult_cond_united_named <- lapply(mult_cond_united_named, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

mult_cond_all_tps <- do.call("rbind",mult_cond_united_named) 

PI5Z_specific_mult_tps <- mult_cond_all_tps %>% 
  filter(!grepl("Single|Same|Level", `05h_2h_8h`)) 



# Single drug driven ---------------------------------------------------

# Select p-sites that are regulated only in one or both single drugs
# 
# SD_regulated = lapply(PI5Z_long, function(d) {
#   if(("PI5Z" %in% d$Drug)) {
#     d = NULL
#   }
#   return(d)
# }) 
# 
# SD_regulated <- SD_regulated %>%
#   discard(is.null) 
# 
# 
# SD_regulated <- lapply(SD_regulated, function(x) long_to_wide(x) %>%
#                          dplyr::select(-contains("adj")) %>% 
#                          delete.na(.,2) %>% 
#                          pivot_longer(-PTM_collapse_key) %>% 
#                          mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
#                          select(-value) %>% 
#                          spread(., name, Sign) %>% 
#                          rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
#                                      .cols = starts_with("logFC")))
# 
# SD_regulated <- do.call("bind_rows", SD_regulated)
# SD_regulated_05h <- SD_regulated %>% 
#   select(PTM_collapse_key, contains("30m")) %>% 
#   delete.na(1) %>% 
#   mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))
# 
# SD_regulated_2h <- SD_regulated %>% 
#   select(PTM_collapse_key, contains("2h")) %>% 
#   delete.na(1)%>% 
#   mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))
# 
# SD_regulated_2h <- SD_regulated_2h[,c(1,3,2)]
# SD_regulated_8h <- SD_regulated %>% 
#   select(PTM_collapse_key, contains("8h")) %>% 
#   delete.na(1)%>% 
#   mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))
# 
# SD_regulated_named <- list("05h" = SD_regulated_05h, "2h" = SD_regulated_2h, "8h" = SD_regulated_8h)
# 
# SD_regulated_named <- lapply(SD_regulated_named, function(x) {
#   colnames(x) <- sub("_.*", "", colnames(x))
#   x[paste(names(x)[2:3], collapse = '_')] <- do.call(paste, c(x[2:3], sep = '_'))
#   x <- x[,c(1,4)]
#   x
# })
# 
# SD_regulated_named <- map2(SD_regulated_named, names(SD_regulated_named), ~ mutate(.x, Timepoint = .y))
# 
# SD_regulated_named <- do.call("rbind", SD_regulated_named)
# 
# SD_regulated_named_split <-  SD_regulated_named %>% 
#   mutate(`5Z_PI_sign` = ifelse(`5Z_PI` %in% c("Up_Unaffected","Down_Unaffected"),"5Z_specific",`5Z_PI`)) %>% 
#   mutate(`5Z_PI_sign` = ifelse(`5Z_PI_sign` %in% c("Unaffected_Up","Unaffected_Down"),"PI_specific",`5Z_PI_sign`)) %>% 
# spread(., Timepoint, `5Z_PI_sign`) %>% 
# split(~PTM)
# 
# SD_regulated_named_split <- lapply(SD_regulated_named_split, function(x) {
#   colnames(x) <- sub("_.*", "", colnames(x))
#   x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
#   x <- x[,c(1,5)]
#   x
# })
# 
# SD_regulated_named_split_all <- do.call("rbind",SD_regulated_named_split) 
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
  split(~PTM)

SD_driven_mult <- lapply(SD_driven_mult, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

SD_all_tps <- do.call("rbind",SD_driven_mult)   

PI_driven_mult <- SD_all_tps %>% 
  filter(!grepl("5Z_|Same|Single|Synergy|specific|Level", `05h_2h_8h`)) %>% 
  filter(!PTM %in% PI_driven$PTM)

X5Z_driven_mult <- SD_all_tps %>% 
  filter(!grepl("PI_|Same|Single|Synergy|specific|Level", `05h_2h_8h`))%>% 
  filter(!PTM %in% X5Z_driven$PTM)

# Coordinated changes

Same_in_all <- SD_all_tps %>% 
  filter(!grepl("PI_|5Z_|Single|Synergy|driven|Level", `05h_2h_8h`))

PI5Z_specific_all <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% PI5Z_specific$PTM_collapse_key | PTM_collapse_key %in% PI5Z_specific_mult_tps$PTM | PTM_collapse_key %in% PI5Z_specific_opposite$PTM) %>% 
  filter(Drug == "PI5Z") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PI5Z_specific")

PI_driven_all <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% PI_driven$PTM | PTM_collapse_key %in% PI_driven_mult$PTM) %>% 
  filter(Drug == "PI") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PI_driven")

X5Z_driven_all <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% X5Z_driven$PTM | PTM_collapse_key %in% X5Z_driven_mult$PTM) %>% 
  filter(Drug == "5Z") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "5Z_driven")

Same_all <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% Same_in_all$PTM) %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "Coordinated")

PI5Z_all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(PI5Z_specific_all, PI_driven_all, X5Z_driven_all,Same_all))
PI5Z_all$Type = factor(PI5Z_all$Type, levels = c("Coordinated","PI_driven","5Z_driven","PI5Z_specific"))
cls <- c("5Z_driven" = "#CC79A7", "PI_driven" = "#009E73","PI5Z_specific" = "#D55E00","Coordinated" = "gray" )
PI5Z_all <- read_csv("Phosphoproteomics/Output/PI5Z_all_DPPs_quantified.csv")
ggplot(PI5Z_all, aes(x=Time, y=count, fill = Type)) + 
  scale_fill_manual(values= cls) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal() + theme(text = element_text(size=18))

ggsave("Figures/ph_PI5Z_specific_barplots.png",
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 300)

write.csv(PI5Z_all, file = "Phosphoproteomics/Output/PI5Z_all_DPPs_quantified.csv", row.names = F)

PI5Z_specific_all_names <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% PI5Z_specific$PTM_collapse_key | PTM_collapse_key %in% PI5Z_specific_mult_tps$PTM | PTM_collapse_key %in% PI5Z_specific_opposite$PTM) %>% 
  filter(Drug == "PI5Z") 

write.csv(PI5Z_specific_all_names, file = "Phosphoproteomics/Output/PI5Z_specific_DPPs.csv", row.names = F)

PI_driven_all_names <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% PI_driven$PTM | PTM_collapse_key %in% PI_driven_mult$PTM) %>% 
  filter(Drug != "5Z")
write.csv(PI_driven_all_names, file = "Phosphoproteomics/Output/PI_driven_DPPs_PI5Z", row.names = F)

X5Z_driven_all_names <- do.call("rbind", PI5Z_long) %>% 
  filter(PTM_collapse_key %in% X5Z_driven$PTM | PTM_collapse_key %in% X5Z_driven_mult$PTM) %>% 
  filter(Drug != "PI") 

write.csv(X5Z_driven_all_names, file = "Phosphoproteomics/Output/5Z_driven_DPPs_PI5Z", row.names = F)


# Leveled-out -------------------------------------------------------------

leveled_out = lapply(PI5Z_long, function(d) {
  if(("PI5Z" %in% d$Drug)) {
    d = NULL
  }
  return(d)
}) 

leveled_out <- leveled_out %>%
  discard(is.null)

## Select p-sites that are only regulated in the synergy in at least one time point.
leveled_out = lapply(leveled_out, function(d) {
  if(length(unique(d$Drug)) <= 1) {
    d = NULL
  }
  return(d)
}) 

leveled_out <- leveled_out %>%
  discard(is.null) %>% 
  do.call("rbind", .) 

leveled_out <- leveled_out %>%
  mutate(Change = ifelse(logFC < 0, "Down","Up"))  %>% 
  split(~ Time)

filter_PTMs <- function(df) {
  PTM_filtered <- df %>% 
    group_by(PTM_collapse_key) %>% 
    summarise(PI_change = Change[Drug == "PI"],
              X5Z_change = Change[Drug == "5Z"]) %>% 
    filter((PI_change == "Up" & X5Z_change == "Down") | 
             (PI_change == "Down" & X5Z_change == "Up")) %>% 
    inner_join(df, by = "PTM_collapse_key")
  return(PTM_filtered)
}

leveled_PI5Z <- lapply(leveled_out, filter_PTMs)

write.csv(X5Z_driven_all_names, file = "Phosphoproteomics/Output/leveledout_DPPs_PI5Z", row.names = F)

