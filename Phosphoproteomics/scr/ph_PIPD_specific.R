library(tidyverse)
library(gridExtra)
library(data.table)
library(purrr)

long_to_wide  <- function(df) {
  dcast(setDT(df),
        PTM_collapse_key ~ Condition,
        value.var = c('logFC', 'adj.P.Val'))
}

delete.na <- function(DF, n = 0) {
  DF[rowSums(is.na(DF)) <= n, ]
}

All_ratios_Limma <- read.csv("Phosphoproteomics/Data/All_ratios_Limma.txt")

PIPD_long = All_ratios_Limma %>% 
  pivot_longer(-PTM_collapse_key) %>% 
  separate(name, into = c("name", "Condition"), sep = "_", extra = "merge")  %>% 
  filter(!is.nan(value)) %>% 
  pivot_wider(names_from = "name", values_from = "value") %>% 
  filter(adj.P.Val < 0.05) %>% 
  separate(Condition, into = c("Drug","Time"), remove = F,sep = "_") %>% 
  filter(!grepl("5Z",Drug)) %>% 
  split(~PTM_collapse_key) 

# Remove p-sites that are not regulated in the synergy
PIPD_regulated = lapply(PIPD_long, function(d) {
  if(!("PIPD" %in% d$Drug)) {
    d = NULL
  }
  return(d)
}) 

PIPD_regulated <- PIPD_regulated %>%
  discard(is.null)

## Select p-sites that are only regulated in the synergy in at least one time point.
PIPD_specific = lapply(PIPD_regulated, function(d) {
  if(length(unique(d$Drug)) > 1) {
    d = NULL
  }
  return(d)
}) 

PIPD_specific <- PIPD_specific %>%
  discard(is.null) %>% 
  do.call("rbind", .)

# Select p-sites that are regulated in more than 2 conditions (synergy and at least one more drug).

mult_cond = lapply(PIPD_regulated, function(d) {
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

PIPD_wide <- lapply(single_tp, function(x) long_to_wide(x) %>%
                      dplyr::select(-contains("adj")) %>% 
                      delete.na(.,2) %>% 
                      pivot_longer(-PTM_collapse_key) %>% 
                      mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
                      dplyr::select(-value) %>% 
                      spread(., name, Sign) %>% 
                      rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
                                  .cols = starts_with("logFC")))


PIPD_wide <- lapply(PIPD_wide, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

PIPD_wide <- map2(PIPD_wide, names(PIPD_wide), ~ mutate(.x, Timepoint = .y))

PIPD_wide <- do.call("rbind", PIPD_wide)

synergy_specific <- c("Up_Unaffected_Down","Unaffected_Up_Down",
                      "Down_Unaffected_Up","Unaffected_Down_Up",
                      "Up_Up_Down","Down_Down_Up")


PIPD_specific_opposite <- subset(PIPD_wide, PD_PI_PIPD %in% synergy_specific)

# Qualitative -------------------------------------------------------------
# 
# higher_than_both_PIPD <- lapply(single_tp, function(x){
#   x %>%
#     group_by(PTM_collapse_key) %>% 
#     filter(n() >= 3) %>% 
#     filter(abs(logFC[Drug == "PIPD"]) > 2 * abs(logFC[Drug == "PI"]) & abs(logFC[Drug == "PIPD"]) > 2 * abs(logFC[Drug == "PD"])) 
# })
# 
# higher_than_PI <- lapply(single_tp, function(x){
#   x %>% 
#     group_by(PTM_collapse_key) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PIPD", "PI")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PIPD"]) > 2 * abs(logFC[Drug == "PI"])) 
# })
# 
# higher_than_PD <- lapply(single_tp, function(x){
#   x %>%
#     group_by(PTM_collapse_key) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(Drug %in% c("PIPD", "PD")) %>% 
#     filter(n_distinct(Drug) == 2) %>% 
#     filter(abs(logFC[Drug == "PIPD"]) > 2 * abs(logFC[Drug == "PD"])) 
# })
# 

# Multiple timepoints  ----------------------------------------------------

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
                      dplyr::select(-value) %>% 
                      spread(., name, Sign) %>% 
                      rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
                                  .cols = starts_with("logFC")))

mult_cond_named <- do.call("bind_rows",mult_cond_split)

mult_cond_named_05h <- mult_cond_named %>% 
  dplyr::select(PTM_collapse_key, contains("30m")) %>% 
  delete.na(2) %>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_2h <- mult_cond_named %>% 
  dplyr::select(PTM_collapse_key, contains("2h")) %>% 
  delete.na(2)%>% 
  mutate(across(everything(), ~replace(., is.na(.) , "Unaffected")))

mult_cond_named_8h <- mult_cond_named %>% 
  dplyr::select(PTM_collapse_key, contains("8h")) %>% 
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

changes_PIPD <- read.csv("Phosphoproteomics/Output/changes_PIPD.csv")
mult_cond_united <- left_join(mult_cond_united, changes_PIPD, by = c("PD_PI_PIPD" = "Change_type"))

mult_cond_united_named <-  mult_cond_united %>% 
  dplyr::select(-PD_PI_PIPD) %>% 
  spread(., Timepoint, Grouped_changes_I)  %>% 
  split(~PTM)

mult_cond_united_named <- lapply(mult_cond_united_named, function(x) {
  colnames(x) <- sub("_.*", "", colnames(x))
  x[paste(names(x)[2:4], collapse = '_')] <- do.call(paste, c(x[2:4], sep = '_'))
  x <- x[,c(1,5)]
  x
})

mult_cond_all_tps <- do.call("rbind",mult_cond_united_named) 

PIPD_specific_mult_tps <- mult_cond_all_tps %>% 
  filter(!grepl("Single|Same|Level", `05h_2h_8h`)) 

all_PIPD_specific <- data.frame("PTM" = c(PIPD_specific$PTM_collapse_key, PIPD_specific_opposite$PTM, PIPD_specific_mult_tps$PTM))

# Single drug driven ---------------------------------------------------

# Select p-sites that are regulated only in one or both single drugs

# SD_regulated = lapply(PIPD_long, function(d) {
#   if(("PIPD" %in% d$Drug)) {
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
#                       dplyr::select(-contains("adj")) %>% 
#                       delete.na(.,2) %>% 
#                       pivot_longer(-PTM_collapse_key) %>% 
#                       mutate(Sign = ifelse(is.na(value),"Unaffected",ifelse(value <0, "Down","Up"))) %>% 
#                       dplyr::select(-value) %>% 
#                       spread(., name, Sign) %>% 
#                       rename_with(.fn = ~ str_replace(.x, "logFC_", ""),
#                                   .cols = starts_with("logFC")))
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

  
# Single drug-driven in a single time point

PI_driven_names <- c("Unaffected_Up_Up","Unaffected_Down_Down","Down_Up_Up","Up_Down_Down")
PI_specific_names <- c("Unaffected_Up_Unaffected","Unaffected_Down_Unaffected")
PI_driven <- subset(PIPD_wide, PD_PI_PIPD %in% PI_driven_names)

PD_driven_names <- c("Up_Unaffected_Up","Down_Unaffected_Down","Up_Down_Up","Down_Up_Down")
PD_specific_names <- c("Up_Unaffected_Unaffected","Down_Unaffected_Unaffected")
PD_driven <- subset(PIPD_wide, PD_PI_PIPD %in% PD_driven_names)

# Single drug-driven in multiple time points

SD_driven_mult <- mult_cond_united
SD_driven_mult <- SD_driven_mult %>%
  mutate(Grouped_changes_I = ifelse(PD_PI_PIPD %in% PI_driven_names, "PI_driven", Grouped_changes_I)) %>%
  mutate(Grouped_changes_I = ifelse(PD_PI_PIPD %in% PD_driven_names, "PD_driven", Grouped_changes_I)) %>%
  mutate(Grouped_changes_I = ifelse(PD_PI_PIPD %in% PI_specific_names, "PI_specific", Grouped_changes_I)) %>%
  mutate(Grouped_changes_I = ifelse(PD_PI_PIPD %in% PD_specific_names, "PI_specific", Grouped_changes_I))

SD_driven_mult <-  SD_driven_mult %>%
  dplyr::select(-PD_PI_PIPD) %>%
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
  filter(!grepl("PD_|Same|Single|Synergy|specific|Level", `05h_2h_8h`)) %>% 
  filter(!PTM %in% PI_driven$PTM)

PD_driven_mult <- SD_all_tps %>% 
  filter(!grepl("PI_|Same|Single|Synergy|specific|Level", `05h_2h_8h`))%>% 
  filter(!PTM %in% PD_driven$PTM)

PI_specific_mult <- SD_all_tps %>% 
  filter(!grepl("PD_|Same|Single|Synergy|driven|Level", `05h_2h_8h`))

PD_specific_mult <- SD_all_tps %>% 
  filter(!grepl("PI_|Same|Single|Synergy|driven|Level", `05h_2h_8h`))

# Coordinated changes

Same_in_all <- SD_all_tps %>% 
  filter(!grepl("PI_|PD_|Single|Synergy|driven|Level", `05h_2h_8h`))

# Combine all results -----------------------------------------------------

PIPD_specific_all <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PIPD_specific$PTM_collapse_key | PTM_collapse_key %in% PIPD_specific_mult_tps$PTM | PTM_collapse_key %in% PIPD_specific_opposite$PTM) %>% 
  filter(Drug == "PIPD") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PIPD_specific")

PI_driven_all <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PI_driven$PTM | PTM_collapse_key %in% PI_driven_mult$PTM) %>% 
  filter(Drug == "PI") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PI_driven")

PD_driven_all <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PD_driven$PTM | PTM_collapse_key %in% PD_driven_mult$PTM) %>% 
  filter(Drug == "PD") %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "PD_driven")

Same_all <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% Same_in_all$PTM) %>% 
  group_by(Time) %>% 
  summarize(count = n()) %>% 
  mutate(Time = gsub("30m","0.5h",Time)) %>% 
  mutate(Type = "Coordinated")

PIPD_all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(PIPD_specific_all, PI_driven_all, PD_driven_all, Same_all))

cls <- c("PD_driven" = "#E69F00", "PI_driven" = "#009E73","PIPD_specific" = "#0072B2","Coordinated" = "gray" )
ggplot(PIPD_all, aes(x=Time, y=count, fill = Type)) + 
  scale_fill_manual(values= cls) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal()+ theme(text = element_text(size=18))

ggsave("Figures/ph_PIPD_specific_barplots.png",
       width = 15, 
       height = 15,
       units = "cm",
       dpi = 300)

write.csv(PIPD_all, file = "Phosphoproteomics/Output/PIPD_all_DPPs_quantified.csv", row.names = F)

PIPD_specific_all_names <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PIPD_specific$PTM_collapse_key | PTM_collapse_key %in% PIPD_specific_mult_tps$PTM | PTM_collapse_key %in% PIPD_specific_opposite$PTM) %>% 
  filter(Drug == "PIPD") 

write.csv(PIPD_specific_all_names, file = "Phosphoproteomics/Output/PIPD_specific_DPPs.csv", row.names = F)

PI_driven_all_names <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PI_driven$PTM | PTM_collapse_key %in% PI_driven_mult$PTM) %>% 
  filter(Drug != "PD")
write.csv(PI_driven_all_names, file = "Phosphoproteomics/Output/PI_driven_DPPs_PIPD", row.names = F)

PD_driven_all_names <- do.call("rbind", PIPD_long) %>% 
  filter(PTM_collapse_key %in% PD_driven$PTM | PTM_collapse_key %in% PD_driven_mult$PTM) %>% 
  filter(Drug != "PI") 

write.csv(PD_driven_all_names, file = "Phosphoproteomics/Output/PD_driven_DPPs_PIPD", row.names = F)


# Leveled-out -------------------------------------------------------------

leveled_out = lapply(PIPD_long, function(d) {
  if(("PIPD" %in% d$Drug)) {
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
              PD_change = Change[Drug == "PD"]) %>% 
    filter((PI_change == "Up" & PD_change == "Down") | 
             (PI_change == "Down" & PD_change == "Up")) %>% 
    inner_join(df, by = "PTM_collapse_key")
  return(PTM_filtered)
}

leveled_PIPD <- lapply(leveled_out, filter_PTMs)

