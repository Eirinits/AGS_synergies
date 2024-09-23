library(decoupleR)
library(tidyverse)
library(here)
library(OmnipathR)
library(dorothea)
library(decoupleR)
library(workflowr)
library(rmarkdown)
library(org.Hs.eg.db)

uniprot_kinases <- OmnipathR::import_omnipath_annotations(resources = "UniProt_keyword") %>%
  dplyr::filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
  distinct() %>%
  pull(genesymbol) %>%
  unique()

omnipath_ptm <- OmnipathR::get_signed_ptms() %>%
  dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
  dplyr::filter(!(stringr::str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
  dplyr::mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
  dplyr::transmute(p_site, enzyme_genesymbol, mor) %>%
  dplyr::filter(enzyme_genesymbol %in% uniprot_kinases)

omnipath_ptm$likelihood <- 1

#we remove ambiguous modes of regulations
omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
omnipath_ptm <- omnipath_ptm[,-5]

names(omnipath_ptm)[c(1,2)] <- c("target","tf")

omnipath_ptm_flt <- omnipath_ptm %>% 
  separate(target, c("target","site"), sep = "_") %>% 
  dplyr::select(-site) %>% 
  distinct()

res_decoupler_df <- bind_rows(res_decoupler, .id = "column_label")

kin_activity <- run_wmean(
  mat = as.matrix(phospho_differential_analysis), 
  network = omnipath_ptm, 
  .source = "tf",
  times = 1000
)

kin_activity <- kin_activity[kin_activity$statistic == "norm_wmean",c(2,4)] %>%
  tibble::column_to_rownames(var = "source")

dictionary <- read.csv("~/Documents/Git/cell_fate_modules/Module_topology/20240104_cell_fate_dictionary.csv")
dic <- dictionary %>% 
  separate_rows(HGNC_symbol, sep = ",") %>% 
  filter(!grepl("rna", HGNC_symbol)) %>% 
  mutate(across(where(is.character), str_trim)) %>% 
  distinct(HGNC_symbol)

clipr::write_clip(dic$HGNC_symbol)
