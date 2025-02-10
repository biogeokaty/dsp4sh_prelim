# Preliminary data exploration of DSP4SH data v2
# Katy Dynarski, October 2023

#### Housekeeping  - Libraries and Data Import ####
library(here)
library(tidyverse)
library(DBI)
library(RSQLite)
library(janitor)

# Read in DB
dsp4sh4 <- dbConnect(SQLite(), here("data", "dsp4sh4.db")) 

# Convert DB tables into dataframes
data <- lapply(setNames(nm = dbListTables(dsp4sh4)), dbReadTable, conn = dsp4sh4)
dspplotmgt <- data$dspplotmgt %>% 
  clean_names()
kssl_lab <- data$kssllabmst %>% 
  clean_names()
coop_lab <- data$cooplabmst %>% 
  clean_names()
pedon <- data$pedon %>% 
  clean_names()
layer_descr <- data$layerdescription %>% 
  clean_names()
layer_desig <- data$layerdesignation %>% 
  clean_names()

# Check data for duplicates
# Cooperator lab data
coop_lab %>%
  add_count(dsp_sample_id) %>%
  select(cooplabmst_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct() # no duplicates

# Cooperator field layer descriptions
layer_descr %>%
  add_count(dsp_sample_id) %>%
  select(layerdes_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct() # no duplicates

# Layer designation
layer_desig %>%
  add_count(dsp_sample_id) %>%
  select(laydesg_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct()
# 182 duplicates here

# pick out the sample IDs so we can filter for the duplicated sample IDs
dupe_ids <- layer_desig %>%
  add_count(dsp_sample_id) %>%
  filter(n>1) %>%
  select(dsp_sample_id)

# filter for duplicated sample IDs
dupe_hzdesig <- layer_desig %>%
  filter(dsp_sample_id %in% dupe_ids$dsp_sample_id) %>%
  arrange(dsp_sample_id)
# issue appears to be that for some sample IDs (i.e. horizon within a pedon), there are multiple horizon designations for the same depth, which is creating extra rows 

# It also looks like there are more samples with lab measurements than field observations
nrow(coop_lab) # 1685
nrow(layer_descr) # 1374
# there are 311 more lab rows than field rows

# Filter for samples with lab data that have no field data
field_samps <- layer_descr %>%
  select(dsp_sample_id)
lab_no_field <- coop_lab %>%
  filter(!dsp_sample_id %in% field_samps$dsp_sample_id)