# 01 - Initial data cleaning for DSP4SH data
# Katy Dynarski, October 2023

#### Housekeeping  - Libraries and Data Import ####
# Read in DB
dsp4sh4 <- dbConnect(SQLite(), here("data_raw", "dsp4sh4.db")) 

# Convert tables into dataframes
data <- lapply(setNames(nm = dbListTables(dsp4sh4)), dbReadTable, conn = dsp4sh4)
dspplotmgt <- data$dspplotmgt %>% 
  clean_names() # contains information about the sampling plots and treatments
kssl_lab <- data$kssllabmst %>% 
  clean_names() # lab data analyzed at KSSL
coop_lab <- data$cooplabmst %>% 
  clean_names() # lab data analyzed at cooperator labs, this is the data we will focus on
pedon <- data$pedon %>% 
  clean_names() # information about the each pedon sampled (location, date, etc)
layer_descr <- data$layerdescription %>% 
  clean_names() # field data for each profile sampled 
layer_desig <- data$layerdesignation %>% 
  clean_names() # designation of the layers in each soil profile sampled - this is where depths and horizon names come from
project_design <- data$projectdesign %>%
  clean_names() # information about the overall project design for individual cooperator projects

# Check data tables for duplication and errors ####
# The goal is to be able to join the lab data to the field description data in order to eventually calculate things like SOC stocks in a soil profile (which requires both lab data and field data). The depths of the soil horizons are stored in layer_desig, so that will also need to be joined in. And, to analyze data we want information about the plots and treatments.

# Start by checking tables for any duplicate values or anything that's "off" and might effect data joining and downstream analysis

# Lab data check
coop_lab %>%
  add_count(dsp_sample_id) %>%
  select(cooplabmst_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct() # No duplicates

# Layer designation table check
layer_desig %>%
  add_count(dsp_sample_id) %>%
  select(laydesg_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct()
# 182 duplicates here, pick out the sample IDs so we can filter for the duplicated sample IDs

dupe_ids <- layer_desig %>%
  add_count(dsp_sample_id) %>%
  filter(n>1) %>%
  select(dsp_sample_id)

# Filter for duplicated sample IDs
dupe_hzdesig <- layer_desig %>%
  filter(dsp_sample_id %in% dupe_ids$dsp_sample_id) %>%
  arrange(dsp_sample_id)
write.csv(dupe_hzdesig, here("data", "dupe_hzdesig.csv"), row.names=FALSE)
# issue appears to be that for some sample IDs (i.e. horizon within a pedon), there are multiple horizon designations for the same depth, which is creating extra rows 
# without fully knowing how horizons data was combined, I would wager a guess that this happened if horizon designations differed between the layer description and the layer lab measurement data sheets
# looking at original Excel sheets - yes, I do think that this is the problem
# differences between lab and field sheets in: horizon name, horizon number, or horizon depths can all cause duplication with the join

# Field layer description data check
layer_descr %>%
  add_count(dsp_sample_id) %>%
  select(layerdes_id, dsp_sample_id, dsp_pedon_id, n) %>%
  filter(n>1) %>%
  distinct() # Also check

# Also want to check - based on field data, how many lab samples would I expect to see? do all the row numbers seem correct?
nrow(coop_lab) # 1685
nrow(layer_descr) # 1374
# there are 311 more lab rows than field rows

# everything in the lab row that's not in the field row
field_samps <- layer_descr %>%
  select(dsp_sample_id)
lab_no_field <- coop_lab %>%
  filter(!dsp_sample_id %in% field_samps$dsp_sample_id)

# Talked to Ekundayo - for some projects, they only did a full field characterization for the central pedon, and then collected soil samples from satellite pedons for lab analysis, so these pedons will have lab data but no field data. Maybe can use generalized horizon labels to map data from the central pedon to the other pedons? Cross that bridge when we get to it.

# Make a table of horizon designations that doesn't have duplicates - from lab data only ####
lab_orig <- read.csv(here("data", "DSP4SH_CL_coop_lab_phase1-2v2.csv")) %>%
  clean_names()
# verify that data are the same size as data from DB
nrow(lab_orig)
nrow(coop_lab)
# yes, good to go

# pull out a new horizon designation dataframe containing same columns as one in DB - but this one will have no duplicates :)
names(layer_desig)
layer_desig_lab <- lab_orig %>%
  select(dsp_sample_id, dsp_pedon_id, layer_no, hzdesg, hrzdep_t, hrzdep_b)

# Make dataframe joining lab data to field data ####
# Join dsp_pedon_id to dsp_plot_id so that site-level data can be attached to pedon-level data
pedon_plot <- pedon %>%
  select(dsp_pedon_id, dsp_plot_id)

# Initially thought about doing a full join, but I think the better better approach is to do a left join starting with the lab data frame so that I am only pulling field data for which lab data exists

layer_join <- coop_lab %>%
  left_join(layer_desig_lab, by=c("dsp_sample_id", "dsp_pedon_id", "layer_no")) %>%
  left_join(layer_descr, by=c("dsp_sample_id", "dsp_pedon_id", "dsp_pedon"), 
            suffix=c("_lab", "_field")) %>%
  left_join(pedon_plot, by="dsp_pedon_id") %>%
  left_join(dspplotmgt, by="dsp_plot_id")

# Checking that all data is assigned to a project ####
# Extract names of projects
projects <- layer_join %>%
  distinct(project)
# there is an NA project name in here so check where that is coming from
proj_na <- layer_join %>%
  filter(is.na(project))
# "NTR -  -1, 2  and 3" pedons - part of Texas A&M Pt 1 project on Amarillo soil series, this should be no-till rye. I suspect problem stemmed fro the weird naming for this plot (no number?) But there appears to be data for it in the original Excel files

# Take a closer look

tx_ntr <- layer_join %>%
  filter((str_detect(dsp_pedon_id, "^NTR")))

# based on comparison with original Excel files and project table, it appears that the "3" got dropped from the "NTR - 3" treatment, so now it isn't joining properly

# how to fix - should fix in original tables and then re-join
# need to detect anytime the string "NTR - -" exists and replace with "NTR - 3 -" and leave the rest of the string. Do this in dsp_pedon_id and dsp_sample_id columns
# test in copy table
tx_ntr2 <- tx_ntr %>%
  mutate(dsp_sample_id = str_replace_all(dsp_sample_id, "NTR -  -", "NTR - 3 -"),
         dsp_pedon_id = str_replace_all(dsp_pedon_id, "NTR -  -", "NTR - 3 -"))
# that works

# replace in dataframe
coop_lab2 <- coop_lab %>%
  mutate(dsp_sample_id = str_replace_all(dsp_sample_id, "NTR -  -", "NTR - 3 -"),
         dsp_pedon_id = str_replace_all(dsp_pedon_id, "NTR -  -", "NTR - 3 -"))
layer_descr2 <- layer_descr %>%
  mutate(dsp_sample_id = str_replace_all(dsp_sample_id, "NTR -  -", "NTR - 3 -"),
         dsp_pedon_id = str_replace_all(dsp_pedon_id, "NTR -  -", "NTR - 3 -"))
layer_desig_lab2 <- layer_desig_lab %>%
  mutate(dsp_sample_id = str_replace_all(dsp_sample_id, "NTR -  -", "NTR - 3 -"),
         dsp_pedon_id = str_replace_all(dsp_pedon_id, "NTR -  -", "NTR - 3 -"))

# Make final join with corrected data ####
layer_join2 <- coop_lab2 %>%
  left_join(layer_desig_lab2, by=c("dsp_sample_id", "dsp_pedon_id", "layer_no")) %>%
  left_join(layer_descr2, by=c("dsp_sample_id", "dsp_pedon_id", "dsp_pedon"), 
            suffix=c("_lab", "_field")) %>%
  left_join(pedon_plot, by="dsp_pedon_id") %>%
  left_join(dspplotmgt, by="dsp_plot_id")

# Check that it worked and all data has assigned project
projects2 <- layer_join2 %>%
  distinct(project)
# yay!! no NA values!

# it looks like there are also a number of pedons that don't have associated plot IDs
# for now - it's hard to backfill this information, just leave it and note

# Write CSV ####
write_csv(layer_join2, here("data_processed", "coop_data_full.csv"), na="NA")
