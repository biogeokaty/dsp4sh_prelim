# 03 - Identifying and filling missing SOC, BD, and CF measurements for DSP4SH data
# Katy Dynarski, March 2024

# This script identifies horizons with missing SOC, bulk density, or coarse fragments data (which are all necessary to calculate SOC stocks) and fills that data in to minimize data loss (i.e. having to throw out entire soil profiles because one horizon is missing a measurement). Missing values for a horizon are typically filled with the average value for a generalized horizon within that particular combination of treatment and project

# Working off of a dataframe Ekundayo put together that fills in missing cooperator SOC values with KSSL SOC - the combined values are in a column called "soc_pct_mod"

# 0 - Load libraries and import data ----
library(aqp)
library(zoo)
library(tidyverse)

coop_data <- read.csv(here("data_processed","02_coop_data_horizons_valid.csv"))

# 1 - SOC% - Identify missing data ----
no_soc <- coop_data %>%
  filter(is.na(soc_pct_mod))
# 27 samples missing SOC% (most of them are from the Texas A&M project where SOC was only measured at KSSL - it looks like neither the cooperators nor KSSL measured SOC in deeper horizons)

# 2 - SOC% - Fill missing values based on average SOC for horizon/treatment ####

# Calculate grouped means for generalized horizons grouped by treatment
# Use na.aggregate() function from zoo pkg, which will fill NA values with grouped means
# Group data by project, treatment, and generalized horizon label for na.aggregate() so that each filled value represents the average value for that horizon within a treatment for a particular project
coop_data_soc_filled <- coop_data %>%
  group_by(project, label, genhz) %>%
  mutate(soc_fill = ifelse(is.na(genhz), soc_pct_mod, na.aggregate(soc_pct_mod))) # ifelse statement fills in original values if there is no genhz (and will leave NAs if there are no SOC measurements and no genhz labels), and fills NA data with mean SOC% according to treatment and genhz if there is a genhz label

# Check for remaining NAs
coop_data_soc_filled %>%
  filter(is.na(soc_fill))
# just Texas A&M Pt 2 project - we don't have a way to fill this so just have to leave it

# 3 - Bulk Density - Identify missing data ####
no_bd <- coop_data %>%
  filter(is.na(bulk_density))
# Many samples are missing bulk density measurements - 221!

# In some projects, it looks like only the deepest soils in a pedon are missing BD measurements. In other cases, it looks like whole pedons did not have BD measurements collected.

# Fill BD values separately for each project and then join together, then join with filled SOC data to make completed df

# 4 - Bulk Density - Fill in missing data for each project ####
# Going to split data into two separate dataframes - one is the Illinois project, which doesn't have horizon labels for most pedons and will require extra math
# All the other data can have bulk density data filled via the same na.aggregate() method

no_ill_bd_fill <- coop_data %>%
  filter(project !="Illinois") %>%
  group_by(project, label, genhz) %>% # group by generalized horizons, don't need project as a grouping variable because this subset is all one project
  mutate(bd_fill = na.aggregate(bulk_density))

# Check for remaining NA values
no_bd_still <- no_ill_bd_fill %>%
  filter(is.na(bd_fill)) # just some deep soils in the UTRGV project

# Fill in data for Illinois project
ill <- coop_data %>%
  filter(project=="Illinois")

# Promote to SPC object so we can use the aqp::dice() function
ill_bd <- ill
depths(ill_bd) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill_bd) <- 'hzdesg'

plotSPC(ill_bd, color="bulk_density")
# none of the samples that had full profiles characterized had bulk density collected, plus a whole site appears to not have BD

# Need to calculate mean bulk density by depth class, since there are no generalized horizon labels for this data

# Step 1 - Add a depth class identifier column
horizons(ill_bd) <- horizons(ill_bd) %>%
  mutate(hrzdep_t_chr = as.character(hrzdep_t),
         hrzdep_b_chr = as.character(hrzdep_b)) %>%
  unite(depth_class, hrzdep_t_chr:hrzdep_b_chr, sep='-')

# Step 2 - not all the depth classes are the same - use AQP dice() to slice by depth, calculate average BD for each depth
ill_bd_dice <- aqp::dice(ill_bd, fm=0:max(ill_bd) ~ bulk_density + dsp_pedon_id + depth_class)
plotSPC(ill_bd_dice, color='bulk_density') # this looks horrible but it will work

# Step 3 - fill in the missing data with average bd for that 1 cm slice and calculate mean BD for original depth class
horizons(ill_bd_dice) <- horizons(ill_bd_dice) %>%
  group_by(hrzdep_b) %>% # group by horizon depth
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(ill_bd_dice, color='bd_fill')

# Step 4 - calculate mean BD for each depth class and add BD data back in
ill_bd_to_join <- horizons(ill_bd_dice) %>%
  group_by(dsp_pedon_id, depth_class) %>%
  summarize(bd_class_avg = mean(bd_fill))

# join calculated BD in to dataframe
ill_bd_filled <- horizons(ill_bd) %>%
  left_join(ill_bd_to_join, by=c("dsp_pedon_id", "depth_class")) %>%
  rename(bd_fill = bd_class_avg)

# make a copy to promote to SPC and plot profiles to check that everything looks right
ill_filled_spc <- ill_bd_filled 
depths(ill_filled_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill_filled_spc) <- 'hzdesg'
plotSPC(ill_filled_spc, color="bd_fill") # looks as good as it will get :)

# 4 - Bulk Density - Join project data back together ####

coop_data_bd_filled <- ill_bd_filled %>%
  select(!c('hzID', 'depth_class')) %>%
  bind_rows(no_ill_bd_fill)

# 5 - Coarse Fragments - Identify and fill missing data ####
# Coarse fragments are a little different than the other two variables we have been looking at - they are needed to calculate SOC stock, but values are often NA in dataframe not because they were overlooked, but because they weren't observed

# Assign 0 to NA values so that calculations are possible
coop_data_coarse_filled <- coop_data %>%
  mutate(coarse_frag_fill = ifelse(is.na(coarse_frag_volume), 0, coarse_frag_volume))

# 6 - Join everything together and write CSV ####
# We have built three different dataframes containing filled values for SOC, bulk density, and coarse fragments that need to be joined back together

soc_fill <- coop_data_soc_filled %>%
  ungroup() %>%
  select(dsp_sample_id, soc_fill)

cf_fill <- coop_data_coarse_filled %>%
  ungroup() %>%
  select(dsp_sample_id, coarse_frag_fill)

coop_data_filled <- coop_data_bd_filled %>%
  left_join(soc_fill, by="dsp_sample_id") %>%
  left_join(cf_fill, by="dsp_sample_id")

# check for NA values
na_check <- coop_data_filled %>%
  filter(is.na(soc_fill) | is.na(bd_fill) | is.na(coarse_frag_fill)) # only NAs left are som SOC values in Texas A&M pt 2 project, bulk density in UTRGV pits below 30 cm, and bulk density in the deep Illinois pits - can't do anything about that now :)

write_csv(coop_data_filled, here("data_processed", "03_coop_data_filled.csv"))
