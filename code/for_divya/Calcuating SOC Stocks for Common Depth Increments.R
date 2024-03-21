# Calculating SOC stocks at specific, common depth increments across many soil profiles
# Katy Dynarski, March 2024

# Load packages
library(aqp)
library(tidyverse)

# Import data
coop_data <- read.csv(here("data_processed", "03_coop_data_filled.csv"))

# For this example, I'm subsetting the full dataset to just look at just Kansas State data to make it easier to plot profiles and check math
sub <- coop_data %>%
  filter(project=="KansasState")

# Promote data to SoilProfileCollection object so you can use AQP functions
sub_spc <- sub # I like to assign my original dataframe a new namebefore promoting to a SoilProfileCollection so I don't lose the dataframe
depths(sub_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b #make data a SoilProfileCollection, get depths from hrzdep_t and hrzdep_b columns
hzdesgnname(sub_spc) <- 'hzdesg' #get horizon names from hzdesg column

# Use dice() function to resample profile into 1 cm increments, keeping the variables needed for SOC stock calculation
# your columns might have different names, but make sure to keep the SOC%, bulk density, and coarse fragment columns
sub_dice <- aqp::dice(sub_spc, fm=0:100 ~ soc_fill + bd_fill + coarse_frag_fill)
# you can plot the soil profiles just to check that everything looks right
plotSPC(sub_dice, color="soc_fill") 

# Extract data from the SoilProfileCollection back into a dataframe and calculate SOC stock in each 1-cm increment
# Information that is needed:
# For each horizon: soil carbon %, Db, horizon depth, and %CF
# then sum all horizons in each pedon

# SOC stock (horizon) = SOC% * BD * depth * (1-%CFrag/100)
# SOC stock (pedon) = sum(SOC horizons)

sub_soc <- horizons(sub_dice) %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)

# Calculate sums for 0-30 cm, 30-50 cm, and 50-100 cm depth increments
# I subsetted the data into 0-30 cm, 30-50 cm, and 50-100 cm dataframes, calculated the sum of SOC stock for each, and then re-joined the data back together
# There may be a more efficient way to do this :)

# 0-30 cm depth sum
soc_tot_30 <- sub_soc %>%
  filter(hrzdep_t < 30) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_0_30cm = sum(soc_stock_hrz))

# 30-50 cm depth sum
soc_tot_50 <- sub_soc %>%
  filter(hrzdep_t >=30 & hrzdep_t < 50) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_30_50cm = sum(soc_stock_hrz))

# 50-100cm depth sum
soc_tot_100 <- sub_soc %>%
  filter(hrzdep_t >=50 & hrzdep_t < 100) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_50_100cm = sum(soc_stock_hrz))

# Join together into one dataframe and calculate total stock from 0-100 cm
soc_stocks <- soc_tot_30 %>%
  left_join(soc_tot_50, by="dsp_pedon_id") %>%
  left_join(soc_tot_100, by="dsp_pedon_id") %>%
  mutate(soc_stock_100cm = soc_stock_0_30cm + soc_stock_30_50cm + soc_stock_50_100cm)

# Also calculate total 0-100 cm stock from original (no dicing) horizon depths as a check
sub_tot <- sub %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_100cm_alt = sum(soc_stock_hrz))

# Join this together to be able to check that math worked - "soc_stock_100cm" column and "soc_stock_100cm_alt" column should have the same values
soc_stocks %>% left_join(sub_tot, by="dsp_pedon_id")

