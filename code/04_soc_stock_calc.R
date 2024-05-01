# Calculating SOC stocks at specific, common depth increments across many soil profiles
# Katy Dynarski, March 2024

# 0 - Load packages and import data ----
library(aqp)
library(tidyverse)

coop_data <- read.csv(here("data_processed", "03_coop_data_filled.csv"))

# 1 - Use aqp::dice() to calculate SOC stocks in each profile at 1-cm increments ----

# Promote data to SoilProfileCollection object so you can use AQP functions

coop_spc <- coop_data # I like to assign my original dataframe a new namebefore promoting to a SoilProfileCollection so I don't lose the dataframe
depths(coop_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b #make data a SoilProfileCollection, get depths from hrzdep_t and hrzdep_b columns
hzdesgnname(coop_spc) <- 'hzdesg' #get horizon names from hzdesg column

# Use dice() function to resample profile into 1 cm increments, keeping the variables needed for SOC stock calculation
# your columns might have different names, but make sure to keep the SOC%, bulk density, and coarse fragment columns
coop_dice <- aqp::dice(coop_spc, fm=0:100 ~ soc_fill + bd_fill + coarse_frag_fill)

# Extract data from the SoilProfileCollection back into a dataframe and calculate SOC stock in each 1-cm increment
# Information that is needed:
# For each horizon: soil carbon %, Db, horizon depth, and %CF
# then sum all horizons in each pedon

# SOC stock (horizon) = SOC% * BD * depth * (1-%CFrag/100)
# SOC stock (pedon) = sum(SOC horizons)

coop_soc <- horizons(coop_dice) %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)

# 2 - Calculate sums for 0-30 cm, 30-50 cm, and 50-100 cm depth increments ----

# I subsetted the data into 0-30 cm, 30-50 cm, and 50-100 cm dataframes, calculated the sum of SOC stock for each, and then re-joined the data back together
# There may be a more efficient way to do this :)

# 0-30 cm depth sum
soc_tot_30 <- coop_soc %>%
  filter(hrzdep_t < 30) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_0_30cm = sum(soc_stock_hrz))

# 30-50 cm depth sum
soc_tot_50 <- coop_soc %>%
  filter(hrzdep_t >=30 & hrzdep_t < 50) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_30_50cm = sum(soc_stock_hrz))

# 50-100cm depth sum
soc_tot_100 <- coop_soc %>%
  filter(hrzdep_t >=50 & hrzdep_t < 100) %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_50_100cm = sum(soc_stock_hrz))

# Join together into one dataframe and calculate total stock from 0-100 cm
soc_stocks <- soc_tot_30 %>%
  left_join(soc_tot_50, by="dsp_pedon_id") %>%
  left_join(soc_tot_100, by="dsp_pedon_id") %>%
  mutate(soc_stock_100cm = soc_stock_0_30cm + soc_stock_30_50cm + soc_stock_50_100cm)

# 4 - Make dataframe with SOC stocks and pedon-level information ----
# what should go into this - SOC stocks for each increment, dsp_plot_id:soil

pedon_data <- coop_data %>%
  select(dsp_pedon_id, dsp_plot_id:soil) %>%
  distinct()

pedon_soc_stock <- pedon_data %>%
  left_join(soc_stocks, by="dsp_pedon_id") 

# Save CSV
write_csv(pedon_soc_stock, here("data_processed", "04_soc_stock_pedon.csv"))

# 5 - Make dataframe with SOC stocks for originally-sampled horizons ----
# Calculate SOC stock in each horizon
coop_soc <-  coop_data %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)

# Save CSV
write_csv(coop_soc, here("data_processed", "04_soc_stock_horizon.csv"))
