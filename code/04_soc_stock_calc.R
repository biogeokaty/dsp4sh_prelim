# 04 - SOC Stock Calculation
# Katy Dynarski, October 2023

# This script takes processed cooperator data that has been checked for NA SOC, bulk density, and coarse fragment values (see script 03_missing_msmts for methods) and calculates SOC stocks for each soil profile

# 0 - Import data ####

coop_data <- read.csv(here("data_processed", "03_coop_data_filled.csv"))

# 1 - Calculate SOC stock in each horizon ####

# Information that is needed:
# For each horizon: soil carbon %, Db, horizon depth, and %CF
# then sum all horizons in each pedon

# SOC stock (horizon) = SOC% * BD * depth * (1-%CFrag/100)
# SOC stock (pedon) = sum(SOC horizons)

# Calculate SOC stock in each horizon
coop_soc <-  coop_data %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
                    cf_mult = 1 - (coarse_frag_fill/100)) %>%
  mutate(soc_stock_hrz = soc_fill * bd_fill * hrzdepth * cf_mult)

# Save CSV
write_csv(coop_soc, here("data_processed", "04_soc_stock_horizon.csv"))

# promote to SPC
coop_spc <- coop_soc
depths(coop_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(coop_spc) <- 'hzdesg'

# 2 - Calculate SOC stocks to 100 cm deep ####

# Identify pedons < 100 cm deep
pedons_less_100 <- coop_soc %>% 
  group_by(dsp_pedon_id) %>%
  summarize(depth = max(hrzdep_b)) %>%
  filter(depth < 100)

# Truncate pedons to 100 cm AND exclude pedons <100 cm 
coop_trunc100 <- trunc(coop_spc, 0, 100)
coop_100cm <- horizons(coop_trunc100) %>%
  filter(!dsp_pedon_id %in% pedons_less_100$dsp_pedon_id)

soc_stock_100cm <- coop_100cm %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_100cm = sum(soc_stock_hrz))
# get NA values for all of the profiles with no BD measurements below 30 cm

# Look at density
ggplot(soc_stock_100cm, aes(x=soc_stock_100cm)) +
  geom_density()

# 3 - Calculate SOC stocks to 30 cm deep ####

# Truncate to 30 cm and calculate those values (this just lets us capture some of the profiles that don't have bulk density below 30 cm)
coop_trunc30 <- trunc(coop_spc, 0, 30)

# Identify any pedons < 30 cm
pedons_less_30 <- coop_soc %>%
  group_by(dsp_pedon_id) %>%
  summarize(depth = max(hrzdep_b)) %>%
  filter(depth < 30)

coop_30cm <- horizons(coop_trunc30) %>%
  filter(!dsp_pedon_id %in% pedons_less_30$dsp_pedon_id)

soc_stock_30cm <- coop_30cm %>%
  group_by(dsp_pedon_id) %>%
  summarize(soc_stock_30cm = sum(soc_stock_hrz))

# 4 - Make dataframe with SOC stocks and pedon-level information ####
# what should go into this - SOC stocks to 30 and 100 cm, dsp_plot_id:soil

pedon_data <- coop_data %>%
  select(dsp_pedon_id, dsp_plot_id:soil) %>%
  distinct()

pedon_soc_stock <- soc_stock_30cm %>%
  left_join(soc_stock_100cm, by="dsp_pedon_id") %>%
  left_join(pedon_data, by="dsp_pedon_id")

# Write CSV
write_csv(pedon_soc_stock, here("data_processed", "04_soc_stock_pedon.csv"))
