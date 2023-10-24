# 03 - SOC Stock Calculation
# Katy Dynarski, October 2023

# Calculate SOC Stocks for each pedon ####

# Information that is needed:
# For each horizon: soil carbon %, Db, horizon depth, and %CF
# then sum all horizons in each pedon
# CF should be in layer description for each layer
# Also need to calculate horizon depth

# SOC stock (pedon) = sum(SOC% * BD * depth * (1-%CFrag/100))

# Calculate SOC stock of each horizon, attach it to profile data SPC, then use profileApply() to calculate SOC stock at each site
# or could use group_by() and summarize()

soc_stock1 <- ksh %>%
  select(dsp_sample_id, dsp_pedon_id, bulk_density, soc_pct, hrzdep_t, hrzdep_b, coarse_frag_volume) %>%
  mutate(hrzdepth = hrzdep_b - hrzdep_t,
         cf_mult = 1 - (coarse_frag_volume/100)) %>%
  mutate(soc_stock = soc_pct * bulk_density * hrzdepth * cf_mult)

soc_stock2 <- soc_stock1 %>%
  select(dsp_sample_id, dsp_pedon_id, hrzdepth, soc_stock)