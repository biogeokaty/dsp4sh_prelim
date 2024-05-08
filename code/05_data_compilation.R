# Data wrangling

# This script contains code to make the different dataframes that will be used in downstream analyses

# 0 - Import data ----
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))
prism <- read.csv(here("data_processed", "dsp4sh_soc_prism3.csv"))
dsp4sh6 <- dbConnect(SQLite(), here("data_raw", "dsp4sh6.db")) 

# 1 - Dataframe of project information ----
# Clean up PRISM data for joining
prism_clean <- prism %>%
  select(DSP_Pedon_ID, pedon_x, pedon_y, PRISM_ppt_30yr_normal_800mM4_annual_bil, PRISM_tmean_30yr_normal_800mM5_annual_bil) %>%
  clean_names() %>%
  distinct() %>%
  rename(map = prism_ppt_30yr_normal_800m_m4_annual_bil,
         mat = prism_tmean_30yr_normal_800m_m5_annual_bil)

# Make dataframe of project overview with climate information from PRISM
project_dat <- soc_pedon %>%
  select(dsp_pedon_id, label_full, label, lu, till, trt, project, soil, site, pedon_x, pedon_y) %>%
  left_join(prism_clean, by=c("dsp_pedon_id", "pedon_x", "pedon_y")) %>%
  mutate(climate = ifelse(mat>15, ifelse(map>900, "warm_wet", "warm_dry"), ifelse(map>900, "cool_wet", "cool_dry"))) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup") %>% # fill in missing climate grouping for rows that don't have lat/long
  ungroup()
write_csv(project_dat, here("data_processed", "05_project_data.csv"))

# Make reduced dataframe of site climate information for later joining
site_clim <- project_dat %>%
  select(project, soil, dsp_pedon_id, pedon_x, pedon_y, mat, map, climate)

# 2 - Pull in clay% data from KSSL ----
# Convert tables into dataframes
db <- lapply(setNames(nm = dbListTables(dsp4sh6)), dbReadTable, conn = dsp4sh6)

# pull out KSSL data
kssl <- db$kssllabmst %>% 
  clean_names()

clay <- kssl %>%
  select(natural_key, clay_tot_psa) %>%
  rename(kssl_labsampnum = natural_key)

# 3 - Attach climate data to pedon and horizon data ----
# Attach climate to pedon data
pedon_clim <- soc_pedon %>%
  left_join(select(site_clim, dsp_pedon_id, project, soil, mat, map, climate), by=c("dsp_pedon_id", "project", "soil")) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup")  # fill in missing climate grouping for rows that don't have lat/long

# Write CSV
write_csv(pedon_clim, here("data_processed", "05_soc_pedon_clim.csv"))

# Attach climate data and clay% to horizon data
horizon_clim <- soc_horizon %>%
  left_join(select(site_clim, dsp_pedon_id, project, soil, mat, map, climate), by=c("dsp_pedon_id", "project", "soil")) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup") %>% # fill in missing climate grouping for rows that don't have lat/long
  left_join(clay, by="kssl_labsampnum") # join in KSSL clay % particle size analysis data

write_csv(horizon_clim, here("data_processed", "05_soc_horizon_clim.csv"))

# 4 - Make filtered pedon and horizon dataframes that exclude soils that aren't consistent between treatments ----

# Check data to see if soil series are consistent between treatments: 
soc_pedon %>% distinct(soil, label, project)
# Moist projects have treatments replicated across soil series, except for University of Minnesota project (Marquis, Readlyn, and Kenyon series)

# For looking at soc stocks to 100 cm -
soils_exclude <- c("Kenyon", "Marquis", "Readlyn")

soc_pedon_filt <- pedon_clim %>%
  filter(!soil %in% soils_exclude)
# What this data should be used for: mixed linear models to look at influence of management on SOC stocks, accounting for variability due to soil series
# Don't apply filtering/exclusion to anything but mixed linear models

write_csv(soc_pedon_filt, here("data_processed", "05_soc_pedon_filt.csv"))

# Make a dataframe of the horizon data excluding the same soils
soc_horizon_filt <- horizon_clim %>%
  filter(!soil %in% soils_exclude)
write_csv(soc_horizon_filt, here("data_processed", "05_soc_horizon_filt.csv"))
# What this data should be used for: depth plots and mixed linear models to look at influence of management on SOC stocks, accounting for variability due to soil series

# 6 - Dataframe of indicator data in surface horizons only ----
# Want average values across 0-10 cm (i.e. average 0-5 and 5-10 for each pedon) for all indicators
surf_all <- horizon_clim %>%
  filter(hrzdep_b == "5" | hrzdep_b=="10") %>%
  group_by(dsp_pedon_id) %>%
  summarize(across(c(soc_pct,bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace, clay_tot_psa), mean)) %>%
  left_join(select(site_clim, dsp_pedon_id, mat, map, climate), by="dsp_pedon_id") %>%
  left_join(select(project_dat, dsp_pedon_id:site), by="dsp_pedon_id")
  
write_csv(surf_all, here("data_processed", "05_surface_horizons.csv"))
# What this data should be used for: analyzing sensitivity of indicators in surface horizons

# 5 - Dataframe for meta-analysis with both SOC stock data (pedon) and surface indicators for matching pedons ----
# want to use this to calculate response ratios :)
surf_sub <- surf_all %>%
  select(dsp_pedon_id, soc_pct, 
         bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace, clay_tot_psa)

meta_df <- pedon_clim %>%
  left_join(surf_sub, by=c("dsp_pedon_id"))

write_csv(meta_df, here("data_processed", "05_meta_df.csv"))