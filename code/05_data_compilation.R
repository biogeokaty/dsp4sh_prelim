# Data wrangling

# This script contains code to make the different dataframes that will be used in downstream analyses

# 0 - Import data ----
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# 1 - Dataframe of project information ----
# Need to pull pedon lat and long
project_dat <- soc_pedon %>%
  select(dsp_pedon_id, label_full, label, lu, till, trt, explanation, project, soil, site, pedon_x, pedon_y)
write_csv(project_dat, here("data_processed", "05_project_data.csv"))

# pull out coordinates to extract WorldClim data
coords <- project_dat %>%
  select(project, soil, dsp_pedon_id, pedon_x, pedon_y)

just_coords <- coords %>%
  select(pedon_x, pedon_y) %>%
  rename(x = pedon_x,
         y=pedon_y) %>%
  relocate(x, .before=y)

climate <- worldclim_global(res=2.5, path=here("data_raw"), var="bio") # download WorldClim data using geodata package, will download the first time you run but is faster after that
# the bioclimatic variables that I want (for now) are BIO1 (MAT) and BIO12 (MAP)

# Make RasterStack of just MAT and MAP
mat_map <- climate[[c(1,12)]] 
names(mat_map) <- c("mat", "map")

site_dat <- terra::extract(mat_map, just_coords)
site_clim <- cbind(coords, site_dat)

# Group sites into warm and wet (MAP > 900 and MAT>15), warm and dry (MAT > 15, MAP < 900), cool and wet (MAT<15, MAP>900), and cool and dry (MAT<15, MAP<900)
site_clim2 <- site_clim %>%
  mutate(climate = ifelse(mat>15, ifelse(map>900, "warm_wet", "warm_dry"), ifelse(map>900, "cool_wet", "cool_dry"))) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup") %>% # fill in missing climate grouping for rows that don't have lat/long
  ungroup()

# Save CSV
write_csv(site_clim2, here("data_processed", "05_site_clim.csv"))

# 2 - Attach climate data to pedon and horizon data ----
# Attach climate to pedon data
pedon_clim <- soc_pedon %>%
  left_join(site_clim2, by=c("dsp_pedon_id", "project", "soil")) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup")  # fill in missing climate grouping for rows that don't have lat/long

# Write CSV
write_csv(pedon_clim_inf, here("data_processed", "05_soc_pedon_clim.csv"))

# Attach climate data to horizon data
horizon_clim <- soc_horizon %>%
  left_join(site_clim2, by=c("dsp_pedon_id", "project", "soil")) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup") # fill in missing climate grouping for rows that don't have lat/long

write_csv(horizon_clim_inf, here("data_processed", "05_soc_horizon_clim.csv"))

# 3 - Make filtered pedon and horizon dataframes that exclude soils that aren't consistent between treatments ----

# Check data to see if soil series are consistent between treatments: 
soc_pedon %>% distinct(soil, label, project)
# Moist projects have treatments replicated across soil series, except for University of Minnesota project (Marquis, Readlyn, and Kenyon series)

# For looking at soc stocks to 100 cm -
soils_exclude <- c("Kenyon", "Marquis", "Readlyn")

soc_pedon_filt <- pedon_clim %>%
  filter(!soil %in% soils_exclude)
# What this data should be used for: mixed linear models to look at influence of management on SOC stocks, accounting for variability due to soil series

write_csv(soc_pedon_filt, here("data_processed", "05_soc_pedon_filt.csv"))

# Make a dataframe of the horizon data excluding the same soils
soc_horizon_filt <- horizon_clim %>%
  filter(!soil %in% soils_exclude)
write_csv(soc_horizon_filt, here("data_processed", "05_soc_horizon_filt.csv"))
# What this data should be used for: depth plots and mixed linear models to look at influence of management on SOC stocks, accounting for variability due to soil series

# 4 - Dataframe of indicator data in surface horizons only ----

# Make version that is filtered (removing projects where soils are not consistent between treatments), and an unfiltered version
surf_all <- horizon_clim %>%
  filter(hrzdep_b == "5" | hrzdep_b=="10") %>%
  group_by(dsp_pedon_id)

write_csv(surf_all, here("data_processed", "05_surface_horizons.csv"))
# What this data should be used for: analyzing sensitivity of indicators in surface horizons, in cases where it doesn't matter that soil series is not consistent between treatments in some projects

surf_filt <- soc_horizon_filt %>%
  filter(hrzdep_b == "5") 
write_csv(surf, here("data_processed", "05_surface_horizons_filt.csv"))
# What this data should be used for: analyzing sensitivity of indicators in surface horizons when soil series needs to be held constant between projects

# 5 - Dataframe with both SOC stock data (pedon) and surface indicators for matching pedons ----
# want to use this to calculate response ratios :)
surf_sub <- surf_all %>%
  select(dsp_pedon_id, soc_pct, 
         bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace)

meta_df <- pedon_clim_inf %>%
  left_join(surf_sub, by=c("dsp_pedon_id"))

write_csv(meta_df, here("data_processed", "05_meta_df.csv"))