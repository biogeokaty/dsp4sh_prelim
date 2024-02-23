# 07 - Big Silly Models
# Katy Dynarski, December 2023/January 2024

# 0 - Import data
soc_pedon100 <- read.csv(here("data_processed", "05_soc_pedon_100cm.csv"))
soc_horizon_filt <- read.csv(here("data_processed", "05_soc_horizon_filtered.csv"))
soc_horizon_all <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))
surf <- read.csv(here("data_processed", "06_surface_horizons.csv"))

# 1 - Extract climate data from WorldClim ----
# Need to pull pedon lat and long
project <- read.csv(here("data_raw", "DSP4SH_pedon_mgtv4.csv")) %>%
  clean_names() %>%
  mutate(pedon_x2 = ifelse(project=="Illinois" | project=="UConn" | project=="UTRGV", pedon_y, pedon_x), # Illinois, UConn, and UTRGV all have coordinates reversed (y coords are in x column)
         pedon_y2 = ifelse(project=="Illinois" | project=="UConn" | project=="UTRGV", pedon_x, pedon_y)) %>%
  select(-pedon_x, -pedon_y) %>%
  rename(lat = pedon_y2,
         long = pedon_x2) %>%
  relocate(lat, .before=long) %>%
  select(dsp_pedon_id, label_full, label, lu, till, trt, explanation, state, project, soil, site, lat, long)
write_csv(project, here("data_processed", "07_project_data.csv"))

coords <- project %>%
  select(project, soil, dsp_pedon_id, lat, long)

just_coords <- coords %>%
  select(lat, long) %>%
  rename(x = long,
         y=lat) %>%
  relocate(x, .before=y)

climate <- getData("worldclim", var="bio", res=2.5)
names(climate) # check naming convention for bioclimatic variables
# the bioclimatic variables that I want (for now) are BIO1 (MAT) and BIO12 (MAP)
plot(climate$bio12, main="MAP") # plot MAP as a check

# Make RasterStack of just MAT and MAP
mat_map <- climate[[c(1,12)]] 
names(mat_map) <- c("mat_10", "map")

site_dat <- raster::extract(mat_map, just_coords)
site_clim <- cbind(coords, site_dat) %>%
  mutate(mat = mat_10/10) %>% # MAT is originally formatted as MAT*10, divide by 10 to get correct number
  select(-mat_10)

soc_pedon_clim <- soc_pedon100 %>%
  left_join(site_clim, by=c("dsp_pedon_id", "project", "soil"))
write_csv(soc_pedon_clim, here("data_processed", "07_soc_pedon_clim.csv"))

# 2 - Calculate average MAP and MAT for each soil series in the project, plot alternative version of stock boxplot for Pedometrics talk ----
avg_clim <- soc_pedon_clim %>%
  group_by(soil) %>%
  summarize(across(map:mat, ~ mean(.x, na.rm = TRUE))) %>%
  rename(avg_map = map,
         avg_mat = mat)

soc_pedon_clim2 <- soc_pedon_clim %>%
  left_join(avg_clim, by="soil")
  

ggplot(soc_pedon_clim2, aes(x=fct_reorder(soil, avg_map), y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  labs(x="Soil Series", y="SOC stock to 100 cm depth (Mg/ha)") +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy()
ggsave(here("figs", "soc_stock_boxplot_soil_mgmt_ordered.png"), width=10, height=7, units="in", dpi=400)

# 2 - Build dataframe for modeling ----
# Want from pedon df: SOC stock to 100 cm, soil, label, lu, till, MAT, MAP
# want from surface df: all the indicators
# want but don't (yet) have: other state factor data - what should this be??

surf_keep <- surf %>%
  select(dsp_pedon_id, soc_stock_hrz, soc_pct, bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace) # same indicators as before

soc_df <- soc_pedon_clim %>%
  left_join(surf_keep, by="dsp_pedon_id") %>%
  mutate(across(where(is.character), as.factor)) %>%
  filter(!is.na(soc_stock_100cm)) # Filter out rows with no SOC stock data - they mess with the random forest model later on


# 3 - Ctree and cforest ----
# Try a ctree - with land use variables and top 5 soil health metrics
soc_tree <- ctree(soc_stock_100cm ~ soc_pct + ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                         label + lu + till + soil + mat + map, data=soc_df)
plot(soc_tree)
# Soil series is the most important determinant, followed by SOC%

# try a version that doesn't include soil series - want to look at the underlying factors integrated by soil series
soc_tree2 <- ctree(soc_stock_100cm ~ soc_pct + ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                    label + lu + till + mat + map, data=soc_df)
plot(soc_tree2)
# SOC% is the most important thing...well duh. SOC stocks are driven by how much soil carbon is there. Should take that out of the tree...

soc_tree3 <- ctree(soc_stock_100cm ~ ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                     label + lu + till + mat + map, data=soc_df)
plot(soc_tree3)
# Now we're getting somewhere...First split is MAP, followed by ACE, then bglucosaminidase, then treatment (!)

# Try as a forest
soc_forest <- cforest(soc_stock_100cm ~ ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                        label + lu + till + mat + map, data=soc_df, ntree=10000)
soc_vi <- vi(soc_forest)
soc_vip <- vip(soc_vi, geom="col")
soc_vip

# Most important variables for predicting SOC stocks to 100 cm are: ACE, MAT, MAP, land use, and KSSL WSA

# Treatment (BAU vs SHM vs Ref) was the second least useful variable - but it is also the most general land-use variable here so that makes sense

# The question with indicators - how much variance is driven by soil series and how much variance is driven by treatment (and does this matter if we can still pick up the effect of treatment) 

# Is it useful to break treatment down into land use and tillage? 

# 4 - Can we predict SOC stocks via the other metrics ----
stock_lm1 <- lm(soc_stock_100cm ~ yoder_agg_stab_mwd + soil_respiration + bglucosidase + bglucosaminidase + 
                 alkaline_phosphatase + acid_phosphatase + phosphodiesterase + arylsulfatase + pox_c + ace, data=soc_df)
summary(stock_lm1) # This is not a very good model. R2 0.19, only ace and bglucosidase are significant
# It's really hard to model this with so many variables because 
stock_sub1 <- regsubsets(soc_stock_100cm ~ soc_pct + yoder_agg_stab_mwd + soil_respiration + bglucosidase + bglucosaminidase + 
                          alkaline_phosphatase + acid_phosphatase + phosphodiesterase + arylsulfatase + pox_c + ace, data=soc_df)
summary(stock_sub1)

stock_lm2 <- lm(ssoc_stock_100cm ~ soc_pct + yoder_agg_stab_mwd + soil_respiration + bglucosidase + bglucosaminidase + 
                  alkaline_phosphatase + acid_phosphatase + phosphodiesterase + arylsulfatase + pox_c + ace, data=soc_df)
summary(stock_lm2)
stock_sub2 <- regsubsets(soc_stock_100cm ~ soc_pct + yoder_agg_stab_mwd + soil_respiration + bglucosidase + bglucosaminidase + 
                           alkaline_phosphatase + acid_phosphatase + phosphodiesterase + arylsulfatase + pox_c + ace, data=soc_df)
summary(stock_sub2)

# There's some really interesting stuff in here - 0-5cm SOC % is by far the best predictor of overall SOC stocks. This makes sense (soils with the most SOC throughout the profile should also have the most SOC in the surface horizon).

# In models with SOC included, metrics like POX-C and ACE don't contribute much extra information (they are basically a proxy for carbon) and get kicked out of the model. Enzyme activity is a much more important variable. But in models without SOC included, ACE is the most important variable (probably because it approximates SOC) followed by bglucosidase

# 5 - separate sites by climatic variables ----
# For this type of analysis, I should bring back in all of the indicator data - doesn't matter if some osils are missing SOC stock data or have different soil series for different reps. 

ggplot(site_clim, aes(x=mat, y=map, colour=project)) +
  geom_point()

# This makes a pretty decent grid - we can separate sites into warm and wet (MAP > 900 and MAT>15), warm and dry (MAT > 15, MAP < 900), cool and wet (MAT<15, MAP>900), and cool and dry (MAT<15, MAP<900)
site_clim2 <- site_clim %>%
  mutate(climate = ifelse(mat>15, ifelse(map>900, "warm_wet", "warm_dry"), ifelse(map>900, "cool_wet", "cool_dry"))) %>%
  group_by(soil) %>%
  fill(climate, .direction = "downup") %>% # fill in missing climate grouping for rows that don't have lat/long
  ungroup()

site_clim_distinct <- site_clim2 %>%
  distinct(soil, project, climate)

# 6 - Grouped analysis of indicator sensitivity in surface soils across all projects ----
# Attach climate data to horizon data
soc_horizon_clim <- soc_horizon_all %>%
  left_join(site_clim_distinct, by="soil")

# Filter for 0-5 cm soils, attach climate data
surface_all <-  soc_horizon_clim %>%
  filter(hrzdep_b == "5") 

surface_indicators_pivot <- surface_all %>%
  select(dsp_pedon_id, project, soil, label, climate, soc_stock_hrz, soc_pct, bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace) %>%
  pivot_longer(soc_stock_hrz:ace, names_to="indicator", values_to="value")

# ANOVA to test significance of relationship between management and indicators within each project (nested within climate)
surface_indicator_sig <- surface_indicators_pivot %>%
  na.omit() %>%
  group_by(climate, project, indicator) %>% 
  nest() %>%
  mutate(anova_obj = map(data, ~anova(aov(value ~ label, data=.x)))) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(climate, project, indicator, anova_tidy) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term=="label") %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant"))

# plot significance of response in each project grouped by climate
# Manually order projects to cooldry, coolwet, warmdry, and warmwet
surface_indicator_sig$project <- factor(surface_indicator_sig$project, 
                                        levels=c("TexasA&MPt-2", "UnivOfMinnesota","WashingtonState","KansasState",
                                                 "Illinois","OregonState", "UConn",
                                                 "UTRGV", "TexasA&MPt-1",
                                                 "NCState"), ordered=TRUE)
  
ggplot(surface_indicator_sig, aes(x=project, y=indicator, fill=sig)) +
  geom_tile() +
  scale_fill_viridis(discrete=TRUE) +
  geom_rect(xmin=c(0), xmax=c(4.5), ymin=c(0), ymax=c(16.5), color="black", fill=NA )+
  geom_rect(xmin=c(4.5), xmax=c(7.5), ymin=c(0), ymax=c(16.5), color="black", fill=NA ) +
  geom_rect(xmin=c(7.5), xmax=c(9.5), ymin=c(0), ymax=c(16.5), color="black", fill=NA )+
  geom_rect(xmin=c(9.5), xmax=c(10.5), ymin=c(0), ymax=c(16.5), color="black", fill=NA )

# This is potentially useful - I think this is a good start - but there are challenges with this analysis:
# not sure that indicator sensitivity is the most useful thing to look at. most indicators are sensitive to management in most projects, it seems that most indicators are sensitive to some extent
# variability mostly seems to exist within individual projects
