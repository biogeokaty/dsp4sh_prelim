# 07 - Big Silly Models
# Katy Dynarski, December 2023/January 2024

# 0 - Import data ----
soc_pedon100 <- read.csv(here("data_processed", "05_soc_pedon_100cm.csv"))
soc_horizon_filt <- read.csv(here("data_processed", "05_soc_horizon_filtered.csv"))
soc_horizon_all <- read.csv(here("data_processed", "05_soc_horizon_clim_infilt.csv"))
surf <- read.csv(here("data_processed", "05_surface_horizons.csv"))

# 1 - Build dataframe for modeling ----
# Want from pedon df: SOC stock to 100 cm, soil, label, lu, till, MAT, MAP, infiltration
# want from surface df: all the indicators
# want but don't (yet) have: other state factor data - what should this be??

surf_keep <- surf %>%
  select(dsp_pedon_id, soc_stock_hrz, soc_pct, bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace) # same indicators as before

soc_df <- soc_pedon100 %>%
  left_join(surf_keep, by="dsp_pedon_id") %>%
  mutate(across(where(is.character), as.factor)) %>%
  filter(!is.na(soc_stock_100cm)) # Filter out rows with no SOC stock data - they mess with the random forest model later on

# 3 - Ctree and cforest ----
# Try a ctree - with land use variables and top 5 soil health metrics
soc_tree <- ctree(soc_stock_100cm ~ soc_pct + ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                         label + lu + till + soil + mat + map + cornell_infiltrometer, data=soc_df)
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


# Try using climate as a factor in pedon data - what is the best way to determine how much additional varaibility is explained by soil vs climate?
# stepwise regression?
soc_pedon100_clim <- soc_pedon100 %>%
  left_join(site_clim_distinct, by=c("soil", "project"))

soc_mixed_clim <- lmer(soc_stock_100cm ~ label + (1|soil) + (1|climate), data = soc_pedon100_clim)
summary(soc_mixed_clim)
# Results of mixed model: soil series and climate each explain about a third of variability not explained by treatment (that seems pretty good to me???)

# What are some other good statistical tools for teasing out the role of climate vs soil?
# I think I made a plot of this for the stoichiometry paper?? try variance partitioning...

