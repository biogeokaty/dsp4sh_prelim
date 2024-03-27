# 07 - Modeling relationships between SOC and indicators/soils/climate
# Katy Dynarski, March 2024
# 0 - Import data ----
meta_df <- read.csv(here("data_processed", "05_meta_df.csv"))

# 1 - Build dataframe for modeling ----
soc_df <- meta_df %>%
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
# Now we're getting somewhere...First split is MAT, followed by ACE/MAP, then land use and treatment

# Try as a forest
# Should include soil series in this analysis to see if this provides more data than MAT/MAP
soc_forest <- cforest(soc_stock_100cm ~ ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                        label + lu + till + mat + map + soil, data=soc_df, ntree=10000)
soc_vi <- vi(soc_forest)
soc_vip <- vip(soc_vi, geom="col")
soc_vip

# Most important variables for predicting SOC stocks to 100 cm are: soil, MAT, ACE, MAP, land use, label, and tillage

# I'm not sure if the indicators (bglucosaminidase, bglucosidase, KSSL WSA, and Yoder agg stab) being less important is due to missing data for some of those measurements in some projects

# 4 - Linear models to predict SOC stocks based on top variables ----
stock_lm1 <- lm(soc_stock_100cm ~ soil + mat + ace + map + lu + till, data=soc_df)
summary(stock_lm1) # This is a decent model, R2 of 0.78, p<0.001

# Try regular subsets to see if model can be reduced
stock_sub1 <- regsubsets(soc_stock_100cm ~ soil + mat + ace + map + lu + till, data=soc_df)
summary(stock_sub1)
# if you're only going to keep three variables, they would be soil series, MAT, and MAP

stock_lm2 <- lm(soc_stock_100cm ~ soil + mat + map, data=soc_df)
summary(stock_lm2)
# this is almost as good of a model as the full one, R2 is 0.76

stock_lm3 <- lm(soc_stock_100cm ~ mat + map, data=soc_df)
summary(stock_lm3)
# Lose a lot of explanatory power when you take out the 

# 5 - Teasing out importance of climate vs soil in driving SOC stocks ----

# try mixed linear model with variance partitioning
soc_mixed_clim <- lmer(soc_stock_100cm ~ (1|soil) + (1|mat) + (1|map), data = soc_df)
summary(soc_mixed_clim)
# Results of mixed model: soil series and climate each explain about a third of variability not explained by treatment, climate actually explains more (that seems pretty good to me???)

soc_mixed_clim_var <- data.frame(VarCorr(soc_mixed_clim)) %>%
  select(grp, vcov) %>%
  rename(group = grp,
         variance = vcov) 

tot_var <- soc_mixed_clim_var %>%
  summarize(sum_var = sum(variance))
temp <- cbind(soc_mixed_clim_var, tot_var)

soc_mixed_clim_var_pct <- temp %>%
  mutate(pct_var = (variance/sum_var)*100) %>%
  select(group, pct_var)

# could also try partial correlations...

# 6 - Climate patterns in other indicators ----
# Want to know if climate is a significant driver in other indicator values (it probably is)
ind_long <- meta_df %>%
  select(project, dsp_pedon_id, soil, label, climate, mat, map, lu, till, soc_stock_0_30cm, soc_stock_100cm, 
         soc_pct:yoder_agg_stab_mwd, soil_respiration:ace) %>%
  pivot_longer(soc_pct:ace, names_to="indicator", values_to="value")

# Make a version of the model with soil, MAT/MAP, land use, and tillage as predictors
ind_lm <- ind_long %>%
  na.omit() %>%
  group_by(indicator) %>% 
  nest() %>%
  mutate(lm_obj = map(data, ~lm(value ~ soil + mat + map + lu + till, data=.x))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  ungroup() %>%
  transmute(indicator, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
flextable(ind_lm)

# try a version with no soil
ind_lm_clim <- ind_long %>%
  na.omit() %>%
  group_by(indicator) %>% 
  nest() %>%
  mutate(lm_obj = map(data, ~lm(value ~ mat + map + lu + till, data=.x))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  ungroup() %>%
  transmute(indicator, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
flextable(ind_lm_clim)

# and a version with only MAT and MAP
ind_lm_clim2 <- ind_long %>%
  na.omit() %>%
  group_by(indicator) %>% 
  nest() %>%
  mutate(lm_obj = map(data, ~lm(value ~ mat + map, data=.x))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  ungroup() %>%
  transmute(indicator, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
flextable(ind_lm_clim2)

# Overall findings from this:
# Soil has the most explanatory power - R2 is always higher when soil is included as a variable in the model. Climate on its own still explains a lot, but without land use and tillage as factors, R2 maxes out at .5 (which is still pretty high to be honest)
# What is the best way to display this information??