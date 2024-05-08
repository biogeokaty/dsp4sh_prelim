# 07 - Modeling relationships between SOC and indicators/soils/climate
# Katy Dynarski, March 2024
# 0 - Import data ----
meta_df <- read.csv(here("data_processed", "05_meta_df.csv"))

# 1 - Build dataframe for modeling ----
soc_df <- meta_df %>%
  mutate(across(where(is.character), as.factor)) %>%
  filter(!is.na(soc_stock_100cm)) # Filter out rows with no SOC stock data - they mess with the random forest model later on

# 2 - Ctree and cforest to identify variables influencing SOC stocks ----
# Try a ctree - with land use variables and top 5 soil health metrics
soc_tree <- ctree(soc_stock_100cm ~ soc_pct + ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                    label + lu + till + soil + mat + map + clay_tot_psa, data=soc_df)
plot(soc_tree)
# Soil series is the most important determinant, followed by SOC% and land use

# try a version that doesn't include soil series - want to look at the underlying factors integrated by soil series
soc_tree2 <- ctree(soc_stock_100cm ~ soc_pct + ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                     label + lu + till + mat + map + clay_tot_psa, data=soc_df)
plot(soc_tree2)
# SOC% is the most important thing...well duh. SOC stocks are driven by how much soil carbon is there. Should take that out of the tree...

soc_tree3 <- ctree(soc_stock_100cm ~ ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                     label + lu + till + mat + map + clay_tot_psa, data=soc_df)
plot(soc_tree3)
# Now we're getting somewhere...First split is MAT, followed by ACE/MAP, then land use and treatment

# Try as a forest
# Should include soil series in this analysis to see if this provides more data than MAT/MAP
soc_forest <- cforest(soc_stock_100cm ~ ace + bglucosidase + yoder_agg_stab_mwd + kssl_wsa + bglucosaminidase +
                        label + lu + till + mat + map + soil + clay_tot_psa, data=soc_df, ntree=30000)
soc_vi <- vi(soc_forest)
soc_vip <- vip(soc_vi, geom="col")
soc_vip

# Most important variables for predicting SOC stocks to 100 cm are: soil, MAT, ACE, MAP, land use, clay content, label, and tillage

# I'm not sure if the indicators (bglucosaminidase, bglucosidase, KSSL WSA, and Yoder agg stab) being less important is due to missing data for some of those measurements in some projects

# 3 - Linear models to predict SOC stocks based on top variables ----
stock_lm1 <- lm(soc_stock_100cm ~ soil + mat + ace + map + lu + clay_tot_psa + till, data=soc_df)
summary(stock_lm1) # This is a decent model, R2 of 0.86, p<0.001

# Try regular subsets to see if model can be reduced
stock_sub1 <- regsubsets(soc_stock_100cm ~ soil + mat + ace + map + lu + clay_tot_psa + till, data=soc_df)
summary(stock_sub1)
# if you're only going to keep three variables, they would be soil series, MAT, and ACE

stock_lm2 <- lm(soc_stock_100cm ~ soil + mat + ace, data=soc_df)
summary(stock_lm2)
# this is almost as good of a model as the full one, R2 is 0.78

stock_lm3 <- lm(soc_stock_100cm ~ mat + ace, data=soc_df)
summary(stock_lm3)
# Lose a lot of explanatory power when you take out the soil series (R2 down to 0.56)

stock_lm4 <- lm(soc_stock_100cm ~ mat + map + clay_tot_psa + ace, data=soc_df)
summary(stock_lm4)
# but can actually get a pretty good model wit MAT, MAP, clay, and ACE (R2=0.77), MAP is not significant

stock_lm5 <- lm(soc_stock_100cm ~ mat + clay_tot_psa + ace, data=soc_df)
summary(stock_lm5)
as.data.frame(partial_r2(stock_lm5))
# model is just as good with no MAP. R2=0.77. Temperature, clay, and ACE are closely correlated with total SOC stocks. 

# 4 - Determine best linear model to predict indicator values using glmulti ----
# Want to know if climate is a significant driver in other indicator values (it probably is)
ind_long <- meta_df %>%
  select(project, dsp_pedon_id, soil, label, climate, mat, map, lu, till, clay_tot_psa, soc_stock_0_30cm, soc_stock_100cm, 
         soc_pct:yoder_agg_stab_mwd, soil_respiration:ace) %>%
  pivot_longer(soc_pct:ace, names_to="indicator", values_to="value")

## 4.1 - glmulti - all sites ----
# use glmulti to find best model for each indicator, pull out the best model from the weight table, fit linear model, and extract fit metrics and coefficients
multi <- ind_long %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(multi_obj = map(data, ~glmulti("value", c("mat", "map", "lu", "till", "clay_tot_psa"), data=.x,
                                        level = 1, maxsize = 3, confsetsize = 128))) %>%
  mutate(weight = map(multi_obj, ~(weightable(.x) %>% slice_head(n=1)))) %>%
  transmute(indicator, data, weight) %>%
  unnest(cols = c(weight)) %>%
  mutate(lm_obj = map(data, ~lm(model, data=.x))) %>%
  mutate(part_r2 = map(lm_obj, ~as.data.frame(sensemakr::partial_r2(.x)) 
                       %>% rownames_to_column() %>% rename(var = "rowname", partial_r2 = "sensemakr::partial_r2(.x)"))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  mutate(lm_coefs = map(lm_obj, broom::tidy)) 

# Make top model for each indicator into a nice table
multi_table <- multi %>%
  ungroup() %>%
  transmute(indicator, model, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  select(indicator, model, adj.r.squared, p.value, AIC, nobs, sig) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) 
flextable(multi_table)
write_csv(multi_table, here("figs", "indicator_lm.csv"))

# Pull out the partial R2 values for all models
multi_partials <- multi %>%
  ungroup() %>%
  transmute(indicator, model, lm_tidy, part_r2) %>%
  unnest(cols=c(lm_tidy, part_r2)) %>%
  select(indicator, model, adj.r.squared, var, partial_r2) %>%
  filter(var!="(Intercept)") %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 
flextable(multi_partials)
write_csv(multi_partials, here("figs", "indicator_partial_r2.csv"))

## 4.2 - plot partial R2 - all sites ----
multi_partials_2plot <- multi_partials %>%
  filter(partial_r2 > 0.1) %>%
  mutate(var = factor(var, levels=c("luRANGE", "luGRASS", "luFOREST", "luFORAGE", "tillPerennial", "clay_tot_psa", "map", "mat")))

var_fill <- c("mat" = "#35978f",
              "map" = "#80cdc1",
              "clay_tot_psa" = "#c7eae5",
              "tillPerennial" = "#f5f5f5",
              "luFORAGE" = "#f6e8c3",
              "luFOREST" = "#dfc27d",
              "luGRASS" = "#bf812d",
              "luRANGE" = "#8c510a")

ggplot(multi_partials_2plot, aes(fill=var, 
                                 y=partial_r2, x=fct_reorder(indicator, adj.r.squared))) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = partial_r2), size=3, position=position_stack(vjust = 0.5)) +
  labs(x="Indicator", y=expression(Partial~R^2)) +
  scale_x_discrete(labels=indicator_labs) +
  scale_fill_manual(name="Predictor variable", 
                    values=var_fill,
                    labels=c("LU - Rangeland", "LU - Grass", "LU - Forest", "LU - Forage", "Tillage - Untilled Perennial",
                             "Soil Clay %", "MAP", "MAT"),
                    guide = guide_legend(reverse = TRUE),) +
  theme_katy() +
  coord_flip()
ggsave(here("figs", "indicator_partial_r2.png"), height=5, width=9, units="in", dpi=400)

## 4.3 - glmulti - reference sites only ----
multi_ref <- ind_long %>%
  filter(label=="Ref") %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(multi_obj = map(data, ~glmulti("value", c("mat", "map", "lu", "clay_tot_psa"), data=.x,
                                      level = 1, maxsize = 3, confsetsize = 128))) %>%
  mutate(weight = map(multi_obj, ~(weightable(.x) %>% slice_head(n=1)))) %>%
  transmute(indicator, data, weight) %>%
  unnest(cols = c(weight)) %>%
  mutate(lm_obj = map(data, ~lm(model, data=.x))) %>%
  mutate(part_r2 = map(lm_obj, ~as.data.frame(sensemakr::partial_r2(.x)) 
                       %>% rownames_to_column() %>% rename(var = "rowname", partial_r2 = "sensemakr::partial_r2(.x)"))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  mutate(lm_coefs = map(lm_obj, broom::tidy)) 

# Make top model for each indicator into a nice table
multi_table_ref <- multi_ref %>%
  ungroup() %>%
  transmute(indicator, model, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  select(indicator, model, adj.r.squared, p.value, AIC, nobs, sig) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) 
flextable(multi_table_ref)
write_csv(multi_table_ref, here("figs", "indicator_lm_ref.csv"))

# Pull out the partial R2 values for all models
multi_partials_ref <- multi_ref %>%
  ungroup() %>%
  transmute(indicator, model, lm_tidy, part_r2) %>%
  unnest(cols=c(lm_tidy, part_r2)) %>%
  select(indicator, model, adj.r.squared, var, partial_r2) %>%
  filter(var!="(Intercept)") %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 
flextable(multi_partials_ref)
write_csv(multi_partials, here("figs", "indicator_partial_r2_ref.csv"))

## 4.4 - plot partial R2 - reference sites only ----
multi_partials_2plot_ref <- multi_partials_ref %>%
  filter(partial_r2 > 0.1) %>%
  mutate(var = factor(var, levels=c("luRANGE", "luGRASS", "luFOREST", "luORCHARD", "clay_tot_psa", "map", "mat")))

var_fill_ref <- c("mat" = "#35978f",
              "map" = "#80cdc1",
              "clay_tot_psa" = "#c7eae5",
              "luORCHARD" = "#f6e8c3",
              "luFOREST" = "#dfc27d",
              "luGRASS" = "#bf812d",
              "luRANGE" = "#8c510a")

ggplot(multi_partials_2plot_ref, aes(fill=var, 
                    y=partial_r2, x=fct_reorder(indicator, adj.r.squared))) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = partial_r2), size=3, position=position_stack(vjust = 0.5)) +
  labs(x="Indicator", y=expression(Partial~R^2)) +
  scale_x_discrete(labels=indicator_labs) +
  scale_fill_manual(name="Predictor variable", 
                    values=var_fill_ref,
                    labels=c("LU - Rangeland", "LU - Grass", "LU - Forest", "LU - Orchard",
                             "Soil Clay %", "MAP", "MAT"),
                    guide = guide_legend(reverse = TRUE),) +
  theme_katy() +
  coord_flip()
ggsave(here("figs", "indicator_partial_r2_ref.png"), height=5, width=9, units="in", dpi=400)

# 5 - ARCHIVE - manually calling lm() for indicators ----
# Make a version of the model with soil, MAT/MAP, land use, and tillage as predictors
ind_lm <- ind_long %>%
  na.omit() %>%
  group_by(indicator) %>% 
  nest() %>%
  mutate(lm_obj = map(data, ~lm(value ~ soil + mat + map + lu + till + clay_tot_psa, data=.x))) %>%
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
  mutate(lm_obj = map(data, ~lm(value ~ mat + map + lu + till + clay_tot_psa, data=.x))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  ungroup() %>%
  transmute(indicator, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
flextable(ind_lm_clim)

# and a version with only MAT and MAP and clay
ind_lm_clim2 <- ind_long %>%
  na.omit() %>%
  group_by(indicator) %>% 
  nest() %>%
  mutate(lm_obj = map(data, ~lm(value ~ mat + map + clay_tot_psa, data=.x))) %>%
  mutate(lm_tidy = map(lm_obj, broom::glance)) %>%
  ungroup() %>%
  transmute(indicator, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
flextable(ind_lm_clim2)

# Overall findings from this:
# Soil has the most explanatory power - R2 is always higher when soil is included as a variable in the model. Climate on its own still explains a lot, but without land use and tillage as factors, R2 maxes out at .69 (which is still pretty high to be honest)
# What is the best way to display this information??