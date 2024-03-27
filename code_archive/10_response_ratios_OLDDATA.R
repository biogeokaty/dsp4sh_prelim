# Response ratios

# 0 - Import data ----
meta_df <- read.csv(here("data_processed", "05_meta_df.csv"))
site_clim <- read.csv(here("data_processed", "05_site_clim.csv"))

mean_sd <- list(
  mean = ~round(mean(.x, na.rm = TRUE), 2), 
  sd = ~round(sd(.x, na.rm = TRUE), 2)
)

# 1 - Calculate response ratios of indicators ----
# need to summarize and get mean, sd, and n for each group in each study
# will need to calculate separate effect sizes for each 
# ....FOR EVERY VARIABLE dun dun dun
# think we can accomplish this with some sort of nesting structure in a map function

# let's test with soc stocks to get the code right

soc_es_in <- meta_df %>%
  group_by(project, label) %>%
  summarize(across(soc_stock_100cm, mean_sd), n=n()) %>%
  pivot_wider(names_from = label, values_from=soc_stock_100cm_mean:n)

soc_es_bau_shm <- escalc(n1i = n_BAU, n2i = n_SHM, m1i = soc_stock_100cm_mean_BAU, m2i = soc_stock_100cm_mean_SHM, 
                  sd1i = soc_stock_100cm_sd_BAU, sd2i = soc_stock_100cm_sd_SHM, data = soc_es_in, measure = "ROM")
hist(soc_es_bau_shm$yi) # histogram of response ratios
forest(soc_es_bau_shm$yi, soc_es_bau_shm$vi) # forest plot of response ratios

fixed <- rma(yi=yi, vi=vi, data = soc_es_bau_shm, method="FE") 
summary(fixed)
forest.rma(fixed)

random <- rma(yi, vi, data=soc_es_bau_shm, method="REML")
summary(random)

# 2 - nested dataframe to calculate effect sizes for many variables ----
# first need to pivot datalonger
meta_long <- meta_df %>%
  select(project, dsp_pedon_id, soil, label, climate, lu, till, trt, soc_stock_30cm:soc_stock_100cm, single_ring_infiltration_1:ace) %>%
  pivot_longer(soc_stock_30cm:ace, names_to="indicator", values_to="value")

# build a nested dataframe and summarize data within each study/label/indicator group with a map() function - this builds the dataframe that can go into escalc() to calculate effect sizes
nested <- meta_long %>%
  group_by(indicator, label, project) %>%
  nest() %>%
  mutate(
    mean = map_dbl(data, ~mean(.$value, na.rm = TRUE)),
    sd = map_dbl(data, ~sd(.$value, na.rm = TRUE)),
    n = map_dbl(data, ~sum(!is.na(.$value)))) %>% 
  select(indicator, label, project, mean, sd, n) %>%
  pivot_wider(names_from = label, values_from=mean:n) %>%
  na.omit()

# calculate effect sizes for all indicators/studies contrasting BAU and SHM groups
es_bau_shm <- escalc(n1i = n_SHM, n2i = n_BAU, m1i = mean_SHM, m2i = mean_BAU, 
                         sd1i = sd_SHM, sd2i = sd_BAU, data = nested, measure = "ROM")

# calculate effect sizes for all indicators/studies contrasting BAU and Ref groups
es_bau_ref <- escalc(n1i = n_Ref, n2i = n_BAU, m1i = mean_Ref, m2i = mean_BAU, 
                     sd1i = sd_Ref, sd2i = sd_BAU, data = nested, measure = "ROM")  

# put two frames together to have one dataframe of all effect sizes
es_all <- es_bau_shm %>%
  left_join(select(es_bau_ref, indicator, project, yi, vi), by=c("indicator", "project"), suffix = c("_shm", "_ref"))

# 3 - Calculate random-effects models for all indicators and plot ----
# need to reshape data once again!
es_long <- es_all %>%
  select(indicator, project, yi_shm:vi_ref) %>%
  group_by(indicator, project) %>%
  pivot_longer(yi_shm:vi_ref,
               names_to = c(".value", "label"), 
               names_sep = "_")
# now should be able to nest and plot

# Test plot of effect sizes for each indicator
ggplot(es_long, aes(x=label, y=yi)) +
  geom_boxplot() +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(vars(indicator)) +
  theme_katy()
# Visually, it looks like the indicators most sensitive to management are ACE, water-stable aggregates (KSSL), arylsulfatase, bglucosaminidase, POXC, SOC%, and aggregate stability (Yoder)

# Can test this out a bit better by fitting random effects models

# nest data by indicators and treatment contrasts, calculate random-effects models, extract data from model and make dataframe needed for plotting forest plot of each indicator and contrast
es_rma <- es_long %>%
  group_by(indicator, label) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab=project, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2)))

# test how this works by plucking the first dataframe (SOC stock to 30cm, SHM vs BAU)
test <- es_rma %>% pluck("plot_df", 1) %>%
  mutate(study = factor(study, levels=rev(study))) # make study name an ordered factor so that summary will go on the bottom of the plot

# test forest plot
ggplot(data=test,aes(x=factor(study, levels=rev(study)),y=es,ymax=es+se,ymin=es-se,size=factor(type),colour=factor(type))) + 
  geom_pointrange() +
  coord_flip() + 
  geom_hline(yintercept=0, lty=2,linewidth=1) +
  scale_size_manual(values=c(0.5,1)) + 
  labs(xlab="Study", ylab="log response ratio") + 
  scale_colour_manual(values=c("grey","black")) + 
  theme_katy() +
  theme(legend.position="none")

# now to make lots and lots of these bwuahaha
es_rma_plot_df <- es_rma %>%
  select(indicator, label, plot_df) %>%
  ungroup() %>%
  transmute(indicator, label, plot_df) %>%
  unnest(cols = c(plot_df))

# make vector for each indicator/label combination
indicator_label <- es_rma_plot_df %>% distinct(indicator, label)
indicators <- indicator_label %>% pull(indicator) %>% as.character
labels <- indicator_label %>% pull(label) %>% as.character

# make forest plot for all combinations of indicators and management contrasts (Ref vs BAU and SH vs BAU)
plot_list <- map2(.x = indicators,
                 .y=labels,
                 .f = ~{
                   es_rma_plot_df %>% 
                     filter(indicator == .x, label == .y) %>%
                     ggplot(aes(x=factor(study, levels=rev(study)),y=es,ymax=es+se,ymin=es-se,size=factor(type),colour=factor(type))) + 
                     geom_pointrange() +
                     coord_flip() + 
                     geom_hline(yintercept=0, lty=2,linewidth=1) +
                     scale_size_manual(values=c(0.5,1)) + 
                     labs(x="Study", y="log response ratio",
                          title=glue::glue({unique(filter(es_rma_plot_df, indicator==.x & label==.y)$indicator)}, 
                                           " - ", {unique(filter(es_rma_plot_df, indicator==.x & label==.y)$label)})) + 
                     scale_colour_manual(values=c("grey","black")) + 
                     theme_katy() +
                     theme(legend.position="none")
                   
                   ggsave(here("figs", "indicator effect sizes",glue::glue("es_", .x, "_", .y, ".png")), 
                          width=8, height=5.5, units="in", dpi=400)
                 })

# Extract summary data from random-effects model, label contrasts as either significant or not significant
es_rma_sig <- es_rma %>%
  ungroup() %>%
  transmute(indicator, label, rma_tidy) %>%
  unnest(cols = c(rma_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(estimate))
# Largest response ratio estimates are for: infiltration rings (though - this is only based on a few projects), arylsulfatase, aggregate stability (both methods), bglucosaminidase, bglucosidase. Interesting that so many enzymes are fairly sensitive!

# Big indicators seem to be enzyme activity and aggregate stability

# 4 - Add moderator variables to dataframe and perform meta-regression ----

# Eventually I want to incorporate land-use variables in this - but they are pretty complex because there are typically multiple land uses and management strategies in each project. Will need to calculate even more comparisons!
# for now - incorporate mean MAT and MAP at site, also climate category - can see if climate grouping adds useful information 

# input data for random effects model is : es_long
# need: mean MAP and MAT data, climate grouping for each project

site_clim_dist <- site_clim %>%
  distinct(project, climate)

site_clim_sum <- site_clim %>%
  group_by(project) %>%
  summarize(across(mat:map, ~ mean(.x, na.rm = TRUE))) %>%
  left_join(site_clim_dist, by="project")

es_long_mod <- es_long %>% 
  ungroup() %>%
  filter(indicator!="cornell_infiltrometer", indicator!="single_ring_infiltration_1", 
         indicator!="single_ring_infiltration_2", indicator!="p_h") %>%
  left_join(site_clim_sum, by="project")
  

es_rma_mod <- es_long_mod %>%
  group_by(indicator, label) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, mod= ~climate, slab=project, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2,
                                           r2 = .x$R2)))

# Compare I^2 value between simpler and more complex model to see if climate categories added any information
# original model summary
summary(es_rma %>% pluck("rma_obj", 1))
# more complex model summary
summary(es_rma_mod %>% pluck("rma_obj", 1))
# I^2 went from 83% to 80% - not a big difference
# The I^2 value reflects how much of the overall variability in data is due to heterogeneity between studies - in this case, climate groupings don't explain much of the heterogeneity between studies

# try with climate as a continuous variable?
es_rma_mod2 <- es_long_mod %>%
  group_by(indicator, label) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, mod= ~mat+map, slab=project, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2,
                                         r2 = .x$R2)))

# pull out a model summary to examine - look at ACE in Ref, which had a large overall estimated effect size but a lot of variability
summary(es_rma %>% pluck("rma_obj", 28)) # model with no moderating variables
summary(es_rma_mod %>% pluck("rma_obj", 28)) # model with climate as a category, actually raises I^2
summary(es_rma_mod2 %>% pluck("rma_obj", 28)) # model with MAT and MAP as continuous vars, accounts for 70% of heterogeneity!

# need to pull out I2 and R2 for these models and put into some kind of table to identify for which indicators is the response to management moderated by climate

es_rma_nomod_comp_df <- es_rma %>%
  select(indicator, label, het_df) %>%
  ungroup() %>%
  transmute(indicator, label, het_df) %>%
  unnest(cols = c(het_df)) %>%
  rename(i2_nomod = i2)

es_rma_cat_comp_df <- es_rma_mod %>%
  select(indicator, label, het_df) %>%
  ungroup() %>%
  transmute(indicator, label, het_df) %>%
  unnest(cols = c(het_df)) %>%
  rename(i2_cat = i2,
         r2_cat = r2)

es_rma_cont_comp_df <- es_rma_mod2 %>%
  select(indicator, label, het_df) %>%
  ungroup() %>%
  transmute(indicator, label, het_df) %>%
  unnest(cols = c(het_df)) %>%
  rename(i2_cont = i2,
         r2_cont = r2)

es_rma_comp_df <- es_rma_nomod_comp_df %>%
  left_join(es_rma_cat_comp_df, by=c("indicator", "label")) %>%
  left_join(es_rma_cont_comp_df, by=c("indicator", "label")) %>%
  arrange(desc(r2_cont))

# Results are super interesting - response ratio variability for certain indicators is really well-explained by climate. For example, including climate as either a continuous or categorical variable explains 100% of the heterogeneity in bulk density response.   

ggplot(filter(es_long_mod, indicator=="bglucosidase", label=="shm"), aes(x=climate, y=yi)) +
  geom_boxplot() +
  theme_katy()

# Looking at some plots, this is a little tricky. Sometimes when the moderating variable explains a lot, it's because there's not many data points and there's only one project in each climate category. So climate looks like it explains all of the heterogeneity, but it's because each project is assigned to a different climate. For that reason, I think it's probably not a good idea to use climate as a categorical variable - better to stick with continuous variables. 

# 5 - Calculate effect size for individual treatments and plot indicator effect sizes ----
# ultimate goal of this section of code is to end up with different land uses and tillage as different rows 
# this actually involves re-calculating effect sizes - instead of having one mean effect size for each SHM and Ref treatment within a project, will separate by the different cropping types and land uses etc within each project
# arguably, this is the correct way to organize this data anyway. If I were doing a traditional meta-analysis, I would definitely record the different land uses as different cases and give each their own row
meta_long2 <- meta_long %>%
  mutate(treat = ifelse(label=="BAU", "BAU", "ASP"))

bau_means <- meta_long2 %>%
  filter(treat=="BAU") %>%
  filter(dsp_pedon_id!="BAU9-1") %>% # exclude pedon BAU9-1 from analysis - it's from WashingtonState study and has implausibly high enzyme activity
  group_by(project, indicator) %>%
  nest() %>%
  mutate(
    mean_bau = map_dbl(data, ~mean(.$value, na.rm = TRUE)),
    sd_bau = map_dbl(data, ~sd(.$value, na.rm = TRUE)),
    n_bau = map_dbl(data, ~sum(!is.na(.$value)))) %>%
  select(indicator, project, mean_bau, sd_bau, n_bau) %>%
  ungroup() 

nested2 <- meta_long2 %>%
  filter(treat=="ASP") %>%
  group_by(indicator, label, project, lu, till, trt) %>%
  nest() %>%
  mutate(
    mean_asp = map_dbl(data, ~mean(.$value, na.rm = TRUE)),
    sd_asp = map_dbl(data, ~sd(.$value, na.rm = TRUE)),
    n_asp = map_dbl(data, ~sum(!is.na(.$value)))) %>%
  select(indicator, label, project, lu, till, trt, mean_asp, sd_asp, n_asp) %>%
  ungroup()

es_in_big <- nested2 %>%
  left_join(bau_means, by=c("indicator", "project")) %>%
  na.omit()

es_asp1 <- escalc(n1i = n_asp, n2i = n_bau, m1i = mean_asp, m2i = mean_bau, 
                     sd1i = sd_asp, sd2i = sd_bau, data = es_in_big, measure = "ROM") 

ggplot(es_asp1 %>% filter(indicator=="soc_pct"), aes(x=till, y=yi, colour=label)) +
  geom_boxplot()

# why is the SOC response so high to "Till" - looked at data, and it's all from Texas A&M project that is coded as "Till" but the treatment name is no-till....should fix

es_asp <- es_asp1 %>%
  na.omit() %>% # omit any rows with NA - this should just be places where NaNs were calculated
  unite("proj_trt", project, trt, remove = FALSE) %>%
  rename(till_orig = till) %>%
  mutate(till = ifelse(till_orig == "Till", "No-till", till_orig)) %>% # fix Texas A&M data that was mislabeled
  relocate(till, .before=till_orig) %>%
  filter(indicator!="cornell_infiltrometer", indicator!="single_ring_infiltration_1", 
         indicator!="single_ring_infiltration_2", indicator!="p_h") %>% # remove indicators that are only in one or two studies
  left_join(site_clim_sum, by="project") # add in climate data


# check by plotting
ggplot(es_asp %>% filter(indicator=="soc_pct"), aes(x=till, y=yi, colour=label)) +
  geom_boxplot()

# Now we have projects properly separated, with tillage and land use data as well as climate!

# Do I get the same results as before

es_asp_rma_nomod <- es_asp %>%
  group_by(indicator, label) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = proj_trt, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2)))

# pull out plotting data
es_asp_rma_plot_df <- es_asp_rma_nomod %>%
  select(indicator, label, plot_df) %>%
  ungroup() %>%
  transmute(indicator, label, plot_df) %>%
  unnest(cols = c(plot_df))

# make vectors of indicator/label combination
indicator_label2 <- es_asp_rma_plot_df %>% distinct(indicator, label)
indicators2 <- indicator_label2 %>% pull(indicator) %>% as.character
labels2 <- indicator_label2 %>% pull(label) %>% as.character

# make plots - should have the same overall effect sizes for each indicator, but more lines to show the different treatments
plot_list <- map2(.x = indicators2,
                  .y=labels2,
                  .f = ~{
                    es_asp_rma_plot_df %>% 
                      filter(indicator == .x, label == .y) %>%
                      ggplot(aes(x=factor(study, levels=rev(study)),y=es,ymax=es+se,ymin=es-se,size=factor(type),colour=factor(type))) + 
                      geom_pointrange() +
                      coord_flip() + 
                      geom_hline(yintercept=0, lty=2,linewidth=1) +
                      scale_size_manual(values=c(0.5,1)) + 
                      labs(x="Study", y="log response ratio",
                           title=glue::glue({unique(filter(es_asp_rma_plot_df, indicator==.x & label==.y)$indicator)}, 
                                            " - ", {unique(filter(es_asp_rma_plot_df, indicator==.x & label==.y)$label)})) + 
                      scale_colour_manual(values=c("grey","black")) + 
                      theme_katy() +
                      theme(legend.position="none")
                    
                    ggsave(here("figs", "indicator effect sizes 2",glue::glue("es_", .x, "_", .y, ".png")), 
                           width=8, height=5.5, units="in", dpi=400)
                  })

# pull out significant indicators and make a table
es_asp_rma_sig <- es_asp_rma_nomod %>%
  ungroup() %>%
  transmute(indicator, label, rma_tidy) %>%
  unnest(cols = c(rma_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  group_by(label) %>%
  arrange(p.value)
gt(es_asp_rma_sig, groupname_col="label")

# 6 - Testing moderator variables in the treatment-specific data ----

# test random effects model with a subset
es_asp_sub <- es_asp %>%
  filter(indicator=="ace")
random1 <- rma(yi, vi, data=es_asp_sub, mod= ~ mat + map + lu + till, method="REML")
summary(random1)
forest.rma(random1)

# One version of the model - group by indicators only (so Ref and SHM are tested in the same model), and include MAT, MAP, land use, and tillage as moderators
es_asp_rma <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ mat + map + lu + till, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2,
                                           r2 = .x$R2))) # pull out I2 and R2

# Pull out dataframe with I2 and R2 values to see for which indicators does adding moderators to the model add significant explanatory power. Then we can pull these models out and look at them individually.
es_asp_rma_het <- es_asp_rma %>%
  select(indicator, het_df) %>%
  ungroup() %>%
  transmute(indicator,  het_df) %>%
  unnest(cols = c(het_df)) %>%
  arrange(desc(r2))
# highest R2 is with KSSL WSA

# Pluck out KSSL WSA model to look at predictors
summary(es_asp_rma %>% pluck("rma_obj", 5))
# Interesting thing is that none of the individual moderators are actually significant (land use range is almost significant, p=0.06). I wonder if this model just explains heterogeneity because it essentially makes each study a unique combination of categorical moderators ?

# 7 - Extract significant moderator variables for each indicator and plot ---- 
# Make a table of significant moderator variables for each indicator
es_asp_rma_sig <- es_asp_rma %>%
  ungroup() %>%
  transmute(indicator, rma_tidy) %>%
  unnest(cols = c(rma_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  filter(sig=="significant") %>%
  filter(term!="intercept") %>%
  select(!type) %>%
  group_by(indicator) %>%
  arrange(p.value, .by_group = TRUE)
gt(es_asp_rma_sig)
# The interesting thing is that some of the indicators where the random-effects model explained the most heterogeneity (e.g. KSSL WSA) do not actually have any significant moderating variables - that to me is a sign the model is overfit. Tillage is never a significant variable, so that really could be removed

# Plot significant moderator variables
# MAT - Significant for ACE, acid phosphatase, arylsulfatase, bglucosaminidase, bglucosidase, SOC%
# Make vector of indicators 
mat_indicators <- c("ace", "acid_phosphatase", "arylsulfatase", "bglucosaminidase", "bglucosidase", "soc_pct")

# Plot with map function
mat_plots <- map(.x = mat_indicators,
                  .f = ~{
                    es_asp %>% 
                      filter(indicator == .x) %>%
                      ggplot(aes(x=mat, y=yi)) +
                      geom_point(aes(colour=lu)) +
                      geom_smooth(method="lm", formula = y~x, color="black") +
                      geom_hline(yintercept=0, linetype="dashed") +
                      stat_cor(aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
                               size=4) +
                      scale_colour_discrete(name="Land use") +
                      labs(x=expression("Mean annual temperature"~(degree*C)), y="Log response ratio",
                           title=glue::glue({unique(filter(es_asp, indicator==.x)$indicator)})) + 
                      theme_katy_grid()
                  })

# Now make into a panel grid with cowplot
mat_ace <- pluck(mat_plots, 1)
mat_acidp <- pluck(mat_plots, 2)
mat_aryl <- pluck(mat_plots, 3)
mat_bgm <- pluck(mat_plots, 4)
mat_bg <- pluck(mat_plots, 5)
mat_soc <- pluck(mat_plots, 6)

mat_grid <- plot_grid(mat_ace + theme(legend.position="none", axis.title.x=element_blank()),
                      mat_acidp + theme(legend.position="none", axis.title.y=element_blank()),
                      mat_aryl + theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank()),
                      mat_bgm + theme(legend.position="none", axis.title.x=element_blank()),
                      mat_bg + theme(legend.position="none", axis.title.y=element_blank()),
                      mat_soc + theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank()),
                      align = 'vh',
                      labels = c("A", "B", "C", "D", "E", "F"),
                      hjust = -1,
                      nrow = 2)

mat_leg <- get_legend(mat_ace + theme(legend.box.margin = margin(0, 0, 0, 12)))

plot_grid(mat_grid, mat_leg, rel_widths = c(3, .4))
ggsave(here("figs", "mat_mod_reg.png"), width=11, height=8, units="in")

# MAP - significant for POX-C and respiration
# Make vector of indicators 
map_indicators <- c("pox_c", "soil_respiration")

# Plot with map function
map_plots <- map(.x = map_indicators,
                 .f = ~{
                   es_asp %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x=map, y=yi)) +
                     geom_point(aes(colour=lu)) +
                     geom_smooth(method="lm", formula = y~x, color="black") +
                     geom_hline(yintercept=0, linetype="dashed") +
                     stat_cor(aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
                              size=4) +
                     scale_colour_discrete(name="Land use") +
                     labs(x=expression("Mean annual precipitation"~(mm~yr^-1)), y="Log response ratio",
                          title=glue::glue({unique(filter(es_asp, indicator==.x)$indicator)})) + 
                     theme_katy_grid()
                 })

# Now make into a panel grid with cowplot
map_poxc <- pluck(map_plots, 1)
map_resp <- pluck(map_plots, 2)

map_grid <- plot_grid(map_poxc + theme(legend.position="none", axis.title.x=element_blank()),
                      map_resp + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()),
                      align = 'vh',
                      labels = c("A", "B"),
                      hjust = -1,
                      nrow = 1)

map_grid_label <- ggdraw(map_grid) +
  draw_label(expression("Mean annual precipitation"~(mm~yr^-1)), color="black", vjust=0, y=0)

map_leg <- get_legend(map_poxc + theme(legend.box.margin = margin(0, 0, 0, 12)))

plot_grid(map_grid_label, map_leg, rel_widths = c(3, .4))
ggsave(here("figs", "map_mod_reg.png"), width=10, height=8, units="in")

# Land use - significant for arylsulfatase, bglucosaminidase, bulk density, phosphodiesterase, soil respiration
# Make vector of indicators
lu_indicators <- c("arylsulfatase", "bglucosaminidase", "bulk_density", "phosphodiesterase","soil_respiration")

# Plot with map function
lu_plots <- map(.x = lu_indicators,
                 .f = ~{
                   es_asp %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x=lu, y=yi, colour=lu)) +
                     geom_boxplot() +
                     geom_hline(yintercept=0, linetype="dashed") +
                     labs(x="Land use", y="Log response ratio",
                          title=glue::glue({unique(filter(es_asp, indicator==.x)$indicator)})) + 
                     theme_katy_grid()
                 })

# Now make into a panel grid with cowplot
lu_aryl <- pluck(lu_plots, 1)
lu_bgm <- pluck(lu_plots, 2)
lu_bd <- pluck(lu_plots, 3)
lu_pd <- pluck(lu_plots, 4)
lu_resp <- pluck(lu_plots, 5)

lu_grid <- plot_grid(lu_aryl + theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     lu_bgm + theme(legend.position="none", axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1)), 
                     lu_bd + theme(legend.position="none", axis.title.y=element_blank(), 
                                   axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     lu_pd + theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     lu_resp + theme(legend.position="none", axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                      align = 'vh',
                      labels = c("A", "B", "C", "D", "E"),
                      hjust = -1,
                      nrow = 2)
lu_grid
ggsave(here("figs", "lu_mod_box.png"), width=11, height=8, units="in")


# 8 - Alternate versions of random-effects models with moderating variables with less explanatory power ----
# try another version of the model with MAT, MAP, and label instead of land use and tillage (more parsimonious to just compare cropping vs perennial)
es_asp_rma2 <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ mat + map + label, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2,
                                           r2 = .x$R2))) # pull out I2 and R2

# Pull out dataframe with I2 and R2 values to see for which indicators does adding moderators to the model add significant explanatory power. Then we can pull these models out and look at them individually.
es_asp_rma_het2 <- es_asp_rma2 %>%
  select(indicator, het_df) %>%
  ungroup() %>%
  transmute(indicator,  het_df) %>%
  unnest(cols = c(het_df)) %>%
  arrange(desc(r2))
# this appears less useful - the highest R2 is .48 for TN% and the I2 for that is still 73%

# try another version of the model with MAt, MAP, and tillage only
es_asp_rma3 <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ mat + map + till, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2,
                                           r2 = .x$R2))) # pull out I2 and R2

# Pull out dataframe with I2 and R2 values to see for which indicators does adding moderators to the model add significant explanatory power. Then we can pull these models out and look at them individually.
es_asp_rma_het3 <- es_asp_rma3 %>%
  select(indicator, het_df) %>%
  ungroup() %>%
  transmute(indicator,  het_df) %>%
  unnest(cols = c(het_df)) %>%
  arrange(desc(r2))
# Also not very good - highest R2 is 0.46 for KSSL WSA, I2 is still 81%

# 9 - Sandbox of making graphs ----
# One question I had - do management decisions affect indicators differently in different climate types
# we can graph this! (we can also model it :))

# first - want to know which indicator has the most heterogeneity
es_asp_rma_nomod2 <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = proj_trt, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) %>%
  # make dataframe with effect sizes and variances for individual studies
  mutate(study_df = map(rma_obj, ~data.frame(es = .x$yi,
                                             se= sqrt(.x$vi),
                                             type = "study",
                                             study=.x$slab))) %>%
  # make dataframe with overall random effects model estimate and standard error
  mutate(summary_df = map(rma_obj, ~data.frame(es = .x$b,
                                               se=.x$se,
                                               type= "summary",
                                               study="summary"))) %>%
  # put two dfs together into one that can be used to make a forest plot
  mutate(plot_df = map2(study_df, summary_df, rbind)) %>%
  mutate(het_df = map(rma_obj, ~data.frame(i2 = .x$I2)))

es_asp_rma_nomod_het <- es_asp_rma_nomod2 %>%
  select(indicator, het_df) %>%
  ungroup() %>%
  transmute(indicator,  het_df) %>%
  unnest(cols = c(het_df)) %>%
  arrange(desc(i2))

# It's POX-C!

ggplot(es_asp %>% filter(indicator=="pox_c"), aes(x=climate, y=yi, colour=label)) +
  geom_boxplot()

ggplot(es_asp %>% filter(indicator=="pox_c") %>% filter(climate=="cool_dry"), aes(x=project, y=yi, colour=label)) +
  geom_boxplot()

# The tough thing here is that many climate categories are represented by one project only so it's really hard to tell 

ggplot(es_asp %>% filter(indicator=="pox_c"), aes(x=map, y=yi)) +
  geom_point(aes(colour=label)) +
  geom_smooth(method="lm", formula = y~x) +
  stat_cor(aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
           size=4) +
  theme_katy()
