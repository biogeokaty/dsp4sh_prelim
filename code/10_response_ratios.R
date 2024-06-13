# Using Meta-Analysis to Examine Soil Health Indicator Sensitivity to Management
# Katy Dynarski, March 2024

# 0 - Import data ----
meta_df <- read.csv(here("data_processed", "05_meta_df.csv"))
project <- read.csv(here("data_processed", "05_project_data.csv"))

mean_sd <- list(
  mean = ~round(mean(.x, na.rm = TRUE), 2), 
  sd = ~round(sd(.x, na.rm = TRUE), 2)
)

min_max <- list(
  min = ~min(.x, na.rm=TRUE), 
  max = ~max(.x, na.rm=TRUE)
)

# 1 - Calculate effect size for individual treatments and plot indicator effect sizes ----

# Pivot data longer
meta_long <- meta_df %>%
  select(project, dsp_pedon_id, soil, label, climate, lu, till, trt, soc_stock_0_30cm, soc_stock_100cm, 
         soc_pct, bulk_density, kssl_wsa:yoder_agg_stab_mwd, soil_respiration:acid_phosphatase, arylsulfatase:ace) %>%
  pivot_longer(soc_stock_0_30cm:ace, names_to="indicator", values_to="value")

# Make a new treatment variable with two values, BAU and ASP. ASP will contain all of the SHM and Ref soils.
meta_long2 <- meta_long %>%
  mutate(treat = ifelse(label=="BAU", "BAU", "ASP"))

# Calculate mean, sd, and n for BAU soils - need this for effect size calculation input
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

# Calculate mean, sd, and n for Ref and SHM
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

# Put data together. This will be input for escalc() function to calculate effect sizes
es_in_big <- nested2 %>%
  left_join(bau_means, by=c("indicator", "project")) %>%
  na.omit()
write_csv(es_in_big, here("figs", "mean_indicator_values.csv"))


# Calculate effect sizes
es_asp1 <- escalc(n1i = n_asp, n2i = n_bau, m1i = mean_asp, m2i = mean_bau, 
                  sd1i = sd_asp, sd2i = sd_bau, data = es_in_big, measure = "ROM")

# Calculate site mean climate data
site_clim_dist <- project %>%
  distinct(project, climate)

site_clim_sum <- project %>%
  group_by(project) %>%
  summarize(across(mat:map, ~ mean(.x, na.rm = TRUE))) %>%
  left_join(site_clim_dist, by="project")

# Calculate site mean clay%
site_clay <- meta_df %>%
  group_by(project, label) %>%
  summarize(clay_mean = mean(clay_tot_psa, na.rm=TRUE))

clay_min_max <- site_clay %>%
  ungroup() %>%
  na.omit() %>%
  summarize(across(where(is.numeric), min_max)) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

es_asp <- es_asp1 %>%
  na.omit() %>% # omit any rows with NA - this should just be places where NaNs were calculated
  unite("proj_trt", project, trt, remove = FALSE) %>%
  rename(till_orig = till) %>%
  mutate(till = ifelse(till_orig == "Till", "No-till", till_orig)) %>% # fix Texas A&M data that was mislabeled - no-till treatment was coded as Till
  relocate(till, .before=till_orig) %>%
  left_join(site_clim_sum, by="project") %>%
  left_join(site_clay, by=c("project", "label"))

# check by plotting
ggplot(es_asp %>% filter(indicator=="soc_pct"), aes(x=till, y=yi, colour=label)) +
  geom_boxplot()
# Now we have projects properly separated, with tillage and land use data as well as climate!

# 2 - Calculate random effects models for all projects with no moderator variables, make forest plots ----

# nest data by indicators and treatment contrasts, calculate random-effects models, extract data from model and make dataframe needed for plotting forest plot of each indicator and contrast
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
  mutate(plot_df = map2(study_df, summary_df, rbind))

# pull out plotting data
es_asp_rma_plot_df <- es_asp_rma_nomod %>%
  select(indicator, label, plot_df, rma_tidy) %>%
  ungroup() %>%
  transmute(indicator, label, plot_df, rma_tidy) %>%
  unnest(cols = c(plot_df, rma_tidy), names_sep = "_") %>%
  select(indicator, label, plot_df_es, plot_df_se, plot_df_type, plot_df_study, rma_tidy_p.value) %>%
  rename(es = plot_df_es,
         se = plot_df_se,
         type = plot_df_type,
         study = plot_df_study,
         summary_pval = rma_tidy_p.value)

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
                      ggplot(aes(x=factor(study, levels=rev(study)),y=es,ymax=es+se,ymin=es-se,size=factor(type),
                                 colour=factor(type))) + 
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

## 2.1 - make plot of summary effect sizes for all indicators ----
es_rma_summary <- es_asp_rma_plot_df %>%
  filter(type=="summary") %>%
  filter(indicator!="soc_stock_100cm") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  mutate(sig = ifelse(summary_pval<0.05, "significant", "not significant")) %>%
  select(-study)

ggplot(es_rma_summary, aes(x=fct_relevel(indicator, rev(indicator_plotting_order)), y=es, ymax=es+se, ymin=es-se, color=label)) + 
  geom_pointrange(size=1, linewidth=1) +
  geom_text(aes(label=ifelse(sig=="significant", "*", "")), 
            color="black", position=position_nudge(x=0.25), size=5) +
  coord_flip() + 
  geom_hline(yintercept=0, lty=2,linewidth=1) +
  labs(x="Indicator", y="Log response ratio") +
  scale_x_discrete(labels=indicator_labs) +
  scale_color_viridis(discrete=TRUE, name="Management") +
  theme_katy()
ggsave(here("figs", "indicator_effect_sizes.png"), width=8, height=5.5, units="in", dpi=400)

# Make table with n for interpretability
es_rma_n_df <- es_asp_rma_nomod %>%
  select(indicator, label, plot_df) %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  ungroup() %>%
  transmute(indicator, label, plot_df) %>%
  unnest(cols = c(plot_df)) %>%
  count(indicator, label, type) %>%
  filter(type=="study") %>%
  select(-type) %>%
  mutate(indicator=factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator, label)
flextable(es_rma_n_df)
write_csv(es_rma_n_df, here("figs", "es_rma_n.csv"))

## 2.2 - make table summarizing overall random-effects model ----
# tau2, I2, H2, Q-test, model estimates

rma_report_table <- es_asp_rma_nomod %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp))) %>%
  select(indicator, label, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, label, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  select(-term, -type) %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator, label)
flextable(rma_report_table)
write_csv(rma_report_table, here("figs", "rma_report_table.csv"))

## 2.3 - ARCHIVE - make plot of summary effect sizes for all indicators, grouped by land use ----
# this would help with identifying which indicators are sensitive to management in a variety of contexts and determine the most useful indicators
# I'm imagining a forest plot grouped by climate, faceted for each indicator
# might also be good to have as a table for me to reference

# Need to calculate another RMA with land use as a grouping factor
lu_rma <- es_asp %>%
  group_by(indicator, label, lu) %>%
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
  mutate(plot_df = map2(study_df, summary_df, rbind))

# table for reporting
lu_rma_report_table <- lu_rma %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp))) %>%
  select(indicator, label, lu, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, label, lu, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant")) %>%
  select(-term, -type) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3))
flextable(lu_rma_report_table)
write_csv(lu_rma_report_table, here("figs", "lu_rma_report_table.csv"))

# pull out plotting data
lu_rma_plot_df <- lu_rma %>%
  select(indicator, label, lu, plot_df, rma_tidy) %>%
  ungroup() %>%
  transmute(indicator, label, lu, plot_df, rma_tidy) %>%
  unnest(cols = c(plot_df, rma_tidy), names_sep = "_") %>%
  select(indicator, label, lu, plot_df_es, plot_df_se, plot_df_type, plot_df_study, rma_tidy_p.value) %>%
  rename(es = plot_df_es,
         se = plot_df_se,
         type = plot_df_type,
         study = plot_df_study,
         summary_pval = rma_tidy_p.value) %>%
  filter(type=="summary") %>%
  filter(indicator!="soc_stock_100cm") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  mutate(sig = ifelse(summary_pval<0.05, "significant", "not significant")) %>%
  select(-study)

# Make plot - land use on x axis, faceted by indicator
lu_rma_plot <- ggplot(lu_rma_plot_df, aes(x=lu, y=es, ymax=es+se, ymin=es-se, color=label)) + 
  geom_pointrange(size=1, linewidth=1) +
  geom_text(aes(label=ifelse(sig=="significant", "*", "")), 
            color="black", position=position_nudge(x=0.25), size=5) +
  coord_flip() + 
  geom_hline(yintercept=0, lty=2,linewidth=1) +
  labs(x="Land use", y="Log response ratio") +
  facet_wrap(vars(indicator), scales="free_x", labeller=labeller(indicator=indicator_labs)) +
  scale_x_discrete(labels=indicator_labs) +
  scale_color_viridis(discrete=TRUE, name="Management") +
  theme_katy()
lu_rma_plot
ggsave(here("figs", "indicator_effect_sizes_lu.png"), width=12.25, height=9.5, units="in", dpi=400)

# Make table with n for interpretability
lu_n_df <- lu_rma %>%
  select(indicator, label, lu, plot_df) %>%
  ungroup() %>%
  transmute(indicator, label, lu, plot_df) %>%
  unnest(cols = c(plot_df)) %>%
  count(indicator, label, lu, type) %>%
  filter(type=="study") %>%
  select(-type)
flextable(lu_n_df)
write_csv(lu_n_df, here("figs", "lu_n.csv"))

## 2.4 - ARCHIVE - make plot of summary effect sizes for all indicators, grouped by tillage ----

# Need to calculate another RMA with tillage as a grouping factor
till_rma <- es_asp %>%
  group_by(indicator, label, till) %>%
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
  mutate(plot_df = map2(study_df, summary_df, rbind))

# table for reporting
till_rma_report_table <- till_rma %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp))) %>%
  select(indicator, label, till, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, label, till, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant")) %>%
  select(-term, -type) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3))
flextable(till_rma_report_table)
write_csv(till_rma_report_table, here("figs", "till_rma_report_table.csv"))

# pull out plotting data
till_rma_plot_df <- till_rma %>%
  select(indicator, label, till, plot_df, rma_tidy) %>%
  ungroup() %>%
  transmute(indicator, label, till, plot_df, rma_tidy) %>%
  unnest(cols = c(plot_df, rma_tidy), names_sep = "_") %>%
  select(indicator, label, till, plot_df_es, plot_df_se, plot_df_type, plot_df_study, rma_tidy_p.value) %>%
  rename(es = plot_df_es,
         se = plot_df_se,
         type = plot_df_type,
         study = plot_df_study,
         summary_pval = rma_tidy_p.value) %>%
  filter(type=="summary") %>%
  filter(indicator!="soc_stock_100cm") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  mutate(sig = ifelse(summary_pval<0.05, "significant", "not significant")) %>%
  select(-study)

# Make plot - tillage on x axis, faceted by indicator
till_rma_plot <- ggplot(till_rma_plot_df, aes(x=till, y=es, ymax=es+se, ymin=es-se, color=label)) + 
  geom_pointrange(size=1, linewidth=1) +
  geom_text(aes(label=ifelse(sig=="significant", "*", "")), 
            color="black", position=position_nudge(x=0.25), size=5) +
  coord_flip() + 
  geom_hline(yintercept=0, lty=2,linewidth=1) +
  labs(x="Tillage", y="Log response ratio") +
  facet_wrap(vars(indicator), scales="free_x", labeller=labeller(indicator=indicator_labs)) +
  scale_x_discrete(labels=indicator_labs) +
  scale_color_viridis(discrete=TRUE, name="Management") +
  theme_katy()
till_rma_plot
ggsave(here("figs", "indicator_effect_sizes_till.png"), width=12.25, height=9.5, units="in", dpi=400)

# Make table with n for interpretability
till_n_df <- till_rma %>%
  select(indicator, label, till, plot_df) %>%
  ungroup() %>%
  transmute(indicator, label, till, plot_df) %>%
  unnest(cols = c(plot_df)) %>%
  count(indicator, label, till, type) %>%
  filter(type=="study") %>%
  select(-type)
flextable(till_n_df)
write_csv(till_n_df, here("figs", "till_n.csv"))

## 2.5 - ARCHIVE - make plot of summary effect sizes for all indicators, grouped by climate ----

# Need to calculate another RMA with climate as a grouping factor
clim_rma <- es_asp %>%
  group_by(indicator, label, climate) %>%
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
  mutate(plot_df = map2(study_df, summary_df, rbind))

# pull out plotting data
clim_rma_plot_df <- clim_rma %>%
  select(indicator, label, climate, plot_df, rma_tidy) %>%
  ungroup() %>%
  transmute(indicator, label, climate, plot_df, rma_tidy) %>%
  unnest(cols = c(plot_df, rma_tidy), names_sep = "_") %>%
  select(indicator, label, climate, plot_df_es, plot_df_se, plot_df_type, plot_df_study, rma_tidy_p.value) %>%
  rename(es = plot_df_es,
         se = plot_df_se,
         type = plot_df_type,
         study = plot_df_study,
         summary_pval = rma_tidy_p.value) %>%
  filter(type=="summary") %>%
  filter(indicator!="soc_stock_100cm") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  mutate(sig = ifelse(summary_pval<0.05, "significant", "not significant")) %>%
  select(-study)

# Make plot - climate on x axis, faceted by indicator
clim_rma_plot <- ggplot(clim_rma_plot_df, aes(x=climate, y=es, ymax=es+se, ymin=es-se, color=label)) + 
  geom_pointrange(size=1, linewidth=1) +
  geom_text(aes(label=ifelse(sig=="significant", "*", "")), 
            color="black", position=position_nudge(x=0.25), size=5) +
  coord_flip() + 
  geom_hline(yintercept=0, lty=2,linewidth=1) +
  labs(x="Climate", y="Log response ratio") +
  facet_wrap(vars(indicator), scales="free_x", labeller=labeller(indicator=indicator_labs)) +
  scale_x_discrete(labels=indicator_labs) +
  scale_color_viridis(discrete=TRUE, name="Management") +
  theme_katy()
clim_rma_plot
ggsave(here("figs", "indicator_effect_sizes_clim.png"), width=12.25, height=9.5, units="in", dpi=400)

# Make table with n for interpretability
clim_n_df <- clim_rma %>%
  select(indicator, label, climate, plot_df) %>%
  ungroup() %>%
  transmute(indicator, label, climate, plot_df) %>%
  unnest(cols = c(plot_df)) %>%
  count(indicator, label, climate, type) %>%
  filter(type=="study") %>%
  select(-type)
flextable(clim_n_df)
write_csv(clim_n_df, here("figs", "clim_n.csv"))

# 3 - Calculate random-effects models with moderating variables ----

## 3.1 - Random-effects model with MAT as a moderator variable ----
# Calculate random-effects model
es_asp_rma_mat <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ mat, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy))

# Make table for reporting, also include QM test of moderators
mat_rma_report_table <- es_asp_rma_mat %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  filter(term != "intercept") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  select(-term, -type) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant"),
         estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator)
flextable(mat_rma_report_table)
write_csv(mat_rma_report_table, here("figs", "mat_rma_report_table_full.csv"))

# Extract significant moderator variables
mat_rma_report_table_sig <- mat_rma_report_table %>%
  filter(q_mod_sig=="significant") 
flextable(mat_rma_report_table_sig)
# Save CSV of significant moderator variables
write_csv(mat_rma_report_table_sig, here("figs", "mat_rma_report_table_sig.csv"))

# Plot significant moderator variables
# MAT - Significant for ACE, acid phosphatase, arylsulfatase, bglucosaminidase, bglucosidase, SOC%

# Make vector of indicators 
mat_indicators <- c("ace", "arylsulfatase", "bglucosaminidase", "bglucosidase", "soc_pct")

# Plot with map function
mat_plots <- map(.x = mat_indicators,
                 .f = ~{
                   es_asp %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x=mat, y=yi)) +
                     geom_point(aes(colour=label)) +
                     geom_smooth(method="lm", formula = y~x, color="black") +
                     geom_hline(yintercept=0, linetype="dashed") +
                     stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
                              size=4) +
                     scale_colour_viridis(discrete=TRUE, name="Management") +
                     labs(x=expression("Mean annual temperature"~(degree*C)), y="Log response ratio",
                          title=glue::glue({filter(indicator_labs_df, indicator==.x)$label})) + 
                     theme_katy_grid()
                 })

# Now make into a panel grid with cowplot
mat_ace <- pluck(mat_plots, 1)
mat_aryl <- pluck(mat_plots, 2)
mat_bgm <- pluck(mat_plots, 3)
mat_bg <- pluck(mat_plots, 4)
mat_soc <- pluck(mat_plots, 5)

mat_grid <- plot_grid(mat_ace + theme(legend.position="none", axis.title.x=element_blank()),
                      mat_aryl + theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank()),
                      mat_bgm + theme(legend.position="none", axis.title.y=element_blank(), axis.title.x=element_blank()),
                      mat_bg + theme(legend.position="none", axis.title.x=element_blank()),
                      mat_soc + theme(legend.position="none", axis.title.y=element_blank()),
                      align = 'vh',
                      labels = c("A", "B", "C", "D", "E"),
                      hjust = -1,
                      nrow = 2)

mat_leg <- get_legend(mat_ace)

plot_grid(mat_grid, mat_leg, rel_widths = c(3, .7))
ggsave(here("figs", "mat_mod_reg.png"), width=11, height=8, units="in")

## 3.2 - Random-effects model with MAP as a moderator variable ----
# Calculate random-effects model
es_asp_rma_map <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ map, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy))

# Make table for reporting, also include QM test of moderators
map_rma_report_table <- es_asp_rma_map %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  filter(term != "intercept") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  select(-term, -type) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant"),
         estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator)
flextable(map_rma_report_table)
write_csv(map_rma_report_table, here("figs", "map_rma_report_table_full.csv"))

# Extract significant moderator variables
map_rma_report_table_sig <- map_rma_report_table %>%
  filter(q_mod_sig=="significant") 
flextable(map_rma_report_table_sig)
# Save CSV of significant moderator variables
write_csv(map_rma_report_table_sig, here("figs", "map_rma_report_table_sig.csv"))

# Plot significant moderator variables
# MAP is significant for alkaline phosphatase, arylsulfatase, POX-C, and SOC stock 0-30cm
# Make vector of indicators 
map_indicators <- c("alkaline_phosphatase", "arylsulfatase", "pox_c")

# Plot with map function
map_plots <- map(.x = map_indicators,
                 .f = ~{
                   es_asp %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x=map, y=yi)) +
                     geom_point(aes(colour=label)) +
                     geom_smooth(method="lm", formula = y~x, color="black") +
                     geom_hline(yintercept=0, linetype="dashed") +
                     stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
                              size=4) +
                     scale_colour_viridis(discrete=TRUE, name="Management") +
                     labs(x=expression("Mean annual precipitation"~(mm~yr^-1)), y="Log response ratio",
                          title=glue::glue({filter(indicator_labs_df, indicator==.x)$label})) + 
                     theme_katy_grid()
                 })

# Now make into a panel grid with cowplot
map_alk <- pluck(map_plots, 1)
map_aryl <- pluck(map_plots, 2)
map_poxc <- pluck(map_plots, 3)

map_grid <- plot_grid(map_alk + theme(legend.position="none", axis.title.x=element_blank()),
                      map_aryl + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()),
                      map_poxc + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()),
                      align = 'vh',
                      labels = c("A", "B", "C"),
                      hjust = -1,
                      nrow = 1)

map_grid_label <- ggdraw(map_grid) +
  draw_label(expression("Mean annual precipitation"~(mm~yr^-1)), color="black", vjust=0, y=0)

map_leg <- get_legend(map_aryl)

plot_grid(map_grid_label, map_leg, rel_widths = c(3, .7))
ggsave(here("figs", "map_mod_reg.png"), width=11.5, height=4, units="in")

## 3.3 - Random-effects model with land use as a moderator variable ----
# Calculate random-effects model
es_asp_rma_lu <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ factor(lu), data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy))

# Make table for reporting, also include QM test of moderators
lu_rma_report_table <- es_asp_rma_lu %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  filter(term != "intercept") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  select(-type) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant"),
         estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator)
flextable(lu_rma_report_table)
write_csv(lu_rma_report_table, here("figs", "lu_rma_report_table_full.csv"))

# Extract significant moderator variables
lu_rma_report_table_sig <- lu_rma_report_table %>%
  filter(q_mod_sig=="significant") 
flextable(lu_rma_report_table_sig)
# Save CSV of significant moderator variables
write_csv(lu_rma_report_table_sig, here("figs", "lu_rma_report_table_sig.csv"))

# Plot significant moderator variables
lu_indicators_df <- lu_rma_report_table_sig %>%
  distinct(indicator)
lu_indicators <- as.vector(lu_indicators_df$indicator)

lu_plots <- map(.x = lu_indicators,
                .f = ~{
                  es_asp %>% 
                    filter(indicator == .x) %>%
                    ggplot(aes(x=lu, y=yi, fill=label)) +
                    geom_boxplot() +
                    geom_hline(yintercept=0, linetype="dashed") +
                    labs(x="Land use", y="Log response ratio",
                         title=glue::glue({filter(indicator_labs_df, indicator==.x)$label})) + 
                    scale_fill_viridis(discrete=TRUE, name="Management") +
                    theme_katy_grid()
                })

# Now make into a panel grid with cowplot
lu_bd <- pluck(lu_plots, 1)
lu_kssl <- pluck(lu_plots, 2)
lu_resp <- pluck(lu_plots, 3)
lu_alk <- pluck(lu_plots, 4)

lu_grid <- plot_grid(lu_bd + theme(legend.position="none", 
                                    axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     lu_kssl + theme(legend.position="none", axis.title.x=element_blank(),
                                    axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1)), 
                     lu_resp + theme(legend.position="none",  
                                   axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     lu_alk + theme(legend.position="none", axis.title.y=element_blank(),
                                   axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     align = 'vh',
                     labels = c("A", "B", "C", "D"),
                     hjust = -1,
                     nrow = 2, ncol = 2)

lu_grid_label <- ggdraw(lu_grid) +
  draw_label("Land use", color="black", vjust=0, y=0)

lu_leg <- get_legend(lu_alk + theme(legend.box.margin = margin(0, 0, 0, 12)))

plot_grid(lu_grid_label, lu_leg, rel_widths = c(3, .7))

ggsave(here("figs", "lu_mod_box.png"), width=13, height=8, units="in")

## 3.4 - Random-effects model with tillage as a moderator variable ----
es_asp_rma_till <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ factor(till), data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy))

# Make table for reporting, also include QM test of moderators
till_rma_report_table <- es_asp_rma_till %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  filter(term != "intercept") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  select(-type) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant"),
         estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator)
flextable(till_rma_report_table)
write_csv(till_rma_report_table, here("figs", "till_rma_report_table_full.csv"))

# Extract significant moderator variables
till_rma_report_table_sig <- till_rma_report_table %>%
  filter(q_mod_sig=="significant") 
flextable(till_rma_report_table_sig)
# Save CSV of significant moderator variables
write_csv(till_rma_report_table_sig, here("figs", "till_rma_report_table_sig.csv"))

# Plot significant moderator variables
# Tillage is significant for bulk density, KSSL WSA, Yoder WSA, and alkaline phosphatase
till_indicators <- c("bulk_density", "kssl_wsa", "yoder_agg_stab_mwd", "alkaline_phosphatase")

till_plots <- map(.x = till_indicators,
                .f = ~{
                  es_asp %>% 
                    filter(indicator == .x) %>%
                    ggplot(aes(x=till, y=yi, fill=label)) +
                    geom_boxplot() +
                    geom_hline(yintercept=0, linetype="dashed") +
                    labs(x="Tillage", y="Log response ratio",
                         title=glue::glue({filter(indicator_labs_df, indicator==.x)$label})) + 
                    scale_fill_viridis(discrete=TRUE, name="Management") +
                    theme_katy_grid()
                })

till_bd <- pluck(till_plots, 1)
till_kssl <- pluck(till_plots, 2)
till_yod <- pluck(till_plots, 3)
till_alk <- pluck(till_plots, 4)

till_grid <- plot_grid(till_bd + theme(legend.position="none", 
                                   axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     till_kssl + theme(legend.position="none", axis.title.x=element_blank(),
                                     axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1)), 
                     till_yod + theme(legend.position="none", 
                                      axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     till_alk + theme(legend.position="none", axis.title.x=element_blank(),
                                      axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1)),
                     align = 'vh',
                     labels = c("A", "B", "C", "D"),
                     hjust = -1,
                     nrow = 2)

till_grid_label <- ggdraw(till_grid) +
  draw_label("Tillage", color="black", vjust=0, y=0)

till_leg <- get_legend(till_alk + theme(legend.box.margin = margin(0, 0, 0, 12)))

plot_grid(till_grid_label, till_leg, rel_widths = c(3, .7))

ggsave(here("figs", "till_mod_box.png"), width=10.5, height=8, units="in")

## 3.5 - Random-effects model with clay% as a moderator variable ----
# Calculate random-effects model
es_asp_rma_clay <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ clay_mean, data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy))
# 7 studies are omitted from model fitting because of missing clay content data, I wouldn't really trust this...

# Make table for reporting, also include QM test of moderators
clay_rma_report_table <- es_asp_rma_clay %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  filter(term != "intercept") %>%
  filter(indicator!="soc_stock_0_30cm") %>%
  filter(indicator!="soc_stock_100cm") %>%
  select(-term, -type) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant"),
         q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant"),
         estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2),
         indicator = factor(indicator, levels=indicator_plotting_order)) %>%
  arrange(indicator)
flextable(clay_rma_report_table)
write_csv(clay_rma_report_table, here("figs", "clay_rma_report_table_full.csv"))

# 4 - Evaluate indicator sensitivity and heterogeneity ----
# want to rank indicators by heterogeneity?

het_ranking <- rma_report_table %>%
  filter(p.value<0.05) %>%
  arrange(i2)
flextable(het_ranking)
# Indicators with the lowest I2 values are alkaline phosphatase, KSSL WSA, SOC...but in general I2 is really high


# 4 - ARCHIVE - Random-effects model with MAT, MAP, land use, and tillage as moderator variables in one model ----

# One version of the model - group by indicators only (so Ref and SHM are tested in the same model), and include MAT, MAP, land use, and tillage as moderators
es_asp_rma_all_mod <- es_asp %>%
  group_by(indicator) %>%
  nest() %>%
  mutate(rma_obj = map(data, ~rma(yi, vi, slab = project, mod = ~ mat + map + factor(lu) + factor(till), data=.x, method="REML"))) %>%
  mutate(rma_tidy = map(rma_obj, broom::tidy)) 

# Make table for reporting, also include QM test of moderators
all_mod_rma_report_table <- es_asp_rma %>%
  mutate(report_df = map(rma_obj, ~data.frame(tau2 = .x$tau2,
                                              i2 = .x$I2,
                                              h2 = .x$H2,
                                              r2 = .x$R2,
                                              q_stat = .x$QE,
                                              q_pval = .x$QEp,
                                              q_mod = .x$QM,
                                              q_mod_pval = .x$QMp))) %>%
  select(indicator, rma_tidy, report_df) %>%
  ungroup() %>%
  transmute(indicator, rma_tidy, report_df) %>%
  unnest(cols = c(rma_tidy, report_df)) %>%
  mutate(q_sig = ifelse(q_pval<0.05, "significant", "not significant")) %>%
  mutate(q_mod_sig = ifelse(q_mod_pval<0.05, "significant", "not significant")) %>%
  filter(term != "intercept") %>%
  select(-type) %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         tau2 = round(tau2, 3),
         i2 = round(i2, 2),
         h2 = round(h2, 2),
         r2 = round(r2, 2),
         q_stat = round(q_stat, 1),
         q_pval = round(q_pval, 3),
         q_mod = round(q_mod, 1),
         q_mod_pval = round(q_mod_pval, 2))
flextable(all_mod_rma_report_table)
write_csv(all_mod_rma_report_table, here("figs", "all_mod_rma_report_table_full.csv"))

# Extract significant moderator variables
all_mod_rma_report_table_sig <- all_mod_rma_report_table %>%
  filter(q_mod_sig=="significant") 
flextable(all_mod_rma_report_table_sig)
# Save CSV of significant moderator variables
write_csv(all_mod_rma_report_table_sig, here("figs", "all_mod_rma_report_table_sig.csv"))

# 5 - ARCHIVE - Alternate versions of random-effects models with moderating variables with less explanatory power ----
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
flextable(es_asp_rma_het2)
# this appears less useful - the highest R2 is .59 for arylsulfatase and the I2 for that is still 80%

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
flextable(es_asp_rma_het3)
# Also not very good - highest R2 is 0.52 for arysulfatase, I2 is still 83%

# Overall lesson from all of these different models is that adding moderator variables tends to not explain a ton of indicator response to management

# 6 - ARCHIVE - Sandbox of making graphs ----
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
flextable(es_asp_rma_nomod_het)
# It's POX-C!

# Plot POX-C effect sizes in different climates
ggplot(es_asp %>% filter(indicator=="pox_c"), aes(x=climate, y=yi, colour=label)) +
  geom_boxplot()

# Plot POX-C effect sizes by project
ggplot(es_asp %>% filter(indicator=="pox_c") %>% filter(climate=="cool_dry"), aes(x=project, y=yi, colour=label)) +
  geom_boxplot()

# The tough thing here is that many climate categories are represented by one project only so it's really hard to tell 

# Try plotting POX-C effect size against precip as a continuous variable
ggplot(es_asp %>% filter(indicator=="pox_c"), aes(x=map, y=yi)) +
  geom_point(aes(colour=label)) +
  geom_smooth(method="lm", formula = y~x) +
  stat_cor(aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
           size=4) +
  theme_katy()
# not sure what this tells us....

# 7 - ARCHIVE - Calculating effect sizes for indicators - first approach ----
# This is a slightly less good approach because it puts Ref and SHM in different columns and doesn't pull specific treatment inforamtion, it's just a lot less useful
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

# 8 - ARCHIVE - Calculate random-effects models for all indicators and plot ----
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

# Reshape dataframe to allow for forest plotting
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

# Make table of most significant indicators 
# Extract summary data from random-effects model, label contrasts as either significant or not significant
es_rma_sig <- es_rma %>%
  ungroup() %>%
  transmute(indicator, label, rma_tidy) %>%
  unnest(cols = c(rma_tidy)) %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant")) %>%
  arrange(desc(estimate))
flextable(es_rma_sig)
# Largest response ratio estimates are for:  arylsulfatase, aggregate stability (both methods), ace, bglucosaminidase, bglucosidase. Interesting that so many enzymes are fairly sensitive!

# Big indicators seem to be enzyme activity and aggregate stability

# 9 - ARCHIVE - Add moderator variables to dataframe and perform meta-regression ----

# Eventually I want to incorporate land-use variables in this - but they are pretty complex because there are typically multiple land uses and management strategies in each project. Will need to calculate even more comparisons!
# for now - incorporate mean MAT and MAP at site, also climate category - can see if climate grouping adds useful information 

# input data for random effects model is : es_long
# need: mean MAP and MAT data, climate grouping for each project

es_long_mod <- es_long %>% 
  ungroup() %>%
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
# I^2 went from 87% to 84.78% - not a big difference
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
summary(es_rma_mod %>% pluck("rma_obj", 28)) # model with climate as a category
summary(es_rma_mod2 %>% pluck("rma_obj", 28)) # model with MAT and MAP as continuous vars, accounts for 4% of heterogeneity - not much

# need to pull out I2 and R2 for all models and put into some kind of table to identify for which indicators is the response to management moderated by climate

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

# Look at table
flextable(es_rma_comp_df)

# Results are super interesting - response ratio variability for certain indicators is really well-explained by climate. For example, including climate as either a continuous or categorical variable explains 100% of the heterogeneity in TN% response.   
# I think this usually reflects some kind of problem in the model

ggplot(filter(es_long_mod, indicator=="arylsulfatase", label=="shm"), aes(x=climate, y=yi)) +
  geom_boxplot() +
  theme_katy()

# Looking at some plots, this is a little tricky. Sometimes when the moderating variable explains a lot, it's because there's not many data points and there's only one project in each climate category. So climate looks like it explains all of the heterogeneity, but it's because each project is assigned to a different climate. For that reason, I think it's probably not a good idea to use climate as a categorical variable - better to stick with continuous variables. 