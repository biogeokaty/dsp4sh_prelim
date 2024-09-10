# 07 - PCA and Significance of Soil Health Metrics/Indicators
# Katy Dynarski, March 2024

# Looking at the other soil health metrics...

# 0 - Import data ----
surf <- read.csv(here("data_processed","05_surface_horizons.csv"))
project <- read.csv(here("data_processed", "05_project_data.csv"))

# 1 - Boxplots of all indicators across soils and managements ----
# 1.1 - Boxplots for each indicator (x axis is project) ----
# can we automate with - list of indicators, plot boxplots across soil series x label, ANOVA and generate letters?
# what are the indicators? soc_pct, soc_stock_hrz, tn_pct:yoder_agg_stab_mwd, p_h:ace
surf_long <- surf %>%
  select(dsp_pedon_id, soil, project, label, lu, climate, 
         soc_pct, bulk_density, kssl_wsa:yoder_agg_stab_mwd, soil_respiration:acid_phosphatase, arylsulfatase:ace) %>%
  pivot_longer(soc_pct:ace, names_to="indicator", values_to="value")

# Make list of indicators for purrr
indicators <- surf_long %>%
  distinct(indicator) %>%
  pull(indicator) %>%
  as.character()

# test function
.x <- indicators[1]

surf_long %>% 
  filter(indicator == .x) %>%
  ggplot(aes(x = project, y = value, fill=label)) +
  geom_boxplot() +
  labs(x="Project", 
       y=glue::glue({filter(indicator_labs_df, indicator==.x)$label}), 
       title=glue::glue("0-10 cm Indicator Values - ", 
                        {filter(indicator_labs_df, indicator==.x)$label})) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# put into map function and iterate
plot_list <- map(.x = indicators,
                 .f = ~{
                   surf_long %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x = project, y = value, fill=label)) +
                     geom_boxplot() +
                     labs(x="Project", 
                          y=glue::glue({filter(indicator_labs_df, indicator==.x)$label}), 
                          title=glue::glue("0-10 cm Indicator Values - ", 
                                           {filter(indicator_labs_df, indicator==.x)$label})) +
                     scale_fill_viridis(discrete=TRUE, name="Management") +
                     theme_katy() +
                     theme(axis.text.x = element_text(angle = 45, hjust=1))
                   
                   ggsave(here("figs/indicator_boxplots", glue::glue("indicator_", .x, ".png")), 
                          width=10, height=7, units="in", dpi=400)
                 })

# 1.2 - Facet wrap of all indicators for each project ----
# Also want to make a version with a facet wrap of all indicators for each project - to see which indicators are measured in which projects
projects <- surf_long %>%
  distinct(project) %>%
  pull(project) %>%
  as.character()

plot_list2 <- map(.x = projects,
                 .f = ~{
                   surf_long %>% 
                     filter(project == .x) %>%
                     ggplot(aes(x=label, y=value, fill=label)) +
                     geom_boxplot() +
                     facet_wrap(vars(indicator), scales="free_y", labeller=labeller(indicator=indicator_labs)) +
                     scale_fill_viridis(discrete=TRUE, name="Management") +
                     labs(x="Management", 
                          title=glue::glue("0-10 cm Indicator Values - ", 
                                           {unique(filter(surf_long, project==.x)$project)})) +
                     theme_katy()
                   
                   ggsave(here("figs/indicator_boxplots_projects", glue::glue("project_", .x, ".png")), 
                          width=12, height=10, units="in", dpi=400)
                 })

# 1.3 - boxplot of two most common indicators with MAP ----
# Make a version that plots two most commonly measured indicators - BG and POXC - alongside MAP
# Calculate site mean climate data
site_clim_dist <- project %>%
  distinct(project, climate)

site_clim_sum <- project %>%
  group_by(project) %>%
  summarize(across(mat:map, ~ mean(.x, na.rm = TRUE))) %>%
  left_join(site_clim_dist, by="project")

# Join to soils data
surf_avg_clim <- surf %>%
  left_join(select(site_clim_sum, project, mat, map), by="project", suffix=c("_pedon", "_site"))

# Plot - make one line plot that shows precipitation across sites, then plot on top of boxplots
map_plot <- ggplot(surf_avg_clim, aes(x=fct_reorder(project, map_site), y=map_site)) +
  geom_point(color="dodgerblue4") +
  geom_line(aes(group=1), color="dodgerblue2") +
  labs(x="Project", y="MAP (mm)") +
  theme_katy_grid() +
  easy_remove_x_axis() 

poxc_box <- ggplot(surf_avg_clim, aes(x=fct_reorder(project, map_site), y=pox_c, fill=label)) +
  geom_boxplot() +
  labs(x="Project", y=expression("POX-C"~(mg~kg^-1))) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy_grid() +
  easy_remove_x_axis()

bg_box <- ggplot(surf_avg_clim, aes(x=fct_reorder(project, map_site), y=bglucosidase, fill=label)) +
  geom_boxplot() +
  labs(x="Project", y=expression(beta*"-glucosidase activity"~(mg~kg^-1~hr-1))) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy_grid() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ind_grid <- plot_grid(map_plot,
                      poxc_box + theme(legend.position="none"), 
                      bg_box + theme(legend.position="none"), 
                      ncol=1, rel_heights=c(.5, 1, 1.25), align="v")

ind_leg <- get_legend(bg_box + theme(legend.box.margin = margin(0, 0, 0, 12)))

ind_full <- plot_grid(ind_grid, ind_leg, rel_widths = c(3, .5))
ind_full
ggsave(here("figs", "indicator_patt_precip.png"), width=11, height=9, units="in", dpi=400)

# 1.4 - Table of mean indicator values in reference, SHM, and BAU systems ----
indicator_summary_wide <- surf_long %>%
  group_by(project, label, indicator) %>%
  summarize(across(value, mean_sd_cv)) %>%
  na.omit() %>%
  unite("summary", value_mean:value_sd)  %>%
  pivot_wider(names_from=indicator, values_from=summary:value_cv)

write_csv(indicator_summary_wide, here("figs", "indicator_summary_wide.csv"))

indicator_summary_ref <- indicator_summary_wide %>%
  filter(label=="Ref")
write_csv(indicator_summary_ref, here("figs", "indicator_summary_ref.csv"))

# Table of CV, calculated as the sd/mean across all sites/treatments
indicator_cv_range <- surf_long %>%
  group_by(indicator) %>%
  summarize(across(value, mean_sd_cv)) %>%
  na.omit()
write_csv(indicator_cv_range, here("figs", "indicator_cv.csv"))

# 1.5 - Table of mean indicator values by land use ----
indicator_summary_lu <- surf_long %>%
  group_by(lu, label, indicator) %>%
  summarize(across(value, mean_sd)) %>%
  na.omit() %>%
  unite("summary", value_mean:value_sd) %>%
  pivot_wider(names_from=indicator, values_from=summary)
write_csv(indicator_summary_lu, here("figs", "indicator_summary_lu.csv"))

# 2 - Correlation matrix to see structure of indicator data ----
indicators_only <- surf %>%
  select(bulk_density, soc_pct, kssl_wsa:yoder_agg_stab_mwd, soil_respiration:acid_phosphatase, arylsulfatase:ace)
indicators_normalized <- scale(indicators_only)

# Plot correlation matrix
corr_matrix <- cor(indicators_normalized, use="pairwise.complete.obs")
corr_pmat <- cor_pmat(indicators_normalized)
ggcorrplot(corr_matrix, p.mat=corr_pmat, hc.order=TRUE, type="lower", lab=TRUE, insig="blank") +
  scale_x_discrete(labels=indicator_labs) +
  scale_y_discrete(labels=indicator_labs) +
  theme_katy() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title=element_blank())
ggsave(here("figs", "surface_indicator_corrplot.png"), height=9, width=9, units="in", dpi=400)

# Make a smaller version of correlation matrix for presentation
indicators_red <- indicators_only %>%
  select(soc_pct, kssl_wsa, ace, pox_c, acid_phosphatase, bglucosidase, bglucosaminidase, arylsulfatase) %>%
  scale()

corr_matrix2 <- cor(indicators_red, use="pairwise.complete.obs")
corr_pmat2 <- cor_pmat(indicators_red)
ggcorrplot(corr_matrix2, p.mat=corr_pmat2, type="lower", lab=TRUE, insig="blank") +
  scale_x_discrete(labels=indicator_labs) +
  scale_y_discrete(labels=indicator_labs) +
  theme_katy() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title=element_blank())
ggsave(here("figs", "surface_indicator_corrplot_reduced.png"), height=6, width=6, units="in", dpi=400)

# 3 - PCA of indicators ----
# Prcomp will not work with missing values - Can impute the missing values with missMDA package!!
# prcomp is preferable to princomp because princomp uses a spectral decomposition approach (examines covariances between variables) and prcomp uses singular value decomposition, which examines covariances between individuals - supposed to have slightly better numerical accuracy
nb <- estim_ncpPCA(indicators_only,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
res.comp <- imputePCA(indicators_only, ncp = nb$ncp)
imp_ind <- res.comp$completeObs
imp_pca <- prcomp(imp_ind, scale.=TRUE)
summary(imp_pca)

# plot screeplot with fviz_eig() from factoextra 
fviz_eig(imp_pca)

# Plot loadings using factoextra
show_col(viridis_pal(option="plasma")(6)) # show hex codes for viridis color palette option plasma so I can pull them into a PCA plot

pca_loadings <- fviz_pca_var(imp_pca,
                             col.var = "contrib", # Color by contributions to the PC
                             gradient.cols = c("#FCA636FF", "#B12A90FF", "#0D0887FF" ),
                             repel = TRUE     # Avoid text overlapping
)
pca_loadings
ggsave(here("figs", "indicators_pca_loadings.png"), width=8, height=6, units="in")

# plot pca biplot with autoplot() function from ggfortify, colored by soil series
pca_soil <- autoplot(imp_pca, data=surf, colour="soil", size=3,
                     loadings.label = TRUE, loadings.label.size = 4,
                     loadings.label.colour = 'black', loadings.color="black", 
                     loadings.label.repel=TRUE) +
  scale_color_viridis(discrete=TRUE, name="Soil Series", option="magma") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_katy()
pca_soil
ggsave(here("figs", "indicators_pca_soil.png"), width=8, height=6, units="in")
    
# plot pca colored by management
pca_mgmt <- autoplot(imp_pca, data=surf, colour="label", size=3,
                     loadings.label = TRUE, loadings.label.size = 4,
                     loadings.label.colour = 'black', loadings.color="black", 
                     loadings.label.repel=TRUE) +
  scale_color_viridis(discrete=TRUE, name="Management") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_katy()
pca_mgmt
ggsave(here("figs", "indicators_pca_mgmt.png"), width=8, height=6, units="in")

# Put two together in a grid
plot_grid(pca_soil, pca_mgmt, ncol=1)
ggsave(here("figs", "indicator_pca_grid.png"), width=7, height=10, units="in", dpi=400)

# try a PCA grouped by climate
pca_clim <- autoplot(imp_pca, data=surf, colour="climate") +
  scale_color_viridis(discrete=TRUE, name="Climate") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_katy()
pca_clim
ggsave(here("figs", "indicators_pca_clim.png"), width=8, height=6, units="in")

# 4 - Redundancy analysis of indicators ?----


# 5 - ARCHIVE - Plot indicators that are the most sensitive to treatment (ANOVA) ----
proj_anova <- surf_long %>%
  na.omit() %>%
  group_by(project, indicator) %>% 
  nest() %>%
  mutate(anova_obj = map(data, ~anova(aov(value ~ label, data=.x)))) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(project, indicator, anova_tidy) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term=="label") %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant"))


ggplot(proj_anova, aes(x=project, y=sig, fill=sig)) +
  geom_bar(stat="identity") +
  facet_wrap(vars(indicator)) +
  theme_katy()
# think a little bit more on how best to present this data - but yes, which indicators are significant varies between projects

# maybe make a top 3?

proj_indicator_top3 <- proj_anova %>%
  group_by(project) %>%
  arrange(p.value, .by_group=TRUE) %>%
  slice_head(n=3) %>%
  select(project, indicator, p.value)
flextable(proj_indicator_top3)

# What I would like is to plot the number of sigs for each indicator and then arrange and plot
indicator_sig_count <- proj_anova %>%
  group_by(indicator, sig) %>%
  count() %>%
  pivot_wider(names_from = "sig", values_from="n") %>%
  select(indicator, significant) %>%
  mutate(significant = replace_na(significant, 0)) %>%
  mutate(half = ifelse(significant > 4, "more than half", "less than half"))

ggplot(filter(indicator_sig_count, indicator!="p_h"), aes(x=fct_reorder(indicator, significant, .desc=TRUE), y=significant)) +
  geom_bar(stat="identity", aes(fill=half)) +
  geom_vline(xintercept=6.5, linetype="dashed", linewidth=2) +
  labs(y="Significance count", x="Indicator") +
  theme_katy() +
  scale_fill_viridis(discrete=TRUE, begin=0, end=0.85) +
  scale_x_discrete(labels=indicator_labs) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none")
ggsave(here("figs", "indicator_sig_count.png"), height=6, width=8, units="in")

# 6 - ARCHIVE - Grouped analysis of indicator sensitivity in surface soils across all projects using ANOVA and climate groupings ----
# Attach climate data to horizon data
# ANOVA to test significance of relationship between management and indicators within each project (nested within climate)
ind_sig_clim <- surf_long %>%
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
ind_sig_clim$project <- factor(ind_sig_clim$project, 
                               levels=c("TexasA&MPt-2", "UnivOfMinnesota","WashingtonState","KansasState",
                                        "Illinois","OregonState", "UConn",
                                        "UTRGV", "TexasA&MPt-1",
                                        "NCState"), ordered=TRUE)

ggplot(ind_sig_clim, aes(x=project, y=indicator, fill=sig)) +
  geom_tile() +
  scale_fill_viridis(discrete=TRUE) +
  geom_rect(xmin=c(0), xmax=c(4.5), ymin=c(0), ymax=c(17.5), color="black", fill=NA )+
  geom_rect(xmin=c(4.5), xmax=c(7.5), ymin=c(0), ymax=c(17.5), color="black", fill=NA ) +
  geom_rect(xmin=c(7.5), xmax=c(9.5), ymin=c(0), ymax=c(17.5), color="black", fill=NA )+
  geom_rect(xmin=c(9.5), xmax=c(10.5), ymin=c(0), ymax=c(17.5), color="black", fill=NA ) +
  scale_y_discrete(labels=indicator_labs)
ggsave(here("figs", "indicators_climate.png"), width=8, height=6, units="in")

# This is potentially useful - I think this is a good start - but there are challenges with this analysis:
# not sure that indicator sensitivity is the most useful thing to look at. most indicators are sensitive to management in most projects, it seems that most indicators are sensitive to some extent
# variability mostly seems to exist within individual projects
# Next Task ... try a meta-analysis approach which will hopefully help crack into these differences!