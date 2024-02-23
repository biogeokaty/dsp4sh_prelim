# 06 - Soil Health Metrics 
# Katy Dynarski, November 2023

# Looking at the other soil health metrics...

# 0 - Import data ----
soc_pedon100 <- read.csv(here("data_processed", "05_soc_pedon_100cm.csv"))
soc_horizon_filt <- read.csv(here("data_processed", "05_soc_horizon_filtered.csv"))
soc_horizon_all <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# 1 - Plot other indicators vs SOC concentrations/stocks ####
# Pull out surface horizons only and just look around
surf <- soc_horizon_filt %>%
  filter(hrzdep_b == "5")
write_csv(surf, here("data_processed", "06_surface_horizons.csv"))

ggplot(surf, aes(x=label, y=soc_pct)) +
  geom_boxplot() +
  facet_wrap(~soil, scales="free_y")

ggplot(surf, aes(x=soil_respiration, y=soc_stock_hrz, color=label)) +
  geom_point() +
  facet_wrap(~soil, scales="free")

# 2 - Boxplots of all indicators across soils and managements ----
# can we automate with - list of indicators, plot boxplots across soil series x label, ANOVA and generate letters?
# what are the indicators? soc_pct, soc_stock_hrz, tn_pct:yoder_agg_stab_mwd, p_h:ace

surf_indicator <- surf %>%
  select(dsp_pedon_id, soil, label, soc_stock_hrz, soc_pct, bulk_density, tn_pct:yoder_agg_stab_mwd, p_h:ace) %>%
  pivot_longer(soc_stock_hrz:ace, names_to="indicator", values_to="value")

# Make list of indicators for purrr
indicators <- surf_indicator %>%
  distinct(indicator) %>%
  pull(indicator) %>%
  as.character()

# test function
.x <- indicators[1]

surf_indicator %>% 
  filter(indicator == .x) %>%
  ggplot(aes(x = soil, y = value, fill=label)) +
  geom_boxplot() +
  labs(x="Soil", 
       y=glue::glue({unique(filter(surf_indicator, indicator==.x)$indicator)}), 
       title=glue::glue("0-5 cm Indicator Values Across Soil and Management - ", 
                        {unique(filter(surf_indicator, indicator==.x)$indicator)})) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy()
ggsave(here("figs/indicator_boxplots", "test.png"), 
       width=10, height=7, units="in", dpi=400)

# put into map function and iterate
plot_list <- map(.x = indicators,
                 .f = ~{
                   surf_indicator %>% 
                     filter(indicator == .x) %>%
                     ggplot(aes(x = soil, y = value, fill=label)) +
                     geom_boxplot() +
                     labs(x="Soil", 
                          y=glue::glue({unique(filter(surf_indicator, indicator==.x)$indicator)}), 
                          title=glue::glue("0-5 cm Indicator Values Across Soil and Management - ", 
                                           {unique(filter(surf_indicator, indicator==.x)$indicator)})) +
                     scale_fill_viridis(discrete=TRUE, name="Management") +
                     theme_katy()
                   
                   ggsave(here("figs/indicator_boxplots", glue::glue("indicator_", .x, ".png")), 
                          width=10, height=7, units="in", dpi=400)
                 })

# 3 - Mixed models to look at which indicators are significantly influenced by management ----

# test out on POX-C - but should build a function and run this in purrr
poxc_mixed <- lmer(pox_c ~ label + (1|soil) + (1|hrzdep_b), data = soc_horizon_all)
summary(poxc_mixed)
# Soil type explains a little less than half the variance in the random effects, depth explains another third
# Label appears significant 

# Plot predictions of surface soil POX-C with treatment
pred_poxc <- ggpredict(poxc_mixed, terms = c("label"))
plot(pred_poxc)
# Label appears significant - test with LR

# Test by comparing full and reduced models, use likelihood ratio test
poxc_mixed_full <- lmer(pox_c ~ label + (1|soil) + (1|hrzdep_b), data = soc_horizon_all, REML = FALSE)
poxc_mixed_reduced <- lmer(pox_c ~ (1|soil)+ (1|hrzdep_b), data = soc_horizon_all, REML = FALSE)
anova(poxc_mixed_full, poxc_mixed_reduced, text="Chisq")
# Models are significantly different - means that label is a significant factor after controlling for soil type

# Test for significantly different groups
poxc_mixed_tukey <- glht(poxc_mixed, linfct = mcp(label = 'Tukey'))
summary(poxc_mixed_tukey)
# Each group is significantly different

# Run this iteratively for all indicators
lmm_sig <- surf_indicator %>% 
  group_by(indicator) %>% 
  nest() %>%
  mutate(anova_obj = map(data, ~{
    lmer_full <- lmer(value ~ label + (1|soil), data = ., REML=FALSE)
    lmer_reduced <- lmer(value ~ (1|soil), data = ., REML=FALSE)
    anova(lmer_full, lmer_reduced, text="Chisq")
    })) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(indicator, anova_tidy) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term == "lmer_full") %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) %>%
  select(indicator, deviance:p.value)

# 4 - Calculate response ratios of indicators ----
# haven't done this yet
response_ratios <- surf_indicator %>% 
  group_by(indicator) %>% 
  nest() %>%
  mutate(anova_obj = map(data, ~{
    lmer_full <- lmer(value ~ label + (1|soil), data = ., REML=FALSE)
    lmer_reduced <- lmer(value ~ (1|soil), data = ., REML=FALSE)
    anova(lmer_full, lmer_reduced, text="Chisq")
  })) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(indicator, anova_tidy) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term == "lmer_full") %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) %>%
  select(indicator, deviance:p.value)

# 5 - Correlation matrix to see structure of indicator data ----
indicators_only <- soc_horizon_filt %>%
  select(bulk_density, soc_pct, tn_pct:yoder_agg_stab_mwd, p_h:ace, clay_pct_field)
indicators_normalized <- scale(indicators_only)

# Make vector of indicator labels so they will print nicely
indicator_labs <- c("soc_pct" = "SOC %",
  "bglucosaminidase" = "B-glucosaminidase",
  "ace" = "ACE",
  "bglucosidase" = "B-glucosidase",
  "kssl_wsa" = "Aggregate stability (KSSL)",
  "yoder_agg_stab_mwd" = "Aggregate stability (MWD)",
  "arylsulfatase" = "Arylsulfatase",
  "pox_c" = "POX-C",
  "soc_stock_hrz" = "SOC stock",
  "tn_pct" = "Total N%",
  "bulk_density" = "Bulk density",
  "soil_respiration" = "Respiration",
  "phosphodiesterase" = "Phosphodiesterase",
  "alkaline_phosphatase" = "Alkaline phosphatase",
  "acid_phosphatase" = "Acid phosphatase",
  "p_h" = "pH")

# Plot correlation matrix
corr_matrix <- cor(indicators_normalized, use="pairwise.complete.obs")
corr_pmat <- cor_pmat(indicators_normalized)
ggcorrplot(corr_matrix, p.mat=corr_pmat, hc.order=TRUE, type="lower", lab=TRUE, insig="blank") +
  scale_x_discrete(labels=indicator_labs) +
  scale_y_discrete(labels=indicator_labs)
ggsave(here("figs", "surface_indicator_corrplot.png"), height=9, width=9, units="in")

# Another option to use corrplot - I like the look of the figure a bit better but unfortunately not compatible with ggsave
png(height=8, width=8, units="in", res=600, file=here("figs", "surface_indicator_corrplot.png"), type = "cairo")
corrplot(corr_matrix, order="FPC", type="lower", diag=FALSE,
         tl.col = 'black', tl.srt = 45, addCoef.col = 'black') 

dev.off() # corrplot doesn't generate ggplot2 compatible objects

# ACE is the strongest covariate with SOC

# Make a smaller version of correlation matrix for presentation
indicators_red <- indicators_only %>%
  select(soc_pct, tn_pct, kssl_wsa, ace, pox_c, acid_phosphatase, phosphodiesterase, bglucosidase) %>%
  scale()

corr_matrix2 <- cor(indicators_red, use="pairwise.complete.obs")
corr_pmat2 <- cor_pmat(indicators_red)
ggcorrplot(corr_matrix2, p.mat=corr_pmat2, type="lower", lab=TRUE, insig="blank") +
  scale_x_discrete(labels=indicator_labs) +
  scale_y_discrete(labels=indicator_labs)
ggsave(here("figs", "surface_indicator_corrplot_reduced.png"), height=6, width=6, units="in")

# 6 - Plot indicators that are the most sensitive to treatment ----
# How i think this might work - do same thing as before but with a map2 to iterate across indicators and projects 
# start with project and then see if there are common threads between projects??

soil_anova <- surf_indicator %>%
  na.omit() %>%
  group_by(soil, indicator) %>% 
  nest() %>%
  mutate(anova_obj = map(data, ~anova(aov(value ~ label, data=.x)))) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(soil, indicator, anova_tidy) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term=="label") %>%
  mutate(sig = ifelse(p.value<0.05, "significant", "not significant"))
  

ggplot(soil_anova, aes(x=soil, y=sig, fill=sig)) +
  geom_bar(stat="identity") +
  facet_wrap(vars(indicator)) +
  theme_katy()
# think a little bit more on how best to present this data - but yes, which indicators are significant varies between soil series

# maybe make a top 3?

soil_indicator_top3 <- soil_anova %>%
  group_by(soil) %>%
  arrange(p.value, .by_group=TRUE) %>%
  slice_head(n=3) %>%
  select(soil, indicator, p.value)
flextable(soil_indicator_top3)

# What I would like is to plot the number of sigs for each indicator and then arrange and plot
indicator_sig_count <- soil_anova %>%
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

# 6 - Plot relationship between SOC and indicators ----
# Mixed model for relationship between SOC and POX-C
poxc_soc_lmer <- lmer(pox_c ~ soc_pct + (1|label) + (1|soil), data = surf)
summary(poxc_soc_lmer)

pred_poxc_soc <- ggpredict(poxc_soc_lmer, terms=c("soc_pct"))
plot(pred_poxc_soc)

ggplot(surf, aes(x=soc_pct, y=pox_c)) +
  labs(x="SOC Concentration (%)", y="POX-C") +
  scale_color_viridis(discrete=TRUE, name="Soil Series") +
  theme_katy()

ggplot(surf, aes(x=soc_pct, y=kssl_wsa)) +
  geom_point(aes(colour=soil, shape=label)) +
  geom_smooth(method=lm, color="black") +
  ylim(c(0,100)) +
  labs(x="SOC Concentration (%)", y="Water-stable aggregates") +
  scale_color_viridis(discrete=TRUE, name="Soil Series") +
  scale_shape_discrete(name="Management") +
  theme_katy()

# 7 - PCA of indicators ----
# Try PCA
indicator_pca <- princomp(corr_matrix) # can't run PCA because there are some comparisons that don't have data
summary(indicator_pca)

# Hard to do this because of missing values...

# Can impute the missing values with missMDA package!!
# try with already scaled data and princomp
# after testing - this is not my preferred method, skip to next section
nb <- estim_ncpPCA(indicators_normalized,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
nb$ncp
res.comp <- imputePCA(indicators_normalized, ncp = nb$ncp)
imp_corr <- cor(res.comp$completeObs)
imp_pca <- princomp(imp_corr)
summary(imp_pca)
imp_pca$loadings[,1:2]

# also try with original data and prcomp, plot with ggfortify
# prcomp is preferable to princomp because princomp uses a spectral decomposition approach (examines covariances between variables) and prcomp uses singular value decomposition, which examines covariances between individuals - supposed to have slightly better numerical accuracy
nb2 <- estim_ncpPCA(indicators_only,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
nb2$ncp
res.comp2 <- imputePCA(indicators_only, ncp = nb2$ncp)
imp_ind <- res.comp2$completeObs
imp_pca2 <- prcomp(imp_ind, scale.=TRUE)
imp_pca3 <- prcomp(res.comp$completeObs) # also run a test with the imputed scaled data. I think it makes more sense to impute the actual data first, THEN scale during the PCA 
summary(imp_pca2)

# plot screeplot with fviz_eig() from factoextra 
fviz_eig(imp_pca2)

# Plot loadings using factoextra
show_col(viridis_pal(option="plasma")(6)) # show hex codes for viridis color palette option plasma

pca_loadings <- fviz_pca_var(imp_pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#FCA636FF", "#B12A90FF", "#0D0887FF" ),
             repel = TRUE     # Avoid text overlapping
)
pca_loadings
ggsave(here("figs", "indicators_pca_loadings.png"), width=8, height=6, units="in")

# plot pca with autoplot() function from ggfortify, colored by soil series
# note - I had originally tried to do this as a biplot but I think that adding the loadings on top of the points really makes things too confusing
pca_soil <- autoplot(imp_pca2, data=soc_horizon_filt, colour="soil") +
  scale_color_viridis(discrete=TRUE, name="Soil Series") +
  xlim(c(-.25, .15)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_katy()
pca_soil
ggsave(here("figs", "indicators_pca_soil.png"), width=8, height=6, units="in")

# plot pca colored by management
pca_mgmt <- autoplot(imp_pca2, data=soc_horizon_filt, colour="label") +
  scale_color_viridis(discrete=TRUE, name="Management") +
  xlim(c(-.25, .15)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_katy()
pca_mgmt
ggsave(here("figs", "indicators_pca_mgmt.png"), width=8, height=6, units="in")

# Put two together in a grid
plot_grid(pca_soil, pca_mgmt)
ggsave(here("figs", "indicator_pca_grid.png"), width=16, height=6, units="in")