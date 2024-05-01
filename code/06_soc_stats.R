# 05 - SOC Stock Summary Statistics and Plots
# Katy Dynarski, March 2024

# ALL THINGS SOC mostly mixed linear models...

# 0 - Import data ----
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))
soc_pedon_filt <- read.csv(here("data_processed", "05_soc_pedon_filt.csv"))
soc_horizon_filt <- read.csv(here("data_processed", "05_soc_horizon_filt.csv"))
project <- read.csv(here("data_processed", "05_project_data.csv"))

# 1 - Summary boxplots and LMER - effect of treatment and project on SOC stocks ####
# First, clean data so there are no NA values
soc_100_clean <- soc_pedon_filt %>%
  select(dsp_pedon_id, project, label, soil, soc_stock_100cm) %>%
  na.omit() 

# Boxplot comparing total SOC stocks (100 cm) between treatments within soil types
ggplot(soc_100_clean, aes(x=soil, y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  labs(x="Soil Series", y="SOC stock to 100 cm depth (Mg/ha)") +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(here("figs", "soc_stock_boxplot_soil_mgmt.png"), width=10, height=7, units="in", dpi=400)
# can see a clear effect of soil type - soil profiles in each project tend to cluster around a similar SOC stock. It's a little easier to see effects of treatment in this plot - reference sites have wildly higher SOC stocks for some projects. Interestingly, we never really see higher SOC stocks under SHM treatment. 

# 2 - Linear mixed model to look at effects of treatment on 100-cm SOC stocks ----
# Random effects: soil
# Fixed effects: label
soc_stock100_mixed <- lmer(soc_stock_100cm ~ label + (1|soil), data = soc_100_clean)
summary(soc_stock100_mixed)
# Results of mixed model: Soil series explains most of the variability in SOC stocks not explained by model. Label does not appear significant, but we can test by making a reduced model.

# Plot predictions of 100cm SOC stocks under each treatment
pred_stock100 <- ggpredict(soc_stock100_mixed, terms = c("label"))
ggplot(pred_stock100, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high)) +
  labs(x="Management", y="Predicted SOC stock to 100 cm (Mg/ha)") +
  theme_katy()
ggsave(here("figs", "predicted_soc_stock_by_treatment.png"))

# Test significance of treatment by comparing full and reduced models, use likelihood ratio test
stock100_mixed_full <- lmer(soc_stock_100cm ~ label + (1|soil), data = soc_100_clean, REML = FALSE)
stock100_mixed_reduced <- lmer(soc_stock_100cm ~ (1|soil), data = soc_100_clean, REML = FALSE)
anova(stock100_mixed_full, stock100_mixed_reduced, text="Chisq")
# Models are not significantly different - treatment is not a significant predictor of SOC stocks to 100 cm, only soil series is

# 3 - Linear mixed model to look at effects of treatment on 30-cm stocks ---- 
soc_30_clean <- soc_pedon_filt %>%
  select(dsp_pedon_id, project, label, soil, soc_stock_0_30cm) %>%
  na.omit()

ggplot(soc_30_clean, aes(x=soil, y=soc_stock_0_30cm, fill=label)) +
  geom_boxplot() +
  labs(x="Soil", y="SOC stock to 30 cm depth (Mg/ha)") +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(here("figs", "soc_30cm_stock_boxplot_soil_mgmt.png"), width=10, height=7, units="in", dpi=400)

# Also test with LME
# Random effects: soil
# Fixed effects: label
soc_stock30_mixed <- lmer(soc_stock_0_30cm ~ label + (1|soil), data = soc_30_clean)
summary(soc_stock30_mixed)
# Results of mixed model: Soil series explains most of the variability in SOC stocks. Label does not appear significant, but we can test by making a reduced model.

# Plot predictions of 30cm SOC stocks under each treatment
pred_stock30 <- ggpredict(soc_stock30_mixed, terms = c("label"))
ggplot(pred_stock30, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high)) +
  labs(x="Management", y="Predicted SOC stock to 30 cm (Mg/ha)") +
  theme_katy()
ggsave(here("figs", "predicted_soc_30cm_stock_by_treatment.png"))

# Test significance of treatment by comparing full and reduced models, use likelihood ratio test
stock30_mixed_full <- lmer(soc_stock_0_30cm ~ label + (1|soil), data = soc_30_clean, REML = FALSE)
stock30_mixed_reduced <- lmer(soc_stock_0_30cm ~ (1|soil), data = soc_30_clean, REML = FALSE)
anova(stock30_mixed_full, stock30_mixed_reduced, text="Chisq")
# Models are significantly different, indicating treatment is a a significant predictor of SOC stocks to 30cm depth! 

# Are there significantly different groups?
stock30_mixed_tukey <- glht(soc_stock30_mixed, linfct = mcp(label = 'Tukey'))
summary(stock30_mixed_tukey)
# Ref is significantly higher than SHM and BAU, SHM and BAU are not significantly different

# 4 - Depth plot of SOC stocks and concentrations by treatment ----
# Need to reshape data so that SOC stocks are long, not wide
stock_depth_df <- soc_pedon_filt %>%
  select(dsp_pedon_id, project, soil, label, soc_stock_0_30cm:soc_stock_50_100cm) %>%
  pivot_longer(soc_stock_0_30cm:soc_stock_50_100cm, names_to="depth", values_to="soc_stock") %>%
  mutate(depth = fct_recode(depth, "0-30 cm" = "soc_stock_0_30cm", "30-50 cm" = "soc_stock_30_50cm", 
                            "50-100 cm" = "soc_stock_50_100cm"))

stock_depth <- ggplot(stock_depth_df, aes(x=depth, y=soc_stock, fill=label)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  labs(title="SOC Stocks by Depth", x="Depth", y="SOC (Mg/ha)") +
  theme_katy_grid()
stock_depth
ggsave(here("figs", "soc_stock_depth.png"))

# Also try with SOC concentrations - stocks is a bit misleading because it's influenced by depth and not all depth increments are equal in size
# for this, we'll need horizon data/AQP::slab
soc_spc <- soc_horizon_filt
depths(soc_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(soc_spc) <- 'hzdesg'
# promote project and label to site-level so they can be used as grouping variables
site(soc_spc) <- ~ project + label + soil

soc_conc_slab_ref <- aqp::slab(subset(soc_spc, label=="Ref"),
                               fm = soil ~ soc_fill,
                               slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="Ref")

soc_conc_slab_shm <- aqp::slab(subset(soc_spc, label=="SHM"),
                               fm = soil ~ soc_fill,
                               slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="SHM")
soc_conc_slab_bau <- aqp::slab(subset(soc_spc, label=="BAU"),
                               fm = soil ~ soc_fill,
                               slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="BAU")

# Put management conditions together
soc_conc_slab_all <- bind_rows(soc_conc_slab_ref, soc_conc_slab_shm, soc_conc_slab_bau) %>%
  mutate(depth_increment = factor(ifelse(bottom=="5", "0-5 cm",
                                         ifelse(bottom=="10", "5-10 cm",
                                                ifelse(bottom=="30", "10-30 cm",
                                                       ifelse(bottom=="60", "30-60 cm",
                                                              "60-100 cm")))),
                                  levels=c("0-5 cm", "5-10 cm", "10-30 cm", "30-60 cm", "60-100 cm"))) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")))

conc_depth <- ggplot(soc_conc_slab_all, aes(x=depth_increment, y=p.q50, fill=label)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  labs(title="SOC Concentration by Depth", x="Depth", y="SOC Concentration (%)") +
  theme_katy_grid()
conc_depth
ggsave(here("figs", "soc_conc_depth.png"))

# Make a ribbon depth plot as well
# Calculate stocks by depth increment for each soil and management condition
slab_ref <- aqp::slab(subset(soc_spc, label=="Ref"),
                      fm = soil ~ soc_stock_hrz,
                      slab.structure = seq(0,100,by=10)) %>%
  mutate(label="Ref")

slab_shm <- aqp::slab(subset(soc_spc, label=="SHM"),
                      fm = soil ~ soc_stock_hrz,
                      slab.structure = seq(0,100,by=10)) %>%
  mutate(label="SHM")
slab_bau <- aqp::slab(subset(soc_spc, label=="BAU"),
                      fm = soil ~ soc_stock_hrz,
                      slab.structure = seq(0,100,by=10)) %>%
  mutate(label="BAU")

# Put management conditions together
slab_mgmt <- bind_rows(slab_ref, slab_shm, slab_bau)

# Plot with all mgmt together
ggplot(slab_mgmt, aes(x=top, y=p.q50)) +
  geom_line(linewidth=1.2, aes(color=label)) +
  geom_ribbon(aes(ymin=p.q25, ymax=p.q75, x=top, fill=label), alpha=0.2) +
  xlim(c(100,0)) +
  coord_flip() +
  labs(title="SOC Stocks by Depth", x="Depth (cm)", y="SOC (Mg/ha)") +
  facet_wrap(~ soil, scales = "free_x") +
  scale_fill_viridis(discrete=TRUE, guide="none") +
  scale_color_viridis(discrete=TRUE, name="Management") +
  theme_katy()
ggsave(here("figs", "soc_stock_depth_ribbon.png"))

# 5 - Linear mixed effects models for effect of treatment on SOC concentration while accounting for soil series and depth ----

# Although there is no overall difference in 100 cm SOC stocks between treatments across all experiments, there are differences in 30 cm SOC stocks, and it does look like there many be significantly higher concentration of SOC in surface horizon (0-5cm) in reference plots. Should do test for soil X depth X treatment effects - mixed linear model to allow for depth as continuous?

# Clean data for modeling - want only what's going into the model, remove any NA values
soc_conc_clean <- soc_horizon_filt %>%
  select(dsp_pedon_id, label, soil, hrzdep_t, hrzdep_b, soc_fill) %>%
  na.omit() 

# Random effects: soil
# Fixed effects: label, depth (? not sure if this should be a random effect instead of a fixed effect)

soc_conc_mixed <- lmer(soc_fill ~ label + (1|soil) + (1|hrzdep_b), data = soc_conc_clean)
summary(soc_conc_mixed)
# Results of mixed model: depth and soil both explain a chunk of variability. Soil explains much more than depth, and there is a lot of residual variability. Reference appears significant after accounting for variability of soil type and depth. 

# Plot predictions of SOC concentration for each treatment
pred_soc_conc <- ggpredict(soc_conc_mixed, terms = c("label"))
plot(pred_soc_conc)

# The question - do we get more explanatory power with label than with only random effects?
# Test by comparing full and reduced models, use likelihood ratio test
soc_conc_mixed_full <- lmer(soc_fill ~ label + (1|soil) + (1|hrzdep_b), data = soc_conc_clean, REML = FALSE)
soc_conc_mixed_reduced <- lmer(soc_fill ~ (1|soil) + (1|hrzdep_b), data = soc_conc_clean, REML = FALSE)
anova(soc_conc_mixed_full, soc_conc_mixed_reduced, text="Chisq")
# Models are significantly different

# Are there significantly different groups?
soc_conc_mixed_tukey <- glht(soc_conc_mixed, linfct = mcp(label = 'Tukey'))
summary(soc_conc_mixed_tukey)
# Ref is significantly different from both SHM and BAU. SHM and BAU are not significantly different from each other. 

# 6 - Depth plot and LMER of horizon thickness-normalized SOC stocks by treatment, soil series, and depth ----

# Calculate stocks normalized to horizon thickness so each horizon can be comparable

# Clean data for modeling - want only what's going into the model, and remove any NA values
soc_hrz_stock_clean <- soc_horizon_filt %>%
  select(dsp_pedon_id, label, soil, hrzdep_t, hrzdep_b, soc_stock_hrz, hrzdepth) %>%
  na.omit() %>%
  mutate(depth_stock = soc_stock_hrz/hrzdepth) # calculate horizon thickness-normalized SOC stocks

# Random effects: soil, pedon
# Fixed effects: label, depth
soc_stock_mixed <- lmer(depth_stock ~ label + (1|soil)  + (1|hrzdep_b), data = soc_hrz_stock_clean)
summary(soc_stock_mixed)
# Soil type and depth both appear to be significant random effects
# Label appears significant 

# Plot predictions of thickness-normalized SOC stocks with treatment
pred_soc_stock <- ggpredict(soc_stock_mixed, terms = c("label"))
plot(pred_soc_stock)

# Again, create reduced model to test if we get more explanatory power with label than with only random effects
# Test by comparing full and reduced models, use likelihood ratio test
soc_stock_mixed_full <- lmer(depth_stock ~ label + (1|soil) + (1|hrzdep_b), data = soc_hrz_stock_clean, REML = FALSE)
soc_stock_mixed_reduced <- lmer(depth_stock ~ (1|soil) + (1|hrzdep_b), data = soc_hrz_stock_clean, REML = FALSE)
anova(soc_stock_mixed_full, soc_stock_mixed_reduced, text="Chisq")
# Models are significantly different - means that label is a significant factor after accounting for depth and soil type

# Test for significantly different groups
soc_stock_mixed_tukey <- glht(soc_stock_mixed, linfct = mcp(label = 'Tukey'))
summary(soc_stock_mixed_tukey)
# Ref is significantly different from both SHM and BAU. SHM and BAU are not significantly different from each other. 

# 7 - Bulk density patterns by treatment, soil and depth ----
# Calculate depth increments for each soil series in different management conditions
bd_slab_ref <- aqp::slab(subset(soc_spc, label=="Ref"),
                         fm = soil ~ bd_fill,
                         slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="Ref")

bd_slab_shm <- aqp::slab(subset(soc_spc, label=="SHM"),
                         fm = soil ~ bd_fill,
                         slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="SHM")
bd_slab_bau <- aqp::slab(subset(soc_spc, label=="BAU"),
                         fm = soil ~ bd_fill,
                         slab.structure = c(0,5,10,30,60,100)) %>%
  mutate(label="BAU")

# Put management conditions together
bd_slab_all <- bind_rows(bd_slab_ref, bd_slab_shm, bd_slab_bau) %>%
  mutate(depth_increment = factor(ifelse(bottom=="5", "0-5 cm",
                                         ifelse(bottom=="10", "5-10 cm",
                                                ifelse(bottom=="30", "10-30 cm",
                                                       ifelse(bottom=="60", "30-60 cm",
                                                              "60-100 cm")))),
                                  levels=c("0-5 cm", "5-10 cm", "10-30 cm", "30-60 cm", "60-100 cm"))) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")))

# Plot bulk density with depth
bd_depth <- ggplot(bd_slab_all, aes(x=depth_increment, y=p.q50, fill=label)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  labs(title="Bulk Density by Depth", x="Depth", y="Bulk Density (g/cm3)") +
  theme_katy_grid()
bd_depth

# Mixed model for bulk density response to soil, depth, and management
bd_clean <- soc_horizon_filt %>%
  select(dsp_pedon_id, label, soil, hrzdep_t, hrzdep_b, bd_fill) %>%
  na.omit() 

# Random effects: soil
# Fixed effects: label, depth

bd_mixed <- lmer(bd_fill ~ label + (1|soil) + (1|hrzdep_b), data = bd_clean)
summary(bd_mixed)
# Results of mixed model: Soil explains the most of the random effects variability, Ref appears significant as a fixed effect

# Plot predictions of SOC concentration for each treatment
pred_bd <- ggpredict(bd_mixed, terms = c("label"))
plot(pred_bd)

# Test by comparing full and reduced models, use likelihood ratio test
bd_mixed_full <- lmer(bd_fill ~ label + (1|soil) + (1|hrzdep_b), data = bd_clean, REML = FALSE)
bd_mixed_reduced <- lmer(bd_fill ~ (1|soil) + (1|hrzdep_b), data = bd_clean, REML = FALSE)
anova(bd_mixed_full, bd_mixed_reduced, text="Chisq")
# Models are significantly different

# Are there significantly different groups?
bd_mixed_tukey <- glht(bd_mixed, linfct = mcp(label = 'Tukey'))
summary(bd_mixed_tukey)
# Ref is significantly different from both SHM and BAU. SHM and BAU are not significantly different from each other. 
# Lower bulk density in reference treatment - this makes sense

# 8 - Depth plot of SOC concentration, stocks, and bulk density with shared legend ---

# remove legend from figs for gridding
stock_depth2 <- stock_depth + theme(legend.position="none")
conc_depth2 <- conc_depth + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank()) #remove y axis from conc plot so axes can be shared
bd_depth2 <- bd_depth + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_blank())
depth_plot_legend <- get_legend(stock_depth)

depth_plot_grid <- plot_grid(stock_depth2, conc_depth2, bd_depth2, depth_plot_legend, labels="", 
                             rel_widths=c(1.4,1,1, .5), nrow=1)
depth_plot_grid
ggsave(here("figs", "soc_stock_conc_depth.png"), height=7, width=14, units="in")

# 9 - Mixed effects model of surface horizon SOC stocks and concentrations by soil series and treatment ----
# this is a prime candidate for automating with purrr!

# Subset just 0-5 and 5-10 cm horizons since these are consistent across all experiments
soc_stock_clean_5 <- soc_hrz_stock_clean %>%
  filter(hrzdep_b == 5)

soc_stock5_mixed <- lmer(soc_stock_hrz ~ label + (1|soil), data = soc_stock_clean_5)
summary(soc_stock5_mixed)

# Plot predictions 
pred_soc_stock5 <- ggpredict(soc_stock5_mixed, terms = c("label"))
ggplot(pred_soc_stock5, aes(x=x, y=predicted)) +
  geom_point() +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high)) +
  labs(x="Management", y="Predicted SOC stock in 0-5 cm horizon (Mg/ha)") +
  theme_katy()
ggsave(here("figs", "predicted_surface_soc_stock_by_treatment.png"))

# Test significance of reduced model
soc_stock5_mixed_full <- lmer(soc_stock_hrz ~ label + (1|soil), data = soc_stock_clean_5, REML = FALSE)
soc_stock5_mixed_reduced <- lmer(soc_stock_hrz ~ (1|soil), data = soc_stock_clean_5, REML = FALSE)
anova(soc_stock5_mixed_full, soc_stock5_mixed_reduced, text="Chisq")
# Two models are significantly different

# test for significantly different groups
soc_stock5_mixed_tukey <- glht(soc_stock5_mixed, linfct = mcp(label = 'Tukey'))
summary(soc_stock5_mixed_tukey)
# All groups are significantly different!!

# Repeat with 5-10 cm depth
soc_stock_clean_10 <- soc_hrz_stock_clean %>%
  filter(hrzdep_b == 10)

soc_stock10_mixed <- lmer(soc_stock_hrz ~ label + (1|soil), data = soc_stock_clean_10)
summary(soc_stock10_mixed)

# Plot predictions 
pred_soc_stock10 <- ggpredict(soc_stock10_mixed, terms = c("label"))
plot(pred_soc_stock10)

# Test significance of reduced model
soc_stock10_mixed_full <- lmer(soc_stock_hrz ~ label + (1|soil), data = soc_stock_clean_10, REML = FALSE)
soc_stock10_mixed_reduced <- lmer(soc_stock_hrz ~ (1|soil), data = soc_stock_clean_10, REML = FALSE)
anova(soc_stock10_mixed_full, soc_stock10_mixed_reduced, text="Chisq")
# Two models are significantly different

soc_stock10_mixed_tukey <- glht(soc_stock10_mixed, linfct = mcp(label = 'Tukey'))
summary(soc_stock10_mixed_tukey)
# Only Ref is significantly different from other treatments