# 05 - SOC Stock Summary Statistics
# Katy Dynarski, November 2023

# This script...

# 0 - Import data ####
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv"))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv"))

# 1 - Visual assessment of normality for pedon SOC stocks ####

# Density functions
ggplot(soc_pedon, aes(x=soc_stock_100cm)) +
  geom_density(aes(color=soil))

ggplot(soc_pedon, aes(x=soc_stock_30cm)) +
  geom_density(aes(color=soil))
# Neither are quite normal, Texas A&M SOC stocks are much lower

# QQ plots
ggplot(soc_pedon, aes(sample=soc_stock_100cm)) +
  geom_qq() +
  geom_qq_line()
ggplot(soc_pedon, aes(sample=soc_stock_30cm)) +
  geom_qq() +
  geom_qq_line()
# Both have a bit of a tail at the low and high end but are overall not terrible

# 2 - Summary statistics table with average pedon SOC stocks for different treatments and projects ####

# make list with mean and sd for summarizing
mean_sd <- list(
  mean = ~round(mean(.x, na.rm = TRUE), 2), 
  sd = ~round(sd(.x, na.rm = TRUE), 2)
)

# make table with summary statistics grouped by project and treatment
soc_summary <- soc_pedon %>%
  group_by(soil, label) %>%
  summarize(across(soc_stock_30cm:soc_stock_100cm, mean_sd), n=n())
flextable(soc_summary)

# 3 - Summary boxplots - effect of treatment and project on SOC stocks ####

# Boxplot comparing total SOC stocks (to 100cm) between treatments - I don't think this is helpful, shouldn't look at treatments without soil types since SOC stocks are so strongly driven by soil type
ggplot(soc_pedon, aes(x=label, y=soc_stock_100cm)) +
  geom_boxplot() +
  labs(x="Treatment", y="SOC stock to 100 cm depth (Mg/ha)") +
  theme_katy()
# Ref sites appear to have slightly higher SOC, SHM and BAU do not look that different - if anything SHM seems to have slightly lower SOC stocks

# Boxplot comparing total SOC stocks (100 cm) between treatments within soil types
# Need to exclude soils that don't have 100-cm stocks OR only have one soil type per treatment
soils_exclude <- c("Canton", "Hidalgo", "Pullman", "Kenyon", "Marquis", "Readlyn")
ggplot(filter(soc_pedon, !soil %in% soils_exclude), aes(x=soil, y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  labs(x="Soil", y="SOC stock to 100 cm depth (Mg/ha)") +
  theme_katy()
ggsave(here("figs", "soc_stock_boxplot_soil_mgmt.png"), width=10, height=7, units="in", dpi=400)
# can see a clear effect of soil type - soil profiles in each project tend to cluster around a similar SOC stock. It's a little easier to see effects of treatment in this plot - reference sites have wildly higher SOC stocks for some projects. Interestingly, we never really see higher SOC stocks under SHM treatment. 

# how well are 30 cm stocks correlated with 100 cm stocks?
ggplot(soc_pedon, aes(x=soc_stock_30cm, y=soc_stock_100cm)) +
  geom_point(aes(color=soil)) +
  geom_smooth(method="lm")

# 4 - Plot SOC stocks with depth ####
# promote horizon data to SPC
soc_spc <- soc_horizon
depths(soc_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(soc_spc) <- 'hzdesg'
# promote project and label to site-level so they can be used as grouping variables
site(soc_spc) <- ~ project + label + soil

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

# Plot with all mgmt together - filter out soils that don't have stocks calculated to 100cm 
ggplot(filter(slab_mgmt, !soil %in% soils_exclude), aes(x=top, y=p.q50, color=label)) +
  geom_line(linewidth=1.2) +
  geom_ribbon(aes(ymin=p.q25, ymax=p.q75, x=top, fill=label), alpha=0.2) +
  xlim(c(100,0)) +
  coord_flip() +
  labs(title="Median SOC Stocks by Depth", x="Depth", y="SOC (Mg/ha)") +
  facet_wrap(~ soil)

# Plot reference conditions only
ggplot(filter(slab_mgmt, !soil %in% soils_exclude & label=="Ref"), aes(x=top, y=p.q50)) +
  geom_line(linewidth=1.2) +
  geom_ribbon(aes(ymin=p.q25, ymax=p.q75, x=top), alpha=0.2) +
  xlim(c(100,0)) +
  coord_flip() +
  labs(title="Median SOC Stocks by Depth - Reference Condition", x="Depth (cm)", y="SOC (Mg/ha)") +
  facet_wrap(~ soil)
ggsave(here('figs', "soc_stocks_by_depth_ref.png"))

# Plot BAU vs SHM
ggplot(filter(slab_mgmt, !soil %in% soils_exclude & !label=="Ref"), aes(x=top, y=p.q50)) +
  geom_line(linewidth=1.2, aes(color=label)) +
  geom_ribbon(aes(ymin=p.q25, ymax=p.q75, x=top, fill=label), alpha=0.2) +
  xlim(c(100,0)) +
  coord_flip() +
  labs(title="Median SOC Stocks by Depth - Alternative Management Conditions", x="Depth (cm)", y="SOC (Mg/ha)") +
  facet_wrap(~ soil)
ggsave(here('figs', "soc_stocks_by_depth_mgmt.png"))

# 5 - SOC depth plot accounting for treatment ####

# Slab into broader depth categories and plot SOC stocks in different 
slab_label <- aqp::slab(soc_spc,
                        fm = label ~ soc_stock_hrz,
                        slab.structure = c(0,5,10,30,100)) %>%
  mutate(depth_increment = factor(ifelse(bottom=="5", "0-5 cm",
                                         ifelse(bottom=="10", "5-10 cm",
                                                ifelse(bottom=="30", "10-30 cm",
                                                       "30-100 cm"))),
                                  levels=c("0-5 cm", "5-10 cm", "10-30 cm", "30-100 cm")))

# Plot
ggplot(slab_label, aes(x=label, y=p.q50, fill=fct_rev(depth_increment))) +
  geom_bar(stat="identity") +
  labs(x="Treatment", y="Median SOC Stock (Mg/ha)") +
  scale_fill_discrete(name="Depth Increment") +
  guides(fill = guide_legend(reverse=TRUE)) +
  ylim(c(105,0)) +
  theme_katy()
# Figure shows that most of variation between treatments in SOC stocks is due to differences in 0-10cm depth (though interestingly, SHM seems to have lower SOC stocks on both 10-30 and 30-100 depth then BAU and Ref treatments)

# 6 - ANOVA for total SOC stocks ####
anova(aov(soc_stock_100cm ~ (project+label)^2, data=soc_pedon))
# significant effect of soil series, label (treatment), and their interaction - not surprising!!

anova(aov(soc_stock_100cm ~ (soil+label)^2, data=soc_pedon))
# if we group by soil series - only significant effect of soil series

# how is that possible
soc_pedon %>%
  group_by(project, soil) %>%
  count()
# University of Minnesota has three different soil series(!)
ggplot(filter(soc_pedon, project=="UnivOfMinnesota"), aes(x=soil, y=soc_stock_100cm, fill=label)) +
  geom_boxplot()
# each treatment is on a different soil series, makes it hard to disentangle soil series from treatment

ggplot(soc_pedon, aes(x=label, y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  facet_wrap(~ soil, scales="free_y")

# it will be really interesting to tease out - are there drivers of SOC stocks across dataset? and do those explain more variability than treatment?

# 7 - Plot SOC stocks in reference only and BAU vs SHM conditions ####

# Reference only
soc_pedon_ref <- soc_pedon %>%
  filter(label=="Ref")

ggplot(soc_pedon_ref, aes(x=soil, y=soc_stock_100cm)) +
  geom_boxplot() +
  theme_katy()

# Management comparison
soc_pedon_mgmt <- soc_pedon %>%
  filter(label=="BAU" | label=="SHM")

ggplot(soc_pedon_mgmt, aes(x=soil, y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  theme_katy()
anova(aov(soc_stock_100cm ~ (soil+label)^2, data=soc_pedon_mgmt))
# may need to remove the Minnesota project for comparison - think it throws things off to have one project where each label is on a different soil series

# 8 - Plot other indicators vs SOC concentrations/stocks ####
surf <- soc_horizon %>%
  filter(hrzdep_b == "5")

ggplot(surf, aes(x=label, y=soil_respiration)) +
  geom_boxplot() +
  facet_wrap(~soil, scales="free_y")

ggplot(surf, aes(x=soil_respiration, y=soc_stock_hrz, color=label)) +
  geom_point() +
  facet_wrap(~soil, scales="free")

# 9 - Boxplots of all indicators across soils and managements ####
# can we automate with - list of indicators, plot boxplots across soil series x label, ANOVA and generate letters?
# what are the indicators? soc_pct, soc_stock_hrz, yoder_agg_stab_mwd, p_h:ace

surf_indicator <- surf %>%
  filter(!project=="UnivOfMinnesota") %>%
  select(dsp_pedon_id, soil, label, soc_stock_hrz, soc_pct, bulk_density, yoder_agg_stab_mwd, p_h:ace) %>%
  pivot_longer(soc_stock_hrz:ace, names_to="indicator", values_to="value")

# Make list of indicators for purrr
indicators <- surf_indicator %>%
  distinct(indicator) %>%
  pull(indicator) %>%
  as.character

# test function
.x <- indicators[1]

surf_indicator %>% 
  filter(indicator == .x) %>%
  ggplot(aes(x = soil, y = value, fill=label)) +
  geom_boxplot() +
  labs(x="Soil", 
       y=glue::glue({unique(filter(surf_indicator, indicator==.y)$indicator)}), 
       title=glue::glue("Indicator Values Across Soil and Management - ", 
                       {unique(filter(surf_indicator, indicator==.y)$indicator)})) +
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
                           title=glue::glue("Indicator Values Across Soil and Management - ", 
                                            {unique(filter(surf_indicator, indicator==.x)$indicator)})) +
                      theme_katy()
                    
                    ggsave(here("figs/indicator_boxplots", glue::glue("indicator_", .x, ".png")), 
                           width=10, height=7, units="in", dpi=400)
                  })

# Nest by indicator, run ANOVA for soil+label interaction for each indicator, make dataframe of ANOVA results for each indicator
surf_nested <- surf_indicator %>%
  group_by(indicator) %>%
  nest()
sig_nested <- surf_nested %>%
  mutate(anova_obj = map(data, ~ anova(aov(value ~ (soil+label)^2, data=.x)))) %>%
  mutate(anova_tidy = map(anova_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(indicator, anova_tidy) %>%
  unnest(cols=c(anova_tidy))
ggplot(filter(sig_nested, term=="label"), aes(x=indicator, y=p.value)) +
  geom_bar(stat="identity") # look at the smallest p-values to see which indicators are most strongly influenced by mgmt
