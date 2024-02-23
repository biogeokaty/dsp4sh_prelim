# 4 - Ribbon plot of SOC stocks with depth in individual soil series ----
# promote horizon data to SPC so we can analyze using slab()

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
ggplot(filter(slab_mgmt, !soil %in% soils_exclude), aes(x=top, y=p.q50)) +
  geom_line(linewidth=1.2, aes(color=label)) +
  geom_ribbon(aes(ymin=p.q25, ymax=p.q75, x=top, fill=label), alpha=0.2) +
  xlim(c(100,0)) +
  coord_flip() +
  labs(title="Median SOC Stocks by Depth", x="Depth", y="SOC (Mg/ha)") +
  facet_wrap(~ soil)
ggsave(here('figs', "soc_stocks_by_depth_all_mgmt.png"))

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