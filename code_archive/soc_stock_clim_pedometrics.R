# alt version of SOC stock figure for pedometrics talk

# 1 - Calculate average MAP and MAT for each soil series in the project, plot alternative version of stock boxplot for Pedometrics talk ----
avg_clim <- soc_pedon_clim %>%
  group_by(soil) %>%
  summarize(across(map:mat, ~ mean(.x, na.rm = TRUE))) %>%
  rename(avg_map = map,
         avg_mat = mat)

soc_pedon_clim2 <- soc_pedon_clim %>%
  left_join(avg_clim, by="soil")

ggplot(soc_pedon_clim2, aes(x=fct_reorder(soil, avg_map), y=soc_stock_100cm, fill=label)) +
  geom_boxplot() +
  labs(x="Soil Series", y="SOC stock to 100 cm depth (Mg/ha)") +
  scale_fill_viridis(discrete=TRUE, name="Management") +
  theme_katy()
ggsave(here("figs", "soc_stock_boxplot_soil_mgmt_ordered.png"), width=10, height=7, units="in", dpi=400)