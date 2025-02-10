# Archive - initial examination of SOC data for all pedons

# 0 - Import data ####
soc_pedon <- read.csv(here("data_processed", "04_soc_stock_pedon.csv")) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")))
soc_horizon <- read.csv(here("data_processed", "04_soc_stock_horizon.csv")) %>%
  mutate(label = factor(label, levels=c("BAU", "SHM", "Ref")))

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
  group_by(project, soil, label) %>%
  summarize(across(soc_stock_30cm:soc_stock_100cm, mean_sd), n=n())
flextable(soc_summary)
# The summary table shows that some some projects didn't have a consistent soil type across the different treatments, or just didn't have sufficient SOC data to calculate 30 cm or 100 cm depths. Data from these soils should be excluded from further analysis. 
# These soils are: Canton (only in SHM), Hidalgo (has 30 cm stocks, not 100 cm), Pullman (only has SOC %), Kenyon (only in Ref), Marquis (only in SHM), Readlyn (only in BAU), and Woodbridge