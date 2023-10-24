# Misc AQP data exploration of DSP4SH data v2
# Katy Dynarski, October 2023

# Import and subset data, initialize SoilProfileCollection
coop_data <- read.csv(here("data_processed", "coop_data_full.csv"))
ks <- coop_data %>%
  filter(project == "KansasState")

# Promote dataframe to SPC object
depths(ks) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ks) <- 'hzdesg'

# Exploratory data analysis ####
# This is based on me working through the "Exploratory Data Analysis" section of Statistics for Soil Survey

# Make dataframe of horizon data (if I were doing this will full dataset, it should end up being the same as the layer_join df)
ksh <- horizons(ks)

# Check text data
sort(unique(ksh$hzdesg))
sort(table(ksh$hzdesg), decreasing=TRUE)
sort(table(ksh$hrzdep_b), decreasing=TRUE)
summary(ksh$hrzdep_b)
ggplot(ksh, aes(x=hrzdep_b)) +
  geom_density()

# Generate summary statistics for some variables of interest
ksh %>%
  select(hzdesg, hrzdep_t, hrzdep_b, soc_pct, bulk_density, pox_c, clay_pct_field, soil_respiration, yoder_agg_stab_mwd) %>%
  summary()

# Assign generalized horizon labels based on visual inspection of profiles
# Sequence: A, Bt, Btk, Bk 
n <- c('A', 'Bt', 'Btk', 'Bk') # generalized horizon label sequence
p <- c('^A',
       '^Bt$|^Bt1|^Bt2|^Bt3|^Bt4',
       '^Btk',
       '^Bk') #regex pattern of horizons to be corresponded to GHL seq

ksh$genhz <- generalize.hz(ksh$hzdesg, n, p) # generate labels

# Check labels
tab <- table(ksh$genhz, ksh$hzdesg)
addmargins(tab)

# Summarize some variables of interest grouped by GHL
ksh %>%
  group_by(genhz) %>%
  summarize(soc_avg = mean(soc_pct, na.rm=TRUE),
            soc_sd = sd(soc_pct, na.rm=TRUE),
            clay_avg = mean(clay_pct_field, na.rm=TRUE),
            clay_sd = sd(clay_pct_field, na.rm=TRUE),
            bd_avg = mean(bulk_density, na.rm=TRUE),
            bd_sd = sd(bulk_density, na.rm=TRUE)
  )

# Make correlation matrix
ksh %>%
  select(soc_pct, bulk_density, pox_c, clay_pct_field, soil_respiration, yoder_agg_stab_mwd) %>%
  cor(use="complete.obs") %>%
  round(2)

# Look at distributions of some continuous variables with density plots
ggplot(ksh, aes(x=soc_pct)) +
  geom_density()
ggplot(ksh, aes(x=pox_c)) +
  geom_density()
ggplot(ksh, aes(x=soil_respiration)) +
  geom_density()
ggplot(ksh, aes(x=clay_pct_field)) +
  geom_density()
ggplot(ksh, aes(x=bulk_density)) +
  geom_density()

# Look at distribution of variables within horizons
ggplot(ksh, aes(x=genhz, y=soc_pct)) +
  geom_boxplot() # unsurprisingly most SOC is in A horizons
ggplot(ksh, aes(x=genhz, y=pox_c)) +
  geom_boxplot() # similar pattern for POXC
ggplot(ksh, aes(x=genhz, y=bulk_density)) +
  geom_boxplot() # bulk density is pretty consistent throughout profile
ggplot(ksh, aes(x=genhz, y=clay_pct_field)) +
  geom_boxplot() # clay content highest in Bt horizon
ggplot(ksh, aes(x=genhz, y=soil_respiration)) +
  geom_boxplot() # respiration highest in A horizon, really thrown off by one high measurement

# QQ plots to check distribution
ggplot(ksh, aes(sample=soc_pct)) +
  geom_qq() +
  geom_qq_line()

ggplot(ksh, aes(sample=clay_pct_field)) +
  geom_qq() +
  geom_qq_line()

ggplot(ksh, aes(sample=bulk_density)) +
  geom_qq() +
  geom_qq_line()

# Plot SOC against POXC
ggplot(ksh, aes(x=soc_pct, y=pox_c, color=genhz)) +
  geom_point()

# Can also look at variable distributions within different horizons
ggplot(ksh, aes(x=soc_pct, color=genhz)) +
  geom_density(linewidth=2)

ggplot(ksh, aes(x=genhz, fill=field_texture)) +
  geom_bar()
