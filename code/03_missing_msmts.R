# 03 - Identifying missing lab and field measurements
# Katy Dynarski, March 2024

# This script identifies horizons with missing SOC, bulk density, or coarse fragments data (which are all necessary to calculate SOC stocks) and fills that data in to minimize data loss (i.e. having to throw out entire soil profiles because one horizon is missing a measurement). Missing values for a horizon are typically filled with the average value for a generalized horizon within that particular combination of treatment and project

# Working off of a dataframe Ekundayo put together that fills in missing cooperator SOC values with KSSL SOC - the combined values are in a column called "soc_pct_mod"

# 0 - Load libraries and import data ----
library(aqp)
library(here)
library(zoo)
library(soilDB)
library(sf)
library(tidyverse)

coop_data <- read.csv(here("data_processed","02_coop_data_horizons_valid.csv"))
coop_spc <- coop_data

# 1 - SOC% - Identify missing data ----
no_soc <- coop_data %>%
  filter(is.na(soc_pct_mod))
# 27 samples missing SOC% (most of them are from the Texas A&M project where SOC was only measured at KSSL)

# Promote to Soil Profile Collection to take a look at the profiles
depths(coop_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(coop_spc) <- 'hzdesg'

# Subset for profiles with missing SOC measurements and plot to see where the SOC measurements are missing
no_soc_spc <- subset(coop_spc, dsp_pedon_id %in% no_soc$dsp_pedon_id)
plotSPC(no_soc_spc, color="soc_pct_mod")

# 2 - SOC% - Fill missing values based on average SOC for horizon/treatment ####

# Calculate grouped means for generalized horizons grouped by treatment
# Use na.aggregate() function from zoo pkg, which will fill NA values with grouped means
# Try this for UnivOfMinnesota subset
minn_soc <- subset(coop_spc, project=="UnivOfMinnesota")

# First, test method in a separate df
soc_fill_minn <- horizons(minn_soc) %>%
  group_by(label, genhz) %>% # group by generalized horizons within treatment to allow us to get means
  mutate(soc_fill = na.aggregate(soc_pct_mod)) # fills any NAs with group mean
is.na(soc_fill_minn$soc_fill) # check for any remaining NAs

# Double check fill values are correct
check <- soc_fill_minn %>%
  group_by(label, genhz) %>%
  summarize(soc_pct_mean = mean(soc_pct_mod, na.rm=TRUE)) %>%
  unite('label_genhz', label:genhz, remove=FALSE)

no_soc_label <- soc_fill_minn %>%
  filter(is.na(soc_pct_mod)) %>%
  select(label, genhz, soc_fill) %>%
  unite('label_genhz', c('label', 'genhz'), remove=FALSE)

check2 <- check %>%
  filter(label_genhz %in% no_soc_label$label_genhz) %>%
  left_join(no_soc_label, by=c('label', 'genhz', 'label_genhz'))
check2 # want values in 'soc_fill' and 'soc_pct_mod' to match for each row
# yes they are the same

# Repeat in SPC subset to check profile plot
horizons(minn_soc) <- horizons(minn_soc) %>%
  group_by(label, genhz) %>% # group by generalized horizons 
  mutate(soc_fill = na.aggregate(soc_pct_mod)) # fills any NAs with group mean
plotSPC(minn_soc, color="soc_fill") # everything is filled

# This works - fill SOC values in actual dataframe
# Will need to add 'project' as a grouping variable for na.aggregate() so that each filled value represents the average value for that horizon within a treatment for a particular project
coop_data_soc_filled <- coop_data %>%
  group_by(project, label, genhz) %>%
  mutate(soc_fill = ifelse(is.na(genhz), soc_pct_mod, na.aggregate(soc_pct_mod))) # ifelse statement fills in original values if there is no genhz (and will leave NAs if there are no SOC measurements and no genhz labels), and fills NA data with mean SOC% according to treatment and genhz if there is a genhz label

# Check that that worked
tam_soc <- coop_data_soc_filled %>%
  filter(project=="TexasA&MPt-2") %>%
  select(dsp_pedon_id, project, label, genhz, hrzdep_t, hrzdep_b, soc_pct, soc_pct_mod, soc_fill)

# Check for remaining NAs
coop_data_soc_filled %>%
  filter(is.na(soc_fill))
# just Texas A&M Pt 2 project

# 3 - Bulk Density - Identify missing data ####
no_bd <- coop_data %>%
  filter(is.na(bulk_density))
# Many samples are missing bulk density measurements - 221!

# In some projects, it looks like only the deepest soils in a pedon are missing BD measurements. In other cases, it looks like whole pedons did not have BD measurements collected.

# Subset by project to understand where and why BD measurements are missing - this will influence method necessary to fix
projects_no_bd <- distinct(no_bd, project)
projects_no_bd$project

# Fill BD values separately for each project and then join together, then join with filled SOC data to make completed df

# 4 - Bulk Density - Fill in missing data for each project ####
# 4.1 - KansasState BD Fill ####
ks_bd <- subset(coop_spc, project=="KansasState")
plotSPC(ks_bd, color='bulk_density')
# Horizons with missing bulk density data are a bit random - KeC1-3 is missing three bottom horizons, KeCF1-3 is missing a measurement in the Bt1 horizon, and KeCF3-3 is missing a measurement in the very bottom of the pedon
# Fill these back in based on ghl and treatment

horizons(ks_bd) <- horizons(ks_bd) %>%
  group_by(label, genhz) %>% # group by generalized horizons, don't need project as a grouping variable because this subset is all one project
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(ks_bd, color="bd_fill") # values all look reasonable

# Make dataframe
ks_bd_filled <- horizons(ks_bd)

# 4.2 - NCState BD Fill ####
ncs_bd <- subset(coop_spc, project=="NCState")
plotSPC(ncs_bd, color='bulk_density')
# Missing BD values are few and random - can use grouped NA aggregation method here

horizons(ncs_bd) <- horizons(ncs_bd) %>%
  group_by(label, genhz) %>% # group by generalized horizons
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(ncs_bd, color="bd_fill")

# Make dataframe
ncs_bd_filled <- horizons(ncs_bd)

# 4.3 - WashingtonState BD Fill ####
wash_bd <- subset(coop_spc, project=="WashingtonState")
plotSPC(wash_bd, color='bulk_density')
# Lots of missing measurements, but they are all pretty random - try grouped NA aggregation

horizons(wash_bd) <- horizons(wash_bd) %>%
  group_by(label, genhz) %>% # group by generalized horizons
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(wash_bd, color="bd_fill") # all looks reasonable

# Make dataframe
wash_bd_filled <- horizons(wash_bd)

# 4.4 - UnivofMinnesota BD Fill ####
minn_bd <- subset(coop_spc, project =="UnivOfMinnesota")
plotSPC(minn_bd, color='bulk_density')
# missing one value, can easily fix this

horizons(minn_bd) <- horizons(minn_bd) %>%
  group_by(label, genhz) %>% # group by generalized horizons
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(minn_bd, color="bd_fill")

# Make dataframe
minn_bd_filled <- horizons(minn_bd)

# 4.5 - Illinois BD Fill ####
ill_bd <- subset(coop_spc, project=="Illinois")
plotSPC(ill_bd, color="bulk_density")
# none of the samples that had full profiles characterized had bulk density collected, plus a whole site appears to not have BD

# Extra twist is that most pedons don't have horizons designated - will need to calculate by depth class?

# try to plot by treatment to get a better sense of what we do and don't have data for
# Make a duplicate for playing with
ill2 <- ill_bd
site(ill2) <- ~ label #promote 'label' (where BAU vs Ref vs SHM lives) to site-level
site(ill2) <- ~ site  #promote site to site-level as well

# Examine by treatment
groupedProfilePlot(ill2, groups='label', color='bulk_density')
# Bulk density measurements only exist for 100-cm pedons in BAU and Ref (not the deep pedons, and no BD for SHM pedons)
groupedProfilePlot(ill2, groups='label', color='soc_pct_mod')
# All pedons have SOC data, though

# Examine by site
groupedProfilePlot(ill2, groups='site', color='bulk_density')
# Each site (except for Arboretum) has a full pit (deeper than 100cm), the full pit does not have bulk density measurements but the 100 cm pits do. At the Rothermel site there are no BD measurements. I won't be able to fill in bulk density for those deep horizons, but should be able to fill everything else...

# Step 1 - Add a depth class identifier column
ill_sub <- horizons(ill_bd) %>%
  mutate(hrzdep_t_chr = as.character(hrzdep_t),
         hrzdep_b_chr = as.character(hrzdep_b)) %>%
  unite(depth_class, hrzdep_t_chr:hrzdep_b_chr, sep='-')

# Step 2 - not all the depth classes are the same - use AQP dice() to slice by depth, calculate average BD for each depth
# make SPC
ill_sub_spc <- ill_sub
depths(ill_sub_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill_sub_spc) <- 'hzdesg'

# dice
ill_bd_dice <- aqp::dice(ill_sub_spc, fm=0:max(ill_sub_spc) ~ bulk_density + dsp_pedon_id + depth_class)
plotSPC(ill_bd_dice, color='bulk_density') # this looks horrible but it will work

# Step 3 - fill in the missing data with average bd for that 1 cm slice and calculate mean BD for original depth class
horizons(ill_bd_dice) <- horizons(ill_bd_dice) %>%
  group_by(hrzdep_b) %>% # group by horizon depth
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(ill_bd_dice, color='bd_fill')

# Step 4 - calculate mean BD for each depth class and add BD data back in
ill_bd_to_join <- horizons(ill_bd_dice) %>%
  group_by(dsp_pedon_id, depth_class) %>%
  summarize(bd_class_avg = mean(bd_fill))

# join calculated BD in to dataframe
ill_sub_bd_filled <- ill_sub %>%
  left_join(ill_bd_to_join, by=c("dsp_pedon_id", "depth_class")) %>%
  rename(bd_fill = bd_class_avg)

# make a copy to promote to SPC and plot profiles to check that everything looks right
ill_sub_filled_spc <- ill_sub_bd_filled 
depths(ill_sub_filled_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill_sub_filled_spc) <- 'hzdesg'
plotSPC(ill_sub_filled_spc, color="bd_fill") # looks good

# 4.6 - UTRGV BD Fill ####
utrgv_data <- coop_data %>%
  filter(project=="UTRGV")

utrgv_bd <- subset(coop_spc, project=="UTRGV")
plotSPC(utrgv_bd, color='bulk_density')
# Bulk density was only collected to 30 cm

utrgv2 <- utrgv_bd
site(utrgv2) <- ~ site
site(utrgv2) <- ~ label

groupedProfilePlot(utrgv2, groups='label', color='bulk_density')
groupedProfilePlot(utrgv2, groups='label', color='soc_pct')
# Missing BD below 30cm but do have SOC

groupedProfilePlot(utrgv2, groups='site', color='bulk_density')

# Are there other projects on the same series?
hidalgo <- coop_data %>%
  filter(soil=="Hidalgo") %>%
  distinct(project) # No, just UTRGV

# Fill this data with representative values for Hidalgo soil series from SSURGO
# get pedon locations
hidalgo_locations <- coop_data %>%
  filter(project=="UTRGV") %>%
  select(project, label, trt, soil, pedon_x, pedon_y) %>%
  group_by(project, label, trt) %>%
  distinct() %>%
  slice(1)

# convert to sf object
hidalgo_pts <- st_as_sf(hidalgo_locations, coords = c("pedon_x", "pedon_y"), crs = 4326)

# SDA query
# First get map unit keys
hidalgo_mu <- SDA_spatialQuery(hidalgo_pts, what = 'mukey', geomIntersection = TRUE, byFeature=TRUE)
# bind to pedons
hidalgo_pts_mu = cbind(hidalgo_pts, hidalgo_mu)

# Get component data
hidalgo_comp_list <- sprintf("mukey IN %s", format_SQL_in_statement(hidalgo_mu$mukey))
# query to return mukey
hidalgo_comp_data <- fetchSDA(WHERE = hidalgo_comp_list, duplicates = TRUE, childs = FALSE, 
                              nullFragsAreZero = TRUE, stringsAsFactors = FALSE, rmHzErrors = FALSE)
plot(hidalgo_comp_data, color = "dbthirdbar_r", label = 'compname', cex.names = 0.75)

hidalgo_comp_db <- as(hidalgo_comp_data, "data.frame") %>%
  filter(compname %in% hidalgo_locations$soil) %>%
  select(compname, comppct_r, hzname, hzdept_r, hzdepb_r, dbthirdbar_r) %>%
  distinct()

# use AQP dice() to slice by depth, calculate average BD for each depth
# make SPC
hidalgo_spc <- hidalgo_comp_db
depths(hidalgo_spc) <- compname ~ hzdept_r + hzdepb_r
hzdesgnname(hidalgo_spc) <- 'hzname'

# dice
hidalgo_dice <- aqp::dice(hidalgo_spc, fm=0:max(hidalgo_spc) ~ dbthirdbar_r + hzname)
plotSPC(hidalgo_dice, color='dbthirdbar_r')

# make diced profile into dataframe, add in depth class of original sampling, calculate mean db for each depth class
hidalgo_mean_bd <- horizons(hidalgo_dice) %>%
  filter(hzdepb_r <= 100) %>% # remove data below 100 cm depth
  mutate(sampled_depth_bottom = case_when(hzdepb_r <=5 ~ 5,
                                          hzdepb_r <=10 ~ 10,
                                          hzdepb_r <=30 ~ 30,
                                          hzdepb_r <=50 ~ 50,
                                          hzdepb_r <=80 ~ 80,
                                          hzdepb_r <=100 ~ 100)) %>%
  group_by(sampled_depth_bottom) %>%
  summarize(ssurgo_mean_bd = mean(dbthirdbar_r))

# fill in the missing data 
utrgv_bd_filled <- utrgv_data %>%
  left_join(hidalgo_mean_bd, by=c("hrzdep_b" = "sampled_depth_bottom")) %>%
  group_by(hrzdep_b) %>% # group by horizon depth
  mutate(bd_fill = case_when(!is.na(bulk_density) ~ bulk_density,
                             is.na(bulk_density) ~ ssurgo_mean_bd)) %>%
  select(-ssurgo_mean_bd)

# 4.7 - UConn BD Fill ####
uconn_bd <- subset(coop_spc, project=="UConn")
plotSPC(uconn_bd, color='bulk_density')
# Missing measurements are few and random - can use grouped NA aggregation

horizons(uconn_bd) <- horizons(uconn_bd) %>%
  group_by(label, genhz) %>% # group by generalized horizons
  mutate(bd_fill = na.aggregate(bulk_density)) # fills any NAs with group mean
plotSPC(uconn_bd, color="bd_fill")

# Make dataframe
uconn_bd_filled <- horizons(uconn_bd)

# 4 - Bulk Density - Join project data back together ####
# need to also add back in the data where no calculations were performed - Texas A&M pt1and2, and OSU
coop_data_bd_fine <- coop_data %>%
  filter(project=="TexasA&MPt-1" | project=="TexasA&MPt-2" | project=="OregonState") %>%
  mutate(bd_fill = bulk_density)

coop_data_bd_filled <- bind_rows(ks_bd_filled, ncs_bd_filled, wash_bd_filled, minn_bd_filled, 
                                 ill_sub_bd_filled, uconn_bd_filled, utrgv_bd_filled) %>%
  select(!c('hzID', 'depth_class')) %>% # all data where BD were missing and needed to be calculated
  bind_rows(coop_data_bd_fine)

# 5 - Coarse Fragments - Identify and fill missing data ####
# Coarse fragments are a little different than the other two variables we have been looking at - they are needed to calculate SOC stock, but values are often NA in dataframe not because they were overlooked, but because they weren't observed

# It might be easier to start by looking at where there ARE coarse fragments
coarse <- coop_data %>%
  filter(coarse_frag_volume > 0)
# Coarse fragments are only identified in three projects - OSU, Illinois and UConn
# In Illinois, they are only found in deepest horizons (well below 100cm) for the full characterization pedons - can disregard these
# Let's look at UConn data
plotSPC(uconn_bd, color='coarse_frag_volume')

# For all other projects, can just assign 0 to NA values so that calculations are possible
coop_data_coarse_filled <- coop_data %>%
  mutate(coarse_frag_fill = ifelse(is.na(coarse_frag_volume), 0, coarse_frag_volume))

# 6 - Join everything together and write CSV ####
# We have built three different dataframes containing filled values for SOC, bulk density, and coarse fragments that need to be joined back together
# Do a left join based on bulk density frame (which has the fewest observations)
# Shrink SOC and CF dataframes so there is less to join in

soc_fill <- coop_data_soc_filled %>%
  ungroup() %>%
  select(dsp_sample_id, soc_fill)

cf_fill <- coop_data_coarse_filled %>%
  ungroup() %>%
  select(dsp_sample_id, coarse_frag_fill)

coop_data_filled <- coop_data_bd_filled %>%
  left_join(soc_fill, by="dsp_sample_id") %>%
  left_join(cf_fill, by="dsp_sample_id")

# check for NA values
na_check <- coop_data_filled %>%
  filter(is.na(soc_fill) | is.na(bd_fill) | is.na(coarse_frag_fill)) # only NAs left are som SOC values in Texas A&M pt 2 project, bulk density in UTRGV pits below 30 cm, and bulk density in the deep Illinois pits - can't do anything about that now :)

write_csv(coop_data_filled, here("data_processed", "03_coop_data_filled.csv"))
