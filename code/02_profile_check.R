# 02 - Construct soil profiles for each project and check data quality
# Katy Dynarski, October 2023

# This script constructs a SoilProfileCollection object for each project in the DSP4SH dataset, checks horizon logic, corrects any issues in the data, builds a corrected dataframe for each project, and re-joins to make a final, corrected dataframe 

# Load data and identify project list ####
coop_data <- read.csv(here("data_processed", "coop_data_full.csv"))
projects <- coop_data %>%
  distinct(project)
projects

# 1 - KansasState ####
ks <- coop_data %>%
  filter(project == "KansasState")

# Promote dataframe to SPC object
depths(ks) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ks) <- 'hzdesg'

# Check that horizon depths make sense
checkHzDepthLogic(ks) # All are valid except for KeC2 - 3
kec2_3 <- horizons(ks) %>%
  filter(str_detect(dsp_pedon_id, "KeC2 - 3"))
# problem is an incorrectly labelled pedon ID - there is a horizon that actually is in pedon KeC2-2 (based on its dsp_sample_id and dsp_pedon values), but the dsp_pedon_id incorrectly places it in pedon KeC2-3

# this is an easy fix ! check in copy dataframe first
kec2_3_fix <- kec2_3 %>%
  mutate(dsp_pedon_id = case_when(dsp_sample_id == "KeC2 - 2-8" ~"KeC2 - 2",
                                  TRUE ~ dsp_pedon_id)) # when dsp_sample_id is the mislabeled one, changes dsp_pedon_id, leaves dsp_pedon_id for all correct data

# fix in SPC object
horizons(ks) <- horizons(ks) %>%
  mutate(dsp_pedon_id = case_when(dsp_sample_id == "KeC2 - 2-8" ~"KeC2 - 2",
                                  TRUE ~ dsp_pedon_id))

# Check horizon logic again
checkHzDepthLogic(ks) # all good now!

# Check visually by plotting
plotSPC(ks, color="soc_pct") # I don't see any gaps in the profiles and all horizon labels look reasonable

# 2 - NCState ####
ncs <- coop_data %>%
  filter(project == "NCState")

depths(ncs) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ncs) <- 'hzdesg'

checkHzDepthLogic(ncs)
plotSPC(ncs, color="soc_pct")

# 3 - Texas A&M pt 1 ####
tam1 <- coop_data %>%
  filter(project=="TexasA&MPt-1")

depths(tam1) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(tam1) <- 'hzdesg'

checkHzDepthLogic(tam1)
plotSPC(tam1, color="soc_pct")

# 4 - Texas A&M pt 2 ####
tam2 <- coop_data %>%
  filter(project=="TexasA&MPt-2")

depths(tam2) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(tam2) <- 'hzdesg'

checkHzDepthLogic(tam2)
plotSPC(tam2)
# Plot appears blank - take a closer look
horizons_tam2 <- horizons(tam2)
# There is no SOC data for this project, additionally no horizon designation given in lab data which is why profiles print blank. Not much sense in worrying about this project for now since there is no SOC data :)

# 5 - WashingtonState ####
wash <- coop_data %>%
  filter(project=="WashingtonState")

depths(wash) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(wash) <- 'hzdesg'

checkHzDepthLogic(wash) #All are valid except for ASP2-3, check that more closely
asp2_3 <- horizons(wash) %>%
  filter(str_detect(dsp_pedon_id, "ASP2-3")) 
# problem is an incorrectly labelled pedon ID - there is a horizon that actually is in pedon ASP2-2 (based on its dsp_sample_id and dsp_pedon value), but the incorrectly places it in pedon ASP2-3

# fix in SPC object
horizons(wash) <- horizons(wash) %>%
  mutate(dsp_pedon_id = case_when(dsp_sample_id == "ASP2-2-6" ~"ASP2-2",
                                  TRUE ~ dsp_pedon_id))

# Check horizon logic again
checkHzDepthLogic(wash) # All good now!
plotSPC(wash, color="soc_pct")

# 6 - UnivofMinnesota ####
minn <- coop_data %>%
  filter(project == "UnivOfMinnesota")

depths(minn) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(minn) <- 'hzdesg'

checkHzDepthLogic(minn)
plotSPC(minn, color="soc_pct") # Looks good!

# 7 - OregonState ####
osu <- coop_data %>%
  filter(project=="OregonState")

depths(osu) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(osu) <- 'hzdesg'

checkHzDepthLogic(osu)
plotSPC(osu, color="soc_pct") #JoV2-3 appears to only go to 20 cm, is this true? Check in original data.
# yes, that is all the data there is
# 8 - Illinois ####
ill <- coop_data %>%
  filter(project=="Illinois")

depths(ill) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill) <- 'hzdesg'

checkHzDepthLogic(ill)
plotSPC(ill, color="soc_pct") #some deep profiles!!

# 9 - UTRGV ####
utrgv <- coop_data %>%
  filter(project=="UTRGV")

depths(utrgv) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(utrgv) <- 'hzdesg'

checkHzDepthLogic(utrgv)
plotSPC(utrgv, color="soc_pct") # such uniform horizon sampling!

# UConn ####
uconn <- coop_data %>%
  filter(project=="UConn")

depths(uconn) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(uconn) <- 'hzdesg'

checkHzDepthLogic(uconn) # more profiles are not valid than are - looks like this is due to overlaps and gaps

uconn_invalid <- checkHzDepthLogic(uconn) %>%
  filter(valid=="FALSE")
uconn_inv_hz <- horizons(uconn) %>%
  filter(dsp_pedon_id %in% uconn_invalid$dsp_pedon_id) # Lots of little gaps here - generally they all happen below the 10cm mark
# This error appears to be in both the lab and the field raw data files - not sure if the data were recorded incorrectly or if there really are horizons missing. Will need to consult with Ekundayo. For now, exclude those profiles from analysis.

uconn_valid <- horizons(uconn) %>%
  filter(!dsp_pedon_id %in% uconn_invalid$dsp_pedon_id)

depths(uconn_valid) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(uconn_valid) <- 'hzdesg'

checkHzDepthLogic(uconn_valid)
plotSPC(uconn_valid, color="soc_pct") # Some of these profiles are shallow, but some are really deep. This project is kind of odd.
