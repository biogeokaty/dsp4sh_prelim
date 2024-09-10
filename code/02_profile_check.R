# Checking profile validity for new DSP4SH data
# Katy DYnarski, March 2024

# This script constructs a SoilProfileCollection object for each project in the DSP4SH dataset, checks horizon logic, corrects any issues in the data, adds generalized horizon labels to the SPC, builds a corrected dataframe for each project, and re-joins to make a final, corrected dataframe 

# 0 - Load data and identify project list ####
coop_data <- read.csv(here("data_raw", "dsp4sh_soc_df_KAD.csv")) %>%
  clean_names()
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

# pull out the incorrect pedon to look at more closely
kec2_3 <- horizons(ks) %>%
  filter(str_detect(dsp_pedon_id, "KeC2-3"))
# problem is an incorrectly labelled pedon ID - there is a horizon that actually is in pedon KeC2-2 (based on its dsp_sample_id and dsp_pedon values), but the dsp_pedon_id incorrectly places it in pedon KeC2-3

# fix in SPC object
horizons(ks) <- horizons(ks) %>%
  mutate(dsp_pedon_id = case_when(dsp_sample_id == "KeC2-2-8" ~"KeC2-2",
                                  TRUE ~ dsp_pedon_id))

# Check horizon logic again
checkHzDepthLogic(ks) # all good now!

# Check visually by plotting
plotSPC(ks, color="soc_pct") # I don't see any gaps in the profiles and all horizon labels look reasonable

# Add generalized horizon labels
# Sequence: A, Bt, Btk, Bk 
n_ks <- c('A', 'Bt', 'Btk', 'Bk') # generalized horizon label sequence
p_ks <- c('^A',
          '^Bt$|^Bt1|^Bt2|^Bt3|^Bt4',
          '^Btk',
          '^Bk')

ks$genhz <- generalize.hz(ks$hzdesg, n_ks, p_ks) # generate labels

# Visually inspect assignment and check that everything looks right
tab_ks <- table(ks$genhz, ks$hzdesg)
addmargins(tab_ks)

# save corrected dataframe
ks_corr <- horizons(ks)

# 2 - NCState ####
ncs <- coop_data %>%
  filter(project == "NCState")

depths(ncs) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ncs) <- 'hzdesg'

checkHzDepthLogic(ncs) 
plotSPC(ncs, color="soc_pct") # all good

# Add generalized horizon labels
# Sequence: A, B, Bt, BC, C
n_ncs <- c('A', 'B', 'Bt', 'BC', 'C') # generalized horizon label sequence
p_ncs <- c('^Ap|^A$',
           '^B$|BA|B/A|A/B|E',
           '^Bt',
           '^BC',
           '^C|\\dC|3Abp|3Bbt')
ncs$genhz <- generalize.hz(ncs$hzdesg, n_ncs, p_ncs) # generate labels

# Visually inspect assignment and check that everything looks right
tab_ncs <- table(ncs$genhz, ncs$hzdesg)
addmargins(tab_ncs)

plotSPC(ncs, color='genhz')

# save corrected dataframe
ncs_corr <- horizons(ncs)

# 3 - Texas A&M pt 1 ####
tam1 <- coop_data %>%
  filter(project=="TexasA&MPt-1")

depths(tam1) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(tam1) <- 'hzdesg'

checkHzDepthLogic(tam1)
plotSPC(tam1, color="soc_pct")

# Add generalized horizon labels
# Sequence: A, Bt
n_tam1 <- c('A', 'Bt') # generalized horizon label sequence
p_tam1 <- c('^Ap|^A$',
            'Bt')
tam1$genhz <- generalize.hz(tam1$hzdesg, n_tam1, p_tam1) # generate labels

# Visually inspect assignment and check that everything looks right
tab_tam1 <- table(tam1$genhz, tam1$hzdesg)
addmargins(tab_tam1)

plotSPC(tam1, color='genhz')

tam1_corr <- horizons(tam1)

# 4 - Texas A&M pt 2 ####
tam2 <- coop_data %>%
  filter(project=="TexasA&MPt-2")

depths(tam2) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(tam2) <- 'hzdesg'

checkHzDepthLogic(tam2)
plotSPC(tam2)
# Plot appears blank - take a closer look

tam2_hz <- horizons(tam2)
# There is no SOC data for this project, additionally no horizon designation given in lab data which is why profiles print blank. 
# Can use KSSL to fill in missing data

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

# Add generalized horizon labels
# Sequence: A, AB, Bw
n_wash <- c('A', 'AB', 'Bw') # generalized horizon label sequence
p_wash <- c('^Ap|^A$',
            'AB',
            '^Bw')
wash$genhz <- generalize.hz(wash$hzdesg, n_wash, p_wash) # generate labels

# Visually inspect assignment and check that everything looks right
tab_wash <- table(wash$genhz, wash$hzdesg)
addmargins(tab_wash)

plotSPC(wash, color='genhz')

# save corrected dataframe
wash_corr <- horizons(wash)

# 6 - UnivofMinnesota ####
minn <- coop_data %>%
  filter(project == "UnivOfMinnesota")

depths(minn) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(minn) <- 'hzdesg'

checkHzDepthLogic(minn)
plotSPC(minn, color="soc_pct") # Looks good!

# Add generalized horizon labels
# Sequence: Ap1, Ap2, A, AB, BA, Bw, 2Bw, Bt
n1 <- c('Ap1', 'Ap2', 'A', 'AB', 'BA', 'Bw', '2Bw', 'Bt') # generalized horizon label sequence
p1 <- c('Ap1',
        'Ap2',
        '^A$|Ap3',
        'AB$',
        'BA$',
        '^Bw|2C1|2BC',
        '^2Bw|2BG',
        'Bt')

minn$genhz <- generalize.hz(minn$hzdesg, n1, p1) # generate labels

# Visually inspect assignment and check that everything looks right
tab_minn <- table(minn$genhz, minn$hzdesg)
addmargins(tab_minn)

# Extract data from SPC object to dataframe for later join
minn_corr <- horizons(minn)

# 7 - OregonState ####
osu <- coop_data %>%
  filter(project=="OregonState")

depths(osu) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(osu) <- 'hzdesg'

checkHzDepthLogic(osu)
plotSPC(osu, color="soc_pct") #JoV2-3 appears to only go to 20 cm, is this true? Check in original data.
# yes, that is all the data there is :()

# Add generalized horizon labels
# Sequence: A, AB, Bt
n_osu <- c('A', 'AB', 'Bt') # generalized horizon label sequence
p_osu <- c('^Ap|^A$',
           'AB|BA',
           '^Bt|\\dBt|2BC|BCt')
osu$genhz <- generalize.hz(osu$hzdesg, n_osu, p_osu) # generate labels

# Visually inspect assignment and check that everything looks right
tab_osu <- table(osu$genhz, osu$hzdesg)
addmargins(tab_osu)

plotSPC(osu, color='genhz')

osu_corr <- horizons(osu)

# 8 - Illinois ####
ill <- coop_data %>%
  filter(project=="Illinois")

depths(ill) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(ill) <- 'hzdesg'

checkHzDepthLogic(ill)
plotSPC(ill, color="soc_pct") #some deep profiles!!

# Add generalized horizon labels
# Sequence: A, Bt, C
n_ill <- c('A', 'Bt', 'C') # generalized horizon label sequence
p_ill <- c('^Ap|^A$|AB',
           'Bt|BAt',
           'C')
ill$genhz <- generalize.hz(ill$hzdesg, n_ill, p_ill) # generate labels

# Visually inspect assignment and check that everything looks right
tab_ill <- table(ill$genhz, ill$hzdesg)
addmargins(tab_ill)

plotSPC(ill, color='genhz') # only deep some pedons have horizons assigned at all - can I assign based on depth?
# bookmark in this

# For now, leave the missing generalized horizon labels

ill_corr <- horizons(ill)

# 9 - UTRGV ####
utrgv <- coop_data %>%
  filter(project=="UTRGV")

depths(utrgv) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(utrgv) <- 'hzdesg'

checkHzDepthLogic(utrgv)
plotSPC(utrgv, color="soc_pct") # such uniform horizon sampling!

# Add generalized horizon labels
# Sequence: A, B, B2ca
n_utrgv <- c('A', 'B', 'BC') # generalized horizon label sequence
p_utrgv <- c('^A',
             'B2$',
             'B2ca')
utrgv$genhz <- generalize.hz(utrgv$hzdesg, n_utrgv, p_utrgv) # generate labels

# Visually inspect assignment and check that everything looks right
tab_utrgv <- table(utrgv$genhz, utrgv$hzdesg)
addmargins(tab_utrgv)

plotSPC(utrgv, color='genhz')

utrgv_corr <- horizons(utrgv)

# 10 - UConn ####
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

# Save CSV to send to Ekundayo
write_csv(uconn_inv_hz, here("data_processed", "uconn_missing_depths.csv"))

uconn_valid <- horizons(uconn) %>%
  filter(!dsp_pedon_id %in% uconn_invalid$dsp_pedon_id)

depths(uconn_valid) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(uconn_valid) <- 'hzdesg'

checkHzDepthLogic(uconn_valid)
plotSPC(uconn_valid, color="soc_pct") # Some of these profiles are shallow, but some are really deep.

# Add generalized horizon labels
# Sequence: A, Bw, BC, C
n_uconn <- c('A', 'Bw', 'BC', 'C') # generalized horizon label sequence
p_uconn <- c('A',
             '^Bw|Bw$',
             '^BC',
             '^C')
uconn_valid$genhz <- generalize.hz(uconn_valid$hzdesg, n_uconn, p_uconn) # generate labels

# Visually inspect assignment and check that everything looks right
tab_uconn <- table(uconn_valid$genhz, uconn_valid$hzdesg)
addmargins(tab_uconn)

plotSPC(uconn_valid, color='genhz')

uconn_corr <- horizons(uconn_valid)

# 11 - Re-join data for corrected project dataset and write CSV ####
coop_corr <- bind_rows(ks_corr, ncs_corr, tam1_corr, tam2_hz, wash_corr, minn_corr, osu_corr, ill_corr, utrgv_corr, uconn_corr) %>%
  select(!hzID)
coop_corr_spc <- coop_corr

# Promote to SPC object
depths(coop_corr_spc) <- dsp_pedon_id ~ hrzdep_t + hrzdep_b
hzdesgnname(coop_corr_spc) <- 'hzdesg'

checkHzDepthLogic(coop_corr_spc)

# Write CSV
write_csv(coop_corr, here("data_processed", "02_coop_data_horizons_valid.csv"), na="NA")
