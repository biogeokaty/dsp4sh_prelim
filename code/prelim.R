# Preliminary data exploration of DSP4SH data v2
# Katy Dynarski, October 2023

#### Housekeeping  - Libraries and Data Import ####
library(here)
library(tidyverse)

pedons <- read.csv(here("data", "DSP4SH_Pedons_mgmt_phase1-2v2.csv"))
