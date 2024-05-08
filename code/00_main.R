# 00 - Main file for DSP4SH analysis
# Katy Dynarski, October 2023

# 0 - Libraries ####
library(here)
library(DBI)
library(RSQLite)
library(janitor)
library(aqp)
library(zoo)
library(flextable)
library(multcompView)
library(multcomp)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(lme4)
library(ggeffects)
library(partykit)
library(vip)
library(leaps)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggeasy)
library(viridis)
library(cowplot)
library(missMDA)
library(ggfortify)
library(scales)
library(metafor)
library(ggpubr)
library(glmulti)
library(sensemakr)
library(tidyverse)

# 1 - ggplot theme function ----
theme_katy <- function(base_size=14) {
  theme_bw(base_size=base_size) %+replace%
    theme(# Legend
      legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(.8)),
      #axes                      
      axis.text=element_text(size=rel(.8)),
      axis.title=element_text(size=rel(1)))
}

theme_katy_grid <- function(base_size=14) {
  theme_bw(base_size=base_size) %+replace%
    theme(# Legend
      legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(.8)),
      #axes                      
      axis.text=element_text(size=rel(.8)),
      axis.title=element_text(size=rel(1)),
      # set margins (necessary for figure to not awkwardly overlap when using plot_grid in cowplot)
      plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
}

# 2 - Make vector of indicator labels so they will print nicely ----
indicator_labs <- c("soc_pct" = "SOC %",
                    "soc_stock_100cm" = "SOC stock (100 cm depth)",
                    "soc_stock_0_30cm" = "SOC stock (30 cm depth)",
                    "bglucosaminidase" = "B-glucosaminidase (NAG)",
                    "ace" = "ACE protein",
                    "bglucosidase" = "B-glucosidase (BG)",
                    "kssl_wsa" = "Aggregate stability (KSSL)",
                    "yoder_agg_stab_mwd" = "Aggregate stability (Yoder MWD)",
                    "arylsulfatase" = "Arylsulfatase",
                    "pox_c" = "POX-C",
                    "tn_pct" = "Total N%",
                    "bulk_density" = "Bulk density",
                    "soil_respiration" = "Respiration",
                    "phosphodiesterase" = "Phosphodiesterase",
                    "alkaline_phosphatase" = "Alkaline phosphatase",
                    "acid_phosphatase" = "Acid phosphatase",
                    "p_h" = "pH")

# Also make into dataframe (need this to use labels within map functions)
indicator_labs_df <- data.frame(indicator_labs) %>%
  rownames_to_column() %>%
  rename(indicator = rowname,
         label = indicator_labs)

# 3 - Run all the scripts :) ----
