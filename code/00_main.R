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
library(gt)
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
library(raster)
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
library(tidyverse)

# 1 - ggplot theme function ####
theme_katy <- function(base_size=14) {
  theme_minimal(base_size=base_size) %+replace%
    theme(# Legend
      legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(.8)),
      #axes                      
      axis.text=element_text(size=rel(.8)),
      axis.title=element_text(size=rel(1)))
}

theme_katy_grid <- function(base_size=14) {
  theme_minimal(base_size=base_size) %+replace%
    theme(# Legend
      legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(.8)),
      #axes                      
      axis.text=element_text(size=rel(.8)),
      axis.title=element_text(size=rel(1)),
      # set margins (necessary for figure to not awkwardly overlap when using plot_grid in cowplot)
      plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
}