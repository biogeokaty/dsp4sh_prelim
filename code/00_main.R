# 00 - Main file for DSP4SH analysis
# Katy Dynarski, October 2023

# 0 - Libraries ####
library(here)
library(tidyverse)
library(DBI)
library(RSQLite)
library(janitor)
library(aqp)
library(zoo)
library(flextable)
library(multcompView)

# 1 - ggplot theme function ####
theme_katy <- function(base_size=14) {
  theme_classic(base_size=base_size) %+replace%
    theme(# Legend
      legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(.8)),
      #axes                      
      axis.text=element_text(size=rel(.8)),
      axis.title=element_text(size=rel(1)))
}

# 2 - Function to make significance letters
hsd_letters<-function(x, y){
  aov_out <- aov(y~x)
  hsd <- TukeyHSD(aov_out)
  hsd_frame <- data.frame(hsd["x"])
  sig <- hsd_frame %>%
    select(contains("p.adj")) %>%
    rename(p.adj = contains("p.adj")) %>%
    na.omit() %>%
    rownames_to_column(var="comparison") %>%
    mutate(is_sig = ifelse(p.adj<0.05, TRUE, FALSE))
  
  mult <- sig[,3]
  names(mult) <- sig$comparison
  
  data.frame(multcompLetters(mult)$Letters) %>%
    rownames_to_column(var="x") %>%
    rename(letters = multcompLetters.mult..Letters)
  
}
