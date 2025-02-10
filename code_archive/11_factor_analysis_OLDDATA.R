# 11 - Factor Analysis 
# Katy Dynarski, March 2024

# 0 - Import data ----
# Data needed for factor analysis is similar to what is needed for PCA - want just the soil health indicators (can even work with just a correlation matrix)

meta_df <- read.csv(here("data_processed", "05_meta_df.csv"))

# 1 - Exploratory Factor Analysis ----
# First step is to prepare data - want only the indicators
surf_ind <- meta_df %>%
select(bulk_density, soc_pct, tn_pct:yoder_agg_stab_mwd, p_h:ace, clay_pct_field)

# Get basic information about each indicator - how many observations are there, and how many NA values are there
surf_ind %>% map_df(~(data.frame(n_distinct = n_distinct(.x),
                                  class = class(.x),
                                 sum(is.na(.x)))),
                     .id = "variable")
# There are 262 observations, should drop any indicators with more than 100 NA values
# Drop tn_pct, ph, phosphodiesterase, and clay_pct - they all have more than 100 observations missing which is going to mess up the factor analysis

surf_ind_for_fa <- surf_ind %>%
  select(-tn_pct, -p_h, -phosphodiesterase, -clay_pct_field)

# Scale indicator data
surf_ind_scale <- scale(surf_ind_for_fa)

# Make correlation matrix to feed into factor analysis
corr_matrix <- cor(surf_ind_scale, use="pairwise.complete.obs")

# Next step is to determine the number of factors to keep
ev <- eigen(corr_matrix, symmetric=TRUE)
ap <- parallel(subject=267,var=12, rep=100,cent=.05) #subject is the number of observations, var is the number of variables
# rep and cent you can keep the same no matter what
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#  Most tests indicate that three factors is the optimal number

# Now, conduct the factor analysis
efa <- fa(corr_matrix, 3, rotate="oblimin", fm="ml")
efa
fa.diagram(efa$loadings)

# 2 - Structural Equation Modeling - do factors describe SOC stocks? ----
# Need a correlation matrix that includes SOC stocks
# Can pull this from the "meta_df" dataframe
# Variables in covariance matrix need to match the order in the factor

# Examine one factor at a time - the first factor includes enzyme activity
# want df with enzyme activities and SOC stocks (use 30 cm stocks because there are more of them)
sem_df1 <- meta_df %>%
select(bglucosaminidase, alkaline_phosphatase, arylsulfatase, bglucosidase, acid_phosphatase, soc_stock_30cm)
sem_scale1 <- scale(sem_df1) %>% na.omit()

# Make correlation matrix 
sem_corr_matrix1 <- cor(sem_scale1, use="pairwise.complete.obs")

# Save equation to be tested as an object
sem_eq1 <- '
F1 =~ bglucosaminidase + alkaline_phosphatase + arylsulfatase + bglucosidase + acid_phosphatase
soc_stock_30cm ~ F1
'

# conduct SEM
sem_fit1 <- sem(sem_eq1, sample.cov = sem_corr_matrix1, sample.nobs = 172, std.lv = TRUE)
summary(sem_fit1, fit.measures = TRUE, standardized = TRUE)
# How to check if the fit is good:
# Looking for comparative fit index (CFI) and Tucker Lewis Index (TLI) >0.9 (ideally >0.95)
# SRMR (standardized root mean square residual) <0.1

# What did we get?
# CFI: 0.928
# TLI: 0.880
# SRMR: 0.059
# Model estimate: -0.156 (coefficient), p=0.847
# Not a significant relationship between enzyme activities and SOC stocks

# Factor 2 is soc, ace, wsa, and pox-c
sem_df2 <- meta_df %>%
  select(soc_pct, ace, kssl_wsa, pox_c, soc_stock_30cm)
sem_scale2 <- scale(sem_df2) %>% na.omit()

# Make correlation matrix 
sem_corr_matrix2 <- cor(sem_scale2, use="pairwise.complete.obs")

# Save equation to be tested as an object
sem_eq2 <- '
F2 =~ soc_pct + ace + kssl_wsa + pox_c
soc_stock_30cm ~ F2
'

# conduct SEM
sem_fit2 <- sem(sem_eq2, sample.cov = sem_corr_matrix2, sample.nobs = 192, std.lv = TRUE)
summary(sem_fit2, fit.measures = TRUE, standardized = TRUE)
# How to check if the fit is good:
# Looking for comparative fit index (CFI) and Tucker Lewis Index (TLI) >0.9 (ideally >0.95)
# SRMR (standardized root mean square residual) <0.1

# What did we get?
# CFI: 0.960
# TLI: 0.919
# SRMR: 0.033
# Model estimate: 0.577 (coefficient), p=0.000
# this is a good fit! Not entirely surprising that the factor with carbon-related traits would explain soil carbon stocks, but still!

# Factor 3: soil respiration and Yoder agg
sem_df3 <- meta_df %>%
  select(soil_respiration, yoder_agg_stab_mwd, soc_stock_30cm)
sem_scale3 <- scale(sem_df3) %>% na.omit()

# Make correlation matrix 
sem_corr_matrix3 <- cor(sem_scale3, use="pairwise.complete.obs")

# Save equation to be tested as an object
sem_eq3 <- '
F3 =~ soil_respiration + yoder_agg_stab_mwd
soc_stock_30cm ~ F3
'

# conduct SEM
sem_fit3 <- sem(sem_eq3, sample.cov = sem_corr_matrix3, sample.nobs = 216, std.lv = TRUE)
summary(sem_fit3, fit.measures = TRUE, standardized = TRUE)
# How to check if the fit is good:
# Looking for comparative fit index (CFI) and Tucker Lewis Index (TLI) >0.9 (ideally >0.95)
# SRMR (standardized root mean square residual) <0.1

# For this model, get an error message that some estimated variances are negative. Usually this means that there is something "wrong" with the model structure - might mean that there isn't enough covariance between the variables to make the model work