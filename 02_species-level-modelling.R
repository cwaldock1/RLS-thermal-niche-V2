# Second script in thermal-niche RLS analysis for species level analysis.
# This section models all the data in qgams and occupancy models and produces the dataset of realized thermal niche estimates for later use. 

# This includes 
# 1. Qgam analysis of abundance
# 2. Occupancy models for estimating niche limits
# 3. Derivation of thermal niche dataset
# 4. Plot of TPCs. 

# Initiated 27/03/2018
# Author: Conor Waldock


# PREAMBLE ----
# Load data from 01_organising-data.R ----
load(file = 'data_derived/01_organising-data_SAVE-IMAGE.RData') 

# Load libraries and packages ----

# Detach packages from script 1. 
detach("package:geosphere", unload=TRUE)
detach("package:rgeos", unload=TRUE)

# Managing and manipulating data. 
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

# Fitting quantile gams. 
library(qgam)

# Fitting occupancy models. 
library(glmmTMB)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MODELLING QGAMS SECTION ----
# Function which fits quantile gams and extracts parameters of interest. ----
FitQuantileGam <- function(TestSpp, k = 4, q = 0.9){
  
  NrowPresences <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences <- TestSpp %>% filter(AbundanceAdult40 == 0) %>% nrow(.)
  
  if(NrowPresences <= NrowAbsences){
    TestSpp <- rbind(TestSpp %>% filter(AbundanceAdult40 > 0), TestSpp %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
  }else{NULL}
  
  # print(nrow(TestSpp))
  
  # Fit generalized additive model
  #Create columns for modelling
  TestSpp$AbundanceAdult40_log <- log(TestSpp$AbundanceAdult40 + 1)
  TestSpp$SamplingIntensity_Scaled <- scale(TestSpp$SamplingIntensity)
  TestSpp$MeanSiteSST_NOAA_Scaled <- scale(TestSpp$MeanSiteSST_NOAA)
  TestSpp$MeanSiteSST_NOAA_scaleBT <- attr(scale(TestSpp$MeanSiteSST_NOAA), "scaled:scale")
  TestSpp$MeanSiteSST_NOAA_centerBT <- attr(scale(TestSpp$MeanSiteSST_NOAA), "scaled:center")
  
  if(is.numeric(k)){
    TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled, k = k) + SamplingIntensity_Scaled, qu = q, data = TestSpp), error = function(e) NA)
    
    i <- 1
    while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
      class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled, k = k) + SamplingIntensity_Scaled, qu = q, data = TestSpp), error = function(e) NA)
      i <-  i + 1}
    
    while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled, k = k) + SamplingIntensity_Scaled, qu = q, data = TestSpp)}
  }else{
    
    TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled) + SamplingIntensity_Scaled, qu = q, data = TestSpp), error = function(e) NA)
    
    i <- 1
    while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
      class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled) + SamplingIntensity_Scaled, qu = q, data = TestSpp), error = function(e) NA)
      i <-  i + 1}
    
    while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanSiteSST_NOAA_Scaled) + SamplingIntensity_Scaled, qu = q, data = TestSpp)}
  }
  
  
  # Extract temperture at maximum abundance, and maximum abundance at optimum temperature. 
  Temps <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanSiteSST_NOAA_Scaled), MaxTemp = max(.$MeanSiteSST_NOAA_Scaled)) %>% unnest()
  TempData <- data.frame(MeanSiteSST_NOAA_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 50))
  TempData$SamplingIntensity_Scaled <- mean(TestSpp$SamplingIntensity_Scaled, na.rm = T)
  AbunPredictions <- predict(TestSpp_gam11, TempData, se = T)
  TempData$Abundance <- AbunPredictions$fit
  TempData$SE <- AbunPredictions$se
  TempData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  TempData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  
  # Estimate optimum temperature based on maximum abundance predicted 
  Topt <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanSiteSST_NOAA_Scaled'] * TestSpp$MeanSiteSST_NOAA_scaleBT[1]) + TestSpp$MeanSiteSST_NOAA_centerBT[1]
  
  # Estimate variation in temperature below thermal optimum. OLD METHOD. 
  # Tsd <- TestSpp %>% filter(MeanSiteSST_NOAA < Topt) %>% filter(AbundanceAdult40 > 0) %>% .$MeanSiteSST_NOAA %>% sd(., na.rm = T)
  # Tsd_V2 <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% .$MeanSiteSST_NOAA %>% sd(., na.rm = T)/2
  
  ### --- Estimate max abundance
  MaxAbundance <- TestSpp %>% .$MaxAbundance %>% unique
  
  ### --- Estimate confidence score criteria 4. 
  T_Opt_Difference_Upper <- max(TestSpp$MeanSiteSST_NOAA) - Topt 
  T_Opt_Difference_Lower <- Topt - min(TestSpp$MeanSiteSST_NOAA)
  
  
  ### --- Estimate confidence score criteria 5. 
  # Example extraction
  T_Gam_pvalue <- summary(TestSpp_gam11)$s.table[4]# < 0.05
  T_Gam_edf <- summary(TestSpp_gam11)$s.table[1]# > 1
  
  ### --- Create output dataframe of predictions from model to save with data (and plot up later)
  PredData <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% 
    do(SpeciesName = unique(.$SpeciesName), 
       MeanSiteSST_NOAA        = seq(min(.$MeanSiteSST_NOAA), max(.$MeanSiteSST_NOAA), length.out = 100), 
       MeanSiteSST_NOAA_Scaled = seq(min(.$MeanSiteSST_NOAA_Scaled), max(.$MeanSiteSST_NOAA_Scaled), length.out = 100)) %>% 
    unnest(SpeciesName) %>% unnest(MeanSiteSST_NOAA, MeanSiteSST_NOAA_Scaled)
  
  PredData$SamplingIntensity_Scaled <- median(TestSpp$SamplingIntensity_Scaled, na.rm = T)
  AbunPredictions <- predict(TestSpp_gam11, PredData, se = T)
  PredData$Abundance <- AbunPredictions$fit
  PredData$SE  <- AbunPredictions$se
  PredData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  PredData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  PredData$Topt <- Topt
  PredData$MaxAbundance <- MaxAbundance
  PredData$q <- q
  PredData$k <- k
  
  
  OutputData <- list(model = TestSpp_gam11, predictions = PredData, rawdata = TestSpp %>% select(SpeciesName, AbundanceAdult40_log,  AbundanceAdult40, MeanSiteSST_NOAA))
  
  ### --- Save qgam outputs for if needed later
  saveRDS(OutputData, file = paste('data_derived/qgam_outputs/', gsub(' ' , '_', unique(TestSpp$SpeciesName)), gsub('-','_',Sys.Date()), 'q_is', gsub('0.','',q),'k_is',k,'.rds', sep = '_'))
  rm(TestSpp_gam11)
  
  return(data_frame(SpeciesName = unique(TestSpp$SpeciesName), 
                    Topt = Topt, 
                    MaxAbundance = MaxAbundance, 
                    T_Opt_Difference_Upper = T_Opt_Difference_Upper,
                    T_Opt_Difference_Lower = T_Opt_Difference_Lower,
                    T_Gam_pvalue = T_Gam_pvalue, 
                    T_Gam_edf = T_Gam_edf
  ))
  
}
# Fit the qgam models to a quantile of 0.8 -----

# Fit qgams in dplyr. Each individual model fits independently 
All_QGams <- RLS_19 %>% group_by(SpeciesName) %>% do(GamThermalNiche = tryCatch(FitQuantileGam(., q = 0.8), error = function(e) NA))

# Check for NAs. 
NA_Models <- which(is.na(All_QGams$GamThermalNiche)) 

# Save the output object. 
saveRDS(All_QGams, file = 'data_derived/All_QGams_2018-02-27.rds')

# Obtain parameters from the quantile gam models. 
Quantile_Parameters <- do.call(rbind, All_QGams$GamThermalNiche)

# How many have non-significant estimates of Topt
Quantile_Parameters[which(Quantile_Parameters$T_Gam_pvalue > 0.05),]
# Step 1. Load object created from the above code and run from 166:255 below ----
All_QGams <- readRDS(file = 'data_derived/All_QGams_2018-02-27.rds')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# DEFINING THERMAL GUILDS ----
# Define thermal guilds based on Topt from the qgam models. Creation of RLS_All, RLS_Trop, and RLS_Temp objects ----

# Redefine thermal guild based on qgam estimates of thermal niche. 23°C from Stuart-Smith et al. 2015 Nature, and Stuart-Smith et al. 2017 Nature Ecology and Evolution. 
ThermalGuilds <- left_join(RLS_19, Quantile_Parameters %>% select(SpeciesName, Topt)) %>% 
  group_by(SpeciesName) %>% 
  do(ThermalGuild = ifelse(.$Topt > 23, 'tropical', 'temperate')) %>% 
  unnest(ThermalGuild) %>% unique()

# Estimate mean temperature across all observations used in future. 
MeanTemperatures <- RLS_19 %>% 
  group_by(SpeciesName) %>% 
  filter(AbundanceAdult40 > 0) %>%  
  do(T_Mean_Obs = mean(.$MeanSiteSST_NOAA, na.rm = T)) %>% 
  unnest(T_Mean_Obs)

# Join in thermal guild character, and mean temperature across ranges. 
RLS_19 <- left_join(RLS_19, ThermalGuilds)
RLS_19 <- left_join(RLS_19, MeanTemperatures)

# Estimate observed thermal midpoints (same as T_Mean_Obs), but also observed upper and lower per species. 
ThermalNicheData_Obs <- RLS_19 %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanSiteSST_NOAA)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanSiteSST_NOAA)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanSiteSST_NOAA, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
  select(-data)

# Define tropical and temperate data. 
RLS_Trop <- RLS_19 %>% filter(ThermalGuild == 'tropical')
RLS_Temp <- RLS_19 %>% filter(ThermalGuild == 'temperate')

# Scale covarites for RLS_Trop and RLS_Temp and create RLS_All ----

# Note that reef area is already log transformed. 

# Scale and centre covariates for tropical species.  
RLS_Trop$MeanSiteSST_NOAA_Scaled <- scale(RLS_Trop$MeanSiteSST_NOAA)
RLS_Trop$HumPop50km2_Scaled <- scale(RLS_Trop$HumPop50km2)
RLS_Trop$npp_mean_Scaled <- scale(RLS_Trop$npp_mean)
RLS_Trop$SamplingIntensity_Scaled <- scale(RLS_Trop$SamplingIntensity)

RLS_Trop$MeanSiteSST_NOAA_centerBT <- attr(scale(RLS_Trop$MeanSiteSST_NOAA), 'scaled:center')
RLS_Trop$HumPop50km2_centerBT <- attr(scale(RLS_Trop$HumPop50km2), 'scaled:center')
RLS_Trop$npp_mean_centerBT <- attr(scale(RLS_Trop$npp_mean), 'scaled:center')
RLS_Trop$SamplingIntensity_centerBT <- attr(scale(RLS_Trop$SamplingIntensity), 'scaled:center')

RLS_Trop$MeanSiteSST_NOAA_scaleBT <- attr(scale(RLS_Trop$MeanSiteSST_NOAA), 'scaled:scale')
RLS_Trop$HumPop50km2_scaleBT <- attr(scale(RLS_Trop$HumPop50km2), 'scaled:scale')
RLS_Trop$npp_mean_scaleBT <- attr(scale(RLS_Trop$npp_mean), 'scaled:scale')
RLS_Trop$SamplingIntensity_scaleBT <- attr(scale(RLS_Trop$SamplingIntensity), 'scaled:scale')

# Scale and centre covariates for temperate species. 
RLS_Temp$MeanSiteSST_NOAA_Scaled <- scale(RLS_Temp$MeanSiteSST_NOAA)
RLS_Temp$HumPop50km2_Scaled <- scale(RLS_Temp$HumPop50km2)
RLS_Temp$npp_mean_Scaled <- scale(RLS_Temp$npp_mean)
RLS_Temp$SamplingIntensity_Scaled <- scale(RLS_Temp$SamplingIntensity)

RLS_Temp$MeanSiteSST_NOAA_centerBT <- attr(scale(RLS_Temp$MeanSiteSST_NOAA), 'scaled:center')
RLS_Temp$HumPop50km2_centerBT <- attr(scale(RLS_Temp$HumPop50km2), 'scaled:center')
RLS_Temp$npp_mean_centerBT <- attr(scale(RLS_Temp$npp_mean), 'scaled:center')
RLS_Temp$SamplingIntensity_centerBT <- attr(scale(RLS_Temp$SamplingIntensity), 'scaled:center')

RLS_Temp$MeanSiteSST_NOAA_scaleBT <- attr(scale(RLS_Temp$MeanSiteSST_NOAA), 'scaled:scale')
RLS_Temp$HumPop50km2_scaleBT <- attr(scale(RLS_Temp$HumPop50km2), 'scaled:scale')
RLS_Temp$npp_mean_scaleBT <- attr(scale(RLS_Temp$npp_mean), 'scaled:scale')
RLS_Temp$SamplingIntensity_scaleBT <- attr(scale(RLS_Temp$SamplingIntensity), 'scaled:scale')

# Check for NAs. 
sapply(RLS_Temp, function(x) sum(is.na(x)))
sapply(RLS_Trop, function(x) sum(is.na(x)))

# Save objects
saveRDS(RLS_Temp, file = 'data_derived/RLS_withCovariates_Temperate_2018-03-28.rds')
saveRDS(RLS_Trop, file = 'data_derived/RLS_withCovariates_Tropical_2018-03-28.rds')

# Create RLS_All and save.
RLS_All <- rbind(RLS_Temp, RLS_Trop)
saveRDS(RLS_All, file = 'data_derived/RLS_withCovariates_All_2018-03-28.rds')







# (step 2) Load objects from above code if needed ----
#RLS_Temp <- readRDS(file = 'data_derived/RLS_withCovariates_Temperate_2018-03-28.rds')
#RLS_Trop <- readRDS(file = 'data_derived/RLS_withCovariates_Tropical_2018-03-28.rds')
#RLS_All  <- readRDS(file = 'data_derived/RLS_withCovariates_All_2018-03-28.rds')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# THERMAL NICHE DATA OBJECT CREATION used throughout to parameterise thermal performance curves ----
# Create thermal niche data object ----

# Select species names and thermal guilds to begin with. 
ThermalNicheData_New <- RLS_All %>% select(SpeciesName, ThermalGuild) %>% unique()

# Merge with Topts created from qgam models 
ThermalNicheData_New <- left_join(ThermalNicheData_New, Quantile_Parameters %>% select(SpeciesName, Topt, MaxAbundance))

# Create object of thermal niche from observations only (for later merge, some redundancy here)
ThermalNicheData_Obs <- RLS_All %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanSiteSST_NOAA)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanSiteSST_NOAA)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanSiteSST_NOAA, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
  select(-data)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MODELLING OCCUPANCY SECTION ----
# Create data for occupancy models  ----

# Remove species for which Tupper and Tlower will be poorly estimated and save in an occupancy model.
# This confidence criteria is defined in script 01_organising-data.R
# Temperate
RLS_Temp_Occ_Tupper <- RLS_Temp %>% filter(Confidence_Occ_Tupper == 1); length(unique(RLS_Temp_Occ_Tupper$SpeciesName)) # n = 169
RLS_Temp_Occ_Tlower <- RLS_Temp %>% filter(Confidence_Occ_Tlower == 1); length(unique(RLS_Temp_Occ_Tlower$SpeciesName)) # n = 147

# Tropical
RLS_Trop_Occ_Tupper <- RLS_Trop %>% filter(Confidence_Occ_Tupper == 1); length(unique(RLS_Trop_Occ_Tupper$SpeciesName)) # n = 328
RLS_Trop_Occ_Tlower <- RLS_Trop %>% filter(Confidence_Occ_Tlower == 1); length(unique(RLS_Trop_Occ_Tlower$SpeciesName)) # n = 494

# Save all for future use. 
# save(RLS_Temp_Occ_Tupper, RLS_Temp_Occ_Tlower, RLS_Trop_Occ_Tupper, RLS_Trop_Occ_Tlower, file = 'data_derived/OccupancyModelData_UpperLowerConfidences_2018-02-28.RData')

# If not loaded already, read in data for occupancy models. 
# load(file = 'data/OccupancyModelData_UpperLowerConfidences_2018-02-28.RData')

# Data to run models for based on only samples above and below T_Mean_Obs. 
# This is not Topt because some Topts are close to sampling limits and so this threshold would provide poor data in occupancy models. 
RLS_Temp_Occ_Tupper2 <- RLS_Temp_Occ_Tupper %>% filter(MeanSiteSST_NOAA > T_Mean_Obs)
RLS_Temp_Occ_Tlower2 <- RLS_Temp_Occ_Tlower %>% filter(MeanSiteSST_NOAA < T_Mean_Obs)
RLS_Trop_Occ_Tupper2 <- RLS_Trop_Occ_Tupper %>% filter(MeanSiteSST_NOAA > T_Mean_Obs)
RLS_Trop_Occ_Tlower2 <- RLS_Trop_Occ_Tlower %>% filter(MeanSiteSST_NOAA < T_Mean_Obs)

# Create functions to handle outputs from occupancy models ----
# Extract global parameter estimates from RLS glmmTMB models looking at temperature
ExtractGlobalTrend <- function(Data, Model, FixedFormula, Error = 'Poisson'){
  
  if(sum(grepl('ThermalGuild', FixedFormula)) < 1){
    SeqTemp <- seq(from = min(Data$MeanSiteSST_NOAA_Scaled), 
                   to = max(Data$MeanSiteSST_NOAA_Scaled), 
                   length = 100)
    
    MyData <- expand.grid(MeanSiteSST_NOAA_Scaled = SeqTemp,
                          HumPop50km2_Scaled = mean(Data$HumPop50km2_Scaled), 
                          npp_mean_Scaled = mean(Data$npp_mean_Scaled),
                          ReefAreaIn15km = 0,
                          ThermalGuild = unique(Data$ThermalGuild))
  }else{
    
    Data$ThermalGuild <- as.factor(Data$ThermalGuild)
    
    MyData <- ddply(Data, 
                    .(ThermalGuild), 
                    summarize,
                    MeanSiteSST_NOAA_Scaled = seq(min(MeanSiteSST_NOAA_Scaled), max(MeanSiteSST_NOAA_Scaled), length = 100))
    
    MyData$HumPop50km2_Scaled <- mean(Data$HumPop50km2_Scaled)
    MyData$npp_mean_Scaled <- mean(Data$npp_mean_Scaled)
    MyData$ReefAreaIn15km = 0
    
  }
  
  
  
  # Create fake relationships
  X  <- model.matrix(FixedFormula, data = MyData)
  
  if(Error == 'binomial'){
    
    MyData$eta <- X %*% fixef(Model)$con
    MyData$mu  <- exp(MyData$eta) / (1 + exp(MyData$eta))
    MyData$SE <- sqrt(diag(X %*% vcov(Model)$con %*% t(X)))
    MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$SE) / (1 + exp(MyData$eta  + 1.96 *MyData$SE))
    MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$SE) / 
      (1 + exp(MyData$eta  - 1.96 *MyData$SE))
    
  }else{NULL}
  
  if(Error == 'Poisson'){
    
    MyData$eta <- X %*% fixef(Model)$con
    MyData$mu <- exp(X %*% fixef(Model)$con)
    MyData$SE <- sqrt(diag(X %*% vcov(Model)$con %*% t(X)))
    MyData$SeUp  <- exp(MyData$eta + 1.96 * MyData$SE) 
    MyData$SeLo  <- exp(MyData$eta - 1.96 * MyData$SE)
    
  }else{NULL}
  
  if(sum(grepl('ThermalGuild', FixedFormula)) < 1){
    
    MyData$ScaledEstimate <- MyData$mu / max(MyData$mu)
    
  }else{
    
    ScaledEstimate <- ddply(MyData, 
                            .(ThermalGuild), 
                            summarize,
                            ScaledEstimate = mu / max(mu))
    MyData$ScaledEstimate <- ScaledEstimate$ScaledEstimate
  }
  
  return(MyData)
  
}


# Extract adjusted random parameter estimates from RLS glmmTMB models looking at temperature.
CreateRandomEffectData <- function(Model, Data, NRandomSlopes, NFixedEffects, 
                                   FixedFormula = formula(~MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +  I(MeanSiteSST_NOAA_Scaled^3) + HumPop50km2_Scaled + npp_mean_Scaled), 
                                   Error = 'Poisson', 
                                   RandomSlopes = c('MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)', 'I(MeanSiteSST_NOAA_Scaled^3)')){
  
  Data$ThermalGuild <- as.factor(Data$ThermalGuild)
  
  MyData_2 <- plyr::ddply(Data, 
                          .(SpeciesName, ThermalGuild), 
                          plyr::summarize,
                          MeanSiteSST_NOAA_Scaled = seq(min(MeanSiteSST_NOAA_Scaled), max(MeanSiteSST_NOAA_Scaled), length = 100))
  
  MyData_2$HumPop50km2_Scaled <- mean(Data$HumPop50km2_Scaled)
  MyData_2$npp_mean_Scaled <- mean(Data$npp_mean_Scaled)
  MyData_2$ReefAreaIn15km = 0
  
  # Do for one species
  SpeciesParameters <- list()
  GlobalFixedEffects <- fixef(Model)$con
  
  for(i in 1:length(unique(Data$SpeciesName))){
    SpeciesParameters[[i]] <- fixef(Model)$con 
    SpeciesParameters[[i]][c('(Intercept)', RandomSlopes)] <- fixef(Model)$con[c('(Intercept)', RandomSlopes)] + ranef(Model)$con$SpeciesName[i,c('(Intercept)', RandomSlopes)]
    #SpeciesParameters[[i]] <- c(fixef(Model)$con[c('(Intercept)', RandomSlopes)] + ranef(Model)$con$SpeciesName[i,c('(Intercept)', RandomSlopes)], fixef(Model)$con)
    SpeciesParameters[[i]] <- data.frame(SpeciesParameters[[i]])
    #colnames(SpeciesParameters[[i]]) <- c('(Intercept)', RandomSlopes)
    colnames(SpeciesParameters[[i]]) <- names(fixef(Model)$con)
  }
  
  if(Error == 'Poisson'){
    MyDataEstimates <- list()
    for(i in 1:length(unique(Data$SpeciesName))){
      # Create dataframe
      MyDataEstimates[[i]] <- data.frame(MyData_2 %>% filter(SpeciesName == unique(Data$SpeciesName)[i]))
      
      # Create model matrix for each species
      X <- model.matrix(FixedFormula, data = MyDataEstimates[[i]])
      
      MyDataEstimates[[i]]$eta <- X %*% as.numeric(SpeciesParameters[[i]])
      #ranef(ZIP_OLRE_M1)[[1]]$SpeciesName
      MyDataEstimates[[i]]$mu <- exp(X %*% as.numeric(SpeciesParameters[[i]]))
      #MyDataEstimates[[i]]$SE <- sqrt(diag(X %*% vcov(ZIP_OLRE_M1)$con %*% t(X))) ### Unsure how to estimate uncertainty around the random effect slopes?!
      #MyDataEstimates[[i]]$SeUp  <- exp(MyData$eta + 1.96 *MyData$SE) 
      #MyDataEstimates[[i]]$SeLo  <- exp(MyData$eta - 1.96 *MyData$SE) 
    }}else{NULL}
  
  if(Error == 'binomial'){
    MyDataEstimates <- list()
    for(i in 1:length(unique(Data$SpeciesName))){
      # Create dataframe
      MyDataEstimates[[i]] <- data.frame(MyData_2 %>% filter(SpeciesName == unique(Data$SpeciesName)[i]))
      
      # Create model matrix for each species
      X <- model.matrix(FixedFormula, data = MyDataEstimates[[i]])
      MyDataEstimates[[i]]$eta <- X %*% as.numeric(SpeciesParameters[[i]])
      MyDataEstimates[[i]]$mu <- exp(X %*% as.numeric(SpeciesParameters[[i]])) / (1 + exp(X %*% as.numeric(SpeciesParameters[[i]])))
    }}else{NULL}
  
  MyDataEstimates <- lapply(MyDataEstimates, function(x){
    ScaledEstimate <- x$mu/max(x$mu) 
    return(data.frame(x, ScaledEstimate))
  })
  
  RandomEffectPlotData <- do.call(rbind, MyDataEstimates)
  
  RandomEffectPlotData$MeanSiteSST_NOAA <- (RandomEffectPlotData$MeanSiteSST_NOAA_Scaled *  Data$MeanSiteSST_NOAA_scaleBT[1]) + Data$MeanSiteSST_NOAA_centerBT[1]
  
  
  return(RandomEffectPlotData)
}


# Define occupancy model formulas ----
# Temperate Species occupancy model formula
Temp_Formula1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + 
                           (1 + MeanSiteSST_NOAA_Scaled| SpeciesName) + (1 | ECOregion/SiteCode))

Temp_Formula2 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2) + 
                           (1 + MeanSiteSST_NOAA_Scaled  | SpeciesName) + (1 | ECOregion/SiteCode))

Temp_Formula3 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2) +
                           (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))

# Tropical Species occupancy model formula
Trop_Formula1 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + 
                           (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))

Trop_Formula2 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) + 
                           (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled| SpeciesName) + (1 | ECOregion/SiteCode))

Trop_Formula3 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                           (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# OCCUPANCY TEMPERATE SPECIES ----
# Fit occupancy models TEMPERATE ----

### Temperate occupancy model and compare AIC between model approaches. 

# Fit model to uppers occupancy limits
Temp_Model_Upper1 <- glmmTMB(Temp_Formula1, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper1, file = 'data_derived/Temp_Model_Upper1_2018-01-03.rds')
Temp_Model_Upper2 <- glmmTMB(Temp_Formula2, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper2, file = 'data_derived/Temp_Model_Upper2_2018-01-03.rds')
Temp_Model_Upper3 <- glmmTMB(Temp_Formula3, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper3, file = 'data_derived/Temp_Model_Upper3_2018-01-03.rds')
AIC(Temp_Model_Upper1, Temp_Model_Upper2, Temp_Model_Upper3)

# Fit model to lowers occupancy limits
Temp_Model_Lower1 <- glmmTMB(Temp_Formula1, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower1, file = 'data_derived/Temp_Model_Lower1_2018-01-03.rds')
Temp_Model_Lower2 <- glmmTMB(Temp_Formula2, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower2, file = 'data_derived/Temp_Model_Lower2_2018-01-03.rds')
Temp_Model_Lower3 <- glmmTMB(Temp_Formula3, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower3, file = 'data_derived/Temp_Model_Lower3_2018-01-03.rds')
AIC(Temp_Model_Lower1, Temp_Model_Lower2, Temp_Model_Lower3)

# Read in occupancy models structure 1 and 3 for plots ----
Temp_Model_Upper1 <- readRDS(file = 'data_derived/Temp_Model_Upper1_2018-01-03.rds')
Temp_Model_Upper3 <- readRDS(file = 'data_derived/Temp_Model_Upper3_2018-01-03.rds')
Temp_Model_Lower1 <- readRDS(file = 'data_derived/Temp_Model_Lower1_2018-01-03.rds')
Temp_Model_Lower3 <- readRDS(file = 'data_derived/Temp_Model_Lower3_2018-01-03.rds')

# Creates plots and predictions from models -----

# Model list is just 1 and 3 for extraction of linear trend (a confidence criteria) and non-linear trend with non-linear random effect. 
Temp_Model_Upper_List <- list(Temp_Model_Upper1, Temp_Model_Upper3) # Create model list for functions to extract predictions
Temp_Model_Lower_List <- list(Temp_Model_Lower1, Temp_Model_Lower3) # Create model list for functions to extract predictions 

# Create list of objects for functions to plot outputs. 
# Extraction of global fixed effects. 
FixedFormula_Temp_List <- list(formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled), 
                               formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2)),
                               formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2)))

# Define random effects for input to function. 
RandomSlopes_Temp_List <- list(c('MeanSiteSST_NOAA_Scaled'), 
                               c('MeanSiteSST_NOAA_Scaled'), 
                               c('MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'))

# Create objects for plot data to be stored in. 
Temp_Model_Upper_Data <- list(); Temp_Model_Upper_RE <- list()
Temp_Model_Lower_Data <- list(); Temp_Model_Lower_RE <- list()

# Run functions for each model in loops. 
for(i in 1:length(unique(Temp_Model_Upper_List))){
  
  Temp_Model_Upper_Data[[i]]       <- ExtractGlobalTrend(Data = RLS_Temp_Occ_Tupper2, 
                                                         Model = Temp_Model_Upper_List[[i]], 
                                                         FixedFormula = FixedFormula_Temp_List[[i]], 
                                                         Error = 'binomial')
  Temp_Model_Upper_Data[[i]]$Model <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Temp_Model_Upper_List[[i]]))), 'AIC = ', round(AIC(Temp_Model_Upper_List[[i]])), i)[[3]]
  
  Temp_Model_Upper_RE[[i]]         <- CreateRandomEffectData(Data = RLS_Temp_Occ_Tupper2, 
                                                             Model = Temp_Model_Upper_List[[i]], 
                                                             FixedFormula = FixedFormula_Temp_List[[i]], 
                                                             RandomSlopes = RandomSlopes_Temp_List[[i]],
                                                             Error = 'binomial')
  Temp_Model_Upper_RE[[i]]$Model   <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Temp_Model_Upper_List[[i]]))), 'AIC = ', round(AIC(Temp_Model_Upper_List[[i]])), i)[[3]]
  
}
for(i in 1:length(unique(Temp_Model_Lower_List))){
  
  Temp_Model_Lower_Data[[i]]       <- ExtractGlobalTrend(Data = RLS_Temp_Occ_Tlower2, 
                                                         Model = Temp_Model_Lower_List[[i]], 
                                                         FixedFormula = FixedFormula_Temp_List[[i]], 
                                                         Error = 'binomial')
  Temp_Model_Lower_Data[[i]]$Model <-paste(as.character(gsub("(.{75})", "\\1 \n", formula(Temp_Model_Lower_List[[i]]))), 'AIC = ', round(AIC(Temp_Model_Lower_List[[i]])), i)[[3]]
  
  Temp_Model_Lower_RE[[i]]         <- CreateRandomEffectData(Data = RLS_Temp_Occ_Tlower2, 
                                                             Model = Temp_Model_Lower_List[[i]], 
                                                             FixedFormula = FixedFormula_Temp_List[[i]], 
                                                             RandomSlopes = RandomSlopes_Temp_List[[i]],
                                                             Error = 'binomial')
  Temp_Model_Lower_RE[[i]]$Model   <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Temp_Model_Lower_List[[i]]))), 'AIC = ', round(AIC(Temp_Model_Lower_List[[i]])), i)[[3]]
  
}

# Bind outputs from functions creating plot data together
Temp_Model_Upper_Data <- do.call(rbind, Temp_Model_Upper_Data)
Temp_Model_Upper_RE   <- do.call(rbind, Temp_Model_Upper_RE)
Temp_Model_Lower_Data <- do.call(rbind, Temp_Model_Lower_Data)
Temp_Model_Lower_RE   <- do.call(rbind, Temp_Model_Lower_RE)

# Plot uppers  
Temp_UpperPlot <- ggplot() + 
  #geom_line(data = Temp_Model_Upper_Data, aes(x = MeanSiteSST_NOAA_Scaled, y = mu)) + 
  geom_line(data = Temp_Model_Upper_RE, aes(x = MeanSiteSST_NOAA_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Temp_Model_Upper_Data, aes(label = Model, x = 0.5, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Plot lowers
Temp_LowerPlot <- ggplot() + 
  #geom_line(data = Temp_Model_Lower_Data, aes(x = MeanSiteSST_NOAA_Scaled, y = mu)) + 
  geom_line(data = Temp_Model_Lower_RE, aes(x = MeanSiteSST_NOAA_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Temp_Model_Lower_Data, aes(label = Model, x = -1, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Save in extra_figures
pdf(file = 'figures_extra/ThermalNicheLimits_TemperateOccModels.pdf')
Temp_LowerPlot
Temp_UpperPlot
dev.off()

# Model selection of fixed effects TEMPERATE ----

# Two term models 
Temp_M1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                # Remove polynomial temperature
Temp_M2 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled  + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                         # Remove all temperature
Temp_M3 <- formula(Presence ~ HumPop50km2_Scaled  + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2) + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode)) # Remove npp
Temp_M4 <- formula(Presence ~ npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2) + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))     # Remove human pop

# Single term models. 
Temp_M5 <- formula(Presence ~ MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2) + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode)) # Temperature only
Temp_M6 <- formula(Presence ~ HumPop50km2_Scaled + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                      # Human only
Temp_M7 <- formula(Presence ~ npp_mean_Scaled + (1 + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                         # npp only. 

# List of models to fit. 
Temp_Model_Selection <- list(Temp_M1, Temp_M2, Temp_M3, Temp_M4, Temp_M5, Temp_M6, Temp_M7)
 
# Empty list to fill 
Temp_Models_Upper_MS <- list()
Temp_Models_Lower_MS <- list()

# For loop to fit all models 
for(i in 1:length(Temp_Model_Selection)){
  print(paste('Model run', i))
  # Fit over each model. 
  Temp_Models_Upper_MS[[i]] <- glmmTMB(Temp_Model_Selection[[i]], family = "binomial", data = RLS_Temp_Occ_Tupper2)
  Temp_Models_Lower_MS[[i]] <- glmmTMB(Temp_Model_Selection[[i]], family = "binomial", data = RLS_Temp_Occ_Tlower2)
}

save(Temp_Models_Upper_MS, Temp_Models_Lower_MS, file = 'data_derived/Temperate-occupancy-model-selection.rdata')
  
# Compare AIC from model lists for upper limits 
Temp_Upper_MS <- c(list(Temp_Model_Upper3), Temp_Models_Upper_MS) # Compare with full model that has 2nd order polynomial too. 
bbmle::AICtab(Temp_Upper_MS)
AIC(Temp_Upper_MS[[1]])
AIC(Temp_Upper_MS[[4]])
FinalTempUpperModel <- Temp_Upper_MS[[4]] # Small difference in AIC. == 2. Parameter does not take away from the model or add. Retain for consistancy across all models. 

# Compare AIC from model lists for lower limits 
Temp_Lower_MS <- c(list(Temp_Model_Lower3), Temp_Models_Lower_MS) # Compare with full model that has 2nd order polynomial too. 
bbmle::AICtab(Temp_Lower_MS)
FinalTempLowerModel <- Temp_Lower_MS[[1]]

save(FinalTempUpperModel, FinalTempLowerModel, file = 'data_derived/Temperate-occupancy-final-models.rdata')

# load in model selections and final models from above code ----
load(file = 'data_derived/Temperate-occupancy-model-selection.rdata')
load(file = 'data_derived/Temperate-occupancy-final-models.rdata')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# OCCUPANCY TROPICAL SPECIES ---- 
# Fit occupancy models TROPICAL -----
### Tropical occupancy model and compare AIC between model approaches. 

# Fit model to uppers occupancy limits
Trop_Model_Upper1 <- glmmTMB(Trop_Formula1, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper1, file = 'data_derived/Trop_Model_Upper1_2018-01-03.rds')
Trop_Model_Upper2 <- glmmTMB(Trop_Formula2, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper2, file = 'data_derived/Trop_Model_Upper2_2018-01-03.rds')
Trop_Model_Upper3 <- glmmTMB(Trop_Formula3, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper3, file = 'data_derived/Trop_Model_Upper3_2018-01-03.rds')
AIC(Trop_Model_Upper1, Trop_Model_Upper2, Trop_Model_Upper3)

# Fit model to lowers occupancy limits
Trop_Model_Lower1 <- glmmTMB(Trop_Formula1, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower1, file = 'data_derived/Trop_Model_Lower1_2018-01-03.rds')
Trop_Model_Lower2 <- glmmTMB(Trop_Formula2, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower2, file = 'data_derived/Trop_Model_Lower2_2018-01-03.rds')
Trop_Model_Lower3 <- glmmTMB(Trop_Formula3, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower3, file = 'data_derived/Trop_Model_Lower3_2018-01-03.rds')
AIC(Trop_Model_Lower1, Trop_Model_Lower2, Trop_Model_Lower3)

# Read in occupancy models structure 1 and 3 for plots ----
Trop_Model_Upper1 <- readRDS(file = 'data_derived/Trop_Model_Upper1_2018-01-03.rds')
Trop_Model_Upper3 <- readRDS(file = 'data_derived/Trop_Model_Upper3_2018-01-03.rds')
Trop_Model_Lower1 <- readRDS(file = 'data_derived/Trop_Model_Lower1_2018-01-03.rds')
Trop_Model_Lower3 <- readRDS(file = 'data_derived/Trop_Model_Lower3_2018-01-03.rds')

# Creates plots and predictions from models ----

# Model list is just 1 and 3 for extraction of linear trend (a confidence criteria) and non-linear trend with non-linear random effect. 
Trop_Model_Upper_List <- list(Trop_Model_Upper1, Trop_Model_Upper3) # Create model list for functions to extract predictions
Trop_Model_Lower_List <- list(Trop_Model_Lower1, Trop_Model_Lower3) # Create model list for functions to extract predictions 

# Create list of objects for functions to plot outputs. 
# Extraction of global fixed effects. 
FixedFormula_Trop_List <- list(formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled), 
                               formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2)),
                               formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled  + I(MeanSiteSST_NOAA_Scaled^2)))

# Define random effects for input to function. 
RandomSlopes_Trop_List <- list(c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled'), 
                               c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled'), 
                               c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'))


# Create objects for plot data to be stored in. 
Trop_Model_Upper_Data <- list(); Trop_Model_Upper_RE <- list()
Trop_Model_Lower_Data <- list(); Trop_Model_Lower_RE <- list()

# Run functions for each model in loops. 
for(i in 1:length(unique(Trop_Model_Upper_List))){
  
  Trop_Model_Upper_Data[[i]]       <- ExtractGlobalTrend(Data = RLS_Trop_Occ_Tupper2, 
                                                         Model = Trop_Model_Upper_List[[i]], 
                                                         FixedFormula = FixedFormula_Trop_List[[i]], 
                                                         Error = 'binomial')
  Trop_Model_Upper_Data[[i]]$Model <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Trop_Model_Upper_List[[i]]))), 'AIC = ', round(AIC(Trop_Model_Upper_List[[i]])), i)[[3]]
  
  Trop_Model_Upper_RE[[i]]         <- CreateRandomEffectData(Data = RLS_Trop_Occ_Tupper2, 
                                                             Model = Trop_Model_Upper_List[[i]], 
                                                             FixedFormula = FixedFormula_Trop_List[[i]], 
                                                             RandomSlopes = RandomSlopes_Trop_List[[i]],
                                                             Error = 'binomial')
  Trop_Model_Upper_RE[[i]]$Model   <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Trop_Model_Upper_List[[i]]))), 'AIC = ', round(AIC(Trop_Model_Upper_List[[i]])), i)[[3]]
  
}
for(i in 1:length(unique(Trop_Model_Lower_List))){
  
  Trop_Model_Lower_Data[[i]]       <- ExtractGlobalTrend(Data = RLS_Trop_Occ_Tlower2, 
                                                         Model = Trop_Model_Lower_List[[i]], 
                                                         FixedFormula = FixedFormula_Trop_List[[i]], 
                                                         Error = 'binomial')
  Trop_Model_Lower_Data[[i]]$Model <-paste(as.character(gsub("(.{75})", "\\1 \n", formula(Trop_Model_Lower_List[[i]]))), 'AIC = ', round(AIC(Trop_Model_Lower_List[[i]])), i)[[3]]
  
  Trop_Model_Lower_RE[[i]]         <- CreateRandomEffectData(Data = RLS_Trop_Occ_Tlower2, 
                                                             Model = Trop_Model_Lower_List[[i]], 
                                                             FixedFormula = FixedFormula_Trop_List[[i]], 
                                                             RandomSlopes = RandomSlopes_Trop_List[[i]],
                                                             Error = 'binomial')
  Trop_Model_Lower_RE[[i]]$Model   <- paste(as.character(gsub("(.{75})", "\\1 \n", formula(Trop_Model_Lower_List[[i]]))), 'AIC = ', round(AIC(Trop_Model_Lower_List[[i]])), i)[[3]]
  
}

# Bind outputs from functions creating plot data together
Trop_Model_Upper_Data <- do.call(rbind, Trop_Model_Upper_Data)
Trop_Model_Upper_RE   <- do.call(rbind, Trop_Model_Upper_RE)
Trop_Model_Lower_Data <- do.call(rbind, Trop_Model_Lower_Data)
Trop_Model_Lower_RE   <- do.call(rbind, Trop_Model_Lower_RE)

# Plot uppers  
Trop_UpperPlot <- ggplot() + 
  #geom_line(data = Trop_Model_Upper_Data, aes(x = MeanSiteSST_NOAA_Scaled, y = mu)) + 
  geom_line(data = Trop_Model_Upper_RE, aes(x = MeanSiteSST_NOAA_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Trop_Model_Upper_Data, aes(label = Model, x = 0.5, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Plot lowers
Trop_LowerPlot <- ggplot() + 
  #geom_line(data = Trop_Model_Lower_Data, aes(x = MeanSiteSST_NOAA_Scaled, y = mu)) + 
  geom_line(data = Trop_Model_Lower_RE, aes(x = MeanSiteSST_NOAA_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Trop_Model_Lower_Data, aes(label = Model, x = -1, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Save in extra_figures
pdf(file = 'figures_extra/ThermalNicheLimits_TropicalOccModels.pdf')
Trop_LowerPlot
Trop_UpperPlot
dev.off()


# Model selection TROPICAL ----

# Remove one covariate
Trop_M1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                           (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))   # Remove reef
Trop_M2 <- formula(Presence ~ ReefAreaIn15km  + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans
Trop_M3 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove NPP
Trop_M4 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove temperature

# Remove two covariates
Trop_M5 <- formula(Presence ~ npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and humans 
Trop_M6 <- formula(Presence ~ HumPop50km2_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and npp 
Trop_M7 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and temperature
Trop_M8 <- formula(Presence ~ ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans and npp
Trop_M9 <- formula(Presence ~ ReefAreaIn15km  + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans and temperature
Trop_M10 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + 
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove npp and temperature

# Single covariates 
Trop_M11 <- formula(Presence ~ HumPop50km2_Scaled + 
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only humans
Trop_M12 <- formula(Presence ~ npp_mean_Scaled + 
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only NPP
Trop_M13 <- formula(Presence ~ ReefAreaIn15km + 
                      (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))        # Only Reef
Trop_M14 <- formula(Presence ~ MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2) + 
                     (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only temperature

# List of models to fit. 
Trop_Model_Selection <- list(Trop_M1, Trop_M2, Trop_M3, Trop_M4, Trop_M5, Trop_M6, Trop_M7, Trop_M8, Trop_M9, Trop_M10, Trop_M11, Trop_M12, Trop_M13, Trop_M14)

# Empty list to fill 
Trop_Models_Upper_MS <- list()
Trop_Models_Lower_MS <- list()

# For loop to fit all models. Takes a while > 12 hours to run.  
for(i in 1:length(Trop_Model_Selection)){
  print(paste('Model run', i))
  # Fit over each model. 
  Trop_Models_Upper_MS[[i]] <- glmmTMB(Trop_Model_Selection[[i]], family = "binomial", data = RLS_Trop_Occ_Tupper2)
  Trop_Models_Lower_MS[[i]] <- glmmTMB(Trop_Model_Selection[[i]], family = "binomial", data = RLS_Trop_Occ_Tlower2)
}

# Save model outputs
save(Trop_Models_Upper_MS, Trop_Models_Lower_MS, file = 'data_derived/Tropical-occupancy-model-selection.rdata')

# Last 4 models for model selection 
Trop_Model_Selection_V2 <- list(Trop_M11, Trop_M12, Trop_M13, Trop_M14)

# Empty list to fill 
Trop_Models_Upper_MS_V2 <- list()
Trop_Models_Lower_MS_V2 <- list()

# For loop to fit all models. Takes a while > 12 hours to run.  
for(i in 1:length(Trop_Model_Selection_V2)){
  print(paste('Model run', i))
  # Fit over each model. 
  Trop_Models_Upper_MS_V2[[i]] <- glmmTMB(Trop_Model_Selection_V2[[i]], family = "binomial", data = RLS_Trop_Occ_Tupper2)
  Trop_Models_Lower_MS_V2[[i]] <- glmmTMB(Trop_Model_Selection_V2[[i]], family = "binomial", data = RLS_Trop_Occ_Tlower2)
}

# Save verson 2. 
save(Trop_Models_Upper_MS_V2, Trop_Models_Lower_MS_V2, file = 'data_derived/Tropical-occupancy-model-selection_V2.rdata')

# Compare AIC from model lists for upper limits 
Trop_Upper_MS <- c(list(Trop_Model_Upper3), Trop_Models_Upper_MS, Trop_Models_Upper_MS_V2) # Compare with full model that has 2nd order polynomial too. 
bbmle::AICtab(Trop_Upper_MS)
Trop_Upper_MS[[3]] %>% fixef()
Trop_Upper_MS[[12]] %>% fixef()
Trop_Upper_MS[[1]] %>% fixef() # Here i'm going to choose model 1 as there is little difference in AIC. 
bbmle::AICtab(Trop_Upper_MS[[1]], Trop_Upper_MS[[3]])
FinalTropUpperModel <- Trop_Upper_MS[[3]]

# Compare AIC from model lists for lower limits 
Trop_Lower_MS <- c(list(Trop_Model_Lower3), Trop_Models_Lower_MS, Trop_Models_Lower_MS_V2) # Compare with full model that has 2nd order polynomial too. 
bbmle::AICtab(Trop_Lower_MS)
Trop_Lower_MS[[1]]  %>% fixef()

# Extract final model
FinalTropLowerModel <- Trop_Lower_MS[[1]]

save(FinalTropUpperModel, FinalTropLowerModel, file = 'data_derived/Tropical-occupancy-final-models.rdata')

# load in model selections and final models from above code ----
# Load these to create the models to select between. But then remove as will clog up memory. 
load(file = 'data_derived/Tropical-occupancy-model-selection.rdata')
load(file = 'data_derived/Tropical-occupancy-model-selection_V2.rdata')
load(file = 'data_derived/Tropical-occupancy-final-models.rdata')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# OCCUPANCY QUALITY CONTROLS ----
# Fit models with same covariates as AIC selected formula apart from fitting with linear slopes to act as a quality controls ----

# Temperate upper
formula(FinalTempUpperModel)
TempUpperModel_Linear <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + (1 + MeanSiteSST_NOAA_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTempUpperModel_Linear <- glmmTMB(TempUpperModel_Linear, family = "binomial", data = RLS_Temp_Occ_Tupper2)

# Temperate lower
formula(FinalTempLowerModel)
TempLowerModel_Linear <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + (1 + MeanSiteSST_NOAA_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTempLowerModel_Linear <- glmmTMB(TempLowerModel_Linear, family = "binomial", data = RLS_Temp_Occ_Tlower2)

# Tropical upper
formula(FinalTropUpperModel)
TropUpperModel_Linear <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTropUpperModel_Linear <- glmmTMB(TropUpperModel_Linear, family = "binomial", data = RLS_Trop_Occ_Tupper2)

# Tropical lower
# Extract final model quality control linear slope. 
formula(FinalTropLowerModel)
TropLowerModel_Linear <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + (1 + ReefAreaIn15km + MeanSiteSST_NOAA_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTropLowerModel_Linear <- glmmTMB(TropLowerModel_Linear, family = "binomial", data = RLS_Trop_Occ_Tlower2)

# Save these
save(FinalTempUpperModel_Linear, FinalTempLowerModel_Linear, FinalTropUpperModel_Linear, FinalTropLowerModel_Linear, file = 'data_derived/occupancy-models-linear-slopes.rdata')


# load in data from code above ----
load(file = 'data_derived/occupancy-models-linear-slopes.rdata')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MODEL SUMMARIES ----
# Summaries of occupancy models for supporting table S1 and table S2 ----
# Tropical models
summary(Trop_Model_Lower1)
summary(Trop_Model_Lower3)
summary(Trop_Model_Upper1)
summary(Trop_Model_Upper3)

# Temperate models
summary(Temp_Model_Lower1)
summary(Temp_Model_Lower3)
summary(Temp_Model_Upper1)
summary(Temp_Model_Upper3)

# Summaries of fixed effects model selections
summary(FinalTempUpperModel)
summary(FinalTempLowerModel)
summary(FinalTropUpperModel)
summary(FinalTropLowerModel)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# EXTRACT NICHE LIMITS FROM RANDOM EFFECTS AND PLOT IN MULTIPANEL ----
# Function to extract niche limits from models ----
ExtractPredictions <- function(Data, Model, ModelNames, RandomSlopes, FixedFormula, Limit, 
                               i = i # This delimits the species to subset by and the random effect
){ 
  
  Data_all <- Data
  
  if(Limit == 'Upper'){
    
    Data <- Data %>% filter(SpeciesName == unique(.$SpeciesName)[i])
    Data$T_Mean_Obs_Scaled <- (Data$T_Mean_Obs - Data$MeanSiteSST_NOAA_centerBT) / Data$MeanSiteSST_NOAA_scaleBT
    Data$T_Upper_Absences_Scaled <- (Data$T_Upper_Absences+5 - Data$MeanSiteSST_NOAA_centerBT) / Data$MeanSiteSST_NOAA_scaleBT
    
    # Extract prediction frame from the data. 
    Preds <- Data %>% 
      select(SpeciesName, ThermalGuild, MeanSiteSST_NOAA_Scaled, T_Mean_Obs_Scaled, T_Upper_Absences_Scaled) %>% 
      group_by(SpeciesName, ThermalGuild) %>% 
      nest() %>% 
      mutate(MeanSiteSST_NOAA_Scaled = purrr::map(data, ~seq(unique(.$T_Mean_Obs_Scaled), unique(.$T_Upper_Absences_Scaled), length.out = 1000))) %>% 
      unnest(MeanSiteSST_NOAA_Scaled)
    
  }else{
    
    Data <- Data %>% filter(SpeciesName == unique(.$SpeciesName)[i])
    Data$T_Mean_Obs_Scaled <- (Data$T_Mean_Obs - Data$MeanSiteSST_NOAA_centerBT) / Data$MeanSiteSST_NOAA_scaleBT
    Data$T_Lower_Absences_Scaled <- (Data$T_Lower_Absences-5 - Data$MeanSiteSST_NOAA_centerBT) / Data$MeanSiteSST_NOAA_scaleBT
    
    # Extract prediction frame from the data. 
    Preds <- Data %>% 
      select(SpeciesName, ThermalGuild, MeanSiteSST_NOAA_Scaled, T_Mean_Obs_Scaled, T_Lower_Absences_Scaled) %>% 
      group_by(SpeciesName, ThermalGuild) %>% 
      nest() %>% 
      mutate(MeanSiteSST_NOAA_Scaled = purrr::map(data, ~seq(unique(.$T_Lower_Absences_Scaled), unique(.$T_Mean_Obs_Scaled), length.out = 1000))) %>% 
      unnest(MeanSiteSST_NOAA_Scaled)
    
  }
  
  Preds$HumPop50km2_Scaled <- mean(Data_all$HumPop50km2_Scaled)
  Preds$npp_mean_Scaled <- mean(Data_all$npp_mean_Scaled)
  if(unique(Data$ThermalGuild) == 'tropical'){Preds$ReefAreaIn15km <- mean(Data_all$ReefAreaIn15km)}else{Preds$ReefAreaIn15km <- NA}
  
  # Extract species parameter estimates from random effects for species of interest
  GlobalFixedEffects <- fixef(Model)$con
  SpeciesParameters <- fixef(Model)$con 
  
  # Check if the random slope appears in the fixed component. If not then assign random slope post-extraction. 
  if(identical(names(fixef(Model)$con[c(RandomSlopes)]) , RandomSlopes)) { 
    SpeciesParameters[c('(Intercept)', RandomSlopes)] <- (fixef(Model)$con[c('(Intercept)', RandomSlopes)] + ranef(Model)$con$SpeciesName[i,c('(Intercept)', RandomSlopes)])
    SpeciesParameters <- data.frame(SpeciesParameters)
    colnames(SpeciesParameters) <- names(fixef(Model)$con)
  }else{
    
    #RandomSlope_RawExtract <- RandomSlopes[which(!names(fixef(Model)$con[c('(Intercept)', RandomSlopes)]) %in% RandomSlopes)[-1]-1]
    #RandomSlope_NumExtract <- which(!names(fixef(Model)$con[c('(Intercept)', RandomSlopes)]) %in% RandomSlopes)[-1]
    #SpeciesParameters[c('(Intercept)', RandomSlopes)][-RandomSlope_NumExtract] <- fixef(Model)$con[c('(Intercept)', RandomSlopes)][-RandomSlope_NumExtract] + ranef(Model)$con$SpeciesName[i,c('(Intercept)', RandomSlopes)][-RandomSlope_NumExtract]
    #SpeciesParameters <- data.frame(SpeciesParameters)
    #SpeciesParameters[which(is.na(SpeciesParameters))] <- ranef(Model)$con$SpeciesName[i,c('(Intercept)', RandomSlopes)][RandomSlope_NumExtract]
    
  }
  
  # Create dataframe
  PredEstimates <- data.frame(Preds %>% filter(SpeciesName == unique(Data$SpeciesName)))
  
  # Create model matrix for each species
  X <- model.matrix(FixedFormula, data = PredEstimates)
  PredEstimates$eta <- as.numeric(X %*% as.numeric(SpeciesParameters))
  PredEstimates$mu <- as.numeric(exp(X %*% as.numeric(SpeciesParameters)) / (1 + exp(X %*% as.numeric(SpeciesParameters))))
  PredEstimates$mu_scaled <- as.numeric(PredEstimates$mu / max(PredEstimates$mu, na.rm = T))
  
  # Handle outputs 
  PredEstimates$MeanSiteSST_NOAA <- (PredEstimates$MeanSiteSST_NOAA_Scaled *  Data$MeanSiteSST_NOAA_scaleBT[1]) + Data$MeanSiteSST_NOAA_centerBT[1]
  
  PredEstimates$Model <- ModelNames
  
  if(Limit == 'Upper'){
    
    ### Estimate Tupper now from this predicted relationship. 
    #Where species have too low confidence above and below sampling limits they are given an NA.
    #min_mu_scaled <- min(PredEstimates$mu_scaled)
    #above_min_mu <- PredEstimates[which(PredEstimates$mu_scaled > min_mu_scaled),]
    #below_min_mu <- PredEstimates[which(PredEstimates$mu_scaled < min_mu_scaled),]
    
    #if(nrow(above_min_mu) > 0 | nrow(below_min_mu) > 0){
    PredEstimates$T_Upper <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$MeanSiteSST_NOAA # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
    PredEstimates$T_Upper_mu_scaled <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$mu_scaled
    # T_SD_Upper <- abs(T_Upper - Data$T_Opt) / 2 # This is a better way of estimateing thermal niche. Assuming that 95% of data fall in 2SD of mean. So Tlower is 0.025. Tupper is 0.975
    #}else{
    # This is where things get tricky. 
    
  }else{NULL}
  if(Limit == 'Lower'){
    ### Estimate TLower now from this predicted relationship. 
    # Where species have too low confidence above and below sampling limits they are given an NA.
    PredEstimates$T_Lower <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$MeanSiteSST_NOAA # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
    PredEstimates$T_Lower_mu_scaled <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$mu_scaled
    # T_SD_Lower <- abs(T_Lower - Data$T_Opt) / 2 # This is a better way of estimateing thermal niche. Assuming that 95% of data fall in 2SD of mean. So Tlower is 0.025. TLower is 0.975
  }else{NULL}
  if(sum(grepl(Limit, c('Upper', 'Lower'))) == 0){return(warning('Error: non-specified upper or lower limit'))}else{NULL}
  
  return(PredEstimates)
  
}

# Extract upper niches estimates from the above occupancy modelling and confidence criteria ----

# Extract predictions: Temperate upper
# Extract from model with linear slope (useful only as confidence criteria)
TempUpper_Output_Linear <- list()
for(i in 1:length(unique(RLS_Temp_Occ_Tupper2$SpeciesName))){
  TempUpper_Output_Linear[[i]] <- ExtractPredictions(Data         = RLS_Temp_Occ_Tupper2,
                                                            Model        = FinalTempUpperModel_Linear,
                                                            ModelNames   = 'Temp_Upper_Occ_Linear',
                                                            RandomSlopes = 'MeanSiteSST_NOAA_Scaled',
                                                            FixedFormula = formula(~HumPop50km2_Scaled + MeanSiteSST_NOAA_Scaled), 
                                                     Limit = 'Upper', 
                                                     i = i)
  print(i / length(unique(RLS_Temp_Occ_Tupper2$SpeciesName)) * 100 )
}

# Extract from model with non-linear quadratic slope + random effect (final model used in analysis )
TempUpper_Output <- list()
for(i in 1:length(unique(RLS_Temp_Occ_Tupper2$SpeciesName))){
  TempUpper_Output[[i]] <- ExtractPredictions(       Data         = RLS_Temp_Occ_Tupper2,
                                                     Model        = FinalTempUpperModel,
                                                     ModelNames   = 'Temp_Upper_Occ',
                                                     RandomSlopes = c('MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'),
                                                     FixedFormula = formula(~HumPop50km2_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)), 
                                                     Limit = 'Upper', 
                                                     i = i)
  print(i / length(unique(RLS_Temp_Occ_Tupper2$SpeciesName)) * 100 )
}

# Extract predictions: Tropical upper
# Extract from model with linear slope (useful only as confidence criteria)
TropUpper_Output_Linear <- list()
TropUpper_Output        <- list()
for(i in 1:length(unique(RLS_Trop_Occ_Tupper2$SpeciesName))){
  TropUpper_Output_Linear[[i]] <- ExtractPredictions(Data         = RLS_Trop_Occ_Tupper2,
                                                            Model        = FinalTropUpperModel_Linear,
                                                            ModelNames   = 'Trop_Upper_Occ_Linear',
                                                            RandomSlopes = c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled), 
                                                     Limit = 'Upper',
                                                     i = i)
  
  TropUpper_Output[[i]] <- ExtractPredictions(       Data         = RLS_Trop_Occ_Tupper2,
                                                            Model        = FinalTropUpperModel,
                                                            ModelNames   = 'Trop_Upper_Occ',
                                                            RandomSlopes = c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)), 
                                                     Limit = 'Upper', 
                                                     i = i)
  
  
  print(i / length(unique(RLS_Trop_Occ_Tupper2$SpeciesName)) * 100 )
}


# Extract non-negative linear slopes. 
# These are used as the bases of confidence in this parameter. 
Ranef_linear_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropUpperModel_Linear)$cond$SpeciesName), 
                                LinearSlope = ranef(FinalTropUpperModel_Linear)$cond$SpeciesName$MeanSiteSST_NOAA_Scaled + fixef(FinalTropUpperModel_Linear)$cond['MeanSiteSST_NOAA_Scaled'])

Ranef_linear_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempUpperModel_Linear)$cond$SpeciesName), 
                                LinearSlope  = ranef(FinalTempUpperModel_Linear)$cond$SpeciesName$MeanSiteSST_NOAA_Scaled + fixef(FinalTempUpperModel_Linear)$cond['MeanSiteSST_NOAA_Scaled'])

Ranef_linear <- rbind(Ranef_linear_Trop, Ranef_linear_Temp)

Ranef_linear$ConfidenceLinearSlope_Upper <- ifelse(Ranef_linear$LinearSlope < -1, 1, 0)

ConfidenceLinearSlope_Upper_Species <- Ranef_linear$SpeciesName[which(Ranef_linear$ConfidenceLinearSlope_Upper == 1)]

# Extract quadratic slopes. 
# These are used to support the confidence in negative linear slopes. Values between -1 and 0 can be weakly negative and overestimate thermal niche.
# We use a cut off of 'unimodal' curves as less than -1. 
# We use a cut-off of 'u-shaped' curves as greater than 0. 
Ranef_quad_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropUpperModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTropUpperModel)$cond$SpeciesName$`I(MeanSiteSST_NOAA_Scaled^2)` + fixef(FinalTropUpperModel)$cond['I(MeanSiteSST_NOAA_Scaled^2)'])

Ranef_quad_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempUpperModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTempUpperModel)$cond$SpeciesName$`I(MeanSiteSST_NOAA_Scaled^2)` + fixef(FinalTempUpperModel)$cond['I(MeanSiteSST_NOAA_Scaled^2)'])

Ranef_quad <- rbind(Ranef_quad_Trop, Ranef_quad_Temp)

# Thresholds for quanratic slopes
Ranef_quad$ConfidenceQuadSlope_Upper_Unimodal <- ifelse(Ranef_quad$QuadSlope < -1, 1, 0)
Ranef_quad$ConfidenceQuadSlope_Upper_U_shaped <- ifelse(Ranef_quad$QuadSlope > 0, 1, 0)

# Assigns high confidence species in quadratic term to objects. 
ConfidenceQuadSlope_Upper_Species_Unimodal <- Ranef_quad$SpeciesName[which(Ranef_quad$ConfidenceQuadSlope_Upper_Unimodal == 1)]
ConfidenceQuadSlope_Upper_Species_UShaped <- Ranef_quad$SpeciesName[which(Ranef_quad$ConfidenceQuadSlope_Upper_U_shaped == 1)]

# Selections high confidence and linear from also high confidence and quadratic.  
Confidence_OccUpperAll_Unimodal <- unique(ConfidenceLinearSlope_Upper_Species[ConfidenceLinearSlope_Upper_Species %in% ConfidenceQuadSlope_Upper_Species_Unimodal])
Confidence_OccUpperAll_U_Shaped <- unique(ConfidenceLinearSlope_Upper_Species[ConfidenceLinearSlope_Upper_Species %in% ConfidenceQuadSlope_Upper_Species_UShaped])

# Bind fits for plotting
Fits_UpperLimits_Linear    <- as_data_frame(rbind(do.call(rbind, TempUpper_Output_Linear), do.call(rbind, TropUpper_Output_Linear)))
Fits_UpperLimits           <- as_data_frame(rbind(do.call(rbind, TempUpper_Output),        do.call(rbind, TropUpper_Output)))

# Bind data for plotting 
AllData_UpperLimits <- rbind(RLS_Temp_Occ_Tupper2, RLS_Trop_Occ_Tupper2)

#### Plot all the thermal upper limits and visually inspect these arbitrary thresholds. 
pdf(file = 'figures_extra/Niche_upper_models_2018-03-05.pdf', width = 10, height = 10)
par(mfrow = c(5,5), mar = c(2,4,2,2))
# Plot those with unimodal curves (i.e., negative coefficients.)
for(i in 1:length(unique(Confidence_OccUpperAll_Unimodal))){
  print(i)
  
  # Model lines extracted from fitted models
  SpeciesFits_Linear <- Fits_UpperLimits_Linear %>% filter(SpeciesName == Confidence_OccUpperAll_Unimodal[i])
  SpeciesFits        <- Fits_UpperLimits %>% filter(SpeciesName == Confidence_OccUpperAll_Unimodal[i])
  
  # Raw data input to fitted model. 
  Species_AllData <- RLS_All %>% filter(SpeciesName == Confidence_OccUpperAll_Unimodal[i])
  
  # Create plots
  plot(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, type = 'l', main = as.character(Confidence_OccUpperAll_Unimodal)[i], ylim = c(0,1), xlim = c(min(SpeciesFits$MeanSiteSST_NOAA),unique(SpeciesFits$T_Upper)+1))
  points(Presence ~ MeanSiteSST_NOAA, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Upper), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits, col = 'dark orange')
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, col = 'Black')
}
# Plot those with u-shaped inverse parabolas (i.e., positive coefficients)
for(i in 1:length(unique(Confidence_OccUpperAll_U_Shaped))){
  print(i)
  
  # Model lines extracted from fitted models
  SpeciesFits_Linear <- Fits_UpperLimits_Linear %>% filter(SpeciesName == Confidence_OccUpperAll_U_Shaped[i])
  SpeciesFits        <- Fits_UpperLimits %>% filter(SpeciesName == Confidence_OccUpperAll_U_Shaped[i])
  
  # Raw data input to fitted model. 
  Species_AllData <- RLS_All %>% filter(SpeciesName == Confidence_OccUpperAll_U_Shaped[i])
  
  # Create plots
  plot(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, type = 'l', main = unique(SpeciesFits$SpeciesName), ylim = c(0,1), xlim = c(min(SpeciesFits$MeanSiteSST_NOAA),unique(SpeciesFits$T_Upper)+1))
  points(Presence ~ MeanSiteSST_NOAA, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Upper), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits, col = 'dark orange', lty = 2)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, col = 'Black')
}
dev.off()
####

# Manually define confidence upper all group from some poor visual fits when the coefficient for temperature quadratic term is > 0. 
# This step is necessary because when the parabola is inverted curves y begin and end at 1, and do not have to reach = 0, meaning that thermal limits may not necessarrily be well defined. 
# See figure with '_ManualQualityControl' for highlighted species.   
VisualPoorFits <- c('Abudefduf sexfasciatus',        # Overestimation of thermal limit to to u-shaped curve
                    'Abudefduf whitleyi',            # Overestimation of thermal limit to to u-shaped curve
                    'Acanthurus auranticavus',       # Very flat curve within range placing thermal limit within temperature range. 
                    'Amblygobius decussatus',        # Very flat curve within range placing thermal limit within temperature range. 
                    'Anampses caeruleopunctatus',    # Very flat curve within range placing thermal limit within temperature range. 
                    'Assessor macneilli',            # Overestimation of thermal limit to to u-shaped curve
                    'Cypho purpurascens',            # Very flat curve within range placing thermal limit within temperature range. 
                    'Chaetodontoplus duboulayi',      # Overestimation of thermal limit to to u-shaped curve
                    'Chaetodontoplus meredithi',      # Very flat curve within range placing thermal limit within temperature range.
                    'Chromis agilis',                # Overestimation of thermal limit to to u-shaped curve
                    'Chromis nitida',                # Very flat curve within range placing thermal limit within temperature range.
                    'Crossosalarias macrospilus',    # Overestimation of thermal limit to to u-shaped curve
                    'Ecsenius mandibularis',         # Very flat curve within range placing thermal limit within temperature range. 
                    'Gerres subfasciatus',          # Overestimation of thermal limit to to u-shaped curve
                    'Gymnothorax javanicus',         # Overestimation of thermal limit to to u-shaped curve
                    'Labropsis australis',           # Overestimation of thermal limit to to u-shaped curve
                    'Macropharyngodon choati',       # Overestimation of thermal limit to to u-shaped curve
                    'Monotaxis heterodon',           # Overestimation of thermal limit to to u-shaped curve
                    'Myripristis kuntee',            # Very flat curve within range placing thermal limit within temperature range. 
                    'Parupeneus barberinoides',      # Very flat curve within range placing thermal limit within temperature range. 
                    'Pomacentrus grammorhynchus',    # Overestimation of thermal limit to to u-shaped curve
                    'Pterocaesio digramma',          # Very flat curve within range placing thermal limit within temperature range. 
                    'Scarus schlegeli',              # Very flat curve within range placing thermal limit within temperature range. 
                    'Sphyraena obtusata',             # Very flat curve within range placing thermal limit within temperature range. 
                    'Stegastes beebei',              # Overestimation of thermal limit to to u-shaped curve 
                    'Thalassoma pavo',               # Overestimation of thermal limit to to u-shaped curve 
                    'Thalassoma purpureum',           # Very flat curve within range placing thermal limit within temperature range. 
                    'Chrysiptera notialis')          # Overestimation of thermal limit to to u-shaped curve 

# Which species have mu > 0.05. For these species is it difficult to define a minimum as the occupancy does not decline monotonically. 
Confidence_mu_scaled_Upper <- Fits_UpperLimits %>% select(SpeciesName, T_Upper_mu_scaled) %>% unique() %>% filter(T_Upper_mu_scaled >= 0.05) %>% .$SpeciesName %>% as.character()

# Combine all confidence estimates
Confidence_OccUpperAll <- c(as.character(Confidence_OccUpperAll_Unimodal), as.character(Confidence_OccUpperAll_U_Shaped))
Confidence_OccUpperAll <- Confidence_OccUpperAll[-which(Confidence_OccUpperAll %in% c(VisualPoorFits, Confidence_mu_scaled_Upper))] 

# Join together sampled niches and observed niches. 
Species_Tuppers <- Fits_UpperLimits %>% 
  filter(SpeciesName %in% as.character(Confidence_OccUpperAll)) %>%
  select(SpeciesName, ThermalGuild, T_Upper) %>% unique()
Species_Tuppers2 <- left_join(Species_Tuppers, ThermalNicheData_Obs)

# Join together niches with parameter values to see how these affect estimates. 
Species_Tuppers3 <- left_join(Species_Tuppers2, left_join(Ranef_quad, Ranef_linear))

# On average there is a xxx degree difference in modelled temperature upper limit and sampled upper limit. 
Temperate_Difference_Upper <- Species_Tuppers3$T_Upper[which(Species_Tuppers3$ThermalGuild == 'temperate')] - Species_Tuppers3$T_Upper_Obs[which(Species_Tuppers3$ThermalGuild == 'temperate')] 
Tropical_Difference_Upper  <- Species_Tuppers3$T_Upper[which(Species_Tuppers3$ThermalGuild == 'tropical')] - Species_Tuppers3$T_Upper_Obs[which(Species_Tuppers3$ThermalGuild == 'tropical')] 

# Median difference used to adjust thermal performance in external analysis. 
median(Temperate_Difference_Upper) # 1.37472
median(Tropical_Difference_Upper)  # 1.43299



# Extract lower niches estimates from the above occupancy modelling and confidence criteria --------------
# Extract predictions: Temperate lower
# Extract from model with linear slope (useful only as confidence criteria)
TempLower_Output_Linear <- list()
for(i in 1:length(unique(RLS_Temp_Occ_Tlower2$SpeciesName))){
  TempLower_Output_Linear[[i]] <- ExtractPredictions(Data         = RLS_Temp_Occ_Tlower2,
                                                            Model        = FinalTempLowerModel_Linear,
                                                            ModelNames   = 'Temp_Lower_Occ_Linear',
                                                            RandomSlopes = 'MeanSiteSST_NOAA_Scaled',
                                                            FixedFormula = formula(~HumPop50km2_Scaled + MeanSiteSST_NOAA_Scaled), 
                                                     Limit = 'Lower', 
                                                     i = i)
  print(i / length(unique(RLS_Temp_Occ_Tlower2$SpeciesName)) * 100 )
}

# Extract from model with non-linear quadratic slope + random effect (final model used in analysis )
TempLower_Output <- list()
for(i in 1:length(unique(RLS_Temp_Occ_Tlower2$SpeciesName))){
  TempLower_Output[[i]] <- ExtractPredictions(       Data         = RLS_Temp_Occ_Tlower2,
                                                            Model        = FinalTempLowerModel,
                                                            ModelNames   = 'Temp_Lower_Occ',
                                                            RandomSlopes = c('MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)), 
                                                            Limit = 'Lower', 
                                                     i = i)
  print(i / length(unique(RLS_Temp_Occ_Tlower2$SpeciesName)) * 100 )
}

# Extract predictions: Tropical lower
# Extract from model with linear slope (useful only as confidence criteria)
TropLower_Output_Linear <- list()
TropLower_Output        <- list()
for(i in 1:length(unique(RLS_Trop_Occ_Tlower2$SpeciesName))){
  TropLower_Output_Linear[[i]] <- ExtractPredictions(Data         = RLS_Trop_Occ_Tlower2,
                                                     Model        = FinalTropLowerModel_Linear,
                                                     ModelNames   = 'Trop_Lower_Occ_Linear',
                                                     RandomSlopes = c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled'),
                                                     FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled), 
                                                     Limit = 'Lower', 
                                                     i = i)
  
  TropLower_Output[[i]] <- ExtractPredictions(       Data         = RLS_Trop_Occ_Tlower2,
                                                     Model        = FinalTropLowerModel,
                                                     ModelNames   = 'Trop_Lower_Occ',
                                                     RandomSlopes = c('ReefAreaIn15km', 'MeanSiteSST_NOAA_Scaled', 'I(MeanSiteSST_NOAA_Scaled^2)'),
                                                     FixedFormula = formula(~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanSiteSST_NOAA_Scaled + I(MeanSiteSST_NOAA_Scaled^2)), 
                                                     Limit = 'Lower', 
                                                     i = i)
  
  print(i / length(unique(RLS_Trop_Occ_Tlower2$SpeciesName)) * 100 )
}

# Extract non-negative linear slopes. 
# These are used as the bases of confidence in this parameter. 
Ranef_linear_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropLowerModel_Linear)$cond$SpeciesName), 
                                LinearSlope = ranef(FinalTropLowerModel_Linear)$cond$SpeciesName$MeanSiteSST_NOAA_Scaled + fixef(FinalTropLowerModel_Linear)$cond['MeanSiteSST_NOAA_Scaled'])

Ranef_linear_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempLowerModel_Linear)$cond$SpeciesName), 
                                LinearSlope  = ranef(FinalTempLowerModel_Linear)$cond$SpeciesName$MeanSiteSST_NOAA_Scaled + fixef(FinalTempLowerModel_Linear)$cond['MeanSiteSST_NOAA_Scaled'])

Ranef_linear <- rbind(Ranef_linear_Trop, Ranef_linear_Temp)

Ranef_linear$ConfidenceLinearSlope_Lower <- ifelse(Ranef_linear$LinearSlope > 1, 1, 0)

ConfidenceLinearSlope_Lower_Species <- Ranef_linear$SpeciesName[which(Ranef_linear$ConfidenceLinearSlope_Lower == 1)]

# Extract quadratic slopes. 
# These are used to support the confidence in negative linear slopes. Values between -1 and 0 can be weakly negative and overestimate thermal niche.
# We use a cut off of 'unimodal' curves as less than -1. 
# We use a cut-off of 'u-shaped' curves as greater than 0. 
Ranef_quad_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropLowerModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTropLowerModel)$cond$SpeciesName$`I(MeanSiteSST_NOAA_Scaled^2)` + fixef(FinalTropLowerModel)$cond['I(MeanSiteSST_NOAA_Scaled^2)'])

Ranef_quad_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempLowerModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTempLowerModel)$cond$SpeciesName$`I(MeanSiteSST_NOAA_Scaled^2)` + fixef(FinalTempLowerModel)$cond['I(MeanSiteSST_NOAA_Scaled^2)'])

Ranef_quad <- rbind(Ranef_quad_Trop, Ranef_quad_Temp)

# Thresholds for quanratic slopes
Ranef_quad$ConfidenceQuadSlope_Lower_Unimodal <- ifelse(Ranef_quad$QuadSlope < -1, 1, 0)
Ranef_quad$ConfidenceQuadSlope_Lower_U_shaped <- ifelse(Ranef_quad$QuadSlope > 0, 1, 0)

# Assigns high confidence species in quadratic term to objects. 
ConfidenceQuadSlope_Lower_Species_Unimodal <- Ranef_quad$SpeciesName[which(Ranef_quad$ConfidenceQuadSlope_Lower_Unimodal == 1)]
ConfidenceQuadSlope_Lower_Species_UShaped <- Ranef_quad$SpeciesName[which(Ranef_quad$ConfidenceQuadSlope_Lower_U_shaped == 1)]

# Selections high confidence and linear from also high confidence and quadratic.  
Confidence_OccLowerAll_Unimodal <- unique(ConfidenceLinearSlope_Lower_Species[ConfidenceLinearSlope_Lower_Species %in% ConfidenceQuadSlope_Lower_Species_Unimodal])
Confidence_OccLowerAll_U_Shaped <- unique(ConfidenceLinearSlope_Lower_Species[ConfidenceLinearSlope_Lower_Species %in% ConfidenceQuadSlope_Lower_Species_UShaped])

# Bind fits for plotting
Fits_LowerLimits_Linear    <- as_data_frame(rbind(do.call(rbind, TempLower_Output_Linear), do.call(rbind, TropLower_Output_Linear)))
Fits_LowerLimits           <- as_data_frame(rbind(do.call(rbind, TempLower_Output),        do.call(rbind, TropLower_Output)))

# Bind data for plotting 
AllData_LowerLimits <- rbind(RLS_Temp_Occ_Tlower2, RLS_Trop_Occ_Tlower2)


#### Plot all the thermal lower limits and visually inspect these arbitrary thresholds. 
pdf(file = 'figures_extra/Niche_lower_models_2018-03-05.pdf', width = 10, height = 10)
par(mfrow = c(5,5), mar = c(2,4,2,2))
# Plot those with unimodal curves (i.e., negative coefficients.)
for(i in 1:length(unique(Confidence_OccLowerAll_Unimodal))){
  print(i)
  
  # Model lines extracted from fitted models
  SpeciesFits_Linear <- Fits_LowerLimits_Linear %>% filter(SpeciesName == Confidence_OccLowerAll_Unimodal[i])
  SpeciesFits        <- Fits_LowerLimits %>% filter(SpeciesName == Confidence_OccLowerAll_Unimodal[i])
  
  # Raw data input to fitted model. 
  Species_AllData <- RLS_All %>% filter(SpeciesName == Confidence_OccLowerAll_Unimodal[i])
  
  # Create plots
  plot(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, type = 'l', main = Confidence_OccLowerAll_Unimodal[i], ylim = c(0,1), xlim = c(unique(SpeciesFits$T_Lower)-1, max(SpeciesFits$MeanSiteSST_NOAA)))
  points(Presence ~ MeanSiteSST_NOAA, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Lower), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits, col = 'dark orange')
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, col = 'Black')
}
# Plot those with u-shaped inverse parabolas (i.e., positive coefficients)
for(i in 1:length(unique(Confidence_OccLowerAll_U_Shaped))){
  print(i)
  
  # Model lines extracted from fitted models
  SpeciesFits_Linear <- Fits_LowerLimits_Linear %>% filter(SpeciesName == Confidence_OccLowerAll_U_Shaped[i])
  SpeciesFits        <- Fits_LowerLimits %>% filter(SpeciesName == Confidence_OccLowerAll_U_Shaped[i])
  
  # Raw data input to fitted model. 
  Species_AllData <- RLS_All %>% filter(SpeciesName == Confidence_OccLowerAll_U_Shaped[i])
  
  # Create plots
  plot(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, type = 'l', main = Confidence_OccLowerAll_U_Shaped[i], ylim = c(0,1), xlim = c(unique(SpeciesFits$T_Lower)-1, max(SpeciesFits$MeanSiteSST_NOAA)))
  points(Presence ~ MeanSiteSST_NOAA, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Lower), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits, col = 'dark orange', lty = 2)
  lines(mu_scaled ~ MeanSiteSST_NOAA, data = SpeciesFits_Linear, col = 'Black')
}
dev.off()
####

# Manually define confidence Lower all group from some poor visual fits when the coefficient for temperature quadratic term is > 0. 
# This step is necessary because when the parabola is inverted curves y begin and end at 1, and do not have to reach = 0, meaning that thermal limits may not necessarrily be well defined. 
# See figure with '_ManualQualityControl' for highlighted species.   
VisualPoorFits_lower <- c('Naso vlamingii',         # Underestimation of thermal limit to to u-shaped curve 
                          'Apogon flavus',                # Underestimation of thermal limit to to u-shaped curve 
                          'Chrysiptera notialis',         # Underestimation of thermal limit to to u-shaped curve
                          'Neoglyphidodon polyacanthus',  # Underestimation of thermal limit to to u-shaped curve
                          'Trachypoma macracanthus'       # Underestimation of thermal limit to to u-shaped curve
)


# Which species have mu > 0.05. For these species is it difficult to define a minimum as the occupancy does not decline monotonically. 
Confidence_mu_scaled_Lower <- Fits_LowerLimits %>% select(SpeciesName, T_Lower_mu_scaled) %>% unique() %>% filter(T_Lower_mu_scaled >= 0.05) %>% .$SpeciesName %>% as.character()

# Combine all confidence estimates
Confidence_OccLowerAll <- c(as.character(Confidence_OccLowerAll_Unimodal), as.character(Confidence_OccLowerAll_U_Shaped))
Confidence_OccLowerAll <- Confidence_OccLowerAll[-which(Confidence_OccLowerAll %in% c(VisualPoorFits_lower, Confidence_mu_scaled_Lower))] 

# Join together sampled niches and observed niches. 
Species_TLowers <- Fits_LowerLimits %>% 
  filter(SpeciesName %in% as.character(Confidence_OccLowerAll)) %>%
  select(SpeciesName, ThermalGuild, T_Lower) %>% unique()
Species_TLowers2 <- left_join(Species_TLowers, ThermalNicheData_Obs)

# Join together niches with parameter values to see how these affect estimates. 
Species_TLowers3 <- left_join(Species_TLowers2, left_join(Ranef_quad, Ranef_linear))

# On average there is a xxx degree difference in modelled temperature Lower limit and sampled Lower limit. 
Temperate_Difference_Lower <- Species_TLowers3$T_Lower[which(Species_TLowers3$ThermalGuild == 'temperate')] - Species_TLowers3$T_Lower_Obs[which(Species_TLowers3$ThermalGuild == 'temperate')] 
Tropical_Difference_Lower  <- Species_TLowers3$T_Lower[which(Species_TLowers3$ThermalGuild == 'tropical')]  -  Species_TLowers3$T_Lower_Obs[which(Species_TLowers3$ThermalGuild == 'tropical')] 

# Median difference used to adjust thermal performance in external analysis. 
median(Temperate_Difference_Lower) # -1.14139
median(Tropical_Difference_Lower)  # -1.04672



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CREATE DATASET OF THERMAL PERFORMANCE PARAMETERS ---- 
# Combine all thermal niche parameters ----

# Thermal maximums from occupancy models. 
ThermalMax <- Species_Tuppers %>% select(SpeciesName, T_Upper)
ThermalMax$T_Upper_Mod <- ThermalMax$T_Upper; ThermalMax$T_Upper <- NULL
ThermalNicheData_New <- left_join(ThermalNicheData_New, ThermalMax)

# Thermal maximums adjusted from occupancy model outputs. 
ThermalNicheData_New <- left_join(ThermalNicheData_New, ThermalNicheData_Obs)
ThermalNicheData_New$T_Upper <- NA
ThermalNicheData_New$T_Upper[which(ThermalNicheData_New$ThermalGuild == 'temperate')] <- median(Temperate_Difference_Upper) + ThermalNicheData_New$T_Upper_Obs[which(ThermalNicheData_New$ThermalGuild == 'temperate')]
ThermalNicheData_New$T_Upper[which(ThermalNicheData_New$ThermalGuild == 'tropical')]  <- median(Tropical_Difference_Upper)   + ThermalNicheData_New$T_Upper_Obs[which(ThermalNicheData_New$ThermalGuild == 'tropical')]
ThermalNicheData_New$T_Upper[!is.na(ThermalNicheData_New$T_Upper_Mod)] <- ThermalNicheData_New$T_Upper_Mod[!is.na(ThermalNicheData_New$T_Upper_Mod)] 


# Thermal mimumums from occupancy models. 
ThermalMin <- Species_TLowers %>% select(SpeciesName, T_Lower)
ThermalMin$T_Lower_Mod <- ThermalMin$T_Lower; ThermalMin$T_Lower <- NULL
ThermalNicheData_New <- left_join(ThermalNicheData_New, ThermalMin)

# Thermal minimums adjusted from occupancy model outputs. 
ThermalNicheData_New$T_Lower <- NA
ThermalNicheData_New$T_Lower[which(ThermalNicheData_New$ThermalGuild == 'temperate')] <- median(Temperate_Difference_Lower)  + ThermalNicheData_New$T_Lower_Obs[which(ThermalNicheData_New$ThermalGuild == 'temperate')]
ThermalNicheData_New$T_Lower[which(ThermalNicheData_New$ThermalGuild == 'tropical')]  <- median(Tropical_Difference_Lower)   + ThermalNicheData_New$T_Lower_Obs[which(ThermalNicheData_New$ThermalGuild == 'tropical')]
ThermalNicheData_New$T_Lower[!is.na(ThermalNicheData_New$T_Lower_Mod)] <- ThermalNicheData_New$T_Lower_Mod[!is.na(ThermalNicheData_New$T_Lower_Mod)] 

pdf(file = 'figures_extra/CompareUppers_ALL_150102018.pdf', width = 6, height = 4)
ggplot() + 
  geom_point(data = ThermalNicheData_New %>% filter(), aes(x = T_Upper_Obs, y = T_Upper)) + 
  geom_point(data = Species_Tuppers3, aes(x = T_Upper_Obs, y = T_Upper), col = 'blue') + 
  geom_abline() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(10, 34, 1)) + 
  scale_x_continuous(breaks = seq(10, 34, 1)) + 
  xlab('Sampling upper limit') + 
  ylab('Modelled upper limit')
dev.off()

pdf(file = 'figures_extra/CompareLowers_ALL_15012018.pdf', width = 6, height = 4)
ggplot() + 
  geom_point(data = ThermalNicheData_New, aes(x = T_Lower_Obs, y = T_Lower)) + 
  geom_point(data = Species_TLowers3, aes(x = T_Lower_Obs, y = T_Lower), col = 'blue') + 
  geom_abline() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(10, 34, 1)) + 
  scale_x_continuous(breaks = seq(10, 34, 1)) + 
  xlab('Sampling Lower limit') + 
  ylab('Modelled Lower limit')
dev.off()


# Confidence criteria for all thermal niche parameters ----

# These match the workflow in the supporting materials 

# 1. Confidence in upper and lower limit

# Create sampling uppers and sampling lowers
AbsenceLimits <- RLS_All %>% select(SpeciesName, MeanSiteSST_NOAA, Presence) %>% group_by(SpeciesName) %>% 
  do(Absence_TUpper = max(.$MeanSiteSST_NOAA), 
     Absence_TLower = min(.$MeanSiteSST_NOAA)) %>% 
  unnest(Absence_TUpper, Absence_TLower)

# Combine with thermal niche data new
ThermalNicheData_New <- left_join(ThermalNicheData_New, AbsenceLimits)

# Create confidence limits (first are they modelled, second are they within 3°C of absence records limits)
ThermalNicheData_New$Conf_T_Upper <- ifelse(is.na(ThermalNicheData_New$T_Upper_Mod), -1, ifelse(ThermalNicheData_New$T_Upper_Mod - ThermalNicheData_New$Absence_TUpper > 3, -1, 0))
ThermalNicheData_New$Conf_T_Lower <- ifelse(is.na(ThermalNicheData_New$T_Lower_Mod), -1, ifelse(ThermalNicheData_New$T_Lower_Mod - ThermalNicheData_New$Absence_TLower < -3, -1, 0))





# 2. Confidence in T_Breadth 
ThermalNicheData_New$T_Range <- ThermalNicheData_New$T_Upper - ThermalNicheData_New$T_Lower

# Estimate ratio between observed range and estimated range. 
ThermalNicheData_New$T_Range_RATIO <- ThermalNicheData_New$T_Range / (ThermalNicheData_New$T_Upper_Obs - ThermalNicheData_New$T_Lower_Obs)

# We see that the observed T_Range is now rarely much bigger than the modelled T_Range.
hist(unique(ThermalNicheData_New$T_Range_RATIO)) 

# Assign confidence scores. 
ThermalNicheData_New$Conf_T_Breadth <- ifelse(ThermalNicheData_New$T_Range_RATIO <= 2, ifelse(ThermalNicheData_New$T_Range_RATIO <= 2 & ThermalNicheData_New$T_Range_RATIO > 0.9, 0, -1), -1)
ThermalNicheData_New %>% select(ThermalGuild, Conf_T_Breadth, SpeciesName) %>% unique(.) %>% filter(Conf_T_Breadth != 0) %>% .$ThermalGuild %>% table(.) # 
# temperate  tropical 
# 16        32

hist(ThermalNicheData_New$Conf_T_Upper + ThermalNicheData_New$Conf_T_Lower + ThermalNicheData_New$Conf_T_Breadth)

# Use realized (sampling) thermal limits to define Tmin or Tmax with previous addition of difference between confidence 3 species sampling limits and modelled thermal limits. 
ThermalNicheData_New$T_Upper[which((ThermalNicheData_New$Conf_T_Upper == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'tropical' )] <- 
  ThermalNicheData_New$T_Upper_Obs[which((ThermalNicheData_New$Conf_T_Upper == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'tropical' )] + median(Tropical_Difference_Upper)
ThermalNicheData_New$T_Upper[which((ThermalNicheData_New$Conf_T_Upper == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'temperate' )] <- 
  ThermalNicheData_New$T_Upper_Obs[which((ThermalNicheData_New$Conf_T_Upper == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'temperate' )] + median(Temperate_Difference_Upper)

ThermalNicheData_New$T_Lower[which((ThermalNicheData_New$Conf_T_Lower == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'tropical' )] <- 
  ThermalNicheData_New$T_Lower_Obs[which((ThermalNicheData_New$Conf_T_Lower == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'tropical' )] + median(Tropical_Difference_Lower)
ThermalNicheData_New$T_Lower[which((ThermalNicheData_New$Conf_T_Lower == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'temperate' )] <- 
  ThermalNicheData_New$T_Lower_Obs[which((ThermalNicheData_New$Conf_T_Lower == -1 | ThermalNicheData_New$Conf_T_Breadth == -1) & ThermalNicheData_New$ThermalGuild == 'temperate' )] + median(Temperate_Difference_Lower)






# 3. Confidence in qgam estimate of Topt. 
ThermalNicheData_New <- left_join(ThermalNicheData_New, Quantile_Parameters %>% select(-MaxAbundance, -Topt, -T_Opt_Difference_Lower))

# Confidence_Qgam
ThermalNicheData_New$Conf_Qgam <- ifelse(ThermalNicheData_New$T_Gam_pvalue < 0.05, 0, -1)

ThermalNicheData_New %>% select(ThermalGuild, Conf_Qgam, SpeciesName) %>% unique(.) %>% filter(Conf_Qgam != 0) %>% .$ThermalGuild %>% table(.) #
# temperate  tropical 
# 26        66  

# Redefine Topt to midpoint when non-significant or no confidence in upper/lower. 
ThermalNicheData_New$Topt[which(ThermalNicheData_New$Conf_Qgam == -1)] <- ThermalNicheData_New$T_Midpoint_Obs[which(ThermalNicheData_New$Conf_Qgam == -1)]

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits. 
ThermalNicheData_New$T_SD_Upper <- abs(ThermalNicheData_New$T_Upper - ThermalNicheData_New$Topt) / 1.96
ThermalNicheData_New$T_SD_Lower <- abs(ThermalNicheData_New$T_Lower - ThermalNicheData_New$Topt) / 1.96

# Create column estimating skew
ThermalNicheData_New$T_Skew <- ThermalNicheData_New$T_SD_Upper - ThermalNicheData_New$T_SD_Lower




# 4. Confidence in Topt (should not buffer exactly Tmin)
ThermalNicheData_New$T_Opt_Difference_Upper <- ThermalNicheData_New$T_Upper - ThermalNicheData_New$Topt
ThermalNicheData_New$T_Opt_Difference_Lower <- ThermalNicheData_New$Topt    - ThermalNicheData_New$T_Lower

# Confidence_T_Opt_Difference_Upper + Confidence_T_Opt_Difference_Lower
ThermalNicheData_New$Confidence_T_Opt_Difference_Upper <- ifelse(ThermalNicheData_New$T_Opt_Difference_Upper < 1, -1, 0)
ThermalNicheData_New$Confidence_T_Opt_Difference_Lower <- ifelse(ThermalNicheData_New$T_Opt_Difference_Lower < 1, -1, 0)

# Combine all confidence scores
ThermalNicheData_New$ConfidenceCombined <- 3 + 
  ThermalNicheData_New$Conf_T_Upper + 
  ThermalNicheData_New$Conf_T_Lower +
  ThermalNicheData_New$Conf_T_Breadth + 
  ThermalNicheData_New$Conf_Qgam + 
  ThermalNicheData_New$Confidence_T_Opt_Difference_Upper + 
  ThermalNicheData_New$Confidence_T_Opt_Difference_Lower

table(ThermalNicheData_New$ConfidenceCombined)
# 0   1   2   3 
# 33 173 361 137 
# ----

# Save objects needed for analyses ----
# Cleared workspace of redundant objects not needed down the line. 
#RLS_All
#RLS_Temp
#RLS_Trop
#ThermalNicheData_New # The main object created from this script.
save.image(file = 'data_derived/objects-from-script-2.RData')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------



