# Fourth script in thermal-niche RLS analysis for net thermal niche shape analysis.
# This section models the new thermal niche parameters across standardised species ecological performance and thermal ranges . 

# This includes 
# 1. Bayesian models fit in JAGS estimating ecological performance across thermal gradients for 
#      a. P-occ-abun
#      b. P-abun
#      b. P-occ
# 2. The above models fit seperately for tropical and temperate species. 

# Load libraries and packages ---- 
library(R2jags)
library(MCMCvis)
library(coda)

# Load RLS_All from script 2 ----
load('data_derived/objects-from-script-2.RData')
# ----------------------------------------------------------------------------

# DEFINE JAGS MODELS ----
# JAGS TPC Model for abundance and occupancy + abundance data ----
# Step 3: Formulate JAGS modelling code
sink("TPC.txt")
cat("
    model{
    
    # Priors for TPC parameters
    T_Opt ~ dnorm(0, 0.001)
    T_SD  ~ dunif(0, 10)
    T_SD_2  ~ dunif(0, 10)
    # T_Upper ~ dnorm(0, 0.001)I(0,10)
    MaxAbun ~ dnorm(1, 0.001)I(0,1)
    Var ~ dnorm(1, 0.001)
    
    # TPC Model
    for(i in 1:nTemp){ 	## loop through sites
    Abundance[i] ~ dnorm(y[i], Var)
    
    # Temp.hi[i]   <-  step(Temp[i]-T_Opt) # Define your breakpoint
    #  y[i]         <-  (1-Temp.hi[i]) *   MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    #                   (Temp.hi[i])   *   MaxAbun*(1-((Temp[i] - T_Opt) / (T_Opt - T_Upper))^2)   # Model above Topt
    
    Temp.hi[i]   <-  step(Temp[i]-T_Opt) # Define your breakpoint
    y[i]         <-  (1-Temp.hi[i])   *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi[i])     *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    
    }
    
    # Solve for Topt
    #    YMax <- MaxAbun*(exp(-(((T_Opt - T_Opt) / (T_SD)))^2))
    
    # Make performance estimates from parameters 
    for(j in 1:nPerformance){
    Temp.hi2[j]   <-  step(Temp_pred[j]-T_Opt) # Define your breakpoint
    Performance[j] <- (1-Temp.hi2[j])   *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi2[j])     *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    
    }
    
    # Estimate skew 
    Skew <- T_SD_2 - T_SD
    
    # Estimate niche breadth in sd 
    T_Breadth <- (((T_SD + T_SD_2))/2) 
    
    #3. Discrepancy measures: 
    #   Pearson residuals PRes and ordinary residuals E
    y_hat <- mean(Abundance)
    
    for (i in 1:nTemp) {
    # VarY[i] <- y[i] * (1 - y[i])
    # PRes[i] <- (Abundance[i] - y[i]) / sqrt(VarY[i])   
    S_res[i]    <- (Abundance[i] - y[i])^2
    S_tot[i]    <- (Abundance[i] - y_hat)^2
    }
    
    R_squared <- 1 - (sum(S_res) / sum(S_tot))
    
    }
    ",fill = TRUE)
sink()


# JAGS TPC Model with logit link for occupancy only model ----
sink("TPC_Occupancy.txt")
cat("
    model{
    
    # Priors for TPC parameters
    T_Opt ~ dnorm(0, 0.001)
    T_SD  ~ dunif(0, 10)
    T_SD_2  ~ dunif(0, 10)
    MaxAbun ~ dnorm(1, 0.001)I(0,1)
    # Var ~ dnorm(1, 0.001)
    
    # TPC Model
    for(i in 1:nTemp){ 	# loop through sites
    Occupancy[i] ~ dbern(y[i])
    Temp.hi[i]  <-  step(Temp[i]-T_Opt) # Define your breakpoint
    logit(y[i]) <-  (1-Temp.hi[i])   *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi[i])     *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    }
    
    # Make performance estimates from parameters 
    for(j in 1:nPerformance){
    Temp.hi2[j]   <-  step(Temp_pred[j]-T_Opt) # Define your breakpoint
    Performance[j] <- (1-Temp.hi2[j])   *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi2[j])     *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    }
    
    # Estimate skew 
    Skew <- T_SD_2 - T_SD
    
    # Estimate niche breadth in sd 
    T_Breadth <- (((T_SD + T_SD_2))/2)
    
    #3. Discrepancy measures: 
    #   Pearson residuals PRes and ordinary residuals E
    for (i in 1:nTemp) {
    E[i]    <- Occupancy[i]  - y[i]
    }
    }
    ",fill = TRUE)
sink()




# ----------------------------------------------------------------------------


# DEFINE FUNCTIONS TO RUN AND HANDLE OUTPUTS OF JAGS MODELS ----
# Function to fit JAGS models to quantiles ----
# Function to fit jags models to defined quantiles on certain data and then save with predefined colour scale for plotting. 
RunJAGSModels <- function(modeldata,
                          params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var"),
                          Quantile = c(0.99),
                          Subset = NULL, # Define which subset we are talking about. LABEL ONLY
                          ThermalGuild = 'temperate',# Define which subset we are talking about
                          Scale = 'Species', # Produces different aggregations of data depending on intput here. 
                          n.thin     = 10, 
                          n.chains   = 3,
                          n.burnin   = 400,
                          n.iter     = 500
){
  
  # Define inputs to for loop.
  QuantileDataList <- list()
  JagsDataList <- list()
  TPC_Model <- list()
  PlotDataJags <- list()
  inits  <- function(){list(T_Opt = 0, T_SD = 1, T_SD_2 = 1, T_Upper = 1)}
  
  # Run for loop to estimate model over each level of quantile. 
  for(i in 1:length(Quantile)){
    
    # Create data. 
    #QuantileDataList[[i]] <- RLS_All_HighConf_Temp %>% group_by(Species_MeanSiteSST_NOAA_V2) %>% do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% unnest(ScaledLogAbundanceAdult40) %>% na.omit()
    
    
    if(Scale == 'Species'){
      # Aggragte by species and temperature groups (species have equal contributions)
      QuantileDataList[[i]] <- modeldata %>% 
        
        # Take quantile for each species and temperature group
        group_by(Species_MeanSiteSST_NOAA_V2, SpeciesName, ThermalGuild) %>% 
        do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% 
        unnest(ScaledLogAbundanceAdult40) %>% 
        group_by(Species_MeanSiteSST_NOAA_V2, ThermalGuild) #%>% 
      
      # Take average quantile across all species
      #do(ScaledLogAbundanceAdult40    = mean(.$ScaledLogAbundanceAdult40, na.rm = T), 
      #   ScaledLogAbundanceAdult40_SD = sd(.$ScaledLogAbundanceAdult40, na.rm = T)) %>% 
      #unnest(ScaledLogAbundanceAdult40, ScaledLogAbundanceAdult40_SD) 
      
    }else{
      
      # Aggragte over all samples (species have different contributions) 
      QuantileDataList[[i]] <- modeldata %>% 
        
        # Take quantile for each species and temperature group
        group_by(Species_MeanSiteSST_NOAA_V2) %>% 
        do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% 
        unnest(ScaledLogAbundanceAdult40)
      
    }
    
    # Create jags data list. 
    JagsDataList[[i]]  <- list(Abundance     = QuantileDataList[[i]]$ScaledLogAbundanceAdult40,
                               Temp          = QuantileDataList[[i]]$Species_MeanSiteSST_NOAA_V2,
                               nTemp         = nrow(QuantileDataList[[i]]), 
                               Temp_pred     = seq(min(QuantileDataList[[i]]$Species_MeanSiteSST_NOAA_V2), max(QuantileDataList[[i]]$Species_MeanSiteSST_NOAA_V2), length.out = 500),
                               nPerformance  = 500)
    
    # Run the jags models
    TPC_Model[[i]]   <- jags(data       = JagsDataList[[i]],
                             inits      = inits,
                             parameters = params,
                             model      = "TPC.txt",
                             n.thin     = n.thin, 
                             n.chains   = n.chains,
                             n.burnin   = n.burnin,
                             n.iter     = n.iter)
    
    # Store the outputs plotlines
    PlotDataJags[[i]] <- data.frame(
      Performance       = TPC_Model[[i]]$BUGSoutput$mean$y,
      Performance_Upper = apply(TPC_Model[[i]]$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.975)),
      Performance_Lower = apply(TPC_Model[[i]]$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.025)), 
      Quantile = as.factor(Quantile[i]), 
      Species_MeanSiteSST_NOAA_V2 = QuantileDataList[[i]]$Species_MeanSiteSST_NOAA_V2, 
      ScaledLogAbundanceAdult40   = QuantileDataList[[i]]$ScaledLogAbundanceAdult40,
      ThermalGuild = ThermalGuild)
    if(is.null(Subset)){NULL}else{PlotDataJags[[i]]$Subset = Subset[i]}
    
  }
  
  ## Aggregate data output for plotting
  # Integrate with data to form ggplot. 
  PlotDataJags <- do.call(rbind, PlotDataJags)
  
  return(list(data = PlotDataJags, models = TPC_Model))
  
}

# Function to fit JAGS mdoel to occupancy patterns using a logit link function for same mathmatical function describing TPC. 
RunJAGSModels_OCC <- function(modeldata,
                              params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var"),
                              Quantile = c(0.99),
                              Subset = NULL, # Define which subset we are talking about. LABEL ONLY
                              ThermalGuild = 'temperate',# Define which subset we are talking about
                              Scale = 'Species', # Produces different aggregations of data depending on intput here. 
                              n.thin     = 10, 
                              n.chains   = 3,
                              n.burnin   = 400,
                              n.iter     = 500
){
  
  # Define inputs to for loop.
  i = 1
  PresentModelData <- data.frame()
  JagsDataList <- list()
  TPC_Model <- data.frame()
  PlotDataJags <- list()
  inits  <- function(){list(T_Opt = 0, T_SD = 1, T_SD_2 = 1, T_Upper = 1)}
  
  # Aggregate by species and temperature groups (species have equal contributions) and take mean of occupancies 
  
  #DataTest <-  modeldata %>% filter(SpeciesName == 'Abudefduf troschelii')
  
  if(Scale == 'Species'){
    PresentModelData <- modeldata %>% 
      # Take mean for each species occupancy rates
      group_by(Species_MeanSiteSST_NOAA_V2, SpeciesName, ThermalGuild) %>% 
      do(Mean_SpeciesPresence = mean(.$Presence, na.rm = T)) %>% 
      unnest(Mean_SpeciesPresence) }else{
        
        PresentModelData <- modeldata %>% 
          # Take mean for each species occupancy rates
          group_by(Species_MeanSiteSST_NOAA_V2, ThermalGuild) %>% 
          do(Mean_SpeciesPresence = mean(.$Presence, na.rm = T)) %>% 
          unnest(Mean_SpeciesPresence) 
      }
  
  
  # Create jags data list. 
  JagsDataList[[i]]  <- list(Abundance = PresentModelData$Mean_SpeciesPresence,
                             Temp          = PresentModelData$Species_MeanSiteSST_NOAA_V2,
                             nTemp         = nrow(PresentModelData), 
                             Temp_pred     = seq(min(PresentModelData$Species_MeanSiteSST_NOAA_V2), max(PresentModelData$Species_MeanSiteSST_NOAA_V2), length.out = 500),
                             nPerformance  = 500)
  
  # Run the jags models
  TPC_Model   <- jags(data       = JagsDataList[[i]],
                      inits      = inits,
                      parameters = params,
                      model      = "TPC.txt",
                      n.thin     = n.thin, 
                      n.chains   = n.chains,
                      n.burnin   = n.burnin,
                      n.iter     = n.iter)
  
  # Store the outputs plotlines
  PlotDataJags[[i]] <- data.frame(
    Performance       = TPC_Model$BUGSoutput$mean$y,
    Performance_Upper = apply(TPC_Model$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.975)),
    Performance_Lower = apply(TPC_Model$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.025)), 
    Species_MeanSiteSST_NOAA_V2 = PresentModelData$Species_MeanSiteSST_NOAA_V2, 
    Mean_SpeciesPresence   = PresentModelData$Mean_SpeciesPresence,
    ThermalGuild = ThermalGuild)
  
  if(is.null(Subset)){NULL}else{PlotDataJags[[i]]$Subset = Subset}
  
  
  
  ## Aggregate data output for plotting
  # Integrate with data to form ggplot. 
  PlotDataJags <- do.call(rbind, PlotDataJags)
  
  return(list(data = PlotDataJags, models = TPC_Model))
  
}


# Function to obtain the posterior distributions of defined parameters. 
ObtainPosterior <- function(BugsList = list(), ModelData = list(), Parameter = character(), Subset = character(), colour = character()){
  
  OutputPosteriorData  <- list()
  for(i in 1:length(BugsList)){
    
    OutputPosteriorData[[i]]  <- data.frame(Value = BugsList[[i]]$sims.list[[Parameter]], Parameter = Parameter, ThermalGuild = unique(ModelData[[i]][['ThermalGuild']]), Subset = Subset[i], colour = colour[i])
    
  }
  
  OutputPosteriorData <- do.call('rbind', OutputPosteriorData)
  return(OutputPosteriorData)
  
}




# ----------------------------------------------------------------------------







# STANDARDISE DATA ACROSS ECOLOGICAL PERFORMANCE AND TEMPERATURE GRADIENTS WITHIN EACH SPECIES ---- 
# Run data standaridisations ----
# Bind in thermal niche parameters to standardised species temperature and abundance ranges. 
RLS_All_JAGS <- RLS_All %>% left_join(., ThermalNicheData_New %>% select(SpeciesName, Topt, T_SD_Upper, T_SD_Lower))

# Centre temperature data within species
RLS_All_JAGS_scales <- RLS_All_JAGS %>% 
  group_by(SpeciesName, ThermalGuild) %>% 
  nest() %>% 
  
  # Scale and transpose data to be comparable across
  mutate(Species_MeanSiteSST_NOAA  = purrr::map(data, ~ (.$MeanSiteSST_NOAA - unique(.$Topt))/mean(c(.$T_SD_Upper , .$T_SD_Lower))), 
         ScaledLogAbundanceAdult40 = purrr::map(data, ~ log(.$AbundanceAdult40+1)/max(log(.$AbundanceAdult40+1)))) %>% 
  unnest(Species_MeanSiteSST_NOAA, ScaledLogAbundanceAdult40)

# Assign columns to data. 
RLS_All_JAGS$Species_MeanSiteSST_NOAA     <- RLS_All_JAGS_scales$Species_MeanSiteSST_NOAA
RLS_All_JAGS$ScaledLogAbundanceAdult40    <- RLS_All_JAGS_scales$ScaledLogAbundanceAdult40
RLS_All_JAGS$T_Opt                        <- RLS_All_JAGS$Topt
RLS_All_JAGS$T_SD_2                       <- RLS_All_JAGS$T_SD_Upper
RLS_All_JAGS$T_SD                         <- RLS_All_JAGS$T_SD_Lower

# Fit one model to each thermal guild. 
TempData_JAGS <- RLS_All_JAGS %>% filter(ThermalGuild == 'temperate') #%>% filter(AbundanceAdult40 > 0)
TropData_JAGS <- RLS_All_JAGS %>% filter(ThermalGuild == 'tropical')  #%>% filter(AbundanceAdult40 > 0)
AllData_JAGS <- rbind(TropData_JAGS, TempData_JAGS)


# ----------------------------------------------------------------------------

# FIT MODELS FOR P-Occ-Abun ---- 
# AllJAGS_0.99  ----
AllJAGS_0.99 <- RunJAGSModels(modeldata = AllData_JAGS,
                              params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                              Quantile = c(0.99),
                              Subset = 'Quantile 0.99', # Define which subset we are talking about
                              ThermalGuild = 'tropical', 
                              Scale = 'Species', 
                              n.thin     = 10, 
                              n.chains   = 4,
                              n.burnin   = 100,
                              n.iter     = 500)
saveRDS(AllJAGS_0.99, 'data_derived/jagsmodels/AllJAGS_0.99.rds')

# Summary statistic for above model 
round(AllJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(AllJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(AllJAGS_0.99$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(AllJAGS_0.99$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(AllJAGS_0.99$models[[1]]$BUGSoutput$summary['R_squared',], 3)
AllJAGS_0.99$models[[1]]$BUGSoutput$DIC
AllJAGS_0.99$models[[1]]$BUGSoutput$pD

# Caterpillar plots for above models
AllJAGS_0.99_mcmc <- AllJAGS_0.99$models[[1]]
AllJAGS_0.99_mcmc <- MCMCchains(AllJAGS_0.99_mcmc, params = c("T_Opt", "T_SD", "T_SD_2", "MaxAbun", "Skew", "T_Breadth"), mcmc.list = T)
MCMCtrace(AllJAGS_0.99_mcmc)


# TropJAGS_0.99 ----
TropJAGS_0.99 <- RunJAGSModels(modeldata = TropData_JAGS,
                               params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                               Quantile = c(0.99),
                               Subset = 'Quantile 0.99', # Define which subset we are talking about
                               ThermalGuild = 'tropical', 
                               Scale = 'Species', 
                               n.thin     = 10, 
                               n.chains   = 4,
                               n.burnin   = 7500,
                               n.iter     = 10000)
saveRDS(TropJAGS_0.99, 'data_derived/jagsmodels/TropJAGS_0.99.rds')

# Summary statistic for above model 
round(TropJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TropJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TropJAGS_0.99$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TropJAGS_0.99$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TropJAGS_0.99$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TropJAGS_0.99$models[[1]]$BUGSoutput$DIC
TropJAGS_0.99$models[[1]]$BUGSoutput$pD

# Caterpillar plots for above models
TropJAGS_0.99_mcmc <- TropJAGS_0.99$models[[1]]
TropJAGS_0.99_mcmc <- MCMCchains(TropJAGS_0.99_mcmc, params = c("T_Opt", "T_SD", "T_SD_2", "MaxAbun", "Skew", "T_Breadth"), mcmc.list = T)
MCMCtrace(TropJAGS_0.99_mcmc)

# TempJAGS_0.99 ----
TempJAGS_0.99 <- RunJAGSModels(modeldata = TempData_JAGS,
                               params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                               Quantile = c(0.99),
                               Subset = 'Quantile 0.99', # Define which subset we are talking about
                               ThermalGuild = 'Temperate', 
                               Scale = 'Species', 
                               n.thin     = 10, 
                               n.chains   = 4,
                               n.burnin   = 7500,
                               n.iter     = 10000)
saveRDS(TempJAGS_0.99, 'data_derived/jagsmodels/TempJAGS_0.99.rds')

# Summary statistic for above model 
round(TempJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TempJAGS_0.99$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TempJAGS_0.99$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TempJAGS_0.99$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TempJAGS_0.99$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TempJAGS_0.99$models[[1]]$BUGSoutput$DIC
TempJAGS_0.99$models[[1]]$BUGSoutput$pD

# Caterpillar plots for above models
TempJAGS_0.99_mcmc <- TempJAGS_0.99$models[[1]]
TempJAGS_0.99_mcmc <- MCMCchains(TempJAGS_0.99_mcmc, params = c("T_Opt", "T_SD", "T_SD_2", "MaxAbun", "Skew", "T_Breadth"), mcmc.list = T)
MCMCtrace(TempJAGS_0.99_mcmc)


# AllJAGS_0.99 aggregation  ----
AllJAGS_0.99_Aggregated <- RunJAGSModels(modeldata = rbind(TropData_JAGS, TempData_JAGS),
                                         params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                         Quantile = c(0.99),
                                         Subset = 'Quantile 0.99', # Define which subset we are talking about
                                         ThermalGuild = 'tropical', 
                                         Scale = 'Aggregated', 
                                         n.thin     = 10, 
                                         n.chains   = 4,
                                         n.burnin   = 7500,
                                         n.iter     = 10000)
saveRDS(AllJAGS_0.99_Aggregated , 'data_derived/jagsmodels/AllJAGS_0.99_Aggregated.rds')
round(AllJAGS_0.99_Aggregated$models[[1]]$BUGSoutput$summary['R_squared',], 3)

# TropJAGS_0.99 aggregation ----
TropJAGS_0.99_Aggregated <- RunJAGSModels(modeldata = TropData_JAGS,
                                          params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                          Quantile = c(0.99),
                                          Subset = 'Quantile 0.99', # Define which subset we are talking about
                                          ThermalGuild = 'tropical', 
                                          Scale = 'Aggregated', 
                                          n.thin     = 10, 
                                          n.chains   = 4,
                                          n.burnin   = 7500,
                                          n.iter     = 10000)

saveRDS(TropJAGS_0.99_Aggregated , 'data_derived/jagsmodels/TropJAGS_0.99_Aggregated.rds')
TropJAGS_0.99_Aggregated$models[[1]]$BUGSoutput$summary['R_squared',]
TempJAGS_0.99_Aggregated$models[[1]]$BUGSoutput$summary['R_squared',]

# TempJAGS_0.99 aggregation ----
TempJAGS_0.99_Aggregated  <- RunJAGSModels(modeldata = TempData_JAGS,
                                           params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                           Quantile = c(0.99),
                                           Subset = 'Quantile 0.99', # Define which subset we are talking about
                                           ThermalGuild = 'temperate', 
                                           Scale = 'Aggregated', 
                                           n.thin     = 10, 
                                           n.chains   = 4,
                                           n.burnin   = 7500,
                                           n.iter     = 10000)

saveRDS(TempJAGS_0.99_Aggregated , 'data_derived/jagsmodels/TempJAGS_0.99_Aggregated.rds')
TropJAGS_0.99_Aggregated$models[[1]]$BUGSoutput$summary['R_squared',]
TempJAGS_0.99_Aggregated$models[[1]]$BUGSoutput$summary['R_squared',]

# ----------------------------------------------------------------------------

# FIT MODELS FOR P-Abun ---- 
# AllJAGS_0.99_ABUN  ---- 
AllJAGS_0.99_ABUN <- RunJAGSModels(modeldata = AllJags_Data %>% filter(AbundanceAdult40 > 0),
                                   params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                   Quantile = c(0.99),
                                   Subset = 'Quantile 0.99', # Define which subset we are talking about
                                   ThermalGuild = 'tropical', 
                                   Scale = 'Species', 
                                   n.thin     = 10, 
                                   n.chains   = 4,
                                   n.burnin   = 7500,
                                   n.iter     = 10000)
saveRDS(AllJAGS_0.99_ABUN, 'data_derived/jagsmodels/AllJAGS_0.99_ABUN.rds')

round(AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['R_squared',], 3)
AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$DIC
AllJAGS_0.99_ABUN$models[[1]]$BUGSoutput$pD

# TropJAGS_0.99_ABUN ----
TropJAGS_0.99_ABUN <- RunJAGSModels(modeldata = TropData_JAGS %>% filter(AbundanceAdult40 > 0),
                                    params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                    Quantile = c(0.99),
                                    Subset = 'Quantile 0.99', # Define which subset we are talking about
                                    ThermalGuild = 'tropical', 
                                    Scale = 'Species', 
                                    n.thin     = 10, 
                                    n.chains   = 4,
                                    n.burnin   = 7500,
                                    n.iter     = 10000)
saveRDS(TropJAGS_0.99_ABUN, 'data_derived/jagsmodels/TropJAGS_0.99_ABUN.rds')

round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$DIC
TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$pD

# TempJAGS_0.99_ABUN ----
TempJAGS_0.99_ABUN <- RunJAGSModels(modeldata = TempData_JAGS %>% filter(AbundanceAdult40 > 0),
                                    params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                    Quantile = c(0.99),
                                    Subset = 'Quantile 0.99', # Define which subset we are talking about
                                    ThermalGuild = 'temperate', 
                                    Scale = 'Species', 
                                    n.thin     = 10, 
                                    n.chains   = 4,
                                    n.burnin   = 7500,
                                    n.iter     = 10000)
saveRDS(TempJAGS_0.99_ABUN, 'data_derived/jagsmodels/TempJAGS_0.99_ABUN.rds')

round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$DIC
TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$pD


# AllJAGS_0.99_ABUN_AGG  ----
AllJAGS_0.99_ABUN_AGG <- RunJAGSModels(modeldata = AllJags_Data %>% filter(AbundanceAdult40 > 0),
                                       params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                       Quantile = c(0.99),
                                       Subset = 'Quantile 0.99', # Define which subset we are talking about
                                       ThermalGuild = 'tropical', 
                                       Scale = 'Aggregated', 
                                       n.thin     = 10, 
                                       n.chains   = 4,
                                       n.burnin   = 7500,
                                       n.iter     = 10000)
saveRDS(AllJAGS_0.99_ABUN_AGG, 'data_derived/jagsmodels/AllJAGS_0.99_ABUN_AGG.rds')

# TropJAGS_0.99_ABUN_AGG ----
TropJAGS_0.99_ABUN_AGG <- RunJAGSModels(modeldata = TropData_JAGS %>% filter(AbundanceAdult40 > 0),
                                        params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                        Quantile = c(0.99),
                                        Subset = 'Quantile 0.99', # Define which subset we are talking about
                                        ThermalGuild = 'tropical', 
                                        Scale = 'Aggregated', 
                                        n.thin     = 10, 
                                        n.chains   = 4,
                                        n.burnin   = 7500,
                                        n.iter     = 10000)
saveRDS(TropJAGS_0.99_ABUN_AGG, 'data_derived/jagsmodels/TropJAGS_0.99_ABUN_AGG.rds')


# TempJAGS_0.99_ABUN_AGG ----
TempJAGS_0.99_ABUN_AGG <- RunJAGSModels(modeldata = TempData_JAGS %>% filter(AbundanceAdult40 > 0),
                                        params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                        Quantile = c(0.99),
                                        Subset = 'Quantile 0.99', # Define which subset we are talking about
                                        ThermalGuild = 'temperate', 
                                        Scale = 'Aggregated', 
                                        n.thin     = 10, 
                                        n.chains   = 4,
                                        n.burnin   = 7500,
                                        n.iter     = 10000)
saveRDS(TempJAGS_0.99_ABUN_AGG, 'data_derived/jagsmodels/TempJAGS_0.99_ABUN_AGG.rds')


# ----------------------------------------------------------------------------

# FIT MODELS FOR P-Occ ----
# AllJAGS_OCC  ----
AllJAGS_OCC <- RunJAGSModels_OCC(modeldata = AllJags_Data,
                                 params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                 Subset = 'Occupancy', # Define which subset we are talking about
                                 ThermalGuild = 'temperate', 
                                 Scale = 'Species', 
                                 n.thin     = 10, 
                                 n.chains   = 4,
                                 n.burnin   = 7500,
                                 n.iter     = 10000)
saveRDS(AllJAGS_OCC, 'data_derived/jagsmodels/AllJAGS_OCC.rds')


round(AllJAGS_OCC$models$BUGSoutput$summary['T_Opt',], 3)
round(AllJAGS_OCC$models$BUGSoutput$summary['T_Breadth',], 3)
round(AllJAGS_OCC$models$BUGSoutput$summary['Skew',], 3)
round(AllJAGS_OCC$models$BUGSoutput$summary['MaxAbun',], 3)
round(AllJAGS_OCC$models$BUGSoutput$summary['R_squared',], 3)
AllJAGS_OCC$models$BUGSoutput$DIC
AllJAGS_OCC$models$BUGSoutput$pD

# TropJAGS_OCC ----
TropJAGS_OCC <- RunJAGSModels_OCC(modeldata = TropData_JAGS,
                                  params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                  Subset = 'Occupancy', # Define which subset we are talking about
                                  ThermalGuild = 'tropical', 
                                  Scale = 'Species', 
                                  n.thin     = 10, 
                                  n.chains   = 4,
                                  n.burnin   = 7500,
                                  n.iter     = 10000)
saveRDS(TropJAGS_OCC, 'data_derived/jagsmodels/TropJAGS_OCC.rds')


round(TropJAGS_OCC$models$BUGSoutput$summary['T_Opt',], 3)
round(TropJAGS_OCC$models$BUGSoutput$summary['T_Breadth',], 3)
round(TropJAGS_OCC$models$BUGSoutput$summary['Skew',], 3)
round(TropJAGS_OCC$models$BUGSoutput$summary['MaxAbun',], 3)
round(TropJAGS_OCC$models$BUGSoutput$summary['R_squared',], 3)
TropJAGS_OCC$models$BUGSoutput$DIC
TropJAGS_OCC$models[[1]]$BUGSoutput$pD

# TempJAGS_OCC ----
TempJAGS_OCC <- RunJAGSModels_OCC(modeldata = TempData_JAGS,
                                  params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                  Subset = 'Occupancy', # Define which subset we are talking about
                                  ThermalGuild = 'temperate', 
                                  Scale = 'Species', 
                                  n.thin     = 10, 
                                  n.chains   = 4,
                                  n.burnin   = 7500,
                                  n.iter     = 10000)
saveRDS(TempJAGS_OCC, 'data_derived/jagsmodels/TempJAGS_OCC.rds')

round(TempJAGS_OCC$models$BUGSoutput$summary['T_Opt',], 3)
round(TempJAGS_OCC$models$BUGSoutput$summary['T_Breadth',], 3)
round(TempJAGS_OCC$models$BUGSoutput$summary['Skew',], 3)
round(TempJAGS_OCC$models$BUGSoutput$summary['MaxAbun',], 3)
round(TempJAGS_OCC$models$BUGSoutput$summary['R_squared',], 3)
TempJAGS_OCC$models$BUGSoutput$DIC
TempJAGS_OCC$models$BUGSoutput$pD


# AllJAGS_OCC_AGG ----
AllJAGS_OCC_AGG <- RunJAGSModels_OCC(modeldata = AllJags_Data,
                                     params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                     Subset = 'Occupancy', # Define which subset we are talking about
                                     ThermalGuild = 'temperate', 
                                     Scale = 'Aggregated', 
                                     n.thin     = 10, 
                                     n.chains   = 4,
                                     n.burnin   = 7500,
                                     n.iter     = 10000)
saveRDS(AllJAGS_OCC_AGG, 'data_derived/jagsmodels/AllJAGS_OCC_AGG.rds')

round(AllJAGS_OCC_AGG$models$BUGSoutput$summary['R_squared',], 4)

# TropJAGS_OCC_AGG ----
TropJAGS_OCC_AGG <- RunJAGSModels_OCC(modeldata = TropData_JAGS,
                                      params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                      Subset = 'Occupancy', # Define which subset we are talking about
                                      ThermalGuild = 'tropical', 
                                      Scale = 'Aggregated', 
                                      n.thin     = 10, 
                                      n.chains   = 4,
                                      n.burnin   = 7500,
                                      n.iter     = 10000)
saveRDS(TropJAGS_OCC_AGG, 'data_derived/jagsmodels/TropJAGS_OCC_AGG.rds')

round(TropJAGS_OCC_AGG$models$BUGSoutput$summary['R_squared',], 4)

# TempJAGS_OCC_AGG ----
TempJAGS_OCC_AGG <- RunJAGSModels_OCC(modeldata = TempData_JAGS,
                                      params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                      Subset = 'Occupancy', # Define which subset we are talking about
                                      ThermalGuild = 'temperate', 
                                      Scale = 'Aggregated', 
                                      n.thin     = 10, 
                                      n.chains   = 4,
                                      n.burnin   = 7500,
                                      n.iter     = 10000)
saveRDS(TempJAGS_OCC_AGG, 'data_derived/jagsmodels/TempJAGS_OCC_AGG.rds')

round(TempJAGS_OCC_AGG$models$BUGSoutput$summary['R_squared',], 4)




# ----------------------------------------------------------------------------

