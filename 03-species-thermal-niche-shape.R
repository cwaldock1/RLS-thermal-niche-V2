# Third script in thermal-niche RLS analysis for species level analysis.
# This section models the thermal niche parameters derived in script 2. 

# This includes 
# 1. Tests of abundance declines at thermal niche edges (X^2 tests)
# 2. Analysis of thermal niche skew vs. Topt (figure 2) 
# 3. 
# 4. 

# Initiated 07/05/2018
# Author: Conor Waldock

# Load libraries ----
library(lme4)
library(MuMIn)
library(remef)
library(stargazer)
# ----------------------------------------------------------------------------


# ANALYSIS OF SHAPE PATTERNS IN QGAMS AND 'ECOLOGICAL PERFORMANCE CURVES' ----
# Extract model runs and predict curves from models ----

# Finds the most recent model runs. 
QuantileModels <- list.files('data_derived/qgam_outputs')

# Which species' have confidence scores of > 3. 
HighConfidenceSpecies <- ThermalNicheData_New$SpeciesName[which(ThermalNicheData_New$ConfidenceCombined == 3)]
HighConfidenceSpecies2 <- gsub(' ', '_', ThermalNicheData_New$SpeciesName[which(ThermalNicheData_New$ConfidenceCombined == 3)])

# Extract the location within the vector for each model of a high-confidence species. 
QuantileModels_Species <- do.call(rbind, lapply(HighConfidenceSpecies2, function(x) grep(x, QuantileModels)))

# Subsets list of all models and returns file names of all relevant models. 
qgamFiles_Conf3 <- QuantileModels[QuantileModels_Species]
qgamFiles_Conf3 <- paste('data_derived/qgam_outputs/', qgamFiles_Conf3, sep = '')

# Create list of model outputs for all species. 
Conf3_Gams <- list()
for(i in 1:length(qgamFiles_Conf3)){ 
  Conf3_Gams[[i]] <- readRDS(file = qgamFiles_Conf3[i])
}

# Extract predicitons from models and combine together (and deviance explained)
Conf3_Gams_Preds <- do.call(rbind, lapply(Conf3_Gams, function(x) x[2][[1]]))






# Taken predictions and define abundance groups for both 25% and 50% thresholds ----

# Define categories
Category_25 <- c()
Abundance_percent_cool <- c()
Abundance_percent_warm <- c()
TG <- c() # Thermal guild for a subset. 
for(i in 1:length(unique(HighConfidenceSpecies))){
  
  QGAMfunction <- Conf3_Gams_Preds %>% filter(SpeciesName == HighConfidenceSpecies[i])
  RLS_Species <- RLS_All %>% filter(SpeciesName == HighConfidenceSpecies[i])
  ThermalNicheData_Spp <- ThermalNicheData_New %>% filter(SpeciesName == HighConfidenceSpecies[i])
  
  # Extract abundance and scale
  Abundance <- QGAMfunction$Abundance
  Temp <- QGAMfunction$MeanSiteSST_NOAA
  scaleAbundance <- Abundance/max(Abundance)
  
  TG[i] <- ThermalNicheData_Spp$ThermalGuild
  
  # No trend: mean percentage decline from peak model performance from 
  MaxRaw <- log(ThermalNicheData_Spp$MaxAbundance + 1)
  Warm_Drop <-1 - ( Abundance[100] / max(Abundance) )
  Cool_Drop <-1 - ( Abundance[1] / max(Abundance) )
  NoTrend_Percentage <- max(Cool_Drop, Warm_Drop)
  NoTrend <- ifelse(NoTrend_Percentage > 0.25, FALSE, TRUE)
  
  # Ramped edges
  WarmEdge <- 1 - scaleAbundance[100] # Amount of decline from max abun
  ColdEdge <- 1 - scaleAbundance[1]   # Amount of decline from max abun
  
  Abundance_percent_warm[i] <- scaleAbundance[100] # Amount of decline from max abun
  Abundance_percent_cool[i] <- scaleAbundance[1]   # Amount of decline from max abun
  
  # Test is edge falls by at least 25% of maximum. 
  ColdEdgeRamp <- ifelse(ColdEdge > 0.25, FALSE, TRUE)
  WarmEdgeRamp <- ifelse(WarmEdge > 0.25, FALSE, TRUE)
  
  # Test if both edges are ramped and therefore are at the centre. 
  AbundanceCentre <- ifelse(ColdEdgeRamp == F & WarmEdgeRamp == F, TRUE, FALSE)
  
  AbundanceEdge   <- ifelse(AbundanceCentre == F & NoTrend == F, TRUE, FALSE)
  
  Categories <- c('No trend', 'Ramped cool-edge', 'Ramped warm-edge', 'Abundance centre', 'Abundance edge')
  
  if(NoTrend == T){ 
    Category_25[i] <- 'No trend'
  }else{
    Category_25[i] <- Categories[which(c(NoTrend, ColdEdgeRamp, WarmEdgeRamp, AbundanceCentre, AbundanceEdge))]
  }
}

# Define categories
Category_50 <- c()
TG_50 <- c() # Thermal guild for a subset. 
for(i in 1:length(unique(HighConfidenceSpecies))){
  
  QGAMfunction <- Conf3_Gams_Preds %>% filter(SpeciesName == HighConfidenceSpecies[i])
  RLS_Species <- RLS_All %>% filter(SpeciesName == HighConfidenceSpecies[i])
  ThermalNicheData_Spp <- ThermalNicheData_New %>% filter(SpeciesName == HighConfidenceSpecies[i])
  
  # Extract abundance and scale
  Abundance <- QGAMfunction$Abundance
  Temp <- QGAMfunction$MeanSiteSST_NOAA
  scaleAbundance <- Abundance/max(Abundance)
  
  TG_50[i] <- ThermalNicheData_Spp$ThermalGuild
  
  # No trend: mean percentage decline from peak model performance from 
  MaxRaw <- log(ThermalNicheData_Spp$MaxAbundance + 1)
  Warm_Drop <-1 - ( Abundance[100] / max(Abundance) )
  Cool_Drop <-1 - ( Abundance[1] / max(Abundance) )
  NoTrend_Percentage <- max(Cool_Drop, Warm_Drop)
  NoTrend <- ifelse(NoTrend_Percentage > 0.5, FALSE, TRUE)
  
  # Ramped edges
  WarmEdge <- 1 - scaleAbundance[100] # Amount of decline from max abun
  ColdEdge <- 1 - scaleAbundance[1]   # Amount of decline from max abun
  
  # Test is edge falls by at least 25% of maximum. 
  ColdEdgeRamp <- ifelse(ColdEdge > 0.25, FALSE, TRUE)
  WarmEdgeRamp <- ifelse(WarmEdge > 0.25, FALSE, TRUE)
  
  # Test if both edges are ramped and therefore are at the centre. 
  AbundanceCentre <- ifelse(ColdEdgeRamp == F & WarmEdgeRamp == F, TRUE, FALSE)
  
  AbundanceEdge   <- ifelse(AbundanceCentre == F & NoTrend == F, TRUE, FALSE)
  
  Categories <- c('No trend', 'Ramped cool-edge', 'Ramped warm-edge', 'Abundance centre', 'Abundance edge')
  
  if(NoTrend == T){ 
    Category_50[i] <- 'No trend'
  }else{
    Category_50[i] <- Categories[which(c(NoTrend, ColdEdgeRamp, WarmEdgeRamp, AbundanceCentre, AbundanceEdge))]
  }
}

# Tests and summaries for 25% thresholds ----

# Proportions of each grouping.
round(table(Category_25) / length(Category_25) * 100, 1)

# Number of species in in group
table(Category_25)

# chisq.test without additional group of no trend. 
chisq.test(x = c(table(Category_25)))
chisq.test(x = c(table(Category_25[Category_25 != 'No trend'])))

# What are the abundance reductions at the edges of species ranges? 
# Estimate abundances at the edges of species ranges as a scaled percentage of maximum from model.

# OVERALL:
round(mean(Abundance_percent_cool), 3)*100
round(mean(Abundance_percent_warm), 3)*100

# ABUNDANCE CENTRE:
round(mean(Abundance_percent_cool[Category_25=='Abundance centre']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre']), 3)*100

# COOL RAMPED:
round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge']), 3)*100

# WARM RAMPED:
round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge']), 3)*100

# Number with peak at extreme of distribution
sum(Abundance_percent_cool == 1)
sum(Abundance_percent_warm == 1)


# How do these patterns vary between thermal guilds? 
Category_25_Temperate <- Category_25[which(ThermalNicheData_New %>% filter(SpeciesName %in% HighConfidenceSpecies) %>% .$ThermalGuild == 'temperate')]
Category_25_Tropical <- Category_25[which(ThermalNicheData_New %>% filter(SpeciesName %in% HighConfidenceSpecies) %>% .$ThermalGuild == 'tropical')]

table(Category_25_Temperate) / length(Category_25_Temperate)
table(Category_25_Tropical) / length(Category_25_Tropical) 

chisq.test(x = c(table(Category_25_Temperate)))
chisq.test(x = c(table(Category_25_Tropical)))

round(mean(Abundance_percent_cool[TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Abundance centre' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Abundance centre' & TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre' & TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge' & TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge' & TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge'& TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge'& TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge'& TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge'& TG == 'temperate']), 3)*100

sum(Abundance_percent_warm == 1) / length(Abundance_percent_warm)
sum(Abundance_percent_warm[TG == 'temperate'] == 1) / length(Abundance_percent_warm)
sum(Abundance_percent_warm[TG == 'tropical'] == 1) / length(Abundance_percent_warm)

sum(Abundance_percent_cool == 1) / length(Abundance_percent_cool)
sum(Abundance_percent_cool[TG == 'temperate'] == 1) / length(Abundance_percent_cool)
sum(Abundance_percent_cool[TG == 'tropical'] == 1) / length(Abundance_percent_cool)

# Tests and summaries for 50% thresholds ----

# Proportions of each grouping.
round(table(Category_50) / length(Category_50) * 100, 1)

# Number of species in in group
table(Category_50)

# chisq.test without additional group of no trend. 
chisq.test(x = c(table(Category_50)))
chisq.test(x = c(table(Category_50[Category_50 != 'No trend'])))

# ----------------------------------------------------------------------------

# ANALYSIS OF THERMAL NICHE SHAPE FOR EACH SPECIES. SKEW VS. TOPT ---- 
# Define species habitat associations ----
SiteCovariates <- read.csv(file = '/Users/caw1g15/Documents/GitHub/RLSThermalResilience/data/Site covariatesV2.csv')
SiteCovariates$HumanPop <- log10(SiteCovariates$HumanPop + 1)
SiteCovariates2 <- left_join(SiteCovariates, RLS_All %>% select(SiteCode, MeanSiteSST_NOAA) %>% unique()) %>%  
  group_by(SiteCode) %>% 
  nest() %>% 
  mutate(MeanSiteSST_NOAA_Site = purrr::map(data, ~mean(.$MeanSiteSST_NOAA, na.rm = T)),
         #SSTmean_Site = purrr::map(data, ~mean(.$SSTmean, na.rm = T)),
         #SSTrange_Site = purrr::map(data, ~mean(.$SSTrange, na.rm = T)),
         #HumanPop_Site = purrr::map(data, ~mean(.$HumanPop, na.rm = T)),
         #Exposure_Site = purrr::map(data, ~mean(.$Exposure, na.rm = T)), 
         #PredatorB_Site = purrr::map(data, ~mean(.$PredatorB, na.rm = T)), 
         LiveCoralCover_Site = purrr::map(data, ~mean(.$LiveCoralCover, na.rm = T)), 
         #ComplexCoralCover_Site = purrr::map(data, ~mean(.$ComplexCoralCover, na.rm = T)), 
         AlgalCover_Site = purrr::map(data, ~mean(.$AlgalCover, na.rm = T))) %>% 
  unnest(MeanSiteSST_NOAA_Site, 
         #SSTmean_Site, 
         #SSTrange_Site, 
         #HumanPop_Site, 
         #Exposure_Site, 
         #PredatorB_Site, 
         LiveCoralCover_Site, 
         #ComplexCoralCover_Site, 
         AlgalCover_Site) %>% 
  select(-data)

# Join RLS_All with the fine scale habitat data from RSS 
RLS_All_WithCovariates <- left_join(RLS_All, SiteCovariates2)

# Summaries covariates per species 
Species_SitePreferences <- RLS_All_WithCovariates %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(#SSTmean_Spp = purrr::map(data, ~mean(.$SSTmean_Site, na.rm = T)),
         #SSTrange_Spp = purrr::map(data, ~mean(.$SSTrange_Site, na.rm = T)),
         #HumanPop_Spp = purrr::map(data, ~mean(.$HumanPop_Site, na.rm = T)),
         #Exposure_Spp = purrr::map(data, ~mean(.$Exposure_Site, na.rm = T)), 
         #PredatorB_Spp = purrr::map(data, ~mean(.$PredatorB_Site, na.rm = T)), 
         LiveCoralCover_Spp = purrr::map(data, ~mean(.$LiveCoralCover_Site, na.rm = T)), 
         #ComplexCoralCover_Spp = purrr::map(data, ~mean(.$ComplexCoralCover_Site, na.rm = T)), 
         AlgalCover_Spp = purrr::map(data, ~mean(.$AlgalCover_Site, na.rm = T))) %>% 
  unnest(#SSTmean_Spp, 
         #SSTrange_Spp, 
         #HumanPop_Spp, 
         #Exposure_Spp, 
         #PredatorB_Spp, 
         LiveCoralCover_Spp, 
         #ComplexCoralCover_Spp, 
         AlgalCover_Spp) %>% 
  select(-data)

ThermalNicheData_New <- left_join(ThermalNicheData_New, Species_SitePreferences)




# SOM figure: Plot correlations between species habitat associations and thermal traits ---- 
pdf(file = 'figures_extra/SOM_Habitat associations vs. Topt.pdf', width = 4, height = 6.5)
gridExtra::grid.arrange(
  ggplot() + 
    geom_point(data = ThermalNicheData_New , aes(x = Topt, y = LiveCoralCover_Spp)) +
    stat_smooth(data = ThermalNicheData_New , aes(x = Topt, y = LiveCoralCover_Spp), method = lm , colour = 'gray80', se = F) + 
    theme_classic() + theme(text = element_text(size = 15)) + 
    ylab('Species mean live coral cover') + 
    xlab('Species thermal optima'),
  
  ggplot() + 
    geom_point(data = ThermalNicheData_New , aes(x = Topt, y = AlgalCover_Spp)) + 
    stat_smooth(data = ThermalNicheData_New , aes(x = Topt, y = AlgalCover_Spp), method = lm , colour = 'gray80', se = F) + 
    theme_classic() + theme(text = element_text(size = 15)) + 
    ylab('Species mean algae cover') + 
    xlab('Species thermal optima')
)
dev.off()
# SOM figure: Plot correlations between site temperature and site habitat ----
# What are the habitat differences across all sites at a larger scale? 
pdf(file = 'figures_extra/SOM_GlobalHabitat_ThermalGuild.pdf', width = 5, height = 4)
ggplot(data  = SiteCovariates2[which(SiteCovariates2$AlgalCover_Site<101),], aes(x = MeanSiteSST_NOAA_Site)) + 
  geom_point(aes(y = LiveCoralCover_Site, colour = "Coral"), alpha = 0.3) + 
  geom_point(aes(y = AlgalCover_Site,  colour = "Algae"), alpha = 0.3) + 
  stat_smooth(aes(y = AlgalCover_Site, colour = "Algae"), method = glm, se = F) +
  stat_smooth(aes(y = LiveCoralCover_Site, colour = "Coral"), method = glm, se = F) + 
  ylab('Live coral cover') + 
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Algae cover"), limits = c(0,100)) + 
  scale_colour_manual(values = c("darkblue", "dark orange")) + 
  labs(colour = "Habitat cover")  + 
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Site sea surface temperature')
dev.off()
# Set up data for further analyses below ----

# Filter to high confidence
ThermalNicheData_New_conf3 <- ThermalNicheData_New %>% filter(ConfidenceCombined == 3)

# Combine in taxonomic information
SpeciesTaxonomy <- read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv') %>% select(SPECIES_NAME, GENUS, FAMILY, ORDER)
SpeciesTaxonomy <- SpeciesTaxonomy %>% select(SPECIES_NAME, GENUS, FAMILY, ORDER) %>% unique()
SpeciesTaxonomy <- SpeciesTaxonomy %>% plyr::rename(., replace = c('SPECIES_NAME' = 'SpeciesName', 'GENUS' = 'Genus', 'FAMILY' = 'Family', 'ORDER' = 'Order'))
ThermalNicheData_New_conf3 <- left_join(ThermalNicheData_New_conf3, SpeciesTaxonomy)

# MAIN ANALYSIS 1. Model t-skew with global model (across all species conf = 3) ----

# Global models for algal cover. 
GlobalModel_algae <- lmer(T_Skew ~  
                            Topt * ThermalGuild + 
                            AlgalCover_Spp * ThermalGuild +
                            (1|Order / Family / Genus), data = ThermalNicheData_New_conf3)

GlobalModel_algae2 <- lmer(T_Skew ~  
                             Topt * ThermalGuild + 
                             AlgalCover_Spp + ThermalGuild +
                             (1|Order / Family / Genus), data = ThermalNicheData_New_conf3)

anova(GlobalModel_algae, GlobalModel_algae2, test = 'LRT') # No interaction between  algae association and Topt
#Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
#GlobalModel_algae2  9 260.46 286.74 -121.23   242.46                        
#GlobalModel_algae  10 262.46 291.66 -121.23   242.46 1e-04      1     0.9932


GlobalModel_algae3 <- lmer(T_Skew ~  
                           Topt + ThermalGuild + 
                           AlgalCover_Spp + ThermalGuild +
                           (1|Order / Family / Genus), data = ThermalNicheData_New_conf3)

anova(GlobalModel_algae2, GlobalModel_algae3, test="LRT") # Significant interaction between thermal niche Topt. 
#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#GlobalModel_algae3  8 291.11 314.47 -137.56   275.11                             
#GlobalModel_algae2  9 260.46 286.74 -121.23   242.46 32.648      1  1.104e-08 ***

GlobalModel_algae4 <- lmer(T_Skew ~  
                             Topt * ThermalGuild + 
                             (1|Order / Family / Genus), data = ThermalNicheData_New_conf3)
anova(GlobalModel_algae2, GlobalModel_algae4, test="LRT") # Significant effect of algae (similar across realms) 
#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
#GlobalModel_algae4  8 267.56 290.92 -125.78   251.56                            
#GlobalModel_algae2  9 260.46 286.74 -121.23   242.46 9.0993      1   0.002557 **

GlobalModel_algae5 <- lmer(T_Skew ~  
                             ThermalGuild + 
                             AlgalCover_Spp +
                             (1|Order / Family / Genus), data = ThermalNicheData_New_conf3)
anova(GlobalModel_algae2, GlobalModel_algae5, test="LRT") # Significant effect of Topt 
#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#GlobalModel_algae5  7 466.15 486.59 -226.08   452.15                             
#GlobalModel_algae2  9 260.46 286.74 -121.23   242.46 209.69      2  < 2.2e-16 ***
  
# Remove taxonomy 
GlobalModel_algae2_ML_nophylo <- lm(T_Skew ~  
                                      Topt * ThermalGuild + 
                                      AlgalCover_Spp + ThermalGuild, data = ThermalNicheData_New_conf3, method = 'ML') # 
GlobalModel_algae2_ML <- lmer(T_Skew ~  
                                Topt * ThermalGuild + 
                                AlgalCover_Spp + ThermalGuild +
                                (1|Order / Family / Genus), data = ThermalNicheData_New_conf3, method = 'ML') # 
# Taxonomic terms are not essential
AICc(GlobalModel_algae2_ML_nophylo, GlobalModel_algae2_ML)



# Test with drop1
drop1(GlobalModel_algae2) # All terms are important. 
#                  Df    AIC
#<none>               260.46
#AlgalCover_Spp     1 267.56
#Topt:ThermalGuild  1 291.11


anova(GlobalModel_algae2, test = 'LRT') 
#Analysis of Variance Table
#                  Df Sum Sq Mean Sq F value
#Topt               1 75.908  75.908 220.622
#ThermalGuild       1 91.663  91.663 266.415
#AlgalCover_Spp     1  5.872   5.872  17.066
#Topt:ThermalGuild  1 12.262  12.262  35.638

summary(GlobalModel_algae2)
#                          Estimate Std. Error t value
#(Intercept)               12.19891    1.40282   8.696
#Topt                      -0.51780    0.05033 -10.288
#ThermalGuildtropical      12.86423    1.56603   8.215
#AlgalCover_Spp            -0.04261    0.01416  -3.009
#Topt:ThermalGuildtropical -0.39104    0.06550  -5.970

r.squaredGLMM(GlobalModel_algae2)
# R2m       R2c 
# 0.7998033 0.8070411 

# Check residuals
plot(GlobalModel_algae2)
plot(fitted(GlobalModel_algae2), resid(GlobalModel_algae2, type = 'pearson'))
qqplot(y = resid(GlobalModel_algae2), x = rnorm(1000)) # Errors are approximately normally distributed


# Predict from 'global' model ----

# Create dataframe to predict over. 
ModelGlobal_Pred <- plyr::ddply(ThermalNicheData_New_conf3, 
                                     .(ThermalGuild), 
                                     plyr::summarize,
                                     Topt = seq(min(Topt), max(Topt), length = 100))

ModelGlobal_Pred$AlgalCover_Spp <- mean(ThermalNicheData_New_conf3 %>%  .$AlgalCover_Spp)

# Predict from model coefficients
X  <- model.matrix(formula(~Topt * ThermalGuild + AlgalCover_Spp + ThermalGuild), data = ModelGlobal_Pred)
ModelGlobal_Pred$y <- as.numeric(X %*% fixef(GlobalModel_algae2))
ModelGlobal_Pred$SE <- as.numeric(sqrt(diag(X %*% vcov(GlobalModel_algae2) %*% t(X))))
ModelGlobal_Pred$SeUp  <- ModelGlobal_Pred$y + 1.96 
ModelGlobal_Pred$SeLo  <- ModelGlobal_Pred$y - 1.96

# Create partial residuals. 
Globalmodel_PartialResiduals <- remef(GlobalModel_algae2, fix = c("AlgalCover_Spp"), ran = list("Genus:(Family:Order)" = 1), grouping = T)
Globalmodel_PartialResiduals = Globalmodel_PartialResiduals + (mean(ThermalNicheData_New_conf3$T_Skew) - mean(Globalmodel_PartialResiduals))


# MAIN RESULTS 'FIGURE 2' ---- 
TND_temperate <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'temperate')
TND_tropical  <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'tropical')

# Select example species near to fitted line
Skew_TempSpp_min <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[3]
Skew_TempSpp_max <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt),] %>% .$SpeciesName %>% .[4]
Skew_TropSpp_min <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'tropical') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[4]
Skew_TropSpp_max <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'tropical') %>% filter(T_Skew == max(.$T_Skew))  %>% .$SpeciesName

pdf(file = 'figure_final/Topt_GlobalModel_V2.pdf', width = 5, height = 5, useDingbats = F)
ggplot() + 
  geom_ribbon(data = ModelGlobal_Pred %>% filter(ThermalGuild == 'temperate'), aes(x = Topt, ymax = SeUp, ymin = SeLo), alpha = 0.1, fill = 'darkblue') + 
  geom_ribbon(data = ModelGlobal_Pred %>% filter(ThermalGuild == 'tropical'), aes(x = Topt,  ymax = SeUp, ymin = SeLo), alpha = 0.1, fill = 'dark orange') + 
  
  
  geom_point(data = TND_temperate %>% filter(ThermalGuild == 'temperate'),  
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$ThermalGuild == 'temperate')]), colour = 'darkblue', alpha = 1, size = 3) + 
  geom_point(data = TND_tropical %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$ThermalGuild == 'tropical')]), colour = 'dark orange', alpha = 1, size = 3) + 
  
  geom_line(data = ModelGlobal_Pred %>% filter(ThermalGuild == 'temperate'), aes(y = y, x = Topt), colour = 'black') + 
  geom_line(data = ModelGlobal_Pred %>% filter(ThermalGuild == 'tropical'), aes(y = y, x = Topt), colour = 'black') + 
  
  geom_point(data = TND_temperate %>% filter(SpeciesName == Skew_TempSpp_min),  
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$SpeciesName == Skew_TempSpp_min)]), col = 'darkblue' , fill = 'white', pch = 21, size = 5) + 
  geom_point(data = TND_tropical %>% filter(SpeciesName == Skew_TropSpp_min), 
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$SpeciesName == Skew_TropSpp_min)]),col = 'dark orange', fill = 'white', pch = 21, size = 5) + 
  
  geom_point(data = TND_temperate %>% filter(SpeciesName == Skew_TempSpp_max),  
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$SpeciesName == Skew_TempSpp_max)]), col = 'darkblue' , fill = 'gray75', pch = 21, size = 5) + 
  geom_point(data = TND_tropical %>% filter(SpeciesName == Skew_TropSpp_max), 
             aes(x = Topt, y = Globalmodel_PartialResiduals[which(ThermalNicheData_New_conf3$SpeciesName == Skew_TropSpp_max)]),col = 'dark orange', fill = 'gray75', pch = 21, size = 5) + 
  
  theme_classic() + theme(text = element_text(size = 15), aspect.ratio = 0.6) + 
  xlab(expression(T["opt"])) + 
  ylab(expression(T["skew"])) + 
  scale_x_continuous(breaks = seq(round(min(ModelGlobal_Pred$Topt)), round(max(ModelGlobal_Pred$Topt)), by = 2))
dev.off()

# ----------------------------------------------------------------------------

# SUPPORTING ANALYSIS: ANALYSIS OF THERMAL NICHE SHAPE FOR EACH SPECIES. SKEW VS. TOPT ----
# Model t-skew with global model (across all species conf = 3) TROPICAL ----

# Select tropical species only 
ThermalNicheData_New_conf3_trop <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'tropical')

# Tropical models for coral cover. 
TropicalModel_coral <- lmer(T_Skew ~  
                              Topt +
                              LiveCoralCover_Spp +
                              (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_trop) # 

TropicalModel_coral2 <- lmer(T_Skew ~  
                               LiveCoralCover_Spp +
                               (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_trop) # 

anova(TropicalModel_coral, TropicalModel_coral2) # Siginificant temperature effect. 
#                     Df    AIC    BIC   logLik deviance Chisq Chi Df Pr(>Chisq)    
#TropicalModel_coral2  6 284.44 298.50 -136.219   272.44                            
#TropicalModel_coral   7 138.54 154.95  -62.271   124.54 147.9      1  < 2.2e-16 ***

TropicalModel_coral3 <- lmer(T_Skew ~  
                               Topt + 
                               (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_trop) # 

anova(TropicalModel_coral, TropicalModel_coral3) # Significant coral effect.  
#                     Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#TropicalModel_coral3  6 141.83 155.90 -64.916   129.83                           
#TropicalModel_coral   7 138.54 154.95 -62.271   124.54 5.2906      1    0.02144 *

# Remove phylogeny 
TropicalModel_coral_ML_nophylo <- lm(T_Skew ~  
                                       Topt +
                                       LiveCoralCover_Spp, data = ThermalNicheData_New_conf3_trop, REML = 'ML') # 
TropicalModel_coral_ML <- lmer(T_Skew ~  
                                 Topt +
                                 LiveCoralCover_Spp +
                                 (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_trop, method = 'ML') # 
AICc(TropicalModel_coral_ML_nophylo, TropicalModel_coral_ML)

# Define final model
TropicalModel_coral

# Keep all remaining terms. 
drop1(TropicalModel_coral)
#                   Df    AIC
#<none>                138.54
#Topt                1 284.44
#LiveCoralCover_Spp  1 141.83

summary(TropicalModel_coral)
#Fixed effects:
#                   Estimate Std. Error t value
#(Intercept)        22.54465    1.11280  20.259
#Topt               -0.91133    0.04320 -21.093
#LiveCoralCover_Spp  0.05915    0.02326   2.543

r.squaredGLMM(TropicalModel_coral)
#R2m       R2c 
#0.8440504 0.8748064 

# Check residuals
plot(TropicalModel_coral)
plot(fitted(TropicalModel_coral), resid(TropicalModel_coral, type = 'pearson'))
qqplot(y = resid(TropicalModel_coral), x = rnorm(1000)) # Errors are normally distributed

# Do we get the same effect when fitted to < 26°C
TropicalModel_coral_Final_UpperTruncated      <- lmer(formula(TropicalModel_coral), data = ThermalNicheData_New_conf3_trop %>% filter(Topt < 26))
# Fixed Effects:
#   (Intercept)                Topt  LiveCoralCover_Spp  
# 28.41255            -1.17050             0.07255  

# Model t-skew with global model (across all species conf = 3) TEMPERATE ----

# Select temperate species only 
ThermalNicheData_New_conf3_temp <- ThermalNicheData_New_conf3 %>% filter(ThermalGuild == 'temperate')

# Temperate models for algae cover. 
TemperateModel_algae <- lmer(T_Skew ~  
                               Topt +
                               AlgalCover_Spp +
                               (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_temp) # 

TemperateModel_algae2 <- lmer(T_Skew ~  
                                AlgalCover_Spp +
                                (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_temp) # 

anova(TemperateModel_algae, TemperateModel_algae2) # Siginificant temperature effect. 
#                      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#TemperateModel_algae2  6 180.72 193.29 -84.361   168.72                             
#TemperateModel_algae   7 126.92 141.58 -56.460   112.92 55.803      1  8.013e-14 ***
  
TemperateModel_algae3 <- lmer(T_Skew ~  
                                Topt + 
                                (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_temp) # 

anova(TemperateModel_algae, TemperateModel_algae3) # Significant algae effect  
#                      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#TemperateModel_algae3  6 130.65 143.22 -59.327   118.65                           
#TemperateModel_algae   7 126.92 141.58 -56.460   112.92 5.7339      1    0.01664 *
  
# Remove phylogeny 
TemperateModel_algae_ML_nophylo <- lm(T_Skew ~  
                                        Topt +
                                        AlgalCover_Spp, data = ThermalNicheData_New_conf3_temp, REML = F) # 
TemperateModel_algae_ML <- lmer(T_Skew ~  
                                  Topt +
                                  AlgalCover_Spp +
                                  (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_temp, REML = F) # 
AICc(TemperateModel_algae_ML_nophylo, TemperateModel_algae_ML)
#                                df     AICc
#TemperateModel_algae_ML_nophylo  4 121.9886
#TemperateModel_algae_ML          7 129.0734

# Define final model
TemperateModel_algae

summary(TemperateModel_algae)
#                Estimate Std. Error t value
# (Intercept)    12.24605    1.63936   7.470
# Topt           -0.52579    0.05443  -9.660
# AlgalCover_Spp -0.04288    0.01765  -2.429

r.squaredGLMM(TemperateModel_algae)
# R2m       R2c 
# 0.6013954 0.6439650 

# Check residuals
plot(TemperateModel_algae)
plot(fitted(TemperateModel_algae), resid(TemperateModel_algae, type = 'pearson'))
qqplot(y = resid(TemperateModel_algae), x = rnorm(1000)) # Errors are normally distributed


# Model t-skew with global model (across all species conf = 4) < 26°C Topt ----

# Subset data for < 26°C global data 
ThermalNicheData_New_conf3_26 <- ThermalNicheData_New_conf3 %>% filter(Topt < 26)

# Global models for algal cover. 
GlobalModel_algae_conservative <- lmer(T_Skew ~  
                                         Topt * ThermalGuild + 
                                         AlgalCover_Spp * ThermalGuild +
                                         (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_26)

GlobalModel_algae_conservative2 <- lmer(T_Skew ~  
                                          Topt * ThermalGuild + 
                                          AlgalCover_Spp + ThermalGuild +
                                          (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_26)

anova(GlobalModel_algae_conservative, GlobalModel_algae_conservative2) # No interaction between algae association and Topt
#                                Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#GlobalModel_algae_conservative2  9 163.09 184.19 -72.546   145.09                         
#GlobalModel_algae_conservative  10 164.24 187.68 -72.121   144.24 0.8507      1     0.3564


GlobalModel_algae_conservative3 <- lmer(T_Skew ~  
                                          Topt + 
                                          AlgalCover_Spp + ThermalGuild +
                                          (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_26)

anova(GlobalModel_algae_conservative2, GlobalModel_algae_conservative3, test="LRT") # Significant interaction between T-opt and thermal guild 
#                                 Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# GlobalModel_algae_conservative3  8 166.37 185.12 -75.184   150.37                           
# GlobalModel_algae_conservative2  9 163.09 184.19 -72.546   145.09 5.2754      1    0.02163 *

GlobalModel_algae_conservative4 <- lmer(T_Skew ~  
                                          Topt*ThermalGuild + 
                                          (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_26)
anova(GlobalModel_algae_conservative2, GlobalModel_algae_conservative4, test="LRT") # Significant effect of algae (similar across realms) 
#                                Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
#GlobalModel_algae_conservative4  8 168.51 187.26 -76.257   152.51                            
#GlobalModel_algae_conservative2  9 163.09 184.19 -72.546   145.09 7.4213      1   0.006445 **

GlobalModel_algae_conservative5 <- lmer(T_Skew ~  
                                          ThermalGuild + 
                                          AlgalCover_Spp +
                                          (1|Order / Family / Genus), data = ThermalNicheData_New_conf3_26)
anova(GlobalModel_algae_conservative2, GlobalModel_algae_conservative5, test="LRT") # Significant effect of Topt 
#                                Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
#GlobalModel_algae_conservative5  7 227.47 243.88 -106.736   213.47                             
#GlobalModel_algae_conservative2  9 163.09 184.19  -72.546   145.09 68.379      2  1.418e-15 ***

drop1(GlobalModel_algae_conservative2)               # All terms are important. 
# T_Skew ~ Topt * ThermalGuild + AlgalCover_Spp + (1 | Order/Family/Genus)
# Df    AIC
# <none>               163.09
# AlgalCover_Spp     1 168.51
# Topt:ThermalGuild  1 166.37

anova(GlobalModel_algae_conservative2, test = 'LRT') 
#                  Df Sum Sq Mean Sq  F value
#Topt               1  3.239   3.239   8.4516
#ThermalGuild       1 39.641  39.641 103.4501
#AlgalCover_Spp     1  2.729   2.729   7.1210
#Topt:ThermalGuild  1  2.200   2.200   5.7426

summary(GlobalModel_algae_conservative2)
#                          Estimate Std. Error t value
#(Intercept)               12.70974    1.62998   7.797
#Topt                      -0.53167    0.05515  -9.641
#ThermalGuildtropical      17.50139    6.11408   2.862
#AlgalCover_Spp            -0.04913    0.01728  -2.843
#Topt:ThermalGuildtropical -0.58931    0.24592  -2.396

r.squaredGLMM(GlobalModel_algae_conservative2)
# R2m       R2c 
# 0.6103523 0.6514894 

# Check residuals
plot(GlobalModel_algae_conservative2)
plot(fitted(GlobalModel_algae_conservative2), resid(GlobalModel_algae_conservative2, type = 'pearson'))
qqplot(y = resid(GlobalModel_algae_conservative2), x = rnorm(1000)) # Errors are normally distributed

  
  
  

# ----------------------------------------------------------------------------

# TABLE FOR SOM: ----
stargazer(TropicalModel_coral, TemperateModel_algae, GlobalModel_algae2, GlobalModel_algae_conservative2,
          type = 'html', out = 'figure_final/Tskew analysis outputs.htm', 
          dep.var.labels=c("T-skew"), 
          covariate.labels = c('T-opt', 
                               'Coral association', 
                               'Thermal guild (tropical)', 
                               'Algae association', 
                               'T-opt : Thermal guild (tropical)', 
                               'Intercept'), 
          column.labels = c('Tropical', 'Temperate', 'Global', 'Global < 26°C'))
# ----------------------------------------------------------------------------





