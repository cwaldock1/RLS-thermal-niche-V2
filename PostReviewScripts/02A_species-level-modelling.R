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
library(summarytools)

# Fitting PCAs 
library(pcaMethods)
library(factoextra)

# Fitting quantile gams. 
library(qgam)

# Fitting occupancy models. 
library(glmmTMB)

# Doing parallel stuff
library(doParallel)

# Load data from 01_organising-data.R ----

# Old file upload with more objects
# load(file = 'data_derived/01_organising-data_SAVE-IMAGE.RData') 
# rm(MaximumAbundance, RLS_HighConfidence_Species, RLS_LowConfidence_Species, RLS11, RLS12, RLS12_DT, RLS13_DT, RLS14, RLS14_DT, RLS15, RLS_Site_Covariates)

# Saved only RLS19 which is necessary for this script. 
#save.image(file = 'data_derived/01-organsiing-data05072018.RData')
#load(file = 'data_derived/01-organsiing-data05072018.RData')
load(file = 'data_derived/01-organsiing-data24092018.RData')

#view(dfSummary(RLS_19))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------





# CREATE COVARIATE PCA SCORES ----
# Read in covariate data and produce pairs plot and explore structure and correlations ----
Covariates <- as_data_frame(read.csv(file =  'data_raw/RLS_Site_Covariates_V2_2018-09-26.csv'))
view(dfSummary(Covariates))

# Filter to only relevant sites for modelling. 
Covariates <- Covariates %>% filter(SiteCode %in% unique(RLS_19$SiteCode))

# Transform those which need it. 
Covariates$HumPop50km2 <- log10(Covariates$HumPop50km2 + 1)
Covariates$ReefAreaIn15km <- log10(Covariates$ReefAreaIn15km + 1)
Covariates$npp_mean <- log(Covariates$npp_mean + 1)
Covariates$npp_max <- log(Covariates$npp_max + 1)
Covariates$Silicate.Mean <- log(Covariates$Silicate.Mean + 1)
Covariates$Iron.Mean <- log10(Covariates$Iron.Mean*1000 + 1)
Covariates$Current.Velocity.Mean <- log10(Covariates$Current.Velocity.Mean + 1)
Covariates$Current.Velocity.Lt.max <- log10(Covariates$Current.Velocity.Lt.max + 1)

# Salinity mean vs. max is highly correalted 
# NPP mean vs. max is highly correlated 
# Oxygen mean and max are highly correlated
# Nitrate shows high number of extreme points that will swamp variation in PCA. It is highly correlated with phosphate anyway. 
# If highly correlated values are input then this will swamp the PCA without and combined information getting through. 

Covariates_V2 <- Covariates %>% dplyr::select(SiteCode, HumPop50km2, ReefAreaIn15km, 
                                       npp_mean, Silicate.Mean, Salinity.Mean, 
                                       Phosphate.Mean, Iron.Mean,
                                       Dissolved.oxygen.Mean, Current.Velocity.Mean, Current.Velocity.Lt.max, 
                                       pH)

pdf('figure_final/CovariatePairsPlot.pdf', width = 15, height = 15)
psych::pairs.panels(Covariates_V2[, -1])
dev.off()

# Perform PCA analysis ----

PCA_data.input <- as.matrix(Covariates_V2[, -1])

# Compute PCA
#PCA_results <- prcomp(na.omit(PCA_data.input), scale = TRUE)
PCA_results <- prcomp(na.omit(PCA_data.input), scale = TRUE) # Just testing a subset. 

# Visualise variance explained by each component
pdf('figure_final/PCA-variance.pdf', width = 4, height = 3)
fviz_eig(PCA_results); dev.off()


# Plotting PCA axis ----
# How much is explained by first 3 components? 
round(summary(PCA_results)$importance, 2)# First 3 components explain 83% variation. 

# Extract results for variables
res.var <- get_pca_var(PCA_results)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
round(res.var$contrib, 2)
res.var$cos2           # Quality of representation 

VariablePCA_Coords <- cbind(melt(res.var$coord), melt(res.var$contrib)[,3])
levels(VariablePCA_Coords$Var1) <- gsub('.Mean','', VariablePCA_Coords$Var1)
names(VariablePCA_Coords)[c(3,4)] <- c('coord', 'contrib')

res.ind <- get_pca_ind(PCA_results)
res.ind$coord          # Coordinates ####################### this is what I want to model with. 
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

res.var1 <- res.var$coord %>% data.frame
res.var1$vars <- rownames(res.var1)

# Produce boxplots of coordinates from indepedent axis. 
ggplot(VariablePCA_Coords) + 
  geom_bar(aes(x = Var1, y = coord, fill = contrib), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 1') + xlab(NULL) +
  facet_wrap(~Var2, scale = 'free_y')


BarAxis1 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.1, fill = Dim.1), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 1') + xlab(NULL) 

BarAxis2 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.2, fill = Dim.2), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 2') + xlab(NULL)

BarAxis3 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.3, fill = Dim.3), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 3') + xlab(NULL)

pdf('figure_final/PCAbar-directions.pdf', width = 5, height = 15)
gridExtra::grid.arrange(BarAxis1, BarAxis2, BarAxis3) 
dev.off() 


# This shows the correlation between variables. Positively correlated variables point to the same sides of the plots. 
Arrow1 <- fviz_pca_var(PCA_results, axes = c(1, 2), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 1-2') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Arrow2 <- fviz_pca_var(PCA_results, axes = c(2, 3), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 2-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Arrow3 <- fviz_pca_var(PCA_results, axes = c(1, 3), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 1-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds1 <- fviz_pca_ind(PCA_results, axes = c(1, 2), col.ind = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 1-2') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds2 <- fviz_pca_ind(PCA_results, axes = c(2, 3), col.ind  = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 2-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds3 <- fviz_pca_ind(PCA_results, axes = c(1, 3), col.ind  = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 1-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)


pdf('figure_final/PCA-multiplot.pdf', width = 15, height = 15)
gridExtra::grid.arrange(
  BarAxis1 , Arrow1  , Arrow3, 
  Inds1   , BarAxis2, Arrow2, 
  Inds3   , Inds2  , BarAxis3, 
  nrow = 3, ncol = 3)
dev.off()

# Combind in site level coordinates
#res.ind <- get_pca_ind(PCA_results)
#Covariates_V2$PCA_1 <- res.ind$coord[,1]
#Covariates_V2$PCA_2 <- res.ind$coord[,2]
#Covariates_V2$PCA_3 <- res.ind$coord[,3]
#Covariates_V2$PCA_4 <- res.ind$coord[,4]


# Relationship between axis and temperature ----

# Match with temperature data
Covariates_WithTemp <- RLS_19 %>% dplyr::select(SiteCode, MeanTemp_CoralWatch) %>% unique() %>% left_join(., Covariates_V2 %>% dplyr::select(SiteCode, PCA_1, PCA_2, PCA_3, PCA_4))

pdf(file = 'figure_final/PairsWithTemp.pdf', width = 8, height = 8)
psych::pairs.panels(Covariates_WithTemp[-1], scale = T)
dev.off()

# Test PCA on subset of species ----
LabroidesSites <- RLS_19 %>%  filter(SpeciesName == unique(.$SpeciesName)[1]) %>% dplyr::select(SiteCode) %>% unique() %>% .$SiteCode %>% as.character()
                              #filter(SpeciesName == 'Labroides dimidiatus') %>% dplyr::select(SiteCode) %>% unique() %>% .$SiteCode %>% as.character()
Covariates_V2_Labroides <- Covariates_V2 %>% filter(SiteCode %in% LabroidesSites) %>% dplyr::select(-PCA_1, -PCA_2, -PCA_3, -PCA_4)

PCA_Labroides <- prcomp(na.omit(as.matrix(Covariates_V2_Labroides[,-1])), scale = TRUE) # Just testing a subset. 
fviz_eig(PCA_Labroides)
round(summary(PCA_Labroides)$importance, 2)# First 3 components explain 83% variation. 

res.var_Lab <- get_pca_var(PCA_Labroides)
res.var_Lab <- res.var_Lab$coord %>% data.frame; 
res.var_Lab$vars <- rownames(res.var_Lab)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.1, fill = Dim.1), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.2, fill = Dim.2), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.3, fill = Dim.3), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.4, fill = Dim.4), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.5, fill = Dim.5), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.6, fill = Dim.6), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.7, fill = Dim.7), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)

res.ind_Lab <- get_pca_ind(PCA_Labroides)
Covariates_V2_Labroides$PCA_1 <- res.ind_Lab$coord[,1]
Covariates_V2_Labroides$PCA_2 <- res.ind_Lab$coord[,2]
Covariates_V2_Labroides$PCA_3 <- res.ind_Lab$coord[,3]
Covariates_V2_Labroides$PCA_4 <- res.ind_Lab$coord[,4]

Covariates_WithTemp_SPP <- RLS_19 %>% filter(SiteCode %in% LabroidesSites) %>% dplyr::select(SiteCode, MeanTemp_CoralWatch) %>% unique() %>% left_join(., Covariates_V2_Labroides %>% dplyr::select(SiteCode, PCA_1, PCA_2, PCA_3, PCA_4))

psych::pairs.panels(Covariates_WithTemp_SPP[-1], scale = T)

# Estimate PCA scores for all species within sites ----

# Example input data 
Species_Sites_Temp <- RLS_19 %>% filter(SpeciesName == unique(.$SpeciesName)[1]) %>% dplyr::select(SpeciesName, SiteCode, MeanTemp_CoralWatch) %>% unique()
Species_Sites_Temp <- TestSpecies %>% dplyr::select(SpeciesName, SiteCode, MeanTemp_CoralWatch) %>% unique()

SpeciesPCAScore <- function(Species_Sites_Temp,          # Input a dataframe of temperature by sitecodes
                            covariates = Covariates_V2){
  
  # Convert species subset to a matrix 
  Covariates_Species <- covariates %>% filter(SiteCode %in% Species_Sites_Temp$SiteCode)
  
  if(sum(apply(Covariates_Species[,-1], 2, sd) == 0) != 0){
  #If a column is singluar it breaks the pc. Remove all singluar columns
  Covariates_Species <- Covariates_Species[, -which(apply(Covariates_Species, 2, sd) == 0)]
  }else{NULL}
  PCA_Species <- prcomp(na.omit(as.matrix(Covariates_Species[,-1])), scale = TRUE) # Just testing a subset. 
  
  # Join up temperature with site codes to ensure consistent sitecode order in single object used for PCA
  Species_Sites_Temp <- left_join(Covariates_Species %>% dplyr::select(SiteCode), Species_Sites_Temp, by = 'SiteCode')
  
  # Create output of PCA for variables to decide how many. 
  PCA_Variances <- round(summary(PCA_Species)$importance, 2)
  PCAs <- colnames( PCA_Variances[, -which(PCA_Variances[2,] < 0.1) ] )
  
  # Create columns in temperature by site code based on PCAs. 
  IndividualPCA <- get_pca_ind(PCA_Species)
  
  ToFill <- matrix(NA, nrow(Covariates_Species), length(PCAs)+3)
  ToFill[,1] <- Species_Sites_Temp$SpeciesName
  ToFill[,2] <- Species_Sites_Temp$SiteCode
  ToFill[,3] <- Species_Sites_Temp$MeanTemp_CoralWatch
  for(i in 1:length(PCAs)){
    ToFill[,i+3] <- IndividualPCA$coord[,i]
  }
  
  # Create dataframe to link with species covariates. 
  PCA_Covarites <- data.frame(ToFill, stringsAsFactors = F)
  names(PCA_Covarites) <- c('SpeciesName', 'SiteCode','MeanTemp_CoralWatch', PCAs) 
  PCA_Covarites[,-c(1:2)] <- apply(PCA_Covarites[,-c(1:2)], 2, as.numeric)
  # psych::pairs.panels(PCA_Covarites[,-1])
  # str(PCA_Covarites)
  
  # Determine correlation between PCAs and temperature 
  TempPCAcor <- matrix(NA, 1, length(PCAs) + 1)
  TempPCAcor[,1] <- unique(Species_Sites_Temp$SpeciesName)
  for(i in 1:length(PCAs)){
    TempPCAcor[,i+1] <- cor(PCA_Covarites$MeanTemp_CoralWatch, 
                            PCA_Covarites[,i+3])
  }
  TempPCAcor <- data.frame(TempPCAcor)
  names(TempPCAcor) <- c('SpeciesName', paste0(PCAs, 'Temp_cor'))
  
  # Bundle up for output
  PCA_Covarites_DF <- as_data_frame(PCA_Covarites) %>% dplyr::select(-MeanTemp_CoralWatch) %>% nest(.key = "PCA_Covarites")
  PCA_Covarites_DF$PCA_Temp_Cors <- list(TempPCAcor) #%>% nest(.key = "PCA_Temp_Cors")
  return(PCA_Covarites_DF) 
}

SpeciesPCAScore(Species_Sites_Temp = Species_Sites_Temp, 
                covariates = Covariates_V2)

# Test in dplyr loop 
PCA_Outputs <- list()
for(i in 1:length(unique(RLS_19$SpeciesName))){
  PCA_Outputs[[i]] <- tryCatch(SpeciesPCAScore(Species_Sites_Temp = RLS_19 %>% filter(SpeciesName == unique(.$SpeciesName)[i]) %>% dplyr::select(SiteCode, MeanTemp_CoralWatch, SpeciesName) %>% unique(), 
                               covariates = Covariates_V2), error = function(e) NA)
}

PCA_Fails <- which(cbind(lapply(PCA_Outputs, function(x) length(class(x)) == 1)) == T)

# Species which run PCAs
PCA_Outputs_Run <- do.call(rbind, PCA_Outputs)

# Combbine. 
PCA_Covariates_Final  <- PCA_Outputs_Run %>% unnest(PCA_Covarites)
PCA_Correlations_Final <- PCA_Outputs_Run %>% unnest(PCA_Temp_Cors)
hist(as.numeric(PCA_Correlations_Final$PC1Temp_cor))
hist(as.numeric(PCA_Correlations_Final$PC2Temp_cor))
hist(as.numeric(PCA_Correlations_Final$PC3Temp_cor))

# These are the species which produce NAs for PCAs. 
MissingSpecies <- unique(RLS_19$SpeciesName)[!unique(RLS_19$SpeciesName) %in% unique(PCA_Correlations_Final$SpeciesName)]
TestSpecies <- RLS_19 %>% filter(SpeciesName == MissingSpecies[1])


# Combine species x site PCA scores with RLS data  ----
RLS_20 <- left_join(RLS_19, PCA_Covariates_Final, by = c('SpeciesName', 'SiteCode'))
sum(is.na(RLS_20$PC1)) # These are the species' removed. 
identical(round(RLS_20$MeanTemp_CoralWatch.x) , round(RLS_20$MeanTemp_CoralWatch.y))
hist(RLS_20$MeanTemp_CoralWatch.x - RLS_20$MeanTemp_CoralWatch.y)

# saveRDS(RLS_20, file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')

# ----


# NEW MODELLING QGAMS SECTION ----
ConstrainAbsences <- function(Data){
  
  NrowPresences <- Data %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences <-  Data %>% filter(AbundanceAdult40 == 0) %>% nrow(.)
  
  if(NrowPresences <= NrowAbsences){
    Data <- rbind(Data %>% filter(AbundanceAdult40 > 0), Data %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
  }else{NULL}
  
  return(Data)}

# For a subset of species (choose by N) fit the qgam model and explore the inflight of PCA on temperature effect
Top10 <- RLS_20 %>% filter(Presence == 1) %>%  group_by(SpeciesName) %>% do(N_Sites = nrow(.)) %>% unnest(N_Sites) %>% .[order(.$N_Sites, decreasing = T),] %>% .$SpeciesName %>% .[1:10]
Bottom10 <- RLS_20 %>% filter(Presence == 1) %>%  group_by(SpeciesName) %>% do(N_Sites = nrow(.)) %>% unnest(N_Sites) %>% .[order(.$N_Sites, decreasing = F),] %>% .$SpeciesName %>% .[1:10]

Species1 <- RLS_20 %>% filter(SpeciesName == Top10[9])# %>% ConstrainAbsences(.)
Species1 <- RLS_20 %>% filter(SpeciesName == Bottom10[9])# %>% ConstrainAbsences(.)
sum(is.na(Species1$PC1)) # These are the species' removed. 

# transform, scale and centre columns for modelling
Species1$AbundanceAdult40_log           <- as.numeric(log(Species1$AbundanceAdult40 + 1))
Species1$SamplingIntensity_Scaled       <- as.numeric(scale(Species1$SamplingIntensity))
Species1$MeanTemp_CoralWatch_Scaled     <- as.numeric(scale(Species1$MeanTemp_CoralWatch))
Species1$Depth_Site                     <- as.numeric(scale(Species1$Depth_Site))
Species1$NEOLI[is.na(Species1$NEOLI)]   <- 0
Species1$PC1                            <- as.numeric(Species1$PC1)
Species1$PC2                            <- as.numeric(Species1$PC2)
Species1$PC3                            <- as.numeric(Species1$PC3)
Species1 <- Species1[,-which(apply(apply(Species1, 2, is.na), 2, sum) > 1)]

# Fit qgam model 
summary(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
        data = Species1, family = gaussian, zi = ~ 1))
GLMM_Resid <- resid(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
        data = Species1, family = gaussian, zi = ~ 1))
hist(GLMM_Resid)

Species1$GLMM_Resid <- GLMM_Resid
Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = 4), qu = 0.8, data = Species1)
Species1_GAM_V2 <- qgam(AbundanceAdult40_log ~ 
                          s(MeanTemp_CoralWatch_Scaled, k = 4) + 
                          s(Depth_Site, k = 4) + 
                          s(PC1, k = 4) + 
                          s(PC2, k = 4) + 
                          s(PC3, k = 4) + 
                          SamplingIntensity_Scaled + 
                          NEOLI
                        , qu = 0.8, data = Species1)

summary(Species1_GAM)
summary(Species1_GAM_V2)

plot(Species1_GAM, pages = 1)
plot(Species1_GAM_V2, pages = 1)

#ResidQgam <- resid(Species1_GAM)

par(mfrow = c(1,2))
plot(AbundanceAdult40_log ~ MeanTemp_CoralWatch_Scaled, data = Species1)
#plot(ResidQgam ~ MeanTemp_CoralWatch_Scaled, data = Species1)
plot(GLMM_Resid ~ MeanTemp_CoralWatch_Scaled, data = Species1)



# Function to fit quantile gams accounting for covariates ----
FitModels_27092015 <- function(Species1, k = 4, q = 0.8, NRUNS = 1, Folder = NULL){
  
  # Check before with a particular dataset - but at present there are 2 NAs which cause problems with 1 species. 
  Species1 <- Species1[!is.na(Species1$MeanTemp_CoralWatch),]
  
  # Ensure balanced design of positive abundance values and absences. 
  NrowPresences <- Species1 %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences  <- Species1 %>% filter(AbundanceAdult40 == 0) %>% nrow(.)

  # Scale and transform variables as appropriate
  Species1$AbundanceAdult40_log           <- as.numeric(log(Species1$AbundanceAdult40 + 1))
  Species1$SamplingIntensity_Scaled       <- as.numeric(scale(Species1$SamplingIntensity))
  Species1$MeanTemp_CoralWatch_Scaled     <- as.numeric(scale(Species1$MeanTemp_CoralWatch))
  Species1$Depth_Site                     <- as.numeric(scale(Species1$Depth_Site))
  Species1$NEOLI[is.na(Species1$NEOLI)]   <- 0
  
  # Remove and columsn with NA values as this malfunctions qgam. 
  Species1 <- Species1[,-which(apply(apply(Species1, 2, is.na), 2, sum) > 1)]
  
  # Convert varying number of PCs to numerics  
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  PC_Columns <- apply(PC_Columns, 2, as.numeric)
  for(i in 1:ncol(PC_Columns)){
    Species1[,colnames(Species1[,grepl('PC', colnames(Species1))])[i]] <- PC_Columns[,i]
  }
  
  # Create a model structure that can be flexible to number of PCs. 
  BaselineFormula <- 'AbundanceAdult40_log ~ poly(Depth_Site, 2) + SamplingIntensity_Scaled + NEOLI'
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  for(i in 1:ncol(PC_Columns)){
    BaselineFormula <- paste0(BaselineFormula, ' + ', 'poly(', colnames(Species1[,grepl('PC', colnames(Species1))])[i], ',2)')
  }
  BaselineFormula <- formula(BaselineFormula)
  
  # Columns for backtransformation of temperature data
  Species1$MeanTemp_CoralWatch_scaleBT  <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:scale")
  Species1$MeanTemp_CoralWatch_centerBT <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:center")
  
  # Bootstrapping procedure for estimation of Topt and its SD.
  Temps <- Species1 %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanTemp_CoralWatch_Scaled), MaxTemp = max(.$MeanTemp_CoralWatch_Scaled)) %>% unnest()
  TempData <- data.frame(MeanTemp_CoralWatch_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 1000))
  
  # Objects to assign inside bootstrap 
  Topt <- c()
  T_Gam_pvalue <- c()
  T_Gam_edf <- c()

  # If there is subsampling then need to bootstrap to get more reliable value for Topt. 
  if(NrowPresences < NrowAbsences){
  for(i in 1:NRUNS){
    
  Species1_Subset <- rbind(Species1 %>% filter(AbundanceAdult40 > 0), Species1 %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
  
  # Extract residuals from fitted values. 
  GLMM_Resid <- resid(glmmTMB(BaselineFormula,
                              data = Species1_Subset, family = gaussian, zi = ~ 1))
  # Fit model to residuals
  Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = k), qu = q, data = Species1_Subset)
  TempData$Abundance <- predict(Species1_GAM, TempData, se = T)$fit
  Topt[i]         <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  T_Gam_pvalue[i] <- summary(Species1_GAM)$s.table[4]
  T_Gam_edf[i]   <- summary(Species1_GAM)$s.table[1]
  }
    
    
    }else{ # If not then use single model. 
   
     # Extract residuals from fitted values. 
    GLMM_Resid <- resid(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
                                data = Species1, family = gaussian, zi = ~ 1))
    # Fit model to residuals
    Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = k), qu = q, data = Species1)
    TempData$Abundance <- predict(Species1_GAM, TempData, se = T)$fit
    Topt         <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
    T_Gam_pvalue <- summary(Species1_GAM)$s.table[4]
    T_Gam_edf    <- summary(Species1_GAM)$s.table[1]
    
    }
  rm(TempData)
  
  # Extract temperture at maximum abundance, and maximum abundance at optimum temperature. 
  #Temps <- Species1 %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanTemp_CoralWatch_Scaled), MaxTemp = max(.$MeanTemp_CoralWatch_Scaled)) %>% unnest()
  #TempData <- data.frame(MeanTemp_CoralWatch_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 1000))
  #AbunPredictions <- predict(Species1_GAM, TempData, se = T)
  #TempData$Abundance <- AbunPredictions$fit
  #TempData$SE <- AbunPredictions$se
  #TempData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  #TempData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  
  # Estimate optimum temperature based on maximum abundance predicted 
  #Topt_Raw <- TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled']# * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  #Topt <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  
  # Take SE at Topt
  #SE <- TempData$SE[which.min(abs((TempData$MeanTemp_CoralWatch_Scaled*Species1$MeanTemp_CoralWatch_scaleBT[1] + Species1$MeanTemp_CoralWatch_centerBT[1]) - Topt))]
  #(quantile(rnorm(mean = Topt_Raw, sd = SE, n = 10000), c(0.025, 0.975)) * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  
  # Estimate variation in temperature below thermal optimum. OLD METHOD. 
  # Tsd <- Species1 %>% filter(MeanTemp_CoralWatch < Topt) %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)
  # Tsd_V2 <- Species1 %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)/2
  
  ### --- Estimate max abundance
  MaxAbundance <- Species1 %>% .$MaxAbundance %>% unique
  
  ### --- Estimate confidence score criteria 4. 
  T_Opt_Difference_Upper <- max(Species1$MeanTemp_CoralWatch) - mean(Topt)
  T_Opt_Difference_Lower <- mean(Topt) - min(Species1$MeanTemp_CoralWatch)
  
  
  ### --- Estimate confidence score criteria 5. 
  # Example extraction
  #T_Gam_pvalue <- summary(Species1_GAM)$s.table[4] # < 0.05
  #T_Gam_edf    <- summary(Species1_GAM)$s.table[1] # > 1
  
  ### --- Create output dataframe of predictions from model to save with data (and plot up later)
  PredData <- Species1 %>% filter(AbundanceAdult40 > 0) %>% 
    do(SpeciesName = unique(.$SpeciesName), 
       MeanTemp_CoralWatch        = seq(min(.$MeanTemp_CoralWatch), max(.$MeanTemp_CoralWatch), length.out = 100), 
       MeanTemp_CoralWatch_Scaled = seq(min(.$MeanTemp_CoralWatch_Scaled), max(.$MeanTemp_CoralWatch_Scaled), length.out = 100)) %>% 
    unnest(SpeciesName) %>% unnest(MeanTemp_CoralWatch, MeanTemp_CoralWatch_Scaled)
  
  AbunPredictions <- predict(Species1_GAM, PredData, se = T)
  PredData$Abundance <- AbunPredictions$fit
  PredData$SE  <- AbunPredictions$se
  PredData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  PredData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  PredData$Topt <- mean(Topt)
  PredData$Topt_SD <- sd(Topt)
  PredData$MaxAbundance <- MaxAbundance
  PredData$q <- q
  PredData$k <- k
  
  OutputData <- list(model = Species1_GAM, predictions = PredData)
  
  ### --- Save qgam outputs for if needed later
  dir.create(paste0('data_derived2/', Folder))
  saveRDS(OutputData, file = paste(paste0('data_derived2/', Folder, '/'), gsub(' ' , '_', unique(Species1$SpeciesName)), gsub('-','_',Sys.Date()), 'q_is', gsub('0.','',q),'k_is',k,'.rds', sep = '_'))
  rm(Species1_GAM)
  
  Topt2 <- mean(Topt)
  Topt_SD <- sd(Topt)
  return(data_frame(SpeciesName = unique(Species1$SpeciesName), 
                    Topt = Topt2, 
                    Topt_SD = Topt_SD,
                    MaxAbundance = MaxAbundance, 
                    T_Opt_Difference_Upper = T_Opt_Difference_Upper,
                    T_Opt_Difference_Lower = T_Opt_Difference_Lower,
                    T_Gam_pvalue = mean(T_Gam_pvalue), 
                    T_Gam_edf = mean(T_Gam_edf)))
}

TestSpecies <- RLS_20 %>% filter(SpeciesName == unique(.$SpeciesName)[349]) # %>% ConstrainAbsences(.)
TestSpecies <- RLS_20 %>% filter(SpeciesName == Bottom10[9])# %>% ConstrainAbsences(.)
FitModels_27092015(Species1 = TestSpecies, Folder = 'TEST')

# save.image('data_derived/SaveObjectsScript2-2018-27-09.RData')
# load('data_derived/SaveObjectsScript2-2018-27-09.RData')

# Fitting models for all species of the 'high confidence species' from previous version ----
# Load in old data 
#load(file = 'data_derived/objects-from-script-2-05062018.RData')

#TestSetSpecies <- ThermalNicheData_New %>% filter(ConfidenceCombined == 3) %>% .$SpeciesName #%>% sample(., 50)
#TestSetSpecies <- RLS_20 %>% filter(SpeciesName %in% TestSetSpecies)

# Runs function in parallel. 
cl <- makeCluster(4)
registerDoParallel(cl)
ModelOutputs <- foreach(i=1:length(unique(RLS_20$SpeciesName)), .packages=c('tidyr', 'dplyr', 'glmmTMB', 'qgam')) %dopar% {
  tryCatch(FitModels_27092015(RLS_20[which(RLS_20$SpeciesName == unique(RLS_20$SpeciesName)[i]),], 
                     Folder = 'AllModelsNew', NRUNS = 5), error = function(e) NA)
}
stopCluster()

QgamModels <- do.call(rbind, ModelOutputs)
#saveRDS(QgamModels, file = 'data_derived/qGamModelOutputs-2019-09-29.rds')

# Match up new and old Topts 
OldTops <- ThermalNicheData_New #%>% filter(SpeciesName %in% unique(TestSetSpecies$SpeciesName))
OldTops$Topt_OLD <- OldTops$Topt
NewTops <- left_join(NewTops, OldTops %>% dplyr::select(SpeciesName, Topt_OLD))

# Plot models and see relationship. 
pdf(file = 'figures_extra/ComparisonTopt-NewVsOld.pdf', height = 5, width = 5)
plot(NewTops$Topt, NewTops$Topt_OLD, xlab = 'New Topt', ylab = 'Old Topt')
abline(0, 1)
dev.off()
summary(lm(NewTops$Topt ~ NewTops$Topt_OLD))

# Dealing with species for which models do not run ----

# Work through function above to figure out why these species' fail and if this can be managed. 
# These are the species with errors in the PCA scores. 
# Not all species are missing PC data. # Kyphosus sectatrix had 2 NAs which caused problem. Now fixed. 

# Missing species in previous run didn't work because of error in PCA code which threw an NA when all data were same
# in a single column. This was all 0s in reef column. This column is now removed for these species. 

# Subset data for species which did no run. 
MissingSpeciesQgam <- unique(RLS_20$SpeciesName)[!QgamModels$SpeciesName %in% unique(RLS_20$SpeciesName)]
RLS_20_MissingSpecies <- RLS_20 %>% filter(SpeciesName %in% MissingSpeciesQgam)

# Run these models in parallel
cl <- makeCluster(4)
registerDoParallel(cl)
ModelOutputs_MissingSpecies <- foreach(i=1:length(unique(RLS_20_MissingSpecies$SpeciesName)), .packages=c('tidyr', 'dplyr', 'glmmTMB', 'qgam')) %dopar% {
  tryCatch(FitModels_27092015(RLS_20_MissingSpecies[which(RLS_20_MissingSpecies$SpeciesName == unique(RLS_20_MissingSpecies$SpeciesName)[i]),], 
                              Folder = 'AllModelsNew', NRUNS = 50), error = function(e) NA)
}
stopCluster()

QgamModels_MissingSpecies <- do.call(rbind, ModelOutputs_MissingSpecies)
saveRDS(QgamModels_MissingSpecies, file = 'data_derived/QgamModelsModelOutputs_MissingSpecies-2019-09-28.rds')

# Bind together the two outputs ----
AllQgams <- rbind(QgamModels[!is.na(QgamModels$Topt),], QgamModels_MissingSpecies)
Quantile_Parameters <- AllQgams
#save(Quantile_Parameters, file = 'data_derived2/AllQgams_2018-09-28.RData')
# ---- 



# REDUNDANT SECTION NOW: ***MODELLING IS DONE ABOVE AS OF 2019-09-28*** MODELLING QGAMS SECTION ----
# Function which fits quantile gams and extracts parameters of interest. ----

# Function to fit model with different covariate structures depending on N
FitQuantileGam_Covariates <- function(TestSpp, k = 4, q = 0.9, Folder = NULL){
  
  NrowPresences <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences <- TestSpp %>% filter(AbundanceAdult40 == 0) %>% nrow(.)
  
  if(NrowPresences <= NrowAbsences){
    TestSpp <- rbind(TestSpp %>% filter(AbundanceAdult40 > 0), TestSpp %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
  }else{NULL}
  
  # print(nrow(TestSpp))
  
  # Fit generalized additive model
  #Create columns for modelling
  TestSpp$AbundanceAdult40_log <- log(TestSpp$AbundanceAdult40 + 1)
  TestSpp$SamplingIntensity_Scaled <- as.numeric(scale(TestSpp$SamplingIntensity))
  TestSpp$MeanTemp_CoralWatch_Scaled <- as.numeric(scale(TestSpp$MeanTemp_CoralWatch))
  TestSpp$ReefAreaIn15km_Scaled <- ifelse(sum(TestSpp$ReefAreaIn15km) > 0, as.numeric(scale(TestSpp$ReefAreaIn15km)), TestSpp$ReefAreaIn15km)
  TestSpp$HumPop50km2_Scaled <- as.numeric(scale(TestSpp$HumPop50km2))
  TestSpp$Depth_Site <- as.numeric(scale(TestSpp$Depth_Site))
  TestSpp$npp_mean_Scaled <- as.numeric(scale(TestSpp$npp_mean))
  TestSpp$MeanTemp_CoralWatch_scaleBT <- attr(scale(TestSpp$MeanTemp_CoralWatch), "scaled:scale")
  TestSpp$MeanTemp_CoralWatch_centerBT <- attr(scale(TestSpp$MeanTemp_CoralWatch), "scaled:center")
  TestSpp$NEOLI[is.na(TestSpp$NEOLI)] <- 0
  
  
  # If the number of data is > 100 then fit the below models, else fit a simpler set
  if(NrowPresences > 100){
    
  if(is.numeric(k)){
    TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp), error = function(e) NA)
    
    i <- 1
    while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
      class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp), error = function(e) NA)
      i <-  i + 1}
    
    if(!is.na(TestSpp_gam11)){
    while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp)}
    }
    }else{
    
    TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled  + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp), error = function(e) NA)
    
    i <- 1
    while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
      class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp), error = function(e) NA)
      i <-  i + 1}
    
    if(!is.na(TestSpp_gam11)){
    while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled + NEOLI + npp_mean_Scaled + HumPop50km2_Scaled + ReefAreaIn15km_Scaled + Depth_Site, qu = q, data = TestSpp)}
    }
    }
    
  }else{
    
    if(is.numeric(k)){
      TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
      
      i <- 1
      while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
        class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
        i <-  i + 1}
      
      if(!is.na(TestSpp_gam11)){
      while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp)}
      }
        }else{
      
      TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
      
      i <- 1
      while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
        class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
        i <-  i + 1}
    
    }
    } # End of if statement based on number of presences and which model to fit. 
  
  if(is.na(TestSpp_gam11)){
    
    if(is.numeric(k)){
      TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
      
      i <- 1
      while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
        class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
        i <-  i + 1}
      
      while(TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)){ TestSpp_gam11 <- qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled, k = k) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp)}
    }else{
      
      TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled  + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
      
      i <- 1
      while(#TestSpp_gam11$converged == F | TestSpp_gam11$outer.info$conv == F | is.logical(TestSpp_gam11)
        class(TestSpp_gam11)[1] != 'qgam' && i < 10){TestSpp_gam11 <- tryCatch(qgam(AbundanceAdult40_log ~ s(MeanTemp_CoralWatch_Scaled) + SamplingIntensity_Scaled + NEOLI, qu = q, data = TestSpp), error = function(e) NA)
        i <-  i + 1}
      
    }
  }else{NULL} # End of if statement based on non-fitting of above models but 100 sites. 

  # Extract temperture at maximum abundance, and maximum abundance at optimum temperature. 
  Temps <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanTemp_CoralWatch_Scaled), MaxTemp = max(.$MeanTemp_CoralWatch_Scaled)) %>% unnest()
  TempData <- data.frame(MeanTemp_CoralWatch_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 50))
  TempData$SamplingIntensity_Scaled <- mean(TestSpp$SamplingIntensity_Scaled, na.rm = T)
  TempData$npp_mean_Scaled <- mean(TestSpp$npp_mean_Scaled, na.rm = T) 
  TempData$HumPop50km2_Scaled <- mean(TestSpp$HumPop50km2_Scaled, na.rm = T)
  TempData$ReefAreaIn15km_Scaled <- mean(TestSpp$ReefAreaIn15km_Scaled, na.rm = T)
  TempData$Depth_Site <- mean(TestSpp$Depth_Site, na.rm = T)
  TempData$NEOLI <- mean(TestSpp$NEOLI, na.rm = T)
  AbunPredictions <- predict(TestSpp_gam11, TempData, se = T)
  TempData$Abundance <- AbunPredictions$fit
  TempData$SE <- AbunPredictions$se
  TempData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  TempData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  
  # Estimate optimum temperature based on maximum abundance predicted 
  Topt <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * TestSpp$MeanTemp_CoralWatch_scaleBT[1]) + TestSpp$MeanTemp_CoralWatch_centerBT[1]
  
  # Estimate variation in temperature below thermal optimum. OLD METHOD. 
  # Tsd <- TestSpp %>% filter(MeanTemp_CoralWatch < Topt) %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)
  # Tsd_V2 <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)/2
  
  ### --- Estimate max abundance
  MaxAbundance <- TestSpp %>% .$MaxAbundance %>% unique
  
  ### --- Estimate confidence score criteria 4. 
  T_Opt_Difference_Upper <- max(TestSpp$MeanTemp_CoralWatch) - Topt 
  T_Opt_Difference_Lower <- Topt - min(TestSpp$MeanTemp_CoralWatch)
  
  
  ### --- Estimate confidence score criteria 5. 
  # Example extraction
  T_Gam_pvalue <- summary(TestSpp_gam11)$s.table[4]# < 0.05
  T_Gam_edf <- summary(TestSpp_gam11)$s.table[1]# > 1
  
  ### --- Create output dataframe of predictions from model to save with data (and plot up later)
  PredData <- TestSpp %>% filter(AbundanceAdult40 > 0) %>% 
    do(SpeciesName = unique(.$SpeciesName), 
       MeanTemp_CoralWatch        = seq(min(.$MeanTemp_CoralWatch), max(.$MeanTemp_CoralWatch), length.out = 100), 
       MeanTemp_CoralWatch_Scaled = seq(min(.$MeanTemp_CoralWatch_Scaled), max(.$MeanTemp_CoralWatch_Scaled), length.out = 100)) %>% 
    unnest(SpeciesName) %>% unnest(MeanTemp_CoralWatch, MeanTemp_CoralWatch_Scaled)
  
  PredData$SamplingIntensity_Scaled <- median(TestSpp$SamplingIntensity_Scaled, na.rm = T)
  PredData$npp_mean_Scaled <- mean(PredData$npp_mean_Scaled, na.rm = T) 
  PredData$HumPop50km2_Scaled <- mean(PredData$HumPop50km2_Scaled, na.rm = T)
  PredData$ReefAreaIn15km_Scaled <- mean(PredData$ReefAreaIn15km_Scaled, na.rm = T)
  PredData$Depth_Site <- mean(PredData$Depth_Site, na.rm = T)
  PredData$NEOLI <- mean(TestSpp$NEOLI, na.rm = T)
  
  AbunPredictions <- predict(TestSpp_gam11, PredData, se = T)
  PredData$Abundance <- AbunPredictions$fit
  PredData$SE  <- AbunPredictions$se
  PredData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  PredData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  PredData$Topt <- Topt
  PredData$MaxAbundance <- MaxAbundance
  PredData$q <- q
  PredData$k <- k
  
  OutputData <- list(model = TestSpp_gam11, predictions = PredData, rawdata = TestSpp %>% select(SpeciesName, AbundanceAdult40_log,  AbundanceAdult40, MeanTemp_CoralWatch))
  
  ### --- Save qgam outputs for if needed later
  dir.create(paste0('data_derived2/', Folder))
  saveRDS(OutputData, file = paste(paste0('data_derived2/', Folder, '/'), gsub(' ' , '_', unique(TestSpp$SpeciesName)), gsub('-','_',Sys.Date()), 'q_is', gsub('0.','',q),'k_is',k,'.rds', sep = '_'))
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

# test on 1 species
TestSpp <- RLS_19 %>% filter(SpeciesName == 'Labroides dimidiatus')


TestSpp <- RLS_19 %>% filter(SpeciesName %in% unique(.$SpeciesName)[NA_Models[4]])
FitQuantileGam_Covariates(TestSpp, k = 4, q = 0.8, Folder = 'Test')
TestSppMod <- readRDS(file = 'data_derived2/Run1/_Labroides_dimidiatus_2018_09_25_q_is_8_k_is_4_.rds')
summary(TestSppMod$model)
plot(TestSppMod$model)

# Fit qgams in dplyr. Each individual model fits independently 
All_QGams <- RLS_19 %>% group_by(SpeciesName) %>% do(GamThermalNiche = tryCatch(FitQuantileGam_Covariates(., q = 0.8, Folder = 'RunMultipleCovariates'), error = function(e) NA))

# Check for NAs. 
NA_Models <- which(is.na(All_QGams$GamThermalNiche)) # Some models produce errors when fit with the above code? Too complex? 
RLS_19 %>% filter(SpeciesName %in% unique(.$SpeciesName)[NA_Models], Presence == 1) %>% nrow(.) / length(NA_Models) # Could be due to low sample size. 

# Refit models to the previous models which were NAs. 
NA_QGams <- RLS_19 %>% 
  filter(SpeciesName %in% unique(.$SpeciesName)[NA_Models]) %>%  
  group_by(SpeciesName) %>% 
  do(GamThermalNiche = tryCatch(FitQuantileGam_Covariates(., q = 0.8, Folder = 'RunMultipleCovariates'), error = function(e) NA))

# Combine lists 
All_QGams <- rbind(All_QGams[-which(is.na(All_QGams$GamThermalNiche)),], NA_QGams[-which(is.na(NA_QGams$GamThermalNiche)),])

# 2 are NAs without enough data to model. 

# Save the output object. 
saveRDS(All_QGams, file = 'data_derived2/All_QGams_2018-08-25.rds')

# Obtain parameters from the quantile gam models. 
Quantile_Parameters <- do.call(rbind, All_QGams$GamThermalNiche)

# How many have non-significant estimates of Topt = 185. Report in manuscript. 
Quantile_Parameters[which(Quantile_Parameters$T_Gam_pvalue > 0.05),]
185/704 #26% show non-significant temperature patterns. 

# FOR REVIEW compare the Topt from the old thermal niche estimates to the new thermal niche estimates when modelled with covariates ----
# Read in the old version of ThermalNicheData and plot. 
load(file = 'data_derived/objects-from-script-2-05062018.RData')

Quantile_Parameters$Topt2 <- Quantile_Parameters$Topt
ThermalNicheData_New_Comparison <- left_join(ThermalNicheData_New, Quantile_Parameters %>% dplyr::select(SpeciesName, Topt2))
ThermalNicheData_New_Comparison2 <- ThermalNicheData_New_Comparison %>% filter(ConfidenceCombined == 3)

pdf('figure_final/review-figures/Topt-comparisons.pdf', width = 3, height = 3)
ggplot(ThermalNicheData_New_Comparison2) + 
  geom_point(aes(x = Topt, y = Topt2)) + 
  geom_abline() + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank()) + 
  xlab('Topt - Temperature only') + ylab('Topt - Including covariates')
dev.off()

ToptModel <- lm(Topt2 ~ Topt, data = ThermalNicheData_New_Comparison2)
summary(ToptModel)

# Fit the qgam models to a quantile of 0.95 -----

# Fit qgams in dplyr. Each individual model fits independently 
All_QGams_0.95 <- RLS_19 %>% group_by(SpeciesName) %>% do(GamThermalNiche = tryCatch(FitQuantileGam(., q = 0.95), error = function(e) NA))

# Check for NAs. 
NA_Models <- which(is.na(All_QGams_0.95$GamThermalNiche)) 
All_QGams_0.95 <- All_QGams_0.95[-NA_Models, ]

# Save the output object. 
saveRDS(All_QGams_0.95, file = 'data_derived2/All_QGams_0.95_2018-02-27.rds')

# Obtain parameters from the quantile gam models. 
Quantile_Parameters_0.95 <- do.call(rbind, All_QGams_0.95$GamThermalNiche)

# How many have non-significant estimates of Topt
Quantile_Parameters_0.95[which(Quantile_Parameters_0.95$T_Gam_pvalue > 0.05),]

# Step 1. Load object created from the above code and run from 166:255 below ----
All_QGams <- readRDS(file = 'data_derived2/All_QGams_2018-08-25.rds')
# Combine and plot estimates of Topt for SOM ----
All_QGams_0.95 
All_QGams
Quantile_Parameters_0.95 <- do.call(rbind, All_QGams_0.95$GamThermalNiche)
names(Quantile_Parameters_0.95) <- paste(names(Quantile_Parameters_0.95), '_95', sep = '')
Quantile_Parameters <- do.call(rbind, All_QGams$GamThermalNiche)
Quantile_Parameters <- Quantile_Parameters[which(Quantile_Parameters$SpeciesName %in% Quantile_Parameters_0.95$SpeciesName_95),]
Quantile_Parameters_combined <- cbind(Quantile_Parameters_0.95, Quantile_Parameters)

pdf('figure_final/SOM-comparison of quantiles.pdf', width = 5, height = 5)
ggplot(Quantile_Parameters_combined) + 
  geom_point(aes(x = Topt, y = Topt_95)) + 
  theme_bw() + 
  xlab(expression(paste(T["opt"], "  (0.8 quantile)"))) + 
  ylab(expression(paste(T["opt"], "  (0.95 quantile)"))) + 
  theme(aspect.ratio = 0.75, panel.grid.major = element_blank() , panel.grid.minor = element_blank())
dev.off()

cor.test(Quantile_Parameters_combined$Topt_95, Quantile_Parameters_combined$Topt, method = 'pearson')

# Summarise model parameters for coefficient terms: NEOLI ----
QgamModels <- paste0('data_derived2/RunMultipleCovariates/', list.files('data_derived2/RunMultipleCovariates'))
AllModels <- list()
for(i in 1:length(QgamModels)){ 
  AllModels[[i]] <- readRDS(QgamModels[i])
}

as.numeric(coef(AllModels[[1]]$model)[c('NEOLI', 'npp_mean_Scaled', 'HumPop50km2_Scaled', 'ReefAreaIn15km_Scaled')])
as.numeric(summary(AllModels[[1]]$model)$pTerms.table['NEOLI','p-value'])

NEOLI_Coefs <- lapply(AllModels, function(x) data.frame(SpeciesName = x$predictions$SpeciesName[1],
                                                        coefs = as.numeric(coef(x$model)['NEOLI']), 
                                                        p_values = as.numeric(summary(x$model)$pTerms.table['NEOLI','p-value'])))
NEOLI_Coefs <- do.call('rbind', NEOLI_Coefs)
hist(NEOLI_Coefs$coefs)
hist(NEOLI_Coefs$p_values)
sum(NEOLI_Coefs$p_values < 0.05) / length(AllModels) # just 24% are significant. 
NEOLI_P0.05 <- NEOLI_Coefs$coefs[which(NEOLI_Coefs$p_values  < 0.05 & NEOLI_Coefs$p_values  > 0.01)]
NEOLI_P0.01 <- NEOLI_Coefs$coefs[which(NEOLI_Coefs$p_values  < 0.01 & NEOLI_Coefs$p_values  > 0.001)]
NEOLI_P0.001 <- NEOLI_Coefs$coefs[which(NEOLI_Coefs$p_values < 0.001)]

pdf('figure_final/review-figures/NEOLI-coefs.pdf', width = 5, height = 12)
gridExtra::grid.arrange(
ggplot() + 
  geom_histogram(data = data.frame(x = NEOLI_P0.001), aes(x = NEOLI_P0.001), alpha = 0.8, fill = 'gray10', position = 'identity') + 
  theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
  xlab('NEOLI coefficient') + 
  ylab(NULL) + 
  xlim(range(c(NEOLI_P0.05, NEOLI_P0.01, NEOLI_P0.001))), 

ggplot() + 
  geom_histogram(data = data.frame(x = NEOLI_P0.01), aes(x = NEOLI_P0.01), alpha = 0.5, fill = 'gray20', position = 'identity') + 
  theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
  xlab('NEOLI coefficient') + 
  ylab(NULL) + 
  xlim(range(c(NEOLI_P0.05, NEOLI_P0.01, NEOLI_P0.001))), 

ggplot() + 
  geom_histogram(data = data.frame(x = NEOLI_P0.05), aes(x = NEOLI_P0.05), alpha = 0.2, fill = 'gray30', position = 'identity') + 
  theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
  xlab('NEOLI coefficient') + 
  ylab(NULL) + 
  xlim(range(c(NEOLI_P0.05, NEOLI_P0.01, NEOLI_P0.001))), 
nrow = 3)
dev.off()

t.test(NEOLI_Coefs[,1])
t.test(NEOLI_P0.001) # Those that are significant are positive. 
t.test(NEOLI_P0.01) # Those that are significant are positive. 
t.test(NEOLI_P0.05) # Those that are significant are positive. 


# Summarise model parameters for coefficient terms: NPP ----

grepl('npp_mean_Scaled', formula(AllModels[[1]]$model)[3])
as.numeric(coef(AllModels[[1]]$model)[c('NEOLI', 'npp_mean_Scaled', 'HumPop50km2_Scaled', 'ReefAreaIn15km_Scaled')])
as.numeric(summary(AllModels[[1]]$model)$pTerms.table['NEOLI','p-value'])

# Grab those that have the covariate of interest
NPP_Models <- do.call(c, lapply(AllModels, function(x) grepl('npp_mean_Scaled', formula(x$model)[3])))
NPP_Coefs <- lapply(AllModels[which(NPP_Models==T)], function(x) data.frame(SpeciesName = x$predictions$SpeciesName[1],
                                                                            coefs = as.numeric(coef(x$model)['npp_mean_Scaled']),
                                                                            p_values = as.numeric(summary(x$model)$pTerms.table['npp_mean_Scaled','p-value'])))
NPP_Coefs <- do.call('rbind', NPP_Coefs)
hist(NPP_Coefs$coefs)
hist(NPP_Coefs$p_values)
sum(NPP_Coefs$p_values < 0.05) / length(AllModels) # just 28% are significant. 
NPP_P0.05 <- NPP_Coefs$coefs[which(NPP_Coefs$p_values  < 0.05 & NPP_Coefs$p_values  > 0.01)]
NPP_P0.01 <- NPP_Coefs$coefs[which(NPP_Coefs$p_values  < 0.01 & NPP_Coefs$p_values  > 0.001)]
NPP_P0.001 <- NPP_Coefs$coefs[which(NPP_Coefs$p_values < 0.001)]

pdf('figure_final/review-figures/NPP-coefs.pdf', width = 5, height = 12)
gridExtra::grid.arrange(
  ggplot() + 
    geom_histogram(data = data.frame(x = NPP_P0.001), aes(x = NPP_P0.001), alpha = 0.8, fill = 'gray10', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('NPP coefficient') + 
    ylab(NULL) + 
    xlim(range(c(NPP_P0.05, NPP_P0.01, NPP_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = NPP_P0.01), aes(x = NPP_P0.01), alpha = 0.5, fill = 'gray20', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('NPP coefficient') + 
    ylab(NULL) + 
    xlim(range(c(NPP_P0.05, NPP_P0.01, NPP_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = NPP_P0.05), aes(x = NPP_P0.05), alpha = 0.2, fill = 'gray30', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('NPP coefficient') + 
    ylab(NULL) + 
    xlim(range(c(NPP_P0.05, NPP_P0.01, NPP_P0.001))), 
  nrow = 3)
dev.off()

t.test(NPP_Coefs[,1])
t.test(NPP_P0.001)
t.test(NPP_P0.01)
t.test(NPP_P0.05) 



# Summarise model parameters for coefficient terms: HUMAN ----

grepl('HumPop50km2_Scaled', formula(AllModels[[1]]$model)[3])
as.numeric(coef(AllModels[[1]]$model)[c('NEOLI', 'npp_mean_Scaled', 'HumPop50km2_Scaled', 'ReefAreaIn15km_Scaled')])
as.numeric(summary(AllModels[[1]]$model)$pTerms.table['NEOLI','p-value'])

# Grab those that have the covariate of interest
HUMAN_Models <- do.call(c, lapply(AllModels, function(x) grepl('HumPop50km2_Scaled', formula(x$model)[3])))
HUMAN_Coefs <- lapply(AllModels[which(HUMAN_Models==T)], function(x) data.frame(SpeciesName = x$predictions$SpeciesName[1],
                                                                                coefs = as.numeric(coef(x$model)['HumPop50km2_Scaled']), 
                                                                                p_values = as.numeric(summary(x$model)$pTerms.table['HumPop50km2_Scaled','p-value'])))
HUMAN_Coefs <- do.call('rbind', HUMAN_Coefs)
hist(HUMAN_Coefs$coefs)
hist(HUMAN_Coefs$p_values)
sum(HUMAN_Coefs$p_values < 0.05) / length(AllModels) # just 14% are significant. 
HUMAN_P0.05 <- HUMAN_Coefs$coefs[which(HUMAN_Coefs$p_values  < 0.05 & HUMAN_Coefs$p_values  > 0.01)]
HUMAN_P0.01 <- HUMAN_Coefs$coefs[which(HUMAN_Coefs$p_values  < 0.01 & HUMAN_Coefs$p_values  > 0.001)]
HUMAN_P0.001 <- HUMAN_Coefs$coefs[which(HUMAN_Coefs$p_values < 0.001)]

pdf('figure_final/review-figures/HUMAN-coefs.pdf', width = 5, height = 12)
gridExtra::grid.arrange(
  ggplot() + 
    geom_histogram(data = data.frame(x = HUMAN_P0.001), aes(x = HUMAN_P0.001), alpha = 0.8, fill = 'gray10', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('HUMAN coefficient') + 
    ylab(NULL) + 
    xlim(range(c(HUMAN_P0.05, HUMAN_P0.01, HUMAN_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = HUMAN_P0.01), aes(x = HUMAN_P0.01), alpha = 0.5, fill = 'gray20', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('HUMAN coefficient') + 
    ylab(NULL) + 
    xlim(range(c(HUMAN_P0.05, HUMAN_P0.01, HUMAN_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = HUMAN_P0.05), aes(x = HUMAN_P0.05), alpha = 0.2, fill = 'gray30', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('HUMAN coefficient') + 
    ylab(NULL) + 
    xlim(range(c(HUMAN_P0.05, HUMAN_P0.01, HUMAN_P0.001))), 
  nrow = 3)
dev.off()

t.test(HUMAN_Coefs[,1])
t.test(HUMAN_P0.001)
t.test(HUMAN_P0.01)
t.test(HUMAN_P0.05) 

# Summarise model parameters for coefficient terms: REEF_AREA ----

grepl('ReefAreaIn15km_Scaled', formula(AllModels[[1]]$model)[3])
as.numeric(coef(AllModels[[1]]$model)[c('NEOLI', 'npp_mean_Scaled', 'HumPop50km2_Scaled', 'ReefAreaIn15km_Scaled')])
as.numeric(summary(AllModels[[1]]$model)$pTerms.table['NEOLI','p-value'])

# Grab those that have the covariate of interest
REEF_AREA_Models <- do.call(c, lapply(AllModels, function(x) grepl('ReefAreaIn15km_Scaled', formula(x$model)[3])))
REEF_AREA_Coefs <- lapply(AllModels[which(REEF_AREA_Models==T)], function(x) data.frame(SpeciesName = x$predictions$SpeciesName[1],
                                                                                        coefs = as.numeric(coef(x$model)['ReefAreaIn15km_Scaled']),
                                                                                        p_values = as.numeric(summary(x$model)$pTerms.table['ReefAreaIn15km_Scaled','p-value'])))
REEF_AREA_Coefs <- do.call('rbind', REEF_AREA_Coefs)
hist(REEF_AREA_Coefs$coefs)
hist(REEF_AREA_Coefs$p_values)
sum(REEF_AREA_Coefs$p_values < 0.05) / length(AllModels) # just 23% are significant. 
REEF_AREA_P0.05 <- REEF_AREA_Coefs$coefs[which(REEF_AREA_Coefs$p_values  < 0.05 & REEF_AREA_Coefs$p_values  > 0.01)]
REEF_AREA_P0.01 <- REEF_AREA_Coefs$coefs[which(REEF_AREA_Coefs$p_values  < 0.01 & REEF_AREA_Coefs$p_values  > 0.001)]
REEF_AREA_P0.001 <- REEF_AREA_Coefs$coefs[which(REEF_AREA_Coefs$p_values < 0.001)]

pdf('figure_final/review-figures/REEF_AREA-coefs.pdf', width = 5, height = 12)
gridExtra::grid.arrange(
  ggplot() + 
    geom_histogram(data = data.frame(x = REEF_AREA_P0.001), aes(x = REEF_AREA_P0.001), alpha = 0.8, fill = 'gray10', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('REEF_AREA coefficient') + 
    ylab(NULL) + 
    xlim(range(c(REEF_AREA_P0.05, REEF_AREA_P0.01, REEF_AREA_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = REEF_AREA_P0.01), aes(x = REEF_AREA_P0.01), alpha = 0.5, fill = 'gray20', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('REEF_AREA coefficient') + 
    ylab(NULL) + 
    xlim(range(c(REEF_AREA_P0.05, REEF_AREA_P0.01, REEF_AREA_P0.001))), 
  
  ggplot() + 
    geom_histogram(data = data.frame(x = REEF_AREA_P0.05), aes(x = REEF_AREA_P0.05), alpha = 0.2, fill = 'gray30', position = 'identity') + 
    theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()) + 
    xlab('REEF_AREA coefficient') + 
    ylab(NULL) + 
    xlim(range(c(REEF_AREA_P0.05, REEF_AREA_P0.01, REEF_AREA_P0.001))), 
  nrow = 3)
dev.off()

t.test(REEF_AREA_Coefs[,1])
t.test(REEF_AREA_P0.001)
t.test(REEF_AREA_P0.01)
t.test(REEF_AREA_P0.05) 
# Combine coefs and see if there are correlations ----
AllCoefs <- join_all(list(NEOLI_Coefs, REEF_AREA_Coefs, HUMAN_Coefs, NPP_Coefs), by='SpeciesName', type='left')

# Not really any strong correlations here
pairs(AllCoefs[,c(2, 4, 6, 8)])

# REDUNDANT SECTION NOW: ***MODELLING IS DONE ABOVE AS OF 2019-09-28*** MODELLING QGAMS SECTION ----




# DEFINING THERMAL GUILDS ----
# Define thermal guilds based on Topt from the qgam models. Creation of RLS_All, RLS_Trop, and RLS_Temp objects ----

# Redefine thermal guild based on qgam estimates of thermal niche. 23°C from Stuart-Smith et al. 2015 Nature, and Stuart-Smith et al. 2017 Nature Ecology and Evolution. 
ThermalGuilds <- left_join(RLS_20, Quantile_Parameters %>% select(SpeciesName, Topt)) %>% 
  group_by(SpeciesName) %>% 
  do(ThermalGuild = ifelse(.$Topt > 23, 'tropical', 'temperate')) %>% 
  unnest(ThermalGuild) %>% unique()

# Estimate mean temperature across all observations used in future. 
MeanTemperatures <- RLS_20 %>% 
  group_by(SpeciesName) %>% 
  filter(AbundanceAdult40 > 0) %>%  
  do(T_Mean_Obs = mean(.$MeanTemp_CoralWatch, na.rm = T)) %>% 
  unnest(T_Mean_Obs)

# Join in thermal guild character, and mean temperature across ranges. 
RLS_20 <- left_join(RLS_20, ThermalGuilds)
RLS_20 <- left_join(RLS_20, MeanTemperatures)

# Estimate observed thermal midpoints (same as T_Mean_Obs), but also observed upper and lower per species. 
ThermalNicheData_Obs <- RLS_20 %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanTemp_CoralWatch)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanTemp_CoralWatch)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanTemp_CoralWatch, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
  select(-data)

# Define tropical and temperate data. 
RLS_Trop <- RLS_20 %>% filter(ThermalGuild == 'tropical')
RLS_Temp <- RLS_20 %>% filter(ThermalGuild == 'temperate')

# Scale covarites for RLS_Trop and RLS_Temp and create RLS_All ----

# Note that reef area is already log transformed. 

# Scale and centre covariates for tropical species.  
RLS_Trop$MeanTemp_CoralWatch_Scaled <- scale(RLS_Trop$MeanTemp_CoralWatch)
RLS_Trop$HumPop50km2_Scaled <- scale(RLS_Trop$HumPop50km2)
RLS_Trop$npp_mean_Scaled <- scale(RLS_Trop$npp_mean)
RLS_Trop$SamplingIntensity_Scaled <- scale(RLS_Trop$SamplingIntensity)
RLS_Trop$Depth_Site_Scaled <- scale(RLS_Trop$Depth_Site)

RLS_Trop$MeanTemp_CoralWatch_centerBT <- attr(scale(RLS_Trop$MeanTemp_CoralWatch), 'scaled:center')
RLS_Trop$HumPop50km2_centerBT <- attr(scale(RLS_Trop$HumPop50km2), 'scaled:center')
RLS_Trop$npp_mean_centerBT <- attr(scale(RLS_Trop$npp_mean), 'scaled:center')
RLS_Trop$SamplingIntensity_centerBT <- attr(scale(RLS_Trop$SamplingIntensity), 'scaled:center')
RLS_Trop$Depth_Site_centerBT <- attr(scale(RLS_Trop$Depth_Site), 'scaled:center')

RLS_Trop$MeanTemp_CoralWatch_scaleBT <- attr(scale(RLS_Trop$MeanTemp_CoralWatch), 'scaled:scale')
RLS_Trop$HumPop50km2_scaleBT <- attr(scale(RLS_Trop$HumPop50km2), 'scaled:scale')
RLS_Trop$npp_mean_scaleBT <- attr(scale(RLS_Trop$npp_mean), 'scaled:scale')
RLS_Trop$SamplingIntensity_scaleBT <- attr(scale(RLS_Trop$SamplingIntensity), 'scaled:scale')
RLS_Trop$Depth_Site_scaleBT <- attr(scale(RLS_Trop$Depth_Site), 'scaled:scale')


# Scale and centre covariates for temperate species. 
RLS_Temp$MeanTemp_CoralWatch_Scaled <- scale(RLS_Temp$MeanTemp_CoralWatch)
RLS_Temp$HumPop50km2_Scaled <- scale(RLS_Temp$HumPop50km2)
RLS_Temp$npp_mean_Scaled <- scale(RLS_Temp$npp_mean)
RLS_Temp$SamplingIntensity_Scaled <- scale(RLS_Temp$SamplingIntensity)
RLS_Temp$Depth_Site_Scaled <- scale(RLS_Temp$Depth_Site)

RLS_Temp$MeanTemp_CoralWatch_centerBT <- attr(scale(RLS_Temp$MeanTemp_CoralWatch), 'scaled:center')
RLS_Temp$HumPop50km2_centerBT <- attr(scale(RLS_Temp$HumPop50km2), 'scaled:center')
RLS_Temp$npp_mean_centerBT <- attr(scale(RLS_Temp$npp_mean), 'scaled:center')
RLS_Temp$SamplingIntensity_centerBT <- attr(scale(RLS_Temp$SamplingIntensity), 'scaled:center')
RLS_Temp$Depth_Site_centerBT <- attr(scale(RLS_Temp$Depth_Site), 'scaled:center')

RLS_Temp$MeanTemp_CoralWatch_scaleBT <- attr(scale(RLS_Temp$MeanTemp_CoralWatch), 'scaled:scale')
RLS_Temp$HumPop50km2_scaleBT <- attr(scale(RLS_Temp$HumPop50km2), 'scaled:scale')
RLS_Temp$npp_mean_scaleBT <- attr(scale(RLS_Temp$npp_mean), 'scaled:scale')
RLS_Temp$SamplingIntensity_scaleBT <- attr(scale(RLS_Temp$SamplingIntensity), 'scaled:scale')
RLS_Temp$Depth_Site_scaleBT <- attr(scale(RLS_Temp$Depth_Site), 'scaled:scale')


# Check for NAs. 
sapply(RLS_Temp, function(x) sum(is.na(x))) # Check where the NAs are being introduced
sapply(RLS_Trop, function(x) sum(is.na(x)))

# Save objects
#saveRDS(RLS_Temp, file = 'data_derived2/RLS_withCovariates_Temperate_2018-09-26.rds')
#saveRDS(RLS_Trop, file = 'data_derived2/RLS_withCovariates_Tropical_2018-09-26.rds')

# Create RLS_All and save.
RLS_All <- rbind(RLS_Temp, RLS_Trop)
#saveRDS(RLS_All, file = 'data_derived2/RLS_withCovariates_All_2018-09-26.rds')







# (step 2) Load objects from above code if needed ----
#RLS_Temp <- readRDS(file = 'data_derived2/RLS_withCovariates_Temperate_2018-09-26.rds')
#RLS_Trop <- readRDS(file = 'data_derived2/RLS_withCovariates_Tropical_2018-09-26.rds')
#RLS_All  <- readRDS(file = 'data_derived2/RLS_withCovariates_All_2018-09-26.rds')
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
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanTemp_CoralWatch)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanTemp_CoralWatch)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanTemp_CoralWatch, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
  select(-data)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------



# MODELLING OCCUPANCY SECTION V2 ----
# Function to fit species' specific occupancy models ----
Species1 <- RLS_20 %>% filter(SpeciesName == Top10[2]) # unique(.$SpeciesName)['Labroides dimidiatus'])

SpeciesSpecific_OccupancyModel <- function(Species1, FOLDER = 'TESTFOLDER'){
  
  # HERE I SHOULD DO SPECIES' SPECIFIC SCALING OF PARAMETERS AGAIN AS IN THE QGAM FUNCTION.
  # Check before with a particular dataset - but at present there are 2 NAs which cause problems with 1 species. 
  Species1 <- Species1[!is.na(Species1$MeanTemp_CoralWatch),]
  
  # Scale and transform variables as appropriate
  Species1$AbundanceAdult40_log           <- as.numeric(log(Species1$AbundanceAdult40 + 1))
  Species1$SamplingIntensity_Scaled       <- as.numeric(scale(Species1$SamplingIntensity))
  Species1$MeanTemp_CoralWatch_Scaled     <- as.numeric(scale(Species1$MeanTemp_CoralWatch))
  Species1$Depth_Site_Scaled              <- as.numeric(scale(Species1$Depth_Site))
  Species1$NEOLI[is.na(Species1$NEOLI)]   <- 0
  
  # Columns for back transformation of temperature 
  # Columns for backtransformation of temperature data
  Species1$MeanTemp_CoralWatch_scaleBT  <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:scale")
  Species1$MeanTemp_CoralWatch_centerBT <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:center")
  
  # Remove NA columns as this crashes the gams. 
  Species1 <- Species1[,-which(apply(apply(Species1, 2, is.na), 2, sum) > 1)]
  
  # Check that species has adequate sampling above and below abundance limits to model occupancy.
  #if(!Species1$Confidence_Occ_Tupper[1] == 1){T_Upper <- NA}#else{NULL}
  #if(!Species1$Confidence_Occ_Tlower[1] == 1){T_Lower <- NA}#else{NULL}
  
  # Convert varying number of PCs to numerics  
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  PC_Columns <- apply(PC_Columns, 2, as.numeric)
  for(i in 1:ncol(PC_Columns)){
    Species1[,colnames(Species1[,grepl('PC', colnames(Species1))])[i]] <- PC_Columns[,i]
  }
  
  # Create a model structure that can be flexible to number of PCs. 
  BaselineFormula <- 'Presence ~ poly(MeanTemp_CoralWatch_Scaled, 2) + s(Depth_Site_Scaled, k = 4) + NEOLI'
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  for(i in 1:ncol(PC_Columns)){
    BaselineFormula <- paste0(BaselineFormula, ' + ', 's(', colnames(Species1[,grepl('PC', colnames(Species1))])[i], ',k=4)')
  }
  BaselineFormula <- formula(BaselineFormula)
  
  # Create data to predict upper and lower limits 
  Species1_Upper <- Species1 %>% filter(MeanTemp_CoralWatch_Scaled > mean(Species1$MeanTemp_CoralWatch_Scaled[which(Species1$Presence==1)]))
  Species1_Lower <- Species1 %>% filter(MeanTemp_CoralWatch_Scaled < mean(Species1$MeanTemp_CoralWatch_Scaled[which(Species1$Presence==1)]))
  
  #plot(Species1_Upper$Presence ~ Species1_Upper$MeanTemp_CoralWatch_Scaled)
  #plot(Species1_Lower$Presence ~ Species1_Lower$MeanTemp_CoralWatch_Scaled)
  
  # Fit the occupancy model. 
  OccupancyModel_Upper <- gam(BaselineFormula, data = Species1_Upper, family = binomial(link = "logit"))
  OccupancyModel_Lower <- gam(BaselineFormula, data = Species1_Lower, family = binomial(link = "logit"))
  
  #plot(OccupancyModel_Upper, pages = 1, all.terms = T)
  #plot(OccupancyModel_Lower, pages = 1, all.terms = T)
  
  # Predict from the occupancy model 
  PredictionFrame <- data.frame(SpeciesName = unique(Species1$SpeciesName),
                                MeanTemp_CoralWatch_Scaled = seq(from = min(Species1$MeanTemp_CoralWatch_Scaled)- 1, to = max(Species1$MeanTemp_CoralWatch_Scaled) + 1, length.out = 100), 
                                PC1 = mean(Species1$PC1), 
                                PC2 = mean(Species1$PC2), 
                                PC3 = mean(Species1$PC3),
                                PC4 = mean(Species1$PC4), 
                                #PC5 = mean(Species1$PC5), 
                                NEOLI = 0, 
                                Depth_Site_Scaled = mean(Species1$Depth_Site_Scaled))
  
  PredictionFrame$Preds_Upper <- boot::inv.logit(predict(OccupancyModel_Upper, PredictionFrame))
  PredictionFrame$Preds_Lower <- boot::inv.logit(predict(OccupancyModel_Lower, PredictionFrame))
  
  # Extract the values from T_Upper and T_Lower based on this occupancy model. 
  
  # Handle outputs 
  PredictionFrame$MeanTemp_CoralWatch <- (PredictionFrame$MeanTemp_CoralWatch_Scaled *  Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  #PredictionFrame$Topt_Occupancy <- PredictionFrame$MeanTemp_CoralWatch[which.max(PredictionFrame$Preds)]
  
  # Create T_Upper and T_Lower depending on the shape of the curve and treat as NA when there is no decline at one thermal limit, for that thermal limit. 
  #if(PredictionFrame$Topt_Occupancy[1] > max(Species1$MeanTemp_CoralWatch)){PredictionFrame$T_Upper <- NA}else{
    PredictionFrame$T_Upper <- PredictionFrame %>% filter(MeanTemp_CoralWatch > mean(Species1$MeanTemp_CoralWatch[which(Species1$Presence==1)])) %>%  .[which.min(abs(.$Preds_Upper-0.01)), ] %>% .$MeanTemp_CoralWatch # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
  #}
  #if(PredictionFrame$Topt_Occupancy[1] < min(Species1$MeanTemp_CoralWatch)){PredictionFrame$T_Lower <- NA}else{
    PredictionFrame$T_Lower <- PredictionFrame %>% filter(MeanTemp_CoralWatch < mean(Species1$MeanTemp_CoralWatch[which(Species1$Presence==1)])) %>%  .[which.min(abs(.$Preds_Lower-0.01)), ] %>% .$MeanTemp_CoralWatch # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
  #}

  # At the end need to make sure I return all the same objects (and more) to build the ThermalNicheData_New object for the
  # rest of the scripts and analysis/quality controls.

  # Create plot and save in a folder output. 
  PlotOccupancyModel <- 
    ggplot() + 
    geom_line(data = PredictionFrame, aes(y = Preds_Lower, x = MeanTemp_CoralWatch), col = 'blue3') + 
    geom_line(data = PredictionFrame, aes(y = Preds_Upper, x = MeanTemp_CoralWatch), col = 'red3') + 
    geom_point(data = Species1, aes(y = Presence, x = MeanTemp_CoralWatch)) + 
    geom_point(data = PredictionFrame, aes(x = T_Lower, y = 0.01), col = 'blue3', size = 3) + 
    geom_point(data = PredictionFrame, aes(x = T_Upper, y = 0.01), col = 'red3', size = 3)
  dir.create(paste0('figures_extra/', FOLDER))
  
  png(file = paste0('figures_extra/', FOLDER,'/', unique(Species1$SpeciesName), '.png'), res = 200, height = 1000, width = 1000)
  print(PlotOccupancyModel)
  dev.off()
  
  # Return the object with parameters estimated from the model above. 
  return(as_data_frame(PredictionFrame %>% dplyr::select(SpeciesName, T_Lower, T_Upper) %>% unique()))
  
}

OccupancyModels_All <- RLS_All %>% group_by(SpeciesName) %>% do(OccupancyModels = SpeciesSpecific_OccupancyModel(., FOLDER = 'TESTMODELS2'))
OccupancyModels_Combined <- do.call(rbind, OccupancyModels_All$OccupancyModels)
hist(OccupancyModels_Combined$Topt_Occupancy)
hist(OccupancyModels_Combined$T_Lower)
hist(OccupancyModels_Combined$T_Upper)

# Create data for occupancy models ----

RLS_Occ_Tupper <- RLS_20 %>% filter(Confidence_Occ_Tupper == 1); length(unique(RLS_Occ_Tupper$SpeciesName)) # n = 447
RLS_Occ_Tlower <- RLS_20 %>% filter(Confidence_Occ_Tlower == 1); length(unique(RLS_Occ_Tlower$SpeciesName)) # n = 637
RLS_20 %>% filter(Confidence_Occ_Tlower == 1 & Confidence_Occ_Tupper == 1) %>% .$SpeciesName %>% unique %>% length # n = 406

# Take a test species to test function
RLS_Occ_Conf <- RLS_20 %>% filter(N_Absences_T_Lower > 50 & N_Absences_T_Upper > 50)
Species1_test <- RLS_20 %>% filter(Confidence_Occ_Tlower == 1 & Confidence_Occ_Tupper == 1, SpeciesName == unique(.$SpeciesName)[1])

SpeciesSpecific_OccupancyModel(Species1 = Species1_test, FOLDER = 'test')

# Run this function over all the high confidence species 
OccupancyModels_All <- RLS_Occ_Conf %>% group_by(SpeciesName) %>% do(OccupancyModels = SpeciesSpecific_OccupancyModel(., FOLDER = 'High-confidence_V2'))
OccupancyModels_Combined <- do.call(rbind, OccupancyModels_All$OccupancyModels)

hist(RLS_Occ_Tupper$N_Absences_T_Upper[RLS_Occ_Tupper$N_Absences_T_Upper > 50])
length(unique(RLS_Occ_Tupper$SpeciesName[RLS_Occ_Tupper$N_Absences_T_Upper > 20]))
hist(RLS_Occ_Tlower$N_Absences_T_Lower)
length(unique(RLS_Occ_Tlower$SpeciesName[RLS_Occ_Tlower$N_Absences_T_Lower > 20]))

# Save test csv for species distribution modelling. 
# Site by lat long code 
Site_Location  <- as_data_frame(read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv')) %>% dplyr::select(SiteCode, SiteLat, SiteLong) %>% unique() # nrow = 587,840
Species1_test <- RLS_20 %>% filter(Confidence_Occ_Tlower == 1 & Confidence_Occ_Tupper == 1, SpeciesName == unique(.$SpeciesName)[2])
OccupancyDataTest <- left_join(Species1_test %>% dplyr::select(SpeciesName, SiteCode, Presence) %>% rename(., Occurrence = Presence), Site_Location)
write.csv(OccupancyDataTest, file = 'data_derived2/OccupancyDataTest.csv')

# Join in occupancy models with observed limits ----

ThermalNicheData_Obs_Occ <- left_join(ThermalNicheData_Obs, OccupancyModels_Combined)

# Only plot species which have occupancy thermal peaks within sampled range. 
ThermalNicheData_Obs_Occ_HighQuality <- 
  ThermalNicheData_Obs_Occ[-which(ThermalNicheData_Obs_Occ$Topt_Occupancy < ThermalNicheData_Obs_Occ$T_Lower_Obs | 
                               ThermalNicheData_Obs_Occ$Topt_Occupancy > ThermalNicheData_Obs_Occ$T_Upper_Obs),]

# What is the relationship between new thermal niche limit estimates and old thermal niche limit estimates? 
ThermalNicheData_Old <- readRDS(file = 'data_derived2/ThermalNicheData_Old.rds')

# Convert to NA when breaks criteria
#UpperModelled <- ThermalNicheData_Obs_Occ_HighQuality[,c('SpeciesName', 'T_Upper_Obs', 'T_Upper')]

T_Lower_Data <- ThermalNicheData_Obs_Occ_HighQuality %>% filter(T_Lower > T_Lower_Obs)
T_Upper_Data <- ThermalNicheData_Obs_Occ_HighQuality %>% filter(T_Upper < T_Upper_Obs)

# ThermalNicheData_Obs_Occ_HighQuality$T_Upper[which(ThermalNicheData_Obs_Occ_HighQuality$T_Upper_Obs >= ThermalNicheData_Obs_Occ_HighQuality$T_Upper)] <- NA
# ThermalNicheData_Obs_Occ_HighQuality$T_Lower[which(ThermalNicheData_Obs_Occ_HighQuality$T_Lower > ThermalNicheData_Obs_Occ_HighQuality$T_Lower_Obs)] <- NA
ThermalNicheData_New

ThermalNicheData_Old_Lower <- ThermalNicheData_Old %>% filter(SpeciesName %in% T_Lower_Data$SpeciesName)
ThermalNicheData_Old_Upper <- ThermalNicheData_Old %>% filter(SpeciesName %in% T_Upper_Data$SpeciesName)

ggplot() + 
  geom_point(aes(ThermalNicheData_Old_Lower$T_Lower, T_Lower_Data$T_Lower), col = 'blue') +
  geom_point(aes(ThermalNicheData_Old_Upper$T_Upper, T_Upper_Data$T_Upper), col = 'red') +
  geom_abline()

# ----



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
# save(RLS_Temp_Occ_Tupper, RLS_Temp_Occ_Tlower, RLS_Trop_Occ_Tupper, RLS_Trop_Occ_Tlower, file = 'data_derived2/OccupancyModelData_UpperLowerConfidences_2018-09-26.RData')

# If not loaded already, read in data for occupancy models. 
# load(file = 'data/OccupancyModelData_UpperLowerConfidences_2018-09-26.RData')

# Data to run models for based on only samples above and below T_Mean_Obs. 
# This is not Topt because some Topts are close to sampling limits and so this threshold would provide poor data in occupancy models. 
RLS_Temp_Occ_Tupper2 <- RLS_Temp_Occ_Tupper %>% filter(MeanTemp_CoralWatch > T_Mean_Obs)
RLS_Temp_Occ_Tlower2 <- RLS_Temp_Occ_Tlower %>% filter(MeanTemp_CoralWatch < T_Mean_Obs)
RLS_Trop_Occ_Tupper2 <- RLS_Trop_Occ_Tupper %>% filter(MeanTemp_CoralWatch > T_Mean_Obs)
RLS_Trop_Occ_Tlower2 <- RLS_Trop_Occ_Tlower %>% filter(MeanTemp_CoralWatch < T_Mean_Obs)

# Create functions to handle outputs from occupancy models ----
# Extract global parameter estimates from RLS glmmTMB models looking at temperature
ExtractGlobalTrend <- function(Data, Model, FixedFormula, Error = 'Poisson'){
  
  if(sum(grepl('ThermalGuild', FixedFormula)) < 1){
    SeqTemp <- seq(from = min(Data$MeanTemp_CoralWatch_Scaled), 
                   to = max(Data$MeanTemp_CoralWatch_Scaled), 
                   length = 100)
    
    MyData <- expand.grid(MeanTemp_CoralWatch_Scaled = SeqTemp,
                          HumPop50km2_Scaled = mean(Data$HumPop50km2_Scaled), 
                          npp_mean_Scaled = mean(Data$npp_mean_Scaled),
                          ReefAreaIn15km = 0,
                          Depth_Site_Scaled = mean(Data$Depth_Site_Scaled, na.rm = T),
                          ThermalGuild = unique(Data$ThermalGuild))
  }else{
    
    Data$ThermalGuild <- as.factor(Data$ThermalGuild)
    
    MyData <- ddply(Data, 
                    .(ThermalGuild), 
                    summarize,
                    MeanTemp_CoralWatch_Scaled = seq(min(MeanTemp_CoralWatch_Scaled), max(MeanTemp_CoralWatch_Scaled), length = 100))
    
    MyData$HumPop50km2_Scaled <- mean(Data$HumPop50km2_Scaled)
    MyData$npp_mean_Scaled <- mean(Data$npp_mean_Scaled)
    MyData$ReefAreaIn15km = 0
    MyData$Depth_Site_Scaled <- mean(Data$Depth_Site_Scaled, na.rm = T)
    
    
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
                                   FixedFormula = formula(~MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +  I(MeanTemp_CoralWatch_Scaled^3) + HumPop50km2_Scaled + npp_mean_Scaled), 
                                   Error = 'Poisson', 
                                   RandomSlopes = c('MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)', 'I(MeanTemp_CoralWatch_Scaled^3)')){
  
  Data$ThermalGuild <- as.factor(Data$ThermalGuild)
  
  MyData_2 <- plyr::ddply(Data, 
                          .(SpeciesName, ThermalGuild), 
                          plyr::summarize,
                          MeanTemp_CoralWatch_Scaled = seq(min(MeanTemp_CoralWatch_Scaled), max(MeanTemp_CoralWatch_Scaled), length = 100))
  
  MyData_2$HumPop50km2_Scaled <- mean(Data$HumPop50km2_Scaled)
  MyData_2$npp_mean_Scaled <- mean(Data$npp_mean_Scaled)
  MyData_2$ReefAreaIn15km = 0
  MyData_2$Depth_Site_Scaled <- mean(Data$Depth_Site_Scaled,na.rm = T)
  
  
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
  
  RandomEffectPlotData$MeanTemp_CoralWatch <- (RandomEffectPlotData$MeanTemp_CoralWatch_Scaled *  Data$MeanTemp_CoralWatch_scaleBT[1]) + Data$MeanTemp_CoralWatch_centerBT[1]
  
  
  return(RandomEffectPlotData)
}


# Define occupancy model formulas ----
# Temperate Species occupancy model formula
Temp_Formula1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + 
                           (1 + MeanTemp_CoralWatch_Scaled| SpeciesName) + (1 | ECOregion/SiteCode))

Temp_Formula2 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled + 
                           (1 + MeanTemp_CoralWatch_Scaled  | SpeciesName) + (1 | ECOregion/SiteCode))

Temp_Formula3 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled + 
                           (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))

# Tropical Species occupancy model formula
Trop_Formula1 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + 
                           (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))

Trop_Formula2 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled + 
                           (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled| SpeciesName) + (1 | ECOregion/SiteCode))

Trop_Formula3 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled + 
                           (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------



# OCCUPANCY TEMPERATE SPECIES ----
# Fit occupancy models TEMPERATE ----

### Temperate occupancy model and compare AIC between model approaches. 

# Fit model to uppers occupancy limits
Temp_Model_Upper1 <- glmmTMB(Temp_Formula1, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper1, file = 'data_derived2/Temp_Model_Upper1_2018-09-26.rds')
#Temp_Model_Upper1 <- readRDS(file = 'data_derived2/Temp_Model_Upper1_2018-09-26.rds')
Temp_Model_Upper2 <- glmmTMB(Temp_Formula2, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper2, file = 'data_derived2/Temp_Model_Upper2_2018-09-26.rds')
#Temp_Model_Upper2 <- readRDS(file = 'data_derived2/Temp_Model_Upper2_2018-09-26.rds')
Temp_Model_Upper3 <- glmmTMB(Temp_Formula3, family = "binomial", data = RLS_Temp_Occ_Tupper2); saveRDS(Temp_Model_Upper3, file = 'data_derived2/Temp_Model_Upper3_2018-09-26.rds')
AIC(Temp_Model_Upper1, Temp_Model_Upper2, Temp_Model_Upper3)

# Fit model to lowers occupancy limits
Temp_Model_Lower1 <- glmmTMB(Temp_Formula1, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower1, file = 'data_derived2/Temp_Model_Lower1_2018-09-26.rds')
Temp_Model_Lower2 <- glmmTMB(Temp_Formula2, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower2, file = 'data_derived2/Temp_Model_Lower2_2018-09-26.rds')
Temp_Model_Lower3 <- glmmTMB(Temp_Formula3, family = "binomial", data = RLS_Temp_Occ_Tlower2); saveRDS(Temp_Model_Lower3, file = 'data_derived2/Temp_Model_Lower3_2018-09-26.rds')
AIC(Temp_Model_Lower1, Temp_Model_Lower2, Temp_Model_Lower3)

# Read in occupancy models structure 1 and 3 for plots ----
Temp_Model_Upper1 <- readRDS(file = 'data_derived2/Temp_Model_Upper1_2018-09-26.rds')
Temp_Model_Upper3 <- readRDS(file = 'data_derived2/Temp_Model_Upper3_2018-09-26.rds')
Temp_Model_Lower1 <- readRDS(file = 'data_derived2/Temp_Model_Lower1_2018-09-26.rds')
Temp_Model_Lower3 <- readRDS(file = 'data_derived2/Temp_Model_Lower3_2018-09-26.rds')

# Creates plots and predictions from models -----

# Model list is just 1 and 3 for extraction of linear trend (a confidence criteria) and non-linear trend with non-linear random effect. 
Temp_Model_Upper_List <- list(Temp_Model_Upper1, Temp_Model_Upper3) # Create model list for functions to extract predictions
Temp_Model_Lower_List <- list(Temp_Model_Lower1, Temp_Model_Lower3) # Create model list for functions to extract predictions 

# Create list of objects for functions to plot outputs. 
# Extraction of global fixed effects. 
FixedFormula_Temp_List <- list(formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled), 
                               formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled),
                               formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled))

# Define random effects for input to function. 
RandomSlopes_Temp_List <- list(c('MeanTemp_CoralWatch_Scaled'), 
                               c('MeanTemp_CoralWatch_Scaled'), 
                               c('MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'))

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
  #geom_line(data = Temp_Model_Upper_Data, aes(x = MeanTemp_CoralWatch_Scaled, y = mu)) + 
  geom_line(data = Temp_Model_Upper_RE, aes(x = MeanTemp_CoralWatch_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Temp_Model_Upper_Data, aes(label = Model, x = 0.5, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Plot lowers
Temp_LowerPlot <- ggplot() + 
  #geom_line(data = Temp_Model_Lower_Data, aes(x = MeanTemp_CoralWatch_Scaled, y = mu)) + 
  geom_line(data = Temp_Model_Lower_RE, aes(x = MeanTemp_CoralWatch_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
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
Temp_M1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                # Remove polynomial temperature
Temp_M2 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled  + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                         # Remove all temperature
Temp_M3 <- formula(Presence ~ HumPop50km2_Scaled  + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode)) # Remove npp
Temp_M4 <- formula(Presence ~ npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))     # Remove human pop

# Single term models. 
Temp_M5 <- formula(Presence ~ MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode)) # Temperature only
Temp_M6 <- formula(Presence ~ HumPop50km2_Scaled + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                      # Human only
Temp_M7 <- formula(Presence ~ npp_mean_Scaled + (1 + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)  | SpeciesName) + (1 | ECOregion/SiteCode))                                         # npp only. 

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

save(Temp_Models_Upper_MS, Temp_Models_Lower_MS, file = 'data_derived2/Temperate-occupancy-model-selection.rdata')
  
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

save(FinalTempUpperModel, FinalTempLowerModel, file = 'data_derived2/Temperate-occupancy-final-models.rdata')

# load in model selections and final models from above code ----
load(file = 'data_derived2/Temperate-occupancy-model-selection.rdata')
load(file = 'data_derived2/Temperate-occupancy-final-models.rdata')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------



# OCCUPANCY TROPICAL SPECIES ---- 
# Fit occupancy models TROPICAL -----
### Tropical occupancy model and compare AIC between model approaches. 

# Fit model to uppers occupancy limits
Trop_Model_Upper1 <- glmmTMB(Trop_Formula1, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper1, file = 'data_derived2/Trop_Model_Upper1_2018-09-26.rds')
Trop_Model_Upper2 <- glmmTMB(Trop_Formula2, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper2, file = 'data_derived2/Trop_Model_Upper2_2018-09-26.rds')
Trop_Model_Upper3 <- glmmTMB(Trop_Formula3, family = "binomial", data = RLS_Trop_Occ_Tupper2); saveRDS(Trop_Model_Upper3, file = 'data_derived2/Trop_Model_Upper3_2018-09-26.rds')
AIC(Trop_Model_Upper1, Trop_Model_Upper2, Trop_Model_Upper3)

# Fit model to lowers occupancy limits
Trop_Model_Lower1 <- glmmTMB(Trop_Formula1, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower1, file = 'data_derived2/Trop_Model_Lower1_2018-09-26.rds')
Trop_Model_Lower2 <- glmmTMB(Trop_Formula2, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower2, file = 'data_derived2/Trop_Model_Lower2_2018-09-26.rds')
Trop_Model_Lower3 <- glmmTMB(Trop_Formula3, family = "binomial", data = RLS_Trop_Occ_Tlower2); saveRDS(Trop_Model_Lower3, file = 'data_derived2/Trop_Model_Lower3_2018-09-26.rds')
AIC(Trop_Model_Lower1, Trop_Model_Lower2, Trop_Model_Lower3)

# Read in occupancy models structure 1 and 3 for plots ----
Trop_Model_Upper1 <- readRDS(file = 'data_derived2/Trop_Model_Upper1_2018-09-26.rds')
Trop_Model_Upper3 <- readRDS(file = 'data_derived2/Trop_Model_Upper3_2018-09-26.rds')
Trop_Model_Lower1 <- readRDS(file = 'data_derived2/Trop_Model_Lower1_2018-09-26.rds')
Trop_Model_Lower3 <- readRDS(file = 'data_derived2/Trop_Model_Lower3_2018-09-26.rds')

# Creates plots and predictions from models ----

# Model list is just 1 and 3 for extraction of linear trend (a confidence criteria) and non-linear trend with non-linear random effect. 
Trop_Model_Upper_List <- list(Trop_Model_Upper1, Trop_Model_Upper3) # Create model list for functions to extract predictions
Trop_Model_Lower_List <- list(Trop_Model_Lower1, Trop_Model_Lower3) # Create model list for functions to extract predictions 

# Create list of objects for functions to plot outputs. 
# Extraction of global fixed effects. 
FixedFormula_Trop_List <- list(formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled), 
                               formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled),
                               formula(~ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled  + I(MeanTemp_CoralWatch_Scaled^2) + Depth_Site_Scaled))

# Define random effects for input to function. 
RandomSlopes_Trop_List <- list(c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled'), 
                               c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled'), 
                               c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'))


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
  #geom_line(data = Trop_Model_Upper_Data, aes(x = MeanTemp_CoralWatch_Scaled, y = mu)) + 
  geom_line(data = Trop_Model_Upper_RE, aes(x = MeanTemp_CoralWatch_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
  facet_wrap(~Model) + 
  geom_text(data = Trop_Model_Upper_Data, aes(label = Model, x = 0.5, y = 0.9), size = 3) + 
  theme(strip.text = element_blank(), aspect.ratio = 1)

# Plot lowers
Trop_LowerPlot <- ggplot() + 
  #geom_line(data = Trop_Model_Lower_Data, aes(x = MeanTemp_CoralWatch_Scaled, y = mu)) + 
  geom_line(data = Trop_Model_Lower_RE, aes(x = MeanTemp_CoralWatch_Scaled, y = mu, group = SpeciesName), alpha = 0.1) + 
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
Trop_M1 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                           (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))   # Remove reef
Trop_M2 <- formula(Presence ~ ReefAreaIn15km  + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans
Trop_M3 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove NPP
Trop_M4 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove temperature

# Remove two covariates
Trop_M5 <- formula(Presence ~ npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and humans 
Trop_M6 <- formula(Presence ~ HumPop50km2_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and npp 
Trop_M7 <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove reef and temperature
Trop_M8 <- formula(Presence ~ ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans and npp
Trop_M9 <- formula(Presence ~ ReefAreaIn15km  + npp_mean_Scaled +
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove humans and temperature
Trop_M10 <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + 
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Remove npp and temperature

# Single covariates 
Trop_M11 <- formula(Presence ~ HumPop50km2_Scaled + 
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only humans
Trop_M12 <- formula(Presence ~ npp_mean_Scaled + 
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only NPP
Trop_M13 <- formula(Presence ~ ReefAreaIn15km + 
                      (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))        # Only Reef
Trop_M14 <- formula(Presence ~ MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2) + 
                     (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)| SpeciesName) + (1 | ECOregion/SiteCode))         # Only temperature

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
save(Trop_Models_Upper_MS, Trop_Models_Lower_MS, file = 'data_derived2/Tropical-occupancy-model-selection.rdata')

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
save(Trop_Models_Upper_MS_V2, Trop_Models_Lower_MS_V2, file = 'data_derived2/Tropical-occupancy-model-selection_V2.rdata')

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

save(FinalTropUpperModel, FinalTropLowerModel, file = 'data_derived2/Tropical-occupancy-final-models.rdata')

# load in model selections and final models from above code ----
# Load these to create the models to select between. But then remove as will clog up memory. 
load(file = 'data_derived2/Tropical-occupancy-model-selection.rdata')
load(file = 'data_derived2/Tropical-occupancy-model-selection_V2.rdata')
load(file = 'data_derived2/Tropical-occupancy-final-models.rdata')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------




# OCCUPANCY QUALITY CONTROLS ----
# Fit models with same covariates as AIC selected formula apart from fitting with linear slopes to act as a quality controls ----

# Temperate upper
formula(Temp_Model_Upper3)
TempUpperModel_Linear <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + (1 + MeanTemp_CoralWatch_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTempUpperModel_Linear <- glmmTMB(TempUpperModel_Linear, family = "binomial", data = RLS_Temp_Occ_Tupper2)

# Temperate lower
formula(Temp_Model_Lower3)
TempLowerModel_Linear <- formula(Presence ~ HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + (1 + MeanTemp_CoralWatch_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTempLowerModel_Linear <- glmmTMB(TempLowerModel_Linear, family = "binomial", data = RLS_Temp_Occ_Tlower2)

# Tropical upper
formula(Trop_Model_Upper3)
TropUpperModel_Linear <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTropUpperModel_Linear <- glmmTMB(TropUpperModel_Linear, family = "binomial", data = RLS_Trop_Occ_Tupper2)

# Tropical lower
# Extract final model quality control linear slope. 
formula(Trop_Model_Lower3)
TropLowerModel_Linear <- formula(Presence ~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + Depth_Site_Scaled + (1 + ReefAreaIn15km + MeanTemp_CoralWatch_Scaled | SpeciesName) + (1 | ECOregion/SiteCode))
FinalTropLowerModel_Linear <- glmmTMB(TropLowerModel_Linear, family = "binomial", data = RLS_Trop_Occ_Tlower2)

# Save these
save(FinalTempUpperModel_Linear, FinalTempLowerModel_Linear, FinalTropUpperModel_Linear, FinalTropLowerModel_Linear, file = 'data_derived2/occupancy-models-linear-slopes.rdata')

# load in data from code above ----
load(file = 'data_derived2/occupancy-models-linear-slopes.rdata')
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
summary(FinalTempLowerModel)
summary(FinalTempUpperModel)
summary(FinalTropLowerModel)
summary(FinalTropUpperModel)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------




# EXTRACT NICHE LIMITS FROM RANDOM EFFECTS AND PLOT IN MULTIPANEL ----
# Function to extract niche limits from models ----
ExtractPredictions <- function(Data, Model, ModelNames, RandomSlopes, FixedFormula, Limit, 
                               i = i # This delimits the species to subset by and the random effect
){ 
  
  Data_all <- Data
  
  if(Limit == 'Upper'){
    
    Data <- Data %>% filter(SpeciesName == unique(.$SpeciesName)[i])
    Data$T_Mean_Obs_Scaled <- (Data$T_Mean_Obs - Data$MeanTemp_CoralWatch_centerBT) / Data$MeanTemp_CoralWatch_scaleBT
    Data$T_Upper_Absences_Scaled <- (Data$T_Upper_Absences+5 - Data$MeanTemp_CoralWatch_centerBT) / Data$MeanTemp_CoralWatch_scaleBT
    
    # Extract prediction frame from the data. 
    Preds <- Data %>% 
      select(SpeciesName, ThermalGuild, MeanTemp_CoralWatch_Scaled, T_Mean_Obs_Scaled, T_Upper_Absences_Scaled) %>% 
      group_by(SpeciesName, ThermalGuild) %>% 
      nest() %>% 
      mutate(MeanTemp_CoralWatch_Scaled = purrr::map(data, ~seq(unique(.$T_Mean_Obs_Scaled), unique(.$T_Upper_Absences_Scaled), length.out = 1000))) %>% 
      unnest(MeanTemp_CoralWatch_Scaled)
    
  }else{
    
    Data <- Data %>% filter(SpeciesName == unique(.$SpeciesName)[i])
    Data$T_Mean_Obs_Scaled <- (Data$T_Mean_Obs - Data$MeanTemp_CoralWatch_centerBT) / Data$MeanTemp_CoralWatch_scaleBT
    Data$T_Lower_Absences_Scaled <- (Data$T_Lower_Absences-5 - Data$MeanTemp_CoralWatch_centerBT) / Data$MeanTemp_CoralWatch_scaleBT
    
    # Extract prediction frame from the data. 
    Preds <- Data %>% 
      select(SpeciesName, ThermalGuild, MeanTemp_CoralWatch_Scaled, T_Mean_Obs_Scaled, T_Lower_Absences_Scaled) %>% 
      group_by(SpeciesName, ThermalGuild) %>% 
      nest() %>% 
      mutate(MeanTemp_CoralWatch_Scaled = purrr::map(data, ~seq(unique(.$T_Lower_Absences_Scaled), unique(.$T_Mean_Obs_Scaled), length.out = 1000))) %>% 
      unnest(MeanTemp_CoralWatch_Scaled)
    
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
  PredEstimates$MeanTemp_CoralWatch <- (PredEstimates$MeanTemp_CoralWatch_Scaled *  Data$MeanTemp_CoralWatch_scaleBT[1]) + Data$MeanTemp_CoralWatch_centerBT[1]
  
  PredEstimates$Model <- ModelNames
  
  if(Limit == 'Upper'){
    
    ### Estimate Tupper now from this predicted relationship. 
    #Where species have too low confidence above and below sampling limits they are given an NA.
    #min_mu_scaled <- min(PredEstimates$mu_scaled)
    #above_min_mu <- PredEstimates[which(PredEstimates$mu_scaled > min_mu_scaled),]
    #below_min_mu <- PredEstimates[which(PredEstimates$mu_scaled < min_mu_scaled),]
    
    #if(nrow(above_min_mu) > 0 | nrow(below_min_mu) > 0){
    PredEstimates$T_Upper <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$MeanTemp_CoralWatch # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
    PredEstimates$T_Upper_mu_scaled <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$mu_scaled
    # T_SD_Upper <- abs(T_Upper - Data$T_Opt) / 2 # This is a better way of estimateing thermal niche. Assuming that 95% of data fall in 2SD of mean. So Tlower is 0.025. Tupper is 0.975
    #}else{
    # This is where things get tricky. 
    
  }else{NULL}
  if(Limit == 'Lower'){
    ### Estimate TLower now from this predicted relationship. 
    # Where species have too low confidence above and below sampling limits they are given an NA.
    PredEstimates$T_Lower <- PredEstimates %>%  .[which.min(abs(.$mu_scaled-0.01)), ] %>% .$MeanTemp_CoralWatch # This should be relative to the maximum occupancy rate per species rather than overall occupanccy rate
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
                                                            RandomSlopes = 'MeanTemp_CoralWatch_Scaled',
                                                            FixedFormula = formula(~HumPop50km2_Scaled + MeanTemp_CoralWatch_Scaled), 
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
                                                     RandomSlopes = c('MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'),
                                                     FixedFormula = formula(~HumPop50km2_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)), 
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
                                                            RandomSlopes = c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled), 
                                                     Limit = 'Upper',
                                                     i = i)
  
  TropUpper_Output[[i]] <- ExtractPredictions(       Data         = RLS_Trop_Occ_Tupper2,
                                                            Model        = FinalTropUpperModel,
                                                            ModelNames   = 'Trop_Upper_Occ',
                                                            RandomSlopes = c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)), 
                                                     Limit = 'Upper', 
                                                     i = i)
  
  
  print(i / length(unique(RLS_Trop_Occ_Tupper2$SpeciesName)) * 100 )
}


# Extract non-negative linear slopes. 
# These are used as the bases of confidence in this parameter. 
Ranef_linear_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropUpperModel_Linear)$cond$SpeciesName), 
                                LinearSlope = ranef(FinalTropUpperModel_Linear)$cond$SpeciesName$MeanTemp_CoralWatch_Scaled + fixef(FinalTropUpperModel_Linear)$cond['MeanTemp_CoralWatch_Scaled'])

Ranef_linear_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempUpperModel_Linear)$cond$SpeciesName), 
                                LinearSlope  = ranef(FinalTempUpperModel_Linear)$cond$SpeciesName$MeanTemp_CoralWatch_Scaled + fixef(FinalTempUpperModel_Linear)$cond['MeanTemp_CoralWatch_Scaled'])

Ranef_linear <- rbind(Ranef_linear_Trop, Ranef_linear_Temp)

Ranef_linear$ConfidenceLinearSlope_Upper <- ifelse(Ranef_linear$LinearSlope < -1, 1, 0)

ConfidenceLinearSlope_Upper_Species <- Ranef_linear$SpeciesName[which(Ranef_linear$ConfidenceLinearSlope_Upper == 1)]

# Extract quadratic slopes. 
# These are used to support the confidence in negative linear slopes. Values between -1 and 0 can be weakly negative and overestimate thermal niche.
# We use a cut off of 'unimodal' curves as less than -1. 
# We use a cut-off of 'u-shaped' curves as greater than 0. 
Ranef_quad_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropUpperModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTropUpperModel)$cond$SpeciesName$`I(MeanTemp_CoralWatch_Scaled^2)` + fixef(FinalTropUpperModel)$cond['I(MeanTemp_CoralWatch_Scaled^2)'])

Ranef_quad_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempUpperModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTempUpperModel)$cond$SpeciesName$`I(MeanTemp_CoralWatch_Scaled^2)` + fixef(FinalTempUpperModel)$cond['I(MeanTemp_CoralWatch_Scaled^2)'])

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
  plot(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, type = 'l', main = as.character(Confidence_OccUpperAll_Unimodal)[i], ylim = c(0,1), xlim = c(min(SpeciesFits$MeanTemp_CoralWatch),unique(SpeciesFits$T_Upper)+1))
  points(Presence ~ MeanTemp_CoralWatch, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Upper), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits, col = 'dark orange')
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, col = 'Black')
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
  plot(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, type = 'l', main = unique(SpeciesFits$SpeciesName), ylim = c(0,1), xlim = c(min(SpeciesFits$MeanTemp_CoralWatch),unique(SpeciesFits$T_Upper)+1))
  points(Presence ~ MeanTemp_CoralWatch, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Upper), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits, col = 'dark orange', lty = 2)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, col = 'Black')
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
                                                            RandomSlopes = 'MeanTemp_CoralWatch_Scaled',
                                                            FixedFormula = formula(~HumPop50km2_Scaled + MeanTemp_CoralWatch_Scaled), 
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
                                                            RandomSlopes = c('MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'),
                                                            FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)), 
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
                                                     RandomSlopes = c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled'),
                                                     FixedFormula = formula(~HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled), 
                                                     Limit = 'Lower', 
                                                     i = i)
  
  TropLower_Output[[i]] <- ExtractPredictions(       Data         = RLS_Trop_Occ_Tlower2,
                                                     Model        = FinalTropLowerModel,
                                                     ModelNames   = 'Trop_Lower_Occ',
                                                     RandomSlopes = c('ReefAreaIn15km', 'MeanTemp_CoralWatch_Scaled', 'I(MeanTemp_CoralWatch_Scaled^2)'),
                                                     FixedFormula = formula(~ ReefAreaIn15km + HumPop50km2_Scaled + npp_mean_Scaled + MeanTemp_CoralWatch_Scaled + I(MeanTemp_CoralWatch_Scaled^2)), 
                                                     Limit = 'Lower', 
                                                     i = i)
  
  print(i / length(unique(RLS_Trop_Occ_Tlower2$SpeciesName)) * 100 )
}

# Extract non-negative linear slopes. 
# These are used as the bases of confidence in this parameter. 
Ranef_linear_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropLowerModel_Linear)$cond$SpeciesName), 
                                LinearSlope = ranef(FinalTropLowerModel_Linear)$cond$SpeciesName$MeanTemp_CoralWatch_Scaled + fixef(FinalTropLowerModel_Linear)$cond['MeanTemp_CoralWatch_Scaled'])

Ranef_linear_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempLowerModel_Linear)$cond$SpeciesName), 
                                LinearSlope  = ranef(FinalTempLowerModel_Linear)$cond$SpeciesName$MeanTemp_CoralWatch_Scaled + fixef(FinalTempLowerModel_Linear)$cond['MeanTemp_CoralWatch_Scaled'])

Ranef_linear <- rbind(Ranef_linear_Trop, Ranef_linear_Temp)

Ranef_linear$ConfidenceLinearSlope_Lower <- ifelse(Ranef_linear$LinearSlope > 1, 1, 0)

ConfidenceLinearSlope_Lower_Species <- Ranef_linear$SpeciesName[which(Ranef_linear$ConfidenceLinearSlope_Lower == 1)]

# Extract quadratic slopes. 
# These are used to support the confidence in negative linear slopes. Values between -1 and 0 can be weakly negative and overestimate thermal niche.
# We use a cut off of 'unimodal' curves as less than -1. 
# We use a cut-off of 'u-shaped' curves as greater than 0. 
Ranef_quad_Trop <- data.frame(SpeciesName = rownames(ranef(FinalTropLowerModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTropLowerModel)$cond$SpeciesName$`I(MeanTemp_CoralWatch_Scaled^2)` + fixef(FinalTropLowerModel)$cond['I(MeanTemp_CoralWatch_Scaled^2)'])

Ranef_quad_Temp <- data.frame(SpeciesName = rownames(ranef(FinalTempLowerModel)$cond$SpeciesName), 
                              QuadSlope  = ranef(FinalTempLowerModel)$cond$SpeciesName$`I(MeanTemp_CoralWatch_Scaled^2)` + fixef(FinalTempLowerModel)$cond['I(MeanTemp_CoralWatch_Scaled^2)'])

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
  plot(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, type = 'l', main = Confidence_OccLowerAll_Unimodal[i], ylim = c(0,1), xlim = c(unique(SpeciesFits$T_Lower)-1, max(SpeciesFits$MeanTemp_CoralWatch)))
  points(Presence ~ MeanTemp_CoralWatch, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Lower), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits, col = 'dark orange')
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, col = 'Black')
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
  plot(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, type = 'l', main = Confidence_OccLowerAll_U_Shaped[i], ylim = c(0,1), xlim = c(unique(SpeciesFits$T_Lower)-1, max(SpeciesFits$MeanTemp_CoralWatch)))
  points(Presence ~ MeanTemp_CoralWatch, data = Species_AllData, ylim = c(0,1))
  points(unique(SpeciesFits$T_Lower), 0.0, cex = 2, col = 'orange', pch = 19)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits, col = 'dark orange', lty = 2)
  lines(mu_scaled ~ MeanTemp_CoralWatch, data = SpeciesFits_Linear, col = 'Black')
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
AbsenceLimits <- RLS_All %>% select(SpeciesName, MeanTemp_CoralWatch, Presence) %>% group_by(SpeciesName) %>% 
  do(Absence_TUpper = max(.$MeanTemp_CoralWatch), 
     Absence_TLower = min(.$MeanTemp_CoralWatch)) %>% 
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
#save.image(file = 'data_derived2/objects-from-script-2.RData')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
