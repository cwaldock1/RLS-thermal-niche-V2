# Load libraries for spatial analysis ----
library(sp)
library(raster)
library(rgdal)
library(sdm)
library(rgeos)
library(RStoolbox)
library(tryCatchLog)
library(sdmvspecies)
library(doParallel)
library(filesstrings)

# SET UP ENVIRONMENTAL DATA ----
# Get rasters
path  <- '/Volumes/Untitled/Raster_files/BioOricleV2/'
lst   <- list.files(path=path,pattern='asc$',full.names = T)
path2 <- '/Volumes/Untitled/Raster_files/MSECData_Yeager2017'
lst2  <- list.files(path=path2,pattern='nc$',full.names = T)[c(3)]#,4,8,12,13,14)]
lst3  <- list.files(path=path2,pattern='nc$',full.names = T)[c(5,8,12,13,14)]

path3 <- '/Volumes/Untitled/Raster_files/GEBCO-30Sec'
lst4  <- list.files(path=path3,pattern='nc$',full.names = T)

# World ocean mask. 
install.packages('spaMM')
library(spaMM)
data("oceanmask")
plot(oceanmask)

preds_BIOORCLE  <- stack(lst)
preds_MSEChuman <- stack(lst2)
preds_MSEChuman <- preds_MSEChuman$X2015
preds_MSEC      <- stack(lst3)
preds_MSEChuman <- rotate(preds_MSEChuman)
preds_MSEC      <- rotate(preds_MSEC)
extent(preds_MSEC)       <- extent(preds_BIOORCLE)
extent(preds_MSEChuman)  <- extent(preds_BIOORCLE)
crs(preds_BIOORCLE)      <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
crs(preds_MSEC)          <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
crs(preds_MSEChuman)     <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

# Handle depth data
preds_DEPTH     <- raster(lst4)
preds_DEPTH <- mask(preds_DEPTH, oceanmask)
preds_DEPTH_reprojected <- projectRaster(preds_DEPTH, preds_BIOORCLE[[1]])
preds_DEPTH_reprojected[preds_DEPTH_reprojected > 100] <- NA
plot(preds_DEPTH_reprojected)


# Stack all the environmental data in one place for use in loop. 
preds_MSEC_Crop_reprojected      <- projectRaster(preds_MSEC, preds_BIOORCLE[[1]])
preds_MSEChuman_Crop_reprojected <- projectRaster(preds_MSEChuman, preds_BIOORCLE[[1]])

# Stack the re-scaled rasters together. 
preds_Stack1 <- stack(preds_MSEC_Crop_reprojected, preds_BIOORCLE)
preds_Stack2 <- stack(preds_MSEChuman_Crop_reprojected, preds_Stack1)
preds_Stack3 <- stack(preds_DEPTH_reprojected, preds_Stack2)
plot(preds_Stack3)

# Drop correlated layers
preds_Stack3 <- dropLayer(preds_Stack3, c(6, 9, 11, 16))
names(preds_Stack3) <- c('Depth', 'HumanPop', 'LandArea', 'NPP', 'ReefArea', 'WaveEnergy', 'CurrentVelocity',
                         'O2', 'Iron', 'Nitrate', 'pH', 'Phosphate', 'Salinity', 'Silicate', 'Temperature')

hist(preds_Stack3$Depth)
hist(log10(preds_Stack3$HumanPop+1))
hist(log10(preds_Stack3$LandArea + 1))
hist(log10(preds_Stack3$NPP + 1))
hist(log10(preds_Stack3$ReefArea + 1))
hist(log10(preds_Stack3$WaveEnergy + 1))
hist(log10(preds_Stack3$CurrentVelocity+1))
hist(preds_Stack3$O2)
hist(log10(preds_Stack3$Iron+0.000001))
hist(preds_Stack3$Nitrate)
hist(preds_Stack3$pH)
hist(preds_Stack3$Phosphate)
hist(preds_Stack3$Salinity)
hist(preds_Stack3$Silicate)
hist(preds_Stack3$Temperature)

# Apply transformations to data layers 
?calc
preds_Stack3$HumanPop <- calc(preds_Stack3$HumanPop, function(x) log10(x+1))
preds_Stack3$LandArea <- calc(preds_Stack3$LandArea, function(x) log10(x+1))
preds_Stack3$NPP <- calc(preds_Stack3$NPP, function(x) log10(x+1))
preds_Stack3$ReefArea <- calc(preds_Stack3$ReefArea, function(x) log10(x+1))
preds_Stack3$WaveEnergy <- calc(preds_Stack3$WaveEnergy, function(x) log10(x+1))
preds_Stack3$CurrentVelocity <- calc(preds_Stack3$CurrentVelocity, function(x) log10(x+1))
preds_Stack3$Iron <- calc(preds_Stack3$Iron, function(x) log10(x+0.000001))

#saveRDS(preds_Stack3, file = '/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack.nc')
preds_Stack3 <- readRDS(file = '/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack.nc')


# SET UP OCCUPANCY DATA ----
OccDat <- read.csv('data_derived2/OccupancyDataTest.csv')
OccDat <- OccDat[sample(1:nrow(OccDat), size = 100),]
Range <- SpatialPointsDataFrame(coords = cbind(OccDat$SiteLong, OccDat$SiteLat), 
                                data = data.frame(OccDat %>% dplyr::select(Occurrence)))
plot(Range)
plot(Range[Range$Occurrence == 0,],col='red',pch=16)
points(Range[Range$Occurrence == 1,],col='blue',pch=16) 

# This needs gBuffer
crs(Range) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
Range2 <- gBuffer(Range, width = 10)
plot(Range2)

# CROP ENVIRONMENT BY OCCUPANCY POLYGONS ---- 
# Crop all covariates by range.
Ext <- extent(bbox(Range2)[1,1]-5,bbox(Range2)[1,2]+5,bbox(Range2)[2,1]-5,bbox(Range2)[2,2]+5)
preds_Crop  <- crop(mask(preds_Stack3, Range2), Ext)

# Apply reprojections
#preds2_Crop_reprojected <- projectRaster(preds2_Crop, preds_Crop[[1]])

#preds3 <- stack(preds_Crop, preds2_Crop_reprojected)

# Crop to species range.
# Ext <- extent(bbox(Range)[1,1]-5,bbox(Range)[1,2]+5,bbox(Range)[2,1]-5,bbox(Range)[2,2]+5)
# preds_Crop <- crop(mask(preds3, Range2), Ext)

plot(preds_Crop[[1]])
points(Range)

# FIT SPECIES' DISTRIBUTION MODELS ----

d <- sdmData(formula = Occurrence~., train = Range, predictors = preds_Crop)

m1 <- sdm(Occurrence~., data=d, methods=c('brt', 'gam', 'svm', 'rf'), 
          replication='cv',cv.folds=5)

boxplot(m1)

getEvaluation(m1, stat = c('AUC','specificity','sensitivity','prevalence'))

rocPlot <- roc(m1,smooth=TRUE)

p1 <- ensemble(m1, newdata = preds_Crop, filename = "p1.img", 
               overwrite = T, mean = T, setting = list(method='weighted',stat='AUC'))

p1_predict <- predict(m1, newdata=preds_Crop, filename='preds.img')
p1_predict[p1_predict$id_1.sp_1.m_brt.re_cros]
EvaluateM1 <- evaluates(x = m1, p = p1_predict$id_1.sp_1.m_brt.re_cros[p1_predict$id_1.sp_1.m_brt.re_cros])


p1[p1<=0.1]<-NA
#p1[p1>=0.1]<-NA
plot(p1)

# Upper and lower. 
quantile(preds_Crop$Present.Surface.Temperature.Mean[p1], 0.01, na.rm = T)
quantile(preds_Crop$Present.Surface.Temperature.Mean[p1], 0.99, na.rm = T)

hist(preds_Crop$Present.Surface.Temperature.Mean[p1])
hist(Species1_test$MeanTemp_CoralWatch)
hist(Species1_test$MeanTemp_CoralWatch[Species1_test$Presence==1])


# Make function that works over a raster stack of all environmental data ----
OccurrenceData    <- OccDat
EnvironmentalData <- preds_Stack3
FolderName        <- 'TestSDM'

# Function fits SDMs, and saves RDS file on hard-drive. 
ExtractSDMs <- function(OccurrenceData, EnvironmentalData){ 
  
  # Convert lats and longs into spatial points data frame. 
  Range <- SpatialPointsDataFrame(coords = cbind(OccurrenceData$SiteLong, OccurrenceData$SiteLat), 
                                  data = data.frame(OccurrenceData %>% dplyr::select(Occurrence)))
  
  # Assign spatial projection 
  crs(Range) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
  
  # Buffer 10latlong away from sites. 
  Range2 <- gBuffer(Range, width = 5)

  # Crop the environmental layers to bounding box (and some) 
  Ext <- extent(bbox(Range)[1,1]-5,bbox(Range)[1,2]+5,bbox(Range)[2,1]-5,bbox(Range)[2,2]+5)
  print('Cropping environmental data to species range')
  EnvironmentalData_Cropped  <- crop(mask(EnvironmentalData, Range2), Ext)
  
  # Ensure there are no singluar raster layers
  Check_SD <- cellStats(EnvironmentalData_Cropped, "sd")
  if(sum(Check_SD == 0) != 0){
    #If a column is singluar it breaks the pc. Remove all singluar columns
    EnvironmentalData_Cropped <- dropLayer(EnvironmentalData_Cropped, which(Check_SD == 0))
  }else{NULL}

  # Perform PCA on selected raster stack variables using RStoolbox::rasterPCA
  print('Fitting PCA to environmental layers')
  EnvironmentalData_Cropped_PCA <- rasterPCA(dropLayer(EnvironmentalData_Cropped, c(1, 15)), nComp = 5, spca = T, maskCheck = TRUE)
  vars <- summary(EnvironmentalData_Cropped_PCA$model)$sdev^2 
  vars <- vars/sum(vars) 
  EnvironmentalData_Cropped_PCA2 <- dropLayer(EnvironmentalData_Cropped_PCA$map, which(vars < 0.1))

  # Stack PCA onto temperature and depth. 
  EnvironmentalData_Cropped_FINAL <- stack(EnvironmentalData_Cropped_PCA2, dropLayer(EnvironmentalData_Cropped, which(!1:15 %in% c(1, 15))))
  
  # Convert to a SDM model object 
  SDM_Data <- sdmData(formula = Occurrence~., train = Range, predictors = EnvironmentalData_Cropped_FINAL)
  
  # Create SDM model
  print('Fitting SDM')
  SDM_Model <- sdm(Occurrence~., data=SDM_Data, 
            methods=c('brt', 'gam', 'svm', 'rf'), 
            replication='cv',cv.folds=5) # These are the 4 most flexible models. 
  
  # Create_ROC plots
  #rocPlot <- roc(SDM_Model, smooth=TRUE)
  
  # Aggregate models 
  print('Creating model ensemble')
  SDM_Ensemble <- ensemble(SDM_Model, newdata = EnvironmentalData_Cropped_FINAL, 
                           #filename = 'test.img',
                 filename = paste0(
                   #'/Volumes/Untitled/Raster_files/RLS_SDMs/', 
                   OccurrenceData$SpeciesName[1], '.grd'), 
                 overwrite = T, 
                 mean = T, 
                 setting=list(method='weighted',stat='AUC'))
  
  # Copy the files across to an external device. 
  file.copy(from = list.files()[grep(OccurrenceData$SpeciesName[1], list.files())], 
            to = "/Volumes/Untitled/Raster_files/RLS_SDMs",
            recursive = FALSE,
            copy.mode = TRUE, 
            copy.date = FALSE)
  

  # Remove the created file to save disk memory
  #unlink(paste0('data_derived2' , '/', FolderName, '/', OccurrenceData$SpeciesName[1], '.img'), recursive = FALSE)
  
  # Create new raster to save
  #print('Saving rds file of in-memory raster')
  #SDM_Ensemble2 <- SDM_Ensemble
  #saveRDS(SDM_Ensemble2, paste0('/Volumes/Untitled/Raster_files/RLS_SDMs/', as.character(OccurrenceData$SpeciesName[1]), '.rds'))
  #SDM_Ensemble2 <- readRDS('/Volumes/Untitled/Raster_files/RLS_SDMs/Abudefduf sexfasciatus.rds')
  
  # Create new raster to aggregate from
  SDM_Ensemble <- sdmvspecies::rescale(SDM_Ensemble)
  T_Lower_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.01, na.rm = T), 3)
  T_Lower_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.01, na.rm = T), 3)
  T_Lower_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.01, na.rm = T), 3)
  T_Upper_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.99, na.rm = T), 3)
  T_Upper_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.99, na.rm = T), 3)
  T_Upper_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.99, na.rm = T), 3)
  
  ModelSummary <- getEvaluation(SDM_Model, stat = c('AUC', 'specificity', 'sensitivity', 'TSS'))
  ModelSummaries <- colMeans(ModelSummary)[2:5]
  ModelSummaries <- t(data.frame(ModelSummaries))

  # Delete the files from local device. 
  unlink(list.files()[grep(OccurrenceData$SpeciesName[1], list.files())]) # Was used to delete the created ensemble model. Now saving. 
  
  # Return summaries of temperature
  return(as_data_frame(cbind(SpeciesName = as.character(OccurrenceData$SpeciesName[1]),
                             T_Lower_0.1, T_Lower_0.25,T_Lower_0.5, 
                             T_Upper_0.1, T_Upper_0.25, T_Upper_0.5, 
                             round(ModelSummaries,3))))
  
}

# Function to extract metrics from an ensemble data layer
EnsembleModelList <- list.files('/Volumes/Untitled/Raster_files/RLS_SDMs', full.names = TRUE)
EnvironmentalData <- preds_Stack3
OccurrenceData <- HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[2])

# Finish later if needed... 
function(OccurrenceData, EnvironmentalData, EnsembleModel){
  
  SDM_Ensemble <- readRDS(EnsembleModelList[grep(OccurrenceData$SpeciesName[2], EnsembleModelList)])
  
SDM_Ensemble <- sdmvspecies::rescale(SDM_Ensemble)
T_Lower_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.01, na.rm = T), 3)
T_Lower_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.01, na.rm = T), 3)
T_Lower_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.01, na.rm = T), 3)
T_Upper_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.99, na.rm = T), 3)
T_Upper_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.99, na.rm = T), 3)
T_Upper_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.99, na.rm = T), 3)

#ModelSummary <- getEvaluation(SDM_Model, stat = c('AUC', 'specificity', 'sensitivity', 'TSS'))
#ModelSummaries <- colMeans(ModelSummary)[2:5]
#ModelSummaries <- t(data.frame(ModelSummaries))
T_Lower <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble], 0.01, na.rm = T), 3)
T_Upper <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble], 0.99, na.rm = T), 3)
unlink(list.files()[grep(OccurrenceData$SpeciesName[1], list.files())])

return(as_data_frame(cbind(SpeciesName = as.character(OccurrenceData$SpeciesName[1]),
                           T_Lower_0.1, T_Lower_0.25,T_Lower_0.5, 
                           T_Upper_0.1, T_Upper_0.25, T_Upper_0.5#, 
                           #round(ModelSummaries,3)
                           )))
}

OccurrenceData    <- OccDat
EnvironmentalData <- preds_Stack3

ExtractSDMs(OccDat[sample(1:nrow(OccDat), 100),], preds_Stack3)
beepr::beep()



# Run function over 50 high influence species from ThermalNicheData_Old ----

ThermalNicheData_Old <- readRDS(file = 'data_derived2/ThermalNicheData_Old.rds')
ThermalNicheData_Old <- ThermalNicheData_Old %>% filter(ConfidenceCombined == 3)
#ThermalNicheData_Old <- ThermalNicheData_Old %>% filter(SpeciesName %in% unique(.$SpeciesName)[1:50])

Site_Location  <- as_data_frame(read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv')) %>% dplyr::select(SiteCode, SiteLat, SiteLong) %>% unique() # nrow = 587,840

RLS_20 <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')
HighQualityOccupancies <- RLS_20 %>% filter(SpeciesName %in% ThermalNicheData_Old$SpeciesName)
HighQualityOccupancies <- left_join(HighQualityOccupancies %>% dplyr::select(SpeciesName, SiteCode, Presence) %>% rename(., Occurrence = Presence), Site_Location)


# Test on full environmental dataset
ExtractSDMs(OccurrenceData = HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[4]), 
            EnvironmentalData = preds_Stack3)
beep()
cl <- makeCluster(4)
registerDoParallel(cl)
sdmModels <- 
  foreach(i=1:length(unique(HighQualityOccupancies$SpeciesName)), 
          .packages=c('tidyr', 'dplyr', 'sdm', 'sp', 'raster', 'rgeos', 'sdmvspecies', 'tryCatchLog', 'RStoolbox', 'filesstrings')) %dopar% 
          {
            tryCatch(ExtractSDMs(OccurrenceData = HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[i]), 
                               EnvironmentalData = preds_Stack3), 
                     error = function(e) NA)
          }
#beepr::beep(4)
stopCluster(cl)

# Read in predicted rasters. 
sdmModels_AllCovariates <- do.call(rbind, sdmModels) # THIS IS RUN FROM LAST 30/09/2018 NIGHT 
saveRDS(sdmModels_AllCovariates, file = 'data_derived2/sdmModels_AllCovariates.rds')
sdmModels_AllCovariates[,-1] <- apply(sdmModels_AllCovariates[,-1], 2, as.numeric)

# Fit SDMs that are simpled in covariate structure to obtain Tupper and Tlower ----

# Test on smaller environmental dataset 
preds_Stack_Subset <- dropLayer(preds_Stack3, which(!1:15 %in% c(1, 2, 4, 5, 15)))
ExtractSDMs(OccurrenceData = HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[4]), 
            EnvironmentalData = preds_Stack_Subset)

beep()


cl <- makeCluster(4)
registerDoParallel(cl)
sdmModels3 <- 
  foreach(i=1:length(unique(HighQualityOccupancies$SpeciesName)), 
          .packages=c('tidyr', 'dplyr', 'sdm', 'sp', 'raster', 'rgeos', 'sdmvspecies', 'tryCatchLog')) %dopar% 
          {
            #tryCatch(
              tryLog(ExtractSDMs(OccurrenceData = HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[i]), 
                          EnvironmentalData = preds_Stack_Subset))
             # , error = function(e) NA)
          }
stopCluster(cl)
beep()

sdmModels_SimpleSDMs <- do.call(rbind, sdmModels3)
hist(as.numeric(sdmModels_SimpleSDMs$AUC))

# Plot up the SDM outputs for high confidence species' ----

# Combine these SEMs with the old thermal niche data
sdmModels_SimpleSDMs$T_Lower_NEW <- as.numeric(sdmModels_SimpleSDMs$T_Lower)
sdmModels_SimpleSDMs$T_Upper_NEW <- as.numeric(sdmModels_SimpleSDMs$T_Upper)

# Join in 
ThermalNicheData_Old2 <- left_join(ThermalNicheData_Old, sdmModels_SimpleSDMs[,c('SpeciesName', 'T_Lower_NEW', 'T_Upper_NEW', 'AUC')])

# Plot comparison of old niche limit with new niche limit. 
ggplot(ThermalNicheData_Old2) + 
  geom_point(aes(x = T_Lower, y = T_Lower_NEW), col = 'blue') + 
  geom_point(aes(x = T_Upper, y = T_Upper_NEW), col = 'red') + 
  geom_abline() + 
  xlab('Old realised thermal limit') + 
  ylab('New realised thermal limit')

Quantile_Parameters$Topt_NEW <- Quantile_Parameters$Topt
ThermalNicheData_Old2 <- left_join(ThermalNicheData_Old2, Quantile_Parameters %>% dplyr::select(SpeciesName, Topt_NEW))

# Plots and test models 
ggplot(ThermalNicheData_Old2) + 
  geom_point(aes(y = T_Skew, x = Topt), col = 'blue') + 
  geom_point(aes(y = (T_Upper_NEW-Topt_NEW) - (Topt_NEW - T_Lower_NEW), x = Topt_NEW), col = 'red') + 
  stat_smooth(aes(y = T_Skew, x = Topt), col = 'blue', method = 'lm') + 
  stat_smooth(aes(y = (T_Upper_NEW-Topt_NEW) - (Topt_NEW - T_Lower_NEW), x = Topt_NEW), col = 'red', method = 'lm') + 
  facet_wrap(~ThermalGuild, scales = 'free_x')

ThermalNicheData_Old2$T_Skew_NEW <- with(ThermalNicheData_Old2, (T_Upper_NEW-Topt_NEW) - (Topt_NEW - T_Lower_NEW))
T_Skew_NEW_m1 <- lme4::lmer(T_Skew_NEW ~ ThermalGuild * Topt_NEW + (1|Order/Family), data = ThermalNicheData_Old2) 
summary(T_Skew_NEW_m1)
T_Skew_NEW_m2 <- lme4::lmer(T_Skew_NEW ~ ThermalGuild + Topt_NEW + (1|Order/Family), data = ThermalNicheData_Old2) 
anova(T_Skew_NEW_m1, T_Skew_NEW_m2)

MuMIn::r.squaredGLMM(T_Skew_NEW_m1)

T_Skew_NEW_m2 <- lme4::lmer(T_Skew_NEW ~ Topt_NEW + (1|Order/Family), data = ThermalNicheData_Old2 %>% filter(ThermalGuild != 'Tropical')) 

pdf(file = 'figures_extra/Comparison-of-skews.pdf', width = 5, height = 10)
gridExtra::grid.arrange(ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_Obs >= Topt_NEW | T_Upper_Obs <= Topt_NEW)) + 
  geom_point(aes(y = T_Lower_NEW, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
  geom_point(aes(y = T_Upper_NEW, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
  geom_linerange(aes(ymax = T_Lower_NEW, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
  geom_linerange(aes(ymax = T_Upper_NEW, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
  geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
  theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
  ylab('SDM thermal niche limits'), 

ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
  geom_point(aes(y = T_Lower, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
  geom_point(aes(y = T_Upper, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
  geom_linerange(aes(ymax = T_Lower, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
  geom_linerange(aes(ymax = T_Upper, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
  geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
  theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
  ylab('Occupancy mixed model niche limit'),

ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
  geom_point(aes(y = T_Lower_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
  geom_point(aes(y = T_Upper_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
  geom_linerange(aes(ymax = T_Lower_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
  geom_linerange(aes(ymax = T_Upper_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
  geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
  ylab('Sampling thermal limit') + xlab('Topt') +
  theme_bw() + theme(aspect.ratio = 0.75), 
nrow = 3)
dev.off()

  
ThermalNicheData_Old2 %>% filter(T_Lower_NEW > Topt_NEW) # 5 cases where topt is < tmin 
ThermalNicheData_Old2 %>% filter(T_Upper_NEW < Topt_NEW) # 13 cases where topt is > tmin 

# Combine and compare across old and new ----
ThermalNicheData_Old2 <- left_join(ThermalNicheData_Old, sdmModels_AllCovariates)

ggplot(ThermalNicheData_Old2) + 
  geom_point(aes(x = T_Lower, y = T_Lower_0.25), col = 'blue') + 
  geom_point(aes(x = T_Upper, y = T_Upper_0.25), col = 'red')

# load
load(file = 'data_derived2/AllQgams_2018-09-28.RData')

# Whats the relation between Tskew and Topt in the new and old data
Quantile_Parameters$Topt_NEW <- Quantile_Parameters$Topt
ThermalNicheData_Old2 <- left_join(ThermalNicheData_Old2, Quantile_Parameters %>% dplyr::select(SpeciesName, Topt_NEW))

ggplot(ThermalNicheData_Old2) + 
  geom_point(aes(y = T_Skew, x = Topt), col = 'blue') + 
  geom_point(aes(y = (T_Upper_0.5-Topt_NEW) - (Topt_NEW - T_Lower_0.5), x = Topt_NEW), col = 'red') + 
  stat_smooth(aes(y = T_Skew, x = Topt), col = 'blue', method = 'lm') + 
  stat_smooth(aes(y = (T_Upper_0.5-Topt_NEW) - (Topt_NEW - T_Lower_0.5), x = Topt_NEW), col = 'red', method = 'lm')

# Whats the relation between observed limits and modelled occupancy thermal limits. 
ThermalNicheData_Obs <- RLS_20 %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanTemp_CoralWatch)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanTemp_CoralWatch)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanTemp_CoralWatch, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
  dplyr::select(-data)

ThermalNicheData_Obs2 <- left_join(ThermalNicheData_Obs, sdmModels_AllCovariates)
ThermalNicheData_Obs2
ggplot(ThermalNicheData_Obs2) + 
  geom_point(aes(y = as.numeric(T_Lower_0.5), x = T_Lower_Obs), col = 'blue') + 
  geom_point(aes(y = as.numeric(T_Upper_0.5), x = T_Upper_Obs), col = 'red') + 
  geom_abline()


pdf(file = 'figures_extra/Comparison-of-skews.pdf', width = 5, height = 10)
gridExtra::grid.arrange(ggplot(ThermalNicheData_Old2)+ #%>% filter(T_Lower_0.1 >= Topt_NEW | T_Upper_0.1 <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower_0.5, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_0.5, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_0.5, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_0.5, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
                          geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('SDM thermal niche limits'), 
                        
                        ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
                          geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('Occupancy mixed model niche limit'),
                        
                        ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
                          geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
                          ylab('Sampling thermal limit') + xlab('Topt') +
                          theme_bw() + theme(aspect.ratio = 0.75), 
                        nrow = 3)
dev.off()

# Same as above but for mid points. 
gridExtra::grid.arrange(ggplot(ThermalNicheData_Old2)+ #%>% filter(T_Lower_0.1 >= Topt_NEW | T_Upper_0.1 <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower_0.5, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_0.5, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_0.5, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_0.5, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'red') + 
                          geom_point(aes(y = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('SDM thermal niche limits'), 
                        
                        ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'red') + 
                          geom_point(aes(y = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('Occupancy mixed model niche limit'),
                        
                        ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower_Obs, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_Obs, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_Obs, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_Obs, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'red') + 
                          geom_point(aes(y = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'black', size = 0.5) + 
                          ylab('Sampling thermal limit') + xlab('Topt') +
                          theme_bw() + theme(aspect.ratio = 0.75), 
                        nrow = 3)


# Compare plots based on midpoints to plots based on topts. 

gridExtra::grid.arrange(ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
  geom_point(aes(y = T_Lower_Obs, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'blue', size = 2) + 
  geom_point(aes(y = T_Upper_Obs, x = T_Midpoint_Obs, pch = ThermalGuild), col = 'red', size = 2) + 
  geom_linerange(aes(ymax = T_Lower_Obs, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'blue') + 
  geom_linerange(aes(ymax = T_Upper_Obs, ymin = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'red') + 
  geom_point(aes(y = T_Midpoint_Obs, x = T_Midpoint_Obs), col = 'black', size = 0.5) + 
  ylab('Sampling thermal limit') + xlab('Topt') +
  theme_bw() + theme(aspect.ratio = 0.75), 
ggplot(ThermalNicheData_Old2) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
  geom_point(aes(y = T_Lower_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
  geom_point(aes(y = T_Upper_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
  geom_linerange(aes(ymax = T_Lower_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
  geom_linerange(aes(ymax = T_Upper_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
  geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
  ylab('Sampling thermal limit') + xlab('Topt') +
  theme_bw() + theme(aspect.ratio = 0.75), nrow = 2)


