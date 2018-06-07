# First script in thermal-niche RLS analysis. 
# Initiated 27/03/2018
# Author: Conor Waldock

# The following raw data are required prior to running these scripts. 
# RLS_GLOBAL_ABUND_BIOMASS_20170507.csv   # Acquiered from Antonia Cooper, UTAS, Taroon on 5th June 2017.  
# species_traits.csv                      # Acquired from Rick Stuart-Smith January 2017
# SiteTemperatureValuesRLS_04_09_2017.csv # Produced on 4th September in external script. 
# RLSWarmingRate_04_09_2017.csv           # Produced on 4th September in external script. 
# RLS_Site_Covariates_V2.csv              # Site level covariates acquired from Rick Stuart-Smith on 5th November 2017. 

###  ORGANISING DATA SECTION ----
### Load relevant packages for data manipulation ----
### 
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(geosphere)
library(rgeos)
library(data.table)
#library(qgam)
#library(glmmTMB)
#library(R2jags)
#library(gridExtra)
#library(MuMIn)
#library(remef) # Creation of partial residual plots
#source(file = "/Users/caw1g15/Dropbox/PhD Plans and Courses/Courses/Highland Statistics Mixed Modelling and GLMM/Highland Statistics Mixed Modelling and GLMM/MCMCSupportHighstatV4.R")
#source('/Users/caw1g15/Dropbox/PhD Plans and Courses/Courses/Highland Statistics Mixed Modelling and GLMM/Highland Statistics Mixed Modelling and GLMM/HighstatLibV10.R')
### 




### Read in raw RLS data and preliminary cleaning ---- 
RLS  <- as_data_frame(read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv')) # nrow = 587,840

# Remove Method 0 data. 
RLS <- RLS %>% filter(Method != 0) %>% filter(Block != 0) %>% filter(Method != 2) # nrow = 433,215: Remove method 0, method 2 and unknown blocks. 

# Check blocks are correct in structure. 
hist(RLS$Block)           # Blocks are equal in number now. 
length(unique(RLS$Block)) # Remove blocks which are 0.
unique(RLS$Method)        # Remove Method which are 0 or 2 

# Remove only vertebrates (i.e., fish)
RLS <- RLS %>% filter(PHYLUM == 'Chordata') # nrow = 432,432

# Total number of SurveyIDs 
length(unique(RLS$SurveyID))

# Total number of Sites
length(unique(RLS$SiteCode))

# Remove those surveys without 2 blocks. 
# Work out number of blocks. 
TestNBlocks <- RLS %>% 
  group_by(SiteCode, Site.name, SiteLat, SiteLong, SurveyDate, Depth, SurveyID, Method) %>% 
  do(NBlocks = length(unique(.$Block))) %>% 
  unnest(NBlocks) %>%
  ungroup() 

# Count number of blocks per survey. 
TestNBlocks %>% unique() %>% filter(NBlocks == 1) # 547
TestNBlocks %>% unique() %>% filter(NBlocks == 2) # 19,968

# Ensure number of blocks = 2 for sampling consistency. 
RLS2 <- left_join(TestNBlocks %>% filter(NBlocks > 1) %>% select(-NBlocks), RLS) # 430,191
RemovedData <- left_join(TestNBlocks %>% filter(NBlocks == 1) %>% select(-NBlocks), RLS) # 2,241
nrow(RLS2) - nrow(RLS) # -2,241

# Remove species with .sp and .spp
RLS2 <- RLS2[which(grepl("sp.", as.character(RLS2$SPECIES_NAME), fixed = T) == F), ]
RLS2 <- RLS2[which(grepl("spp.", as.character(RLS2$SPECIES_NAME), fixed = T) == F), ]
nrow(RLS2) # 425,164
430191 - 425164 # 5027 records which don't have species fully resolved. 

# Rename species name for convienience
RLS2 <- RLS2 %>% dplyr::rename(SpeciesName = SPECIES_NAME)

# Total number of SurveyIDs 
length(unique(RLS2$SurveyID))

# Total number of Sites
length(unique(RLS2$SiteCode))




### Remove recruits and work with adult abundances only.  ----   
# Read in species traits to obtain maximum lengths to define adults
SpeciesTraits <- read.csv('data_raw/species_traits.csv')

SpeciesTraits %>% filter(!.$MaxLength == '') %>% nrow(.) - nrow(SpeciesTraits)

# Remove species without maximum length data
SpeciesTraits <- SpeciesTraits %>% filter(!.$MaxLength == '') %>% dplyr::rename(CurrentSpeciesName = SpeciesName)
SpeciesTraits$SpeciesName <- SpeciesTraits$CurrentSpeciesName; SpeciesTraits$CurrentSpeciesName <- NULL

# Extract only columns of interest
MaximumLengths <- SpeciesTraits %>% select(SpeciesName, MaxLength)
MaximumLengths$SpeciesName <- as.character(MaximumLengths$SpeciesName)

# Merge in species maximum size. 
RLSDT <- data.table(RLS2)

# Remove duplicates from maximum length data. 
Duplicates <- data.table(MaximumLengths)[,count := .N,key = SpeciesName][count>1, which=F] %>% unique()
NonDuplicates <- data.table(MaximumLengths)[,count := .N,key = SpeciesName][count==1, which=F] %>% unique()
AllMaximumLengths <- rbind(Duplicates, NonDuplicates)
AllMaximumLengths <- AllMaximumLengths[which(AllMaximumLengths$SpeciesName %in% RLS2$SpeciesName),] 
AllMaximumLengths <- AllMaximumLengths %>% 
  group_by(SpeciesName) %>% 
  do(MaxLength = mean(.$MaxLength, na.rm = T)) %>% 
  unnest(MaxLength)

# Join RLS to maximum length data. 
RLS3 <- left_join(RLS2,  AllMaximumLengths) # 425,164
identical(nrow(RLS3), nrow(RLS2))           # Check merge.

# Code to identify abundances based on 25% threshold 
RLS3$AdultThreshold25 <- RLS3$MaxLength*0.25 # Using 0.25. 
RLS3$AdultThreshold40 <- RLS3$MaxLength*0.4  # Using 0.4 based on Rick Stuart-Smith's suggestion. 

# Converting function to vectorised form. Estimating number of individuals in difference size classes. 
AbundanceData <- RLS3 %>% select(SpeciesName, SiteCode, Site.name, SiteLat, SiteLong, SurveyDate, Depth, SurveyID, Method, Block, 
                                 AdultThreshold25, AdultThreshold40, TotalAbundance, Unsized,
                                 X2.5, X5.0, X7.5, X10.0, X12.5, X15.0, 
                                 X20.0, X25.0, X30.0, X35.0, X40.0, X50.0, X62.5, 
                                 X75.0, X87.5, X100.0, X112.5, X125.0, X137.5,
                                 X150.0, X162.5, X175.0, X187.5, X200.0, X250.0, 
                                 X300.0, X350.0, X400.0)

# Create logical classes as columns identifying whether adult or recruit based on 25% and 40% of maximum length thresholds. 
SizeData <- gather(AbundanceData %>% select(-Unsized), SizeClass, SizeClassAbundance, X2.5:X400.0)
SizeData$SizeClass <- as.numeric(gsub('X', '', SizeData$SizeClass))
SizeData$AdultClass25 <- SizeData$SizeClass > SizeData$AdultThreshold25
SizeData$AdultClass40 <- SizeData$SizeClass > SizeData$AdultThreshold40
SizeData$RecruitClass25 <- SizeData$SizeClass < SizeData$AdultThreshold25
SizeData$RecruitClass40 <- SizeData$SizeClass < SizeData$AdultThreshold40

# Fill columns with abundances
SizeData$AdultClass25 <- ifelse(SizeData$AdultClass25 == T, SizeData$SizeClassAbundance, NA)
SizeData$AdultClass40 <- ifelse(SizeData$AdultClass40 == T, SizeData$SizeClassAbundance, NA)
SizeData$RecruitClass25 <- ifelse(SizeData$RecruitClass25 == T, SizeData$SizeClassAbundance, NA)
SizeData$RecruitClass40 <- ifelse(SizeData$RecruitClass40 == T, SizeData$SizeClassAbundance, NA)

# Sum across species to get adult and recruit abundances
SizeDataDT <- data.table(SizeData)
setkey(SizeDataDT, SpeciesName, SiteCode, Site.name, SiteLat, SiteLong, SurveyDate, Depth, SurveyID, Method, Block)
SizeAbundanceData <- SizeDataDT[, list(AbundanceAdult25 = sum(AdultClass25, na.rm = T),
                                       AbundanceAdult40 = sum(AdultClass40, na.rm = T),
                                       AbundanceRecruit25 = sum(RecruitClass25, na.rm = T),
                                       AbundanceRecruit40 = sum(RecruitClass40, na.rm = T)), by = key(SizeDataDT)]

# Join back in to dataframe to work with new abundances. 
RLS4 <- as_data_frame(left_join(SizeAbundanceData, RLS3 %>% select(SpeciesName, SiteCode, Site.name, SiteLat, SiteLong, SurveyDate, Depth, SurveyID, Method, Block, 
                                                                   Country, State, Location, ECOregion, province, realm , Visibility, Rugosity..max., Rugosity..mean., 
                                                                   Diver, RecordedSpecies, SpeciesID, PHYLUM, CLASS, ORDER, FAMILY, GENUS, MaxLength, AdultThreshold25, AdultThreshold40)))
nrow(RLS4) - nrow(RLS3)   # 0 different in length
# 425,164 matches the row number expected. 

# Removes all 0s that are from the size class aggregation. 
RLS4 <- RLS4 %>% filter(AbundanceAdult40 > 0) %>% select(-AbundanceAdult25, -AbundanceRecruit25, -AbundanceRecruit40)





### Create geobuffer for absences in occupancy modelling. ----

# Remove hanging objects from workspace 
rm(RLS, RLS2, RLS3, RLSDT, AbundanceData, AllMaximumLengths, Duplicates, NonDuplicates, SizeAbundanceData, SizeData, SizeDataDT, SpeciesTraits, TestNBlocks, RemovedData, MaximumLengths)

# INITIAL TEST OF METHOD 
#TestSpecies <- RLS4 %>% filter(SpeciesName == 'Labroides dimidiatus')
#CoordTest <- coordinates(cbind(TestSpecies$SiteLong, TestSpecies$SiteLat))

# Compare buffer units. Difference is marginal.
#SpeciesPointData <- SpatialPoints(CoordTest, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84  +no_defs +towgs84=0,0,0"))
#SpeciesPointData_meters <- spTransform(SpeciesPointData, CRS("+proj=cea +units=km"))

#bufferedPoints <- gBuffer(SpeciesPointData, width=10, byid=T)
#bufferedPoints_meters <- gBuffer(SpeciesPointData_meters, width=1000, byid=T)

#bufferedPoints_SPDF <- SpatialPolygonsDataFrame(bufferedPoints, data = data.frame(TestSpecies))
#bufferedPoints_meters_SPDF <- SpatialPolygonsDataFrame(bufferedPoints_meters, data = data.frame(TestSpecies))

#bufferedPoints_SPDF@data$id <- rownames(bufferedPoints_SPDF@data)
#bufferedPoints_meters_SPDF@data$id <- rownames(bufferedPoints_meters_SPDF@data)
#bufferedPoints_meters_SPDF <- spTransform(bufferedPoints_meters_SPDF, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84  +no_defs +towgs84=0,0,0"))

# create a data.frame from our spatial object
#bufferedPoints_DF <- fortify(bufferedPoints_SPDF, region = "id")
#bufferedPoints_meters_DF <- fortify(bufferedPoints_meters_SPDF, region = "id")

# Test plots on map to ensure scale expected
#world <- map_data(map('world', fill = T, bg = 'red', colour = 'blue', plot = F))
#ggplot() + 
#  geom_map(data=world, map=world,
#           aes(x=long, y=lat, map_id=region), col = 'gray50', fill="gray50", size = 0.1) + 
#  theme_minimal() +
#  theme(axis.title = element_blank(),
#        axis.text = element_blank())+
#  geom_polygon(data = bufferedPoints_DF, aes(x = long, y = lat, group = id), col = 'red', fill = NA) +
#  geom_polygon(data = bufferedPoints_meters_DF, aes(x = long, y = lat, group = id), col = 'blue', fill = NA) +
#  
#  geom_point(data = TestSpecies %>% .[which(.$AbundanceAdult40 != 0),], aes(y = SiteLat, x = SiteLong), pch = 21) +
#  coord_fixed() + 
#  xlim(c(min(TestSpecies %>% .[which(.$AbundanceAdult40 != 0),] %>% .$SiteLong), 
#         max(TestSpecies %>% .[which(.$AbundanceAdult40 != 0),] %>% .$SiteLong))) + 
#  ylim(c(min(TestSpecies %>% .[which(.$AbundanceAdult40 != 0),] %>% .$SiteLat), 
#         max(TestSpecies%>% .[which(.$AbundanceAdult40 != 0),] %>% .$SiteLat))) + 
#  ggtitle(unique(TestSpecies$SpeciesName))


# Convert all RLS data to spatial points object
AllDataPoint <- SpatialPoints(coordinates(unique(cbind(RLS4$SiteLong, RLS4$SiteLat))), 
                              proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +units=m'))

# Extract all sites
All_sites <- RLS4 %>% dplyr::select(realm, ECOregion, Block, Method, SurveyID, SurveyDate, SiteCode, Site.name, SiteLat, SiteLong) %>% unique()

### Write geospatial 10 degree buffer into for loop for all species -> then turn all to 0, bind with abundance data and sum. 
SpeciesBuffered0s <- list()
for(i in 1:length(unique(RLS4$SpeciesName))){
  
  # Keep count
  if(i %% 10==0) {
    # Print on the screen some message
    cat(paste0("iteration: ", i, "\n"))
  }
  Sys.sleep(0.1) # Just for waiting a bit in this example
  
  # Obtain species presences and turn into spatial data. 
  SpeciesPresences <- RLS4 %>% filter(SpeciesName == unique(RLS4$SpeciesName)[i])
  SpeciesPresences_coord <- coordinates(cbind(SpeciesPresences$SiteLong, SpeciesPresences$SiteLat))
  SpeciesPresences_sp <- SpatialPoints(SpeciesPresences_coord, proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +units=m'))
  
  # Convert whole database into spatial data to subset absences from
  #  AllDataPoint <- SpatialPoints(coordinates(unique(cbind(RLS4$SiteLong, RLS4$SiteLat))), 
  #                                proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'))
  
  # Create polygon buffer around all species' presences
  bufferedPoints <- gBuffer(SpeciesPresences_sp, width=10, byid=T)
  
  # Extract all
  SpeciesBuffered0s[[i]] <- as_data_frame(left_join(data.frame(AllDataPoint[bufferedPoints]), All_sites %>% unique(), by = c('coords.x1' = 'SiteLong', 'coords.x2' = 'SiteLat'))) %>% na.omit()
  
  SpeciesBuffered0s[[i]]$AbundanceAdult40 <- 0
  
  colnames(SpeciesBuffered0s[[i]])[1:2] <- c('SiteLong', 'SiteLat')
  
  SpeciesBuffered0s[[i]]$SpeciesName <- unique(RLS4$SpeciesName)[i]
  
}

# Create dataframe from list
SpeciesBuffered0s_V2 <- as_data_frame(do.call(rbind, SpeciesBuffered0s))

# Bind in abundance data and absences data.
RLS5 <- rbind(SpeciesBuffered0s_V2, RLS4 %>% dplyr::select(colnames(SpeciesBuffered0s_V2)))

# Sum abundances to remove double 0s
RLS5_DT <- data.table(RLS5)
setkey(RLS5_DT, SpeciesName, ECOregion, realm, Block, Method, SurveyID, SurveyDate, SiteCode, Site.name, SiteLat, SiteLong)
RLS5_DT_2 <- RLS5_DT[, j = list(AbundanceAdult40 = sum(AbundanceAdult40, na.rm = T)), 
                     by = key(RLS5_DT)]

RLS6 <- as_data_frame(RLS5_DT_2)



### Estimate first set of quality controls (Max abundance at block level is > 3)  ---- 

# Remove hanging objects
rm(SpeciesBuffered0s, SpeciesBuffered0s_V2, RLS5_DT_2, RLS5_DT, RLS5, RLS4, All_sites, 
   SpeciesPresences, SpeciesPresences_coord, AllDataPoint, bufferedPoints, i, SpeciesPresences_sp)

### Maximum abundance when present (N > 3)
RLS7 <- RLS6 %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(Confidence_MaxAbun = purrr::map(data, ~max(.$AbundanceAdult40, na.rm = T))) %>% 
  unnest(Confidence_MaxAbun)

# Create confidence score. 
RLS7$Confidence_MaxAbun <- ifelse(RLS7$Confidence_MaxAbun > 3, 1, 0)
hist(RLS7$Confidence_MaxAbun)
sum(RLS7$Confidence_MaxAbun>0)  # 1644 species maximum abundances > 3
sum(!RLS7$Confidence_MaxAbun>0) # 711 species maximum abundances < 3

# Unnest data. 
RLS7 <- RLS7 %>% unnest(data)




### Match up data with covariates for modelling  ---- 

# Create columns to match with site temperature data scale
RLS7$Latitude.1 <- round(RLS7$SiteLat)
RLS7$Longitude.1 <- round(RLS7$SiteLong)

# Read in site temperature data from NOAA
SiteTemperature <- as_data_frame(read.csv(file = 'data_raw/SiteTemperatureValuesRLS_04_09_2017.csv')) %>% select(-X)
SiteTemperature$Latitude.1 <- SiteTemperature$Latitude;SiteTemperature$Latitude<-NULL
SiteTemperature$Longitude.1 <- SiteTemperature$Longitude;SiteTemperature$Longitude<-NULL

# Join in site temperature data and RLS data. 
RLS8 <- left_join(RLS7, SiteTemperature)

# Rename for mergeing
RLS8 <- dplyr::rename(RLS8, Latitude = SiteLat, Longitude = SiteLong)

# Read in covariates and warming rate data
WarmingRates <- read.csv(file = 'data_raw/RLSWarmingRate_04_09_2017.csv') %>% select(-X, -Latitude, -Longitude)
RLS_Site_Covariates <- as_data_frame(read.csv(file = 'data_raw/RLS_Site_Covariates_V2.csv')) %>% dplyr::select(-X)
RLS_Site_Covariates <- dplyr::rename(RLS_Site_Covariates, Latitude = lat, Longitude = long)

# Rescale to log 
RLS_Site_Covariates$HumPop50km2 <- log10(RLS_Site_Covariates$HumPop50km2 + 1) 
RLS_Site_Covariates$ReefAreaIn15km <- log10(RLS_Site_Covariates$ReefAreaIn15km + 1)

# Join with RLS data
RLS8 <- left_join(RLS8, RLS_Site_Covariates)
RLS8 <- left_join(RLS8, na.omit(WarmingRates))


### Create sampling intensity column ---- 

# Round temperature data and assess number of samples per °C for each species. 
RLS8$MeanSiteSST_NOAA_ROUNDED <- round(RLS8$MeanSiteSST_NOAA)

# Count number of samples per temperature band. 
SamplingIntensityTemps  <- RLS8 %>% group_by(SpeciesName) %>% nest() %>% mutate(TemperatureTable = purrr::map(data, ~data.frame(table(.$MeanSiteSST_NOAA_ROUNDED))))

# Handle sampling intensity data. 
SamplingIntensityTemps2 <- SamplingIntensityTemps %>% unnest(TemperatureTable) %>% dplyr::rename(., 
                                                                                                 MeanSiteSST_NOAA_ROUNDED = Var1, 
                                                                                                 SamplingIntensity = Freq)

# Convert to numeric
SamplingIntensityTemps2$MeanSiteSST_NOAA_ROUNDED <- SamplingIntensityTemps2$MeanSiteSST_NOAA_ROUNDED %>% as.numeric()

RLS11 <- left_join(RLS8, SamplingIntensityTemps2)

### Aggregate data from block to site and scale covariates ----

# Remove hanging objects
rm(RLS10, RLS8, RLS6, RLS7, SpeciesMeanTemps, WarmingRates, SiteTemperature, SamplingIntensityTemps, SamplingIntensityTemps2)

# Take columns of interest. Scaled columns are actually removed here so there is some redundancy in the above code. 
RLS12 <- RLS11 %>% dplyr::select(SpeciesName, Block, SurveyID, SiteCode, ECOregion, AbundanceAdult40, SamplingIntensity,
                                 ReefAreaIn15km, HumPop50km2, npp_mean, MeanSiteSST_NOAA)

# Sum abundance data across sampling blocks. 
RLS12_DT <- data.table(RLS12)
setkey(RLS12_DT, SpeciesName, SurveyID, SiteCode)
RLS13_DT <- RLS12_DT[ ,.(AbundanceAdult40 = sum(AbundanceAdult40),
                         SamplingIntensity = mean(SamplingIntensity, na.rm = T),
                         ReefAreaIn15km = mean(ReefAreaIn15km, na.rm = T),
                         HumPop50km2 = mean(HumPop50km2, na.rm = T),
                         npp_mean = mean(npp_mean, na.rm = T),
                         MeanSiteSST_NOAA = mean(MeanSiteSST_NOAA, na.rm = T)), key = key(RLS12_DT)]

# Average abundance data across sampling surveyIDs
setkey(RLS13_DT, SpeciesName, SiteCode)
RLS14_DT <- RLS13_DT[ ,.(AbundanceAdult40 = ceiling(mean(AbundanceAdult40, na.rm = T)),
                         SamplingIntensity = mean(SamplingIntensity, na.rm = T),
                         ReefAreaIn15km = mean(ReefAreaIn15km, na.rm = T),
                         HumPop50km2 = mean(HumPop50km2, na.rm = T),
                         npp_mean = mean(npp_mean, na.rm = T),
                         MeanSiteSST_NOAA = mean(MeanSiteSST_NOAA, na.rm = T)), key = key(RLS13_DT)]

RLS14 <- as_data_frame(RLS14_DT)

# Data is now SiteCode level. 
RLS15 <- left_join(RLS14, RLS12 %>% dplyr::select(SpeciesName, SiteCode, ECOregion) %>% unique()) # 1,409,932 rows. 

# Table number of observations per random effect
#hist(table(RLS15$SpeciesName))
#hist(table(RLS15$SiteCode))

# Check and remove NAs. 
colSums(is.na(RLS15))
RLS15 <- RLS15[!is.na(RLS15$ReefAreaIn15km), ]
RLS15 <- RLS15[!is.na(RLS15$HumPop50km2), ]
RLS15 <- RLS15[!is.na(RLS15$npp_mean), ]
colSums(is.na(RLS15))

# Explore site level covariate correlations 
#MyVar <- c('ReefAreaIn15km', 'HumPop50km2', 'npp_mean', 'MeanSiteSST_NOAA')
#Covs <- unique(RLS15[,MyVar])
#Mydotplot(log(Covs)) # Function sourced from Zurr code from course (load entire script)
#corvif(Covs)  # Not much problem here. 
#Mypairs(Covs) # Reef and temperature are problematic. Model temperate and tropical species seperately. 

# Create presence column for modelling with binomial glmm.  
RLS15$Presence <- ifelse(RLS15$AbundanceAdult40 > 0, 1, 0)
RLS15$OLRE <- as.factor(1:nrow(RLS15))

# Define maximum abundance here. This will be used as a scaling factor in physiological models later
MaximumAbundance <- RLS15 %>% filter(AbundanceAdult40 > 0) %>%   
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(MaxAbundance = purrr::map(data, ~round(quantile(.$AbundanceAdult40, 0.95)))) %>% ungroup() %>% 
  unnest(MaxAbundance) %>% select(-data)

# Join in maximum abundance with RLSData. 
RLS16 <- left_join(RLS15, MaximumAbundance)




### Apply stage 1 and 2 of confidence scoring pre-modelling phase ----

# Remove hanging objects
rm(RLS15, RLS14, RLS14_DT, RLS13_DT, RLS11, RLS12, RLS12_DT, RLS13_DT)

# Check memory usage 
print(object.size(x=lapply(ls(), get)), units="Mb") # 189.5Mb

# 1. Confidence score based on number of observations. 

# >= 30 = 3
#  < 30 = 0

# 1. Confidence score based on number of observations. 
RLS17 <- RLS16 %>% 
  #filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(N_Obs = purrr::map(data, ~sum(.$AbundanceAdult40 > 0))) %>% 
  unnest(N_Obs) %>% select(-data) %>% 
  left_join(RLS16, .)

RLS17$Confidence_NObs <- ifelse(RLS17$N_Obs >= 30, 3, 0)

table(RLS17$Confidence_NObs)/nrow(RLS17)
# 0         3 
# 0.40      0.59 


# 2. Confidence score modified based on thermal range of data. 

# > 3°C range in temperature = keep score from 1. 
# < 3°C range in temperature and score = 3, demote to 2. 
# < 3°C range in temperature and score = 2, demote to 1. 

# Estimate known range in temperatures for each species. 
RLS18 <- RLS17 %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs = purrr::map(data, ~quantile(.$MeanSiteSST_NOAA, 0.99)),
         T_Lower_Obs = purrr::map(data, ~quantile(.$MeanSiteSST_NOAA, 0.01))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Range_Obs = purrr::map(data, ~.$T_Upper_Obs - .$T_Lower_Obs)) %>% 
  unnest(T_Range_Obs) %>% select(-data) %>% 
  left_join(RLS17, .)

# Assign thermal randge confidence scores
RLS18$Confidence_TRange_Obs <- NA
RLS18$Confidence_TRange_Obs[which(RLS18$Confidence_NObs == 3)] <- ifelse(RLS18$T_Range_Obs[which(RLS18$Confidence_NObs == 3)] >= 3, 3, 0)
RLS18$Confidence_TRange_Obs[which(RLS18$Confidence_NObs < 3)]  <- ifelse(RLS18$T_Range_Obs[which(RLS18$Confidence_NObs < 3)] >= 3, 0, 0)
table(RLS18$Confidence_TRange_Obs)/nrow(RLS18)
# 0        3 
# 0.43     0.56 


# Confidence score based on number of absences above and below observed maximum and minimum sampling temperatures. 
RLS18_V2 <- RLS18 %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  
  # Estimate upper and lower limits
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanSiteSST_NOAA)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanSiteSST_NOAA)), 
         T_Mean_Obs  = purrr::map(data, ~ mean(.$MeanSiteSST_NOAA, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Mean_Obs) %>% ungroup() %>% select(-data) %>% left_join(RLS18, ., by = 'SpeciesName') %>% 
  
  # Count number of rows > 0 above and below these limits. 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(N_Absences_T_Upper = purrr::map(data, ~nrow(.[which(.$MeanSiteSST_NOAA >= .$T_Upper_Obs), ] %>% filter(AbundanceAdult40 == 0))),
         N_Absences_T_Lower = purrr::map(data, ~nrow(.[which(.$MeanSiteSST_NOAA <= .$T_Lower_Obs), ] %>% filter(AbundanceAdult40 == 0)))) %>% 
  unnest(N_Absences_T_Upper, N_Absences_T_Lower) %>% ungroup() %>% select(-data) %>% left_join(RLS18, ., by = 'SpeciesName')

hist(log(RLS18_V2$N_Absences_T_Upper + 1))
hist(log(RLS18_V2$N_Absences_T_Lower + 1)) # On average there are far fewer species with poorly estimated cold thermal limits than upper limits. 

# Mean temperature of species with upper limits that cannot be estimated. 
RLS18_V2 %>% filter(N_Absences_T_Upper == 0) %>% .$MeanSiteSST_NOAA %>% mean

# Create confidence scores for these parameters. 
RLS18_V2$Confidence_Occ_Tupper  <-NA
RLS18_V2$Confidence_Occ_Tupper  <- ifelse(RLS18_V2$N_Absences_T_Upper > 10, 1, 0)
RLS18_V2$Confidence_Occ_Tlower  <- NA
RLS18_V2$Confidence_Occ_Tlower  <- ifelse(RLS18_V2$N_Absences_T_Lower > 10, 1, 0)


# Estimate maximum temperature in sampling range INCLUDING PRESENCES
RLS18_V2 <- RLS18_V2 %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  # Estimate upper and lower limits
  mutate(T_Upper_Absences = purrr::map(data, ~max(.$MeanSiteSST_NOAA)),
         T_Lower_Absences = purrr::map(data, ~min(.$MeanSiteSST_NOAA)), 
         T_Mean_Absences  = purrr::map(data, ~ mean(.$MeanSiteSST_NOAA))) %>% 
  unnest(T_Upper_Absences, T_Lower_Absences, T_Mean_Absences) %>% ungroup() %>% select(-data) %>% left_join(RLS18_V2, ., by = 'SpeciesName')

# 1650 species with low confidence scores. 
RLS_LowConfidence_Species  <- unique(RLS18_V2[which(RLS18_V2$Confidence_NObs  != 3 | RLS18_V2$Confidence_TRange_Obs != 3), 'SpeciesName'])
# 705 species with high confidence scores. 
RLS_HighConfidence_Species <- unique(RLS18_V2[-which(RLS18_V2$Confidence_NObs != 3 | RLS18_V2$Confidence_TRange_Obs != 3), 'SpeciesName'])

# Filter RLS data into two streams. Low confidence and high confidence. 
RLS_LowConfidence <- RLS18 %>% filter(SpeciesName %in% RLS_LowConfidence_Species$SpeciesName) # 616,949 observations
saveRDS(RLS_LowConfidence, file = 'data_derived/RLS_LowConfidence_Species_2017-03-27.rds')

RLS_19 <- RLS18_V2 %>% filter(SpeciesName %in% RLS_HighConfidence_Species$SpeciesName)         # 787,893 observations


nrow(RLS_19)
RLS_19 %>% filter(AbundanceAdult40 != 0) %>% nrow
length(unique(as.character(RLS_19$SiteCode)))
length(unique(as.character(RLS_19$SpeciesName)))
RLS_19 %>% do()

Survey_no <- RLS %>% select(SiteCode, SurveyID) %>% unique() %>% left_join(RLS_19 %>% select(SiteCode) %>% unique(), .)
Survey_no %>% group_by(SiteCode) %>% do(No_unique = length(unique(.$SurveyID))) %>% unnest(No_unique) %>% .$No_unique %>% mean()


### SAVE PROGRESS ----

# Remove hanging objects
rm(RLS16, RLS17, RLS18, RLS18_V2, RLS_LowConfidence)

# Save the whole workspace 
print(object.size(x=lapply(ls(), get)), units="Mb")
save.image(file = 'data_derived/01_organising-data_SAVE-IMAGE.RData')
#load(      file = 'data_derived/FinalTPC_Workflow_12_01_2018.RData')

### ----







