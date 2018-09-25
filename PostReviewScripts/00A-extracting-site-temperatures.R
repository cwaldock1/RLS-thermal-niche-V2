##############################################################
###
### Script matching species observations to temperature data
### to estimate niche upper and lower realized limits. 
###
### Extract covariates for all sites too based on data in Yeager, L.A., Marchand, P., Gill, D.A., Baum, J.K., and McPherson, J.M. In review. Queryable global layers of environmental and anthropogenic variables for marine ecosystem studies.
###
### Author: Conor A Waldock.
###
### Date of script initiation: 04/09/2017
### 
##############################################################

# --------------------------
# Load libraries and packages.
library(plyr)
library(dplyr)
library(pbapply)
library(data.table)
library(raster)
library(parallel)
library(mgcv)
library(MuMIn)
library(lubridate)
library(rgdal)


# Handling temperature data  ----------------------------
# Read in temperature datasets

# Hadley centre is at 1°x 1° grids
#had <-"/Users/caw1g15/Dropbox/PhD Work/Temperature driven community change (s-Change)/driver_data/HADISST/HadISST_sst-2.nc"
#had <- brick(had, varname="sst")
#crs(had)

# NOAA data https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html
#NOAA <- "/Users/caw1g15/Dropbox/PhD Work/Temperature driven community change (s-Change)/driver_data/NOAA_OI_SST_V2/Weekly/sst.wkmean.1990-present.nc"
#NOAA <- brick(NOAA)
#NOAA <- rotate(NOAA)
#crs(NOAA)

# READ IN FINE-SCALE TEMPERATURE DATA AND AGGREGATE TO YEAR AVERAGE AND MEAN ----

Read1 <- list.files('/Volumes/Untitled/Raster_files/CoralWatch/2017')[100]
TestRaster <- raster(paste('/Volumes/Untitled/Raster_files/CoralWatch/2017/', Read1, sep = ''))
plot(TestRaster)
crs(TestRaster) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

RLS_Sites <- SpatialPointsDataFrame(RLSTemp2[,c('Longitude', 'Latitude')], data = RLSTemp2, proj4string = crs(TestRaster))

# Check match up 
plot(TestRaster)
plot(RLS_Sites, add = T)


# Extract from 1996 (2006-10 years) to 2017 to get temperature 10 years prior to each survey
lubridate::year(dmy_hm(RLSTemp$SurveyDate)) %>% min
lubridate::year(dmy_hm(RLSTemp$SurveyDate)) %>% max
10284 * 11 * 365 # would create 41,290,260 temperature values.  

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

# START OF LOOP FOR DAILY DATA EXTRACTION FROM RLS SITES
# Aggregate to weekly temperature for each RLS site. 
Years <- 1993:2018#:2018
# Loop through years to extract from each yearly folder
#foreach(year=1:length(Years), .packages=c('raster', 'data.table', 'dplyr', 'rgdal', 'lubridate')) %dopar% { # (year in 1:length(Years)){
TimeStart <- Sys.time()
for(year in 1:length(Years)){
YearFolder <- paste('/Volumes/Untitled/Raster_files/CoralWatch/', Years[year], sep = '')
YearFiles <- list.files(YearFolder)

RLS_Sites_List <- list()

for(i in 1:length(YearFiles)){
  
  print(paste(Years[year], 'day', i, sep = '-'))
  print(Sys.time() - TimeStart)
  # Pull data 
  RasterLocation <- paste(YearFolder, YearFiles[i], sep = '/')
  TempRaster <- raster(RasterLocation)
  crs(TempRaster) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

  # Find sites that need buffering 
TestExtraction <- extract(TempRaster, RLS_Sites)
RLS_Sites_NAs <- RLS_Sites[is.na(TestExtraction),]
TestExtraction2 <- extract(TempRaster, RLS_Sites_NAs, buffer = 7500, fun = mean)
TestExtraction[is.na(TestExtraction)] <- TestExtraction2
RLS_Sites$Temperature <- TestExtraction
RLS_Sites$Date <- as.Date(sub('\\..*', '', strsplit(YearFiles[i], '_')[[1]][3]),format = "%Y%m%d")
RLS_Sites_List[[i]] <- as_data_frame(RLS_Sites)
RLS_Sites_List[[i]]$Week <- lubridate::week(RLS_Sites_List[[i]]$Date)
}
RLS_Sites_Rbind <- data.table(do.call(rbind, RLS_Sites_List))
saveRDS(RLS_Sites_Rbind, file = paste('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/', paste('RLS_Daily_Temp', as.character(Years[year]), '.rds', sep = ''), sep = ''))
rm(RLS_Sites_Rbind, RLS_Sites_List)
}

#stopCluster(cl)

# End of all yearly loops
# Now aggregate by week

# ---- 



# Correct dates from previous loop error (now fixed) ----
list.files('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction')

# Assign corrected weeks and dates
Years <- 1985:2018
for(i in 1:length(Years)){
  print(Years[i])
DailyTemps <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/RLS_Daily_Temp',Years[i],'.rds'))
Read1 <- list.files(paste0('/Volumes/Untitled/Raster_files/CoralWatch/',Years[i]))
Dates <- melt(lapply(lapply(strsplit(Read1, '_'), function(x) sub('\\..*', '', x[3])), function(x) as.Date(x, format = "%Y%m%d")), value = 'Date')
DailyTemps <- DailyTemps[order(DailyTemps$RowID),]
DailyTemps$Date <- Dates$value
DailyTemps$Week <- lubridate::week(DailyTemps$Date)
saveRDS(DailyTemps, file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp',Years[i],'.rds'))
}

DailyTemps_2001 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2001','.rds'))
DailyTemps_2002 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2002','.rds'))
DailyTemps_2003 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2003','.rds'))
DailyTemps_2004 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2004','.rds'))
DailyTemps_2005 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2005','.rds'))
DailyTemps_2006 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2006','.rds'))
DailyTemps_2007 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2007','.rds'))
DailyTemps_2008 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2008','.rds'))
DailyTemps_2009 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2009','.rds'))
DailyTemps_2010 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2010','.rds'))
DailyTemps_2011 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2011','.rds'))
DailyTemps_2012 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2012','.rds'))
DailyTemps_2013 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2013','.rds'))
DailyTemps_2014 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2014','.rds'))
DailyTemps_2015 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2015','.rds'))
DailyTemps_2016 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2016','.rds'))
DailyTemps_2017 <- readRDS(file = paste0('/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/RLS_Daily_Temp','2017','.rds'))


Daily_TempsAll <- data.table(rbind.fill(DailyTemps_2003, DailyTemps_2004, DailyTemps_2005, 
                                        DailyTemps_2006, DailyTemps_2007, DailyTemps_2008, DailyTemps_2009, DailyTemps_2010,
                                        DailyTemps_2011, DailyTemps_2012, DailyTemps_2013, DailyTemps_2014, DailyTemps_2015,
                                        DailyTemps_2016, DailyTemps_2017))
rm(DailyTemps_2001, DailyTemps_2002, DailyTemps_2003, DailyTemps_2004, DailyTemps_2005, 
   DailyTemps_2006, DailyTemps_2007, DailyTemps_2008, DailyTemps_2009, DailyTemps_2010,
   DailyTemps_2011, DailyTemps_2012, DailyTemps_2013, DailyTemps_2014, DailyTemps_2015,
   DailyTemps_2016, DailyTemps_2017)

# Read in RLS data ----
RLSTemp <- as_data_frame(read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv'))[,1:33]
RLSTemp2 <- RLSTemp %>% ungroup() %>%  dplyr::select(SiteLong, SiteLat, SurveyID, SurveyDate) %>% unique()
RLSTemp2$SurveyDate <- lubridate::dmy_hm(RLSTemp2$SurveyDate)
RLSTemp2$RowID <- as.character(1:nrow(RLSTemp2))
names(RLSTemp2) <- c('Longitude', 'Latitude', 'SurveyID', 'SurveyDate','RowID')
min(RLSTemp2$SurveyDate); max(RLSTemp2$SurveyDate)

# Create 10 year, 5 year 1 year and 0.5 years temperature. 
# Create single example to loop through 
YearTime  <- as.period(2, unit = 'year')

# Create a list of extraction dates for each survey ID 
ExtractionDate_List <- list()
SurveyIDs <- unique(RLSTemp2$SurveyID)
for(i in 1:length(unique(RLSTemp2$SurveyID))){
print(i)
  if(as_date(RLSTemp2$SurveyDate[i]) == as_date('2012-02-29 UTC')){ RLSTemp2$SurveyDate[i] <- RLSTemp2$SurveyDate[i] - as.period(1, unit = 'day')}else{
  ExtractionDate_List[[i]] <- data.frame(Date = seq((RLSTemp2$SurveyDate[i] - YearTime), RLSTemp2$SurveyDate[i] , by = 'days'), SurveyID = SurveyIDs[i])
  }
}

# Turn the extraction into a dataframe 
ExtractionDate_DT <- data.table(do.call(rbind.fill, ExtractionDate_List))
ExtractionDate_DT[, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run on the bigger object below so run sequentially. 

# Create similar date column in temperature data
Daily_TempsAll[1:10000000, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 
Daily_TempsAll[1000001:20000000, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 
Daily_TempsAll[2000001:30000000, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 
Daily_TempsAll[3000001:40000000, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 
Daily_TempsAll[4000001:50000000, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 
Daily_TempsAll[5000001:56346036, ExtractionCode:=paste(Date, SurveyID)] # This takes a while to run. 

# Save the above object 
#saveRDS(Daily_TempsAll, file = '/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/Daily_TempsAll.rds')
Daily_TempsAll <- readRDS(file = '/Volumes/Untitled/Raster_files/CoralWatch/RLS_Daily_Temp_Extraction/EDIT/Daily_TempsAll.rds')

# Perform extraction of only survey ID and date ranges specified 
Daily_TempsAll_V2 <- Daily_TempsAll[ExtractionCode %in% ExtractionDate_DT$ExtractionCode]

# Aggregate temperatures for each survey ID 
setkey(Daily_TempsAll_V2, SurveyID)
Daily_TempsAll_V2[, `:=`(MeanTemp=mean(Temperature), MinTemp=quantile(Temperature,0.05), MaxTemp=quantile(Temperature, 0.95), SDTemp=sd(Temperature)), by = SurveyID]

# Shorten object. 
Daily_TempsAll_V3 <- unique(Daily_TempsAll_V2[,c('SurveyID', 'Longitude', 'Latitude', 'MeanTemp', 'MinTemp', 'MaxTemp', 'SDTemp')])
write.csv(Daily_TempsAll_V3, '/Users/caw1g15/Documents/GitHub/RLS-thermal-niche-V2/data_derived/RLS_Temperature_CoralWatch0.05.csv')

ggplot(Daily_TempsAll_V3) + 
  geom_point(aes(x = MeanTemp, y = MinTemp))

ggplot(Daily_TempsAll_V3) + 
  geom_point(aes(x = MeanTemp, y = MaxTemp))

ggplot(Daily_TempsAll_V3) + 
  geom_point(aes(x = MeanTemp, y = SDTemp))


# ----

# --- 
# Function to extract 10 years data across each site. 
Temp_Trend <- function(env_dat = env_dat, occ_dat,  origin = as.Date("1800-1-1")){
  
  # calculate starting julian day for each month in env_dat
  month_intervals <-  as.numeric(env_dat@z[["Date"]] - origin)
  
  # extract environmental variable (SST here) for this point
  env_val <- raster::extract(
    env_dat,
    coordinates(occ_dat),
    layer = findInterval(month_intervals[length(month_intervals)-520], month_intervals),
    nl = 521
    #df = F
  )
  
  Date <- as.Date(gsub('X','',dimnames(env_val)[[2]]), format = '%Y.%m.%d')
  
  RowID <- occ_dat@data$RowID
  
  dimnames(env_val) <- NULL
  
  env_val <- data.frame(Temp = c(env_val), Date = Date, RowID = RowID)
  
  # Return the value
  return(env_val)
  
}

# ---
# Extract site level temperature data for past 10 years
# Use summer and winter temperatures to estimate species minimum and maximum temperature limits. 

# Read in RLS data and convert to site, lat and long codes. 
RLSTemp <- as_data_frame(read.csv('data_raw/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv'))[,1:33]
RLSTemp$SiteLong.1 <- round(RLSTemp$SiteLong) 
RLSTemp$SiteLat.1 <- round(RLSTemp$SiteLat)
RLSTemp2 <- RLSTemp %>% ungroup() %>%  dplyr::select(SiteLong, SiteLat, SurveyID) %>% unique()
RLSTemp2$RowID <- as.character(1:nrow(RLSTemp2))
names(RLSTemp2) <- c('Longitude', 'Latitude', 'SurveyID', 'RowID')

# Writing csv of only latitudes and longitudes for extraction of NPP and market access from https://shiny.sesync.org/apps/msec/
#RLSTemp2 <- RLSTemp %>% ungroup() %>%  dplyr::select(SiteLong, SiteLat) %>% unique()
#names(RLSTemp2) <- c('long', 'lat')
#write.csv(RLSTemp2, file = 'data/SiteLongitudeLatitudes_08_07_2015.csv')

# Convert each site to a unique spatial points dataframe.
# Create month extraction dataframes
RLSTemp2List <- split(RLSTemp2, RLSTemp2$RowID) # Split into a list
RLSTemp2SPDF <- mclapply(RLSTemp2List, function(x) SpatialPointsDataFrame(x[,c('Longitude','Latitude')], x, proj4string = crs(NOAA)), mc.cores = getOption('mc.cores', 4L)) # Convert all to spatial data frames

# Extract temperature data for past 10 years from each site. 
RLSTemp2Extraction <- pblapply(RLSTemp2SPDF, FUN = function(x) Temp_Trend(occ_dat = x, env_dat = NOAA)) 

# Convert extractions into dataframe
RLSTemp2ExtractionsDF <- do.call('rbind.fill', RLSTemp2Extraction)

# Merge back in with latitudes and longitudes. 
RLSTempExtract <- left_join(RLSTemp2, RLSTemp2ExtractionsDF)


# ---
# Mean site temperatures 
MeanSiteTemps <- RLSTempExtract %>% 
  group_by(RowID, Latitude, Longitude) %>% 
  do(MeanSiteSST_NOAA = mean(.$Temp, na.rm = T)) %>% 
  unnest(MeanSiteSST_NOAA)

# ---
# Determine north or southern for summer and winter sst estimates
RLSTempExtract$Month <- lubridate::month(RLSTempExtract$Date)

RLSTempExtract$Hemisphere <- ifelse(RLSTempExtract$Latitude >= 0, 
                                    'Northern', 
                                    'Southern')

# Mean summer temperature
RLSTempExtractSummer <- rbind(RLSTempExtract %>% filter(Hemisphere == 'Northern') %>% filter(Month %in% c(5,6,7,8)),
                              RLSTempExtract %>% filter(Hemisphere == 'Southern') %>% filter(Month %in% c(11,12,1,2)))
                              
MeanSiteSummerSST <- RLSTempExtractSummer %>% 
  group_by(RowID, Latitude, Longitude) %>% 
  do(MeanSiteSummerSST_NOAA = mean(.$Temp, na.rm = T)) %>% 
  unnest(MeanSiteSummerSST_NOAA)


# Mean winter temperature
RLSTempExtractWinter <- rbind(RLSTempExtract %>% filter(Hemisphere == 'Southern') %>% filter(Month %in% c(5,6,7,8)),
                              RLSTempExtract %>% filter(Hemisphere == 'Northern') %>% filter(Month %in% c(11,12,1,2)))

MeanSiteWinterSST <- RLSTempExtractWinter %>% 
  group_by(RowID, Latitude, Longitude) %>% 
  do(MeanSiteWinterSST_NOAA = mean(.$Temp, na.rm = T)) %>% 
  unnest(MeanSiteWinterSST_NOAA)

# Merge winter and summer site estimates in to data.
SiteTemperatures <- left_join(left_join(MeanSiteTemps, MeanSiteSummerSST), MeanSiteWinterSST)


# --- 
# Determine minimum and maximum temperatures per year for each site. 

# Extract years to group by
RLSTempExtract$Year <- lubridate::year(RLSTempExtract$Date)

# Extract only extreme values
SiteYearExtremes <- RLSTempExtract %>% 
  group_by(Latitude, Longitude, Year) %>% 
  do(MinimumSiteSST_NOAA = min(.$Temp, na.rm = T),
     MaximumSiteSST_NOAA = max(.$Temp, na.rm = T)) %>% 
  unnest(MinimumSiteSST_NOAA, MaximumSiteSST_NOAA)

# Match up sites
SiteYearExtremes2 <- left_join(SiteYearExtremes, 
          RLSTemp %>% 
            dplyr::select(SiteLong.1, SiteLat.1, SiteCode) %>% 
            dplyr::rename(Latitude = SiteLat.1, 
                          Longitude = SiteLong.1) %>% 
            unique())

SiteYearExtremes2$YearPlus1 <- SiteYearExtremes2$Year + 1

# Write as csv.
write.csv(SiteYearExtremes2, file = 'data/SiteExtremeSST_NOAA.csv', row.names = F)

# 



# ------
# Extract species site 12 weeks prior to survey for species x survey temperatures. 
Species_Temp_Trend <- function(env_dat = env_dat, occ_dat,  origin = as.Date("1800-1-1")){
  
  # calculate starting julian day for each month in env_dat
  month_intervals <-  as.numeric(env_dat@z[["Date"]] - origin)
  
  # calculate julian day for occurence data
  focal_date <- as.numeric(occ_dat@data$Date - origin)
  
  # extract environmental variable (SST here) for this point
  env_val <- raster::extract(
    env_dat,
    coordinates(occ_dat),
    layer = findInterval(focal_date, month_intervals)-11,
    nl = 12
    #df = F
  )
  
  # Date <- as.Date(gsub('X','', dimnames(env_val)[[2]]), format = '%Y.%m.%d')
  
  # RowID <- occ_dat@data$RowID
  
  dimnames(env_val) <- NULL
  
  env_val <- data.frame(SurveyID = occ_dat@data$SurveyID, 
                        Temp = mean(env_val, na.rm = T), 
                        SurveyDate = occ_dat@data$Date, 
                        RowID = occ_dat@data$RowID, 
                        stringsAsFactors = F)
  
  # Return the value
  return(env_val)
  
}

# Read in RLS data and convert to site, lat and long with species names preseved.
# Read in RLS data and convert to site, lat and long codes.
RLSTempSpp <- read.csv('data/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv')[,1:29]
RLSTempSpp$SiteLong.1 <- round(RLSTempSpp$SiteLong) 
RLSTempSpp$SiteLat.1 <- round(RLSTempSpp$SiteLat)
RLSTempSpp <- RLSTempSpp %>% ungroup() %>%  dplyr::select(SurveyID, SiteLong.1, SiteLat.1, SurveyDate) %>% unique()
RLSTempSpp$RowID <- as.character(1:nrow(RLSTempSpp))
names(RLSTempSpp) <- c('SurveyID', 'Longitude', 'Latitude', 'SurveyDate', 'RowID')

# Convert date to proper format
RLSTempSpp$Date <- as.Date(RLSTempSpp$SurveyDate, format = '%d/%m/%Y')

# Convert columns from characters. 
RLSTempSpp$SurveyID <- as.character(RLSTempSpp$SurveyID)
RLSTempSpp$SurveyDate <- as.character(RLSTempSpp$SurveyDate)

# Convert each site to a unique spatial points dataframe.
# Create month extraction dataframes
RLSTempSppList <- split(RLSTempSpp, RLSTempSpp$RowID) # Split into a list

# Check minimum date isn't less than the minimum date of the raster.
RLSTempSppSPDF <- pblapply(RLSTempSppList, function(x) SpatialPointsDataFrame(x[,c('Longitude','Latitude')], x, proj4string = crs(NOAA))) #, mc.cores = getOption('mc.cores', 4L)) # Convert all to spatial data frames

rbind(Species_Temp_Trend(RLSTempSppSPDF[[1]], env_dat = NOAA),
Species_Temp_Trend(RLSTempSppSPDF[[2]], env_dat = NOAA))

# Extract temperature data for past 10 years from each site. 
#RLSTempSppExtraction2 <- pblapply(RLSTempSppSPDF[1:100], FUN = function(x) Species_Temp_Trend(occ_dat = x, env_dat = NOAA)) 

RLSTempSppExtraction <- pblapply(RLSTempSppSPDF, FUN = function(x) Species_Temp_Trend(occ_dat = x, env_dat = NOAA)) #,mc.cores = getOption('mc.cores', 4L))

# Convert extractions into dataframe
RLSTempSppExtractionDF <- do.call('rbind.fill', RLSTempSppExtraction)

# Rename date for merge 
RLSTempSppExtractionDF$Date <- RLSTempSppExtractionDF$SurveyDate; RLSTempSppExtractionDF$SurveyDate <- NULL

# Merge back in with latitudes and longitudes. 
RLSTempSppExtract <- left_join(RLSTempSppExtractionDF, 
                               RLSTempSpp %>% dplyr::select(SurveyID, Latitude, Longitude, Date, RowID))

save(RLSTempSppExtract, file = 'data/RLSTempSppExtractionDF_backup3.RData')

# Change column name
RLSTempSppExtract$TimeStampSST_NOAA <- RLSTempSppExtract$Temp; RLSTempSppExtract$Temp <- NULL


### --------
### Outputs: 

# ----------------
# SiteTemperatures
SiteTemperatures # Site level mean, summer and winter temperatures for 10 years 
write.csv(SiteTemperatures, file = 'data/SiteTemperatureValuesRLS_04_09_2017.csv')

# ----------------------------------
# Species x site x date temperatures
# Gives the temperature 12 weeks prior to sampling for each species site date combination. 
RLSTempSppExtract
write.csv(RLSTempSppExtract, file = 'data/SpeciesTimeStampTemperatureValuesRLS_04_09_2017.csv')

### --------------------------------
### Extract site level warming rates

### IPCC warming rates available from https://www.esrl.noaa.gov/psd/ipcc/ocn/

IPCCWarmingRates <- "/Users/caw1g15/Dropbox/PhD Work/Tasmania RLSTemp/Data/SSTanomaly_RCP85.nc"
setwd('/Users/caw1g15/Dropbox/PhD Work/Tasmania RLSTemp/Data')
IPCCWarmingRates2 <- raster::raster('SSTanomaly_RCP85.nc', level = 1, varname = 'anomaly')
IPCCWarmingRates3 <- raster::rotate(IPCCWarmingRates2)
crs(IPCCWarmingRates3)

# Match up unique latitudes and longitudes to sites in new RLS data. 
RLSTemp <- read.csv('data/RLS_GLOBAL_ABUND_BIOMASS_20170507.csv')[,1:29]
RLSTemp$SiteLong.1 <- round(RLSTemp$SiteLong) 
RLSTemp$SiteLat.1 <- round(RLSTemp$SiteLat)
RLSTemp2 <- RLSTemp %>% ungroup() %>%  dplyr::select(SiteLong.1, SiteLat.1) %>% unique()
names(RLSTemp2) <- c('Longitude', 'Latitude')

# Conver to coordiantes 
coordinates(RLSTemp2) <- RLSTemp2
crs(RLSTemp2) <-  '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

# Extract values from IPCC raster
RLSTemp2$WarmingRates <- extract(IPCCWarmingRates3, RLSTemp2)

# Extract values using buffer for those too close to coast. 
RLSTempNAs <- RLSTemp2[is.na(RLSTemp2$WarmingRates),]
RLSTempNAs$WarmingRates <- extract(IPCCWarmingRates3, RLSTempNAs, buffer = 0.5, df = T)[,2]

# Bind together results
RLSWarmingRates <- rbind(RLSTemp2, RLSTempNAs)

# Merge warming rates in with SiteCode
RLSTempMerge <- dplyr::rename(RLSTemp, Longitude = SiteLong.1, Latitude = SiteLat.1)
RLSWarmingRateFinal <- left_join(as_data_frame(RLSTempMerge) %>% dplyr::select(Longitude, Latitude, SiteCode),
                                 as_data_frame(RLSWarmingRates) %>% dplyr::select(Longitude, Latitude, WarmingRates)) %>% unique()


RLSWarmingRateFinal$WarmingRates <- RLSWarmingRateFinal$WarmingRates / 10

# Re-set wd and save
setwd('/Users/caw1g15/Documents/GitHub/RLSThermalResilience')
write.csv(RLSWarmingRateFinal, file = 'data/RLSWarmingRate_04_09_2017.csv')


### --- Extraction of covariate information from MSEC layers.

# - Human population density
HP50km2 <- "data/MSECData_Yeager2017/msec_humanpop_50km.nc"
HP50km2 <- brick(HP50km2)
HP50km2 <- rotate(HP50km2)
crs(HP50km2)

HP50km2 <- HP50km2$X2015

RLS_SP <- SpatialPoints(RLSTemp2)
crs(RLS_SP) <- crs(HP50km2)
RLSTemp2$HumPop50km2 <- extract(x = HP50km2, y = RLS_SP, na.rm = T)


# - Reef area
ReefArea <- "data/MSECData_Yeager2017/msec_reefarea_15km.nc"
ReefArea <- raster(ReefArea)
ReefArea <- rotate(ReefArea)
crs(ReefArea)

RLSTemp2 <- RLSTemp %>% ungroup() %>%  dplyr::select(SiteLong, SiteLat) %>% unique()
names(RLSTemp2) <- c('long', 'lat')

RLS_SP <- SpatialPoints(RLSTemp2)
crs(RLS_SP) <- crs(ReefArea)
RLSTemp2$ReefAreaIn15km <- extract(x = ReefArea, y = RLS_SP, na.rm = T)


# - NPP
NPP_Mean <- "data/MSECData_Yeager2017/msec_npp_mean.nc"
NPP_Mean <- raster(NPP_Mean)
NPP_Mean <- rotate(NPP_Mean)
RLSTemp2$npp_mean <- extract(x = NPP_Mean, y = RLS_SP, na.rm = T)

NPP_Min <- "data/MSECData_Yeager2017/msec_npp_min.nc"
NPP_Min <- raster(NPP_Min)
NPP_Min <- rotate(NPP_Min)
RLSTemp2$npp_min <- extract(x = NPP_Min, y = RLS_SP, na.rm = T)

NPP_Max <- "data/MSECData_Yeager2017/msec_npp_max.nc"
NPP_Max <- raster(NPP_Max)
NPP_Max <- rotate(NPP_Max)
RLSTemp2$npp_max <- extract(x = NPP_Max, y = RLS_SP, na.rm = T)

NPP_SD <- "data/MSECData_Yeager2017/msec_npp_sd.nc"
NPP_SD <- raster(NPP_SD)
NPP_SD <- rotate(NPP_SD)
RLSTemp2$npp_sd <- extract(x = NPP_SD, y = RLS_SP, na.rm = T)

NPP_SD_inter <- "data/MSECData_Yeager2017/msec_npp_sdinter.nc"
NPP_SD_inter <- raster(NPP_SD_inter)
NPP_SD_inter <- rotate(NPP_SD_inter)
RLSTemp2$npp_interann_sd <- extract(x = NPP_SD_inter, y = RLS_SP, na.rm = T)

NPP_flag <- "data/MSECData_Yeager2017/msec_npp_flag.nc"
NPP_flag <- raster(NPP_flag)
NPP_flag <- rotate(NPP_flag)
RLSTemp2$npp_flag <- extract(x = NPP_flag, y = RLS_SP, na.rm = T)


# - Distance to market
Market <- "data/MSECData_Yeager2017/msec_distmarket.nc"
Market <- raster(Market)
Market <- rotate(Market)
crs(Market)
RLSTemp2$dist_market <- extract(x = Market, y = RLS_SP, na.rm = T)

write.csv(RLSTemp2, file =  'data/RLS_Site_Covariates_V2.csv')
