# This script was used to tune Maxent for focal species using a bunch of values of hyperparameters based on some references:

# Blair, M. E., Sterling, E. J., Dusch, M., Raxworthy, C. J., & Pearson, R. G. (2013). Ecological divergence and speciation between lemur (Eulemur) sister species in 
# Madagascar. Journal of evolutionary biology, 26(8), 1790-1801.
# Note: The approach results in regularization parameters = 1 for all species based on AUC for 25% held-out test data.

# Elith, J., Kearney, M., & Phillips, S. (2010). The art of modelling range‐shifting species. Methods in ecology and evolution, 1(4), 330-342.
# Similar source: http://josephgrinnell.weebly.com/exercise-tuning-maxent-using-beta-and-aic.html
# Note: This paper uses beta = 10 and hinge feature only. 

# https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix1_case_study.md
# Xiao tested the difference between betamultiplier = 1 and betamultiplier = 100

# Kass, J. M., Muscarella, R., Galante, P. J., Bohl, C. L., Pinilla‐Buitrago, G. E., Boria, R. A., ... & Anderson, R. P. (2021). ENMeval 2.0: Redesigned for customizable 
# and reproducible modeling of species’ niches and distributions. Methods in Ecology and Evolution, 12(9), 1602-1608.

# Cobos, M. E., Peterson, A. T., Barve, N., & Osorio-Olvera, L. (2019). kuenm: an R package for detailed development of ecological niche models using Maxent. PeerJ, 7, e6281.
# Note: Jamie's package can work with 'kuenm' package for evaluating different choices of hyperparameters and paritioning methods using AICc
# Source: https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html

# Give Java Virtual machine 45GB RAM
# Give Java Virtual machine 45GB RAM
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx45g"))

# Try 20GB RAM
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx20480m"))
gc()

# Install package
# install.packages("ENMeval")

# Load package
library(ENMeval)
library(rJava)
library(dismo)

library(sp)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(raster)

# Date: 05/14/2024
# This will need the install of Rtools version 4.4.6104
library(kuenm)

# Date: 08/12/2024
# This will need the install of Rtools version 4.4.6104
library(ecospat)


# Load environmental data first for getting crs and removing occurrences with NA in those layers
# Load data
feature_list <- list.files(path = "processed_raster_layers/", pattern = '\\.tif$', full.names = TRUE)

# Regular expression: https://stackoverflow.com/questions/4876813/using-r-to-list-all-files-with-a-specified-extension

# Combine all predictors
feature_all <- stack(feature_list)

# Make sure to declare the categorical variable as a factor
feature_all$dominant_land_cover <- raster::as.factor(feature_all$dominant_land_cover)

# Create a name list without character
spplist <- c()
for (i in 1:length(occlist)) {
  message(i)
  # Get the file name after the folder location
  spplist_i <- sub(".*Output/occurrence/", "", occlist[i])
  
  # Remove .csv after species name
  spplist[i] <- sub("*.csv$", "", spplist_i)
}

# list all the csv files with occurrences
occlist <- list.files("Output/occurrence/", 
                      pattern = ".csv", full.names = TRUE)
bg_list <- list.files("Output/bg", 
                      pattern = ".csv", full.names = TRUE)

crs_list <- readRDS("crs_list.rds")

calibration_area_list <- list.files("Output/calibration_area/", pattern = ".shp", full.names = TRUE)

# Date: 08/12/2024
# Change: Crop all features using the buffered IUCN range maps because global scale rasters eat a lot of RAM
feature_list <- c()

Sys.time()
for (i in 1:length(occlist)) {
  calibration_i <- read_sf(calibration_area_list[i])
  feature_list[[i]] <- crop(feature_all, calibration_i)
  gc()
  feature_list[[i]] <- mask(feature_list[[i]], calibration_i)
  gc()
  writeRaster(feature_list[[i]],
              # a series of names for output files
              filename = paste0("Output/calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), names(feature_list[[i]]), ".asc"), 
              format = "ascii", ## the output format
              bylayer = TRUE, ## this will save a series of layers
              overwrite = T)
  }
Sys.time()

# Date: 08/14/2024
# Change: Convert dominant land cover to factors!!!!


# Create the list to save occurrence and bg with environmental data extracted
extraction_start <- Sys.time()
p_list <- c()
bg_thin_model <- c()

# Extract environmental data for occurrences
for (i in 1:length(occlist)) {
  temp_occ <- read.csv(occlist[i])
  coordinates(temp_occ) <- ~x + y
  crs(temp_occ) <- crs_list[[i]]
  temp_occ <- extract(feature_all, temp_occ) %>% as.data.frame()
  p_list[[i]] <- temp_occ
}

# Extract environmental data for bg points
for (i in 1:length(bg_list)) {
  temp_bg_thin <- read.csv(bg_list[i])
  coordinates(temp_bg_thin) <- ~x + y
  crs(temp_bg_thin) <- crs_list[[i]]
  temp_bg_thin <- extract(feature_all, temp_bg_thin) %>% as.data.frame()
  bg_thin_model[[i]] <- temp_bg_thin
}
extraction_end <- Sys.time()

# Create the vector of occ and bg for each species
pa_list <- c()

for (i in 1:length(occlist)) {
  pa_list[[i]] <- c(rep(1, nrow(p_list[[i]])), rep(0, nrow(bg_thin_model[[i]])))
}

# Create dataframes of occ and background points for each species
p_feature <- c()
for (i in 1:length(occlist)) {
  p_feature[[i]] <- rbind(p_list[[i]], bg_thin_model[[i]]) %>% as.data.frame()
}


# Create data frames for ENMevaluate function
occs_list <- c()
bgs_list <-c()
IUCN_list <- c("D:/xinchen/mee_sdm_experiment/0_downloaded_data/4_processed_range_maps/Prionailurus_bengalensis_ESRI102025.shp",
               "D:/xinchen/mee_sdm_experiment/0_downloaded_data/4_processed_range_maps/Zamia_prasina_NAD83_Albers.shp")
IUCN_WGS84_list <- c()

for (i in 1:length(occlist)) {
  
  # Transform projected coordinates under each IUCN range map to WGS84 to meet the requirement of ENMevaluate function
  
  # Presences
  # Read IUCN range map
  IUCN_i <- st_read(IUCN_list[i])
  # Read occurrences with projected (equal area projection) coordinates 
  temp_occ <- read.csv(occlist[i])
  # Convert data frame to sf object
  temp_occ <- st_as_sf(temp_occ, coords = c("x", "y"), na.fail = FALSE, crs = crs(IUCN_i))
  # Project occurrences from equal area projection to WGS84 
  temp_occ <- st_transform(temp_occ, crs(feature_all))
  # Extract lat and long
  temp_lat_long <- as.data.frame(st_coordinates(temp_occ)) %>% dplyr::select(X, Y) %>% rename(c(longitude = "X", latitude = "Y"))
  IUCN_WGS84_list[[i]] <- st_transform(IUCN_i, crs(feature_all))
  # Combine lat, long with extracted environmental values
  occs_list[[i]] <- cbind(temp_lat_long, p_list[[i]]) %>% as.data.frame()
  
  # Background points
  # Read background points with projected (equal area projection) coordinates 
  temp_bg_thin <- read.csv(bg_list[i])
  # Convert data frame to sf object
  temp_bg_thin <- st_as_sf(temp_bg_thin, coords = c("x", "y"), na.fail = FALSE, crs = crs(IUCN_i))
  # Project occurrences from equal area projection to WGS84 
  temp_bg_thin <- st_transform(temp_bg_thin, crs(feature_all))
  # Extract lat and long
  temp_bg_lat_long <- as.data.frame(st_coordinates(temp_bg_thin)) %>% dplyr::select(X, Y) %>% rename(c(longitude = "X", latitude = "Y"))
  # Combine lat, long with extracted environmental values
  bgs_list[[i]] <- cbind(temp_bg_lat_long, bg_thin_model[[i]]) %>% as.data.frame()
}

#######################################################################################################################################
# Part1. Implement ENMevaluate function from 'ENMeval' package
#######################################################################################################################################
# Train Maxent
# Create a list to save Maxent object for extracting coefficients
maxent_list_1st <- vector("list", length(spplist))
names(maxent_list_1st) <- spplist

set.seed(1)
maxent_start <- Sys.time()
for (i in 1:length(spplist)) {
  
  me_i <-  ENMevaluate(occs = occs_list[[i]], bg = bgs_list[[i]], 
                       algorithm = 'maxent.jar', partitions = 'randomkfold',
                       categoricals = "dominant_land_cover", 
                       tune.args = list(fc = c("LQ"), rm = 1:10)
                       )
  
  # Save ith maxent object for making figures or extracting coefficients
  maxent_list_1st[[i]] <- me_i
  
}
maxent_end <- Sys.time()


# Test on the speed of the 2nd species because it has much smaller range
#set.seed(4)
#maxent_start <- Sys.time()

  
  #me_i <-  ENMevaluate(occs = occs_list[[2]][, 1:2], bg = bgs_list[[2]][, 1:2], envs = feature_list[[2]],
  #                     algorithm = 'maxent.jar', partitions = 'randomkfold',
  #                     categoricals = "dominant_land_cover", 
  #                     tune.args = list(fc = c("LQ"), rm = 1:10)
  #)
  
  ## Save ith maxent object for making figures or extracting coefficients
  #maxent_list_1st[[2]] <- me_i
  
#maxent_end <- Sys.time()

#######################################################################################################################################
# Record of the most running logs to track how much time this process needs
# > set.seed(3)
# > maxent_start <- Sys.time()
# > for (i in 1:length(spplist)) {
#   +     
#     +     me_i <-  ENMevaluate(occs = occs_list[[i]][, 1:2], bg = bgs_list[[i]][, 1:2], envs = feature_list[[i]], 
#                                +                          algorithm = 'maxent.jar', partitions = 'checkerboard1',
#                                +                          categoricals = "dominant_land_cover", 
#                                +                          tune.args = list(fc = c("LQ"), rm = 1:10)
#                                +     )
#     +     
#       +     # Save ith maxent object for making figures or extracting coefficients
#       +     maxent_list_1st[[i]] <- me_i
#       +     
#         + }
# *** Running initial checks... ***
#   
#   * Found 609551 raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables.
# * Removed 30 occurrence points with NA predictor variable values.
# * Removed 625 background points with NA predictor variable values.
# * Assigning variable dominant_land_cover to categorical ...
# * Clamping predictor variable rasters...
# * Model evaluations with checkerboard (2-fold) cross validation...
# 
# *** Running ENMeval v2.0.4 with maxent.jar v3.4.3 from dismo package v1.3.14 ***
#   
#   |==============================================================================================================================================================================================| 100%
# ENMevaluate completed in 89 minutes 22.6 seconds.
# *** Running initial checks... ***
#   
#   * Found 14604 raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables.
# * Removed 2 occurrence localities that shared the same grid cell.
# * Removed 2 occurrence points with NA predictor variable values.
# * Removed 108 background points with NA predictor variable values.
# * Assigning variable dominant_land_cover to categorical ...
# * Clamping predictor variable rasters...
# * Model evaluations with checkerboard (2-fold) cross validation...
# 
# *** Running ENMeval v2.0.4 with maxent.jar v3.4.3 from dismo package v1.3.14 ***
#   
#   |==============================================================================================================================================================================================| 100%
# ENMevaluate completed in 4 minutes 11.4 seconds.


# 08/29/2024
# Note: Part2 can be ignored because the best model was determined by best continous Boyce Index, best AUC, and the lowest AICc
#######################################################################################################################################
# Part2. Extract 'MaxEnt' objects with the 0 delta.AICc, best auc.val.avg, and lowest AICc
#######################################################################################################################################
# Date: 08/12/2024
# Change: Crop future bioclimatic rasters using the buffered IUCN range maps (calibration area) because global scale rasters eat a lot of RAM


feature_bio_list <- list.files(path = "D:/xinchen/mee_sdm_experiment/5_R_script", pattern = '\\.tif$', full.names = TRUE)
feature_bio_all <- stack(feature_bio_list)

# Change the extent of each raster layer based on human impact index layer
# Estimated processing time was ~ 2-4 hours
# Estimated RAM required: ~40GB

options(scipen = 99) # reset scipen if ENMevaluate was run 
# Source: https://stackoverflow.com/questions/76889909/support-with-r-error-code-error-in-if-getoptionscipen-mindigits

feature_bio_list <- c()
Sys.time()
for (i in 1:length(occlist)) {
  options(scipen = 999)
  calibration_i <- read_sf(calibration_area_list[i])
  feature_bio_list[[i]] <- crop(feature_bio_all, calibration_i)
  gc()
  feature_bio_list[[i]] <- mask(feature_bio_list[[i]], calibration_i)
  gc()
  writeRaster(feature_bio_list[[i]],
              # a series of names for output files
              filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names(feature_bio_list[[i]]), ".asc"), 
              format = "ascii", ## the output format
              bylayer = TRUE, ## this will save a series of layers
              overwrite = T)
}
Sys.time()

# Model selection
# Species 1
lowest_aic <- eval.results(maxent_list_1st[[1]]) %>%
  filter(delta.AICc == 0)

highest_cbi <- eval.results(maxent_list_1st[[1]]) %>%
  filter(cbi.val.avg == max(cbi.val.avg))

maxent_lowest_aic <- eval.models(maxent_list_1st[[1]])[[lowest_aic$tune.args]]


maxent_highest_cbi <- eval.models(maxent_list_1st[[1]])[[highest_cbi$tune.args]]

# Predict suitability using best AICc
ped1 <- predict(maxent_lowest_aic, feature_list[[1]])  # studyArea is the clipped rasters 
# Plot base map
plot(feature_all$bio1, col = "grey", legend=FALSE, ext=extent(feature_list[[1]]$bio1), main = "Best AICc")
# plot the continuous prediction
plot(ped1, add = TRUE)

# Predict suitability using best CBI
ped2 <- predict(maxent_highest_cbi, feature_list[[1]])  # studyArea is the clipped rasters 
# Plot base map
plot(feature_all$bio1, col = "grey", legend=FALSE, ext=extent(feature_list[[1]]$bio1), main = "Best CBI")
# plot the continuous prediction
plot(ped2, add = TRUE)

# Species 2
# Model selection
lowest_aic_2 <- eval.results(maxent_list_1st[[2]]) %>%
  filter(delta.AICc == 0)

highest_cbi_2 <- eval.results(maxent_list_1st[[2]]) %>%
  filter(cbi.val.avg == max(cbi.val.avg))

maxent_lowest_aic_2 <- eval.models(maxent_list_1st[[2]])[[lowest_aic$tune.args]]


maxent_highest_cbi_2 <- eval.models(maxent_list_1st[[2]])[[highest_cbi$tune.args]]

# Predict suitability using best AICc
ped1_2 <- predict(maxent_lowest_aic, feature_list[[2]])  # studyArea is the clipped rasters 
# Plot base map
plot(feature_all$bio1, col = "grey", legend=FALSE, ext=extent(feature_list[[2]]$bio1), main = "Best AUC")
# plot the continuous prediction
plot(ped1_2, add = TRUE)

# Predict suitability using best CBI
ped2_2 <- predict(maxent_highest_cbi, feature_list[[2]])  # studyArea is the clipped rasters 
# Plot base map
plot(feature_all$bio1, col = "grey", legend=FALSE, ext=extent(feature_list[[2]]$bio1), main = "Best CBI")
# plot the continuous prediction
plot(ped2_2, add = TRUE)
#######################################################################################################################################

# Prepare future layers by cropping future layers with calibration area raster
stack_1 <- c()
stack_2 <- c()
stack_3 <- c()
stack_4 <- c()

# Start
# > Sys.time()
# [1] "2024-08-21 10:24:10 EDT"

# End
# > Sys.time()
# [1] "2024-08-21 15:33:06 EDT"
# This process only makes sure the number of pixels are the same with the mask or calibration area but the spatial extent is bit off
Sys.time()
for (i in 1:length(occlist)) {
  options(scipen = 999)
  calibration_i <- feature_list[[i]]$hii # Use the processed raster for cropping future layers to make sure they have the same spatial extent?

  # ssp585 2041-2070
  feature_list_ssp585_2041_2070 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2041-2070",
                                              pattern = '\\.tif$', full.names = TRUE)
  for (j in 1:length(feature_list_ssp585_2041_2070)) {
    #message(names(feature_list_ssp585_2041_2070[j]))
    message(j)
    raster_j <- raster(feature_list_ssp585_2041_2070[j])
    message(names(feature_list_ssp126_2071_2100[j]))
    names_j <- names(raster_j)
    cells_j <- cellsFromExtent(raster_j, extent(calibration_i))
    raster_j <- raster_j[cells_j, drop = FALSE]
    writeRaster(raster_j,
                # a series of names for output files
                filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names_j, ".asc"), 
                format = "ascii", ## the output format
                bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
    
  }
  # stack_1 <- raster::stack(feature_list_ssp585_2041_2070)
  # stack_1 <- crop(stack_1, calibration_i)
  # gc()
  # stack_1 <- mask(stack_1, calibration_i)
  # writeRaster(stack_1,
  #             # a series of names for output files
  #             filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names(stack_1), ".asc"), 
  #             format = "ascii", ## the output format
  #             bylayer = TRUE, ## this will save a series of layers
  #             overwrite = T)
  # gc()
  
  
  # ssp126 2041-2070
  options(scipen = 999)
  feature_list_ssp126_2041_2070 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2041-2070",
                                              pattern = '\\.tif$', full.names = TRUE)
  for (j in 1:length(feature_list_ssp126_2041_2070)) {
    #message(names(feature_list_ssp126_2041_2070[j]))
    message(j)
    raster_j <- raster(feature_list_ssp126_2041_2070[j])
    message(names(feature_list_ssp126_2071_2100[j]))
    names_j <- names(raster_j)
    cells_j <- cellsFromExtent(raster_j, extent(calibration_i))
    raster_j <- raster_j[cells_j, drop = FALSE]
    writeRaster(raster_j,
                # a series of names for output files
                filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names_j, ".asc"), 
                format = "ascii", ## the output format
                bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
  }
  # stack_2 <- raster::stack(feature_list_ssp126_2041_2070)
  # stack_2 <- crop(stack_2, calibration_i)
  # gc()
  # stack_2 <- mask(stack_2, calibration_i)
  # writeRaster(stack_2,
  #             # a series of names for output files
  #             filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names(stack_2), ".asc"), 
  #             format = "ascii", ## the output format
  #             bylayer = TRUE, ## this will save a series of layers
  #             overwrite = T)
  # gc()
  
  
  # ssp585 2071-2100
  options(scipen = 999)
  feature_list_ssp585_2071_2100 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2071-2100",
                                              pattern = '\\.tif$', full.names = TRUE)
  for (j in 1:length(feature_list_ssp585_2071_2100)) {
    #message(names(feature_list_ssp585_2071_2100[j]))
    message(j)
    raster_j <- raster(feature_list_ssp585_2071_2100[j])
    message(names(feature_list_ssp126_2071_2100[j]))
    names_j <- names(raster_j)
    cells_j <- cellsFromExtent(raster_j, extent(calibration_i))
    raster_j <- raster_j[cells_j, drop = FALSE]
    writeRaster(raster_j,
                # a series of names for output files
                filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names_j, ".asc"), 
                format = "ascii", ## the output format
                bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
  }
  # stack_3 <- raster::stack(feature_list_ssp585_2071_2100)
  # stack_3 <- crop(stack_3, calibration_i)
  # gc()
  # stack_3 <- mask(stack_3, calibration_i)
  # writeRaster(stack_3,
  #             # a series of names for output files
  #             filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names(stack_3), ".asc"), 
  #             format = "ascii", ## the output format
  #             bylayer = TRUE, ## this will save a series of layers
  #             overwrite = T)
  # gc()
  
  # ssp126 2071-2100
  options(scipen = 999)
  feature_list_ssp126_2071_2100 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2071-2100",
                                              pattern = '\\.tif$', full.names = TRUE)
  for (j in 1:length(feature_list_ssp126_2071_2100)) {
    #message(names(feature_list_ssp126_2071_2100[j]))
    message(j)
    raster_j <- raster(feature_list_ssp126_2071_2100[j])
    message(names(feature_list_ssp126_2071_2100[j]))
    names_j <- names(raster_j)
    cells_j <- cellsFromExtent(raster_j, extent(calibration_i))
    raster_j <- raster_j[cells_j, drop = FALSE]
    writeRaster(raster_j,
                # a series of names for output files
                filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names_j, ".asc"), 
                format = "ascii", ## the output format
                bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
  }
  # stack_4 <- raster::stack(feature_list_ssp126_2071_2100)
  # stack_4 <- crop(stack_4, calibration_i)
  # gc()
  # stack_4 <- mask(stack_4, calibration_i)
  # writeRaster(stack_4,
  #             # a series of names for output files
  #             filename = paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names(stack_4), ".asc"), 
  #             format = "ascii", ## the output format
  #             bylayer = TRUE, ## this will save a series of layers
  #             overwrite = T)
  # gc()
}
Sys.time()

# Modify climate layers because the future climate layers have tiny different spatial extent
Sys.time()
for (i in 1:length(occlist)) {
  # Track the loop
  message(i)
  print(spplist[i])
  
  # List the raster layers
  feature_list_i <- list.files(paste0("Output/future_calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/"), pattern = ".asc$", full.names = TRUE)
  
  # Read one raster with correct spatial extent of all other non-climate rater layers
  calibration_i <- feature_list[[i]]$hii
  
  for (j in 1:length(feature_list_i)) {
    message(j)
    # Read raster layer
    raster_j <- raster(feature_list_i[j])
    names_j <- names(raster_j)
    raster_j <- crop(raster_j, calibration_i)
    raster_j <- raster(vals=values(raster_j), ext = extent(calibration_i), crs=crs(calibration_i),
                       nrows=dim(calibration_i)[1],ncols=dim(calibration_i)[2])
    raster_j <- mask(raster_j, calibration_i) # Unfortunately, raster::mask function requires same spatial extent if the mask is in raster object format
    writeRaster(raster_j,
                # a series of names for output files
                filename = paste0("Output/future_calibration_raster/Adjusted_1/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", names_j, ".asc"),
                format = "ascii", ## the output format
                #bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
  }
}
Sys.time()


# Predict suitability using the final Maxent for each species
# Model selection based on continuous Boyce Index (CBI), AICc, and AUC
Maxent_final_list <- c()
highest_cbi_list <- c()

for (i in 1:length(occlist)) {
  # Select the evaluation metrics based on the highest continuous Boyce Index (CBI) first
  highest_cbi_list[[i]] <- eval.results(maxent_list_1st[[i]]) %>% filter(cbi.val.avg == max(cbi.val.avg))
  
  # Subset the corresponding Maxent with the highest CBI
  Maxent_final_list[[i]] <- eval.models(maxent_list_1st[[i]])[[highest_cbi_list[[i]]$tune.args]]
  
}

# Create the index for each GCM
gcm_list <- c()
for (i in 1:5) {
  gcm_list[[i]] <- sapply(1:6, FUN = function(X)(i+5*(X-1)))
}

# Create a list of GCM names
gcm_name_list <- c("gfdl.esm4", "ipsl.cm6a.lr", "mpi.esm1.2.hr", "mri.esm2.0", "ukesm1.0.ll")

# Making prediction for present
Sys.time()
for (i in 1:length(occlist)) {
  feature_list_present_i <- feature_list[[i]]
  feature_list_present_i$dominant_land_cover <- as.factor(feature_list_present_i$dominant_land_cover)
  pred_present_i <- predict(Maxent_final_list[[i]], feature_list_present_i)
  writeRaster(pred_present_i,
              # a series of names for output files
              filename = paste0("Output/present_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", "suitability_present", ".tif"),
              format = "GTiff", ## the output format
              #bylayer = TRUE, ## this will save a series of layers
              options=c("COMPRESS=NONE", "TFW=YES"),
              overwrite = T)
}
Sys.time()

# Making prediction for future scenarios
for (i in 1:length(occlist)) {
  # List non-climate and climate layers from feature list for species i 
  feature_list_i <- list.files(paste0("Output/calibration_raster/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i]))), pattern = ".asc$", full.names = TRUE)[-c(1:6)]
  feature_list_bio_i <- list.files(paste0("Output/future_calibration_raster/Adjusted/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/"), pattern = ".asc$", full.names = TRUE)
  
  # Subset climate layers for each future scenario
  bio_stack_2041_2070_ssp126_index_i <- grep("2041.+2070.+ssp126", feature_list_bio_i)
  bio_stack_2041_2070_ssp585_index_i <- grep("2041.+2070.+ssp585", feature_list_bio_i)
  bio_stack_2071_2100_ssp126_index_i <- grep("2071.+2100.+ssp126", feature_list_bio_i)
  bio_stack_2071_2100_ssp585_index_i <- grep("2071.+2100.+ssp585", feature_list_bio_i)
    
  # Loop all future scenarios: 2 SSP x 2 time period
  for (j in 1:5) {
    # 2041-2070 ssp126
    feature_list_1_j <- raster::stack(c(feature_list_bio_i[bio_stack_2041_2070_ssp126_index_i[gcm_list[[j]]]], feature_list_i))
    names(feature_list_1_j) <- names(feature_list[[i]])
    feature_list_1_j$dominant_land_cover <- raster::as.factor(feature_list_1_j$dominant_land_cover)
    pred_1_j <- predict(Maxent_final_list[[i]], feature_list_1_j)
    
    writeRaster(pred_1_j,
                # a series of names for output files
                filename = paste0("Output/future_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", gcm_name_list[j], "2041_2070_ssp126_pred", ".asc"),
                format = "ascii", ## the output format
                #bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
    
    # 2041-2070 ssp585
    feature_list_2_j <- raster::stack(c(feature_list_bio_i[bio_stack_2041_2070_ssp585_index_i[gcm_list[[j]]]], feature_list_i))
    names(feature_list_2_j) <- names(feature_list[[i]])
    feature_list_2_j$dominant_land_cover <- raster::as.factor(feature_list_2_j$dominant_land_cover)
    pred_2_j <- predict(Maxent_final_list[[i]], feature_list_2_j)
    
    writeRaster(pred_2_j,
                # a series of names for output files
                filename = paste0("Output/future_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", gcm_name_list[j], "2041_2070_ssp585_pred", ".asc"),
                format = "ascii", ## the output format
                #bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
    
    # 2071-2100 ssp126
    feature_list_3_j <- raster::stack(c(feature_list_bio_i[bio_stack_2071_2100_ssp126_index_i[gcm_list[[j]]]], feature_list_i))
    names(feature_list_3_j) <- names(feature_list[[i]])
    feature_list_3_j$dominant_land_cover <- raster::as.factor(feature_list_3_j$dominant_land_cover)
    pred_3_j <- predict(Maxent_final_list[[i]], feature_list_3_j)
    
    writeRaster(pred_3_j,
                # a series of names for output files
                filename = paste0("Output/future_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", gcm_name_list[j], "2071_2100_ssp126_pred", ".asc"),
                format = "ascii", ## the output format
                #bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
    
    # 2071-2100 ssp585
    feature_list_4_j <- raster::stack(c(feature_list_bio_i[bio_stack_2071_2100_ssp585_index_i[gcm_list[[j]]]], feature_list_i))
    names(feature_list_4_j) <- names(feature_list[[i]])
    feature_list_4_j$dominant_land_cover <- raster::as.factor(feature_list_4_j$dominant_land_cover)
    pred_4_j <- predict(Maxent_final_list[[i]], feature_list_4_j)
    
    writeRaster(pred_4_j,
                # a series of names for output files
                filename = paste0("Output/future_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/", gcm_name_list[j], "2071_2100_ssp585_pred", ".asc"),
                format = "ascii", ## the output format
                #bylayer = TRUE, ## this will save a series of layers
                overwrite = T)
  }
  
}



# Species 1
#lowest_aic <- eval.results(maxent_list_1st[[1]]) %>%
#  filter(delta.AICc == 0)

#highest_cbi <- eval.results(maxent_list_1st[[1]]) %>%
#  filter(cbi.val.avg == max(cbi.val.avg))

#maxent_lowest_aic <- eval.models(maxent_list_1st[[1]])[[lowest_aic$tune.args]]


#maxent_highest_cbi <- eval.models(maxent_list_1st[[1]])[[highest_cbi$tune.args]]




