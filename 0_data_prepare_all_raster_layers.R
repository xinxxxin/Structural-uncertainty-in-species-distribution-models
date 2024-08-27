## This script was used to model all species using Maxent and investigate response curve
## The occurrence data was downloaded from GBIF and then only the occurrences within 1.1 radius of IUCN range maps between 1980-2020 were selected
## The background points were generated within within 1.1 radius of IUCN range maps and then thinned to make sure one bg point per pixel
## The bioclimatic variables were downloaded from CHELSA 2.1


# There are two main parts of this script:
# part1
# 1. Making shapefiles based on csv and save them in a vector
# 2. Generating 1.1 radius of IUCN maps
# 3. Spatial thinning and then export thinned occurrence as csv files
# 4. Generate background points within 1.1 radius of IUCN maps and spatial thinning for background points

# part2
# 5. Extracting environmental factors for both occurrence and background points
# 6. Running Maxent for each species
# 7. Extracting information using rmaxent package

# Load packages
library(data.table)
library(sp)
library(mapview) # This is very nice for quick visualization

# library(rangeBuilder)
library(spThin)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(raster)
library(dismo)

#####################################################################################################################################################
# Part1. Making shapefiles based on csv files, project the shapefile to ESRI 102025 Albers Northern Asia, remove occurrences with NA in environmental 
# layer, and then save each shapefile in a vector
# Part2. Generate 1.1 radius of IUCN range maps
# Part3. Spatial thinning for occurrence data
# Part4. Generate background points within 1.1 radius of IUCN range maps and spatial thinning for background points
#####################################################################################################################################################
# Load environmental data first for getting crs and removing occurrences with NA in those layers
# Load data
feature_list <- list.files(path = "D:/xinchen/mee_sdm_experiment/4_processed_raster_layers",
                           pattern = '\\.tif$', full.names = TRUE)

# Regular expression: https://stackoverflow.com/questions/4876813/using-r-to-list-all-files-with-a-specified-extension

# There is a little off for the extents of some raster layers: 17042x43200, 17040x43200, and 17041x43200
# Cut the land cover and bioclimatic variables based on human impact index layer since it has the smallest dimension
hii <- raster(feature_list[12])
feature_all <- stack(hii)

# Change the extent of each raster layer based on human impact index layer
# Estimated processing time was ~ 2-4 hours
# Estimated RAM required: ~40GB

Sys.time()
for (i in 1:13) {
  # Track the loop
  print(i)
  
  # Read raster layer
  raster_i <- raster(feature_list[-12][i])
  
  # Change the extent of each rater based on hii
  raster_i <- crop(raster_i, extent(hii))
  raster_i <- raster(vals=values(raster_i), ext=extent(hii), crs=crs(hii),
                     nrows=dim(hii)[1],ncols=dim(hii)[2])
  
  # Stack each raster to feature_all
  feature_all <- stack(feature_all, raster_i)
  
  # Check the number of layers in the stack
  nlayers(feature_all)
}
Sys.time()

# Name the raster stack
names(feature_all) <- c("hii", "bio1", "bio10", "bio11", "bio12", "bio16", "bio17", "cropland", "dominant_land_cover", "elevation", "grassland", "haerbaceous_covered_area", 
                        "shrubs_covbered_area", "trees_covered_area")

# Export each raster in raster stack 
Sys.time()
for (i in 1:14) {
  
  # If the raster layers need to be visualized in ArcGIS, options of TFW should be set to YES
  # Source: https://terpconnect.umd.edu/~egurarie/research/NWT/Step02b_ExportingRastersToArcGIS.html
  # Example: writeRaster(PEM60RAT, filename = "PEM60RAT.tif", options=c('TFW=YES'))
  writeRaster(feature_all[[i]], paste0("processed_raster_layers/", names(feature_all)[i] ,".tif"))
  gc()
}
Sys.time()

#saveRDS(feature_all, "1_feature_all.rds")

# Create CRS for IUCN range maps
crs_list <- list(Prionailurus_bengalensis = rep(NA, 1), Zamia_prasina = rep(NA, 1))
crs_list$Prionailurus_bengalensis <-  crs(read_sf("D:/xinchen/mee_sdm_experiment/0_downloaded_data/4_processed_range_maps/Prionailurus_bengalensis_ESRI102025.shp"))
crs_list$Zamia_prasina <- crs(read_sf("D:/xinchen/mee_sdm_experiment/0_downloaded_data/4_processed_range_maps/Zamia_prasina_NAD83_Albers.shp"))

# Export crs information
saveRDS(crs_list, "crs_list.rds")
