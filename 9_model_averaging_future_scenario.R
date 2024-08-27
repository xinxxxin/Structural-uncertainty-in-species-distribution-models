# This script was used to average the suitablity maps of future scenarios for each species

# Load packages
library(terra)

# Create lists for different naming purpose
occlist <- list.files("Output/occurrence/", 
                      pattern = ".csv", full.names = TRUE)

# Average model predictions of future scenarios for each species
Sys.time()
for (i in 1:length(occlist)) {
  message(i)
  # list suitability maps under each scenarios
  map_i <- list.files(paste0("Output/future_prediction/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])), "/"), pattern = ".asc$", full.names = TRUE)
  output_path_i <- paste0("Output/future_prediction/suitability_maps/", sub("*.csv$", "", sub(".*Output/occurrence/", "", occlist[i])))
  
  # 2041-2070 ssp126
  map_stack_2041_2070_ssp126_index_i <- grep("2041.+2070.+ssp126", map_i)
  map_stack_2041_2070_ssp126_i <- terra::rast(map_i[map_stack_2041_2070_ssp126_index_i])
  mean_2041_2070_ssp126_i <- terra::mean(map_stack_2041_2070_ssp126_i)
  
  # 2041-2070 ssp585
  map_stack_2041_2070_ssp585_index_i <- grep("2041.+2070.+ssp585", map_i)
  map_stack_2041_2070_ssp585_i <- terra::rast(map_i[map_stack_2041_2070_ssp585_index_i])
  mean_2041_2070_ssp585_i <- terra::mean(map_stack_2041_2070_ssp585_i)
  
  # 2071-2100 ssp126
  map_stack_2071_2100_ssp126_index_i <- grep("2071.+2100.+ssp126", map_i)
  map_stack_2071_2100_ssp126_i <- terra::rast(map_i[map_stack_2071_2100_ssp126_index_i])
  mean_2071_2100_ssp126_i <- terra::mean(map_stack_2071_2100_ssp126_i)
  
  # 2071-2100 ssp585
  map_stack_2071_2100_ssp585_index_i <- grep("2071.+2100.+ssp585", map_i)
  map_stack_2071_2100_ssp585_i <- terra::rast(map_i[map_stack_2071_2100_ssp585_index_i])
  mean_2071_2100_ssp585_i <- terra::mean(map_stack_2071_2100_ssp585_i)

  # Average different ssp for each time period
  mean_2041_2070_i <- terra::mean(c(mean_2041_2070_ssp126_i, mean_2041_2070_ssp585_i))
  mean_2071_2100_i <- terra::mean(c(mean_2071_2100_ssp126_i, mean_2071_2100_ssp585_i))
  
  # Export each raster layers
  # 2041-2070 ssp126
  terra::writeRaster(mean_2041_2070_ssp126_i, 
                     file.path(output_path_i, "suitability_2041_2070_ssp126.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
  # 2041-2070 ssp585
  terra::writeRaster(mean_2041_2070_ssp585_i, 
                     file.path(output_path_i, "suitability_2041_2070_ssp585.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
  # 2071-2100 ssp126
  terra::writeRaster(mean_2071_2100_ssp126_i, 
                     file.path(output_path_i, "suitability_2071_2100_ssp126.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
  # 2071-2100 ssp585
  terra::writeRaster(mean_2071_2100_ssp585_i, 
                     file.path(output_path_i, "suitability_2071_2100_ssp585.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
  # 2041-2070 mean
  terra::writeRaster(mean_2041_2070_i, 
                     file.path(output_path_i, "suitability_2041_2070_mean.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
  # 2071-2100 mean
  terra::writeRaster(mean_2071_2100_i, 
                     file.path(output_path_i, "suitability_2071_2100_mean.tif"), 
                     overwrite = TRUE, 
                     gdal=c("COMPRESS=NONE", "TFW=YES") # They require GeoTiff file format
                     )
  
}
Sys.time()