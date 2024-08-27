# This script was used to change the extent of 6 bioclimatic variables of 5 different GCMs by 2 SSPs from CHELSA

# Load packages
library(raster)

# There is a little off for the extents of some raster layers: 17042x43200, 17040x43200, and 17041x43200
# Cut the land cover and bioclimatic variables based on human impact index layer since it has the smallest dimension
hii <- raster(list.files(path = "D:/xinchen/mee_sdm_experiment/4_processed_raster_layers",
                         pattern = '\\.tif$', full.names = TRUE)[12])

feature_list <- list.files(path = "D:/xinchen/mee_sdm_experiment/5_R_script", pattern = '\\.tif$', full.names = TRUE)
feature_all <- stack(feature_list)

# Change the extent of each raster layer based on human impact index layer
# Estimated processing time was ~ 2-4 hours
# Estimated RAM required: ~40GB

Sys.time()
for (i in 1:length(names(feature_all))) {
  # Track the loop
  message(names(feature_all)[i])
  
  # Read raster layer
  raster_i <- feature_all[[i]]
  
  # Change the extent of each rater based on hii
  raster_i <- crop(raster_i, extent(hii))
  raster_i <- raster(vals=values(raster_i), ext=extent(hii), crs=crs(hii),
                     nrows=dim(hii)[1],ncols=dim(hii)[2])
  
  # If the raster layers need to be visualized in ArcGIS, options of TFW should be set to YES
  # Source: https://terpconnect.umd.edu/~egurarie/research/NWT/Step02b_ExportingRastersToArcGIS.html
  # Example: writeRaster(PEM60RAT, filename = "PEM60RAT.tif", options=c('TFW=YES'))
  writeRaster(raster_i, paste0("processed_future_raster_layers/", names(feature_all)[i], ".tif"))
  gc()
}
Sys.time()