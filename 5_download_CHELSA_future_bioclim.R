# This script was used to download 6 bioclimatic variables of 5 different GCMs by 2 SSPs from CHELSA

# Load package
library(raster)
library(terra)
library(RCurl)

# Adjust timeout in case downloading was interrupted
options(timeout=10000000000000)


# Read URLs
future_ssp585_2041_2070 <- read.table("5_envidatS3paths_ssp585_2041_2070.txt", quote="\"", comment.char="")
future_ssp585_2071_2100 <- read.table("5_envidatS3paths_ssp585_2071_2100.txt", quote="\"", comment.char="")

future_ssp126_2041_2070 <- read.table("5_envidatS3paths_ssp126_2041_2070.txt", quote="\"", comment.char="")
future_ssp126_2071_2100 <- read.table("5_envidatS3paths_ssp126_2071_2100.txt", quote="\"", comment.char="")

# Download data
# 2041-2070
Sys.time()
# ssp585
for (i in 1:30) {
  message(i)
  download.file(future_ssp585_2041_2070[i,], 
                destfile = paste0("D://xinchen//mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2041-2070/", 
                                  sub(".*/ssp585/bio/", "", future_ssp585_2041_2070[i, ])),
                method = "libcurl",
                mode = "wb" # For Windows OS
                )
  gc()
  }

# ssp126
for (i in 1:30) {
  message(i)
  download.file(future_ssp126_2041_2070[i,], 
                destfile = paste0("D://xinchen//mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2041-2070/",
                                  sub(".*/ssp126/bio/", "", future_ssp126_2041_2070[i, ]), sep = ""),
                method = "libcurl",
                mode = "wb" # For Windows OS
                )
  gc()
  }

Sys.time()


# 2071-2100
Sys.time()
# ssp585
for (i in 1:30) {
  message(i)
  download.file(future_ssp585_2071_2100[i,], 
                destfile = paste0("D://xinchen//mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2071-2100/", 
                                  sub(".*/ssp585/bio/", "", future_ssp585_2071_2100[i, ])),
                method = "libcurl",
                mode = "wb" # For Windows OS
                )
  gc()
  }

# ssp126
for (i in 1:30) {
  message(i)
  download.file(future_ssp126_2071_2100[i,], 
                destfile = paste0("D://xinchen//mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2071-2100/", 
                                  sub(".*/ssp126/bio/", "", future_ssp126_2071_2100[i, ])),
                method = "libcurl",
                mode = "wb" # For Windows OS
  )
  gc()
  }
Sys.time()

# Create a list for bioclmatic names
name_list <- c("bio1", "bio10", "bio11", "bio12", "bio16", "bio17")

# Average 6 bioclimatic variables for each SSP
# ssp585 2041-2070
feature_list_ssp585_2041_2070 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2041-2070",
                                            pattern = '\\.tif$', full.names = TRUE)

for (i in 1:6) {
  message(paste0(name_list[i], "_2041_2070_ssp585.tif"))
  message(i)
  file_list_i <- feature_list_ssp585_2041_2070[(1+5*(i-1)):(5*i)]
  stack_i <- raster::stack(file_list_i)
  r_avg_i <- terra::mean(stack_i)
  writeRaster(r_avg_i, paste0(name_list[i], "_2041_2070_ssp585.tif"))
  gc()
  }

# ssp126 2041-2070
feature_list_ssp126_2041_2070 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2041-2070",
                                            pattern = '\\.tif$', full.names = TRUE)

for (i in 1:6) {
  message(paste0(name_list[i], "_2041_2070_ssp126.tif"))
  file_list_i <- feature_list_ssp126_2041_2070[(1+5*(i-1)):(5*i)]
  stack_i <- raster::stack(file_list_i)
  r_avg_i <- terra::mean(stack_i)
  writeRaster(r_avg_i, paste0(name_list[i], "_2041_2070_ssp126.tif"))
  gc()
  }

# ssp585 2071-2100
feature_list_ssp585_2071_2100 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp585/2071-2100",
                                            pattern = '\\.tif$', full.names = TRUE)

for (i in 1:6) {
  message(paste0(name_list[i], "_2071_2100_ssp585.tif"))
  file_list_i <- feature_list_ssp585_2041_2070[(1+5*(i-1)):(5*i)]
  stack_i <- raster::stack(file_list_i)
  r_avg_i <- terra::mean(stack_i)
  writeRaster(r_avg_i, paste0(name_list[i], "_2071_2100_ssp585.tif"))
  gc()
  }

# ssp126 2071-2100
feature_list_ssp126_2071_2100 <- list.files(path = "D:/xinchen/mee_sdm_experiment/0_downloaded_data/2_future_ssp_scenarios/ssp126/2071-2100",
                                            pattern = '\\.tif$', full.names = TRUE)

for (i in 1:6) {
  message(paste0(name_list[i], "_2071_2100_ssp126.tif"))
  file_list_i <- feature_list_ssp126_2041_2070[(1+5*(i-1)):(5*i)]
  stack_i <- raster::stack(file_list_i)
  r_avg_i <- terra::mean(stack_i)
  writeRaster(r_avg_i, paste0(name_list[i], "_2071_2100_ssp126.tif"))
  gc()
}



