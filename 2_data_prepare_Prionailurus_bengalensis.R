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

library(sm)

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
  writeRaster(feature_all[[i]], paste0("processed_raster_layers/", names(feature_all)[i] ,".tif"))
  gc()
}
Sys.time()

saveRDS(feature_all, "1_feature_all.rds")

# Load GBIF data for Prionailurus bengalensis
Sys.time()
occurrence <- fread("D:/xinchen/mee_sdm_experiment/0_downloaded_data/7_GBIF_data/0054460-240626123714530/occurrence.txt",
                    select = c("decimalLatitude", "decimalLongitude", "year",
                               "genericName", "specificEpithet", "order", "basisOfRecord"))
Sys.time()

# # Split occurrence for 1990 to 2020
# occurrence <- occurrence[occurrence$year >= 1980& occurrence$year <= 2020, ]
# 
# # Create a column using genericName+specificEpithet
# occurrence$gbif_name <- paste0(occurrence$genericName, paste0(" ", occurrence$specificEpithet))

# Keep basisOfRecord of "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE", "MATERIAL_CITATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE".
# For "PRESERVED_SPECIMEN", even if a record could be outside of species native range, https://data-blog.gbif.org/post/living-specimen-to-preserved-specimen-understanding-basis-of-record/
# However, our filtering method using IUCN range maps will remove the specimen that fall outside of species range.

# Source: https://docs.gbif.org/course-data-use/en/basis-of-record.html
occurrence <- occurrence[occurrence$basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE", 
                                                         "MATERIAL_CITATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE"), ]

# Check valid lat and long for occurrence
nrow(occurrence[is.na(occurrence$decimalLatitude)])

# Remove the occurrences without valid coordinates
occurrence <- occurrence[!is.na(occurrence$decimalLatitude), ]

# Double check genericName and specificEpithet
table(occurrence$specificEpithet)
table(occurrence$genericName)

# Based on GBIF species name: https://www.gbif.org/species/2434903
# The genericName can only be either Prionailurus or Felis and the specificEpithet can only be bengalensis
# There are some wield unranked species names such as BOLD:AAH7627, BOLD:ACV6153, BOLD:ADR3171, BOLD:ADR9017.
occurrence <- occurrence[occurrence$genericName %in% c("Felis", "Prionailurus"), ]
occurrence <- occurrence[occurrence$specificEpithet %in% c("bengalensis"), ]

# Create a column with checked GBIF name
occurrence$gbif_name <- "Prionailurus_bengalensis"

# Load dissolved IUCN range maps
IUCN_all <- read_sf("D:/xinchen/mee_sdm_experiment/0_downloaded_data/4_processed_range_maps/Prionailurus_bengalensis_ESRI102025.shp")

# Load new matching table for connecting GBIF name with sci_name in IUCN maps
IUCN_species_list_9_matching <- data.frame(sci_name = "Prionailurus bengalensis",
                                           gbif_name = "Prionailurus_bengalensis")

# Create a column for naming output files
# Date: 04/08/2024
# Code review: sort(as.vector(IUCN_species_list_9_matching$gbif_name) caused mismatch between gbif_name and created GBIF_name after check 
# WRONG CODE EXAMPLE: IUCN_species_list_9_matching$GBIF_name <- gsub(" ", "_", sort(as.vector(IUCN_species_list_9_matching$gbif_name)))
IUCN_species_list_9_matching$GBIF_name <- gsub(" ", "_", as.vector(IUCN_species_list_9_matching$gbif_name))

# Create spplist using the GBIF_name column in IUCN_species_list_9_matching because they represent >= 90% range intersecting North America
spplist <- IUCN_species_list_9_matching$GBIF_name

# Create the dataframe to save the data
summary_mcp_ratio <- data.frame(sci_name = IUCN_species_list_9_matching$sci_name,
                                number_of_occurrence = rep(NA, nrow(IUCN_species_list_9_matching)),
                                GBIF_name = IUCN_species_list_9_matching$gbif_name,
                                area_IUCN = rep(NA, nrow(IUCN_species_list_9_matching)),
                                area_IUCN_1.1buffer = rep(NA, nrow(IUCN_species_list_9_matching)),
                                area_mcp = rep(NA, nrow(IUCN_species_list_9_matching)),
                                ratio_mcp = rep(NA, nrow(IUCN_species_list_9_matching)))


# Define the location of projected, NA removed, and spatially thinned occurrence
occ_output <- "D:/xinchen/mee_sdm_experiment/5_R_script/Output/Prionailurus_bengalensis/occurrence/"
bg_output <- "D:/xinchen/mee_sdm_experiment/5_R_script/Output/Prionailurus_bengalensis/bg"
cali_area <- "D:/xinchen/mee_sdm_experiment/5_R_script/Output/Prionailurus_bengalensis/calibration_area/"

# list_run is used to control the looped species and locate the species that may take a lot of RAM for generating buffers
list_run <- c(1:length(spplist))

# Loop thru spplist
part1_start <- Sys.time()
for (i in list_run) {
  # Dynamically track species in the spplist
  message(i)
  
  # Find the matching GBIF species name
  temp_occ_GBIF_species_name <- spplist[i]
  
  # Extract occurrence from downloaded occurrences data based on a time window
  temp_occ <- as.data.frame(occurrence[occurrence$gbif_name == temp_occ_GBIF_species_name&
                                         occurrence$year >= 1980& occurrence$year <= 2024
                                       , c("gbif_name", "decimalLongitude", "decimalLatitude")])
  
  # Remove NAs
  temp_occ <- na.omit(temp_occ)
  
  temp_occ$decimalLongitude <- as.numeric(temp_occ$decimalLongitude)
  temp_occ$decimalLatitude <- as.numeric(temp_occ$decimalLatitude)
  
  # Remove duplicated records because they affect making alpha shape, it is surprising that GBIF data has so many duplicated coordinates
  temp_occ <- temp_occ[!duplicated(temp_occ), ]
  
  # Calculate the number of occurrence
  summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i], "number_of_occurrence"] <- nrow(temp_occ)
  
  # Create shapefile for occurrence
  # Find the corresponding IUCN range map for occurrence
  # Create 1.1 radius of IUCN range map
  if (nrow(temp_occ) > 0){
    # Create spatial data frame
    coordinates(temp_occ) <- ~decimalLongitude + decimalLatitude
    
    # Define CRS of occurrence
    proj4string(temp_occ) <- CRS("+init=EPSG:4326")
    
    # Project occurrence onto the coordinate system of IUCN range maps
    temp_occ <- st_transform(st_as_sf(temp_occ), crs(IUCN_all))
    
    # Find the range map of the species
    IUCN_temp_occ <- st_make_valid(IUCN_all)
    
    # Date: 08/07/2024
    # Dissolve self-provided range maps if the range maps are not from the dissolved IUCN range mpas dataset
    IUCN_temp_occ$gbif_name <- spplist[i]
    # sf dissolve using dplyr pipepline approach, source: https://stackoverflow.com/questions/69351968/arcgis-dissolve-equivalent-in-r
    # Or, https://github.com/r-spatial/sf/issues/290
    
    # Dynamically visuallize the processed range map
    IUCN_temp_occ %>% 
      group_by(gbif_name) %>% 
      summarise() %>% 
      st_cast() %>% 
      mapview()
    
    IUCN_temp_occ <- IUCN_temp_occ %>% group_by(gbif_name) %>% summarize() %>% st_cast()
    
    # Calculate the radius. This part has been updated because the self-provided range maps may not be dissolved for each species.
    # Sum over the area returned by st_area function to obtain the total area of a species range
    # IUCN_temp_occ_r <- sqrt((st_area(IUCN_temp_occ)/pi))
    IUCN_temp_occ_r <- sqrt((sum(st_area(IUCN_temp_occ))/pi))
    
    # Buffer IUCN range map by 1.1 radius
    IUCN_temp_occ_1.1_buffer <- st_buffer(IUCN_temp_occ, IUCN_temp_occ_r*0.1)
    IUCN_temp_occ_1.1_buffer <- st_transform(IUCN_temp_occ_1.1_buffer, crs = crs(IUCN_all)) # Make sure it under the same CRS of IUCN_all
    #summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i], "area_IUCN_1.1buffer"] <- st_area(IUCN_temp_occ_1.1_buffer)
    summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i], "area_IUCN_1.1buffer"] <- sum(st_area(IUCN_temp_occ_1.1_buffer))
    
    # Date: 08/12/2024
    # Change: Export buffered IUCN range map for cropping environmental layers during modeling
    # Project the buffered IUCN range map to the CRS of feature_all
    IUCN_temp_occ_1.1_buffer_out <- st_transform(IUCN_temp_occ_1.1_buffer, crs = crs(feature_all))
    st_write(IUCN_temp_occ_1.1_buffer_out, dsn = cali_area, layer = paste0(spplist[i], "_calibration_area.shp", sep = ""), driver="ESRI Shapefile")
    
    # Keep the occurrence only in 1.1 buffered IUCN
    temp_occ_in_1.1_buffer <- st_as_sf(temp_occ)[(st_intersects(st_as_sf(temp_occ), IUCN_temp_occ_1.1_buffer, sparse = F))[, 1],]
    
    if (nrow(temp_occ_in_1.1_buffer) != 0){
      # Spatial thinning for occurrence
      # Extract raster values from environmental layers
      # This step makes sure that the occurrences would not have NAs in all feature layers
      # Note: The vector data will be projected onto the crs of raster data, which is convenient
      temp_occ_feature_i <- raster::extract(feature_all, temp_occ_in_1.1_buffer, cellnumbers = TRUE)
      
      # Remove the occurrences with NA in environmental layers
      # Remember 1st column is cellnumber
      temp_occ <- temp_occ_in_1.1_buffer[temp_occ_feature_i[, 1] %in% na.omit(temp_occ_feature_i)[, 1], ]
      
      # nrow(temp_occ) can be > 1, =1, or = 0
      # Only nrow(temp_occ) > 1 do spatial thinning and output occurrence
      
      # Check if there is only one record after NA removed
      if (nrow(temp_occ) > 1) {
        
        # Convert occ back to spatialdataframe to make gridSample function work
        temp_occ <- sf::as_Spatial(temp_occ)
        
        # Create a RasterLayer with the extent of occ
        r_10km_i <- raster(temp_occ)
        
        # Set the resolution of the cells to 1km
        res(r_10km_i) <- 50000
        
        # Expand (extend) the extent of the RasterLayer a little
        # extend the boundary by 10 times of resolution in case some occurrences on the boundary influence gridSample function
        r_10km_i <- extend(r_10km_i, extent(r_10km_i) + 50000*10)
        
        # Sample
        temp_occ <- gridSample(temp_occ, r_10km_i, n = 1) %>% as.data.frame()
        names(temp_occ) <- c("x", "y")
        
        # # Check the spatially thinned occurrences
        # p <- rasterToPolygons(r_10km_i)
        # plot(p, border='gray')
        # points(temp_occ, cex=1, col='red', pch='x')
        
      } else if (nrow(temp_occ) == 1){
        # For the species with only one occurrence within IUCN range map
        # Still send the one occurrence to temp_occ but it will be filtered out by the threshold (e.g. >= 10) for the number of occurrence later
        temp_occ <- as.data.frame(st_coordinates(temp_occ)) %>% dplyr::select(X, Y)
        names(temp_occ) <- c("x", "y")
      }
      
      # Decide the number of occurrences after spatial thinning
      # Source: https://onlinelibrary.wiley.com/doi/10.1111/ecog.01509
      # " Model performance is known to rapidly decrease for sample sizes smaller than 20 (Stockwell and Peterson 2002) or 
      # 15 (Papeş and Gaubert 2007), and is dramatically poor for samples sizes smaller than 5 records (Pearson et al. 2007). 
      if (nrow(temp_occ) >= 10) {
        # Export occ as csv files for background points sampling
        write.csv(temp_occ, paste0(occ_output, "/", spplist[i], ".csv", sep = ""), row.names = FALSE)
        
        # Generate bg points only for the species with occurrence >= 10
        # Generate 10000 background points
        # set.seed(1)
        
        # Date: 1/30/2024
        # Change: Use terra::spatSample or sf::st_sample to draw bg points using IUCN map because sampleRandom needs IUCN map to be converted to raster
        # Note: Current SDM tutorial uses terra:spatSample https://rspatial.org/sdm/3_sdm_absence-background.html
        #       Compared with old version using 'raster' package: https://rspatial.org/raster/sdm/2_sdm_occdata.html#sampling-bias
        
        # # Convert sf feature to raster to make sampleRandom function run
        # # st_sample sometimes cannot fully create random samples across the given spatial data
        # IUCN_temp_occ_1.1_buffer_raster <- raster::raster(crs = st_crs(IUCN_temp_occ_1.1_buffer), 
        #                                                   vals = 0, resolution = c(1000, 1000), 
        #                                                   ext = extent(IUCN_temp_occ_1.1_buffer)) %>% terra::rasterize(IUCN_temp_occ_1.1_buffer, .)
        # 
        # # sampleRandom function cannot generate random points >= ncell, creating the number of background points = ncell does not make sense
        # # because it would be obsence everywhere. Therefore, thinning background points using 10km as same as that for thinning occurrences.
        # if (ncell(IUCN_temp_occ_1.1_buffer_raster) >= 10000){
        #   
        #   # Create 10000 bg points if ncell >= 10000
        #   temp_bg <- sampleRandom(x = IUCN_temp_occ_1.1_buffer_raster,
        #                           size = 10000,
        #                           na.rm = T, # removes the 'Not Applicable' points
        #                           sp = T) # return spatial points 
        # } else {
        #   
        #   # Create bg points = ncell
        #   temp_bg <- sampleRandom(x = IUCN_temp_occ_1.1_buffer_raster,
        #                           size = ncell(IUCN_temp_occ_1.1_buffer_raster),
        #                           na.rm = T, # removes the 'Not Applicable' points
        #                           sp = T) # return spatial points 
        # }
        
        # Spatial thinning for background points using 10km because some species have bg points in every pixel
        # # Project background points onto NDA 83 Albers
        # temp_bg <- st_transform(st_as_sf(temp_bg), crs(feature_all))
        
        # Date: 08/13/2024
        # Change: 10,000 locations selected using probabilistic target-group sampling based on 
        # Fitzpatrick et al. (2013) MAXENT vs. MAXLIKE: Empirical Comparisons with 
        # Ant Species Distributions. Ecosphere
        set.seed(1)
        # Try using spatially thinned occurrences
        temp_occ_i <- st_as_sf(temp_occ, coords = c("x", "y"), na.fail = FALSE, crs = crs(IUCN_all))
        # Project occurrences from equal area projection to WGS84 
        temp_occ_i <- st_transform(temp_occ_i, crs(feature_all))
        
        # # Try using occurrences before spatial thinning
        # temp_occ_i <- temp_occ_in_1.1_buffer[temp_occ_feature_i[, 1] %in% na.omit(temp_occ_feature_i)[, 1], ]
        # temp_occ_i <- st_transform(temp_occ_i, crs = crs(feature_all))
        mask_i <- raster(paste0("D:/xinchen/mee_sdm_experiment/5_R_script/Output/calibration_raster/", spplist[i], "/bio1.asc"))
        mask_i <- mask_i > -10000 # Make sure the cell values are all converted to 1
        temp_occ_i <- sf::as_Spatial(temp_occ_i)
        
        # To deal with spatial sampling bias, we will create
        # a raster that reflects the sampling density
        # using kernel density estimation (KDE)
        bias_i <- cellFromXY(mask_i, temp_occ_i)
        cells_i <- unique(sort(bias_i))
        kernelXY_i <- xyFromCell(mask_i, cells_i)
        samps_i <- as.numeric(table(bias_i)) 
        # number of samples in each grid cell
        head(samps_i)
        max(samps_i) #6 is the maximum because we do not have many observations within 1km pixel
        
        # code to make KDE raster
        KDEsur <- sm.density(kernelXY_i, 
                             #weights=samps, 
                             display="none", 
                             ngrid=8726, 
                             ylim=c(-0.3068711,51.66813), 
                             xlim=c(67.20751,139.9242), 
                             nbins=0)
        
        KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
        KDErast <-  SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                                           length(KDEsur$estimate))))
        KDErast <- raster(KDErast)
        plot(KDErast)
        KDErast <- resample(KDErast, mask_i)
        KDErast <- KDErast*mask_i
        plot(KDErast)
        #KDEpts <- rasterToPoints(KDErast)
        
        temp_bg <- randomPoints(KDErast, n = 10000, ext = extent(mask_i), prob = TRUE) %>% as.data.frame()
        
        coordinates(temp_bg) <- ~x + y
        
        # Define CRS of bg points
        proj4string(temp_bg) <- CRS("+init=EPSG:4326")
        
        # Project bg points to equal area projection the same as IUCN range map
        temp_bg <- spTransform(temp_bg, crs(IUCN_all))
        
        # temp_bg <-terra::spatSample(x = terra::vect(IUCN_temp_occ_1.1_buffer), size = 10000, method = "random") # return SpatVector object
        # 
        # # Convert SpatVector object to sf to sp for spatial thinning using gridSample function
        # # Reference: https://rspatial.org/raster/sdm/2_sdm_occdata.html#data-cleaning
        # temp_bg <- sf::as_Spatial(sf::st_as_sf(temp_bg))
        
        # Spatial thinning for background points
        # Create a RasterLayer with the extent of occ
        r_10km_i <- raster(temp_bg)
        
        # set the resolution of the cells to 1km
        res(r_10km_i) <- 1000
        
        # expand (extend) the extent of the RasterLayer a little
        r_10km_i <- extend(r_10km_i, extent(r_10km_i) + 1000)
        
        # sample
        temp_bg <- gridSample(temp_bg, r_10km_i, n = 1) %>% as.data.frame()
        names(temp_bg) <- c("x", "y")
        
        # Export occ as csv files for background points sampling
        write.csv(temp_bg, paste0(bg_output, "/", spplist[i], "_background", ".csv", sep = ""), row.names = FALSE)
        
        # # Generate the mcp
        # if (nrow(temp_occ_in_1.1_buffer) > 5) {
        #   # Create concave hulls
        #   IUCN_temp_occ_mcp <- concaveman(temp_occ_in_1.1_buffer, concavity = 0.001)
        #   
        #   # Estimate mcp using alpha shape
        #   # Calculate mcp_ratio
        #   summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i], "area_mcp"] <- st_area(IUCN_temp_occ_mcp)
        #   summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i], "ratio_mcp"] <- st_area(IUCN_temp_occ_mcp)/st_area(IUCN_temp_occ_1.1_buffer)
        #   
        #   # Plot occurrence on range map using ggplot
        #   # Save the plot object in the list because of the memory issue
        #   plot_list[[i]] <- local({
        #     i <- i
        #     plot_i <-   ggplot() + geom_sf(data = IUCN_temp_occ_1.1_buffer, fill = "transparent") + 
        #       geom_sf(data = temp_occ_in_1.1_buffer, colour = "blue") +
        #       geom_sf(data = IUCN_temp_occ_mcp, fill = "transparent") +
        #       ggtitle (paste(spplist[i], "\nNumber of occurrence = ", nrow(temp_occ), 
        #                      "\nRatio of mcp = ", summary_mcp_ratio[summary_mcp_ratio$GBIF_name == spplist[i],
        #                                                             "ratio_mcp"])) + 
        #       theme_bw() + theme(plot.title = element_text(hjust = .5))
        #     print(plot_i)
        #   })
        # } 
      }
    }
  }
}
part1_end <- Sys.time()
