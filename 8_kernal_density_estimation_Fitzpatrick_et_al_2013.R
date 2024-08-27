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
        max(samps_i) #3 is the maximum because we do not have many observations within 1km pixel
        
        # code to make KDE raster
        KDEsur <- sm.density(kernelXY_i, 
                             #weights=samps, 
                             display="none", 
                             ngrid=963, 
                             ylim=c(14.2598,21.90146), 
                             xlim=c(-94.42582,-86.40082), 
                             nbins=0)
        
        KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
        KDErast <-  SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                                           length(KDEsur$estimate))))
        KDErast <- raster(KDErast)
        plot(KDErast)
        KDErast <- resample(KDErast, mask_i)
        KDErast <- KDErast*mask_i
        plot(KDErast)
        KDEpts <- rasterToPoints(KDErast)
        
        temp_bg <- randomPoints(KDErast, n = 10000, ext = extent(mask_i), prob = TRUE) %>% as.data.frame()
        
        coordinates(temp_bg) <- ~x + y
        
        # Define CRS of bg points
        proj4string(temp_bg) <- CRS("+init=EPSG:4326")
        
        # Project bg points to equal area projection the same as IUCN range map
        temp_bg <- spTransform(temp_bg, crs(IUCN_all))
