################################################################################
# script to model ant distributions in New England using maxent
# 
#
# Results described in:
# Fitzpatrick et al. (2013) MAXENT vs. MAXLIKE: Empirical Comparisons with 
# Ant Species Distributions. Ecosphere
#
# DESCRIPTION & NOTES
# The code needs the following to run:
# (1) climate & elevation rasters (provided as neClim)
# (2) ant distribution data (provided as antData_use4R_snappedToClim.csv). Note
# that as the file name suggests these data have been snapped to cell centriods.
#
# Code is provided as is, without support 
################################################################################

################################################################################
# CHUNK 1: Read in & prep data
################################################################################
library(dismo)
library(raster)
library(foreach)
library(doParallel)
library(dismo)
library(raster)
library(sp)
library(sm)

setwd("D:/xinchen/mee_sdm_experiment/8_taget_group_backgrounds/Week-7-20240813T221121Z-001/Week-7/data") # setwd to where data files are stored

# read in climate & elev data
antsGeoXY <- read.csv("antData.csv")
neClim <- stack("neClim.grd") # three variables
antsSp <- antsGeoXY
coordinates(antsSp) <- c("x", "y")
projection(antsSp) <- projection(neClim)

plot(neClim$bio_1)
points(antsSp, pch=20, cex=0.5, col=rgb(0,0,0,0.25))

# make a mask raster to use for processing
mask <- neClim[[1]]>-1000

# To deal with spatial sampling bias, we will create
# a raster that reflects the sampling density
# using kernel density estimation (KDE)
bias <- cellFromXY(mask, antsSp)
cells <- unique(sort(bias))
kernelXY <- xyFromCell(mask, cells)
samps <- as.numeric(table(bias)) 
# number of samples in each grid cell
head(samps)
max(samps) #255!

# code to make KDE raster
KDEsur <- sm.density(kernelXY, 
                     #weights=samps, 
                     display="none", 
                     ngrid=812, 
                     ylim=c(40,48), 
                     xlim=c(-75,-65), 
                     nbins=0)
KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
KDErast <-  SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                                   length(KDEsur$estimate))))
KDErast <- raster(KDErast)
plot(KDErast)
KDErast <- resample(KDErast, mask)
KDErast <- KDErast*mask
plot(KDErast)
KDEpts <- rasterToPoints(KDErast)
################################################################################
