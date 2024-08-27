## This script was used to model all species using Maxent and investigate response curve
## The occurrence data was downloaded from GBIF and then only the occurrences within 1.1 radius of IUCN range maps between 1990-2020 were selected
## The background points were generated within within 1.1 radius of IUCN range maps and then thinned to make sure one bg point per pixel
## The bioclimatic variables were calculated for 1991-2020 period

# part2
# 5. Extracting environmental factors for both occurrence and background points
# 6. Running Maxent for each species
# 7. Extracting information using rmaxent package

# Install 'rJava' package without compilation
# install.packages("rJava")

# Give Java Virtual machine 45GB RAM
# Give Java Virtual machine 45GB RAM
options(java.parameters = "-Xmx45g")

# Load packages
library(data.table)
library(sp)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(raster)

library(ggpubr)
library(glue)

# For R 4.4.0, installing 'dismo' automatically brings 'raster' package back
library(dismo)
library(ggplot2)

# Install rmaxent for extracting information from maxent object
library(devtools)
# devtools::install_github('johnbaums/rmaxent')

library(rmaxent)

# # Install kuenm package from Github
# devtools::install_github("marlonecobos/kuenm")

# Date: 05/14/2024
# This will need the install of Rtools version 4.4.6104
library(kuenm)

#######################################################################################################################################
# Part5. Extracting environmental factors for both occurrence and background points
#######################################################################################################################################

# Load environmental data first for getting crs and removing occurrences with NA in those layers
# Load data
feature_list <- list.files(path = "processed_raster_layers/", pattern = '\\.tif$', full.names = TRUE)

# Regular expression: https://stackoverflow.com/questions/4876813/using-r-to-list-all-files-with-a-specified-extension

# Combine all predictors
feature_all <- stack(feature_list)

# list all the csv files with occurrences
occlist <- list.files("Output/occurrence/", 
                      pattern = ".csv", full.names = TRUE)
bg_list <- list.files("Output/bg", 
                      pattern = ".csv", full.names = TRUE)

crs_list <- readRDS("crs_list.rds")

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

#######################################################################################################################################
# Part5. Modeling species distribution using Maxent
#######################################################################################################################################
# Maxent arguments
# Source: https://groups.google.com/g/maxent/c/yRBlvZ1_9rQ?pli=1
# https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix3_maxentParameters_v2.pdf

# Source: https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix2_prepPara.R

# For more argument: https://groups.google.com/g/maxent/c/yRBlvZ1_9rQ

prepPara <- function(userfeatures=NULL, #NULL=autofeature, could be any combination of # c("L", "Q", "H", "H", "P")
                     responsecurves=TRUE,
                     jackknife=TRUE,      
                     outputformat="logistic",
                     outputfiletype="asc", 
                     projectionlayers=NULL,
                     randomseed=FALSE,
                     removeduplicates=TRUE,
                     betamultiplier=NULL,
                     biasfile=NULL,
                     testsamplesfile=NULL,
                     replicates=1,
                     replicatetype="crossvalidate",
                     writeplotdata=TRUE,
                     extrapolate=TRUE,
                     doclamp=TRUE,
                     beta_threshold=NULL,
                     beta_categorical=NULL,
                     beta_lqp=NULL,
                     beta_hinge=NULL,
                     applythresholdrule=NULL
){
  #20 & 29-33 features, default is autofeature
  if(is.null(userfeatures)){
    args_out <- c("autofeature")
  } else {
    args_out <- c("noautofeature")
    if(grepl("L",userfeatures)) args_out <- c(args_out,"linear") else args_out <- c(args_out,"nolinear")
    if(grepl("Q",userfeatures)) args_out <- c(args_out,"quadratic") else args_out <- c(args_out,"noquadratic")
    if(grepl("H",userfeatures)) args_out <- c(args_out,"hinge") else args_out <- c(args_out,"nohinge")
    if(grepl("P",userfeatures)) args_out <- c(args_out,"product") else args_out <- c(args_out,"noproduct")
    if(grepl("T",userfeatures)) args_out <- c(args_out,"threshold") else args_out <- c(args_out,"nothreshold")
  }
  
  #1 
  if(responsecurves) args_out <- c(args_out,"responsecurves") else args_out <- c(args_out,"noresponsecurves")
  #2
  #if(picture) args_out <- c(args_out,"pictures") else args_out <- c(args_out,"nopictures")
  #3
  if(jackknife) args_out <- c(args_out,"jackknife") else args_out <- c(args_out,"nojackknife")
  #4
  args_out <- c(args_out,paste0("outputformat=",outputformat))
  #5
  args_out <- c(args_out,paste0("outputfiletype=",outputfiletype))
  #7
  if(!is.null(projectionlayers))    args_out <- c(args_out,paste0("projectionlayers=",projectionlayers))
  #10
  if(randomseed) args_out <- c(args_out,"randomseed") else args_out <- c(args_out,"norandomseed")
  #16
  if(removeduplicates) args_out <- c(args_out,"removeduplicates") else args_out <- c(args_out,"noremoveduplicates")
  #20 & 53-56
  # check if negative
  betas <- c( betamultiplier,beta_threshold,beta_categorical,beta_lqp,beta_hinge)
  if(! is.null(betas) ){
    for(i in 1:length(betas)){
      if(betas[i] <0) stop("betamultiplier has to be positive")
    }
  }
  if (  !is.null(betamultiplier)  ){
    args_out <- c(args_out,paste0("betamultiplier=",betamultiplier))
  } else {
    if(!is.null(beta_threshold)) args_out <- c(args_out,paste0("beta_threshold=",beta_threshold))
    if(!is.null(beta_categorical)) args_out <- c(args_out,paste0("beta_categorical=",beta_categorical))
    if(!is.null(beta_lqp)) args_out <- c(args_out,paste0("beta_lqp=",beta_lqp))
    if(!is.null(beta_hinge)) args_out <- c(args_out,paste0("beta_hinge=",beta_hinge))
  }
  #22
  if(!is.null(biasfile))    args_out <- c(args_out,paste0("biasfile=",biasfile))
  #23
  if(!is.null(testsamplesfile))    args_out <- c(args_out,paste0("testsamplesfile=",testsamplesfile))
  #24&25
  replicates <- as.integer(replicates)
  if(replicates>1 ){
    args_out <- c(args_out,
                  paste0("replicates=",replicates),
                  paste0("replicatetype=",replicatetype) )
  }
  #37
  if(writeplotdata) args_out <- c(args_out,"writeplotdata") else args_out <- c(args_out,"nowriteplotdata")
  #39
  if(extrapolate) args_out <- c(args_out,"extrapolate") else args_out <- c(args_out,"noextrapolate")
  #42
  if(doclamp) args_out <- c(args_out,"doclamp") else args_out <- c(args_out,"nodoclamp")
  #60
  if(!is.null(applythresholdrule))    args_out <- c(args_out,paste0("applythresholdrule=",applythresholdrule))
  
  return(args_out)
}

# Create a name list without character
spplist <- c()
for (i in 1:length(occlist)) {
  message(i)
  # Get the file name after the folder location
  spplist_i <- sub(".*Output/occurrence/", "", occlist[i])
  
  # Remove .csv after species name
  spplist[i] <- sub("*.csv$", "", spplist_i)
}

# Train Maxent
# Create a list to save Maxent object for extracting coefficients
maxent_list_1st <- vector("list", length(spplist))
names(maxent_list_1st) <- spplist

set.seed(3)
maxent_start <- Sys.time()
for (i in 1:length(spplist)) {
  
  me_i <- maxent(x = p_feature[[i]],
                 p = pa_list[[i]],
                 path = paste0("Output/Maxent/", spplist[i]),
                 args = prepPara(userfeatures = "LQ",
                                 doclamp = FALSE)
                 # For more argument: https://groups.google.com/g/maxent/c/yRBlvZ1_9rQ
  )
  
  # Save ith maxent object for making figures or extracting coefficients
  maxent_list_1st[[i]] <- me_i
  
}
maxent_end <- Sys.time()

# Save all raster layers in case the temporary file fails!!!!!!!!!!!!!!!!!!!!!
save.image("4_modeling_Maxent.RData")

#######################################################################################################################################
# Part6. Extracting lambda coefficients and response curves using rmaxent
#######################################################################################################################################

# When the lambda of linear feature is positive while the lambda of quadratic feature is negative
# The response curve looks like a inverted U-shape with one peak 

# Take the 1st species in the list as an example
# Showing all lambdas of the maxent object
parse_lambdas(maxent_list_1st[[1]])

# Extracting the lambda of cropland from the data frame returned by pasrse_lambda function
parse_lambdas(maxent_list_1st[[1]])$lambdas[6, c(2,3)]

# Extract response curve to a data frame and plot it
a <- as.data.frame(response(maxent_list_1st[[1]], "cropland"))

ggplot(a, aes(V1, p)) + geom_line() + ylim(c(0,1))

# Compare to the original response function in dismo
response(maxent_list_1st[[1]], "cropland")


# Create response curves for each variable
# Human impact index
response_hii <- data.frame()

for (i in 1:length(spplist)) {
  hii_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_human_impact_index_2010_only.dat", sep = ""))
  hii_i$species <- spplist[i]
  response_hii <- rbind(response_hii, hii_i)
}

# Plot the response curve using spp as a group
response_hii$species <- as.factor(response_hii$species)
ggplot(response_hii, aes(x, y, color = species)) + 
  geom_line(cex = 1.1) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

ggplot(response_hii, aes(x, y, alpha = species)) + 
  geom_line(cex = 1) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

# Road density
response_highway <- data.frame()

for (i in 1:length(spplist)) {
  highway_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_density_highways_only.dat", sep = ""))
  highway_i$species <- spplist[i]
  response_highway <- rbind(response_highway, highway_i)
}

response_highway$species <- as.factor(response_highway$species)
ggplot(response_highway, aes(x, y, alpha = species)) + 
  geom_line(cex = 1) +
  scale_x_continuous(limits = c(0, 5000), "Highway Density (m/sq km)") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")


# Create a data frame for saving permutation variable importance
var_names <- c("Nighttime Lights",
               "Mean Temperature of Warmest Quarter", 
               "Mean Temperature of Coldest Quarter", 
               "Annual Precipitation", 
               "Precipitation of Wettest Quarter", 
               "Precipitation of Driest Quarter", 
               "Annual Mean Temperature", 
               "Cropland", 
               "Highway Density", 
               #"Food Waste Index",
               "Human Impact", 
               "Pasture", 
               "Population")

# Create a data frame that has the same number row as spplist
perm_var_imp <- data.frame(variable = rep(var_names, 351),
                           value = rep(NA, 12*351),
                           species = rep(NA, 12*351))

for (i in 1:length(spplist)) {
  
  perm_var_imp$species[(12*(i-1) + 1):(12*i)] <- spplist[i]
  perm_var_imp$value[(12*(i-1) + 1):(12*i)] <- maxent_list_1st[[i]]@results[19:30, 1]
  
}

perm_var_imp$variable <- as.factor(perm_var_imp$variable)

library(ggplot2)
ggplot(perm_var_imp, aes(x = variable, y = value, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot() + ylab("Permutation variable importance (%)") + scale_x_discrete(guide = guide_axis(angle = 45), "") +
  theme(legend.position = "none")



# Sorting variable importance and add the rank of variables
perm_var_imp_rank <- perm_var_imp
perm_var_imp_rank$rank <- 1

for (i in 1:length(spplist)) {
  perm_var_imp_rank_i <- perm_var_imp_rank[perm_var_imp_rank$species == spplist[i], ]
  
  # Order variable importance
  perm_var_imp_rank_i <- perm_var_imp_rank_i[order(-perm_var_imp_rank_i$value), ]
  perm_var_imp_rank_i$rank <- 1:12
  
  # Replace the original dataframe
  perm_var_imp_rank[perm_var_imp_rank$species == spplist[i], ] <- perm_var_imp_rank_i
  
}

library(ggplot2)
ggplot(perm_var_imp_rank, aes(x = variable, y = rank, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot() + ylab("Rank of permutation variable importance") + scale_x_discrete(guide = guide_axis(angle = 0), "") +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  coord_flip() +
  theme(legend.position = "none")

# Compare the number of species between annual mean temperature and human impact across all ranks of these two variables
ggplot(perm_var_imp_rank[perm_var_imp_rank$variable == "Human Impact"|perm_var_imp_rank$variable == "Annual Mean Temperature", ], aes(x = rank, fill = variable)) + 
  geom_histogram(position = "dodge") + scale_x_continuous(name = "Rank of permutation variable importance", breaks = c(1:12)) +
  theme(legend.position = "bottom")

# Compare the number of species between Mean Temperature of Coldest Quarter and human impact across all ranks of these two variables
ggplot(perm_var_imp_rank[perm_var_imp_rank$variable == "Human Impact"|perm_var_imp_rank$variable == "Mean Temperature of Coldest Quarter", ], aes(x = rank, fill = variable)) + 
  geom_histogram(position = "dodge") + scale_x_continuous(name = "Rank of permutation variable importance", breaks = c(1:12)) +
  theme(legend.position = "bottom")

# Look at how Human Impact versus Mean Temperature of Coldest Quarter at each value of variable importance rank
# This tells us which variable is more important when they are at higher ranks and which variable is more important when they are at lower ranks
ggplot(perm_var_imp_rank[perm_var_imp_rank$variable == "Human Impact"|perm_var_imp_rank$variable == "Annual Mean Temperature", ], aes(x = rank, fill = variable)) + 
  geom_histogram(position = "dodge", stat = "count") + xlab("Rank of permutation variable importance") + 
  scale_x_continuous(guide = guide_axis(angle = 0), name = "", breaks = c(1:12)) + 
  theme(legend.position = "bottom")

ggplot(perm_var_imp_rank[perm_var_imp_rank$variable == "Human Impact"|perm_var_imp_rank$variable == "Mean Temperature of Coldest Quarter", ], aes(x = rank, fill = variable)) + 
  geom_histogram(position = "dodge", stat = "count") + xlab("Rank of permutation variable importance") + 
  scale_x_continuous(guide = guide_axis(angle = 0), name = "", breaks = c(1:12)) + 
  theme(legend.position = "bottom")

###########################################################################################################################################
# Date: 03/24/2023
# Change: Add index for constraining the values to the limit of that variable found in the training range (clampping)
#         Split the values of predictor based on the quantile of the variable
# Example: Add [!duplicated(response_hii$y)&response_hii$x >= 0, ] in the response of human index
# Note: Some species did not have response$x >= 0 indicating that human impact was not important
###########################################################################################################################################

# For positive response, which.max >= 4000
# For negative response, which.max <10
# For positive and negative response, which.max >=10 but < 4000
# Create a new object
response_hii$category <- "not important"

for (i in 1:length(spplist)) {
  # Remove duplicated values and those less than 0 for ith species
  response_hii_i <- response_hii[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$species == spplist[i], ]
  # Check the maximum probabilities for ith species
  if (nrow(response_hii_i) != 0) {
    
    max_position <- which.max(response_hii_i$y)
    
    if (response_hii_i$x[max_position] < 100) {
      
      response_hii$category[response_hii$species == spplist[i]] <- "decreasing"
      
    } else if (response_hii_i$x[max_position] >= 100&response_hii_i$x[max_position] < 4000) {
      
      response_hii$category[response_hii$species == spplist[i]] <- "concave"
      
    } else {
      
      response_hii$category[response_hii$species == spplist[i]] <- "increasing"
      
    }
  }
  
}

# Date: 02/23/2024
# Update: Try to create a continuous variable in the data frame of response curve to visualize different species using a continuous color palette
# Solution: Coercing species name to a numerical value
# response_hii$color <- as.numeric(response_hii$species)

# Plot the response curves for each category
ggplot(response_hii[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "increasing", ], aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  # scale_color_paletteer_c("grDevices::Plasma",
  #                         name = "Human Impact Index",
  #                         direction = 1, 
  #                         limits=c(0, 8), 
  #                         breaks = c(seq(0, 8, by = 1))
  #                         ) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")
# Calculate how many species have increasing trend
length(unique(response_hii$species[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "increasing"]))
# 141

ggplot(response_hii[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "concave", ], aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")
# Calculate how many species have concave trend
length(unique(response_hii$species[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "concave"]))
# 151

ggplot(response_hii[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "decreasing", ], aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")
# Calculate how many species have decreasing trend
length(unique(response_hii$species[!duplicated(response_hii$y)&response_hii$x >= 0&response_hii$category == "decreasing"]))
# 31

#######################################################################################################################################################
# Date: 04/08/2024
# Change: Match species name with their GBIF order to visualize the response curve for different orders
# Logic: The response curve has GBIF name instead of IUCN name but the order is connected with IUCN name
#        Therefore, using GBIF name and IUCN name matching table can retrieve order for each species.

# Load GBIF data with year, genericName, specificEpithet, order, and basisOfRecord only
Sys.time()
occurrence <- fread("D:/Xin/Data_Back_Up/FSU/D_drive/Projects/SDM_human/Raw_data/IUCN_mammals/IUCN_intersection_ratio/IUCN_intersection_calculation/occurrences_GBIF_all_mammals/0065184-230530130749713/occurrence.txt",
                    select = c("year", "genericName", "specificEpithet", "order", "basisOfRecord"))
Sys.time()

# Split occurrence for 1990 to 2020
occurrence <- occurrence[occurrence$year >= 1990& occurrence$year <= 2020, ]

# Create a column using genericName+specificEpithet
occurrence$gbif_name <- paste0(occurrence$genericName, paste0(" ", occurrence$specificEpithet))

# Keep basisOfRecord of "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", and "OCCURRENCE"
# Source: https://docs.gbif.org/course-data-use/en/basis-of-record.html
occurrence <- occurrence[occurrence$basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE", 
                                                         "MATERIAL_CITATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE"), ]

# Remove duplicated records
occurrence <- occurrence[!duplicated(occurrence$gbif_name), ]

# Create a column for matching GBIF species names with the modeled spplist
occurrence$GBIF_name <- gsub(" ", "_", as.vector(occurrence$gbif_name))


# Merge order to the matching table
species_order <- data.frame(GBIF_name = spplist)
species_order <- merge(species_order, occurrence[, c("order", "GBIF_name")], by = "GBIF_name")

############################################################################################################################################################
# Not run
# # Load Synonymy_table_valid_species_only.csv from PHYLACINE data set that has IUCN name and Order for each species
# IUCN_order <- read.csv("Synonymy_table_valid_species_only.csv")

# Date: 04/08/2024
# Note: Let's not use order from PHYLACINE now
# Not run
# Note: species_searchedname in IUCN_species_list_9_matching should be the same as Binomial.1.2 in IUCN_order
# Test: Bionomial.1.2 in IUCN_order is equal to IUCN.2016.3.Genus + "_" + IUCN.2016.3.Species
# Not run

# # Create a dataframe to connect GBIF_name with order based on IUCN name
# species_order <- merge(IUCN_species_list_9_matching, IUCN_order[, 1:2], 
#                        by.x = "species_searchedname", by.y = "Binomial.1.2", all.x = TRUE)
# 
# # Check NA's in Order.1.2 column
# species_order[is.na(species_order$Order.1.2), ]

# There are some NA's in species order because Binomial.1.2 may use GBIF name
# For example,  Myodes_californicus (GBIF name) is the Synonym of Clethrionomys californicus (IUCN name) but PHYLACINE data set uses its GBIF name instead of
# its IUCN name. Source: https://www.gbif.org/species/2439141

# Export species_order as csv and manually find the order for the species from Synonymy_table_valid_species_only.csv by searching for GBIF_name
# because the IUCN/species_searched name does not match with Binomial.1.2
# write.csv(species_order, "species_order.csv", row.names = FALSE)

# After manually checking the missing order in Synonymy_table_valid_species_only.csv and GBIF
# Load the fixed data back to the environment
# species_order_checked <- read.csv("species_order_manual_checked.csv")
############################################################################################################################################################


# Merge response of human impact with species order
response_hii_order <- merge(response_hii, species_order, by.x = "species", by.y = "GBIF_name", all.x = TRUE)

# For Chiroptera as a quick example
ggplot(response_hii_order[!duplicated(response_hii_order$y)&response_hii_order$x >= 0&response_hii_order$x >= 0&response_hii_order$order == "Chiroptera", ], 
       aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  #paletteer::scale_color_paletteer_d("LaCroixColoR::Apricot")
  paletteer::scale_color_paletteer_d("ggthemes::excel_Ion",
                                     name = "Species",
                                     direction = 1 
                                     #limits=c(0, 8), 
                                     #breaks = c(seq(0, 8, by = 1))
  ) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "bottom")

#######################################################################################################################################################
# Date: 03/02/2024
# Change: Match species name with their order to visualize the variable importance for different orders
# Logic: The species names have been successfully matched with their orders in the previous step so merging perm_var_imp and per_var_imp_rank with 
#        species_order_checked to add order to variable importance for each species

# Merge variable importance and rank of variable importance with species order
perm_var_imp_order <- merge(perm_var_imp, species_order, by.x = "species", by.y = "GBIF_name", all.x = TRUE)

perm_var_imp_rank_order <- merge(perm_var_imp_rank, species_order, by.x = "species", by.y = "GBIF_name", all.x = TRUE)

library(ggplot2)
list_order <- unique(perm_var_imp_order$order)

# Loop over all orders to generate variable importance
for (i in 1:length(list_order)) {
  # Subset the var importance using order i
  perm_var_imp_order_i <- perm_var_imp_order[perm_var_imp_order$order == list_order[i], ]
  
  print(ggplot(perm_var_imp_order_i, aes(x = variable, y = value, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
          geom_boxplot() + ylab("Permutation variable importance (%)") + scale_x_discrete(guide = guide_axis(angle = 45), "") +
          ggtitle(paste("Variable importance for ", list_order[i], sep = "")) +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5)))
}

# Loop over all orders to generate order of variable importance
for (i in 1:length(list_order)) {
  # Subset the var importance using order i
  perm_var_imp_rank_order_i <- perm_var_imp_rank_order[perm_var_imp_rank_order$order == list_order[i], ]
  
  print(ggplot(perm_var_imp_rank_order_i, aes(x = variable, y = rank, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
          geom_boxplot() + ylab("Permutation variable importance (%)") + scale_x_discrete(guide = guide_axis(angle = 0), "") +
          scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
          ggtitle(paste("Variable rank for ", list_order[i], sep = "")) +
          coord_flip() +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5)))
}

ggplot(perm_var_imp, aes(x = variable, y = value, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot() + ylab("Permutation variable importance (%)") + scale_x_discrete(guide = guide_axis(angle = 45), "") +
  theme(legend.position = "none")

ggplot(perm_var_imp_rank, aes(x = variable, y = rank, fill = variable)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot() + ylab("Rank of permutation variable importance") + scale_x_discrete(guide = guide_axis(angle = 0), "") +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  coord_flip() +
  theme(legend.position = "none")

#######################################################################################################################################################
# Date: 03/08/2024
# Note: Check the correlation between human impact and nighttime lights for two bat species Antrozous pallidus and Eptesicus fuscus
# Antrozous pallidus
cor(na.omit(p_feature[[6]]))

# Eptesicus fuscus
cor(na.omit(p_feature[[73]]))

#######################################################################################################################################################
# Date: 03/15/2024
# Note: Split response curve of human impact using mid point where the probability of presence is at the maximum
#       The range of human impact is [0, 5865]

# Remove all the values beyond the range of human impact because they are extrapolation
response_hii_order_sub <- response_hii_order[response_hii_order$x >= 0 & response_hii_order$x <= 5865, ]

# For positive response, max position >= median of range of human impact
# For negative response, max position < median of range of human impact
# Create a new object
response_hii_order_sub$category <- "loser"

for (i in 1:length(spplist)) {
  # Remove duplicated values and those less than 0 for ith species
  response_hii_order_sub_i <- response_hii_order_sub[response_hii_order_sub$species == spplist[i], ]
  
  # Check the maximum probabilities for ith species
  if (nrow(response_hii_order_sub_i) != 0) {
    
    max_position <- which.max(response_hii_order_sub_i$y)
    #median_position <- floor(median(response_hii_order_sub_i$x))
    median_position <- floor(quantile(response_hii_order_sub_i$x, probs = 0.75))
    
    if (response_hii_order_sub_i$x[max_position] < median_position) {
      
      response_hii_order_sub$category[response_hii_order_sub$species == spplist[i]] <- "loser"
      
    } else if (response_hii_order_sub_i$x[max_position] >= median_position) {
      
      response_hii_order_sub$category[response_hii_order_sub$species == spplist[i]] <- "winner"
      
      # } else {
      #   
      #   response_hii_order_sub$category[response_hii_order_sub$species == spplist[i]] <- "loser"
      
    }
  }
}

# Visualize it
ggplot(response_hii_order_sub[response_hii_order_sub$order == "Rodentia"&response_hii_order_sub$category == "winner", ], 
       aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

ggplot(response_hii_order_sub[response_hii_order_sub$order == "Rodentia"&response_hii_order_sub$category == "loser", ], 
       aes(x, y, color = species)) + 
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
# Date: 03/21/2024
# Calculate the slope of each pair of points 
# Create a function to calculate the slope of every two points
# slope_function <- function(data){
#   data_save <- data.frame()
#   # Loop through all species j
#   for (j in 1:length(unique(data$species))) {
#     data_j <- data[data$species == unique(data$species)[j], ]
#     
#     # Loop through every row of data_j
#     for (i in 2:nrow(data_j)) {
#       
#       # Extract the ith and i-1 row from the data
#       data_i <- data_j[(i-1):i, ] # Date: 03/28/2023 Change i:(i-1) to (i-1):i because the logic of diff function
#       slope_i <- diff(data_i[, "y"])/diff(data_i[, "x"])
#       data_j[i, "slope"] <- slope_i
#     }
#     # Update data_j to the blank data frame
#     data_save <- rbind(data_save, data_j)
#   }
#   return(data_save)
# }
# 
# # Calculate slope of every pair of points for each species
# response_slope <- slope_function(response_hii_order_sub)
# 
# # Use sum of slope to define the winner and loser
# sum_slope <- response_slope[!is.na(response_slope$slope), ]
# sum_slope <- aggregate(slope ~ species, data = sum_slope, sum)
# 
# hist(sum_slope$slope)

# Thresholding the response of human impact for each species using their corresponding training data
response_hii_order_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_hii_i <- max(p_feature[[i]][9], na.rm = TRUE)
  lower_hii_i <- min(p_feature[[i]][9], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_hii_order_no_extra_i <- response_hii_order_sub[response_hii_order_sub$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_hii_order_no_extra_i <- response_hii_order_no_extra_i[response_hii_order_no_extra_i$x >= lower_hii_i, ]
  response_hii_order_no_extra_i <- response_hii_order_no_extra_i[response_hii_order_no_extra_i$x <= upper_hii_i, ]
  
  # Merge ith data frame back to the final data frame
  response_hii_order_no_extra <- rbind(response_hii_order_no_extra, response_hii_order_no_extra_i)
}

################################################################################################################################
# Calculate the slope of each pair of points 
# Create a function to calculate the slope of every two points
diff_function <- function(data){
  data_save <- data.frame()
  # Loop through all species j
  for (j in 1:length(unique(data$species))) {
    data_j <- data[data$species == unique(data$species)[j], ]
    
    # Loop through every row of data_j
    for (i in 2:nrow(data_j)) {
      
      # Extract the ith and i-1 row from the data
      data_i <- data_j[(i-1):i, ]
      if (nrow(data_i) == 1) {
        
        data_j[i, "slope"] <- NA
        
      } else {
        
        slope_i <- diff(data_i[, "y"])/1 # Convert diff(data_i[, "x"]) to unit 1 so that each species would be fair to be compared
        # The logic is that no matter how fast or slow the probability of presence change, we only care about increase or decrease
        data_j[i, "slope"] <- slope_i
        
      }
    }
    # Update data_j to the blank data frame
    data_save <- rbind(data_save, data_j)
  }
  return(data_save)
}
################################################################################################################################

# Calculate slope of every pair of points for each species
response_diff <- diff_function(response_hii_order_no_extra)

# Use sum of slope to define the winner and loser
sum_diff <- response_diff[!is.na(response_diff$slope), ]


################################################################################################################################


##################################################################################################################################################################
# Date: 05/30/2024
# Change: Win: >90% of slopes are +    &   sum(slope) > 0.1
#         Lose: >90% of slopes are -   &   sum(slope) <  -0.1


##################################################################################################################################################################

##############################################################################################################################################
# Date: 05/30/2024
# Change: Try visualizing response curve using presences
#         Thresholding the response curve based on species presence using p_list instead of p_feature

# Thresholding the response of human impact for each species using their corresponding training data
response_hii_order_presence_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_hii_i <- max(p_list[[i]][9], na.rm = TRUE)
  lower_hii_i <- min(p_list[[i]][9], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_hii_order_presence_no_extra_i <- response_hii_order_sub[response_hii_order_sub$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_hii_order_presence_no_extra_i <- response_hii_order_presence_no_extra_i[response_hii_order_presence_no_extra_i$x >= lower_hii_i, ]
  response_hii_order_presence_no_extra_i <- response_hii_order_presence_no_extra_i[response_hii_order_presence_no_extra_i$x <= upper_hii_i, ]
  
  # Merge ith data frame back to the final data frame
  response_hii_order_presence_no_extra <- rbind(response_hii_order_presence_no_extra, response_hii_order_presence_no_extra_i)
}
##############################################################################################################################################

# For response curves of human impact index based on presence
# Calculate slope of every pair of points for each species
response_diff <- diff_function(response_hii_order_no_extra)

# Use sum of slope to define the winner and loser
sum_diff <- response_diff[!is.na(response_diff$slope), ] # Note: This sum_diff is NOTTTTT aggregated one with only one summed slope per species!!!

#xf3 check which sp has all positive slopes
s="Ammospermophilus_harrisii"
group_full = as.character(unique(sum_diff$species))
group_win_abs = vector()
group_lose_abs = vector()
for(s in group_full){
  slope_total  = sum(subset(sum_diff,species==s)$slope)
  one_check = subset(sum_diff,species==s)$slope>0
  #if( all(one_check) )   group_win_abs = c(group_win_abs,s)    # all positive
  #if( all(!one_check) )   group_lose_abs = c(group_lose_abs,s) # all negative
  #if( sum(one_check)/length(one_check)>=0.9     )   group_win_abs = c(group_win_abs,s)
  #if( sum(!one_check)/length(one_check)>=0.9    )   group_lose_abs = c(group_lose_abs,s)
  if( sum(one_check)/length(one_check)>=0.9   & slope_total>0.1   )   group_win_abs = c(group_win_abs,s)
  if( sum(!one_check)/length(one_check)>=0.9  & slope_total< (-0.1)  )   group_lose_abs = c(group_lose_abs,s)
}
group_other = group_full[!group_full %in%  c(group_win_abs,group_lose_abs) ]
win_lose_df = data.frame( species=c(group_win_abs,
                                    group_lose_abs,
                                    group_other))
win_lose_df$win_group =c(rep("win",length(group_win_abs)),
                         rep("lose",length(group_lose_abs)),
                         rep("other",length(group_other)))
table(win_lose_df$win_group)
df_curve_abs1 = merge(response_hii_order_no_extra, win_lose_df, by="species")
ggplot(df_curve_abs1,
       aes(x, y, group=species,color = win_group)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Relative probability of presence") +
  ylim(c(0,1)) +facet_grid(order~win_group)
#theme(legend.position = "none")
ggsave("09_debug/response curve v1.tiff",dpi = 600,width = 20,height = 20,units = "cm")

# Try visualizing response curve by winner/Mid/Loser x Order with gradients of variable importance of human impact index variable index
perm_var_imp_order_human_impact <- perm_var_imp_order[perm_var_imp_order$variable == "Human Impact", ]
perm_var_imp_order_human_impact$variable <- "human_impact_index_2010"

df_curve_abs1_plot <- merge(df_curve_abs1, perm_var_imp_order_human_impact, 
                            by = c("species", "variable", "order"))

# Vertical 600x750
ggplot(df_curve_abs1_plot, aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")

# Date: 06/13/2024
# Change: Rename loser-other-winners with losers-neutrals-winners
# New facet label names for win_group
win_group_labs <- c("losers", "neutrals", "winners")
names(win_group_labs) <- c("lose", "other", "win")

# Vertical 800x900
ggplot(df_curve_abs1_plot, aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")

# For response curves of human impact index based on presence
# Calculate slope of every pair of points for each species
response_presence_diff <- diff_function(response_hii_order_presence_no_extra)

# Use sum of slope to define the winner and loser
sum_presence_diff <- response_presence_diff[!is.na(response_presence_diff$slope), ] # Note: This sum_presence_diff is NOTTTTT aggregated one with only one summed slope per species!!!

#xf3 check which sp has all positive slopes
s="Ammospermophilus_harrisii"
group_full = as.character(unique(sum_presence_diff$species))
group_win_abs = vector()
group_lose_abs = vector()
for(s in group_full){
  slope_total  = sum(subset(sum_presence_diff,species==s)$slope)
  one_check = subset(sum_presence_diff,species==s)$slope>0
  
  # # Add one condition that is the endpoint > first point, indicating that the response curve is overall increasing and probability of presence at hii > 0 is greater than that at hii = 0
  # first_hii <- dplyr::first(subset(sum_presence_diff,species==s)$slope)
  # last_hii <- dplyr::last(subset(sum_presence_diff,species==s)$slope)
  
  #if( all(one_check) )   group_win_abs = c(group_win_abs,s)    # all positive
  #if( all(!one_check) )   group_lose_abs = c(group_lose_abs,s) # all negative
  #if( sum(one_check)/length(one_check)>=0.9     )   group_win_abs = c(group_win_abs,s)
  #if( sum(!one_check)/length(one_check)>=0.9    )   group_lose_abs = c(group_lose_abs,s)
  if( sum(one_check)/length(one_check)>=0.9  & slope_total> 0.1)   group_win_abs = c(group_win_abs,s)
  if( sum(!one_check)/length(one_check)>=0.9 & slope_total< (-0.1))   group_lose_abs = c(group_lose_abs,s)
}

group_other = group_full[!group_full %in%  c(group_win_abs,group_lose_abs) ]
win_lose_df = data.frame( species=c(group_win_abs,
                                    group_lose_abs,
                                    group_other))
win_lose_df$win_group =c(rep("win",length(group_win_abs)),
                         rep("lose",length(group_lose_abs)),
                         rep("other",length(group_other)))
table(win_lose_df$win_group)
df_curve_abs2 <- merge(response_hii_order_presence_no_extra, win_lose_df,by="species")
ggplot(df_curve_abs2,
       aes(x, y, group=species,color = win_group)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Relative probability of presence") +
  ylim(c(0,1)) +facet_grid(order~win_group)
#theme(legend.position = "none")
ggsave("09_debug/response curve v1.tiff",dpi = 600,width = 20,height = 20,units = "cm")

# Try visualizing response curve by winner/Mid/Loser x Order with gradients of variable importance of human impact index variable index
perm_var_imp_order_human_impact <- perm_var_imp_order[perm_var_imp_order$variable == "Human Impact", ]
perm_var_imp_order_human_impact$variable <- "human_impact_index_2010"

df_curve_abs2_plot <- merge(df_curve_abs2, perm_var_imp_order_human_impact, 
                            by = c("species", "variable", "order"))

ggplot(df_curve_abs2_plot, aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free") +
  theme(legend.position = "none")

# Export defined winner, losers, and other for matching with trait data
write.csv(win_lose_df, "win_lose_df.csv", row.names = FALSE)

# Sum over the slope (difference) for the response curves
sum_diff_presence <- aggregate(slope ~ species, data = sum_presence_diff, sum)
write.csv(sum_diff_presence[, c("species", "slope")], "sum_diff_presence.csv", row.names = FALSE)

##################################################################################################################################################################































sum_diff <- aggregate(slope ~ species, data = sum_diff, sum)

hist(sum_diff$slope, main = "Histogram of slopes for human impact")
round(quantile(sum_diff$slope), digits = 2)

# 0%    25%    50%   75%   100% 
# -0.72 -0.16  0.26  0.51  0.90 

# Divide all species into different group
positive_diff <- sum_diff[sum_diff$slope > 0, "species"]
negative_diff <- sum_diff[sum_diff$slope <= 0, "species"]

# Visualize the species with positive trend and negative trend
ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% positive_diff, ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% negative_diff, ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

# Visualize the species with different quantiles of slopes
# upper 75% quantile
ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% sum_diff$species[sum_diff$slope >= 0.51], ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

# lower 25% quantile
ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% sum_diff$species[sum_diff$slope < -0.16], ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")


# Split winner from loser using 50% quantile of slopes
ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% sum_diff$species[sum_diff$slope >= 0.24], ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

ggplot(response_hii_order_no_extra[response_hii_order_no_extra$species %in% sum_diff$species[sum_diff$slope < 0.24], ],
       aes(x, y, color = species)) +
  geom_line(cex = 0.8, alpha = 0.5) +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  theme(legend.position = "none")

# Create winner versus loser species list
spp_winner_loser <- sum_diff
spp_winner_loser$group <- "Mid"
spp_winner_loser$group[sum_diff$slope >= 0.26] <- "Winner"
spp_winner_loser$group[sum_diff$slope < -0.16] <- "Loser"

# Merge winner versus loser back to permutation importance of human impact
perm_var_imp_order_spp_group <- merge(perm_var_imp_order, perm_var_imp_rank, by = c("variable", "value", "species"))

perm_var_imp_order_spp_group <- merge(perm_var_imp_order_spp_group, spp_winner_loser[, c(1,3)], by = "species")

# Export the file for investigating trait
write.csv(perm_var_imp_order_spp_group, "perm_var_imp_order_spp_group.csv", row.names = FALSE)

#############################################################################################################################
## Date: 04/21/2024
## Note: Calculate overall variable importance of human factors versus that of climatic factors
##       Calculate the sum of slope (difference) of the response curves of all human factors
#############################################################################################################################

# Create response curves for each human factor
# Human impact index
response_hii <- data.frame()

for (i in 1:length(spplist)) {
  hii_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_human_impact_index_2010_only.dat", sep = ""))
  hii_i$species <- spplist[i]
  response_hii <- rbind(response_hii, hii_i)
}

# Cropland
response_crop <- data.frame()

for (i in 1:length(spplist)) {
  crop_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_cropland_only.dat", sep = ""))
  crop_i$species <- spplist[i]
  response_crop  <- rbind(response_crop, crop_i)
}

# Density highways
response_highways <- data.frame()

for (i in 1:length(spplist)) {
  highways_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_density_highways_only.dat", sep = ""))
  highways_i$species <- spplist[i]
  response_highways  <- rbind(response_highways, highways_i)
}


# Nighttimelights
response_nighttimelights <- data.frame()

for (i in 1:length(spplist)) {
  nighttimelights_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_Nighttimelights_2012_only.dat", sep = ""))
  nighttimelights_i$species <- spplist[i]
  response_nighttimelights  <- rbind(response_nighttimelights, nighttimelights_i)
}

# Pasture
response_pasture <- data.frame()

for (i in 1:length(spplist)) {
  pasture_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_pasture_only.dat", sep = ""))
  pasture_i$species <- spplist[i]
  response_pasture  <- rbind(response_pasture, pasture_i)
}

# Population
response_population <- data.frame()

for (i in 1:length(spplist)) {
  population_i <- read.csv(paste("TestRun_NA_1990_2020/output/", spplist[i], "/plots/species_population_2010_only.dat", sep = ""))
  population_i$species <- spplist[i]
  response_population  <- rbind(response_population, population_i)
}


# Thresholding the response of human impact for each species using their corresponding training data
# Human impact index
response_hii_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_hii_i <- max(p_feature[[i]]["human_impact_index_2010"], na.rm = TRUE)
  lower_hii_i <- min(p_feature[[i]]["human_impact_index_2010"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_hii_no_extra_i <- response_hii[response_hii$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_hii_no_extra_i <- response_hii_no_extra_i[response_hii_no_extra_i$x >= lower_hii_i, ]
  response_hii_no_extra_i <- response_hii_no_extra_i[response_hii_no_extra_i$x <= upper_hii_i, ]
  
  # Merge ith data frame back to the final data frame
  response_hii_no_extra <- rbind(response_hii_no_extra, response_hii_no_extra_i)
}

# Cropland
response_crop_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_cropland_i <- max(p_feature[[i]]["cropland"], na.rm = TRUE)
  lower_cropland_i <- min(p_feature[[i]]["cropland"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_crop_no_extra_i <- response_crop[response_crop$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_crop_no_extra_i <- response_crop_no_extra_i[response_crop_no_extra_i$x >= lower_cropland_i, ]
  response_crop_no_extra_i <- response_crop_no_extra_i[response_crop_no_extra_i$x <= upper_cropland_i, ]
  
  # Merge ith data frame back to the final data frame
  response_crop_no_extra <- rbind(response_crop_no_extra, response_crop_no_extra_i)
}

# Density highways
response_highways_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_highways_i <- max(p_feature[[i]]["density_highways"], na.rm = TRUE)
  lower_highways_i <- min(p_feature[[i]]["density_highways"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_highways_no_extra_i <- response_highways[response_highways$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_highways_no_extra_i <- response_highways_no_extra_i[response_highways_no_extra_i$x >= lower_highways_i, ]
  response_highways_no_extra_i <- response_highways_no_extra_i[response_highways_no_extra_i$x <= upper_highways_i, ]
  
  # Merge ith data frame back to the final data frame
  response_highways_no_extra <- rbind(response_highways_no_extra, response_highways_no_extra_i)
}

# Nighttimelights
response_nighttimelights_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_nighttimelights_i <- max(p_feature[[i]]["Nighttimelights_2012"], na.rm = TRUE)
  lower_nighttimelights_i <- min(p_feature[[i]]["Nighttimelights_2012"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_nighttimelights_no_extra_i <- response_nighttimelights[response_nighttimelights$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_nighttimelights_no_extra_i <- response_nighttimelights_no_extra_i[response_nighttimelights_no_extra_i$x >= lower_nighttimelights_i, ]
  response_nighttimelights_no_extra_i <- response_nighttimelights_no_extra_i[response_nighttimelights_no_extra_i$x <= upper_nighttimelights_i, ]
  
  # Merge ith data frame back to the final data frame
  response_nighttimelights_no_extra <- rbind(response_nighttimelights_no_extra, response_nighttimelights_no_extra_i)
}

# Pasture
response_pasture_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_pasture_i <- max(p_feature[[i]]["pasture"], na.rm = TRUE)
  lower_pasture_i <- min(p_feature[[i]]["pasture"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_pasture_no_extra_i <- response_pasture[response_pasture$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_pasture_no_extra_i <- response_pasture_no_extra_i[response_pasture_no_extra_i$x >= lower_pasture_i, ]
  response_pasture_no_extra_i <- response_pasture_no_extra_i[response_pasture_no_extra_i$x <= upper_pasture_i, ]
  
  # Merge ith data frame back to the final data frame
  response_pasture_no_extra <- rbind(response_pasture_no_extra, response_pasture_no_extra_i)
}

# Population
response_population_no_extra <- data.frame()

for (i in 1:length(spplist)) {
  # Load species name
  spp_i <- spplist[i]
  # Calculate the bounds for human impact
  upper_population_i <- max(p_feature[[i]]["population_2010"], na.rm = TRUE)
  lower_population_i <- min(p_feature[[i]]["population_2010"], na.rm = TRUE)
  # Pass the variable importance to data frame i
  response_population_no_extra_i <- response_population[response_population$species == spp_i, ]
  # Thresholding ith data frame using the upper and lower bound
  response_population_no_extra_i <- response_population_no_extra_i[response_population_no_extra_i$x >= lower_population_i, ]
  response_population_no_extra_i <- response_population_no_extra_i[response_population_no_extra_i$x <= upper_population_i, ]
  
  # Merge ith data frame back to the final data frame
  response_population_no_extra <- rbind(response_population_no_extra, response_population_no_extra_i)
}

# Calculate the sum of slope (difference) for the response curves of all human factors
# Human impact index
response_diff_hii <- diff_function(response_hii_no_extra)
response_diff_hii <- response_diff_hii[!is.na(response_diff_hii$slope), ]
sum_diff_hii <- aggregate(slope ~ species, data = response_diff_hii, sum)
names(sum_diff_hii)[2] <- "slope_of_hii"

# Cropland
response_diff_crop <- diff_function(response_crop_no_extra)
response_diff_crop <- response_diff_crop[!is.na(response_diff_crop$slope), ]
sum_diff_crop <- aggregate(slope ~ species, data = response_diff_crop, sum)
names(sum_diff_crop)[2] <- "slope_of_crop"

# Density highways
response_diff_highways <- diff_function(response_highways_no_extra)
response_diff_highways <- response_diff_highways[!is.na(response_diff_highways$slope), ]
sum_diff_highways <- aggregate(slope ~ species, data = response_diff_highways, sum)
names(sum_diff_highways)[2] <- "slope_of_highways"

# Nighttimelights
response_diff_nighttimelights <- diff_function(response_nighttimelights_no_extra)
response_diff_nighttimelights <- response_diff_nighttimelights[!is.na(response_diff_nighttimelights$slope), ]
sum_diff_nighttimelights <- aggregate(slope ~ species, data = response_diff_nighttimelights, sum)
names(sum_diff_nighttimelights)[2] <- "slope_of_nighttimelights"

# Pasture
response_diff_pasture <- diff_function(response_pasture_no_extra)
response_diff_pasture <- response_diff_pasture[!is.na(response_diff_pasture$slope), ]
sum_diff_pasture <- aggregate(slope ~ species, data = response_diff_pasture, sum)
names(sum_diff_pasture)[2] <- "slope_of_pasture"

# Population
response_diff_population <- diff_function(response_population_no_extra)
response_diff_population <- response_diff_population[!is.na(response_diff_population$slope), ]
sum_diff_population <- aggregate(slope ~ species, data = response_diff_population, sum)
names(sum_diff_population)[2] <- "slope_of_population"

# Date: 04/18/2024
# Note: Based on the number of rows in sum_diff_crop, sum_diff_highways, sum_diff_pasture, 
#       there are some species that do not respond to cropland, highways, and pasture

# Manually check the species with no response to cropland, highways, and pasture
# The slope should be consistently zero for the response curves of these species
spplist[!spplist %in% sum_diff_pasture$species]
spplist[!spplist %in% sum_diff_crop$species]

spplist[!spplist %in% sum_diff_highways$species]

# Merge the slopes (difference) for the response curves of all human factors for 351 species
sum_diff_all <- merge(sum_diff_hii, sum_diff_crop, by = "species", all.x = TRUE)
sum_diff_all <- merge(sum_diff_all, sum_diff_highways, by = "species", all.x = TRUE)
sum_diff_all <- merge(sum_diff_all, sum_diff_nighttimelights, by = "species", all.x = TRUE)
sum_diff_all <- merge(sum_diff_all, sum_diff_pasture, by = "species", all.x = TRUE)
sum_diff_all <- merge(sum_diff_all, sum_diff_population, by = "species", all.x = TRUE)

# Sum over variable importance of all human factors versus climatic factors
var_imp_all <- perm_var_imp_rank_order
var_imp_all$category <- "Climatic"
var_imp_all$category[var_imp_all$variable %in% c("Human Impact", "Cropland", "Highway Density", "Nighttime Lights",
                                                 "Pasture", "Population")] <- "Human"

var_imp_all <- aggregate(value ~ species + category, data = var_imp_all, sum)

var_imp_clim <- var_imp_all[var_imp_all$category == "Climatic", ]
names(var_imp_clim)[3] <- "clim_varimp"
var_imp_human <- var_imp_all[var_imp_all$category == "Human", ]
names(var_imp_human)[3] <- "human_varimp"

# Merge species, order, variable importance of human impact index, 
# sum of slope for human factors, sum of variable importance of human and climatic factors
# Date: 05/02/2024
# Change: Only consider the sum of slope for human impact index
#         Since the number of x values varies among species, so the sum of slope for human impact index has to be adjusted by the number of
#         x values along the response curve of human impact index
spp_var_imp_slope <- merge(sum_diff_all[, c(1,2)], species_order, by.x = "species", by.y = "GBIF_name")
spp_var_imp_slope <- merge(spp_var_imp_slope, var_imp_clim[, c(1,3)], by = "species")
spp_var_imp_slope <- merge(spp_var_imp_slope, var_imp_human[, c(1,3)], by = "species")
spp_var_imp_slope <- merge(spp_var_imp_slope, perm_var_imp_order[perm_var_imp_order$variable == "Human Impact", 1:3], 
                           by = "species")
# Remove variable name of human impact index and rename its variable importance
spp_var_imp_slope <- spp_var_imp_slope[, -6]
names(spp_var_imp_slope)[6] <- "human_impact_varimp"

# Create a winner/loser category based on the sum of slopes of response curves for all human factors
# sum_of_slope_q <- quantile(spp_var_imp_slope$slope_of_hii)

# spp_var_imp_slope$category <- "Mid"
# spp_var_imp_slope$category[spp_var_imp_slope$sum_of_slope < sum_of_slope_q[2]] <- "Loser"
# spp_var_imp_slope$category[spp_var_imp_slope$sum_of_slope >= sum_of_slope_q[4]] <- "Winner"


# Date: 05/23/2024
# Change: Find the species with all positive slopes and all negative slopes
#         The species with all positive slopes --> minimum slope is > 0
#         The species with all negative slopes --> maximum slope is < 0

# Calculate the sum of slope (difference) for the response curves of all human factors
# Human impact index
# response_diff_hii <- diff_function(response_hii_no_extra)
# response_diff_hii <- response_diff_hii[!is.na(response_diff_hii$slope), ]
# sum_diff_hii <- aggregate(slope ~ species, data = response_diff_hii, sum)
# names(sum_diff_hii)[2] <- "slope_of_hii"
min_diff_hii <- aggregate(slope ~ species, data = response_diff_hii, min)
min_diff_hii$category_positive <- "0"
min_diff_hii$category_positive[min_diff_hii$slope > 0] <- "1"

max_diff_hii <- aggregate(slope ~ species, data = response_diff_hii, max)
max_diff_hii$category_negative <- "0"
max_diff_hii$category_negative[max_diff_hii$slope < 0] <- "1"

seg_spp <- as.data.frame(table(response_hii_order_no_extra$species))
spp_var_imp_slope <- merge(spp_var_imp_slope, seg_spp, by.x = "species", by.y = "Var1")

spp_var_imp_slope$slope_adjusted <- spp_var_imp_slope$slope_of_hii/spp_var_imp_slope$Freq

# Merge the category of all positive and negative slopes
spp_var_imp_slope <- merge(spp_var_imp_slope, min_diff_hii[, c(1,3)], by = "species")
spp_var_imp_slope <- merge(spp_var_imp_slope, max_diff_hii[, c(1,3)], by = "species")


sum_of_slope_q <- quantile(spp_var_imp_slope$slope_of_hii)
spp_var_imp_slope$category <- "Mid"
spp_var_imp_slope$category[spp_var_imp_slope$slope_of_hii < 0 | spp_var_imp_slope$category_negative == "1"] <- "Loser"
spp_var_imp_slope$category[spp_var_imp_slope$slope_of_hii >= 0 | spp_var_imp_slope$category_positive == "1"] <- "Winner"

# Date: 05/23/2024
# Change: Redefine the category by Loser = negative slope only and Winner = all positive slopes

# Export the data for trait analysis
write.csv(spp_var_imp_slope[, 1:9], "spp_var_imp_slope.csv", row.names = FALSE)

#############################################################################################################################
## Date: 04/21/2024
## Note: Create some plots for the results
##       1. Variable importance color coded by human versus climatic variables
##       2. Response curve of different human factors for all species
#############################################################################################################################
library(ggplot2)
# Variable importance
perm_var_imp_order_plot <- perm_var_imp_order
perm_var_imp_order_plot$Predictors <- "Human factor"
perm_var_imp_order_plot$Predictors[perm_var_imp_order_plot$variable %in% c("Mean Temperature of Warmest Quarter",
                                                                           "Mean Temperature of Coldest Quarter",
                                                                           "Annual Precipitation",
                                                                           "Precipitation of Wettest Quarter",
                                                                           "Precipitation of Driest Quarter",
                                                                           "Annual Mean Temperature")] <- "Climatic factor"
# Sort variable importance for human and climatic factors
clim_mean <- aggregate(value~variable + Predictors, 
                       data = perm_var_imp_order_plot[perm_var_imp_order_plot$Predictors == "Climatic factor", ], 
                       mean)
clim_mean <- clim_mean[order(clim_mean$value, decreasing = FALSE), ]

human_mean <- aggregate(value~variable + Predictors, 
                        data = perm_var_imp_order_plot[perm_var_imp_order_plot$Predictors == "Human factor", ], 
                        mean)
human_mean <- human_mean[order(human_mean$value, decreasing = FALSE), ]


perm_var_imp_order_plot$variable <- factor(perm_var_imp_order_plot$variable, levels = c(as.character(clim_mean$variable), as.character(human_mean$variable)))


# Stats summary helper function
# Function to produce summary statistics for lower and upper 5%
data_summary_sd <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x) #quantile(x, 0.95)
  ymax <- m + sd(x) #quantile(x, 0.05)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}  

data_summary <- function(x) {
  m <- median(x)
  ymin <- quantile(x, 0)
  ymax <- quantile(x, 1)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}

# Function to produce summary statistics for 1st and 3rd quantile
data_summary_1 <- function(x) {
  m <- median(x)
  ymin <- quantile(x, 0.25)
  ymax <- quantile(x, 0.75)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}

# Variable importance
# Overall
# Order by factors and variable importance

# Date: 05/23/2024
# Change: Rename population and nighttime lights as human population and artificial nighttime lights
predictor_label <- levels(perm_var_imp_order_plot$variable)
predictor_label <- c("bio17", "bio16", "bio12", "bio1", "bio10", "bio11", "POP", "NTL", "HD", "PAS", "CROP", "HII")

ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point",
               size = 2,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Permutation variable importance (%)") + 
  scale_x_discrete(guide = guide_axis(angle = 0), "", 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  #scale_color_manual(values = c("#7570B3", "#fb7b09")) +
  scale_color_manual(values = c("#4A7BB7", "#A50026")) +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_flip() +
  theme(legend.position = "bottom")

ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Permutation importance (%)") + 
  scale_x_discrete(name = " ",
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  # theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0.2, linewidth = 0.07) +
#   stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.6, color = "black") +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Permutation variable importance (%)") + 
#   scale_x_discrete(guide = guide_axis(angle = 0), "") +
#   scale_fill_manual(values = c("#7570B3", "#fb7b09")) +
#   #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   coord_flip() +
#   theme(legend.position = "bottom")

# By order
ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point",
               size = 1.5,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Permutation variable importance (%)") + 
  scale_x_discrete(guide = guide_axis(angle = 0), "", 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  #scale_color_manual(values = c("#7570B3", "#fb7b09")) +
  scale_color_manual(values = c("#4A7BB7", "#A50026")) +
  facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_flip() +
  theme(legend.position = "bottom")

# By order 800x600
ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Permutation importance (%)") + 
  scale_x_discrete(#guide = guide_axis(angle = 0), 
    "", 
    breaks = levels(perm_var_imp_order_plot$variable),
    labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  facet_wrap(. ~ order, ncol = 2, scales = "free") +
  # coord_flip() +
  #theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        strip.background = element_rect(color = "white", fill = "grey"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# By order 800x600 for human > clim
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Permutation importance (%)") + 
  scale_x_discrete(#guide = guide_axis(angle = 0), 
    "", 
    breaks = levels(perm_var_imp_order_plot$variable),
    labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
                     symnum.args = list(cutpoints = c(0, 0.001, Inf), 
                                        symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
                     vjust = 0.2) +
  facet_wrap(. ~ order, ncol = 2, scales = "free") +
  # coord_flip() +
  #theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        strip.background = element_rect(color = "white", fill = "grey"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors)) +
#   stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.2, linewidth = 0.07) +
#   stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.6, color = "black") +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Permutation variable importance (%)") + 
#   scale_x_discrete(guide = guide_axis(angle = 0), "") +
#   scale_fill_manual(values = c("#7570B3", "#fb7b09")) +
#   facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   coord_flip() +
#   theme(legend.position = "bottom")

# Rank of variable importance
perm_var_imp_rank_order_plot <- perm_var_imp_rank_order
perm_var_imp_rank_order_plot$Predictors <- "Human factor"
perm_var_imp_rank_order_plot$Predictors[perm_var_imp_rank_order_plot$variable %in% c("Mean Temperature of Warmest Quarter",
                                                                                     "Mean Temperature of Coldest Quarter",
                                                                                     "Annual Precipitation",
                                                                                     "Precipitation of Wettest Quarter",
                                                                                     "Precipitation of Driest Quarter",
                                                                                     "Annual Mean Temperature")] <- "Climatic factor"          


# Sort variable importance for human and climatic factors
clim_mean <- aggregate(value~variable + Predictors, 
                       data = perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$Predictors == "Climatic factor", ], 
                       mean)
clim_mean <- clim_mean[order(clim_mean$value, decreasing = FALSE), ]

human_mean <- aggregate(value~variable + Predictors, 
                        data = perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$Predictors == "Human factor", ], 
                        mean)
human_mean <- human_mean[order(human_mean$value, decreasing = FALSE), ]


perm_var_imp_rank_order_plot$variable <- factor(perm_var_imp_rank_order_plot$variable, levels = c(as.character(clim_mean$variable), as.character(human_mean$variable)))


# Overall
# Order by factors and variable importance
# Date: 05/23/2024
# Change: Rename population and nighttime lights as human population and artificial nighttime lights
predictor_rank_label <- levels(perm_var_imp_rank_order_plot$variable)
predictor_rank_label <-  c("bio17", "bio16", "bio12", "bio1", "bio10", "bio11", "POP", "NTL", "HD", "PAS", "CROP", "HII")

ggplot(perm_var_imp_rank_order_plot, aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Rank of permutation variable importance") + 
  scale_x_discrete(#guide = guide_axis(angle = 0), 
    name = "",
    breaks = levels(perm_var_imp_rank_order_plot$variable),
    labels = predictor_rank_label) +
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:1)+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# ggplot(perm_var_imp_rank_order_plot, aes(x = variable, y = rank, fill = Predictors)) + 
#   stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.2, linewidth = 0.07) +
#   stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.6, color = "black") +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Rank of permutation variable importance") + 
#   scale_x_discrete(guide = guide_axis(angle = 0), "") +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
#   scale_fill_manual(values = c("#7570B3", "#fb7b09")) +
#   coord_flip() +
#   theme(legend.position = "none")

# By order 800x600
ggplot(perm_var_imp_rank_order_plot, aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  #ylab("Rank of permutation variable importance") + 
  scale_x_discrete(name = "",
                   breaks = levels(perm_var_imp_rank_order_plot$variable),
                   labels = predictor_rank_label) +
  scale_y_reverse(name = "Rank of permutation importance", breaks = c(11, 9, 7, 5, 3, 1))+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  facet_wrap(. ~ order, ncol = 2, scales = "free") +
  # coord_flip() +
  #theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        strip.background = element_rect(color = "white", fill = "grey"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# ggplot(perm_var_imp_rank_order_plot, aes(x = variable, y = rank, color = Predictors)) + 
#   stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.2, linewidth = 0.07) +
#   stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.6, color = "black") +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Rank of permutation variable importance") + 
#   scale_x_discrete(guide = guide_axis(angle = 0), "") +
#   scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
#   scale_color_manual(values = c("#7570B3", "#fb7b09")) +
#   facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   coord_flip() +
#   theme(legend.position = "none")

# Visualize human impact index by order first but this may be improved by a good definition of winner versus losers
response_hii_order_no_extra_plot <- response_hii_order_no_extra
# Merge winner and losers based on sum of slopes
response_hii_order_no_extra_plot <- merge(response_hii_order_no_extra_plot, spp_var_imp_slope[, c("species", "category")])

# Wit winners/losers
ggplot(response_hii_order_no_extra_plot, aes(x, y, group = species)) + 
  geom_line(aes(color = order), alpha = 0.6, cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Probability of presence", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  facet_grid(order ~ category, scales = "free") +
  theme(legend.position = "none")

# Try visualizing response curve by winner/Mid/Loser x Order with gradients of variable importance of human impact index variable index
# Note: See line 770-832

perm_var_imp_order_human_impact <- perm_var_imp_order[perm_var_imp_order$variable == "Human Impact", ]
perm_var_imp_order_human_impact$variable <- "human_impact_index_2010"

response_hii_order_no_extra_plot_imp <- merge(response_hii_order_no_extra_plot, perm_var_imp_order_human_impact, 
                                              by = c("species", "variable", "order"))

ggplot(response_hii_order_no_extra_plot_imp, aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = order), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  facet_grid(order ~ category, scales = "free") +
  theme(legend.position = "none")

# With winners/losers defined by overall slope of response curves for all human factors
ggplot(response_hii_order_no_extra_plot, aes(x, y, group = species)) + 
  geom_line(aes(color = order, linetype = category, linewidth = category), alpha = 0.6) +
  scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous("Human Impact Index") +
  ylab("Probability of presence") +
  ylim(c(0,1)) +
  facet_wrap(. ~ order, ncol = 3, scales = "free") +
  theme(legend.position = "none")

#############################################################################################################################
## Date: 05/18/2024
## Note: Create bar plot to show the % of spp with winner, mid, losers by orders
## Change: Updated win, other, and lose based on new criteria
#############################################################################################################################
# figure2_ratio <- read.csv("D:/xinchen/SDM_human/R_script/22_figure2_ratio.csv")

figure2_ratio <- merge(perm_var_imp_order_human_impact[, c("species", "order")], win_lose_df, by = "species")
figure2_ratio <- aggregate(cbind(species) ~ order + win_group, data = figure2_ratio, FUN = length)
figure2_order_win_group <- data.frame(order = rep(sort(unique(figure2_ratio$order)), 3),
                                      win_group = rep(c("lose", "other", "win"), each = 7),
                                      n = rep(c(9, 20, 44, 2, 19, 223, 34), 3),
                                      ratio = rep(NA, 21))
figure2_ratio_plot <- merge(figure2_order_win_group, figure2_ratio, by = c("order", "win_group"))
# figure2_ratio_plot$species[is.na(figure2_ratio_plot$species)] <- 0
figure2_ratio_plot$ratio <- figure2_ratio_plot$species/figure2_ratio_plot$n * 100

# 700x400
ggplot(data = figure2_ratio_plot, aes(x = order, y = ratio, fill = win_group)) + 
  geom_col() +
  scale_fill_brewer(name = "Species response:", palette = "Dark2", 
                    labels = c("losers", "neutrals", "winners")) + 
  scale_x_discrete(labels = c("Artiodactyla\nn = 9", 
                              "Carnivora\nn = 20", 
                              "Chiroptera\nn = 44", 
                              "Didelphimorphia\nn = 2", 
                              "Lagomorpha\nn = 19", 
                              "Rodentia\nn = 223", 
                              "Soricomorpha\nn = 34")) +
  geom_text(aes(label = paste0(sprintf("%1.1f", ratio), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(x  = NULL, 
       y = "Proportion of species with different\nresponses to human impact (%)") +
  theme_classic() +
  theme(legend.position = "bottom")

#############################################################################################################################

#############################################################################################################################
## Date: 06/06/2024
## Note: Combine variable importance and rank of variable importance using cowplot library
#############################################################################################################################
library(cowplot)

# 820x280
plot_grid(NULL,
          ggplot(perm_var_imp_order_plot, aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            ylab("Permutation importance (%)") + 
            scale_x_discrete(name = " ",
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                             breaks = levels(perm_var_imp_order_plot$variable),
                             labels = predictor_label) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            # theme_classic() +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          ggplot(perm_var_imp_rank_order_plot, aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            ylab("Rank of permutation variable importance") + 
            scale_x_discrete(#guide = guide_axis(angle = 0), 
              name = "",
              breaks = levels(perm_var_imp_rank_order_plot$variable),
              labels = predictor_rank_label) +
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:1)+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()), 
          labels = c('a', '', 'b', ''), 
          ncol = 4, 
          rel_widths = c(0.05, 1, 0.05, 1),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)
#############################################################################################################################

# # Export the data for trait analysis
# write.csv(spp_var_imp_slope[, 1:9], "spp_var_imp_slope.csv", row.names = FALSE)

# Date: 06/06/2024
# Change: Merge newly defined winner, other, losers with the variable importance data frame
spp_var_imp_slope_new <- merge(spp_var_imp_slope[, 1:8], win_lose_df, by = "species")

# For sum of climatic factors > human factors
# Calculate the number by species
sum(spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp)
sum(spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp)

# Calculate the ratio by species
round(rbind(sum(spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp),
            sum(spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp))/rep(351, 2), 2)

# Calculate the number by order
rbind(table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp]), 
      table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp]))

# Calculate the ratio by order
round(rbind(table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp]), 
            table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp]))/rbind(table(spp_var_imp_slope_new$order), table(spp_var_imp_slope_new$order)), 2)

# Date: 06/13/2023
# Change: Test whether all climate is different from all human factors
# Visualize the sum of permutation importance for climate versus human factors
sum_clim_var_imp <- spp_var_imp_slope_new[, c("species", "order", "clim_varimp")]
names(sum_clim_var_imp)[3] <- "varimp"
sum_clim_var_imp$predictors <- "Sum of bioclimatic factors"

sum_human_var_imp <- spp_var_imp_slope_new[, c("species", "order", "human_varimp")]
names(sum_human_var_imp)[3] <- "varimp"
sum_human_var_imp$predictors <- "Sum of human factors"

sum_clim_human_var_imp <- rbind(sum_clim_var_imp, sum_human_var_imp)

library(ggpubr)
# For stat_summary do not use ylim or limits because this will throw data away
# Instead, use coord_cartesian(ylim = c(15, 30))
# Source: https://ggplot2.tidyverse.org/reference/stat_summary.html

# ggplot(sum_clim_human_var_imp, aes(x = predictors, y = varimp, fill = predictors, color = predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) +
#   #geom_boxplot() +
#   scale_y_continuous("Permutation importance (%)", breaks = c(0, 20, 40, 60, 80)) +
#   scale_x_discrete(guide = guide_axis(angle = 0), 
#                    name = "",
#                    labels = c("CLIM", "HUMAN")) +
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Sum of bioclimatic factors", "Sum of human factors")),
#                      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")),
#                      vjust = -0.2,
#                      label.y = 88) +
#   coord_cartesian(ylim = c(-5, 98)) +
#   #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none")

# Permutation importance 400x300
# Only one comparison between HII and bio11
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40, 60, 80)) + 
  scale_x_discrete(name = " ",
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")),
                     vjust = -0.5,
                     label.y = 36) +
  coord_cartesian(ylim = c(-5, 46)) +
  # theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Three comparisons between HII and bio11, HII and bio10, HII and bio1
# Run Wilcoxon signed ranked test with bonferroni adjustment for multiple comparisons
# Overall
stat.test.adj <- perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", 
                                                                                 "Mean Temperature of Coldest Quarter",
                                                                                 "Mean Temperature of Warmest Quarter",
                                                                                 "Annual Mean Temperature",
                                                                                 "Annual Precipitation"), ] %>% 
  pairwise_wilcox_test(value ~ variable, paired = TRUE, 
                       comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"), 
                                          c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                          c("Human Impact", "Annual Mean Temperature"),
                                          c("Human Impact", "Annual Precipitation"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
stat.test.adj

# Double check the test results
stat.test.adj.all <- pairwise_wilcox_test(perm_var_imp_order_plot, value ~ variable, paired = TRUE, p.adjust.method = "bonferroni")
stat.test.adj.all[stat.test.adj.all$group2 == "Human Impact", ]

#  .y.   group1                              group2          n1    n2 statistic        p    p.adj p.adj.signif
# 3 value Annual Precipitation                Human Impact   351   351    18444  2.73e- 9 1.80e- 7 ****        
# 4 value Annual Mean Temperature             Human Impact   351   351    23061  3   e- 3 2.14e- 1 ns          
# 5 value Mean Temperature of Warmest Quarter Human Impact   351   351    23298  3   e- 3 1.94e- 1 ns          
# 6 value Mean Temperature of Coldest Quarter Human Impact   351   351    32166. 1.76e- 1 1   e+ 0 ns  

# 600 x 400
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40, 60, 80), expand = expansion(mult = 0, 4)) + 
  scale_x_discrete(name = " ",
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                     p.adjust.method = "bonferroni",
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"),
                                        c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                        c("Human Impact", "Annual Mean Temperature"),
                                        c("Human Impact", "Annual Precipitation")),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                        symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                     vjust = -0.25,
                     label.y = c(49, 44, 39, 34), hide.ns = FALSE, tip.length = 0.01, bracket.size = .5, size = 3) +
  geom_text(x = 2.8, y = 56, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
  # stat_compare_means(method = "kruskal.test",
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
  #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
  #                    label.x = 2, label.y = 50) +
  coord_cartesian(ylim = c(-5, 55)) +
  # theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Order
stat.test.order.adj <- perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", 
                                                                                       "Mean Temperature of Coldest Quarter",
                                                                                       "Mean Temperature of Warmest Quarter",
                                                                                       "Annual Mean Temperature",
                                                                                       "Annual Precipitation"), ] %>% group_by(order) %>% 
  wilcox_test(value ~ variable, paired = TRUE,
              comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"), 
                                 c("Human Impact", "Annual Mean Temperature"),
                                 c("Human Impact", "Annual Precipitation"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
stat.test.order.adj %>% print(n = 28)


# Rank of permutation importance
# Overall
stat.test.rank.adj <- perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", 
                                                                                                "Mean Temperature of Coldest Quarter",
                                                                                                "Mean Temperature of Warmest Quarter",
                                                                                                "Annual Mean Temperature",
                                                                                                "Annual Precipitation"), ] %>% 
  pairwise_wilcox_test(rank ~ variable, paired = TRUE, 
                       comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"), 
                                          c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                          c("Human Impact", "Annual Mean Temperature"),
                                          c("Human Impact", "Annual Precipitation"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
stat.test.rank.adj


# Double check the test results
stat.test.rank.adj.all <- pairwise_wilcox_test(perm_var_imp_rank_order_plot, rank ~ variable, paired = TRUE, p.adjust.method = "bonferroni")
stat.test.rank.adj.all[stat.test.rank.adj.all$group2 == "Human Impact", ]

#  .y.   group1                              group2          n1    n2 statistic        p    p.adj p.adj.signif
# 3 rank  Annual Precipitation                Human Impact   351   351    43437  3.90e-11 2.57e- 9 ****        
# 4 rank  Annual Mean Temperature             Human Impact   351   351    42488. 9.99e-10 6.59e- 8 ****        
# 5 rank  Mean Temperature of Warmest Quarter Human Impact   351   351    38187  1.19e- 4 8   e- 3 **          
# 6 rank  Mean Temperature of Coldest Quarter Human Impact   351   351    29377  4.26e- 1 1   e+ 0 ns    

# Rank of permutation importance 400x300
ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species,], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  ylab("Rank of permutation variable importance") + 
  scale_x_discrete(#guide = guide_axis(angle = 0), 
    name = "",
    breaks = levels(perm_var_imp_rank_order_plot$variable),
    labels = predictor_rank_label) +
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0, expand = expansion(mult = 0, 1))+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  
  stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                     p.adjust.method = "bonferroni",
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"),
                                        c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                        c("Human Impact", "Annual Mean Temperature"),
                                        c("Human Impact", "Annual Precipitation")),
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf),
                                        symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                     vjust = -0.25,
                     label.y = c(2.6, 1.3, 0, -1.3), hide.ns = FALSE, tip.length = 0.01, bracket.size = .5, size = 3) +
  geom_text(x = 2.8, y = 3.7, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
  # stat_compare_means(method = "kruskal.test",
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
  #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
  #                    label.x = 2, label.y = 50) +
  coord_cartesian(ylim = c(12.5, -3.5)) +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())


# # Permutation importance by order 800x600
# ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species&perm_var_imp_order_plot$order != "Didelphimorphia", ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   # stat_summary(fun = mean, geom = "pointrange",
#   #              fun.max = function(x) mean(x) + sd(x),
#   #              fun.min = function(x) mean(x) - sd(x)) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50)) + 
#   scale_x_discrete(name = " ",
#                    #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
#                    breaks = levels(perm_var_imp_order_plot$variable),
#                    labels = predictor_label) +
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")),
#                      vjust = -0.5,
#                      label.y = 44) +
#   coord_cartesian(ylim = c(-12, 55)) +
#   # theme_classic() +
#   facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 9, face = "bold"),
#         axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -1),
#         panel.grid.minor = element_blank())

# # Rank of permutation importance by order 600x800
# ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species&perm_var_imp_rank_order_plot$order != "Didelphimorphia", ], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   # stat_summary(fun = mean, geom = "pointrange",
#   #              fun.max = function(x) mean(x) + sd(x),
#   #              fun.min = function(x) mean(x) - sd(x)) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Rank of permutation variable importance") + 
#   scale_x_discrete(#guide = guide_axis(angle = 0), 
#     name = "",
#     breaks = levels(perm_var_imp_rank_order_plot$variable),
#     labels = predictor_rank_label) +
#   scale_y_reverse(name = "Rank of permutation importance", breaks = 12:1)+
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")),
#                      vjust = -0.5,
#                      label.y = 0) +
#   coord_cartesian(ylim = c(12.5, -2)) +
#   # theme_classic() +
#   facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 9, face = "bold"),
#         axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -1),
#         panel.grid.minor = element_blank())


# Create a panel for permutation importance and its rank 1000x400
plot_grid(NULL,
          ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40, 60, 80), expand = expansion(mult = 0, 4)) + 
            scale_x_discrete(name = " ",
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                             breaks = levels(perm_var_imp_order_plot$variable),
                             labels = predictor_label) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                               p.adjust.method = "bonferroni",
                               comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"),
                                                  c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                                  c("Human Impact", "Annual Mean Temperature"),
                                                  c("Human Impact", "Annual Precipitation")),
                               symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                                  symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                               vjust = -0.25,
                               label.y = c(50-0.5, 45-0.5, 40-0.5, 35-0.5), hide.ns = FALSE, tip.length = 0.004, bracket.size = .5, size = 3) +
            geom_text(x = 2.8, y = 56, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
            # stat_compare_means(method = "kruskal.test",
            #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
            #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
            #                    label.x = 2, label.y = 50) +
            coord_cartesian(ylim = c(-5, 55)) +
            # theme_classic() +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species,], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            ylab("Rank of permutation variable importance") + 
            scale_x_discrete(#guide = guide_axis(angle = 0), 
              name = "",
              breaks = levels(perm_var_imp_rank_order_plot$variable),
              labels = predictor_rank_label) +
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0, expand = expansion(mult = 0, 1))+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            
            stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                               p.adjust.method = "bonferroni",
                               comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"),
                                                  c("Human Impact", "Mean Temperature of Warmest Quarter"),
                                                  c("Human Impact", "Annual Mean Temperature"),
                                                  c("Human Impact", "Annual Precipitation")),
                               symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf),
                                                  symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p < 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                               vjust = -0.25,
                               label.y = c(2.6, 1.3, 0, -1.3), hide.ns = FALSE, tip.length = 0.01, bracket.size = .5, size = 3) +
            geom_text(x = 2.8, y = 3.7, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
            # stat_compare_means(method = "kruskal.test",
            #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
            #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
            #                    label.x = 2, label.y = 50) +
            coord_cartesian(ylim = c(12.5, -3.5)) +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          labels = c('a', '', 'b', ''), 
          ncol = 4, 
          rel_widths = c(0.2, 2),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)


# This is wrong because this is not pairwise test!!!!!!!!!!!!!!!!!!
# Permutation importance by order
# stat.test.order.adj <- perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", 
#                                                                                        "Mean Temperature of Coldest Quarter",
#                                                                                        "Mean Temperature of Warmest Quarter",
#                                                                                        "Annual Mean Temperature",
#                                                                                        "Annual Precipitation"), ] %>% group_by(order) %>%
#   pairwise_wilcox_test(value ~ variable, paired = TRUE,
#                        comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"), 
#                                           c("Human Impact", "Mean Temperature of Warmest Quarter"),
#                                           c("Human Impact", "Annual Mean Temperature"),
#                                           c("Human Impact", "Annual Precipitation"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
# stat.test.order.adj %>% print(n = 28)

# Permutation importance by order 800x300
# Human impact versus Mean Temperature of Coldest Quarter
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Mean Temperature of Coldest Quarter"), ], 
       aes(x = order, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
           y.position = 100) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Mean Temperature of Warmest Quarter
# Permutation importance by order 800x300
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Mean Temperature of Warmest Quarter"), ], 
       aes(x = order, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
           y.position = 100) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Annual Mean Temperature
# Permutation importance by order 800x300
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Annual Mean Temperature"), ], 
       aes(x = order, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
           y.position = 100) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Annual Precipitation
# Permutation importance by order 800x300
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Annual Precipitation"), ], 
       aes(x = order, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
           y.position = 100) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())


# This is wrong because this is not pairwise test!!!!!!!!!!!!!!!!!!
# # Rank of permutation importance by order
# stat.test.rank.order.adj <- perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", 
#                                                                                                       "Mean Temperature of Coldest Quarter",
#                                                                                                       "Mean Temperature of Warmest Quarter",
#                                                                                                       "Annual Mean Temperature",
#                                                                                                       "Annual Precipitation"), ] %>% group_by(order) %>%
#   pairwise_wilcox_test(rank ~ variable, paired = TRUE,
#                        comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter"), 
#                                           c("Human Impact", "Mean Temperature of Warmest Quarter"),
#                                           c("Human Impact", "Annual Mean Temperature"),
#                                           c("Human Impact", "Annual Precipitation"))) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
# stat.test.rank.order.adj %>% print(n = 28)

# Human impact versus Mean Temperature of Coldest Quarter
ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Mean Temperature of Coldest Quarter"), ], 
       aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_cartesian(ylim = c(8, -1.6)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Mean Temperature of Warmest Quarter
ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Mean Temperature of Warmest Quarter"), ], 
       aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_cartesian(ylim = c(10, -2.3)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Annual Mean Temperature
ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Annual Mean Temperature"), ], 
       aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_cartesian(ylim = c(11, -2.3)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Human impact versus Annual Precipitation
ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Annual Precipitation"), ], 
       aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.8)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  scale_x_discrete(name = " "
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
  ) +
  geom_pwc(aes(group = variable), tip.length = 0,
           hide.ns = FALSE,
           method = "wilcox.test", 
           p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
           label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
  # theme_classic() +
  #facet_wrap(. ~ order, ncol = 2, scales = "free") +
  coord_cartesian(ylim = c(12, -2.3)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())


# Create a panel for permutation importance and its rank by order for HII versus bio11
plot_grid(NULL,
          ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Mean Temperature of Coldest Quarter"), ], 
                 aes(x = order, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
                     y.position = 100) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Mean Temperature of Coldest Quarter"), ], 
                 aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            coord_cartesian(ylim = c(8, -1.6)) +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          labels = c('a', '', 'b', ''), 
          ncol = 2, 
          rel_widths = c(0.2, 2),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)

# Create a panel for permutation importance and its rank by order for HII versus bio10
plot_grid(NULL,
          # Human impact versus Mean Temperature of Warmest Quarter
          # Permutation importance by order 800x300
          ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Mean Temperature of Warmest Quarter"), ], 
                 aes(x = order, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
                     y.position = 100) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          # Human impact versus Mean Temperature of Warmest Quarter
          ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Mean Temperature of Warmest Quarter"), ], 
                 aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            coord_cartesian(ylim = c(10, -2.3)) +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          labels = c('a', '', 'b', ''), 
          ncol = 2, 
          rel_widths = c(0.2, 2),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)

# Create a panel for permutation importance and its rank by order for HII versus bio1
plot_grid(NULL,
          # Human impact versus Annual Mean Temperature
          # Permutation importance by order 800x300
          ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Annual Mean Temperature"), ], 
                 aes(x = order, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
                     y.position = 100) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          # Human impact versus Annual Mean Temperature
          ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Annual Mean Temperature"), ], 
                 aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            coord_cartesian(ylim = c(11, -2.3)) +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          labels = c('a', '', 'b', ''), 
          ncol = 2, 
          rel_widths = c(0.2, 2),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)

# Create a panel for permutation importance and its rank by order for HII versus bio12
plot_grid(NULL,
          # Human impact versus Annual Precipitation
          # Permutation importance by order 800x300
          ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$variable %in% c("Human Impact", "Annual Precipitation"), ], 
                 aes(x = order, y = value, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.05, 0.25))) + 
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5,
                     y.position = 100) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          NULL,
          # Human impact versus Annual Precipitation
          ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$variable %in% c("Human Impact", "Annual Precipitation"), ], 
                 aes(x = order, y = rank, fill = Predictors, color = Predictors)) +
            stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8, position = position_dodge(0.8)) +
            # stat_summary(fun = mean, geom = "pointrange",
            #              fun.max = function(x) mean(x) + sd(x),
            #              fun.min = function(x) mean(x) - sd(x)) +
            stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
                         fun.max = function(x) mean(x) + sd(x),
                         fun.min = function(x) mean(x) - sd(x),
                         position = position_dodge(0.8)) +
            #stat_boxplot(geom = "errorbar", width = 0.25) + 
            #geom_boxplot() + 
            scale_y_reverse(name = "Rank of permutation importance", breaks = 12:0)+
            scale_color_manual(values = c("#2e4057", "#d1495b")) +
            scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
            scale_x_discrete(name = " "
                             #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5)
            ) +
            geom_pwc(aes(group = variable), tip.length = 0,
                     hide.ns = FALSE,
                     method = "wilcox.test", 
                     p.adjust.method = "bonferroni", method.args = list(paired = TRUE),
                     label = "Wilcoxon, \np = {p.adj.format}", label.size = 3.5, y.position = -0.5) +
            # theme_classic() +
            #facet_wrap(. ~ order, ncol = 2, scales = "free") +
            coord_cartesian(ylim = c(12, -2.3)) +
            theme(legend.position = "none",
                  legend.title = element_text(size = 9, face = "bold"),
                  axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
                  axis.title.y = element_text(vjust = +3),
                  axis.title.x = element_text(vjust = -1),
                  panel.grid.minor = element_blank()),
          labels = c('a', '', 'b', ''), 
          ncol = 2, 
          rel_widths = c(0.2, 2),
          hjust = -0.1, vjust = 1,
          label_fontface = "bold",
          label_size = 14)


# # For human factors < clim factors
# sum_clim_var_imp_1 <- spp_var_imp_slope_new[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp, c("species", "order", "clim_varimp")]
# names(sum_clim_var_imp_1)[3] <- "varimp"
# sum_clim_var_imp_1$predictors <- "Sum of bioclimatic factors"
# 
# sum_human_var_imp_1 <- spp_var_imp_slope_new[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp, c("species", "order", "human_varimp")]
# names(sum_human_var_imp_1)[3] <- "varimp"
# sum_human_var_imp_1$predictors <- "Sum of human factors"
# 
# sum_clim_human_var_imp_1 <- rbind(sum_clim_var_imp_1, sum_human_var_imp_1)
# 
# ggplot(sum_clim_human_var_imp_1, aes(x = predictors, y = varimp, fill = predictors, color = predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) +
#   #geom_boxplot() +
#   scale_y_continuous("Permutation importance (%)", breaks = c(20, 40, 60, 80)) +
#   scale_x_discrete(guide = guide_axis(angle = 0), 
#                    name = "",
#                    labels = c("CLIM", "HUMAN")) +
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Sum of bioclimatic factors", "Sum of human factors")),
#                      symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                      vjust = -0.2,
#                      label.y = 90) +
#   coord_cartesian(ylim = c(10, 100)) +
#   #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none")
# 
# # Permutation importance 400x300
# ggplot(perm_var_imp_order_plot[!perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   # stat_summary(fun = mean, geom = "pointrange",
#   #              fun.max = function(x) mean(x) + sd(x),
#   #              fun.min = function(x) mean(x) - sd(x)) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40)) + 
#   scale_x_discrete(name = " ",
#                    #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
#                    breaks = levels(perm_var_imp_order_plot$variable),
#                    labels = predictor_label) +
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                      symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                      vjust = -0.2,
#                      label.y = 38) +
#   coord_cartesian(ylim = c(-5, 45)) +
#   # theme_classic() +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 9, face = "bold"),
#         axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -1),
#         panel.grid.minor = element_blank())
# 
# # Permutation importance by order 800x600
# ggplot(perm_var_imp_order_plot[!perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   # stat_summary(fun = mean, geom = "pointrange",
#   #              fun.max = function(x) mean(x) + sd(x),
#   #              fun.min = function(x) mean(x) - sd(x)) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   scale_y_continuous(name = "Permutation importance (%)", breaks = c(20, 40, 60, 80, 100)) + 
#   scale_x_discrete(name = " ",
#                    #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
#                    breaks = levels(perm_var_imp_order_plot$variable),
#                    labels = predictor_label) +
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#   coord_cartesian(ylim = c(-5, 76)) +
#   # theme_classic() +
#   facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 9, face = "bold"),
#         axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -1),
#         panel.grid.minor = element_blank())
# 
# # Rank of permutation importance 400x300
# ggplot(perm_var_imp_rank_order_plot[!perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species,], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
#   stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#   # stat_summary(fun = mean, geom = "pointrange",
#   #              fun.max = function(x) mean(x) + sd(x),
#   #              fun.min = function(x) mean(x) - sd(x)) +
#   stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                fun.max = function(x) mean(x) + sd(x),
#                fun.min = function(x) mean(x) - sd(x)) +
#   #stat_boxplot(geom = "errorbar", width = 0.25) + 
#   #geom_boxplot() + 
#   ylab("Rank of permutation variable importance") + 
#   scale_x_discrete(#guide = guide_axis(angle = 0), 
#     name = "",
#     breaks = levels(perm_var_imp_rank_order_plot$variable),
#     labels = predictor_rank_label) +
#   scale_y_reverse(name = "Rank of permutation importance", breaks = 12:1)+
#   scale_color_manual(values = c("#2e4057", "#d1495b")) +
#   scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#   stat_compare_means(method = "wilcox.test", 
#                      comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                      symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                         symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                      vjust = -0.2,
#                      label.y = -.6) +
#   coord_cartesian(ylim = c(12.5, -.4)) +
#   #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 9, face = "bold"),
#         axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -1),
#         panel.grid.minor = element_blank())
# 
# 
# # Create a 2x3 panel
# plot_grid(NULL,
#           ggplot(sum_clim_human_var_imp, aes(x = predictors, y = varimp, fill = predictors, color = predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) +
#             #geom_boxplot() +
#             scale_y_continuous("Permutation importance (%)", breaks = c(20, 40, 60, 80)) +
#             scale_x_discrete(guide = guide_axis(angle = 0), 
#                              name = "",
#                              labels = c("CLIM", "HUMAN")) +
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Sum of bioclimatic factors", "Sum of human factors")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = 85) +
#             coord_cartesian(ylim = c(10, 95)) +
#             #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#             theme(legend.position = "none"),
#           NULL,
#           # Permutation importance 400x300
#           ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             # stat_summary(fun = mean, geom = "pointrange",
#             #              fun.max = function(x) mean(x) + sd(x),
#             #              fun.min = function(x) mean(x) - sd(x)) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) + 
#             #geom_boxplot() + 
#             scale_y_continuous(name = "Permutation importance (%)", breaks = c(20, 40, 60, 80)) + 
#             scale_x_discrete(name = " ",
#                              #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
#                              breaks = levels(perm_var_imp_order_plot$variable),
#                              labels = predictor_label) +
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = 65) +
#             coord_cartesian(ylim = c(-5, 76)) +
#             # theme_classic() +
#             theme(legend.position = "none",
#                   legend.title = element_text(size = 9, face = "bold"),
#                   axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#                   axis.title.y = element_text(vjust = +3),
#                   axis.title.x = element_text(vjust = -1),
#                   panel.grid.minor = element_blank()),
#           NULL,
#           ggplot(perm_var_imp_rank_order_plot[perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species,], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             # stat_summary(fun = mean, geom = "pointrange",
#             #              fun.max = function(x) mean(x) + sd(x),
#             #              fun.min = function(x) mean(x) - sd(x)) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) + 
#             #geom_boxplot() + 
#             ylab("Rank of permutation variable importance") + 
#             scale_x_discrete(#guide = guide_axis(angle = 0), 
#               name = "",
#               breaks = levels(perm_var_imp_rank_order_plot$variable),
#               labels = predictor_rank_label) +
#             scale_y_reverse(name = "Rank of permutation\nimportance", breaks = 12:1)+
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = 0) +
#             coord_cartesian(ylim = c(12.5, -1.6)) +
#             #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#             theme(legend.position = "none",
#                   legend.title = element_text(size = 9, face = "bold"),
#                   axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#                   axis.title.y = element_text(vjust = +3),
#                   axis.title.x = element_text(vjust = -1),
#                   panel.grid.minor = element_blank()),
#           NULL,
#           ggplot(sum_clim_human_var_imp_1, aes(x = predictors, y = varimp, fill = predictors, color = predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) +
#             #geom_boxplot() +
#             scale_y_continuous("Permutation importance (%)", breaks = c(20, 40, 60, 80)) +
#             scale_x_discrete(guide = guide_axis(angle = 0), 
#                              name = "",
#                              labels = c("CLIM", "HUMAN")) +
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Sum of bioclimatic factors", "Sum of human factors")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = 90) +
#             coord_cartesian(ylim = c(10, 100)) +
#             #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#             theme(legend.position = "none"),
#           NULL,
#           ggplot(perm_var_imp_order_plot[!perm_var_imp_order_plot$species %in% sum_human_var_imp$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             # stat_summary(fun = mean, geom = "pointrange",
#             #              fun.max = function(x) mean(x) + sd(x),
#             #              fun.min = function(x) mean(x) - sd(x)) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) + 
#             #geom_boxplot() + 
#             scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40)) + 
#             scale_x_discrete(name = " ",
#                              #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
#                              breaks = levels(perm_var_imp_order_plot$variable),
#                              labels = predictor_label) +
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = 38) +
#             coord_cartesian(ylim = c(-5, 47)) +
#             # theme_classic() +
#             theme(legend.position = "none",
#                   legend.title = element_text(size = 9, face = "bold"),
#                   axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#                   axis.title.y = element_text(vjust = +3),
#                   axis.title.x = element_text(vjust = -1),
#                   panel.grid.minor = element_blank()),
#           NULL,
#           ggplot(perm_var_imp_rank_order_plot[!perm_var_imp_rank_order_plot$species %in% sum_human_var_imp$species,], aes(x = variable, y = rank, fill = Predictors, color = Predictors)) +
#             stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
#             # stat_summary(fun = mean, geom = "pointrange",
#             #              fun.max = function(x) mean(x) + sd(x),
#             #              fun.min = function(x) mean(x) - sd(x)) +
#             stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
#                          fun.max = function(x) mean(x) + sd(x),
#                          fun.min = function(x) mean(x) - sd(x)) +
#             #stat_boxplot(geom = "errorbar", width = 0.25) + 
#             #geom_boxplot() + 
#             ylab("Rank of permutation variable importance") + 
#             scale_x_discrete(#guide = guide_axis(angle = 0), 
#               name = "",
#               breaks = levels(perm_var_imp_rank_order_plot$variable),
#               labels = predictor_rank_label) +
#             scale_y_reverse(name = "Rank of permutation\nimportance", breaks = 12:1)+
#             scale_color_manual(values = c("#2e4057", "#d1495b")) +
#             scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
#             stat_compare_means(method = "wilcox.test", 
#                                comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
#                                symnum.args = list(cutpoints = c(0, 0.001, Inf), 
#                                                   symbols = c("Wilcoxon, p < 0.001", "Wilcoxon, p > 0.001")),
#                                tip.length = 0.02,
#                                vjust = -0.2,
#                                label.y = -0.6) +
#             coord_cartesian(ylim = c(12.5, -1)) +
#             #facet_wrap(. ~ order, ncol = 2, scales = "free") +
#             theme(legend.position = "none",
#                   legend.title = element_text(size = 9, face = "bold"),
#                   axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
#                   axis.title.y = element_text(vjust = +3),
#                   axis.title.x = element_text(vjust = -1),
#                   panel.grid.minor = element_blank()), 
#           labels = c('a', '', 'b', '', 'c', '', 'd', '', 'e', '', 'f', ''), 
#           ncol = 6, 
#           rel_widths = c(0.09, 1.1, 0.09, 2, 0.09, 2),
#           hjust = -0.1, vjust = 1,
#           label_fontface = "bold",
#           label_size = 14)


# Calculate the ratio by species
round(rbind(sum(spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp),
            sum(spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp))/rep(351, 2), 2)

# Calculate the number by order
rbind(table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp]), 
      table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp]))

# Calculate the ratio by order
round(rbind(table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp]), 
            table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp]))/rbind(table(spp_var_imp_slope_new$order), table(spp_var_imp_slope_new$order)), 2)

# Crate a dataframe to record 
order_human_greater_clim <- rbind(t(round(rbind(sum(spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp),
                                                sum(spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp)), 2)),
                                  t(rbind(table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp > spp_var_imp_slope_new$human_varimp]), 
                                          table(spp_var_imp_slope_new$order[spp_var_imp_slope_new$clim_varimp <= spp_var_imp_slope_new$human_varimp]))))
row.names(order_human_greater_clim)[1] <- "Total"
order_human_greater_clim <- as.data.frame(order_human_greater_clim)
order_human_greater_clim <- order_human_greater_clim[order(order_human_greater_clim$V2, decreasing = FALSE), ]
#order_human_greater_clim$perc <- paste0(sprintf("%4.1f", round(order_human_greater_clim$V2/(order_human_greater_clim$V1 + order_human_greater_clim$V2)*100, 1)), "%")
order_human_greater_clim$perc <- round(order_human_greater_clim$V2/(order_human_greater_clim$V1 + order_human_greater_clim$V2)*100, 1)
names(order_human_greater_clim)[2] <- "n"
# order_human_greater_clim$col <- factor(row.names(order_human_greater_clim), levels = c("Total\nn = 351",
#                                                                                        "Rodentia\nn = 223",
#                                                                                        "Chiroptera\nn = 44", 
#                                                                                        "Carnivora\nn = 20", 
#                                                                                        "Soricomorpha\nn = 34",
#                                                                                        "Lagomorpha\nn = 19", 
#                                                                                        "Artiodactyla\nn = 9", 
#                                                                                        "Didelphimorphia\nn = 2"))

order_human_greater_clim$col <- factor(row.names(order_human_greater_clim), levels = rev(c("Artiodactyla",
                                                                                           "Carnivora", 
                                                                                           "Chiroptera",
                                                                                           "Didelphimorphia",
                                                                                           "Lagomorpha",
                                                                                           "Rodentia",
                                                                                           "Soricomorpha",
                                                                                           "Total")))

# Source: Bar chart for numbers with ratios
# Bar chart 400x300
#label_perc <- paste0(sprintf("%4.1f", round(order_human_greater_clim$perc, 1)), "%")
label_perc <- sprintf("%4.1f", round(order_human_greater_clim$perc, 1))
ggplot(order_human_greater_clim, aes(x = perc, y = col)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = label_perc), 
            hjust = 0, nudge_x = -4,
            fontface = "bold", size = 4) +
  scale_x_continuous("\nProportion of species with sum of permutation\nimportance for HUMAN > CLIM factors (%)",
                     breaks = c(0, 10, 20, 30, 40, 50),
                     labels = c(0, 10, 20, 30, 40, 50)) +
  scale_y_discrete("") +
  # Make sure labels do not get cut, Part 1
  coord_cartesian(xlim = c(0, 50), clip = "off") +
  theme_minimal() +
  # Make sure labels do not get cut, Part 2
  theme(plot.margin = margin(r = 25))


# Find all the species with the rank of permutation importance = 1
perm_var_imp_rank_one_order <- perm_var_imp_rank_order[perm_var_imp_rank_order$rank == 1, ]

# Create a column to distinguish human factors versus climatic factors
perm_var_imp_rank_one_order$Predictors <- "Human factor"
perm_var_imp_rank_one_order$Predictors[perm_var_imp_rank_one_order$variable %in% c("Mean Temperature of Warmest Quarter",
                                                                                   "Mean Temperature of Coldest Quarter",
                                                                                   "Annual Precipitation",
                                                                                   "Precipitation of Wettest Quarter",
                                                                                   "Precipitation of Driest Quarter",
                                                                                   "Annual Mean Temperature")] <- "Climatic factor"
# Calculate the number by species
table(perm_var_imp_rank_one_order$Predictors)

# Calculate the ratio by species
round(table(perm_var_imp_rank_one_order$Predictors)/351, 2)

# Calculate the number by order
table(perm_var_imp_rank_one_order$Predictors, perm_var_imp_rank_one_order$order)

# Calculate the ratio by order
round(table(perm_var_imp_rank_one_order$Predictors, perm_var_imp_rank_one_order$order)/rbind(colSums(table(perm_var_imp_rank_one_order$Predictors, perm_var_imp_rank_one_order$order)), 
                                                                                             colSums(table(perm_var_imp_rank_one_order$Predictors, perm_var_imp_rank_one_order$order))), 2)

# Export a table showing the top 1 predictor in each order
write.csv(table(perm_var_imp_rank_one_order$variable, perm_var_imp_rank_one_order$order), "top_pred_by_order.csv", row.names = TRUE)


# Crate a dataframe to record 
order_human_greater_clim_rank <- rbind(table(perm_var_imp_rank_one_order$Predictors),
                                       t(table(perm_var_imp_rank_one_order$Predictors, perm_var_imp_rank_one_order$order)))

row.names(order_human_greater_clim_rank)[1] <- "All species"
order_human_greater_clim_rank <- as.data.frame(order_human_greater_clim_rank)
names(order_human_greater_clim_rank) <- c("V1", "n")
order_human_greater_clim_rank <- order_human_greater_clim_rank[order(order_human_greater_clim_rank$n, decreasing = FALSE), ]
# order_human_greater_clim_rank$perc <- paste0(sprintf("%4.1f", round(order_human_greater_clim_rank$n/(order_human_greater_clim_rank$V1 + order_human_greater_clim_rank$n)*100, 1)), "%")
order_human_greater_clim_rank$perc <- round(order_human_greater_clim_rank$n/(order_human_greater_clim_rank$V1 + order_human_greater_clim_rank$n)*100, 1)
# order_human_greater_clim_rank$col <- factor(row.names(order_human_greater_clim_rank), levels = c("Total\nn = 351",
#                                                                                                  "Rodentia\nn = 223",
#                                                                                                  "Chiroptera\nn = 44",
#                                                                                                  "Carnivora\nn = 20",
#                                                                                                  "Soricomorpha\nn = 34",
#                                                                                                  "Lagomorpha\nn = 19",
#                                                                                                  "Artiodactyla\nn = 9",
#                                                                                                  "Didelphimorphia\nn = 2"))

order_human_greater_clim_rank$col <- factor(row.names(order_human_greater_clim_rank), levels = rev(c("Didelphimorphia",
                                                                                                     "Chiroptera",
                                                                                                     "Carnivora",
                                                                                                     "Rodentia",
                                                                                                     "Artiodactyla",
                                                                                                     "Lagomorpha",
                                                                                                     "Soricomorpha",
                                                                                                     "All species")))

# Bar chart 400x300
label_perc_rank <- sprintf("%4.1f", round(order_human_greater_clim_rank$perc, 1))
ggplot(order_human_greater_clim_rank, aes(x = perc, y = col)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = label_perc_rank), 
            hjust = 0, nudge_x = -4,
            fontface = "bold", size = 4) +
  scale_x_continuous("\nProportion of species with HUMAN factor as the\nmost important predictor",
                     breaks = c(0, 10, 20, 30, 40, 50),
                     labels = c(0, 10, 20, 30, 40, 50)) +
  scale_y_discrete("") +
  # Make sure labels do not get cut, Part 1
  coord_cartesian(xlim = c(0, 50), clip = "off") +
  theme_minimal() +
  # Make sure labels do not get cut, Part 2
  theme(plot.margin = margin(r = 25))


# Date: 06/18/2024
# Change: Summarize the percent of most important human factors and each most important factor in each order because we will conclude human impact based on the most important predictor
# Note: This logic is similar to that in The Importance of the Human Footprint in Shaping the Global Distribution of Terrestrial, Freshwater and Marine Invaders
# Source: https://doi.org/10.1371/journal.pone.0125801

# Load manually decast data frame
top_pred_by_order_ratio_plot <- read.csv("Preliminary_results/2024/3rd/3_top_pred_by_order_ratio_plot.csv")
top_pred_by_order_ratio_plot$Ratio <- 100*top_pred_by_order_ratio_plot$Ratio

top_pred_by_order_ratio_plot_overall <- top_pred_by_order_ratio_plot[top_pred_by_order_ratio_plot$Order == "Overall", ]
label_perc_top_pred <- c(sprintf("%4.1f", round(top_pred_by_order_ratio_plot_overall$Ratio, 1))[1:9], rep("", 3))
top_pred_by_order_ratio_plot_overall$Abbv <- factor(top_pred_by_order_ratio_plot_overall$Abbv, levels = c(top_pred_by_order_ratio_plot_overall$Abbv))
top_pred_by_order_ratio_plot$Group <- as.factor(top_pred_by_order_ratio_plot$Group)

ggplot(top_pred_by_order_ratio_plot_overall, aes(x = Ratio, y = rev(Abbv), fill = Group)) +
  geom_col() +
  geom_text(aes(label = label_perc_top_pred), 
            hjust = 0, nudge_x = -1.5,
            fontface = "bold", size = 4, color = "white") +
  scale_x_continuous("\nProportion of species for each\nmost important predictor",
                     breaks = c(0, 10, 20, 30),
                     labels = c(0, 10, 20, 30)) +
  scale_y_discrete("",
                   label = rev(top_pred_by_order_ratio_plot_overall$Abbv)) +
  scale_fill_manual(values = c("#2e4057", "#d1495b")) +
  # Make sure labels do not get cut, Part 1
  coord_cartesian(xlim = c(0, 30), clip = "off") +
  #facet_wrap(. ~ Order, row = 1) +
  theme_minimal() +
  # Make sure labels do not get cut, Part 2
  theme(legend.position = "none",
        plot.margin = margin(r = 25))

# Create a 1x2 panel 800x300
plot_grid(NULL,
          ggplot(order_human_greater_clim_rank, aes(x = perc, y = col)) +
            geom_col(fill = "grey70") +
            geom_text(aes(label = label_perc_rank), 
                      hjust = 0, nudge_x = -5,
                      fontface = "bold", size = 3) +
            scale_x_continuous("\nProportion of species with HUMAN factor as the\nmost important predictor",
                               breaks = c(0, 10, 20, 30, 40, 50),
                               labels = c(0, 10, 20, 30, 40, 50)) +
            scale_y_discrete("") +
            # Make sure labels do not get cut, Part 1
            coord_cartesian(xlim = c(0, 50), clip = "off") +
            theme_minimal() +
            # Make sure labels do not get cut, Part 2
            theme(plot.margin = margin(r = 25)),
          NULL,
          ggplot(top_pred_by_order_ratio_plot_overall, aes(x = Ratio, y = rev(Abbv), fill = Group)) +
            geom_col() +
            geom_text(aes(label = label_perc_top_pred), 
                      hjust = 0, nudge_x = -2.5,
                      fontface = "bold", size = 3, color = "white") +
            scale_x_continuous("\nProportion of species for each\nmost important predictor",
                               breaks = c(0, 10, 20, 30),
                               labels = c(0, 10, 20, 30)) +
            scale_y_discrete("",
                             label = rev(top_pred_by_order_ratio_plot_overall$Abbv)) +
            scale_fill_manual(values = c("#2e4057", "#d1495b")) +
            # Make sure labels do not get cut, Part 1
            coord_cartesian(xlim = c(0, 30), clip = "off") +
            #facet_wrap(. ~ Order, row = 1) +
            theme_minimal() +
            # Make sure labels do not get cut, Part 2
            theme(legend.position = "none",
                  plot.margin = margin(r = 25)),
          labels = c('a', '', 'b', ''), 
          ncol = 4, 
          rel_widths = c(0.02, 2, 0.02, 2),
          hjust = 0, vjust = 0.9,
          label_fontface = "bold",
          label_size = 14)

top_pred_by_order_ratio_plot_overall$Abbv <- factor(top_pred_by_order_ratio_plot_overall$Abbv, levels = c(top_pred_by_order_ratio_plot_overall$Abbv))

# By order
top_pred_by_order_ratio_plot$Abbv <- factor(top_pred_by_order_ratio_plot$Abbv, levels = c(top_pred_by_order_ratio_plot_overall$Abbv))

ggplot(top_pred_by_order_ratio_plot[top_pred_by_order_ratio_plot$Order != "Overall", ], aes(x = Ratio, y = rev(Abbv), fill = Group)) +
  geom_col() +
  scale_x_continuous("\nPercentage of species for each\nmost important predictor (%)") +
  scale_y_discrete("",
                   label = rev(top_pred_by_order_ratio_plot_overall$Abbv)) +
  scale_fill_manual(values = c("#2e4057", "#d1495b")) +
  # Make sure labels do not get cut, Part 1
  facet_wrap(. ~ Order, scales = "free") +
  theme_minimal() +
  # Make sure labels do not get cut, Part 2
  theme(legend.position = "none",
        plot.margin = margin(r = 10)) +
  geom_text(aes(label = c(sprintf("%4.1f", round(top_pred_by_order_ratio_plot$Ratio[top_pred_by_order_ratio_plot$Order != "Overall"], 1)))), 
            hjust = 1, nudge_x = 0,
            fontface = "bold", size = 2.5, color = "white") 



# Date: 06/24/2024
# Change: Summarize AUC values and some other metrics

auc_list <- data.frame()

for (i in 1:length(spplist)) {
  auc_i <- data.frame(AUC = maxent_list_1st[[i]]@results[5])
  auc_i$species <- spplist[i]
  auc_list <- rbind(auc_list, auc_i)
}

# Calculate mean and sd of AUC values
round(mean(auc_list$AUC), 2)
round(sd(auc_list$AUC), 2)


sensitivity_list <- data.frame()

for (i in 1:length(spplist)) {
  sensitivity_i <- data.frame(Sens = maxent_list_1st[[i]]@results[84])
  sensitivity_i$species <- spplist[i]
  sensitivity_list <- rbind(sensitivity_list, sensitivity_i)
}

# Calculate mean and sd of AUC values
round(mean(1-sensitivity_list$Sens), 2)
round(sd(1-sensitivity_list$Sens), 2)


# Date: 06/27/2024
# Note: Find some species based on the slope of their response curves

# Find some indicator species for winner, neutrals, and losers groups
# Loser in Rodentia order
win_lose_df[win_lose_df$win_group == "lose", ]$species[win_lose_df[win_lose_df$win_group == "lose", ]$species %in% df_curve_abs1_plot[df_curve_abs1_plot$order == "Rodentia", "species"]]

# Marmota_olympus
ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Marmota_olympus", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")

# Arborimus longicaudus 
ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Arborimus_longicaudus", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")

# Dipodomys_ingens
ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Dipodomys_ingens", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")


# Loser in Soricomorpha order
win_lose_df[win_lose_df$win_group == "lose", ]$species[win_lose_df[win_lose_df$win_group == "lose", ]$species %in% df_curve_abs1_plot[df_curve_abs1_plot$order == "Soricomorpha", "species"]]

# Sorex dispar 

# Sorex_longirostris

# Neutrals
# Need to create a slope function similar to diff_function 

################################################################################################################################
# Calculate the slope based on the lowest point and highest point 
# Create a function to calculate the slope between highest and lowest point
slop_function <- function(data){
  data_save <- data.frame()
  # Loop through all species j
  for (j in 1:length(unique(data$species))) {
    data_j <- data[data$species == unique(data$species)[j], ]
    
    # Decide if the slope is positive or negative based on the relative position of highest and lowest point
    # If the highest point is on the left of the lowest point, slope is negative
    # If the highest point is on the right of the lowest point, slope is positive
    
    max_x_i <- data_j$x[which.max(data_j$y)]
    min_x_i <- data_j$x[which.min(data_j$y)]
    max_i <- max(data_j$y)
    min_i <- min(data_j$y)
    
    # Recale the slope from the range of its corresponding HII to 0-100
    slope_i <- ((max_i - min_i)/(max_x_i - min_x_i))*(abs(max(data_j$x) - min(data_j$x)))
    
    # Save the calculated slope
    data_save[j, "slope"] <- slope_i
    
    # Save the species name
    data_save[j, "species"] <- unique(data$species)[j]
    
    # # Update data_j to the blank data frame
    # data_save <- rbind(data_save, data_j)
  }
  return(data_save)
}
################################################################################################################################
# Calculate rescaled slope for each species
slope_all <- slop_function(response_hii_order_no_extra)
hist(slope_all$slope)

slope_all <- merge(slope_all, spp_var_imp_slope_new[, c("species", "order", "win_group")], by = "species", all.x = TRUE)

# Find the species with the smallest slope in Soricomorpha
slope_all_soricomorpha <- slope_all[slope_all$order == "Soricomorpha", ]

slope_all_soricomorpha[which.min(abs(slope_all_soricomorpha$slope)), ]

# Sorex preblei

ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Sorex_preblei", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")

# Find the species with the smallest slope in Rodentia
slope_all_rodentia <- slope_all[slope_all$order == "Rodentia", ]

slope_all_rodentia[which.min(abs(slope_all_rodentia$slope)), ]

# Zapus trinotatus
ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Zapus_trinotatus", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")

slope_all_rodentia[which.min(abs(slope_all_rodentia$slope[slope_all_rodentia$species != "Zapus_trinotatus"])), ]

# Peromyscus pectoralis
ggplot(df_curve_abs1_plot[df_curve_abs1_plot$species == "Peromyscus_pectoralis", ], aes(x, y, group = species, alpha = value)) + 
  geom_line(aes(color = win_group), cex = 0.8) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_brewer(palette = "Dark2", 
                      labels = c("losers", "neutrals", "winners")) +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "none")


# Merge all species occurrences by order to create a map showing high coverage of species occurrence overlaid with human impact
# Merge all species occurrences
p_data_frame <- data.frame()
occlist_updated <- list.files("D:/xinchen/SDM_human/R_script/TestRun_NA_1990_2020/occurrences",
                              pattern = ".csv", full.names = TRUE)

for (i in 1:length(occlist_updated)) {
  temp_occ <- read.csv(occlist_updated[i])
  temp_occ$species <- spplist[i]
  p_data_frame <- rbind(p_data_frame, temp_occ)
}

nrow(p_data_frame)

p_data_frame <- merge(p_data_frame, spp_var_imp_slope_new[, c("species", "order")], by = "species")

# Calculate the number of thinned background points for each species
bg_list_updated <- list.files("D:/xinchen/SDM_human/R_script/TestRun_NA_1990_2020/bg",
                              pattern = ".csv", full.names = TRUE)
bg_data_frame <- spp_var_imp_slope_new[, c("species", "order")]
bg_data_frame$bg <- NA

for (i in 1:length(bg_list_updated)) {
  temp_bg_thin <- read.csv(bg_list_updated[i])
  bg_data_frame$bg[i] <- nrow(temp_bg_thin)
}

# Summarize the number of occurrences for each species
p_data_frame_sum <- aggregate(cbind(x) ~ species + order, data = p_data_frame, FUN = length)

# Merge the number of thinned background points with the number of presences
p_data_frame_sum <- merge(p_data_frame_sum, bg_data_frame, by = c("species", "order"))

# Merge AUC with the number of thinned presences
p_data_frame_sum <- merge(p_data_frame_sum, auc_list, by = "species")
p_data_frame_sum$AUC <- round(p_data_frame_sum$AUC, 2)


# Export data for making maps
write.csv(p_data_frame, "p_data_frame.csv", row.names = FALSE)

write.csv(p_data_frame_sum[, -4], "p_data_number_species_by_order.csv", row.names = FALSE)

######################################################################################################################################################################################
# Date: 07/03/2024
# Note: Add some additional comparison for Discussion
# Compare the results of Roland's paper (https://onlinelibrary.wiley.com/doi/10.1111/ddi.13900)
# Load the species list
species_list_Roland <- read.csv("D:/xinchen/SDM_human/R_script/Preliminary_results/2024/3rd/Roland_result/species_list_Roland.csv")
sum(species_list_Roland$species %in% spplist)

# Double check the test results for the species in Roland's paper
stat.test.Roland.adj.all <- pairwise_wilcox_test(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% species_list_Roland$species, ], 
                                                 value ~ variable, paired = TRUE, p.adjust.method = "bonferroni")
stat.test.Roland.adj.all[stat.test.Roland.adj.all$group2 == "Human Impact", ]


# Variable importance
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% species_list_Roland$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40, 60, 80), expand = expansion(mult = 0, 10)) + 
  scale_x_discrete(name = " ",
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                     p.adjust.method = "bonferroni",
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                        symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                     vjust = -0.25,
                     label.y = c(66, 44, 39, 34), hide.ns = FALSE, tip.length = 0.01, bracket.size = .5, size = 3) +
  geom_text(x = 2.8, y = 75, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
  # stat_compare_means(method = "kruskal.test",
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
  #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
  #                    label.x = 2, label.y = 72) +
  coord_cartesian(ylim = c(-5, 68)) +
  # theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Response curve
# Manual palettee source: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
ggplot(df_curve_abs1_plot_Suraci[df_curve_abs1_plot_Suraci$response == 1, ], aes(x, y, group = species), alpha = 0.5) + 
  geom_line(aes(color = common_name), size = 1.1, alpha = 1) +
  #geom_textline(aes(label = common_name), linecolour = "NA", textcolour = "black", hjust = 1) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  #scale_colour_manual(values = c("#D95F02", "#7570B3"),
  #                    labels = c("neutrals", "winners")) +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#7570B3", "#52854C","#D95F02","#56B4E9"), name = "Species") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "right")


# Compare the results of Suraci's paper (https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.15650)
# Load the species list
species_list_Suraci <- read.csv("D:/xinchen/SDM_human/R_script/Preliminary_results/2024/3rd/comparable_studies/Suraci_et_al_results/species_list_Suraci.csv")
sum(species_list_Suraci$species %in% spplist)


# Double check the test results for the species in Suraci's paper
stat.test.Suraci.adj.all <- pairwise_wilcox_test(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% species_list_Suraci$species, ], 
                                                 value ~ variable, paired = TRUE, p.adjust.method = "bonferroni")
stat.test.Suraci.adj.all[stat.test.Suraci.adj.all$group2 == "Human Impact", ]

# Variable importance
ggplot(perm_var_imp_order_plot[perm_var_imp_order_plot$species %in% species_list_Suraci$species, ], aes(x = variable, y = value, fill = Predictors, color = Predictors)) +
  stat_summary(fun.data = data_summary_sd, geom = "errorbar", width = 0, linewidth = .8) +
  # stat_summary(fun = mean, geom = "pointrange",
  #              fun.max = function(x) mean(x) + sd(x),
  #              fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, stroke = 0.8,
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  #stat_boxplot(geom = "errorbar", width = 0.25) + 
  #geom_boxplot() + 
  scale_y_continuous(name = "Permutation importance (%)", breaks = c(0, 20, 40, 60, 80), expand = expansion(mult = 0, 10)) + 
  scale_x_discrete(name = " ",
                   #guide = guide_axis(angle = 30, vjust = 0, hjust = 0.5), 
                   breaks = levels(perm_var_imp_order_plot$variable),
                   labels = predictor_label) +
  scale_color_manual(values = c("#2e4057", "#d1495b")) +
  scale_fill_manual(name = "Group", values = c("#2e4057", "#d1495b")) +
  stat_compare_means(method = c("wilcox.test"), paired = TRUE,
                     p.adjust.method = "bonferroni",
                     comparisons = list(c("Human Impact", "Mean Temperature of Coldest Quarter")),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                        symbols = c("Wilcoxon, p < 0.01", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05", "Wilcoxon, p > 0.05")), # Make sure the range of p.adj matches with the added labels for p
                     vjust = -0.25,
                     label.y = c(66, 44, 39, 34), hide.ns = FALSE, tip.length = 0.01, bracket.size = .5, size = 3) +
  geom_text(x = 2.8, y = 75, label = "Kruskal-Wallis, p < 0.01", color = "black", size = 3.5) +
  # stat_compare_means(method = "kruskal.test",
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
  #                                       symbols = c("Kruskal-Wallis, p < 0.001", "Kruskal-Wallis, p < 0.01", "Kruskal-Wallis, p < 0.05", "Kruskal-Wallis, p > 0.05")),
  #                    label.x = 2, label.y = 72) +
  coord_cartesian(ylim = c(-5, 68)) +
  # theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -1),
        panel.grid.minor = element_blank())

# Response curve
library(geomtextpath)
# Source: https://stackoverflow.com/questions/29357612/plot-labels-at-ends-of-lines
df_curve_abs1_plot_Suraci <- merge(df_curve_abs1_plot[df_curve_abs1_plot$species %in% species_list_Suraci$species, ], species_list_Suraci, by = c("species", "order"))

# Plot the 6 common species in Figure 3b
ggplot(df_curve_abs1_plot_Suraci[df_curve_abs1_plot_Suraci$response == 1, ], aes(x, y, group = species)) + 
  geom_line(aes(color = common_name)) +
  #geom_textline(aes(label = common_name), linecolour = "NA", textcolour = "black", hjust = 1) +
  #scale_linewidth_manual(values = c(0.6, 1)) +
  #scale_linetype_manual(values = c("dashed", "solid")) +
  #scale_colour_manual(values = c("#D95F02", "#7570B3"),
  #                    labels = c("neutrals", "winners")) +
  scale_colour_brewer(palette = "Dark2", name = "Species") +
  scale_x_continuous("Human Impact Index", breaks = c(0, 2500, 5000)) +
  scale_y_continuous("Relative probability of presence", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(order ~ win_group, scales = "free", 
             labeller = labeller(win_group = win_group_labs)) +
  theme(legend.position = "right")

## end of line
ggplot(temp.dat) +
  geom_textline(aes(
    x = Year, y = Capex, group = State, colour = State, label = State
  ),
  hjust = 1
  ) +
  theme(legend.position = "none")













# Run some quick aggregations

# Number of species in winners and losers by order
aggregate(cbind(species) ~ order + category, data = spp_var_imp_slope, length)

# Mean variable importance for all climatic vars and human factors and sum of slopes for all human factors by order
aggregate(cbind(clim_varimp, human_varimp, sum_of_slope) ~ order, data = spp_var_imp_slope, mean)

# Find the most important human factor and climatic factor for each order
# aggregate(cbind(value) ~ order, data = perm_var_imp_order, max)
# 
# 
# perm_var_imp_order[perm_var_imp_order$order == "Artiodactyla" & perm_var_imp_order$value == 53.9045, ]
# perm_var_imp_order[perm_var_imp_order$order == "Carnivora" & perm_var_imp_order$value == 92.5638, ]
# perm_var_imp_order[perm_var_imp_order$order == "Chiroptera" & perm_var_imp_order$value == 65.3399, ]
# perm_var_imp_order[perm_var_imp_order$order == "Didelphimorphia" & perm_var_imp_order$value == 85.0542, ]
# perm_var_imp_order[perm_var_imp_order$order == "Lagomorpha" & perm_var_imp_order$value == 86.5552, ]
# perm_var_imp_order[perm_var_imp_order$order == "Rodentia" & perm_var_imp_order$value == 98.2340, ]
# perm_var_imp_order[perm_var_imp_order$order == "Soricomorpha" & perm_var_imp_order$value == 90.7049, ]


mean_var_imp_order <- aggregate(cbind(value) ~ order + variable, data = perm_var_imp_order, mean)
mean_var_imp_order <- reshape2::melt(mean_var_imp_order, id.vars = "order", 
                                     variable.name = "variable",
                                     value.name = "value")

# Human factors
mean_var_imp_order_human <- mean_var_imp_order[!mean_var_imp_order$variable %in% c("Mean Temperature of Warmest Quarter",
                                                                                   "Mean Temperature of Coldest Quarter",
                                                                                   "Annual Precipitation",
                                                                                   "Precipitation of Wettest Quarter",
                                                                                   "Precipitation of Driest Quarter",
                                                                                   "Annual Mean Temperature"), ]

mean_var_imp_order_human <- reshape2::dcast(mean_var_imp_order_human, order~variable)
mean_var_imp_order_human[, c(2:7)] <- round(mean_var_imp_order_human[, c(2:7)] , 2)

# Climate factors
mean_var_imp_order_clim <- mean_var_imp_order[mean_var_imp_order$variable %in% c("Mean Temperature of Warmest Quarter",
                                                                                 "Mean Temperature of Coldest Quarter",
                                                                                 "Annual Precipitation",
                                                                                 "Precipitation of Wettest Quarter",
                                                                                 "Precipitation of Driest Quarter",
                                                                                 "Annual Mean Temperature"), ]

mean_var_imp_order_clim <- reshape2::dcast(mean_var_imp_order_clim, order~variable)
mean_var_imp_order_clim[, c(2:7)] <- round(mean_var_imp_order_clim[, c(2:7)] , 2)

# Export these two tables
write.csv(mean_var_imp_order_human, "Preliminary_results/2024/Table2_var_imp_human.csv", row.names = FALSE)
write.csv(mean_var_imp_order_clim, "Preliminary_results/2024/Table2_var_imp_clim.csv", row.names = FALSE)



# Use sum of slope to define the winner and loser
sum_diff <- response_diff[!is.na(response_diff$slope), ]
sum_diff <- aggregate(slope ~ species, data = sum_diff, sum)

hist(sum_diff$slope, main = "Histogram of slopes for human impact")
round(quantile(sum_diff$slope), digits = 2)

# Merge response of human impact with species order
response_hii_order <- merge(response_hii, species_order, by.x = "species", by.y = "GBIF_name", all.x = TRUE)
