### Explore Blocking Options
## Jan. 28, 2021


library(tidyverse)
library(blockCV)
library(raster)
library(sf)

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) 

# n_blocks (shared among all methods)
n_blocks <- 5


## SPATIAL BLOCKING ------------------------------------------------------------

coast <- readRDS(here::here("data", "coast_sf.RDS"))
bath_list <- readRDS(here::here("data", "bathy_lowres_rasters.RDS"))

bath <- merge(bath_list[[1]], bath_list[[2]]) %>% 
  merge(., bath_list[[3]])
crs(bath) <- "+proj=longlat +datum=WGS84" 

# reproject raster
bath_utm <- projectRaster(bath, crs = "+proj=utm +zone=10 +units=m")

depth_sf <- st_as_sf(depth_dat_raw, coords = c("longitude", "latitude"), 
                     crs = "+proj=longlat +datum=WGS84") %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))

plot(bath_utm[[1]])
plot(depth_sf, pch = 16, col = "red", add = TRUE)


# spatial correlation in predictors
# TODO: currently only uses bathymetry (others were non-sensical)
sac <- spatialAutoRange(rasterLayer = bath_list[[1]],
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)

# test spatial blockings
sb <- spatialBlock(speciesData = depth_sf,
                   species = NULL,
                   rasterLayer = bath_utm,
                   theRange = 55000, # size of the blocks
                   k = n_blocks,
                   selection = "systematic",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)

# visualize location of points
sb$plots + geom_sf(data = depth_sf, alpha = 0.5)

sp_folds <- sb$foldID


## TEMPORAL BLOCKING -----------------------------------------------------------

# previous analyses indicate temporal autocorrelation on the span of several 
# hours, bin accordingly within an individual (units = minutes)
depth_dat_raw$timestamp_f = cut_width(depth_dat_raw$timestamp_n, 
                                      width = 180, boundary = -0.1)
depth_dat_raw$id = paste(depth_dat_raw$vemco_code, depth_dat_raw$timestamp_f,
                         sep = "_")
                         
# make vector of 5 blocks
time_folds = data.frame(
  id = unique(depth_dat_raw$id),
  time_block = sample.int(n_blocks, length(unique(depth_dat_raw$id)), 
                     replace = T)) %>% 
  right_join(., depth_dat_raw %>% dplyr::select(id), by = "id") 


## INDIVIDUAL BLOCKING ---------------------------------------------------------

# make vector of 5 blocks
ind_folds = data.frame(
  vemco_code = unique(depth_dat_raw$vemco_code),
  block = sample.int(n_blocks, length(unique(depth_dat_raw$vemco_code)), 
                     replace = T)) 



## LIST OF BLOCKING IDS --------------------------------------------------------

block_list <- list(space = sp_folds, time = time_folds, individual = ind_folds)
saveRDS(block_list, 
        here::here("data", "5block_ids.RDS"))


# blockCV Example --------------------------------------------------------------

# import raster data
awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))

# import presence-absence species data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pa_data <- st_as_sf(PA, coords = c("x", "y"), crs = crs(awt))
# see the first few rows
pa_data

plot(awt[[1]]) # plot raster data
plot(pa_data[which(pa_data$Species==1), ], pch = 16, col="red", add=TRUE) # add presence points
plot(pa_data[which(pa_data$Species==0), ], pch = 16, col="blue", add=TRUE) # add absence points
legend(x=500000, y=8250000, legend=c("Presence","Absence"), col=c(2, 4), pch=c(16,16), bty="n")

# import presence-background species data
PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pb_data <- st_as_sf(PB, coords = c("x", "y"), crs = crs(awt))
# number of presence and background records
table(pb_data$Species)

## spatial blocks
# spatial blocking by specified range with systematic assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species",
                   rasterLayer = awt,
                   theRange = 70000, # size of the blocks
                   k = 5,
                   selection = "systematic",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)

# visualize location of points
sb$plots + geom_sf(data = pa_data, alpha = 0.5)

# spatial blocking by rows and columns with checkerboard assignment
sb2 <- spatialBlock(speciesData = pb_data, # presence-background data
                    species = "Species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 6,
                    k = 5,
                    selection = "systematic",
                    biomod2Format = TRUE)


## buffering
# buffering with presence-absence data (all data used)
bf1 <- buffering(speciesData = pa_data,
                 theRange = 70000,
                 species = "Species", # to count the number of presences and absences/backgrounds
                 spDataType = "PA", # presence-absence  data type
                 progress = TRUE)
bf1$plots


## environmental block
eb <- envBlock(rasterLayer = awt,
               speciesData = pa_data,
               species = "Species",
               k = 5,
               standardization = "standard", # rescale variables between 0 and 1
               rasterBlock = FALSE,
               numLimit = 50)


## spatial AC in predictors
sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)


## visualization tools
foldExplorer(blocks = sb, 
             rasterLayer = awt, 
             speciesData = pa_data)

# explore the block size
rangeExplorer(rasterLayer = awt) # the only mandatory input

# add species data to add them on the map
rangeExplorer(rasterLayer = awt,
              speciesData = pa_data,
              species = "Species",
              rangeTable = NULL,
              minRange = 30000, # limit the search domain
              maxRange = 100000)


## example model
# extract the raster values for the species points as a dataframe
mydata <- raster::extract(awt, pb_data)
mydata <- as.data.frame(mydata)
# create a vector of 1 (for presence) and 0 (for background samples)
pb <- pb_data$Species

# extract the folds in spatialBlock object created 
# in the previous section (with presence-background data)
# the foldID only works for spatialBlock and envBlock folds
folds <- sb2$foldID

# create an empty vector to store the AUC of each fold
AUCs <- vector(mode = "numeric")
for(k in seq_len(5)){
  # extracting the training and testing indices
  # this way only works with foldID
  trainSet <- which(folds != k) # training set indices
  testSet <- which(folds == k) # testing set indices
  # fitting a maxent model using linear, quadratic and hinge features
  mx <- maxnet(p = pb[trainSet], 
               data = mydata[trainSet, ], 
               maxnet.formula(p = pb[trainSet], 
                              data = mydata[trainSet, ], 
                              classes = "default"))
  testTable <- pb_data[testSet, ] # a table for testing predictions and reference data
  testTable$pred <- predict(mx, mydata[testSet, ], type = "cloglog") # predict the test set
  # calculate area under the ROC curve
  precrec_obj <- evalmod(scores = testTable$pred, labels = testTable$Species)
  AUCs[k] <- auc(precrec_obj)[1,4] # extract AUC-ROC
}

# print the mean of AUCs
print(mean(AUCs))


## BLOCKING BY SUBJECT IN CARET ------------------------------------------------


library(caret)
library(nlme)

data(Orthodont)
head(Orthodont)
subjects <- as.character(unique(Orthodont$Subject))

## figure out folds at the subject level

set.seed(134)
sub_folds <- createFolds(y = subjects, list = TRUE, returnTrain = TRUE)
sub_folds2 <- groupKFold(group = subjects, k = length(unique(subjects)))

## now create the mappings to which *rows* are in the training set
## based on which subjects are left in or out

in_train <- holdout <- vector(mode = "list", length = length(sub_folds))

row_index <- 1:nrow(Orthodont)

for(i in seq(along = sub_folds)) {
  ## Which subjects are in fold i
  sub_in <- subjects[sub_folds[[i]]]
  ## which rows of the data correspond to those subjects
  in_train[[i]] <- row_index[Orthodont$Subject %in% sub_in]
  holdout[[i]]  <- row_index[!(Orthodont$Subject %in% sub_in)]  
}

names(in_train) <- names(holdout) <- names(sub_folds)

ctrl <- trainControl(method = "cv",
                     savePredictions = TRUE,
                     index = in_train,
                     indexOut = holdout)

mod <- train(distance ~ (age+Sex)^2, data = Orthodont,
             method = "lm", 
             trControl = ctrl)

first_fold <- subset(mod$pred, Resample == "Fold01")

## These were used to fit the model
table(Orthodont$Subject[-first_fold$rowIndex])
## These were heldout:
table(Orthodont$Subject[first_fold$rowIndex])



