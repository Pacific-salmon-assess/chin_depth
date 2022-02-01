### Explore Blocking Options
## Jan. 28, 2021

library(tidyverse)
library(blockCV)
library(raster)
library(sf)

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) 

# number of blocks shared among all methods
n_blocks <- 10

set.seed(456)


## SPATIAL BLOCKING ------------------------------------------------------------

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

# check
depth_dat_raw$space_block <- as.factor(sb$foldID)

coast <- readRDS(here::here("data", "crop_coast_sf.RDS"))
ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = depth_dat_raw, 
             aes(x = longitude, y = latitude, fill = space_block),
             shape = 21) +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "black"))


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

# check
depth_dat_raw$time_block <- as.factor(time_folds$time_block)
#subsample tags to make visualization easier
sub_tags <- sample(depth_dat_raw$vemco_code, size = 50, replace = FALSE)

ggplot(depth_dat_raw %>% filter(vemco_code %in% sub_tags)) +
  geom_point(aes(x = vemco_code, y = date_time_local, fill = time_block),
             shape = 21) +
  facet_wrap(~year, scales = "free")


## INDIVIDUAL BLOCKING ---------------------------------------------------------

# make vector of 5 blocks
ind_folds = data.frame(
  vemco_code = unique(depth_dat_raw$vemco_code),
  ind_block = sample.int(n_blocks, length(unique(depth_dat_raw$vemco_code)), 
                     replace = T) %>% as.factor()) 

# check
left_join(depth_dat_raw, ind_folds, by = "vemco_code")  %>% 
  dplyr::select(vemco_code, ind_block) %>% 
  distinct()


## EXPORT BLOCKING IDS ---------------------------------------------------------

block_list <- list(space = sp_folds, time = time_folds, individual = ind_folds)
saveRDS(block_list, 
        here::here("data", "10block_ids.RDS"))


## COMPARE BLOCKING STRUCTURES -------------------------------------------------

library(caret)
library(recipes)
library(gbm)


# subset to relevant data and split by block type
depth_list <- depth_dat_raw %>% 
  mutate(logit_rel_depth = qlogis(rel_depth)) %>% 
  left_join(., ind_folds, by = "vemco_code") %>%
  dplyr::select(
    logit_rel_depth, stage, latitude, longitude, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    space = space_block, time = time_block, ind = ind_block
  ) %>% 
  pivot_longer(., cols = space:ind, names_to = "block_type") %>% 
  split(., .$block_type)

# save as tibble
block_tbl <- tibble(
  block = c("space", "time", "ind"),
  full_dat = depth_list
) %>% 
  mutate(
    train_dat = purrr::map(full_dat, function (x) 
      x %>% filter(!value == "5") %>% droplevels
      ),
    test_dat = purrr::map(full_dat, function (x) 
      x %>% filter(value == "5") %>% droplevels
    ),
    depth_recipe = purrr::map(train_dat, function (x) 
      recipe(logit_rel_depth ~ ., 
             data = x %>% dplyr::select(-block_type, -value)) %>% 
        step_nzv(all_predictors()) %>% 
        step_dummy(all_predictors(), -all_numeric())
      ),
    # identify blocks in training data
    depth_ctrl = purrr::map(train_dat, function (x) {
      train_folds <- groupKFold(x$value, k = length(unique(x$value)))
      trainControl(
        method="repeatedcv",
        index = train_folds
      )
    })
    ) 

# check recipe
prep(block_tbl$depth_recipe[[1]]) %>%
  bake(., 
       new_data = block_tbl$train_dat[[1]] %>% 
         dplyr::select(-block_type, -value)
       ) %>%
  glimpse()


# parallel based on OS
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}


# fit GBMs (random forest substantially slower and only moderately better 
# performance)
gbmGrid <-  expand.grid(interaction.depth = c(5, 10), #c(3, 5, 9),
                        n.trees = seq(20, 350, by = 10),
                        shrinkage = 0.1,
                        n.minobsinnode = c(5, 10, 20))

tictoc::tic()
gbm_list <- purrr::pmap(
  list(block_tbl$depth_recipe,
       block_tbl$train_dat,
       block_tbl$depth_ctrl),
  function (x, y, z) {
    train(x, y, method = "gbm", metric = "RMSE", maximize = FALSE, 
          trControl = z, 
          # tuneLength = 3
          tuneGrid = gbmGrid
          )  
  }
)
names(gbm_list) <- block_tbl$block
tictoc::toc()


# visualize 
trellis.par.set(caretTheme())
purrr::map(gbm_list, plot) 

pred_foo <- function(mod, dat) {
  preds <- predict(mod, newdata = dat)
  dat$logit_preds <- preds
  
  par(mfrow = c(2, 1))
  plot(logit_preds ~ logit_rel_depth, data = dat)
  abline(0, 1, col = "red")
  plot(plogis(logit_preds) ~ plogis(logit_rel_depth), data = dat)
  abline(0, 1, col = "red")
}

# predictions, even with the training data, not very promising
pred_foo(gbm_list[[1]], dat = block_tbl$train_dat[[1]])
pred_foo(gbm_list[[2]], dat = block_tbl$train_dat[[2]])
pred_foo(gbm_list[[3]], dat = block_tbl$train_dat[[3]])

pred_foo(gbm_list[[3]], dat = block_tbl$test_dat[[3]])


pred_foo(depth_rf, dat = train_depth)
pred_foo(depth_rf, dat = test_depth)
