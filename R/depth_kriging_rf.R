## Fit Random Forest Regression Kriging

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Fit equivalent model with interpolated training data and including RF
#regression kriging to better account for spatial covariance
#following Fox et al. 2020 PloS One.
#RF model structure


library(tidyverse)
library(randomForest)
library(quantregForest)
library(slmrf)
library(caret)
library(recipes)


depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(stage = as.factor(stage))


## PRELIMINARY CLEANING --------------------------------------------------------

## block and interpolate data as in depth_quantreg_rf.R

# add individual block
set.seed(1234)
ind_folds <- data.frame(
  vemco_code = unique(depth_dat_raw$vemco_code),
  ind_block =  sample.int(
    8, length(unique(depth_dat_raw$vemco_code)), replace = T
  ) %>% 
    as.factor()
)

depth_dat <- depth_dat_raw %>% 
  left_join(., ind_folds, by = "vemco_code") %>% 
  dplyr::select(
    depth = pos_depth, fl, mean_log_e, stage, utm_x, utm_y, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  #impute missing ROMS values
  step_impute_knn(all_predictors(), neighbors = 3) %>%
  # step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric())


train_folds <- groupKFold(train_depth$ind_block,
                          k = length(unique(train_depth$ind_block)))
depth_ctrl <-   trainControl(
  method="repeatedcv",
  index = train_folds
)


# refit top model in random forest, using quantregforest to generate uncertainty
# intervals

# apply recipe to dataframe to interpolate values
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))
saveRDS(train_depth_baked, 
        here::here("data", "baked_training_data.RDS"))


# fit standard RF version of quant-reg model using hyperpars based on exp
# analyses
# rf_refit <- readRDS(here::here("data", "model_fits", "depth_quantrf.rds"))
rf_refit <- randomForest::randomForest(depth ~ ., 
                                       data = train_depth_baked, 
                                       mtry = 6, 
                                       ntree = 1000)


# generate distance matrix for coordinates
coord_mat <- train_depth %>% 
  select(utm_x, utm_y) %>%
  distinct() %>% 
  as.matrix()
dist_mat <- compute_distance(coord_mat)

fit <- fit_rfrk(rf_refit, train_depth$depth, dist_mat)


# check computation time with subset
dum <- sample_n(train_depth_baked, 4000, replace = FALSE) %>% 
  select(-(agg_Col:agg_WA_OR))

coord_mat <- dum %>% 
  distinct() %>% 
  as.matrix()
dist_mat <- compute_distance(coord_mat)

dum_fit <- randomForest::randomForest(depth ~ ., 
                                       data = dum, 
                                       mtry = 6, 
                                       ntree = 500)

fit <- fit_rfrk(dum_fit, dum$depth, dist_mat)


