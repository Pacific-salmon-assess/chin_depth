### Depth Model Caret -- Model Comparison
## April 8, 2022
## Compete different model structure to identify best predictor
## Multiple dimensions: 
## 1) Model type (GBM vs. RF)
## 2) Hyperparameters
## 3) Response variable (logit transformed, relative depth, absolute depth)
## 4) Model structure (separate models for stages or common model)

library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)


# use 15 minute bins for fitting due to large number of observations
depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_15min.RDS")) %>% 
  mutate(logit_rel_depth = qlogis(rel_depth),
         stage = as.factor(stage))


# check no infinite values in transformed relative depth
nrow(depth_dat_raw[is.infinite(depth_dat_raw$logit_rel_depth), ])
depth_dat_raw <- depth_dat_raw[is.finite(depth_dat_raw$logit_rel_depth), ]


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
    depth = pos_depth, rel_depth, logit_rel_depth, latitude, longitude, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, stage, ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels() 
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels() 


## MODELS ----------------------------------------------------------------------

## define hyperparameters
gbm_grid <-  expand.grid(interaction.depth = c(2, 5, 10), #c(3, 5, 9),
                         n.trees = c(seq(10, 100, by = 10), 150, 200, 250, 300),
                         shrinkage = 0.1,
                         n.minobsinnode = c(5, 10, 20))

# for RF two dimensions
rf_grid <- expand.grid(tune_length = 8,
                       n.trees = seq(100, 1000, by = 100))


## define stage-specific datasets and consolidate all into tbls
stage_foo <- function(dat, stage_in) {
  dat %>% 
    filter(stage == stage_in) %>% 
    rename(stage_dat = stage)
}

train_mat <- stage_foo(train_depth, stage = "mature")
test_mat <- stage_foo(test_depth, stage = "mature")
train_imm <- stage_foo(train_depth, stage = "immature")
test_imm <- stage_foo(test_depth, stage = "immature")


## create dataframe of models to compete
model_tbl <- expand.grid(
  model_type = c("gbm", "rf"),
  response = c("logit_rel_depth", "rel_depth", "depth"),
  stage_dat = c("integrated", "mature", "immature")
  ) %>% 
  as_tibble() %>% 
  mutate(
    hyper_pars = ifelse(model_type == "gbm", gbm_grid, rf_grid)
  )

# function to add data based on input variables
add_data <- function(response, stage_dat) {
  # select dataset
  if (stage_dat == "integrated") {
    train_dum <- train_depth %>% 
      dplyr::select(-ind_block)
    test_dum <- test_depth %>% 
      dplyr::select(-ind_block)
  } 
  if (stage_dat %in% c("mature", "immature")) {
    train_dum <- stage_foo(train_depth, stage = stage_dat) %>% 
      dplyr::select(-ind_block)
    test_dum <- stage_foo(test_depth, stage = stage_dat) %>% 
      dplyr::select(-ind_block)
  }
  
  # subset based on response
  list(
    train = train_dum %>% select(depth_var = response, latitude:ind_block),
    test = test_dum %>% select(depth_var = response, latitude:ind_block)
  )
}

dum <- map2(model_tbl$response %>% as.character(), 
            model_tbl$stage_dat %>% as.character(), 
            .f = add_data)
model_tbl$train_data <- map(dum, function (x) x$train)
model_tbl$test_data <- map(dum, function (x) x$test)

## add recipe 
model_tbl$recipe <- map2(
  model_tbl$train_data,
  function (x) {
    dum_in <- x %>% dplyr::select(-ind_block)
    
    recipe(depth_var ~ ., data = x) %>% 
      step_impute_knn(all_predictors(), neighbors = 3) %>%
      step_nzv(all_predictors()) %>% 
      step_dummy(all_predictors(), -all_numeric())
})
depth_recipe <- recipe(logit_rel_depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  step_impute_knn(all_predictors(), neighbors = 3) %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric())


