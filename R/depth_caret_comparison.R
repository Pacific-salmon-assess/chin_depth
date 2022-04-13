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


# parallelize based on operating system
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}


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

## create dataframe of models to compete
model_tbl <- expand.grid(
  model_type = c("gbm", "rf"),
  response = c("logit_rel_depth", "rel_depth", "depth"),
  stage_dat = c("integrated", "mature", "immature")
  ) %>% 
  as_tibble() 

# function to add data based on input variables
add_data <- function(response, stage_dat) {
  # select dataset
  if (stage_dat == "integrated") {
    train_dum <- train_depth %>% 
      dplyr::select(depth_var = response, latitude:stage, ind_block)
    test_dum <- test_depth %>% 
      dplyr::select(depth_var = response, latitude:stage, ind_block)
  } 
  if (stage_dat %in% c("mature", "immature")) {
    train_dum <- stage_foo(train_depth, stage = stage_dat) %>% 
      dplyr::select(depth_var = response, latitude:roms_temp, ind_block)
    test_dum <- stage_foo(test_depth, stage = stage_dat) %>% 
      dplyr::select(depth_var = response, latitude:roms_temp, ind_block)
  }
  
  # subset based on response
  list(train = train_dum, test = test_dum)
}

dum <- map2(model_tbl$response %>% as.character(), 
            model_tbl$stage_dat %>% as.character(), 
            .f = add_data)
model_tbl$train_data <- map(dum, function (x) x$train)
model_tbl$test_data <- map(dum, function (x) x$test)


## add recipe 
# NOTE ind_block needs to be retained to identify subsequent blocking but must
# be removed from recipe and prior to fitting
model_tbl$recipe <- purrr::map(
  model_tbl$train_data,
  function (x) {
    dum <- x %>% 
      dplyr::select(-ind_block)
    
    recipe(depth_var ~ ., data = dum) %>% 
      step_impute_knn(all_predictors(), neighbors = 3) %>%
      step_nzv(all_predictors()) %>% 
      step_dummy(all_predictors(), -all_numeric())
})


## shared settings
train_folds <- groupKFold(model_tbl$train_data[[1]]$ind_block,
                          k = length(unique(model_tbl$train_data[[1]]$ind_block)))
ctrl <-   trainControl(
  method="repeatedcv",
  index = train_folds
)

## define hyperparameters
gbm_grid <-  expand.grid(interaction.depth = c(2, 5, 10),
                         n.trees = c(seq(10, 100, by = 10), 150, 200, 250, 300),
                         shrinkage = 0.1,
                         n.minobsinnode = c(5, 10, 20))

rf_grid <- expand.grid(tune_length = 8,
                       n.trees = seq(100, 500, by = 100))


model = sub_tbl$model_type[[1]]
recipe = sub_tbl$recipe[[1]]
train_data = sub_tbl$train_data[[1]]

fit_foo <- function(model, recipe, train_data) {
  if (model == "gbm") {
    fit <- train(
      recipe,
      train_data %>% dplyr::select(-ind_block),
      method = "gbm", 
      metric = "RMSE",
      maximize = FALSE,
      trControl = ctrl,
      tuneGrid = gbm_grid[1:2, ]
    )
  }
  if (model == "rf") {
    for (i in 1:2) {
      fit <- train(
        recipe,
        train_data %>% dplyr::select(-ind_block),
        method = "ranger", 
        metric = "RMSE",
        maximize = FALSE,
        tuneLength = 2,#unique(rf_grid$tune_length),
        trControl = ctrl,
        num.trees = rf_grid$n.trees[i]
      ) 
    }
  }
  list(results = fit$results,
       top_model = fit$finalModel)
}


tt <- fit_foo(sub_tbl$model_type[[1]],
        sub_tbl$recipe[[1]],
        sub_tbl$train_data[[1]])


sub_tbl <- model_tbl[13:14,]

fit_list <- pmap(list(sub_tbl$model_type,
                      sub_tbl$recipe,
                      sub_tbl$train_data),
                 .f = fit_foo)

depth_gbm <- train(
  model_tbl$recipe[[13]],
  model_tbl$train_data[[13]] %>% dplyr::select(-ind_block),
  method = "gbm", 
  metric = "RMSE",
  maximize = FALSE,
  # tuneLength = 2,
  trControl = ctrl,
  tuneGrid = gbm_grid[1:2, ]
)

depth_rf <- train(
  model_tbl$recipe[[1]],
  model_tbl$train_data[[1]] %>% dplyr::select(-ind_block),
  method = "ranger", 
  metric = "RMSE",
  maximize = FALSE,
  tuneLength = 3,
  trControl = ctrl,
  num.trees = 100
)


min(depth_gbm$results$RMSE)
