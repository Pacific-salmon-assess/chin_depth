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
    depth = pos_depth, rel_depth, logit_rel_depth, utm_y, utm_x, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, stage, ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels() 
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels() 


## MODELS ----------------------------------------------------------------------

# helper function for cleaning up stage-specific data
stage_foo <- function(dat, stage_in) {
  dat %>% 
    filter(stage == stage_in) %>% 
    rename(stage_dat = stage)
}


## create dataframe of models to compete
## NOTE exclude stage-specific models for now because issues with convergence 
# with binned data
model_tbl <- expand.grid(
  model_type = c("gbm", "rf"),
  response = c("logit_rel_depth", "rel_depth", "depth"),
  stage_dat = c("integrated")#, "mature", "immature")
  ) %>% 
  as_tibble() 

# function to add data based on input variables
add_data <- function(response, stage_dat) {
  # select dataset
  if (stage_dat == "integrated") {
    train_dum <- train_depth %>% 
      dplyr::select(depth_var = response, utm_y:stage, ind_block)
    test_dum <- test_depth %>% 
      dplyr::select(depth_var = response, utm_y:stage, ind_block)
  } 
  if (stage_dat %in% c("mature", "immature")) {
    train_dum <- stage_foo(train_depth, stage = stage_dat) %>% 
      dplyr::select(depth_var = response, utm_y:roms_temp, ind_block)
    test_dum <- stage_foo(test_depth, stage = stage_dat) %>% 
      dplyr::select(depth_var = response, utm_y:roms_temp, ind_block)
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


imp_train <- prep(model_tbl$recipe[[5]]) %>%
  bake(.,
       new_data = model_tbl$train_data[[5]] %>%
         dplyr::select(-ind_block)) %>%
  glimpse()


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
                         shrinkage = c(0.01, 0.1),
                         n.minobsinnode = c(5, 10, 20))

rf_grid <- expand.grid(tune_length = 6,
                       n.trees = seq(100, 500, by = 100))

sub_tbl <- model_tbl[c(5, 6),]

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
      tuneGrid = gbm_grid
    )
    out <- list(results = fit$results,
                top_model = fit$finalModel)
  }
  if (model == "rf") {
    # iterate over different number of trees
    fits <- vector(length = length(rf_grid$n.trees), mode = "list")
    names(fit_list) <- paste("trees_", rf_grid$n.trees, sep = "")
    for (i in seq_along(fit_list)) {
      fits[[i]] <- train(
        recipe,
        train_data %>% dplyr::select(-ind_block),
        method = "ranger", 
        metric = "RMSE",
        maximize = FALSE,
        tuneLength = unique(rf_grid$tune_length),
        trControl = ctrl,
        num.trees = rf_grid$n.trees[i]
      ) 
      fits[[i]]$results$n_trees <- rf_grid$n.trees[i]
    }
    fit_results <- purrr::map(fits, function (x) x$results) %>% 
      bind_rows()
    best_trees <- fit_results %>% 
      filter(RMSE == min(fit_results$RMSE)) %>% 
      pull(n_trees) %>% 
      paste("trees_", ., sep = "")
    
    out <- list(results = fit_results,
         top_model = fit_list[[best_trees]]$finalModel)
  }
  return(out)
}

tt <- fit_foo(model = model_tbl$model_type[[1]],
              recipe = model_tbl$recipe[[1]],
              train_data = model_tbl$train_data[[1]])


## fit models (separately)\
gbm_tbl <- model_tbl %>% filter(model_type == "gbm")
gbm_list <- pmap(list("gbm",
                      gbm_tbl$recipe,
                      gbm_tbl$train_data),
                 .f = fit_foo)
names(gbm_list) <- gbm_tbl$response
saveRDS(gbm_list, here::here("data", "gbm_model_comparison.rds"))


# performance table
gbm_dat <- map2(names(gbm_list), gbm_list, function (name, x) {
  x$results %>% 
    mutate(response = name)
}) %>% 
  bind_rows() 
ggplot(gbm_dat) +
  geom_point(aes(x = n.trees, y = RMSE, color = as.factor(interaction.depth),
                 shape = as.factor(shrinkage))) +
  facet_grid(response~n.minobsinnode, scales = "free_y")



purrr::map2(gbm_list, gbm_tbl$response, function (x) {
  x$results %>% 
    mutate(resp = y)
}) %>% 
  bind_rows() %>% 
  glimpse()

rf_tbl <- model_tbl %>% filter(model_type == "rf")
rf_list <- pmap(list("rf",
                      rf_tbl$recipe,
                      rf_tbl$train_data),
                 .f = fit_foo)


# predictions
preds <- predict(tt$top_model, 
                 newdata = imp_train %>%
                   select(-depth_var))

plot(preds ~ imp_train$depth_var, alpha = 0.4)


##