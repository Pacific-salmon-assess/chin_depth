### Depth Model Caret -- Model Comparison
## April 8, 2022
## Compete different model structure to identify best predictor
## Multiple dimensions: 
## 1) Model type (GBM vs. RF)
## 2) Hyperparameters
## 3) Response variable (logit transformed, relative depth, absolute depth)


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)
library(future)
library(furrr)


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

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) 


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
  filter(!is.na(roms_temp)) %>%
  dplyr::select(
    depth = pos_depth, rel_depth, logit_rel_depth,
    fl, mean_log_e, stage, utm_x, utm_y, day_night,
    # det_day = local_day,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels() 
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels() 


## PREP MODELS -----------------------------------------------------------------

# helper function for cleaning up stage-specific data
stage_foo <- function(dat, stage_in) {
  dat %>% 
    filter(stage == stage_in) %>% 
    rename(stage_dat = stage)
}


## create dataframe of models to compete
model_tbl <- expand.grid(
  model_type = c("gbm", "rf"),
  response = c("logit_rel_depth", "rel_depth", "depth"),
  stage_dat = c("integrated")
  ) %>% 
  as_tibble() 

# function to add data based on input variables
add_data <- function(response, stage_dat) {
  # select dataset
  if (stage_dat == "integrated") {
    train_dum <- train_depth %>% 
      dplyr::select(depth_var = all_of(response), fl:ind_block)
    test_dum <- test_depth %>% 
      dplyr::select(depth_var = all_of(response), fl:ind_block)
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
      dplyr::select(-ind_block, -max_bathy)
    
    recipe(depth_var ~ ., data = dum) %>% 
      step_dummy(all_predictors(), -all_numeric())
})


## FIT MODELS ------------------------------------------------------------------

## shared settings
train_folds <- groupKFold(
  model_tbl$train_data[[1]]$ind_block,
  k = length(unique(model_tbl$train_data[[1]]$ind_block))
  )
ctrl <-   trainControl(
  method="repeatedcv",
  index = train_folds
)


## define hyperparameters
gbm_grid <-  expand.grid(
  interaction.depth = c(2, 5, 10),
  n.trees = c(seq(10, 100, by = 10), 150, 200, 250, 300, 500),
  shrinkage = c(0.01, 0.1),
  n.minobsinnode = c(5, 10, 20)
)

rf_grid <- expand.grid(tune_length = 10,
                       n.trees = seq(500, 2500, by = 500))


fit_foo <- function(model, recipe, train_data) {
  if (model == "gbm") {
    fit <- train(
      recipe,
      train_data %>% dplyr::select(-ind_block, -max_bathy),
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
    names(fits) <- paste("trees_", rf_grid$n.trees, sep = "")
    for (i in seq_along(fits)) {
      fits[[i]] <- train(
        recipe,
        train_data %>% dplyr::select(-ind_block, -max_bathy),
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
         top_model = fits[[best_trees]]$finalModel)
  }
  return(out)
}


# set up parallel
plan(multisession, workers = 8)


## fit models (separately)
# gbm_tbl <- model_tbl %>% filter(model_type == "gbm")
# gbm_list <- future_pmap(list("gbm",
#                              gbm_tbl$recipe,
#                              gbm_tbl$train_data),
#                         .f = fit_foo,
#                         .options = furrr_options(seed = TRUE))
# names(gbm_list) <- gbm_tbl$response
# saveRDS(gbm_list,
#         here::here("data", "model_fits", "gbm_model_comparison.rds"))
gbm_list <- readRDS(here::here("data", "model_fits", "gbm_model_comparison.rds"))


# rf_tbl <- model_tbl %>% filter(model_type == "rf")
# rf_list <- future_pmap(list("rf",
#                             rf_tbl$recipe,
#                             rf_tbl$train_data),
#                        .f = fit_foo,
#                        .options = furrr_options(seed = TRUE))
# names(rf_list) <- rf_tbl$response
# saveRDS(rf_list,
#         here::here("data", "model_fits", "rf_model_comparison.rds"))
rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))


## COMPARE MODEL STRUCTURES ----------------------------------------------------


## top models for RF cannot be fit to new data so instead extract values from 
# top ranger models then store new ranger objects in tibble
rf_train_list <- model_tbl %>% 
  filter(model_type == "rf") %>% 
  pull(train_data)
top_rangers <- purrr::map2(rf_list, rf_train_list, function (x, y) {
  top_mod <- x$top_model
  baked_dat <- prep(model_tbl$recipe[[1]]) %>%
    bake(.,
         new_data = y)
  
  ranger_rf <- ranger::ranger(
    depth_var ~ .,
    data = baked_dat,
    num.trees = top_mod$param$num.trees,
    mtry = top_mod$tuneValue$mtry
  )
})

# add models to tbl
gbm_tbl <- model_tbl %>% filter(model_type == "gbm")
rf_tbl <- model_tbl %>% filter(model_type == "rf")
rf_tbl$top_model <- top_rangers #map(rf_list, function (x) x$top_model)
gbm_tbl$top_model <- map(gbm_list, function (x) x$top_model) 
model_tbl <- rbind(rf_tbl, gbm_tbl)


# calculate RMSE relative to observations in real space for top models
# bathy vectors to adjust proportional data
train_bathy <- train_depth$max_bathy
test_bathy <- test_depth$max_bathy

real_train <- rf_tbl %>% 
  filter(response == "depth") %>% 
  pull(train_data) %>% 
  as.data.frame() %>% 
  mutate(max_bathy = train_bathy)
real_test <- rf_tbl %>% 
  filter(response == "depth") %>% 
  pull(test_data) %>% 
  as.data.frame() %>% 
  mutate(max_bathy = test_bathy)

rmse_foo <- function(
    mod_in,
    space = c("logit_rel_depth", "rel_depth", "depth"), 
    model_type,
    dat_in
    ) {
  # apply recipe to convert factors to numeric
  baked_dat <- prep(model_tbl$recipe[[1]]) %>%
      bake(.,
           new_data = dat_in %>%
             dplyr::select(-depth_var, -ind_block, -max_bathy))
  
  if (model_type == "rf") {
    preds <- predict(mod_in,
                     data = baked_dat)$predictions
  } else if (model_type == "gbm") {
    preds <- predict(mod_in,
                     newdata = baked_dat)
  }
  if (space == "logit_rel_depth") {
    preds <- plogis(preds)
  }
  if (space %in% c("logit_rel_depth", "rel_depth")) {
    preds <- preds * dat_in$max_bathy
  }
  
  Metrics::rmse(dat_in$depth_var, preds)
}

transformed_rmse_train <- purrr::pmap(
  list(
    model_tbl$top_model, 
    model_tbl$response, 
    model_tbl$model_type
  ),
  rmse_foo,
  dat_in = real_train) %>% 
  unlist()
transformed_rmse_test <- purrr::pmap(
  list(
    model_tbl$top_model, 
    model_tbl$response, 
    model_tbl$model_type
  ),
  rmse_foo,
  dat_in = real_test) %>% 
  unlist()


# summarize RMSE
rmse_out <- model_tbl %>%
  mutate(rmse_train = transformed_rmse_train,
         rmse_test = transformed_rmse_test) %>% 
  select(model_type, response, rmse_train, rmse_test) %>% 
  pivot_longer(cols = c(rmse_train, rmse_test), names_to = "dataset",
               names_prefix = "rmse_",
               values_to = "rmse") %>% 
  mutate(
    dataset = fct_relevel(dataset, "test", after = Inf),
    response = fct_relevel(response, "depth", "rel_depth", "logit_rel_depth")
  )


png(here::here("figs", "model_comp", "rmse_plot.png"), units = "in",
    height = 3.5, width = 6, res = 250)
ggplot(rmse_out) +
  geom_point(aes(x = model_type, y = rmse, fill = response),
             shape = 21,
             position = position_dodge(width=0.75)) +
  facet_wrap(~dataset) +
  ggsidekick::theme_sleek() +
  scale_fill_discrete(name = "Response\nDistribution") +
  labs(y = "Root Mean Square Error", x = "Model Type")
dev.off()

write.csv(rmse_out,
          here::here("figs", "model_comp", "top_model_transformed_rmse.csv"),
          row.names = FALSE)



## HYPER PARAMETER DIAGNOSTICS -------------------------------------------------

# performance tables and figures
gbm_dat <- map2(names(gbm_list), gbm_list, function (name, x) {
  x$results %>% 
    mutate(response = name)
}) %>% 
  bind_rows() 
gbm_hyper_plot <- ggplot(gbm_dat) +
  geom_point(aes(x = n.trees, y = RMSE, color = as.factor(interaction.depth),
                 shape = as.factor(shrinkage))) +
  facet_grid(response~n.minobsinnode, scales = "free_y")

rf_dat <- map2(names(rf_list), rf_list, function (name, x) {
  x$results %>% 
    mutate(response = name)
}) %>% 
  bind_rows() 
rf_hyper_plot <- ggplot(rf_dat) +
  geom_point(aes(x = as.factor(mtry), y = RMSE, color = splitrule,
                 shape = as.factor(min.node.size))) +
  facet_grid(response~n_trees, scales = "free_y")


## export 
pdf(here::here("figs", "model_comp", "hyperpar_tuning.pdf"))
gbm_hyper_plot
rf_hyper_plot
dev.off()
