### Depth Model Caret -- Model Comparison with Weights
## Dec 8, 2023
## Compete different model structure to identify best predictor
## Multiple dimensions: 
## 1) Model type (GBM vs. RF)
## 2) Hyperparameters
## 3) Response variable (logit transformed, relative depth, absolute depth)
## Weights variable by inverse of a) n_det or b) sqrt(n_det)

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
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  filter(!is.na(roms_temp)) %>%
  # sample_n(., size = 5000) %>% 
  group_by(vemco_code) %>% 
  #weight based on number of observations
  mutate(
    n_dets = n(),
    wt = 1 / sqrt(n_dets),
    wt2 = 1 / n_dets
  ) %>% 
  ungroup() %>% 
  filter(
    !grepl("2022", vemco_code)
  ) 


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
    depth = pos_depth, rel_depth, logit_rel_depth,
    fl, lipid, med_stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block, wt, wt2
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels() 
train_weights <- train_depth$wt
train_weights2 <- train_depth$wt2
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels() 


## PREP MODELS -----------------------------------------------------------------

## create dataframe of models to compete
model_tbl <- expand.grid(
  model_type = c("gbm", "rf", "rf_weighted", "rf_weighted2"),
  response = c("logit_rel_depth", "rel_depth", "depth")
  ) %>% 
  as_tibble() 

# function to add data based on input variables
add_data <- function(response) {
  train_dum <- train_depth %>% 
    dplyr::select(depth_var = all_of(response), fl:ind_block)
  test_dum <- test_depth %>% 
    dplyr::select(depth_var = all_of(response), fl:ind_block)
  # subset based on response
  list(train = train_dum, test = test_dum)
}

dum <- map(model_tbl$response %>% as.character(), 
           .f = add_data)
model_tbl$train_data <- map(dum, function (x) x$train)
model_tbl$test_data <- map(dum, function (x) x$test)


## add recipe and bake data
# NOTE ind_block needs to be retained to identify subsequent blocking but must
# be removed from recipe and prior to fitting
model_tbl$baked_train <- purrr::map(
  model_tbl$train_data,
  function (x) {
    dum <- x %>% 
      dplyr::select(-ind_block, -max_bathy)
    
    depth_recipe <- recipe(depth_var ~ ., data = dum) %>% 
      step_dummy(all_predictors(), -all_numeric())
    
    # bake outside of recipe due to issues with weights
    prep(depth_recipe) %>%
      bake(., 
           new_data = dum)
})


## FIT MODELS ------------------------------------------------------------------

## shared settings
train_folds <- groupKFold(
  model_tbl$train_data[[1]]$ind_block,
  k = length(unique(model_tbl$train_data[[1]]$ind_block))
  )
ctrl <-   trainControl(
  method = "repeatedcv",
  index = train_folds
)


## define hyperparameters
gbm_grid <-  expand.grid(
  interaction.depth = c(2, 5, 10),
  n.trees = c(seq(10, 100, by = 10), 150, 200, 250, 300, 500),
  shrinkage = c(0.01, 0.1),
  n.minobsinnode = c(5, 10, 20)
)

rf_grid <- expand.grid(mtry = seq(1, 17, by = 2),
                       min.node.size = 5,
                       splitrule = c("variance", "extratrees"))
rf_n_trees <- seq(1000, 3000, by = 500)



fit_foo <- function(model, baked_train) {
  if (model == "gbm") {
    fit <- train(
      depth_var ~ .,
      baked_train,
      method = "gbm", 
      metric = "RMSE",
      maximize = FALSE,
      trControl = ctrl,
      tuneGrid = gbm_grid
    )
    out <- list(results = fit$results,
                top_model = fit$finalModel)
  }
  
  if (model %in% c("rf_weighted", "rf", "rf_weighted2")) {
    # iterate over different number of trees
    fits <- vector(length = length(rf_n_trees), mode = "list")
    names(fits) <- paste("trees_", rf_n_trees, sep = "")
    wts <- if (model == "rf_weighted") { 
      train_weights
    } else if (model == "rf_weighted2") {
        train_weights2
    } else {
        NULL
      }
    for (i in 1:length(rf_n_trees)) {
      fits[[i]] <- train(
        depth_var ~ .,
        baked_train,
        method = "ranger", 
        num.threads = 4,
        weights = wts,
        metric = "RMSE",
        maximize = FALSE,
        # tuneLength = unique(rf_grid$tune_length),
        tuneGrid = rf_grid,
        trControl = ctrl,
        num.trees = rf_n_trees[i]
      ) 
      fits[[i]]$results$n_trees <- rf_n_trees[i]
    }
    fit_results <- purrr::map(fits, function (x) x$results) %>% 
      bind_rows()
    best_trees <- fit_results %>% 
      filter(RMSE == min(fit_results$RMSE)) %>% 
      pull(n_trees) %>% 
      paste("trees_", ., sep = "")
    
    out <- list(
      results = fit_results,
      top_model = fits[[best_trees]]$finalModel
    )
  }
  return(out)
}


# set up parallel
plan(multisession, workers = ncores)


# fit models (separately)
gbm_tbl <- model_tbl %>% filter(model_type == "gbm")
gbm_list <- future_pmap(list("gbm",
                             # gbm_tbl$recipe,
                             gbm_tbl$baked_train),
                        .f = fit_foo,
                        .options = furrr_options(seed = TRUE))
names(gbm_list) <- gbm_tbl$response
saveRDS(gbm_list,
        here::here("data", "model_fits", "gbm_model_comparison.rds"))
gbm_list <- readRDS(
  here::here("data", "model_fits", "gbm_model_comparison.rds"))


rf_tbl <- model_tbl %>% filter(model_type == "rf")
rf_list <- future_pmap(list("rf",
                            rf_tbl$baked_train),
                       .f = fit_foo,
                       .options = furrr_options(seed = TRUE))
names(rf_list) <- rf_tbl$response
saveRDS(rf_list,
        here::here("data", "model_fits", "rf_model_comparison.rds"))
rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))


rf_weighted_tbl <- model_tbl %>% filter(model_type == "rf_weighted")
rf_weighted_list <- future_pmap(
  list("rf_weighted",
       # rf_weighted_tbl$recipe,
       rf_weighted_tbl$baked_train),
  .f = fit_foo,
  .options = furrr_options(seed = TRUE)
)
names(rf_weighted_list) <- rf_weighted_tbl$response
saveRDS(rf_weighted_list,
        here::here("data", "model_fits", "rf_model_comparison_weighted.rds"))
rf_weighted_list <- readRDS(
  here::here("data", "model_fits", "rf_model_comparison_weighted.rds")
  )

rf_weighted2_tbl <- model_tbl %>% filter(model_type == "rf_weighted2")
rf_weighted2_list <- future_pmap(
  list("rf_weighted2",
       rf_weighted2_tbl$baked_train),
  .f = fit_foo,
  .options = furrr_options(seed = TRUE)
)
names(rf_weighted2_list) <- rf_weighted2_tbl$response
saveRDS(rf_weighted2_list,
        here::here("data", "model_fits", "rf_model_comparison_weighted2.rds"))
rf_weighted2_list <- readRDS(
  here::here("data", "model_fits", "rf_model_comparison_weighted2.rds")
)


## COMPARE MODEL STRUCTURES ----------------------------------------------------


## top models for RF cannot be fit to new data so instead extract values from 
# top ranger models then store new ranger objects in tibble
rf_train_list <- model_tbl %>% 
  filter(model_type == "rf") %>% 
  pull(baked_train)
top_rangers <- purrr::map2(
  rf_list, rf_train_list, function (x, y) {
    top_mod <- x$top_model
    ranger::ranger(
      depth_var ~ .,
      data = y,
      num.trees = top_mod$param$num.trees,
      mtry = top_mod$tuneValue$mtry,
      splitrule = top_mod$splitrule
    )
  }
  )

rf_weighted_train_list <- model_tbl %>% 
  filter(model_type == "rf_weighted") %>% 
  pull(baked_train)
top_rangers_weighted <- purrr::map2(
  rf_weighted_list, rf_weighted_train_list, function (x, y) {
    top_mod <- x$top_model
    
    ranger::ranger(
      depth_var ~ .,
      data = y,
      case.weights = train_weights,
      num.trees = top_mod$param$num.trees,
      mtry = top_mod$tuneValue$mtry,
      splitrule = top_mod$splitrule
    )
    }
  )

rf_weighted2_train_list <- model_tbl %>% 
  filter(model_type == "rf_weighted2") %>% 
  pull(baked_train)
top_rangers_weighted2 <- purrr::map2(
  rf_weighted2_list, rf_weighted2_train_list, function (x, y) {
    top_mod <- x$top_model
    
    ranger::ranger(
      depth_var ~ .,
      data = y,
      case.weights = train_weights2,
      num.trees = top_mod$param$num.trees,
      mtry = top_mod$tuneValue$mtry,
      splitrule = top_mod$splitrule
    )
  }
)


# add models to tbl
gbm_tbl <- model_tbl %>% filter(model_type == "gbm")
rf_tbl <- model_tbl %>% filter(model_type == "rf")
rf_weighted_tbl <- model_tbl %>% filter(model_type == "rf_weighted")
rf_weighted2_tbl <- model_tbl %>% filter(model_type == "rf_weighted2")
rf_tbl$top_model <- top_rangers
rf_weighted_tbl$top_model <- top_rangers_weighted
rf_weighted2_tbl$top_model <- top_rangers_weighted2
gbm_tbl$top_model <- map(gbm_list, function (x) x$top_model) 
model_tbl <- list(gbm_tbl, rf_tbl, rf_weighted_tbl, rf_weighted2_tbl) %>% 
  bind_rows()
saveRDS(model_tbl,
        here::here("data", "model_fits", "comp_top_model_tbl.rds"))
model_tbl <- readRDS(here::here("data", "model_fits", "comp_top_model_tbl.rds"))


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

rmse_foo <- function (
    mod_in, space = c("logit_rel_depth", "rel_depth", "depth"), 
    model_type, dat_in
    ) {
  # apply recipe to convert factors to numeric
  dum <- dat_in %>% 
    dplyr::select(-ind_block, -max_bathy)
  
  depth_recipe <- recipe(depth_var ~ ., data = dum) %>% 
    step_dummy(all_predictors(), -all_numeric())
  
  # bake outside of recipe due to issues with weights
  baked_dat <- prep(depth_recipe) %>%
    bake(., 
         new_data = dum)
  
  if (grepl("rf", model_type)) {
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
    response = fct_relevel(response, "depth", "rel_depth"#, "logit_rel_depth"
                           )
  )


png(here::here("figs", "model_comp", "rmse_plot1.png"), units = "in",
    height = 3.5, width = 6, res = 250)
ggplot(rmse_out %>% filter(model_type %in% c("gbm", "rf_weighted"))) +
  geom_point(aes(x = model_type, y = rmse, fill = response),
             shape = 21,
             position = position_dodge(width=0.75)) +
  facet_wrap(~dataset) +
  ggsidekick::theme_sleek() +
  scale_fill_discrete(name = "Response\nDistribution") +
  labs(y = "Root Mean Square Error", x = "Model Type")
dev.off()

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


write.csv(rmse_out %>% arrange(dataset, desc(rmse)),
          here::here("figs", "model_comp", "top_model_transformed_rmse.csv"),
          row.names = FALSE)



## HYPER PARAMETER DIAGNOSTICS -------------------------------------------------

# performance tables and figures
gbm_dat <- map2(names(gbm_list), gbm_list, function (name, x) {
  x$results %>% 
    mutate(response = name,
           interaction.depth = as.factor(interaction.depth),
           shrinkage = as.factor(shrinkage))
}) %>% 
  bind_rows() 
gbm_hyper_plot <- ggplot(gbm_dat) +
  geom_line(aes(x = n.trees, y = RMSE, color = interaction.depth,
                 lty = shrinkage)) +
  facet_grid(response~n.minobsinnode, scales = "free_y") +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()

# weighted and unweighted
rf_dat <- map2(names(rf_list), rf_list, function (name, x) {
  x$results %>% 
    mutate(response = name,
           min.node.size = as.factor(min.node.size),
           model = "rf")
}) %>% 
  bind_rows() 
rf_weighted_dat <- map2(names(rf_weighted_list), rf_weighted_list, function (name, x) {
  x$results %>% 
    mutate(response = name,
           min.node.size = as.factor(min.node.size),
           model = "rf_weighted")
}) %>% 
  bind_rows()
rf_weighted2_dat <- map2(names(rf_weighted2_list), rf_weighted2_list, function (name, x) {
  x$results %>% 
    mutate(response = name,
           min.node.size = as.factor(min.node.size),
           model = "rf_weighted2")
}) %>% 
  bind_rows() 


rf_hyper_plot <- ggplot(rf_dat) +
  geom_line(aes(x = mtry, y = RMSE, color = splitrule)) +
  facet_grid(response ~ n_trees, scales = "free_y") +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()

rf_hyper_plot_comp <- list(rf_dat, rf_weighted_dat, rf_weighted2_dat) %>% 
  bind_rows() %>% 
  filter(n_trees == "1000" & splitrule == "extratrees") %>% 
  ggplot(.) +
  geom_line(aes(x = mtry, y = RMSE, color = model)) +
  facet_wrap(~response, scales = "free_y") +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()

## export 
png(here::here("figs", "model_comp", "gbm_hyper_plot.png"), units = "in",
    height = 4.5, width = 7.5, res = 250)
gbm_hyper_plot
dev.off()

png(here::here("figs", "model_comp", "rf_hyper_plot_comp.png"), units = "in",
    height = 4.5, width = 7.5, res = 250)
rf_hyper_plot_comp
dev.off()

png(here::here("figs", "model_comp", "rf_hyper_plot.png"), units = "in",
    height = 4.5, width = 7.5, res = 250)
rf_hyper_plot
dev.off()
