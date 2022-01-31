### Depth Models GBM via Caret
## Jan. 25, 2021


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)


# add block IDs (generated in blocking)
# note that these aren't stratified (e.g. individuals are mixed between blocks
# in the time/space blocking structures)
block_list <- readRDS(here::here("data", "5block_ids.RDS"))

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_15min.RDS")) %>% 
  mutate(
    space_block = block_list$space,
    time_block = block_list$time$time_block
  ) %>% 
  left_join(., 
            block_list$individual %>% 
              rename(ind_block = block), 
            by = "vemco_code") 


# check balance of blocks
depth_dat_raw %>% 
  pivot_longer(ends_with("block"), names_to = "block_type") %>%
  ggplot(.) +
  geom_bar(aes(x = as.factor(value)), stat = "count") +
  facet_wrap(~block_type)


# ggplot(depth_dat_raw) +
#   geom_point(aes(x = mean_bathy_c, y = rel_depth))
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = day_c, y = rel_depth), alpha = 0.4) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = hour_c, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_boxplot(aes(x = day_night, y = rel_depth)) +
#   facet_grid(stage~region_f)

# TODO: fit model stratified by tag, stage, space/time(?)

# depth_imm <- depth_dat_raw %>%  
#   filter(stage == "immature") %>%
#   droplevels()

depth_dat <- depth_dat_raw %>% 
  mutate(logit_rel_depth = qlogis(rel_depth)) %>% 
  dplyr::select(
    logit_rel_depth, stage, latitude, longitude, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    ind_block
  ) 

# check no infinite values in transformed relative depth
nrow(depth_dat[is.infinite(depth_dat$logit_rel_depth), ])
depth_dat <- depth_dat[is.finite(depth_dat$logit_rel_depth), ]
# hist(depth_dat$logit_rel_depth)


# subset into training/testing based on 
set.seed(123)
dat_split <- rsample::initial_split(depth_dat %>% dplyr::select(-ind_block),
                                    prop = 0.8, strata = stage)
train_depth <- rsample::training(dat_split)
test_depth <- rsample::testing(dat_split)


# train_depth <- depth_dat %>% 
#   filter(!ind_block == "5") 
# train_folds <- groupKFold(train_depth$ind_block, 
#                           k = length(unique(train_depth$ind_block))) 
# test_depth <- depth_dat %>% 
#   filter(ind_block == "5") 

ggplot(test_depth) +
  geom_point(aes(x = hour, y = logit_rel_depth), alpha = 0.5)

  
## CARET PRE-PROCESSING --------------------------------------------------------

#check correlations with additional bathy variables
corr <- cor(train_depth %>%  dplyr::select(hour:shore_dist, timestamp_n))
ggcorrplot::ggcorrplot(corr)
# shore distance and mean bathy just under 0.75, leave for now


depth_recipe <- recipe(logit_rel_depth ~ ., 
                       data = train_depth) %>% 
  step_nzv(all_predictors()) %>% 
  #consider adding PCA for bathymetric features
  #step_pca(contains("VSA"), prefix = "surf_area_",  threshold = .95) %>% 
  step_dummy(all_predictors(), -all_numeric()) #%>% 
  # step_center(all_predictors()) %>%
  # step_scale(all_predictors())

# check recipe
prep(depth_recipe) %>%
  bake(., new_data = train_depth) %>%
  glimpse()


## FIT CARET -------------------------------------------------------------------

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

# note does not account for different tags among folds
depth_ctrl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
# 4 folds grouped by blocking structure
# depth_ctrl <- trainControl(
#   method="repeatedcv", 
#   index = train_folds
# )

# boosted gradient model 
# adjust grid space for hyperparameter tuning
gbmGrid <-  expand.grid(interaction.depth = c(3, 5, 9),
                        n.trees = (1:25)*50,
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

tictoc::tic()
depth_gbm <- train(depth_recipe, 
                   train_depth, #%>% dplyr::select(-ind_block),
                    method = "gbm", 
                    metric = "RMSE",
                    maximize = FALSE,
                    # tuneLength = 10,
                    trControl = depth_ctrl,
                    tuneGrid = gbmGrid)
tictoc::toc()


# random forest model
tictoc::tic()
depth_rf <- train(depth_recipe, 
                  train_depth, #%>% dplyr::select(-ind_block),
                   method = "ranger", 
                   metric = "RMSE",
                   maximize = FALSE,
                   trControl = depth_ctrl,
                   tuneLength = 6,
                   num.trees = 500)
tictoc::toc()


trellis.par.set(caretTheme())
plot(depth_gbm)  
plot(depth_rf)


## compare
bwplot(resamples(
  list("GBM" = depth_gbm, 
       "RF" = depth_rf)),
       metric = "RMSE")


# predictions
pred_foo <- function(mod, dat = test_depth) {
  preds <- predict(mod, newdata = dat)
  dat$logit_preds <- preds
  
  par(mfrow = c(2, 1))
  plot(logit_preds ~ logit_rel_depth, data = dat)
  abline(0, 1, col = "red")
  plot(plogis(logit_preds) ~ plogis(logit_rel_depth), data = dat)
  abline(0, 1, col = "red")
}

pred_foo(depth_gbm, dat = train_depth)
pred_foo(depth_rf, dat = train_depth)
pred_foo(depth_rf, dat = test_depth)

# evaluate patterns
library(DALEX)
library(DALEXtra)

explainer_gbm <- explain(
  depth_gbm,
  data = dplyr::select(
    train_depth, timestamp_n, hour, det_day, mean_bathy, mean_slope, shore_dist,
    latitude, longitude, stage
  ),
  y = train_depth$logit_rel_depth,
  label = "gbm"
)
explainer_rf <- explain(
  depth_gbm,
  data = dplyr::select(
    train_depth, timestamp_n, hour, det_day, mean_bathy, mean_slope, shore_dist,
    latitude, longitude, stage
  ),
  y = train_depth$logit_rel_depth,
  label = "random forest"
)

# variable importance
plot(feature_importance(explainer_gbm))
plot(feature_importance(explainer_rf))


# function to visualize profiles
make_pdp <- function(param) {
  pdp_gbm <- model_profile(explainer_gbm, N = 400, variables = param)
  pdp_rf <- model_profile(explainer_rf, N = 400, variables = param)
  rbind(as_tibble(pdp_gbm$agr_profiles), as_tibble(pdp_rf$agr_profiles)) %>%
    ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
    geom_line() + xlab(param)
}

# TODO: how to optimize visuals (e.g. smooths?)
# TODO: how to add estimates of uncertainty
make_pdp("hour")
make_pdp("det_day")
make_pdp("mean_bathy")
make_pdp("mean_slope")
make_pdp("shore_dist")
make_pdp("latitude")
make_pdp("longitude")


# generate spatial predictions based on bathymetric data
# TODO: interpolate missing slope estimates (n = 544 of ~20k)

coast <- readRDS(here::here("data", "crop_coast_sf.RDS"))

bath_grid <- readRDS(here::here("data", "pred_bathy_grid.RDS")) %>%
  filter(!is.na(slope)) %>% 
  rename(longitude = X, latitude = Y, mean_bathy = depth, mean_slope = slope)

# stratify predictions by non-spatial covariates
dat <- expand_grid(
  hour = c(0.5, 12.5),
  det_day = c(90, 150, 210, 270)
) %>% 
  mutate(
    # create prediction grid based on fixed covariates
    pred_grid = purrr::map2(hour, det_day, function (x, y) {
      bath_grid %>% 
        mutate(
          hour = x,
          det_day = y
        )
    }),
    # generate predictions
    preds = purrr::map(pred_grid, function (x) {
      dum <- predict(depth_rf, newdata = x)
      # calculate actual depth
      bath_grid %>% 
        mutate(
          logit_pred = dum,
          rel_pred = plogis(logit_pred),
          pred = mean_bathy * rel_pred
        ) %>% 
        filter(
          !(pred > 300)
        )
      }
      )
  )
  


ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat$preds[[1]], 
              aes(x = longitude, y = latitude, fill = pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat$preds[[8]], 
              aes(x = longitude, y = latitude, fill = pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()
