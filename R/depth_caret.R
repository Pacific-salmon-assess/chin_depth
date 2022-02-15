### Depth Models GBM via Caret
## Jan. 25, 2021


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)


# add block IDs (generated in blocking) (IGNORE AND JUST GENERATE IND FACTORS)
# block_list <- readRDS(here::here("data", "10block_ids.RDS"))

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_15min.RDS")) %>% 
  mutate(logit_rel_depth = qlogis(rel_depth),
         stage = as.factor(stage))


# check no infinite values in transformed relative depth
nrow(depth_dat_raw[is.infinite(depth_dat_raw$logit_rel_depth), ])
depth_dat_raw <- depth_dat_raw[is.finite(depth_dat_raw$logit_rel_depth), ]
# hist(depth_dat$logit_rel_depth)


## basic visualizations
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = mean_bathy_c, y = rel_depth))
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = day_c, y = rel_depth), alpha = 0.4) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = hour_c, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#     geom_point(aes(x = v, y = rel_depth), alpha = 0.5) +
#     facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = u, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = w, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = zoo, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_boxplot(aes(x = day_night, y = rel_depth)) +
#   facet_grid(stage~region_f)



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
    logit_rel_depth, stage, latitude, longitude, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, ind_block
  ) 
# imps <- preProcess(depth_dat, method = "knnImpute")
# impute_depth <- predict(imps, depth_dat)

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()


length(train_depth$u[!is.na(train_depth$u)])
length(train_depth$v[!is.na(train_depth$v)])
length(train_depth$zoo[!is.na(train_depth$zoo)])
length(train_depth$w[!is.na(train_depth$w)])


# subset into training/testing based randomly 
# set.seed(123)
# dat_split <- rsample::initial_split(depth_dat %>% dplyr::select(-ind_block),
#                                     prop = 0.8, strata = stage)
# train_depth <- rsample::training(dat_split)
# test_depth <- rsample::testing(dat_split)

  
## CARET PRE-PROCESSING --------------------------------------------------------


depth_recipe <- recipe(logit_rel_depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  step_impute_knn(all_predictors(), neighbors = 3) %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric()) #%>% 
  #impute missing ROMS values
  # step_center(all_predictors()) %>%
  # step_scale(all_predictors())

# check recipe
imp_train <- prep(depth_recipe) %>%
  bake(., new_data = train_depth %>% 
         dplyr::select(-ind_block)) %>%
  glimpse()

#check correlations with additional bathy variables
corr <- cor(imp_train %>%  dplyr::select(latitude:w))
ggcorrplot::ggcorrplot(corr)
# looks pretty good


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

# control if using random blocks
# depth_ctrl <- trainControl(## 10-fold CV
#   method = "repeatedcv",
#   number = 10,
#   ## repeated ten times
#   repeats = 10)

# identify blocks in training data
train_folds <- groupKFold(train_depth$ind_block,
                          k = length(unique(train_depth$ind_block)))
depth_ctrl <-   trainControl(
  method="repeatedcv",
  index = train_folds
)


# boosted gradient model 
# adjust grid space for hyperparameter tuning
gbm_grid <-  expand.grid(interaction.depth = c(2, 5, 10), #c(3, 5, 9),
                        n.trees = seq(10, 400, by = 10),
                        shrinkage = 0.1,
                        n.minobsinnode = c(5, 10, 20))

tictoc::tic()
depth_gbm <- train(
  depth_recipe,
  train_depth %>% dplyr::select(-ind_block),
  method = "gbm", 
  metric = "RMSE",
  maximize = FALSE,
  # tuneLength = 10,
  trControl = depth_ctrl,
  tuneGrid = gbm_grid
)
tictoc::toc()


# random forest model
tictoc::tic()
depth_rf <- train(
  depth_recipe,
  train_depth %>% dplyr::select(-ind_block),
  method = "ranger", 
  metric = "RMSE",
  maximize = FALSE,
  tuneLength = 6,
  trControl = depth_ctrl,
  num.trees = 500
)
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

pred_foo(depth_gbm, dat = test_depth)
pred_foo(depth_rf, dat = test_depth)


# evaluate patterns
library(DALEX)
library(DALEXtra)

explainer_gbm <- explain(
  depth_gbm,
  data = dplyr::select(
    train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w,
    latitude, longitude, stage
  ),
  y = train_depth$logit_rel_depth,
  label = "gbm"
)
explainer_rf <- explain(
  depth_gbm,
  data = dplyr::select(
    train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w,
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
make_pdp("u")
make_pdp("v")
make_pdp("w")
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
