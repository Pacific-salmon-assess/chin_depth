### Depth Models via Caret
## Jan. 25, 2021
## Focus on random forest models with depth as response variable and less than 
## 200 trees after comparison in depth_caret_comparison.R ()

library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)


# add block IDs (generated in blocking) (IGNORE AND JUST GENERATE IND FACTORS)
# block_list <- readRDS(here::here("data", "10block_ids.RDS"))

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(logit_rel_depth = qlogis(rel_depth),
         stage = as.factor(stage))


# check no infinite values in transformed relative depth
nrow(depth_dat_raw[is.infinite(depth_dat_raw$logit_rel_depth), ])
depth_dat_raw <- depth_dat_raw[is.finite(depth_dat_raw$logit_rel_depth), ]
# hist(depth_dat$logit_rel_depth)


## basic visualizations
pdf(here::here("figs", "raw_scatter.pdf"), width = 9, height = 5)
ggplot(depth_dat_raw) +
  geom_point(aes(x = mean_bathy, y = rel_depth), alpha = 0.4) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_point(aes(x = det_day, y = rel_depth), alpha = 0.4) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_point(aes(x = hour, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
    geom_point(aes(x = v, y = rel_depth), alpha = 0.5) +
    facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_point(aes(x = u, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_point(aes(x = w, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_point(aes(x = zoo, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw %>% filter(!region_f == "columbia")) +
  geom_point(aes(x = roms_temp, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat_raw) +
  geom_boxplot(aes(x = day_night, y = rel_depth)) +
  facet_grid(stage~region_f)
dev.off()


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
    depth = pos_depth, stage, utm_x, utm_y, 
    hour, det_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, ind_block
  ) 
# imps <- preProcess(depth_dat, method = "knnImpute")
# impute_depth <- predict(imps, depth_dat)

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

# evaluate prevalence of additional vars
length(train_depth$u[!is.na(train_depth$u)])
length(train_depth$v[!is.na(train_depth$v)])
length(train_depth$mean_slope[!is.na(train_depth$mean_slope)])
length(train_depth$w[!is.na(train_depth$w)])
length(train_depth$roms_temp[!is.na(train_depth$roms_temp)])


# subset into training/testing based randomly 
# set.seed(123)
# dat_split <- rsample::initial_split(depth_dat %>% dplyr::select(-ind_block),
#                                     prop = 0.8, strata = stage)
# train_depth <- rsample::training(dat_split)
# test_depth <- rsample::testing(dat_split)

  
## CARET PRE-PROCESSING --------------------------------------------------------

depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  #impute missing ROMS values
  step_impute_knn(all_predictors(), neighbors = 3) %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric()) #%>% 
  # step_center(all_predictors()) %>%
  # step_scale(all_predictors())

# check recipe
imp_train <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block)) %>%
  glimpse()

#check correlations with additional bathy variables
corr <- cor(imp_train %>% dplyr::select(utm_x:roms_temp))
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
# gbm_grid <-  expand.grid(interaction.depth = c(2, 5, 10), #c(3, 5, 9),
#                         n.trees = seq(10, 400, by = 10),
#                         shrinkage = 0.1,
#                         n.minobsinnode = c(5, 10, 20))
# 
# tictoc::tic()
# depth_gbm <- train(
#   depth_recipe,
#   train_depth %>% dplyr::select(-ind_block),
#   method = "gbm", 
#   metric = "RMSE",
#   maximize = FALSE,
#   # tuneLength = 10,
#   trControl = depth_ctrl,
#   tuneGrid = gbm_grid
# )
# tictoc::toc()

fits <- readRDS(here::here("data", "model_fits", "depth_rf_nobin_list.rds"))


# random forest model
tictoc::tic()
# depth_rf <- train(
#   depth_recipe,
#   train_depth %>% dplyr::select(-ind_block),
#   method = "ranger", 
#   metric = "RMSE",
#   maximize = FALSE,
#   tuneLength = 6,
#   trControl = depth_ctrl,
#   num.trees = 500
# )

tree_seq <- seq(50, 200, by = 50)
# fits <- vector(length = length(tree_seq), mode = "list")
# names(fits) <- paste("trees_", tree_seq, sep = "")
for (i in 3:4#seq_along(fits)
     ) {
  fits[[i]] <- train(
    depth_recipe,
    train_depth %>% dplyr::select(-ind_block),
    method = "ranger", 
    metric = "RMSE",
    maximize = FALSE,
    tuneLength = 10,
    trControl = depth_ctrl,
    num.trees = tree_seq[i]
  ) 
  fits[[i]]$results$n_trees <- tree_seq[i]
}
tictoc::toc()


# save models
# saveRDS(depth_gbm, here::here("data", "model_fits", "depth_gbm_15min.rds"))
# saveRDS(depth_rf, here::here("data", "model_fits", "depth_rf_15min.rds"))
saveRDS(fits, here::here("data", "model_fits", "depth_rf_nobin_list.rds"))


fit_results <- purrr::map(fits, function (x) x$results) %>% 
  bind_rows()


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
  trim_dat <- dat #%>% sample_n(1000)
  preds <- predict(mod, newdata = trim_dat)
  
  trim_dat$logit_preds <- preds
  
  par(mfrow = c(2, 1))
  plot(logit_preds ~ logit_rel_depth, data = trim_dat)
  abline(0, 1, col = "red")
  plot(plogis(logit_preds) ~ plogis(logit_rel_depth), data = trim_dat)
  abline(0, 1, col = "red")
  
  # rmse_out <- sqrt(mean((dat$logit_rel_depth - dat$logit_preds)^2))
  # rmse_out_real <- sqrt(mean(
  #   (plogis(dat$logit_rel_depth) - plogis(dat$logit_preds))^2)
  #   )
  # paste("rmse link =", rmse_out, "rmse real =", rmse_out_real)
}


# one option for quantiles using ranger (ask SA about best options)
# rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
# pred <- predict(rf, mtcars[27:32, ], type = "quantiles", 
#                 quantiles = c(0.1, 0.5, 0.9))
# pred$predictions


pred_foo(depth_gbm, dat = train_depth)
pred_foo(depth_rf, dat = train_depth)

pred_foo(depth_gbm, dat = test_depth)
pred_foo(depth_rf, dat = test_depth)

png(here::here("figs", "depth_ml", "predictive_performance_15min_gbm.png"),
    height = 8, width = 4, res = 200, units = "in")
pred_foo(depth_gbm, dat = test_depth)
dev.off()

png(here::here("figs", "depth_ml", "predictive_performance_15min_rf.png"),
    height = 8, width = 4, res = 200, units = "in")
pred_foo(depth_rf, dat = test_depth)
dev.off()


# evaluate patterns
library(DALEX)
library(DALEXtra)

explainer_gbm <- explain(
  depth_gbm,
  data = dplyr::select(
    train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w, 
    roms_temp,
    latitude, longitude, stage
  ),
  y = train_depth$logit_rel_depth,
  label = "gbm"
)
explainer_rf <- explain(
  depth_rf,
  data = dplyr::select(
    train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w, 
    roms_temp,
    latitude, longitude, stage
  ),
  y = train_depth$logit_rel_depth,
  label = "random forest"
)

# check residuals
gbm_hist <- DALEX::model_performance(explainer_gbm)
rf_hist <- DALEX::model_performance(explainer_rf) 

png(here::here("figs", "depth_ml", "hist_resids_15min.png"))
plot(rf_hist, gbm_hist, geom = "histogram")
dev.off()


# variable importance
png(here::here("figs", "depth_ml", "predictor_importance_15min_gbm.png"))
plot(feature_importance(explainer_gbm))
dev.off()
png(here::here("figs", "depth_ml", "predictor_importance_15min_rf.png"))
plot(feature_importance(explainer_rf))
dev.off()


# function to visualize profiles (replaced with standard plots below)
# make_pdp <- function(param) {
#   pdp_gbm <- model_profile(explainer_gbm, N = 400, variables = param)
#   pdp_rf <- model_profile(explainer_rf, N = 400, variables = param)
#   rbind(as_tibble(pdp_gbm$agr_profiles), as_tibble(pdp_rf$agr_profiles)) %>%
#     ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
#     geom_line() + 
#     xlab(param) + ylab("y_hat (logit rel. depth)") +
#     ggsidekick::theme_sleek()
# }

# TODO: how to optimize visuals (e.g. smooths?)
# TODO: how to add estimates of uncertainty
pdp_rf <- model_profile(explainer_rf_real, type = "partial", groups = "stage")
pdp_gbm <- model_profile(explainer_rf_real, type = "partial", groups = "stage")

png(here::here("figs", "depth_ml", "predictor_profiles_15min_rf.png"))
plot(pdp_rf, geom = "profiles")
dev.off()
png(here::here("figs", "depth_ml", "predictor_profiles_15min_gbm.png"))
plot(pdp_gbm, geom = "profiles")
dev.off()


# generate spatial predictions based on bathymetric data
coast <- readRDS(here::here("data", "crop_coast_sf.RDS"))

ggplot() +
  geom_sf(data = coast)

bath_grid <- readRDS(here::here("data", "pred_bathy_grid.RDS")) %>%
  filter(!is.na(slope), 
         depth < 400) %>% 
  dplyr::rename(mean_bathy = depth, mean_slope = slope)

# calculate mean roms_variables for different seasons (using monthly averages)
roms_month_means <- roms_dat %>% 
  filter(month %in% c(1, 4, 7, 10)) %>% 
  mutate(month = as.factor(as.character(month))) %>% 
  group_by(month) %>% 
  dplyr::summarize(
    u = mean(u, na.rm = T),
    roms_temp = mean(roms_temp, na.rm = T),
    v = mean(v, na.rm = T),
    w = mean(w, na.rm = T)
    )

# stratify predictions by non-spatial covariates
dat_tbl <- expand_grid(
  hour = c(0.5, 12.5),
  # jan 1, apr 1, jul 1, oct 1
  det_day = c(1, 91, 182, 274),
  stage = unique(depth_dat$stage)#,
  # replace with seasonal means 
  # u = 0,
  # v = 0,
  # w = 0
) %>% 
  mutate(
    day = fct_recode(as.factor(hour), "day" = "0.5", "night" = "12.5"),
    season = fct_recode(as.factor(det_day), "winter" = "1", "spring" = "91",
                        "summer" = "182", "fall" = "274"),
    month = factor(as.factor(det_day), labels = c("1", "4", "7", "10")),
    # create prediction grid based on fixed covariates
    pred_grid = purrr::pmap(
      list(month, hour, det_day, stage),
      function (w, x, y, z) {
        bath_grid %>%
          mutate(
            month = w,
            hour = x,
            det_day = y,
            stage = z) %>%
          left_join(., roms_month_means, by = "month")
      }
    ),
    # generate predictions
    preds = purrr::map(pred_grid, function (x) {
      dum <- predict(depth_rf, newdata = x)
      # calculate actual depth
      bath_grid %>%
        mutate(
          logit_pred = dum,
          rel_pred = plogis(logit_pred),
          pred = mean_bathy * rel_pred
        )
      }
      )
  ) 
  
dat <- dat_tbl %>% 
  select(stage, day, season, preds) %>% 
  unnest(cols = preds)

pred_depth <- ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(day == "day"), 
              aes(x = longitude, y = latitude, fill = pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  facet_grid(season ~ stage)
pred_rel_depth <- ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(day == "day"), 
              aes(x = longitude, y = latitude, fill = rel_pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  facet_grid(season ~ stage)

pdf(here::here("figs", "depth_ml", "pred_depth_map.pdf"))
pred_depth
pred_rel_depth
dev.off()

png(here::here("figs", "depth_ml", "pred_depth_mat_imm.png"), res = 250, 
    units = "in", height = 6, width = 8)
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(season %in% c("summer")), 
              aes(x = longitude, y = latitude, fill = pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  facet_grid(stage~day)
dev.off()

png(here::here("figs", "depth_ml", "pred_depth_imm_season.png"), res = 250, 
    units = "in", height = 6, width = 8)
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(stage == "immature", day == "day",
                                    season %in% c("winter", "summer")), 
              aes(x = longitude, y = latitude, fill = pred)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  facet_wrap(~ season)
dev.off()


# bathymetry plot
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = bath_grid, 
              aes(x = longitude, y = latitude, fill = mean_bathy)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()
