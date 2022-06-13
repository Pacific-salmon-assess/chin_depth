library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)

## Explore Random Forest Models 
# (script largely duplicated in rf_model_summary.Rmd)

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Fit various hyperparameters, including 50-200 trees in depth_caret.R and 
#save fits to explore here. Models are stored in a list with different tree 
#lengths so best model for each tree is available.

fits <- readRDS(here::here("data", "model_fits", "depth_rf_nobin_list.rds"))

# bind results together
fit_results <- purrr::map(fits, function (x) x$results) %>% 
  bind_rows()

# Model performance similar among tree sizes but peaks at intermediate mtry and 
# with extratrees split rule. Use 200 trees best model for subsequent 
# exploration. 
ggplot(fit_results) +
  geom_point(aes(x = as.factor(mtry), y = RMSE, color = splitrule)) +
  facet_grid(as.factor(min.node.size)~n_trees, scales = "free_y")

n_trees <- fit_results %>% 
  filter(RMSE == min(RMSE)) %>% 
  pull(n_trees)


## REFIT -----------------------------------------------------------------------

# issues using explainers with original fits so refit here

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(stage = as.factor(stage))


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

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  #impute missing ROMS values
  step_impute_knn(all_predictors(), neighbors = 3) %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric())


train_folds <- groupKFold(train_depth$ind_block,
                          k = length(unique(train_depth$ind_block)))
depth_ctrl <-   trainControl(
  method="repeatedcv",
  index = train_folds
)


# fit with random forest to identify top model
depth_rf <- train(
  depth_recipe,
  train_depth %>% dplyr::select(-ind_block),
  method = "ranger",
  metric = "RMSE",
  maximize = FALSE,
  tuneLength = 6,
  trControl = depth_ctrl,
  num.trees = n_trees
)



# refit top model in random forest, using quantregforest to generate uncertainty
# intervals

# apply recipe to dataframe to interpolate values
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))

rf_refit <- randomForest::randomForest(depth ~ ., 
                                       data = train_depth_baked, 
                                       mtry = depth_rf$bestTune$mtry, 
                                       ntree = n_trees)


# new_dat <- expand.grid(
#   stage_mature = c(0, 1),
#   mean_bathy = seq(0,
#                    round(max(train_depth_baked$mean_bathy), 0),
#                    length.out = 100)
# ) %>% 
#   mutate(utm_x = mean(train_depth_baked$utm_x),
#          hour = mean(train_depth_baked$hour),
#          det_day = mean(train_depth_baked$det_day),
#          utm_y = mean(train_depth_baked$utm_y),
#          mean_slope = mean(train_depth_baked$mean_slope),
#          shore_dist = mean(train_depth_baked$shore_dist),
#          u = mean(train_depth_baked$u),
#          v = mean(train_depth_baked$v),
#          w = mean(train_depth_baked$w),
#          roms_temp = mean(train_depth_baked$roms_temp))
# 
# pred_rf <- predict(rf_refit, new_dat, predict_all = TRUE)


#  explainer
explainer_rf <- explain(
  depth_rf,
  data = dplyr::select(
    train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w, 
    roms_temp,
    utm_x, utm_y, stage
  ),
  y = train_depth$depth,
  label = "random forest"
)
var_imp_plot <- plot(feature_importance(explainer_rf))

# partial dependencies
pdp <- model_profile(explainer_rf, type = "partial", groups = "stage")
pdp_plot <- plot(pdp, geom = "profiles")



# save for use in rmd
saveRDS(list(explainer = explainer_rf,
             var_imp_plot = var_imp_plot,
             pdp = pdp),
        here::here("data", "model_fits", "depth_rf_nobin_explainers.rds")
        )

png(here::here("figs", "depth_ml", "predictor_importance_nobin_rf.png"))
# var_imp_plot
exp_list$var_imp_plot
dev.off()

#Predictions with training data look pretty good although the model does 
#chronically underpredict deepest depths and has an unusual bifurcation at 
#deeper values.
rf_mod <- explainer_rf$finalModel
train_dat <- explainer_rf$trainingData 

plot(rf_mod$predictions ~ train_dat$depth)
abline(0, 1, col = "red")


# histogram of residuals looks pretty good.
rf_hist <- DALEX::model_performance(explainer_rf) 
plot(rf_hist, geom = "histogram")

#Stage and bathymetry strongest predictors of depth distribution.
png(here::here("figs", "depth_ml", "predictor_importance_nobin_rf.png"))
plot(feature_importance(explainer_rf))
dev.off()


#Partial effects plots
pdp_rf <- model_profile(explainer_rf, type = "partial", groups = "stage")

png(here::here("figs", "depth_ml", "predictor_profiles_nobin_rf.png"))
plot(pdp_rf, geom = "profiles")
dev.off()


# Spatial predictions
coast <- readRDS(here::here("data", "crop_coast_sf.RDS")) 

bath_grid <- readRDS(here::here("data", "pred_bathy_grid.RDS")) %>%
  filter(!is.na(slope), 
         depth < 400) %>% 
  dplyr::rename(mean_bathy = depth, mean_slope = slope)

# calculate mean roms_variables for different seasons (using monthly averages)
roms_month_means <- readRDS(here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(month = lubridate::month(date_time)) %>% 
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
  stage = unique(depth_dat$stage)
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
      bath_grid %>%
        mutate(
          pred = predict(fits[[3]], newdata = x)
        )
    }
    )
  ) 

dat <- dat_tbl %>% 
  select(stage, day, season, preds) %>% 
  unnest(cols = preds) 


png(here::here("figs", "depth_ml", "pred_depth_mat_imm.png"), res = 250, 
    units = "in", height = 6, width = 8)
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(season %in% c("winter")), 
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


pdf(here::here("figs", "depth_ml", "pred_depth_pres.pdf"),
    height = 6, width = 8)
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(stage == "immature", day == "day",
                                    season %in% c("winter", "summer")), 
              aes(x = longitude, y = latitude, fill = pred)) +
  coord_sf(ylim = c(46.2, 49.2), expand = FALSE) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~season) +
  theme(axis.title = element_blank())
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(season == "summer", day == "day"), 
              aes(x = longitude, y = latitude, fill = pred)) +
  coord_sf(ylim = c(46.2, 49.2), expand = FALSE) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~stage) +
  theme(axis.text = element_blank())
ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = dat %>% filter(season %in% c("summer"),
                                    stage == "mature"), 
              aes(x = longitude, y = latitude, fill = pred)) +
  coord_sf(ylim = c(46.2, 49.2), expand = FALSE) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~day) +
  theme(axis.text = element_blank())
dev.off()