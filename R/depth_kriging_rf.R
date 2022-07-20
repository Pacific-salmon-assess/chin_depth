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
## TODO: interpolation ideally occurs a priori but necessary until ROMS extractions
# updated

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
test_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = test_depth %>% 
         dplyr::select(-ind_block))
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))
saveRDS(train_depth_baked, 
        here::here("data", "baked_training_data.RDS"))
saveRDS(test_depth_baked, 
        here::here("data", "baked_testing_data.RDS"))


# # fit standard RF version of quant-reg model using hyperpars based on exp
# # analyses
# # rf_refit <- readRDS(here::here("data", "model_fits", "depth_quantrf.rds"))
# rf_refit <- randomForest::randomForest(depth ~ ., 
#                                        data = train_depth_baked, 
#                                        mtry = 6, 
#                                        ntree = 1000)
# 
# 
# # generate distance matrix for coordinates
# coord_mat <- train_depth %>% 
#   select(utm_x, utm_y) %>%
#   as.matrix()
# dist_mat <- compute_distance(coord_mat)
# 
# fit <- fit_rfrk(rf_refit, train_depth$depth, dist_mat)



## fit to subset to evaluate relative performance initially
dum <- sample_n(train_depth_baked, 3000, replace = FALSE) 

coord_mat <- dum %>% 
  select(utm_x, utm_y) %>%
  # distinct() %>%
  as.matrix()
dist_mat <- compute_distance(coord_mat)

rf_fit <- randomForest::randomForest(depth ~ ., 
                                     data = dum, 
                                     mtry = 6, 
                                     ntree = 500)

# as above but excludes spatial coordinates
dum_fit <- randomForest::randomForest(depth ~ ., 
                                       data = dum %>% 
                                        select(!starts_with("utm")), 
                                       mtry = 6, 
                                       ntree = 500)

tictoc::tic()
fit <- fit_rfrk(dum_fit, dum$depth, dist_mat)
tictoc::toc()


## COMPARE PERFORMANCE ---------------------------------------------------------

# calculate RMSE of each model with out of sample data
rf_pred <- predict(rf_fit, newdata = test_depth_baked)
rfrk_pred <- pred_rfrk(
  rf = dum_fit,
  rfrk = fit,
  obsv_coord = dum %>% select(utm_x, utm_y),
  pred_coord = test_depth_baked %>% select(utm_x, utm_y),
  Xp = test_depth_baked %>% select(-c(utm_x, utm_y)),
  computeSE = FALSE
)
plot(rf_pred ~ rfrk_pred[, "pred"])


Metrics::rmse(test_depth_baked$depth, rf_pred)
Metrics::rmse(test_depth_baked$depth, rfrk_pred[, "pred"])

# better predictive performance from vanilla random forest than spatial kriging


## COMPARE SPATIAL PREDICTIONS ------------------------------------------------

# Spatial predictions
bath_grid <- readRDS(here::here("data", "pred_bathy_grid_utm.RDS")) %>%
  dplyr::rename(mean_bathy = depth, mean_slope = slope, utm_x_m = X,
                utm_y_m = Y) %>% 
  mutate(utm_x = utm_x_m / 1000,
         utm_y = utm_y_m / 1000)

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

# calculate mean biological data by stage
stage_means <- depth_dat_raw %>% 
  select(vemco_code, stage, fl, mean_log_e) %>% 
  distinct() %>% 
  group_by(stage) %>% 
  summarize(fl = mean(fl, na.rm = TRUE),
            mean_log_e = mean(mean_log_e, na.rm = TRUE))

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
            stage = z,
            stage_mature = ifelse(stage == "mature", 1, 0)) %>%
          left_join(., roms_month_means, by = "month") %>% 
          left_join(., stage_means, by = "stage")
      }
    ),
    # generate predictions
    rf_preds = purrr::map(pred_grid, function (x) {
      bath_grid %>%
        mutate(
          rf_pred = predict(dum_fit, newdata = x),
          rfrk_pred_mat = pred_rfrk(rf = dum_fit, 
                                rfrk = fit, 
                                obsv_coord = dum %>% 
                                  select(utm_x, utm_y),
                                pred_coord = x %>% 
                                  select(utm_x, utm_y),
                                Xp = x, 
                                computeSE = FALSE),
          rfrk_pred = rfrk_pred_mat[ , "pred"]
        )
    }
    )
  ) 


# test kriging model
rf_pred <- predict(rf_fit, newdata = dat_tbl$pred_grid[[1]])
rfrk_pred <- pred_rfrk(
  rf = dum_fit,
  rfrk = fit,
  obsv_coord = dum %>% select(utm_x, utm_y),
  pred_coord = dat_tbl$pred_grid[[1]] %>% select(utm_x, utm_y),
  Xp = dat_tbl$pred_grid[[1]] %>% select(-c(utm_x, utm_y)),
  computeSE = FALSE
)
plot(rf_pred ~ rfrk_pred[, "pred"])

rfrk_preds <- dat_tbl$pred_grid[[1]] %>%
  mutate(pred = tt[ , "pred"])

# plot 
pred_dat <- dat_tbl %>% 
  select(stage, day, season, rf_preds) %>% 
  unnest(cols = rf_preds) #%>%
  # mutate(utm_x_m = utm_x * 1000,
  #        utm_y_m = utm_y * 1000#,
  #        #pred_int_width = pred_up - pred_lo
  #        ) %>% 
  # filter(!mean_bathy > 400)

pred_dat$rfrk_pred2 <- pred_dat$rfrk_pred[ , 1]
  

# coastline
coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))

  
  
  
ggplot() + 
  geom_sf(data = coast_utm) +
  # geom_raster(data = rfrk_preds,
  #             aes(x = utm_x_m, y = utm_y_m, fill = pred)) +
  geom_raster(data = pred_dat %>% filter(stage == "immature",
                                         day == "day",
                                         season %in% c("winter", "summer")),
              aes(x = utm_x_m, y = utm_y_m, fill = pred)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  # facet_wrap(~season) +
  theme(axis.title = element_blank()) +
  labs(title = "Immature Daytime")


# trim to high contrast subset 
imm_summer_preds <- pred_dat %>% 
  filter(stage == "immature",
         day == "day",
         season == "summer")
