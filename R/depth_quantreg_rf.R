## Fit Quantile Regression Random Forest

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Fit equivalent model with interpolated training data and generate 
#prediction intervals.


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)
library(quantregForest)


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

n_trees_in <- fit_results %>% 
  filter(RMSE == min(RMSE)) %>% 
  pull(n_trees)

mtry_in <- fit_results %>% 
  filter(RMSE == min(RMSE)) %>% 
  pull(mtry)


## REFIT -----------------------------------------------------------------------

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
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))

# rf_refit <- randomForest::randomForest(depth ~ ., 
#                                        data = train_depth_baked, 
#                                        mtry = depth_rf$bestTune$mtry, 
#                                        ntree = n_trees)
rf_refit <- quantregForest::quantregForest(
  x = train_depth_baked %>% dplyr::select(-depth),
  y = train_depth_baked$depth, 
  data = train_depth_baked, 
  mtry = mtry_in, 
  ntree = n_trees_in,
  importance = TRUE
)


# COUNTERFACTUAL PREDICT -------------------------------------------------------

# generate predictions for maturity stage and bathymetric effects
new_dat <- expand.grid(
  stage_mature = c(0, 1),
  mean_bathy = seq(0,
                   round(max(train_depth_baked$mean_bathy), 0),
                   length.out = 100)
) %>%
  mutate(utm_x = mean(train_depth_baked$utm_x),
         hour = mean(train_depth_baked$hour),
         det_day = mean(train_depth_baked$det_day),
         utm_y = mean(train_depth_baked$utm_y),
         mean_slope = mean(train_depth_baked$mean_slope),
         shore_dist = mean(train_depth_baked$shore_dist),
         u = mean(train_depth_baked$u),
         v = mean(train_depth_baked$v),
         w = mean(train_depth_baked$w),
         roms_temp = mean(train_depth_baked$roms_temp))

pred_rf <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                   newdata = new_dat, all = TRUE)
colnames(pred_rf) <- c("lo", "mean", "up")

cbind(new_dat, pred_rf) %>%
  mutate(stage_f = as.factor(stage_mature)) %>% 
  ggplot(., aes(x = mean_bathy)) +
  geom_line(aes(y = mean, color = stage_f), size = 1.5) +
  geom_ribbon(aes(ymin = lo, ymax = up, fill = stage_f), alpha = 0.2) +
  ggsidekick::theme_sleek()


# VARIABLE IMPORTANCE ----------------------------------------------------------

importance(rf_refit, quantiles = c(0.1, 0.5, 0.9))
varImpPlot.qrf(rf_refit, quantiles = c(0.1, 0.5, 0.9), symbols = F, color = T, 
               which.sort = 2)



data(mtcars)
mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000,
                          keep.forest=FALSE, importance=TRUE)
importance(mtcars.rf)
importance(mtcars.rf, type=1)


# SPATIAL PREDICT --------------------------------------------------------------


coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -128, ymin = 45.5, xmax = -122, ymax = 51) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


bath_grid_in <- readRDS(here::here("data", "pred_bathy_grid_utm.RDS")) %>% 
  mutate(id = row_number()) %>% 
  filter(depth < 400)

bath_recipe <- recipe(id ~ ., 
                       data = bath_grid_in) %>% 
  #impute missing ROMS values
  step_impute_knn(all_predictors(), neighbors = 3) %>% 
  prep()
bath_grid <- bake(bath_recipe, new_data = NULL) %>% 
  mutate(utm_x = X / 1000,
         utm_y = Y / 1000) %>% 
  select(utm_x, utm_y, mean_bathy = depth, mean_slope = slope, shore_dist)


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
  stage_mature = c(0, 1)
) %>% 
  mutate(
    day = fct_recode(as.factor(hour), "day" = "0.5", "night" = "12.5"),
    season = fct_recode(as.factor(det_day), "winter" = "1", "spring" = "91",
                        "summer" = "182", "fall" = "274"),
    month = factor(as.factor(det_day), labels = c("1", "4", "7", "10")),
    # create prediction grid based on fixed covariates
    pred_grid = purrr::pmap(
      list(month, hour, det_day, stage_mature),
      function (w, x, y, z) {
        bath_grid %>%
          mutate(
            month = w,
            hour = x,
            det_day = y,
            stage_mature = z) %>%
          left_join(., roms_month_means, by = "month")
      }
    ),
    #generate predictions
    preds = purrr::map(pred_grid, function (x) {
      pred_rf <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                         newdata = x, all = TRUE)
      colnames(pred_rf) <- c("lo", "med", "up")

      bath_grid %>%
        mutate(
          pred_med = pred_rf[, "med"],
          pred_lo = pred_rf[, "lo"],
          pred_up = pred_rf[, "up"]
        )
    }
    )
  ) 

pred_dat <- dat_tbl %>% 
  select(stage_mature, day, season, preds) %>% 
  unnest(cols = preds) %>% 
  mutate(utm_x_m = utm_x * 1000,
         utm_y_m = utm_y * 1000) %>% 
  filter(!mean_bathy > 400)
  

pdf(here::here("figs", "depth_ml", "pred_depth_quantreg.pdf"),
    height = 6, width = 8)
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(stage_mature == "0", day == "day",
                                    season %in% c("winter", "summer")), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~season) +
  theme(axis.title = element_blank())
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(season == "summer", day == "day"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~stage_mature) +
  theme(axis.text = element_blank())
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(season %in% c("summer"),
                                    stage_mature == "1"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~day) +
  theme(axis.text = element_blank())
dev.off()