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


depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(stage = as.factor(stage))


# fits from model comparison (depth best supported response)
fits <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))$depth


# Model performance similar among tree sizes but peaks at intermediate mtry and 
# with extratrees split rule. Use 200 trees best model for subsequent 
# exploration. 
ggplot(fits$results) +
  geom_point(aes(x = as.factor(mtry), y = RMSE, color = splitrule)) +
  facet_grid(as.factor(min.node.size)~n_trees, scales = "free_y")

n_trees_in <- fits$results %>% 
  filter(RMSE == min(RMSE)) %>% 
  pull(n_trees)

mtry_in <- fits$results %>% 
  filter(RMSE == min(RMSE)) %>% 
  pull(mtry)


## REFIT -----------------------------------------------------------------------

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
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))

rf_refit <- quantregForest::quantregForest(
  x = train_depth_baked %>% dplyr::select(-depth),
  y = train_depth_baked$depth,
  data = train_depth_baked,
  mtry = 6,
  ntree = 1000,
  importance = TRUE
)
saveRDS(rf_refit, here::here("data", "model_fits", "depth_quantrf.rds"))
rf_refit <- readRDS(here::here("data", "model_fits", "depth_quantrf.rds"))


# CHECK PREDS ------------------------------------------------------------------

obs_preds <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                     newdata = train_depth_baked, all = TRUE)
colnames(obs_preds) <- c("lo", "mean", "up")

cbind(train_depth_baked, obs_preds) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = depth, y = mean, ymin = lo, ymax = up),
                  alpha = 0.3)


# VARIABLE IMPORTANCE ----------------------------------------------------------

imp_dat <- as.data.frame(rf_refit$importance, row.names = FALSE) %>%
  janitor::clean_names() %>% 
  mutate(var = rownames(rf_refit$importance) %>% 
           fct_reorder(., -percent_inc_mse),
         mse_sd =  rf_refit$importanceSD,
         up = percent_inc_mse + (0.5 * mse_sd),
         lo = percent_inc_mse - (0.5 * mse_sd),
         category = case_when(
           var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y", 
                      "mean_slope") ~ "spatial",
           var %in% c("det_day", "hour") ~ "temporal",
           var %in% c("stage_mature", "fl", "mean_log_e") ~ "biological",
           TRUE ~ "dynamic"
         )) %>% 
  arrange(-percent_inc_mse) 

imp_plot <- ggplot(imp_dat, aes(x = var, y = percent_inc_mse)) +
  geom_point(aes(fill = category), shape = 21, size = 1.5) +
  ggsidekick::theme_sleek()#+
  # scale results in whiskers being not visible
  # geom_pointrange(aes(ymin = lo, ymax = up))

pdf(here::here("figs", "depth_ml", "importance_quantreg.pdf"),
    height = 6, width = 9)
imp_plot
dev.off()


# COUNTERFACTUAL PREDICT -------------------------------------------------------

# generate predictions for maturity stage and different counterfacs (e.g. 
# most important 3 variables)
new_dat <- train_depth_baked %>%
  group_by(stage_mature) %>% 
  summarize(fl = mean(fl),
            mean_log_e = mean(mean_log_e)) %>% 
  mutate(mean_bathy = mean(train_depth_baked$mean_bathy),
         #fl = mean(train_depth_baked$fl),
         utm_x = mean(train_depth_baked$utm_x),
         hour = mean(train_depth_baked$hour),
         det_day = mean(train_depth_baked$det_day),
         utm_y = mean(train_depth_baked$utm_y),
         mean_slope = mean(train_depth_baked$mean_slope),
         shore_dist = mean(train_depth_baked$shore_dist),
         u = mean(train_depth_baked$u),
         v = mean(train_depth_baked$v),
         w = mean(train_depth_baked$w),
         roms_temp = mean(train_depth_baked$roms_temp)#,
         #mean_log_e = mean(train_depth_baked$mean_log_e)
)

# make tibble for different counterfacs
pred_foo <- function(var_in) {
  # necessary to deal with string input
  varname <- ensym(var_in)
  
  # generate stage-specific means for subset of variables
  if (var_in %in% c("det_day", "fl", "mean_log_e")) {
    group_vals <- train_depth_baked %>% 
      group_by(stage_mature) %>% 
      summarize(min_v = min(!!varname),
                max_v = max(!!varname)) 
  } else {
    group_vals <- train_depth_baked %>% 
      summarize(min_v = min(!!varname),
                max_v = max(!!varname))
  }
  
  # change to dataframe
  var_seq <- NULL
  for (i in 1:nrow(group_vals)) {
    var_seq <- c(var_seq, 
                 seq(group_vals$min_v[i], group_vals$max_v[i], 
                     length.out = 100))
  }
  if (var_in %in% c("det_day", "fl", "mean_log_e")) {
    preds_in1 <- data.frame(
      stage_mature = c(rep(0, 100),
                       rep(1, 100)),
      dum = var_seq
    )
  } else {
    preds_in1 <- data.frame(
      stage_mature = c(rep(0, 100),
                       rep(1, 100)),
      dum = rep(var_seq, 2)
    )
  }
  
  preds_in <- preds_in1 %>% 
    rename(!!varname := dum) %>%
    left_join(.,
              # add other variables
              new_dat %>% select(- {{ var_in }}),
              by = "stage_mature")
   
  # make predictions
  preds_out <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                       newdata = preds_in, all = TRUE)
  colnames(preds_out) <- c("lo", "mean", "up")
  cbind(preds_in, preds_out) %>%
    mutate(stage_f = as.factor(stage_mature))
}

counterfac_tbl <- tibble(
  var_in = c("mean_bathy", "hour", "det_day", "shore_dist", "fl", 
             "mean_log_e")) %>% 
  mutate(
    preds = purrr::map(var_in, pred_foo)
  )

plot_list <- purrr::map2(
  counterfac_tbl$var_in, 
  counterfac_tbl$preds, 
  function (var, preds) {
    preds %>%
      ggplot(., aes_string(x = var)) +
      geom_line(aes(y = mean, color = stage_f), size = 1.5) +
      geom_ribbon(aes(ymin = lo, ymax = up, fill = stage_f), alpha = 0.4) +
      ggsidekick::theme_sleek()
  }
)


# similar to above but check for seasonal differences in diurnal cycles
min_hour <- min(train_depth_baked[ , "hour"])
max_hour <- max(train_depth_baked[ , "hour"])

hr_new_dat <- expand_grid(stage_mature = c(0, 1),
                        hour = seq(min_hour, max_hour, length.out = 100),
                        det_day = c(1, 91, 182, 274)) %>% 
  left_join(., 
            # add other variables
            new_dat %>% select(-hour, -det_day),
            by = "stage_mature")

hr_preds <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                     newdata = hr_new_dat, all = TRUE)
colnames(hr_preds) <- c("lo", "mean", "up")
diurnal_seasonal_preds <- cbind(hr_new_dat, hr_preds) %>%
  mutate(stage_f = as.factor(stage_mature)) %>% 
  ggplot(., aes(x = hour)) +
  geom_line(aes(y = mean, color = stage_f), size = 1.5) +
  # geom_ribbon(aes(ymin = lo, ymax = up, fill = stage_f), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  facet_wrap(~as.factor(det_day)) +
  labs(title = "Seasonal Changes in Diurnal Cycle")


pdf(here::here("figs", "depth_ml", "counterfactual_quantreg.pdf"),
    height = 6, width = 8)
plot_list
diurnal_seasonal_preds
dev.off()



# SPATIAL PREDICT --------------------------------------------------------------


coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
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
    w = mean(w, na.rm = T),
    fl = mean(fl, na.rm = T),
    mean_log_e = mean(mean_log_e, na.rm = T)
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
         utm_y_m = utm_y * 1000,
         pred_int_width = pred_up - pred_lo) %>% 
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
  theme(axis.title = element_blank()) +
  labs(title = "Immature Daytime")
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(season == "summer", day == "day"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~stage_mature) +
  theme(axis.text = element_blank()) +
  labs(title = "Summer Daytime")
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(season %in% c("summer"),
                                         day == "day"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~stage_mature) +
  theme(axis.text = element_blank()) +
  labs(title = "Summer Daytime Prediction Interval Width")
ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = pred_dat %>% filter(season %in% c("summer"),
                                    stage_mature == "1"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "depth") +
  ggsidekick::theme_sleek() +
  facet_wrap(~day) +
  theme(axis.text = element_blank()) +
  labs(title = "Mature Summer")
dev.off()