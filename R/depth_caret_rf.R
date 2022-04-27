library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)

## Explore Random Forest Models

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Fit various hyperparameters, including 50-200 trees in depth_caret and 
#save fits to explore here. Models are stored in a list with different tree 
#lengths so best model for each tree is available.

fits <- readRDS(here::here("data", "model_fits", "depth_rf_nobin_list.rds"))

# bind results together
fit_results <- purrr::map(fits, function (x) x$results) %>% 
  bind_rows()

# Model performance similar among tree sizes but peaks at intermediate mtry and 
# with extratrees split rule. Use 100 trees best model for subsequent 
# exploration. 
ggplot(fit_results) +
  geom_point(aes(x = as.factor(mtry), y = RMSE, color = splitrule)) +
  facet_grid(as.factor(min.node.size)~n_trees, scales = "free_y")

#Predictions with training data look pretty good although the model does 
#chronically underpredict deepest depths and has an unusual bifurcation at 
#deeper values.
rf_mod <- fits[[4]]$finalModel
train_dat <- fits[[4]]$trainingData 

plot(rf_mod$predictions ~ train_dat$depth)
abline(0, 1, col = "red")

explainer_rf <- explain(
  fits[[4]],
  data = dplyr::select(
    train_dat, 
    hour, det_day, mean_bathy, mean_slope, shore_dist, u, v, w, 
    roms_temp, utm_x, utm_y, stage
  ),
  y = train_dat$depth,
  label = "random forest"
)

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
          pred = predict(fits[[4]], newdata = x)
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