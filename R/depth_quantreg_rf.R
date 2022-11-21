## Fit Quantile Regression Random Forest

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Fit equivalent model with interpolated training data and generate 
#quantile prediction intervals.
# July 8 (fit model with stock as covariate w/ very minimal improvement in 
# performance; fits in depth_quantrf_stock.rds)


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
  # approximately 4.6k detections have no available ROMS data; exclude for now
  filter(!is.na(roms_temp))


# number of detections
depth_dat_raw %>% 
  group_by(vemco_code) %>% 
  tally() %>% 
  pull(n) %>% 
  range()

# calculate timespan overwhich detections provided
timespan <- depth_dat_raw %>% 
  group_by(vemco_code) %>% 
  summarize(
    min_time = min(date_time_local),
    max_time = max(date_time_local)
  ) %>% 
  mutate(
    timespan = difftime(max_time, min_time, units = "days")
  ) %>% 
  pull(timespan) 
hist(as.numeric(timespan))  


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
    depth = pos_depth, fl, mean_log_e, stage, utm_x, utm_y, day_night,
    # det_day = local_day,
    det_dayx, det_dayy,
    mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block
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


train_folds <- caret::groupKFold(train_depth$ind_block,
                          k = length(unique(train_depth$ind_block)))
depth_ctrl <-   caret::trainControl(
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

dum <- cbind(train_depth_baked, obs_preds) %>% 
  mutate(resid = depth - mean)
ggplot(dum) +
  geom_pointrange(aes(x = depth, y = mean, ymin = lo, ymax = up),
                  alpha = 0.3)

# explore strange bifurcation
dum2 <- dum %>% 
  filter(abs(resid) > 25) %>% 
  glimpse()

ggplot(dum2) +
  geom_pointrange(aes(x = depth, y = mean, ymin = lo, ymax = up,
                      fill = as.factor(stage_mature)), shape = 21,
                  alpha = 0.3)


# VARIABLE IMPORTANCE ----------------------------------------------------------

imp_dat <- rf_refit$importance#ranger::importance(ranger_rf, type = 1, scale = F) 
imp_dat2 <- as.data.frame(imp_dat) %>%
  janitor::clean_names() %>% 
  mutate(
    var = rownames(imp_dat))

imp_dat <- as.data.frame(rf_refit$importance, row.names = FALSE) %>%
  janitor::clean_names() %>% 
  mutate(
    var = rownames(rf_refit$importance) %>% 
      fct_reorder(., -percent_inc_mse),
    # mse_sd =  rf_refit$importanceSD,
    # up = percent_inc_mse + mse_sd,
    # lo = percent_inc_mse - mse_sd,
    category = case_when(
      var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y",
                 "mean_slope") ~ "spatial",
      var %in% c("det_day", "det_dayx", "det_dayy",
                 "day_night_night", "moon_illuminated") ~ "temporal",
      var %in% c("stage_mature", "fl", "mean_log_e") ~ "biological",
      TRUE ~ "dynamic"
    )
  ) %>% 
  arrange(-percent_inc_mse) 
imp_dat$var_f = factor(
  imp_dat$var#, 
  # labels = c("Bottom Depth", "Fork Length", "Maturity", "Year Day", "UTM Y",
  #            "Condition", "UTM X", "Bottom Slope", "Moon Phase", "Zooplankton", 
  #            "Shore Distance", "Temperature", "Oxygen",
  #            "Thermocline Depth", "Day/Night", "H Current 1", "H Current 2", 
  #            "Vertical Current")
)

imp_plot <- ggplot(imp_dat, aes(x = var_f, y = percent_inc_mse)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


png(here::here("figs", "depth_ml", "importance_quantreg.png"),
    height = 4, width = 6, units = "in", res = 250)
imp_plot
dev.off()


# COUNTERFACTUAL PREDICT -------------------------------------------------------

# generate predictions for maturity stage and different counterfacs (e.g. 
# most important 3 variables)
new_dat <- train_depth_baked %>%
  # group_by(stage_mature) %>% 
  # dplyr::summarize(
  #   fl = mean(fl),
  #   mean_log_e = mean(mean_log_e)
  # ) %>% 
  summarize(
    fl = mean(fl),
    mean_log_e = mean(mean_log_e),
    mean_bathy = mean(train_depth_baked$mean_bathy),
    utm_x = mean(train_depth_baked$utm_x),
    local_day = mean(depth_dat_raw$local_day),
    utm_y = mean(train_depth_baked$utm_y),
    mean_slope = mean(train_depth_baked$mean_slope),
    shore_dist = mean(train_depth_baked$shore_dist),
    u = mean(train_depth_baked$u),
    v = mean(train_depth_baked$v),
    w = mean(train_depth_baked$w),
    zoo = mean(train_depth_baked$zoo),
    oxygen = mean(train_depth_baked$oxygen),
    thermo_depth = mean(train_depth_baked$thermo_depth),
    roms_temp = mean(train_depth_baked$roms_temp),
    moon_illuminated = mean(train_depth_baked$moon_illuminated),
    day_night_night = 0.5,
    stage_mature = 0.5
  ) %>% 
  #duplicate 100 times
  as_tibble() %>% 
  slice(rep(1:n(), each = 100))


# make tibble for different counterfacs (commented out sections make 
# stage-specific)
gen_pred_dat <- function(var_in) {
  # necessary to deal with string input
  varname <- ensym(var_in)
  
  # generate stage-specific means for subset of variables
  # if (var_in %in% c("det_day", "fl", "mean_log_e")) {
  #   group_vals <- train_depth_baked %>% 
  #     dplyr::group_by(stage_mature) %>% 
  #     dplyr::summarize(min_v = min(!!varname),
  #                      max_v = max(!!varname)) 
  # } else {
    group_vals <- depth_dat_raw %>% 
      dplyr::summarize(min_v = min(!!varname),
                       max_v = max(!!varname))
  # }
  
  # change to dataframe
  var_seq <- NULL
  for (i in 1:nrow(group_vals)) {
      var_seq <- c(var_seq, 
                   seq(group_vals$min_v[i], group_vals$max_v[i], 
                       length.out = 100))
  }
  
  # if (var_in %in% c("det_day", "fl", "mean_log_e")) {
  #   preds_in1 <- data.frame(
  #     stage_mature = c(rep(0, 100),
  #                      rep(1, 100)),
  #     dum = var_seq
  #   )
  # } else {
    # preds_in1 <- data.frame(
    #   # stage_mature = c(rep(0, 100),
    #   #                  rep(1, 100)),
    #   # dum = rep(var_seq, 2)
    #   dum = var_seq
    # )
  # }
  
  # preds_in <- preds_in1 %>% 
  #   dplyr::rename(!!varname := dum) %>%
  #   left_join(.,
  #             # add other variables
  #             new_dat %>% select(- {{ var_in }}) %>% glimpse(),
  #             by = "stage_mature") 
  new_dat %>% 
    select(- {{ var_in }}) %>% 
    mutate(dum = var_seq) %>% 
    dplyr::rename(!!varname := dum) %>% 
    mutate(
      det_dayx = sin(2 * pi * local_day / 365),
      det_dayy = cos(2 * pi * local_day / 365)
    )
}

pred_foo <- function(preds_in) {
  preds_out <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                       newdata = preds_in, all = TRUE)
  colnames(preds_out) <- c("lo", "mean", "up")
  cbind(preds_in, preds_out)
}


counterfac_tbl <- tibble(
  var_in = c("mean_bathy", "fl", "local_day", "mean_log_e", "thermo_depth", 
             "roms_temp", "zoo", "oxygen", "shore_dist", 
             "moon_illuminated")) %>% 
  mutate(
    pred_dat_in = purrr::map(var_in, gen_pred_dat),
    preds = purrr::map(pred_dat_in, pred_foo)
  )

plot_list <- purrr::map2(
  counterfac_tbl$var_in, 
  counterfac_tbl$preds, 
  function (var, preds) {
    preds %>%
      ggplot(., aes_string(x = var)) +
      geom_line(aes(y = mean#, color = stage_f
                    ), size = 1.5) +
      geom_ribbon(aes(ymin = lo, ymax = up#, fill = stage_f
                      ), alpha = 0.4) +
      ggsidekick::theme_sleek()
  }
)


## categorical predictions for day/night and maturity impacts
new_dat_trim <- new_dat[1:2, ] %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  ) 
mat_plot <- new_dat_trim %>% 
  mutate(
    stage_mature = c(0, 1)
  ) %>% 
  pred_foo(.) %>% 
  ggplot(., aes(x = as.factor(stage_mature))) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
  ggsidekick::theme_sleek() +
  labs(title = "Maturity Counterfactual")
dn_plot <- new_dat_trim %>% 
  mutate(
    day_night_night = c(0, 1)
  ) %>% 
  pred_foo(.) %>% 
  ggplot(., aes(x = as.factor(day_night_night))) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
  ggsidekick::theme_sleek() +
  labs(title = "Day/Night Counterfactual")


pdf(here::here("figs", "depth_ml", "counterfactual_quantreg.pdf"),
    height = 6, width = 8)
plot_list
mat_plot
dn_plot
dev.off()



# SPATIAL PREDICT --------------------------------------------------------------

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


bath_grid_in <- readRDS(here::here("data", "pred_bathy_grid_utm.RDS")) 
bath_grid <- bath_grid_in %>% 
  mutate(utm_x = X / 1000,
         utm_y = Y / 1000) %>% 
  select(utm_x, utm_y, mean_bathy = depth, mean_slope = slope, shore_dist)


# generate grid for predicting residual spatial effects (i.e. all variables)
# at mean value except spatial coords and stage
rand_bath_grid <- bath_grid %>% 
  mutate(mean_bathy = mean(mean_bathy),
         slope = mean(mean_slope),
         shore_dist = mean(shore_dist),
         )


# calculate mean roms_variables for different seasons (using monthly averages)
roms_month_means <- readRDS(here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(month = lubridate::month(date_time_local)) %>% 
  filter(month %in% c(1, 4, 7, 10)) %>% 
  mutate(month = as.factor(as.character(month))) %>% 
  group_by(month) %>% 
  dplyr::summarize(
    u = mean(u, na.rm = T),
    roms_temp = mean(roms_temp, na.rm = T),
    v = mean(v, na.rm = T),
    w = mean(w, na.rm = T),
    zoo = mean(zoo, na.rm = T),
    oxygen = mean(oxygen, na.rm = T),
    thermo_depth = mean(thermo_depth, na.rm = T)) 


## first set are average spatial predictions (i.e. mean biological and temporal 
## attributes)

# biological data
bio_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, mean_log_e, stage) %>% 
  distinct()

# stratify predictions by non-spatial covariates
pred_dat1 <- bath_grid %>% 
  mutate(
    day_night_night = 0.5,
    local_day = 197, #july 1
    moon_illuminated = 0.5,
    stage_mature = 0.5,
    month = "7",
    fl = mean(bio_dat$fl),
    mean_log_e = mean(bio_dat$mean_log_e)
  ) %>% 
  left_join(roms_month_means, by = "month") %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  )

pred_rf1 <- predict(rf_refit, quantiles = c(0.1, 0.5, 0.9),
                   newdata = pred_dat1, all = TRUE)
colnames(pred_rf1) <- c("lo", "med", "up")

pred_dat2 <- pred_dat1 %>%
  mutate(
    pred_med = pred_rf1[, "med"],
    pred_lo = pred_rf1[, "lo"],
    pred_up = pred_rf1[, "up"],
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000,
    pred_int_width = pred_up - pred_lo
  ) %>% 
  filter(mean_bathy < 400)

base_plot <- ggplot() + 
  # geom_sf(data = coast_utm) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) 

mean_depth <- base_plot +
  geom_raster(data = pred_dat2, 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  scale_fill_viridis_c(name = "Mean Depth") +
  theme(legend.position = "top")

var_depth <- base_plot +
  geom_raster(data = pred_dat2, 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  scale_fill_viridis_c(name = "Variation Depth", option = "C")  +
  theme(legend.position = "top")

avg_depth <- cowplot::plot_grid(mean_depth, var_depth, ncol = 2)

png(here::here("figs", "ms_figs", "avg_depth.png"), height = 5, width = 8, 
    units = "in", res = 250)
avg_depth
dev.off()


## generate spatial contrasts
# 1) seasonal effects for immature fish
# 2) maturity effects for avg (July 1)
# 3) moon illumination effects for avg  (July 1)
stage_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, mean_log_e, stage) %>% 
  distinct() %>% 
  group_by(stage) %>% 
  dplyr::summarize(
    fl = mean(fl),
    mean_log_e = mean(mean_log_e)
  ) %>% 
  mutate(stage_mature = ifelse(stage == "mature", 1, 0)) %>% 
  rbind(., 
        data.frame(
          stage = "average",
          fl = mean(depth_dat_raw$fl),
          mean_log_e = mean(depth_dat_raw$mean_log_e),
          stage_mature = 0.5
        ))

pred_tbl <- tibble(
  contrast = rep(c("season", "maturity", "moon light", "dvm"), each = 2),
  local_day = c(16, 197, 197, 197, 197, 197, 197, 197),
  stage_mature = c(0, 0, 0, 1, 0.5, 0.5, 0.5, 0.5),
  moon_illuminated = c(0.5, 0.5, 0.5, 0.5, 0, 1, 0.5, 0.5),
  day_night_night = c(0.5, 0.5, 0.5, 0.5, 1, 1, 0, 1)
) %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365),
    season = fct_recode(as.factor(local_day), 
                        "winter" = "16", "summer" = "197"),
    month = factor(as.factor(local_day), labels = c("1", "7")),
    # create prediction grid based on fixed covariates
    pred_grid = purrr::pmap(
      list(det_dayx, det_dayy, moon_illuminated, month, day_night_night,
           stage_mature),
      function (u, v, w, x, y, z) {
        bath_grid %>%
          mutate(
            det_dayx = u,
            det_dayy = v,
            moon_illuminated = w,
            month = x,
            day_night_night = y,
            stage_mature = z) %>%
          left_join(., stage_dat %>% select(-stage), by = "stage_mature") %>%
          left_join(., roms_month_means, by = "month")
      }
    ),
    # generate predictions
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

pred_dat <- pred_tbl %>% 
  select(contrast, stage_mature, season, moon_illuminated, day_night_night,
         preds) %>% 
  unnest(cols = preds) %>% 
  mutate(utm_x_m = utm_x * 1000,
         utm_y_m = utm_y * 1000,
         pred_int_width = pred_up - pred_lo) %>% 
  filter(!mean_bathy > 400)


# calculate differences for each contrast scenario
season_eff <- pred_dat %>% 
  filter(contrast == "season") %>% 
  mutate(
    mean_depth = mean(pred_med)
  ) %>% 
  select(season, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med,
         mean_depth) %>% 
  pivot_wider(names_from = season, values_from = pred_med) %>%
  mutate(season_diff = (winter - summer) / summer) 
season_map <- base_plot +
  geom_raster(data = season_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = season_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Mean Depth Difference"
    ) +
  labs(title = "Season Effects")


mat_eff <- pred_dat %>% 
  filter(contrast == "maturity") %>% 
  mutate(
    mean_depth = mean(pred_med)
  ) %>% 
  select(stage_mature, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med,
         mean_depth) %>% 
  pivot_wider(names_from = stage_mature, values_from = pred_med) %>%
  mutate(mat_diff = (`0` - `1`) / `1`)
mat_map <- base_plot +
  geom_raster(data = mat_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = mat_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Mean Depth Difference"
  ) +
  labs(title = "Maturity Effects")


moon_eff <- pred_dat %>% 
  filter(contrast == "moon light") %>% 
  mutate(
    mean_depth = mean(pred_med)
  ) %>% 
  select(moon_illuminated, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med,
         mean_depth) %>% 
  pivot_wider(names_from = moon_illuminated, values_from = pred_med) %>%
  mutate(moon_diff = (`0` - `1`) / `1`) 
moonlight_map <- base_plot +
  geom_raster(data = moon_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = moon_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Mean Depth Difference"
  ) +
  labs(title = "Moonlight Effects")


dvm_eff <- pred_dat %>% 
  filter(contrast == "dvm") %>% 
  mutate(
    mean_depth = mean(pred_med)
  ) %>% 
  select(day_night_night, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med,
         mean_depth) %>% 
  pivot_wider(names_from = day_night_night, values_from = pred_med) %>%
  # negative is deeper at night relative to night, pos shallower at night rel to night
  mutate(dvm_diff = (`0` - `1`) / `1`) 
dvm_map <- base_plot +
  geom_raster(data = dvm_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = dvm_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Mean Depth Difference"
  ) +
  labs(title = "DVM Effects")



pdf(here::here("figs", "depth_ml", "spatial_contrasts.pdf"))
season_map
mat_map
moonlight_map
dvm_map
dev.off()


# SPATIAL RESIDUALS ------------------------------------------------------------

dum <- train_depth_baked %>% 
  mutate(
    pred = rf_refit$predicted,
    resid = pred - depth,
    season = case_when(
      det_day >= 1 & det_day < 90 ~ "winter",
      det_day >= 91 & det_day < 182 ~ "spring",
      det_day >= 182 & det_day < 274 ~ "summer",
      det_day >= 274 ~ "fall"
    ),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
    )

base_plot +
  geom_jitter(data = dum, aes(x = utm_x_m, y = utm_y_m, fill = resid),
              shape = 21) +
  scale_fill_viridis_c(name = "Residual")
# no obvious pattern

