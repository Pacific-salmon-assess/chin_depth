## Fit Quantile Regression Random Forest

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (1500) and fit to untransformed depth 
#data. Fit equivalent model with interpolated training data and generate 
#quantile prediction intervals.
# July 8 (fit model with stock as covariate w/ very minimal improvement in 
# performance; fits in depth_quantrf_stock.rds)
# Nov 22 fit using ranger rather than quantregForest to be consistent with
# model comparison
# Nov 23 switch to transformed depth based on updated model comparison


library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)


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
    depth = rel_depth, fl, mean_log_e, stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block, -max_bathy)) %>% 
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
         dplyr::select(-ind_block, -max_bathy))

# pull model attributes from top ranger
rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))
top_mod <- rf_list[[2]]$top_model

ranger_rf <- ranger::ranger(
  depth ~ .,
  data = train_depth_baked,
  num.trees = top_mod$param$num.trees,
  mtry = top_mod$tuneValue$mtry,
  keep.inbag = TRUE,
  quantreg = TRUE,
  importance = "permutation"
)

saveRDS(ranger_rf,
        here::here("data", "model_fits", "relative_rf_ranger.rds"))
ranger_rf <- readRDS(here::here("data", "model_fits", "relative_rf_ranger.rds"))


# CHECK PREDS ------------------------------------------------------------------

obs_preds <- predict(ranger_rf,
                     data = train_depth_baked)

dum <- train_depth %>% 
  mutate(mean_pred = obs_preds$predictions,
         mean_pred_real = mean_pred * max_bathy,
         depth_real = depth * max_bathy,
         split_group = "train")
plot(depth ~ mean_pred, dum)
plot(depth_real ~ mean_pred_real, dum)


# all preds 
test_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = test_depth %>% 
         dplyr::select(-ind_block, -max_bathy))
test_preds <- predict(ranger_rf,
                      data = test_depth_baked)
dum_test <- test_depth %>% 
  mutate(mean_pred = test_preds$predictions,
         mean_pred_real = mean_pred * max_bathy,
         depth_real = depth * max_bathy,
         split_group = "test")

all_preds <- rbind(dum, dum_test) 

fit_obs <- ggplot() +
  geom_point(
    data = all_preds,
    aes(x = depth_real, y = mean_pred_real, fill = split_group),
    shape = 21, alpha = 0.2
    ) +
  labs(
    x = "Observed Depth", y = "Predicted Mean Depth"
  ) +
  scale_fill_discrete(name = "") +
  ggsidekick::theme_sleek() 

png(here::here("figs", "ms_figs_rel", "obs_preds_rel.png"),
    height = 5.5, width = 6, units = "in", res = 200)
fit_obs
dev.off()

# rmse of each group
Metrics::rmse(dum$depth, dum$mean_pred)
Metrics::rmse(dum_test$depth, dum_test$mean_pred)


# VARIABLE IMPORTANCE ----------------------------------------------------------

imp_vals <- ranger::importance(ranger_rf, type = "permutation", scale = F) 
imp_dat <- data.frame(
  var = names(imp_vals),
  val = imp_vals
  ) %>% 
  mutate(
    var = fct_reorder(as.factor(var), -val),
    category = case_when(
      var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y",
                 "mean_slope") ~ "spatial",
      var %in% c("det_day", "det_dayx", "det_dayy",
                 "day_night_night", "moon_illuminated") ~ "temporal",
      var %in% c("stage_mature", "fl", "mean_log_e") ~ "biological",
      TRUE ~ "dynamic"
    )
  ) %>% 
  arrange(-val) 
imp_dat$var_f = factor(
  imp_dat$var#, 
  # labels = c("Bottom Depth", "Year Day 1", "Fork Length", "Maturity", 
  #            "UTM Y", "UTM X", "Condition", "Year Day 2", "Bottom Slope", 
  #            "Moon Phase", "Bottom Slope", "Zooplankton",
  #            "Shore Distance", "Temperature", "Oxygen",
  #            "Thermocline Depth", "Day/Night", "H Current 1", "H Current 2",
  #            "Vertical Current")
)

imp_plot <- ggplot(imp_dat, aes(x = var_f, y = val)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(here::here("figs", "ms_figs_rel", "importance_quantreg.png"),
    height = 4, width = 6, units = "in", res = 250)
imp_plot
dev.off()


# VARIABLE IMPORTANCE WITH PARTY PACKAGE ---------------------------------------

## TAKES TOO LONG TO COMPLETE

# rf_party <- party::cforest(
#   depth ~ .,
#   data = train_depth_baked,
#   control = party::cforest_unbiased(
#     mtry = top_mod$tuneValue$mtry,
#     ntree = top_mod$param$num.trees
#     )
# )
# 
# imp1 <- permimp::permimp(rf_party, conditional = TRUE, progressBar = TRUE)
# 

# imp <- rf_party %>%
#   party::varimp(conditional = TRUE) %>% 
#   as_tibble() %>% 
#   rownames_to_column("Feature") %>% 
#   rename(Importance = value)


# COUNTERFACTUAL PREDICT -------------------------------------------------------

# generate predictions for maturity stage and different counterfacs (e.g. 
# most important 4 variables)
new_dat <- train_depth_baked %>%
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
  group_vals <- depth_dat_raw %>% 
      dplyr::summarize(min_v = min(!!varname),
                       max_v = max(!!varname))
  # change to dataframe
  var_seq <- NULL
  for (i in 1:nrow(group_vals)) {
      var_seq <- c(var_seq, 
                   seq(group_vals$min_v[i], group_vals$max_v[i], 
                       length.out = 100))
  }
  new_dat %>% 
    select(- {{ var_in }}) %>% 
    mutate(dum = var_seq) %>% 
    dplyr::rename(!!varname := dum) %>% 
    mutate(
      det_dayx = sin(2 * pi * local_day / 365),
      det_dayy = cos(2 * pi * local_day / 365)
    )
}


# predictions function
pred_foo <- function(preds_in, ...) {
  preds1 <- predict(
    ranger_rf, 
    data = preds_in,
    ...
  )
  
  if (is.null(preds1$se)) {
    preds_out <- preds1$predictions
    colnames(preds_out) <- c("lo", "mean", "up")
  }
  if (!is.null(preds1$se)) {
    preds_out <- cbind(
      preds1$predictions,
      preds1$se
    )
    colnames(preds_out) <- c("mean", "se")
  }
  
  cbind(preds_in, preds_out)
}


# plotting function
plot_foo <- function (data, var) {
  ggplot(data, aes_string(x = var)) +
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = lo, ymax = up), alpha = 0.4) +
    ggsidekick::theme_sleek() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
                       labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                       limits = c(-1, 0))
} 


counterfac_tbl <- tibble(
  var_in = c("mean_bathy", "fl", "local_day", "mean_log_e", "thermo_depth", 
             "roms_temp", "zoo", "oxygen", "shore_dist", 
             "moon_illuminated", "mean_slope")) %>% 
  mutate(
    pred_dat_in = purrr::map(var_in, 
                             gen_pred_dat),
    preds = purrr::map(pred_dat_in, 
                       pred_foo, 
                       type = "quantiles")
  )

plot_list <- purrr::map2(
  counterfac_tbl$var_in, 
  counterfac_tbl$preds, 
  function (var, preds) {
    dum <- preds %>%
      mutate(
        lo = -1 * lo, #(mean + (qnorm(0.025) * se)),
        up = -1 * up, #(mean + (qnorm(0.975) * se)),
        mean = -1 * mean
      )
    plot_foo(data = dum, var = var)
  }
)


## categorical predictions for day/night and maturity impacts
new_dat_trim <- new_dat[1:2, ] %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  ) 
mat_dat <- new_dat_trim %>%  
  mutate(
    stage_mature = c(0, 1),
    stage = factor(stage_mature, labels = c("immature", "mature"))
  ) %>% 
  pred_foo(., 
           type = "quantiles") %>% 
  mutate(
    lo = -1 * lo, #(mean + (qnorm(0.025) * se)),
    up = -1 * up, #(mean + (qnorm(0.975) * se)),
    mean = -1 * mean
  ) 

mat_plot <- ggplot(mat_dat, aes(x = stage)) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
  ggsidekick::theme_sleek() +
  labs(title = "Maturity Counterfactual") +
  scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                     limits = c(-1, 0))

dn_plot <- new_dat_trim %>% 
  mutate(
    day_night_night = c(0, 1)
  ) %>% 
  pred_foo(., 
           type = "quantiles") %>% 
  mutate(
    lo = -1 * lo, #(mean + (qnorm(0.025) * se)),
    up = -1 * up, #(mean + (qnorm(0.975) * se)),
    mean = -1 * mean
  ) %>% 
  ggplot(., aes(x = as.factor(day_night_night))) +
  geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
  ggsidekick::theme_sleek() +
  labs(title = "Day/Night Counterfactual") +
  scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                     limits = c(-1, 0))


pdf(here::here("figs", "ms_figs_rel", "counterfactual_ranger.pdf"),
    height = 6, width = 8)
plot_list
mat_plot
dn_plot
dev.off()




## cleaned up versions for manuscript main text
# calculate overall depth range for y axis
bathy_cond <- plot_list[[1]] + 
  labs(x = "Mean Bathymetry Depth") +
  theme(
    axis.title.y = element_blank()
  )
yday_cond <- plot_list[[3]] + 
  labs(x = "Year Day") +
  theme(
    axis.title.y = element_blank()
  )
# slope_cond <- plot_list[[11]] + 
#   labs(x = "Mean Slope") +
#   theme(
#     axis.title.y = element_blank()
#   )
# moon_cond <- plot_list[[10]] + 
#  labs(x = "Proportion of Moon Illuminated") +
#   theme(
#     axis.title.y = element_blank()
#   )

panel1 <- cowplot::plot_grid(
  yday_cond,
  bathy_cond,
  # slope_cond,
  # moon_cond, 
  nrow = 1
)

y_grob <- grid::textGrob("Predicted Bathymetric\nDepth Ratio", rot = 90)

png(here::here("figs", "ms_figs_rel", "counterfac_effects.png"),
    height = 2.5, width = 5.25, 
    units = "in", res = 250)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    panel1, left = y_grob)
)
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
  select(utm_x, utm_y, mean_bathy = depth, max_bathy = max_depth,
         mean_slope = slope, shore_dist) %>% 
  filter(!mean_bathy > 400,
         !max_bathy > 500,
         utm_y > 5100)

ggplot(bath_grid) +
  geom_raster(aes(x = utm_x, y = utm_y, fill = mean_bathy))

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

pred_rf1 <- predict(ranger_rf,
                    type = "quantiles",
                    quantiles = c(0.1, 0.5, 0.9),
                    data = pred_dat1)
colnames(pred_rf1$predictions) <- c("lo", "med", "up")

pred_dat2 <- pred_dat1 %>%
  mutate(
    rel_pred_med = pred_rf1$predictions[, "med"],
    rel_pred_lo = pred_rf1$predictions[, "lo"],
    rel_pred_up = pred_rf1$predictions[, "up"],
    pred_med = pred_rf1$predictions[, "med"] * max_bathy,
    pred_lo = pred_rf1$predictions[, "lo"] * max_bathy,
    pred_up = pred_rf1$predictions[, "up"] * max_bathy,
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000,
    pred_int_width = pred_up - pred_lo
  ) %>% 
  filter(mean_bathy < 400)

base_plot <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
  
rel_depth <- base_plot +
  geom_raster(data = pred_dat2, 
              aes(x = utm_x_m, y = utm_y_m, fill = rel_pred_med)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Bathymetric\nDepth Ratio",
                       direction = -1, 
                       option = "A") +
  theme(legend.position = "top")

mean_depth <- base_plot +
  geom_raster(data = pred_dat2, 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Mean Depth",
                       direction = -1) +
  theme(legend.position = "top")

var_depth <- base_plot +
  geom_raster(data = pred_dat2, 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Prediction\nInterval Width", 
                       option = "C",
                       direction = -1)  +
  theme(legend.position = "top")

avg_depth <- cowplot::plot_grid(rel_depth, mean_depth, var_depth, ncol = 3)

png(here::here("figs", "ms_figs_rel", "avg_depth.png"), height = 5, width = 8, 
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
      )
    ,
    # generate predictions
    preds = purrr::map(pred_grid, function (x) {
      pred_rf <- predict(ranger_rf, 
                         type = "quantiles",
                         quantiles = c(0.1, 0.5, 0.9),
                         data = x, 
                         all = TRUE)
      colnames(pred_rf$predictions) <- c("lo", "med", "up")
      
      bath_grid %>%
        mutate(
          rel_pred_med = pred_rf$predictions[, "med"],
          rel_pred_lo = pred_rf$predictions[, "lo"],
          rel_pred_up = pred_rf$predictions[, "up"],
          pred_med = pred_rf$predictions[, "med"] * max_bathy,
          pred_lo = pred_rf$predictions[, "lo"] * max_bathy,
          pred_up = pred_rf$predictions[, "up"] * max_bathy,
          utm_x_m = utm_x * 1000,
          utm_y_m = utm_y * 1000
        )
    }
    )
  )
  
pred_dat <- pred_tbl %>% 
  select(contrast, season, #stage_mature, , moon_illuminated, day_night_night,
         preds) %>% 
  unnest(cols = preds)


# calculate differences for each contrast scenario (NOTE: estimates unaffected
# by whether real or scaled predictions are used)
season_eff <- pred_dat %>% 
  filter(contrast == "season") %>% 
  select(season, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med) %>% 
  pivot_wider(names_from = season, values_from = pred_med) %>%
  mutate(season_diff = (winter - summer) / summer) 
season_map <- base_plot +
  geom_raster(data = season_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = season_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  # scale_fill_gradientn(colours = c("red", "white", "blue"),
  #                      values = scales::rescale(c(-3, 0, 3)),
  #                      guide = "colorbar",
  #                      limits = c(-3, 3)) +
  labs(title = "Season Effects")


mat_eff <- pred_dat %>% 
  filter(contrast == "maturity") %>% 
  select(stage_mature, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med) %>% 
  pivot_wider(names_from = stage_mature, values_from = pred_med) %>%
  mutate(mat_diff = (`0` - `1`) / `1`)
mat_map <- base_plot +
  geom_raster(data = mat_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = mat_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Maturity Effects")


moon_eff <- pred_dat %>% 
  filter(contrast == "moon light") %>% 
  select(moon_illuminated, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med) %>% 
  pivot_wider(names_from = moon_illuminated, values_from = pred_med) %>%
  mutate(moon_diff = (`0` - `1`) / `1`)
moonlight_map <- base_plot +
  geom_raster(data = moon_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = moon_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Moonlight Effects")


dvm_eff <- pred_dat %>% 
  filter(contrast == "dvm") %>% 
  select(day_night_night, mean_bathy:shore_dist, utm_x_m, utm_y_m, pred_med) %>% 
  pivot_wider(names_from = day_night_night, values_from = pred_med) %>%
  # negative is deeper at night relative to night, pos shallower at night rel to night
  mutate(dvm_diff = (`0` - `1`) / `1`) 
dvm_map <- base_plot +
  geom_raster(data = dvm_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = dvm_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "DVM Effects")


pdf(here::here("figs", "depth_ml", "spatial_contrasts_relative.pdf"))
season_map
mat_map
moonlight_map
dvm_map
dev.off()


# combine season and maturity predictions and plot joined version
season_eff2 <- season_eff %>% 
  mutate(comp = "season") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = season_diff)
mat_eff2 <- mat_eff %>% 
  mutate(comp = "maturity") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = mat_diff)
moon_eff2 <- moon_eff %>% 
  mutate(comp = "moonlight") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = moon_diff)
dvm_eff2 <- dvm_eff %>% 
  mutate(comp = "dvm") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = dvm_diff)

comb_preds <- list(season_eff2,
                   mat_eff2,
                   moon_eff2,
                   dvm_eff2) %>% 
  bind_rows()

png(here::here("figs", "ms_figs", "contrast_map.png"), height = 5, width = 8, 
    units = "in", res = 250)
base_plot +
  geom_raster(data = comb_preds, 
              aes(x = utm_x_m, y = utm_y_m, fill = rel_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Relative Depth Difference"
  ) +
  facet_wrap(~comp) +
  theme(legend.position = "top")
dev.off()
