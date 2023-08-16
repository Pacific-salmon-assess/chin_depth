## Fit Regression Random Forest w/ Percentile Intervals

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees and fit to relative depth 
#data. Fit equivalent model with interpolated training data and generate 
#quantile prediction intervals.
# July 8 (fit model with stock as covariate w/ very minimal improvement in 
# performance; fits in depth_quantrf_stock.rds)
# Nov 22 fit using ranger rather than quantregForest to be consistent with
# model comparison
# Nov 23 switch to transformed depth based on updated model comparison
# May 25 add last batch of detections and additional ROMS data

library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)


depth_dat_raw1 <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  # approximately 7k detections have no available ROMS data; exclude 
  filter(!is.na(roms_temp))


# remove 2022 tag releases (~6k dets) for training model
depth_dat_raw <- depth_dat_raw1 %>% 
  filter(!grepl("2022", vemco_code))



## FIT -------------------------------------------------------------------------

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
    depth = rel_depth, fl, lipid, stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block
  )  

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()
test_depth_22 <- depth_dat_raw1 %>% 
  filter(grepl("2022", vemco_code)) %>% 
  dplyr::select(
    depth = rel_depth, fl, lipid, stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated)


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

# apply recipe to dataframe to make dummy variables and 
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block, -max_bathy))

#pull model attributes from top ranger
# rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))
# top_mod <- rf_list[[2]]$top_model
# 
# ranger_rf <- ranger::ranger(
#   depth ~ .,
#   data = train_depth_baked,
#   #hyperpars based on values from top model which is not saved on all locals
#   num.trees = 2000,
#   mtry = 5,
#   keep.inbag = TRUE,
#   quantreg = TRUE,
#   importance = "permutation"
# )

# saveRDS(ranger_rf,
#         here::here("data", "model_fits", "relative_rf_ranger.rds"))

ranger_rf <- readRDS(here::here("data", "model_fits", "relative_rf_ranger.rds"))
# ranger_rf <- readRDS(here::here("data", "model_fits", "relative_rf_ranger_FULL.rds"))


# CHECK PREDS ------------------------------------------------------------------

obs_preds <- predict(ranger_rf,
                     data = train_depth_baked)

dum <- train_depth %>% 
  mutate(mean_pred = obs_preds$predictions,
         mean_pred_real = mean_pred * max_bathy,
         depth_real = depth * max_bathy,
         split_group = "train 2019-21")
plot(depth ~ mean_pred, dum)
plot(depth_real ~ mean_pred_real, dum)


# hold out predictions
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
         split_group = "test 2019-21")


# 2022 predictions
test_depth_baked_22 <- prep(depth_recipe) %>%
  bake(., 
       new_data = test_depth_22 %>% 
         dplyr::select(-max_bathy))
test_preds_22<- predict(ranger_rf,
                      data = test_depth_baked_22)
dum_test_22 <- test_depth_22 %>% 
  mutate(ind_block = NA,
         mean_pred = test_preds_22$predictions,
         mean_pred_real = mean_pred * max_bathy,
         depth_real = depth * max_bathy,
         split_group = "test 2022")

all_preds <- do.call(rbind,
                     list(dum, dum_test, dum_test_22)) %>% 
  mutate(
    split_group = fct_relevel(
      as.factor(split_group), "train 2019-21", after = 0)
  )

fit_obs <- ggplot() +
  geom_point(
    data = all_preds,
    aes(x = depth_real, y = mean_pred_real, fill = split_group),
    shape = 21, alpha = 0.025
    ) +
  labs(
    x = "Observed Depth (m)", y = "Predicted Mean Depth (m)"
  ) +
  scale_fill_discrete(name = "") +
  ggsidekick::theme_sleek() +
  facet_wrap(~split_group) +
  theme(legend.position = "none")

png(here::here("figs", "ms_figs_rel", "obs_preds_rel.png"),
    height = 3, width = 6, units = "in", res = 200)
fit_obs
dev.off()

dum_test_22$resid <- dum_test_22$mean_pred_real - dum_test_22$depth_real
hist(dum_test_22$resid)


# rmse of each group
Metrics::rmse(dum$depth, dum$mean_pred)
Metrics::rmse(dum_test$depth, dum_test$mean_pred)


# VARIABLE IMPORTANCE ----------------------------------------------------------

imp_vals <- ranger::importance(ranger_rf, type = "permutation", scale = F) 

# key for axis labels
var_name_key <- data.frame(
  var = names(imp_vals),
  var_f = c("Fork Length", "Lipid Content", "UTM X", "UTM Y", "Year Day 2", 
      "Year Day 1", "Bottom Depth", "Bottom Slope", "Shore Distance", 
      "Hor. Current 1", "Hor. Current 2", "Vert. Current", "Temperature", 
      "Zooplankton", "Oxygen", "Thermocline Depth", "Lunar Cycle", "Maturity",
      "Day/Night")
)


imp_dat <- data.frame(
  var = names(imp_vals),
  val = imp_vals
  ) %>% 
  left_join(., var_name_key, by = "var") %>% 
  arrange(-val) %>% 
  mutate(
    var = fct_reorder(as.factor(var), -val),
    category = case_when(
      var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y",
                 "mean_slope") ~ "spatial",
      var %in% c("det_day", "det_dayx", "det_dayy",
                 "day_night_night", "moon_illuminated") ~ "temporal",
      var %in% c("stage_mature", "fl", "lipid") ~ "biological",
      TRUE ~ "dynamic"
    )
  ) 
  

imp_plot <- ggplot(imp_dat, aes(x = fct_reorder(var_f, - val), y = val)) +
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


# SPATIAL PREDICT --------------------------------------------------------------

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


# REPLACED BY ESTIMATED VALUES FROM ROMS
# calculate mean roms_variables for different seasons (using monthly averages)
roms_month_means <- readRDS(here::here("data", "depth_dat_nobin.RDS")) %>%
  mutate(month = lubridate::month(date_time_local)) %>%
  filter(month %in% c(1, #4, 
                      7)) %>%
  mutate(season = ifelse(month == "1", "winter", "summer")) %>%
  group_by(season) %>%
  dplyr::summarize(
    thermo_depth = mean(thermo_depth, na.rm = T)
    )


# input grid includes ROMS data for entire study area for specified dates
bath_grid_in <- readRDS(here::here("data", "pred_bathy_grid_roms.RDS")) 
bath_grid <- bath_grid_in %>% 
  filter(!mean_bathy > 400,
         !max_bathy > 500,
         utm_y > 5100) %>% 
  left_join(., roms_month_means, by = "season")


ggplot(bath_grid %>% filter(local_day == "211")) +
  geom_raster(aes(x = utm_x, y = utm_y, fill = mean_bathy))


# base template plot for maps
base_plot <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  # set limitsto avoid boundary effects from UTM projection
  scale_x_continuous(
    limits = c(210000, 560000), 
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(5100000, 5470000), expand = c(0, 0))



## first set are average spatial predictions (i.e. mean biological and temporal 
## attributes)

# biological data
bio_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, lipid, stage) %>% 
  distinct()

# stratify predictions by non-spatial covariates
pred_dat1 <- bath_grid %>% 
  mutate(
    fl = mean(bio_dat$fl),
    lipid = mean(bio_dat$lipid),
    stage_mature = 0.5
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
  ) 

  
rel_depth <- base_plot +
  geom_raster(data = pred_dat2 %>% filter(season == "summer"), 
              aes(x = utm_x_m, y = utm_y_m, fill = rel_pred_med)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Bathymetric\nDepth Ratio",
                       direction = -1, 
                       option = "A") +
  theme(legend.position = "top",
        axis.text = element_blank()) 

mean_depth <- base_plot +
  geom_raster(data = pred_dat2 %>% filter(season == "summer"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Mean Depth (m)",
                       direction = -1) +
  theme(legend.position = "top")

var_depth <- base_plot +
  geom_raster(data = pred_dat2 %>% filter(season == "summer"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Prediction\nInt. Width (m)", 
                       option = "C",
                       direction = -1)  +
  theme(legend.position = "top",
        axis.text = element_blank())

avg_depth1 <- cowplot::plot_grid(plotlist = list(rel_depth, var_depth),
                                 ncol = 1)
avg_depth <- cowplot::plot_grid(plotlist = list(mean_depth, avg_depth1),
                                ncol = 2,
                                rel_widths = c(1.5, 1))

png(here::here("figs", "ms_figs_rel", "avg_depth.png"), height = 5.5, width = 7, 
    units = "in", res = 250)
avg_depth
dev.off()


## generate spatial contrasts
# 1) seasonal effects for immature fish
# 2) maturity effects for avg (July 1)
# 3) moon illumination effects for avg  (July 1)
stage_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, lipid, stage) %>% 
  distinct() %>% 
  group_by(stage) %>% 
  dplyr::summarize(
    fl = mean(fl),
    lipid = mean(lipid)
  ) %>% 
  mutate(stage_mature = ifelse(stage == "mature", 1, 0)) %>% 
  rbind(., 
        data.frame(
          stage = "average",
          fl = mean(depth_dat_raw$fl),
          lipid = mean(depth_dat_raw$lipid),
          stage_mature = 0.5
        ))

# subset bath_grid to remove covariates defined in predictive tibble
bath_grid_trim <- bath_grid %>% 
  select(
    -c(local_day, moon_illuminated, det_dayx, det_dayy, day_night_night)
  ) %>% 
  group_by(season) %>% 
  group_nest()

depth_dat_summer <- depth_dat_raw %>% 
  filter(local_day > 152 & local_day < 243) 
t_range <- c(mean(depth_dat_summer$roms_temp) - sd(depth_dat_summer$roms_temp),
             mean(depth_dat_summer$roms_temp) + sd(depth_dat_summer$roms_temp))
z_range <- c(mean(depth_dat_summer$zoo) - sd(depth_dat_summer$zoo),
             mean(depth_dat_summer$zoo) + sd(depth_dat_summer$zoo))

pred_tbl <- tibble(
  contrast = rep(c("season", "maturity", "moon light", 
                   "day night", "temp", "zoo"), each = 2),
  local_day = c(46, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211),
  stage_mature = c(0, 0, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  moon_illuminated = c(0.5, 0.5, 0.5, 0.5, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  day_night_night = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 1, 0.5, 0.5, 0.5, 0.5),
  #leave NaNs to be replaced by bath_grid_trim except for temp and zoo contrasts
  adj_temp = c(rep(NaN, 8), t_range[1], t_range[2], NaN, NaN),
  adj_zoo = c(rep(NaN, 10), z_range[1], z_range[2])
) %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365),
    season = fct_recode(as.factor(local_day), 
                        "winter" = "46", "summer" = "211"),
    month = factor(as.factor(local_day), labels = c("1", "7"))) %>% 
  left_join(., bath_grid_trim, by = "season") %>%
  left_join(., stage_dat, by = "stage_mature") %>% 
  unnest(cols = c(data)) %>% 
  # replace temp adn zoo data for those specific contrasts
  mutate(
    roms_temp = ifelse(!is.na(adj_temp), adj_temp, roms_temp),
    zoo = ifelse(!is.na(adj_zoo), adj_zoo, zoo)
  ) %>% 
  group_by(contrast, .add = TRUE) %>% 
  group_nest(.key = "pred_grid") %>%
  mutate(
    # generate predictions
    preds = purrr::map(pred_grid, function (x) {
      pred_rf <- predict(ranger_rf,
                         type = "quantiles",
                         quantiles = c(0.1, 0.5, 0.9),
                         data = x,
                         all = TRUE)
      colnames(pred_rf$predictions) <- c("lo", "med", "up")
      
      x %>%
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
  select(-pred_grid) %>%
  unnest(cols = preds)

# ## Show raw differences between plots
# season_map_abs <- base_plot +
#   geom_raster(data = pred_dat %>% 
#                 filter(contrast == "season"), 
#               aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
#   geom_sf(data = coast_utm) +
#   scale_fill_viridis_c(name = "Mean Depth", direction = -1) +
#   labs(title = "Season Effects") +
#   facet_wrap(~season)
# 
# mat_names <- c(
#   `0` = "Immature",
#   `1` = "Mature"
# )
# mat_map_abs <- base_plot +
#   geom_raster(data = pred_dat %>% 
#                 filter(contrast == "maturity"), 
#               aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
#   geom_sf(data = coast_utm) +
#   scale_fill_viridis_c(name = "Mean Depth", direction = -1) +
#   labs(title = "Life Stage Effects") +
#   facet_wrap(~stage_mature,
#              labeller = as_labeller(mat_names))
# 
# moon_map_abs <- base_plot +
#   geom_raster(data = pred_dat %>% 
#                 filter(contrast == "moon light"), 
#               aes(x = utm_x_m, y = utm_y_m, fill = pred_med)) +
#   geom_sf(data = coast_utm) +
#   scale_fill_viridis_c(name = "Mean Depth", direction = -1) +
#   labs(title = "Lunar Illumination") +
#   facet_wrap(~moon_illuminated)
# 
# png(here::here("figs", "ms_figs_rel", "contrast_map_abs.png"), 
#     height = 5.5, width = 5, 
#     units = "in", res = 250)
# cowplot::plot_grid(
#   plotlist = list(season_map_abs, mat_map_abs#, moon_map_abs, dvm_map_abs
#                   ),
#   ncol = 1)
# dev.off()


## Show relative differences 

# calculate differences for each contrast scenario (NOTE: estimates unaffected
# by whether real or scaled predictions are used)
season_eff <- pred_dat %>% 
  filter(contrast == "season") %>% 
  select(season, mean_bathy:shore_dist, utm_x_m, utm_y_m, rel_pred_med) %>% 
  pivot_wider(names_from = season, values_from = rel_pred_med) %>%
  mutate(season_diff = (winter - summer)) 
season_map <- base_plot +
  geom_raster(data = season_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = season_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Season Effects")


mat_eff <- pred_dat %>% 
  filter(contrast == "maturity") %>% 
  select(stage_mature, mean_bathy:shore_dist, utm_x_m, utm_y_m, rel_pred_med) %>% 
  pivot_wider(names_from = stage_mature, values_from = rel_pred_med) %>%
  mutate(mat_diff = (`0` - `1`))
mat_map <- base_plot +
  geom_raster(data = mat_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = mat_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Maturity Effects")


moon_eff <- pred_dat %>% 
  filter(contrast == "moon light") %>% 
  select(moon_illuminated, mean_bathy:shore_dist, utm_x_m, utm_y_m, rel_pred_med) %>%
  pivot_wider(names_from = moon_illuminated, values_from = rel_pred_med) %>%
  # negative is deeper with full moonlight, pos deeper with no moonlight 
  mutate(moon_diff = (`0` - `1`))
moonlight_map <- base_plot +
  geom_raster(data = moon_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = moon_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Moonlight Effects")


sst_eff <- pred_dat %>%
  filter(contrast == "temp") %>%
  mutate(roms_temp_f = factor(roms_temp, labels = c("low", "high"))) %>% 
  select(roms_temp_f, mean_bathy:shore_dist, utm_x_m, utm_y_m,
         rel_pred_med) %>%
  pivot_wider(names_from = roms_temp_f, values_from = rel_pred_med) %>%
  # negative is deeper at higher temps, positive deeper at lower temps
  mutate(temp_diff = (low - high))
sst_map <- base_plot +
  geom_raster(data = sst_eff,
              aes(x = utm_x_m, y = utm_y_m, fill = temp_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "SST Effects")

zoo_eff <- pred_dat %>%
  filter(contrast == "zoo") %>%
  mutate(zoo_f = factor(zoo, labels = c("low", "high"))) %>% 
  select(zoo_f, mean_bathy:shore_dist, utm_x_m, utm_y_m,
         rel_pred_med) %>%
  pivot_wider(names_from = zoo_f, values_from = rel_pred_med) %>%
  # negative is deeper at higher zoo conc, positive deeper at lower zoo
  mutate(zoo_diff = (low - high))
zoo_map <- base_plot +
  geom_raster(data = zoo_eff,
              aes(x = utm_x_m, y = utm_y_m, fill = zoo_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Zooplankton Effects")

dn_eff <- pred_dat %>% 
  filter(contrast == "day night") %>% 
  select(day_night_night, mean_bathy:shore_dist, utm_x_m, utm_y_m, rel_pred_med) %>%
  pivot_wider(names_from = day_night_night, values_from = rel_pred_med) %>%
  # negative is deeper at night, pos deeper during day 
  mutate(dn_diff = (`0` - `1`))
dn_map <- base_plot +
  geom_raster(data = dn_eff, 
              aes(x = utm_x_m, y = utm_y_m, fill = dn_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2() +
  labs(title = "Day/Night Effects")


# combine season and maturity predictions and plot joined version
# season_eff2 <- season_eff %>% 
#   mutate(comp = "season") %>% 
#   select(mean_bathy:utm_y_m, comp, rel_diff = season_diff)
mat_eff2 <- mat_eff %>% 
  mutate(comp = "maturity") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = mat_diff)
moon_eff2 <- moon_eff %>% 
  mutate(comp = "moonlight") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = moon_diff)
sst_eff2 <- sst_eff %>% 
  mutate(comp = "roms_temp") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = temp_diff)
zoo_eff2 <- zoo_eff %>% 
  mutate(comp = "zoo") %>% 
  select(mean_bathy:utm_y_m, comp, rel_diff = zoo_diff)

comb_preds <- list(
  moon_eff2,
  sst_eff2,
  zoo_eff2,
  mat_eff2
) %>% 
  bind_rows() %>% 
  mutate(
    comp = factor(
      comp, levels = c("moonlight", "roms_temp", "zoo", "maturity"),
      labels = c("lunar cycle", "SST", "zooplankton", "maturity/size/lipid")
    )
  )

png(here::here("figs", "ms_figs_rel", "contrast_map.png"), 
    height = 5.5, width = 5, 
    units = "in", res = 250)
base_plot +
  geom_raster(data = comb_preds, 
              aes(x = utm_x_m, y = utm_y_m, fill = rel_diff)) +
  geom_sf(data = coast_utm) +
  scale_fill_gradient2(
    name = "Bathymetric\nDepth Ratio\nDifference"
  ) +
  facet_wrap(~comp) +
  theme(legend.position = "top") +
  theme(
    axis.text = element_text(size = 8)
  )
dev.off()



# day-night contrast
pred_dat %>% 
  filter(contrast == "day night") %>% 
  mutate(rel_int_width = rel_pred_up - rel_pred_lo) %>% 
  select(day_night_night, mean_bathy:shore_dist, utm_x_m, utm_y_m, 
         rel_pred_med, rel_int_width) %>% 
  group_by(day_night_night) %>% 
  summarize(mean_rel_depth = mean(rel_pred_med),
            mean_int_width = mean(rel_int_width))



## LATENT SPATIAL PROCESSES ----------------------------------------------------

# set non-coordinate spatial variables to mean values 
pred_latent <- pred_dat1 %>% 
  filter(season == "summer") %>% 
  mutate(
    mean_bathy = mean(mean_bathy),
    mean_slope = mean(mean_slope),
    shore_dist = mean(shore_dist),
    oxygen = mean(oxygen),
    roms_temp = mean(roms_temp),
    u = mean(u),
    v = mean(v),
    w = mean(w),
    zoo = mean(zoo),
    thermo_depth = mean(thermo_depth)
  )

pred_latent_rf <- predict(ranger_rf,
                    type = "quantiles",
                    quantiles = c(0.1, 0.5, 0.9),
                    data = pred_latent)

colnames(pred_latent_rf$predictions) <- c("lo", "med", "up")

pred_latent2 <- pred_latent %>%
  mutate(
    rel_pred_med = pred_latent_rf$predictions[, "med"],
    rel_pred_lo = pred_latent_rf$predictions[, "lo"],
    rel_pred_up = pred_latent_rf$predictions[, "up"],
    pred_med = pred_latent_rf$predictions[, "med"] * max_bathy,
    pred_lo = pred_latent_rf$predictions[, "lo"] * max_bathy,
    pred_up = pred_latent_rf$predictions[, "up"] * max_bathy,
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000,
    pred_int_width = pred_up - pred_lo
  ) 

# plot that is added to counterfac panel below
rel_latent <- base_plot +
  geom_raster(data = pred_latent2, 
              aes(x = utm_x_m, y = utm_y_m, fill = rel_pred_med)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Predicted\nBathymetric\nDepth Ratio",
                       direction = -1, 
                       option = "A") 


# COUNTERFACTUAL PREDICT -------------------------------------------------------

new_dat <- data.frame(
    utm_x = median(train_depth_baked$utm_x),
    utm_y = median(train_depth_baked$utm_y),
    fl = median(train_depth_baked$fl),
    lipid = median(train_depth_baked$lipid),
    mean_bathy = median(train_depth_baked$mean_bathy),
    local_day = median(depth_dat_raw$local_day),
    mean_slope = median(train_depth_baked$mean_slope),
    shore_dist = median(train_depth_baked$shore_dist),
    u = median(train_depth_baked$u),
    v = median(train_depth_baked$v),
    w = median(train_depth_baked$w),
    zoo = median(train_depth_baked$zoo),
    oxygen = median(train_depth_baked$oxygen),
    thermo_depth = median(train_depth_baked$thermo_depth),
    roms_temp = median(train_depth_baked$roms_temp),
    moon_illuminated = median(train_depth_baked$moon_illuminated),
    day_night_night = 0.5,
    stage_mature = 0.5
  ) %>% 
  #duplicate 100 times
  as_tibble() %>% 
  slice(rep(1:n(), each = 100))


# make tibble for different counterfacs
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
  # var_seq2 <- rep(var_seq, times = length(unique(new_dat$site)))
  var_seq2 <- var_seq
  new_dat %>% 
    select(- {{ var_in }}) %>% 
    mutate(dum = var_seq2) %>% 
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
    colnames(preds_out) <- c("lo", "med", "up")
    dum <- cbind(preds_in, preds_out) %>% 
      mutate(med = -1 * med)
  }
  if (!is.null(preds1$se)) {
    preds_out <- cbind(
      preds1$predictions,
      preds1$se
    ) 
    colnames(preds_out) <- c("mean", "se")
    dum <- cbind(preds_in, preds_out) %>%
      mutate(
        lo = mean + (qnorm(0.025) * se),
        up = mean + (qnorm(0.975) * se), 
        mean = -1 * mean
      )
  }
  
  dum %>%
    mutate(
      lo = -1 * lo, 
      up = -1 * up
    )
}

counterfac_tbl <- tibble(
  var_in = c("mean_bathy", "fl", "utm_x", "utm_y",
             "local_day", "lipid", "thermo_depth",
             "roms_temp", "zoo", "oxygen", "shore_dist",
             "moon_illuminated", "mean_slope"
  )) %>%
  mutate(
    pred_dat_in = purrr::map(var_in,
                             gen_pred_dat),
    preds_ci = purrr::map(pred_dat_in,
                          pred_foo,
                          type = "se",
                          se.method = "infjack")
  )
saveRDS(counterfac_tbl,
        here::here("data", "counterfac_preds_ci.rds"))
counterfac_tbl <- readRDS(here::here("data", "counterfac_preds_ci.rds"))



# plotting function
plot_foo <- function (data, ...) {
  ggplot(data, mapping = aes(!!!ensyms(...))) +
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = lo, ymax = up), alpha = 0.4) +
    ggsidekick::theme_sleek() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
                       labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                       limits = c(-0.75, -0.05))
}


bathy_cond <- counterfac_tbl %>% 
  filter(var_in == "mean_bathy") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "mean_bathy") +  
  labs(x = "Mean Bottom Depth") +
  theme(
    axis.title.y = element_blank()
  ) 

yday_cond <- counterfac_tbl %>% 
  filter(var_in == "local_day") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "local_day") +  
  labs(x = "Year Day") +
  theme(
    axis.title.y = element_blank()
  ) 


panel1 <- cowplot::plot_grid(
  yday_cond,
  bathy_cond,
  nrow = 1
)
pp <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    panel1, 
    left = grid::textGrob("Predicted Bathymetric\nDepth Ratio", rot = 90)
  )
)

png(here::here("figs", "ms_figs_rel", "counterfac_effects.png"),
    height = 5.5, width = 5.5, 
    units = "in", res = 250)
cowplot::plot_grid(
  pp,
  rel_latent,
  nrow = 2
) 
dev.off()


temp_cond <- counterfac_tbl %>% 
  filter(var_in == "roms_temp") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "roms_temp") +  
  theme(
    axis.title.y = element_blank()
  ) 
zoo_cond <- counterfac_tbl %>% 
  filter(var_in == "zoo") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "zoo") +  
  theme(
    axis.title.y = element_blank()
  ) 
moon_cond <- counterfac_tbl %>% 
  filter(var_in == "moon_illuminated") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "moon_illuminated") +  
  theme(
    axis.title.y = element_blank()
  ) 


# pdf for region specific relationships
# pdf(here::here("figs", "ms_figs_rel", "bathy_counterfac_region.pdf"),
#     height = 6, width = 8)
# ggplot(counterfac_tbl$preds[[1]], aes(x = mean_bathy, colour = site)) +
#   geom_line(aes(y = mean)) +
#   geom_ribbon(aes(ymin = lo, ymax = up, fill = site), alpha = 0.2) +
#   ggsidekick::theme_sleek() +
#   scale_x_continuous(expand = c(0, 0)) +
#   facet_wrap(~site)
# dev.off()


## categorical predictions for day/night and maturity impacts
# new_dat_trim <- new_dat[1:2, ] %>% 
#   mutate(
#     det_dayx = sin(2 * pi * local_day / 365),
#     det_dayy = cos(2 * pi * local_day / 365)
#   ) 
# mat_dat <- new_dat_trim %>%  
#   mutate(
#     stage_mature = c(0, 1),
#     stage = factor(stage_mature, labels = c("immature", "mature"))
#   ) %>% 
#   pred_foo(., 
#            type = "quantiles") %>% 
#   mutate(
#     lo = -1 * lo, #(mean + (qnorm(0.025) * se)),
#     up = -1 * up, #(mean + (qnorm(0.975) * se)),
#     mean = -1 * mean
#   ) 
# 
# mat_plot <- ggplot(mat_dat, aes(x = stage)) +
#   geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
#   ggsidekick::theme_sleek() +
#   labs(title = "Maturity Counterfactual") +
#   scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
#                      labels = c("0", "0.25", "0.5", "0.75", "1.0"),
#                      limits = c(-1, 0))
# 
# dn_plot <- new_dat_trim %>% 
#   mutate(
#     day_night_night = c(0, 1)
#   ) %>% 
#   pred_foo(., 
#            type = "quantiles") %>% 
#   mutate(
#     lo = -1 * lo, #(mean + (qnorm(0.025) * se)),
#     up = -1 * up, #(mean + (qnorm(0.975) * se)),
#     mean = -1 * mean
#   ) %>% 
#   ggplot(., aes(x = as.factor(day_night_night))) +
#   geom_pointrange(aes(y = mean, ymin = lo, ymax = up)) +
#   ggsidekick::theme_sleek() +
#   labs(title = "Day/Night Counterfactual") +
#   scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1.0),
#                      labels = c("0", "0.25", "0.5", "0.75", "1.0"),
#                      limits = c(-1, 0))
# 
# 
# pdf(here::here("figs", "ms_figs_rel", "counterfactual_ranger.pdf"),
#     height = 6, width = 8)
# plot_list
# mat_plot
# dn_plot
# dev.off()
# 

