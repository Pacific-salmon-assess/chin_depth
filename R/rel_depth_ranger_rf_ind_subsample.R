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
  filter(!is.na(roms_temp)) %>% 
  group_by(vemco_code) %>% 
  mutate(n_dets = n()) 

high_dets_subset <- depth_dat_raw1 %>% 
  filter(n_dets > 750) %>% 
  slice_sample(n = 750)

depth_dat_raw2 <- depth_dat_raw1 %>% 
  filter(!n_dets < 750) %>% 
  rbind(., high_dets_subset) %>% 
  ungroup()


# remove 2022 tag releases (~6k dets) for training model
depth_dat_raw <- depth_dat_raw2 %>% 
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
test_depth_22 <- depth_dat_raw2 %>% 
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
ranger_rf <- ranger::ranger(
  depth ~ .,
  data = train_depth_baked,
  #hyperpars based on values from top model which is not saved on all locals
  num.trees = 1500,
  mtry = 13,
  # keep.inbag = TRUE for quantile predictions
  keep.inbag = TRUE,
  quantreg = TRUE,
  importance = "permutation",
  splitrule = "extratrees"
)

saveRDS(ranger_rf,
        here::here("data", "model_fits", "relative_rf_ranger_subset.rds"))

ranger_rf <- readRDS(
  here::here("data", "model_fits", "relative_rf_ranger_subset.rds"))


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
  
shape_pal <- c(21, 22, 23, 24)
names(shape_pal) <- levels(imp_dat$category)

imp_plot <- ggplot(imp_dat, aes(x = fct_reorder(var_f, - val), y = val)) +
  geom_point(aes(fill = category, shape = category), size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  scale_shape_manual(values = shape_pal, guide = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 0.08),
    expand = c(0, 0)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = shape_pal))
  )

png(here::here("figs", "ms_figs_rel", "importance_quantreg.png"),
    height = 100, width = 170, units = "mm", res = 300)
imp_plot
dev.off()


# SPATIAL PREDICT --------------------------------------------------------------

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


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
  theme(legend.position = "top",
        axis.text = element_blank()) 
  
var_depth <- base_plot +
  geom_raster(data = pred_dat2 %>% filter(season == "summer"), 
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "80% Prediction\nInt. Width (m)", 
                       option = "C",
                       direction = -1)  +
  theme(legend.position = "top",
        axis.text = element_blank())

avg_depth1 <- cowplot::plot_grid(plotlist = list(rel_depth, mean_depth),
                                 ncol = 2)
# avg_depth <- cowplot::plot_grid(plotlist = list(mean_depth, avg_depth1),
#                                 ncol = 2,
#                                 rel_widths = c(1.5, 1))

png(here::here("figs", "ms_figs_rel", "avg_depth.png"), 
    height = 110, width = 170, units = "mm", res = 300)
avg_depth1
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

pred_latent_rf <- predict(
  ranger_rf,
  type = "quantiles",
  quantiles = c(0.1, 0.5, 0.9),
  data = pred_latent
)

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
                       option = "A") +
  geom_text(aes(x = -Inf, y = Inf, label = "f)"), hjust = -0.5, vjust = 1.5) +
  theme(legend.position = c(0.15, 0.27),#"left",
        legend.key.size = unit(0.75, 'cm'),
        axis.text = element_blank(),
        legend.title = element_blank())


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

# counterfac_tbl <- tibble(
#   var_in = c("mean_bathy", "fl", "utm_x", "utm_y",
#              "local_day", "lipid", "thermo_depth",
#              "roms_temp", "zoo", "oxygen", "shore_dist",
#              "moon_illuminated", "mean_slope"
#   )) %>%
#   mutate(
#     pred_dat_in = purrr::map(var_in,
#                              gen_pred_dat),
#     preds_ci = purrr::map(pred_dat_in,
#                           pred_foo,
#                           type = "se",
#                           se.method = "infjack")
#   )
# saveRDS(counterfac_tbl,
#         here::here("data", "counterfac_preds_ci.rds"))
counterfac_tbl <- readRDS(here::here("data", "counterfac_preds_ci.rds"))


# plotting function
plot_foo <- function (data, ...) {
  ggplot(data, mapping = aes(!!!ensyms(...))) +
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = lo, ymax = up), alpha = 0.4) +
    ggsidekick::theme_sleek() +
    scale_x_continuous(expand = c(0, 0))+
    theme(
      axis.title.y = element_blank()
    )  +
    scale_y_continuous(
      breaks = c(-0.2, -0.4, -0.6),
      labels = c("0.2", "0.4", "0.6")#,
      # limits = c(-0.6, -0.1)
      ) +
    coord_cartesian(ylim = c(-0.6, -0.1))
}


bathy_cond <- counterfac_tbl %>% 
  filter(var_in == "mean_bathy") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "mean_bathy") +  
  labs(x = "Mean Bottom\nDepth") +
  geom_text(x = -Inf, y = Inf, label = "a)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)

yday_cond <- counterfac_tbl %>% 
  filter(var_in == "local_day") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "local_day") +  
  labs(x = "Year\nDay") +
  geom_text(x = -Inf, y = Inf, label = "c)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)


slope_cond <- counterfac_tbl %>% 
  filter(var_in == "mean_slope") %>% 
  pull(preds_ci) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "mean_slope") +  
  labs(x = "Mean\nSlope") +
  geom_text(x = -Inf, y = Inf, label = "b)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)


# dist_cond <- counterfac_tbl %>% 
#   filter(var_in == "shore_dist") %>% 
#   pull(preds_ci) %>% 
#   as.data.frame() %>% 
#   plot_foo(data = .,
#            x = "shore_dist") +  
#   labs(x = "Mean Distance\nto Shore (km)") +
#   scale_x_continuous(breaks = c(15000, 35000, 55000),
#                      labels = c("15", "35", "55"),
#                      limits = c(5, 63000),
#                      expand = c(0, 0))


# maturity predictions
mat_pred_in <- rbind(
  new_dat %>% 
    mutate(
      stage_mature = 0
    ),
  new_dat %>% 
    mutate(
      stage_mature = 1
    )
) %>%
  # replace size and lipid with stage-specific averages
  select(-fl, -lipid) %>% 
  left_join(., stage_dat, by = "stage_mature") %>% 
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  )
mat_preds <- pred_foo(mat_pred_in, type = "se", se.method = "infjack")

mat_cond <- mat_preds %>% 
  mutate(
    lo = mean + (qnorm(0.025) * se),
    up = mean + (qnorm(0.975) * se),
    stage_f = factor(stage_mature, levels = c(0, 1), 
                     labels = c("immature", "mature"))
  ) %>% 
  ggplot(.) +
  geom_pointrange(
    aes(x = stage_f, y = mean, ymin = lo, ymax = up)) +
  labs(x = "Maturity Stage +\nLength + Lipid") +
  ggsidekick::theme_sleek() +
  geom_text(aes(x = -Inf, y = Inf, label = "d)"), hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE) +
  theme(
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(breaks = c(-0.2, -0.4, -0.6),
                     labels = c("0.2", "0.4", "0.6"),
                     limits = c(-0.6, -0.15)
  )


# day-night predictions
dn_pred_in <- rbind(
  new_dat %>% 
    mutate(
      day_night_night = 0
    ),
  new_dat %>% 
    mutate(
      day_night_night = 1
    )
) %>%
  mutate(
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  )
dn_preds <- pred_foo(dn_pred_in, type = "se", se.method = "infjack")

dn_cond <- dn_preds %>% 
  mutate(
    lo = mean + (qnorm(0.025) * se),
    up = mean + (qnorm(0.975) * se),
    dn_f = factor(day_night_night, levels = c(0, 1), 
                     labels = c("day", "night"))
  ) %>% 
  ggplot(.) +
  geom_pointrange(
    aes(x = dn_f, y = mean, ymin = lo, ymax = up)) +
  labs(x = "Diel\nCycle") +
  ggsidekick::theme_sleek() +
  geom_text(x = -Inf, y = Inf, label = "e)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE) +
  theme(
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(breaks = c(-0.2, -0.4, -0.6),
                     labels = c("0.2", "0.4", "0.6"),
                     limits = c(-0.6, -0.15)
  )


panel1 <- cowplot::plot_grid(
  bathy_cond,
  slope_cond,
  yday_cond,
  nrow = 1
)
panel2 <- cowplot::plot_grid(
  mat_cond,
  dn_cond,
  ncol = 1
)
panel3 <- cowplot::plot_grid(
  panel2,
  rel_latent,
  ncol = 2,
  rel_widths = c(0.75, 1.5)
)
pp <- cowplot::plot_grid(
  panel1,
  panel3,
  nrow = 2,
  rel_heights = c(0.5, 1)
) 

png(here::here("figs", "ms_figs_rel", "counterfac_effects.png"),
    height = 190, width = 170, units = "mm", res = 300)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    pp,
    left = grid::textGrob("Predicted Bathymetric Depth Ratio", rot = 90)
  )
)
dev.off()

