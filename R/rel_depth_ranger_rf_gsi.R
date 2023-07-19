## Fit Regression Random Forest w/ Percentile Intervals

## Equivalent to rel_depth_ranger_rf.R except adds population identity
# as covariates

library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)


depth_dat_raw1 <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  # approximately 6k detections have no available ROMS data; exclude for now
  filter(!is.na(roms_temp)) %>% 
  #bin aggregates
  mutate(
    agg = case_when(
      grepl("Fraser", agg) ~ "Fraser",
      grepl("Col", agg) ~ "Col",
      TRUE ~ agg
    )
  )


# remove 2022 tag releases (~6k dets) for training model
depth_dat_raw <- depth_dat_raw1 %>% 
  filter(!grepl("2022", vemco_code))


# number of detections per tag
depth_dat_raw1 %>% 
  group_by(vemco_code) %>% 
  tally() %>% 
  pull(n) %>% 
  range()

# calculate timespan overwhich detections provided
timespan <- depth_dat_raw1 %>% 
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
    agg,
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
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    agg)


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
rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))
top_mod <- rf_list[[2]]$top_model

ranger_rf <- ranger::ranger(
  depth ~ .,
  data = train_depth_baked,
  #hyperpars based on values from top model which is not saved on all locals
  num.trees = 2000,
  mtry = 11,
  keep.inbag = TRUE,
  quantreg = TRUE,
  importance = "permutation"
)

saveRDS(ranger_rf,
        here::here("data", "model_fits", "relative_rf_ranger_gsi.rds"))
ranger_rf <- readRDS(here::here("data", "model_fits", "relative_rf_ranger.rds"))


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
    x = "Observed Depth", y = "Predicted Mean Depth"
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
      var %in% c("stage_mature", "fl", "lipid") ~ "biological",
      grepl("agg", var) ~ "stock",
      TRUE ~ "dynamic"
    )
  ) %>%
  arrange(-val)
imp_dat$var_f = factor(
  imp_dat$var, 
  labels = c("Year Day 1", "Bottom Depth", "UTM X", "UTM Y", "Lunar Cycle", 
             "Temperature", "Zooplankton", "Maturity", "Bottom Slope",
             "Year Day 2", "Oxygen", "Fork Length", "Shore Distance", 
             "Thermocline Depth", "H Current 1", "Lipid Content", "H Current 2",
             "Day/Night", "Vertical Current", "Columbia", "Puget",
             "Fraser", "WA/OR", "ECVI")
)

imp_plot <- ggplot(imp_dat, aes(x = var_f, y = val)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(here::here("figs", "ms_figs_rel", "importance_quantreg_gsi.png"),
    height = 4, width = 6, units = "in", res = 250)
imp_plot
dev.off()

