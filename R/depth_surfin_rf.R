## Fit Random Forest with Variance Estimates 

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Alternative to depth_quantreg_rf to generate uncertainty estimates
# See http://shftan.github.io/surfin/demo.html for details
# Nov. 10, 2022


library(plyr)
library(tidyverse)
library(recipes)
library(surfin)
library(randomForest)


depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>% 
  mutate(stage = as.factor(stage)) %>% 
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
    det_day = local_day, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    ind_block
  ) 

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block)) %>% 
  #impute missing ROMS values (now done externally)
  # step_impute_knn(all_predictors(), neighbors = 3) %>%
  # step_nzv(all_predictors()) %>% 
  step_dummy(all_predictors(), -all_numeric())


train_folds <- caret::groupKFold(
  train_depth$ind_block,
  k = length(unique(train_depth$ind_block))
)
depth_ctrl <-   caret::trainControl(
  method="repeatedcv",
  index = train_folds
)


# refit top model using surfin
# use 1000 trees as in depth_quantreg and recommendation of 10-50 for B 

# apply recipe to dataframe to interpolate values
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block)) %>% 
  select(-depth)


fit <- forest(train_depth_baked, train_depth$depth, 
              var.type="ustat", B = 25, ntree = 1000)
fit2 <- forest(train_depth_baked, train_depth$depth, 
              var.type="infjack", B = 25, ntree = 1000)
rf_fit <- randomForest(train_depth_baked, train_depth$depth, keep.inbag = TRUE, 
                       ntree = 1000) 
saveRDS(rf_fit, here::here("data", "model_fits", "depth_rf.rds"))


# oob predictive performance
plot(fit$predicted ~ train_depth$depth)
plot(rf_fit$predicted ~ train_depth$depth)

# compare to random forest predictions and quantile random forest predictions
qr_fit <- readRDS(here::here("data", "model_fits", "depth_quantrf.rds"))
rf_fit <- readRDS(here::here("data", "model_fits", "depth_rf.rds"))

plot(fit$predicted ~ qr_fit$predicted)
plot(rf_fit$predicted ~ fit$predicted)
plot(fit$predicted ~ fit2$predicted)




## variable importance

imp_dat <- as.data.frame(rf_fit$importance, row.names = FALSE) %>%
  janitor::clean_names() %>% 
  mutate(
    var = rownames(rf_fit$importance) %>% 
      fct_reorder(., -inc_node_purity ),
    category = case_when(
      var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y", 
                 "mean_slope") ~ "spatial",
      var %in% c("det_day", "day_night_night", "moon_illuminated") ~ "temporal",
      var %in% c("stage_mature", "fl", "mean_log_e") ~ "biological",
      TRUE ~ "dynamic"
    )
  ) %>% 
  arrange(-inc_node_purity) 
imp_dat$var_f = factor(
  imp_dat$var, 
  labels = c("Bottom Depth", "Fork Length", "Maturity", "Year Day", "UTM Y",
             "Condition", "UTM X", "Moon Phase", "Bottom Slope", "Zooplankton", 
             "Shore Distance", "Temperature", "Oxygen",
             "Thermocline Depth", "Day/Night", "H Current 1", "H Current 2", 
             "Vertical Current")
)

imp_plot <- ggplot(imp_dat, aes(x = var_f, y = percent_inc_mse)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# scale results in whiskers being not visible
# geom_pointrange(aes(ymin = lo, ymax = up), shape = 21)

png(here::here("figs", "depth_ml", "importance_quantreg.png"),
    height = 4, width = 6, units = "in", res = 250)
imp_plot
dev.off()