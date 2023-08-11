## Fit Random Forest with Variance Estimates 

#Model comparison (depth_caret_comparison.R) indicates top model is random 
#forest with moderate number of trees (<200) and fit to untransformed depth 
#data. Alternative to depth_quantreg_rf to generate uncertainty estimates
# See http://shftan.github.io/surfin/demo.html for details
# Nov. 10, 2022
## NOTE: surfin estimates seem to diverge from other RF packages (ranger,
# quantregforest, randomForest); DO NOT USE


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
saveRDS(fit2, here::here("data", "model_fits", "surfin_depth_rf.rds"))

rf_fit <- randomForest(train_depth_baked, train_depth$depth, keep.inbag = TRUE, 
                       ntree = 1000) 
saveRDS(rf_fit, here::here("data", "model_fits", "depth_rf.rds"))


# oob predictive performance
plot(fit2$predicted ~ train_depth$depth)
plot(rf_fit$predicted ~ train_depth$depth)

# poorly correlated with one another
plot(rf_fit$predicted ~ fit2$predicted)
plot(fit$predicted ~ fit2$predicted)


# compare ranger and qunatile regression
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block))
ranger_rf <- ranger(depth ~ ., data = train_depth_baked, quantreg = TRUE, 
                    importance = "permutation")
qr_rf <- readRDS(here::here("data", "model_fits", "depth_quantrf.rds"))


qr_preds <- predict(qr_rf, quantiles = c(0.1, 0.5, 0.9),
                     newdata = train_depth_baked[1:100, ], all = TRUE)
ranger_preds <- predict(ranger_rf, data = train_depth_baked[1:100, ], 
                        quantiles = c(0.1, 0.5, 0.9), type = "quantiles")



