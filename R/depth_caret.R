### Depth Models GBM via Caret
## Jan. 25, 2021


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)

time_foo <- function(x) {
  lubridate::hour(x) + (lubridate::minute(x) / 60) + 
    (lubridate::second(x) / 3600) 
}

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_60min.RDS")) 


# ggplot(depth_dat_raw) +
#   geom_point(aes(x = mean_bathy_c, y = rel_depth))
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = day_c, y = rel_depth), alpha = 0.4) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_point(aes(x = hour_c, y = rel_depth), alpha = 0.5) +
#   facet_grid(stage~region_f)
# ggplot(depth_dat_raw) +
#   geom_boxplot(aes(x = day_night, y = rel_depth)) +
#   facet_grid(stage~region_f)


depth_imm <- depth_dat_raw %>%  
  filter(stage == "immature") 
depth_dat <- depth_imm %>% 
  mutate(hour = time_foo(date_time_local),
         logit_rel_depth = qlogis(rel_depth)) %>% 
  select(logit_rel_depth,# region_f, 
         hour, det_day, mean_bathy = pos_mean_bathy)

# subset based into training/testing
test_tags <- sample(depth_imm$vemco_code, size = 5, replace = F)
test_depth <- depth_dat[depth_imm$vemco_code %in% test_tags, ]
train_depth <- depth_dat[!depth_imm$vemco_code %in% test_tags, ]

  
## CARET PRE-PROCESSING --------------------------------------------------------

depth_recipe <- recipe(logit_rel_depth ~ ., data = depth_dat) %>% 
  step_nzv(all_predictors()) %>% 
  #consider adding PCA for bathymetric features
  #step_pca(contains("VSA"), prefix = "surf_area_",  threshold = .95) %>% 
  step_dummy(all_nominal(), -all_outcomes()) %>% 
  step_center(all_predictors()) %>%
  step_scale(all_predictors())


# # random splitting 
# # stratified splitting of data
# set.seed(998)
# folds <- groupKFold(group = depth_imm$vemco_code, k = 5)
# training <- Sonar[ inTraining,]
# testing  <- Sonar[-inTraining,]


## FIT CARET -------------------------------------------------------------------


# note does not account for different tags among folds
depth_ctrl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

depth_gbm <- train(depth_recipe, train_depth,
                 method = "gbm", 
                 metric = "RMSE",
                 maximize = FALSE,
                 # tuneLength = 10,
                 trControl = depth_ctrl)

# adjust grid space for hyperparameter tuning
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                        n.trees = (3:30)*50,
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

depth_gbm2 <- train(depth_recipe, train_depth,
                   method = "gbm", 
                   metric = "RMSE",
                   maximize = FALSE,
                   # tuneLength = 10,
                   trControl = depth_ctrl,
                   tuneGrid = gbmGrid)

trellis.par.set(caretTheme())
plot(depth_gbm2)  
