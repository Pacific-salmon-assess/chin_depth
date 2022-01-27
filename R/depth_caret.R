### Depth Models GBM via Caret
## Jan. 25, 2021


library(plyr)
library(tidyverse)
library(caret)
library(recipes)
library(gbm)


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
  filter(stage == "immature") %>% 
  droplevels()

depth_dat <- depth_imm %>% 
  mutate(logit_rel_depth = qlogis(rel_depth)) %>% 
  select(logit_rel_depth, latitude, longitude,
         hour, det_day, mean_bathy, mean_slope, shore_dist,
         vemco_code) 

# check no infinite values in transformed relative depth
nrow(depth_dat[is.infinite(depth_dat$logit_rel_depth), ])
depth_dat <- depth_dat[is.finite(depth_dat$logit_rel_depth), ]
hist(depth_dat$logit_rel_depth)


# subset into training/testing based on 
set.seed(123)
test_tags <- sample(depth_imm$vemco_code, size = 2, replace = F)
train_depth <- depth_dat %>% 
  filter(!vemco_code %in% test_tags) %>%
  select(-vemco_code)
test_depth <- depth_dat %>% 
  filter(vemco_code %in% test_tags) %>%
  select(-vemco_code)

ggplot(test_depth) +
  geom_point(aes(x = hour, y = logit_rel_depth), alpha = 0.5) +
  facet_wrap(~region_f)

  
## CARET PRE-PROCESSING --------------------------------------------------------

#check correlations with additional bathy variables
corr <- cor(test_depth %>%  select(hour:shore_dist))
ggcorrplot::ggcorrplot(corr)
# shore distance and mean bathy just under 0.75, leave for now


depth_recipe <- recipe(logit_rel_depth ~ ., data = train_depth) %>% 
  step_nzv(all_predictors()) %>% 
  #consider adding PCA for bathymetric features
  #step_pca(contains("VSA"), prefix = "surf_area_",  threshold = .95) %>% 
  step_dummy(all_predictors(), -all_numeric()) #%>% 
  # step_center(all_predictors()) %>%
  # step_scale(all_predictors())

# check recipe
prep(depth_recipe) %>% 
  bake(., new_data = train_depth) %>% 
  glimpse()


# # random splitting 
# # stratified splitting of data
# set.seed(998)
# folds <- groupKFold(group = depth_imm$vemco_code, k = 5)
# training <- Sonar[ inTraining,]
# testing  <- Sonar[-inTraining,]


## FIT CARET -------------------------------------------------------------------

# parallelize
ncores <- parallel::detectCores() - 2
doMC::registerDoMC(ncores)


# note does not account for different tags among folds
depth_ctrl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)


# boosted gradient model 
# adjust grid space for hyperparameter tuning
gbmGrid <-  expand.grid(interaction.depth = c(5, 9),
                        n.trees = (10:25)*50,
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

depth_gbm <- train(depth_recipe, train_depth,
                    method = "gbm", 
                    metric = "RMSE",
                    maximize = FALSE,
                    # tuneLength = 10,
                    trControl = depth_ctrl,
                    tuneGrid = gbmGrid)

trellis.par.set(caretTheme())
plot(depth_gbm)  


# random forest model
depth_rf <- train(depth_recipe, train_depth,
                   method = "ranger", 
                   metric = "RMSE",
                   maximize = FALSE,
                   trControl = depth_ctrl,
                   tuneLength = 8,
                   num.trees = 1000)
plot(depth_rf)

## compare
bwplot(resamples(
  list("GBM" = depth_gbm, 
       "RF" = depth_rf)),
       metric = "RMSE")


# predictions
pred_foo <- function(mod, dat = test_depth) {
  preds <- predict(mod, newdata = dat)
  dat$logit_preds <- preds
  
  par(mfrow = c(2, 1))
  plot(logit_preds ~ logit_rel_depth, data = dat)
  abline(0, 1, col = "red")
  plot(plogis(logit_preds) ~ plogis(logit_rel_depth), data = dat)
  abline(0, 1, col = "red")
}

pred_foo(depth_gbm, dat = test_depth)
pred_foo(depth_rf, dat = train_depth)


# variable importance
gbm_imp <- varImp(depth_gbm, scale = F)
plot(gbm_imp)
rf_imp <- varImp(depth_rf)

# evaluate patterns
explainer_gbm <- explain(
  depth_gbm,
  data = select(train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist,
                latitude, longitude),
  y = train_depth$logit_rel_depth,
  label = "gbm"
)
explainer_rf <- explain(
  depth_gbm,
  data = select(train_depth, hour, det_day, mean_bathy, mean_slope, shore_dist,
                latitude, longitude),
  y = train_depth$logit_rel_depth,
  label = "random forest"
)

# function to visualize profiles
make_pdp <- function(param) {
  pdp_gbm <- model_profile(explainer_gbm, N = 400, variables = param)
  pdp_rf <- model_profile(explainer_rf, N = 400, variables = param)
  rbind(as_tibble(pdp_gbm$agr_profiles), as_tibble(pdp_rf$agr_profiles)) %>%
    ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
    geom_line() + xlab(param)
}

make_pdp("hour")
make_pdp("det_day")
make_pdp("mean_bathy")
make_pdp("mean_slope")
make_pdp("shore_dist")
make_pdp("latitude")
make_pdp("longitude")


dum <- model_profile(explainer_gbm, N = 400, variables = "hour")
dum2 <- model_profile(explainer_rf, N = 400, variables = "hour")

rbind(as_tibble(dum$agr_profiles), as_tibble(dum2$agr_profiles)) %>%
  ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
  geom_line() + xlab("hour")

library(DALEX)
library(DALEXtra)

explainer_rf <- explain_tidymodels(
  rf_fit,
  data = select(dat_train, hour_c, day_c, max_bathy_c, latitude, longitude),
  y = dat_train$logit_rel_depth,
  label = "random forest"
)

make_pdp <- function(param) {
  pdp_rf <- model_profile(explainer_rf, N = 400, variables = param)
  as_tibble(pdp_rf$agr_profiles) %>%
    ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
    geom_line() + xlab(param)
}

make_pdp("hour_c")
make_pdp("day_c")
make_pdp("max_bathy_c")
make_pdp("latitude")
make_pdp("longitude")