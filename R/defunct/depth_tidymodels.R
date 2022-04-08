# Most columns in DF are self-explanatory, but depth is tag depth and bathy
# (bathymetry) represents the mean, max, or SD within 800 meters of the moored
# receiver (approximate detection radius)
# Lat/longs refer to location ofreceiver mooring where detection occurred
# Maturation stage is inferred based on location/timing of final detection or
# body size and is relevant because immature fish are at large longer and aren't
# undergoing a directed migration
# Depth estimates are means within a 60 minute timebin; this is partially to
# make the number of datapoints tractable and to help deal with autocorrelation
# which is significant even with this pooling; bins were assigned relative to
# the first detection in the first year
# Region is a fairly arbitrary assignment, hence why I'd like to ultimately move
# to a spatially-explicit model

library(tidymodels)

all_cores <- parallel::detectCores(logical = FALSE)
library(doFuture)
registerDoFuture()
cl <- parallel::makeCluster(all_cores)
plan(cluster, workers = cl)

depth_dat <- readRDS(here::here("data", "depth_dat_60min.RDS")) %>% 
  mutate(
    day_c = as.numeric(scale(det_day, center = TRUE, scale = FALSE)),
    hour_c = as.numeric(scale(hour, center = TRUE, scale = FALSE)),
    max_bathy_c = as.numeric(scale(max_bathy, center = TRUE, scale = FALSE))
  )

ggplot(depth_dat, aes(longitude, latitude, colour = pos_depth)) +
  geom_point()

# try turning 0-1 variable to -Inf, Inf:
depth_dat$logit_rel_depth <- qlogis(depth_dat$rel_depth)
nrow(depth_dat[is.infinite(depth_dat$logit_rel_depth), ])
depth_dat <- depth_dat[is.finite(depth_dat$logit_rel_depth), ]
hist(depth_dat$logit_rel_depth)

# split for final testing later:
set.seed(92819)
dat_split <- initial_split(depth_dat, prop = 0.8, strata = stage)
dat_train <- training(dat_split)
dat_test <- testing(dat_split)

# TODO: split based on individual fish! And possibly based on spatial block.

rec <- recipe(logit_rel_depth ~ hour_c + day_c +
    max_bathy_c + latitude + longitude, data = dat_train)

# rec <- recipe(logit_rel_depth ~ latitude + longitude, data = dat_train)

rf_model <- rand_forest(mode = "regression", mtry = 5, trees = 1000, min_n = 10) %>%
  set_engine("ranger")

rf_wflow <-
  workflow() %>%
  add_model(rf_model) %>%
  add_recipe(rec)

rf_fit <- fit(rf_wflow, dat_train)
p1 <- predict(rf_fit, new_data = dat_train)
plot(dat_train$logit_rel_depth, p1$.pred);abline(0, 1)
plot(plogis(dat_train$logit_rel_depth), plogis(p1$.pred))

p1 <- predict(rf_fit, new_data = dat_test)
plot(dat_test$logit_rel_depth, p1$.pred);abline(0, 1)
plot(plogis(dat_test$logit_rel_depth), plogis(p1$.pred))

hyper_par_grid <- expand_grid(
  mtry = c(1, 3, 5),
  trees = c(500, 1000, 2000),
  min_n = c(5, 10, 20)
)
nrow(hyper_par_grid)

# TODO: split based on fish and possibly spatial block:
cell_folds <- vfold_cv(dat_train, v = 5L, repeats = 1L)

tune_spec <-
  rand_forest(mode = "regression",
    mtry = tune(),
    trees = tune(),
    min_n = tune()
  ) %>%
  set_engine("ranger") %>%
  set_mode("regression")

rf_wflow_tune <-
  workflow() %>%
  add_model(tune_spec) %>%
  add_recipe(rec)

tune_res <-
  rf_wflow_tune %>%
  tune_grid(
    resamples = cell_folds,
    grid = hyper_par_grid,
    control = control_grid(save_pred = TRUE)
  )

plot_tuning <- function(metric) {
  tune_res %>%
    collect_metrics() %>%
    filter(.metric == metric) %>%
    mutate(min_n = as.factor(min_n)) %>%
    ggplot(aes(mtry, mean, color = min_n)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ trees, scales = "free") +
    scale_x_log10(labels = scales::label_number()) +
    scale_color_viridis_d(option = "plasma", begin = .9, end = 0)
}
plot_tuning("rsq")
plot_tuning("rmse")

best_rmse <- select_best(tune_res, "rmse")
best_rmse

rf_final_wf <- finalize_workflow(rf_wflow_tune, best_rmse)
rf_final_wf

rf_fit <- fit(rf_final_wf, dat_train)

p1 <- predict(rf_fit, new_data = dat_train)
plot(dat_train$logit_rel_depth, p1$.pred);abline(0, 1)
plot(plogis(dat_train$logit_rel_depth), plogis(p1$.pred))

p1 <- predict(rf_fit, new_data = dat_test)
plot(dat_test$logit_rel_depth, p1$.pred);abline(0, 1)
plot(plogis(dat_test$logit_rel_depth), plogis(p1$.pred))

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

pdp_rf <- model_profile(explainer_rf, N = 400, variables = c("longitude", "latitude"))
# pdp <- as_tibble(pdp_rf$agr_profiles)
# pdp %>% select(-`_label_`, -`_ids_`) %>% filter(`_vname_` == "longitude")

