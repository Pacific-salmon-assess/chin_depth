## Compare original RF, weighted RF and hierarchical GAM

# Fit approximately equivalent model to ML alternative based on feedback from
# reviewers
# Fit model, test for autocorrelation in residuals, estimate conditional effects
# and evaluate hold out performance on equivalent data
# NOTE: GAMM with explicit autocorrelation structure failed to converge

library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)
library(mgcv)
library(ape)


# parallelize based on operating system
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}



depth_dat_raw1 <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  # approximately 7k detections have no available ROMS data; exclude 
  filter(!is.na(roms_temp)) %>% 
  # define time in hours for binning (all roms variables constant)
  mutate(
    date_time_hour = format(date_time_utm, format = "%Y-%m-%d %H"),
    day_night_dummy = ifelse(day_night == "night", 1, 0)
  ) %>% 
  group_by(vemco_code) %>% 
  #weight based on number of observations
  mutate(n_dets = n(),
         wt = 1 / sqrt(n_dets),
         wt2 = 1 / n_dets,
         # add time steps
         start_time = min(date_time_utm),
         end_time = max(date_time_utm),
         timestamp = difftime(start_time, date_time_utm, units = "mins"),
         timestamp = -1 * round(as.numeric(timestamp)),
         start_time_hour = min(date_time_hour),
         end_time_hour = max(date_time_hour),
         timestamp_hour = difftime(start_time_hour, 
                                   date_time_hour, units = "hours"),
         timestamp_hour = -1 * round(as.numeric(timestamp_hour))) %>% 
  ungroup() 


det_rates <- depth_dat_raw1 %>% 
  group_by(vemco_code) %>% 
  summarize(
    n_dets = n(),
    n_days = length(unique(local_day)),
    n_rec = length(unique(receiver_name))
  ) %>% 
  mutate(
    n_det_bin = cut(n_dets, breaks = c(0, 10, 100, 1000, 10000)),
    n_day_bin = cut(n_days, breaks = c(0, 10, 30, 100, 500)),
    n_rec_bin = cut(n_rec, breaks = c(0, 5, 20, 50, 100))
  )

det_rates %>% 
  group_by(n_det_bin) %>% 
  tally()
det_rates %>% 
  group_by(n_day_bin) %>% 
  tally()
det_rates %>% 
  group_by(n_rec_bin) %>% 
  tally()


# remove 2022 tag releases (~6k dets) for training model
depth_dat_raw <- depth_dat_raw1 %>% 
  filter(!grepl("2022", vemco_code))


# add individual block
set.seed(1234)
ind_folds <- data.frame(
  vemco_code = unique(depth_dat_raw$vemco_code),
  ind_block =  sample.int(
    8, length(unique(depth_dat_raw$vemco_code)), replace = T
  ) %>% 
    as.factor()
)

# binned data
depth_dat_bin <- depth_dat_raw %>%
  group_by(vemco_code,
           timestamp_hour, 
           fl, lipid, med_stage, utm_x, utm_y,
           day_night_dummy, local_day, mean_bathy, mean_slope, shore_dist,
           u, v, w, roms_temp, zoo, oxygen, thermo_depth) %>%
  summarize(
    moon_illuminated = mean(moon_illuminated),
    mean_rel_depth = mean(rel_depth),
    sd_rel_depth = sd(rel_depth),
    .groups = "drop"
  ) %>%
  left_join(., ind_folds, by = "vemco_code") %>%
  ungroup()

depth_dat <- depth_dat_raw %>% 
  left_join(., ind_folds, by = "vemco_code")


# split by individual blocking (identical to random forest)
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()

train_bin <- depth_dat_bin %>% 
  filter(!ind_block == "5") %>% 
  droplevels() 

test_depth_22 <- depth_dat_raw1 %>% filter(grepl("2022", vemco_code))


## FIT GAMS --------------------------------------------------------------------


fit_full <- gam(
  rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k=c(10, 10)) +
    s(mean_bathy, k = 3) + 
    s(local_day, bs = "cc", k = 5) + 
    s(med_stage, k = 4) +
    day_night_dummy + 
    s(vemco_code, bs = "re"),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

fit_sub <- gam(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k = c(10, 10)) +
    s(mean_bathy, k = 5) + s(mean_slope, k = 5) + 
    s(local_day, bs = "cc", k = 5) + 
    s(med_stage, k = 4) +
    day_night_dummy + thermo_depth +
    fl + lipid + s(zoo, k = 5) + shore_dist + moon_illuminated + u + v + w + 
    roms_temp +
    s(vemco_code, bs = "re"),
  data = train_bin,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

fit_full_trim <- gam(
  rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k=c(10, 10)) +
    s(local_day, bs = "cc", k = 5) + 
    s(med_stage, k = 4) +
    day_night_dummy + 
    s(vemco_code, bs = "re"),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

fit_sub_trim <- gam(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k=c(10, 10)) +
    s(local_day, bs = "cc", k = 5) + 
    s(med_stage, k = 4) +
    day_night_dummy + 
    s(vemco_code, bs = "re"),
  
  data = train_bin,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

concurvity(fit_sub)
concurvity(fit_sub_trim)
concurvity(fit_full)
concurvity(fit_full_trim)


saveRDS(fit_sub, here::here("data", "model_fits", "gam_fits_sub.rds"))
saveRDS(fit_full, here::here("data", "model_fits", "gam_fits_full.rds"))
saveRDS(fit_sub_trim, here::here("data", "model_fits", "gam_fits_sub_trim.rds"))
saveRDS(fit_full_trim, here::here("data", "model_fits", "gam_fits_full_trim.rds"))

fit_full <- readRDS(here::here("data", "model_fits", "gam_fits_full.rds"))
fit_sub <- readRDS(here::here("data", "model_fits", "gam_fits_sub.rds"))


## GAM CONVERGENCE -------------------------------------------------------------

hist(resid(fit_full))
hist(resid(fit_sub))

# DOESN'T WORK FOR BETA

# fit_list <- list(fit_full, fit_sub)
# 
# sims_list <- purrr::map2(
#   fit_list, list(train_depth, train_bin),
#   ~ simulate(.x, data = .y, nsim = 50)
# )
# fixed_pred_list <- purrr::map(
#   fit_list, 
#   ~ predict(.x, type = "response")
# )
# qq_list <- purrr::pmap(
#   list(sims_list, fixed_pred_list), 
#   function(y, z) {
#     dharma_res <- DHARMa::createDHARMa(
#       simulatedResponse = y,
#       observedResponse = train_depth$mean_rel_depth,
#       fittedPredictedResponse = z
#     )
#     plot(dharma_res)
#   }
# )
# # looks good


## FIT RANGERS -----------------------------------------------------------------

train_depth_ml <- train_depth %>% 
  select(rel_depth, fl, lipid, med_stage, utm_y, utm_x, mean_bathy, mean_slope,
         shore_dist, u, v, w, roms_temp, zoo, oxygen, thermo_depth, det_dayx,
         det_dayy, moon_illuminated, day_night_dummy)

# hyperpars based on top model from depth_caret_comparison_weighted.R
ranger_rf <- ranger::ranger(
  rel_depth ~ .,
  data = train_depth_ml,
  num.trees =  1000, #$rf_list$rel_depth$top_model$num.trees,
  mtry = 17, #rf_list$rel_depth$top_model$mtry,
  # keep.inbag = TRUE for quantile predictions
  keep.inbag = TRUE,
  quantreg = TRUE,
  importance = "permutation",
  splitrule = "extratrees"
)

ranger_rf_w <- ranger::ranger(
  rel_depth ~ .,
  data = train_depth_ml,
  num.trees =  2500,#rf_weighted_list$rel_depth$top_model$num.trees,
  mtry = 9, #rf_weighted_list$rel_depth$top_model$mtry,
  case.weights = train_depth$wt,
  # keep.inbag = TRUE for quantile predictions
  keep.inbag = TRUE,
  # quantreg = TRUE,
  importance = "permutation",
  splitrule = "extratrees"
)

ranger_rf_w2 <- ranger::ranger(
  rel_depth ~ .,
  data = train_depth_ml,
  num.trees =  2500,#rf_weighted_list$rel_depth$top_model$num.trees,
  mtry = 11, #rf_weighted_list$rel_depth$top_model$mtry,
  case.weights = train_depth$wt2,
  # keep.inbag = TRUE for quantile predictions
  keep.inbag = TRUE,
  # quantreg = TRUE,
  importance = "permutation",
  splitrule = "extratrees"
)

ranger_list <- list(ranger_rf,
                    ranger_rf_w,
                    ranger_rf_w2)
saveRDS(ranger_list, here::here("data", "model_fits", "supp_ranger_fits.rds"))


rf_preds <- predict(ranger_rf, data = train_depth_ml)
rf_w_preds <- predict(ranger_rf_w, data = train_depth_ml)
rf_w2_preds <- predict(ranger_rf_w2, data = train_depth_ml)


## CHECK AUTOCORRELATION -------------------------------------------------------

gam_preds <- predict(fit_full, type = "response") %>% as.numeric()

# calculate pearson residuals from each model
p_resid <- function (pred, obs) {
  (obs - pred) / sqrt(pred * (1 - pred))
}

resid_dat <- train_depth %>% 
  mutate(
    resid_gam = p_resid(gam_preds, rel_depth),
    resid_rf = p_resid(rf_preds$predictions, rel_depth),
    resid_rf_w = p_resid(rf_w_preds$predictions, rel_depth),
    resid_rf_w2 = p_resid(rf_w2_preds$predictions, rel_depth)
  ) %>% 
  arrange(
    date_time_utm
  )


## temporal autocorrelation

# ID tags w/ mult dets
kept_tags <- resid_dat %>% 
  filter(n_dets > 1) %>% 
  pull(vemco_code) %>% 
  unique()
kept_tags_bin <- train_bin %>% 
  group_by(vemco_code) %>% 
  mutate(n_dets = n()) %>% 
  filter(n_dets > 1) %>% 
  pull(vemco_code) %>% 
  unique()

train_bin$resid_gam_bin <- residuals(fit_sub, type = "pearson")

# infill function 
infill_foo <- function(dum, timestep = "timestamp") {
  dd <- dum %>% 
    rename(tt = {{ timestep }})
  expand.grid(
    timestamp = seq(min(dd$tt), max(dd$tt), by = 1)
  ) %>% 
    left_join(., 
              dd %>%
                select(vemco_code, timestamp = tt, starts_with("resid")),
              by = "timestamp")
}

# calculate acf for each tag
resid_ts <- resid_dat %>% 
  droplevels() %>% 
  split(., .$vemco_code) %>% 
  furrr::future_map(
    ., ~ infill_foo(., timestep = "timestamp")
  )
resid_ts_bin <- train_bin %>% 
  droplevels() %>% 
  split(., .$vemco_code) %>% 
  furrr::future_map(
    ., ~ infill_foo(., timestep = "timestamp_hour")
  )

ar1_est <- resid_ts %>% 
  bind_rows() %>% 
  filter(vemco_code %in% kept_tags) %>% 
  group_by(vemco_code) %>% 
  summarize(
    gam = acf(resid_gam, plot = FALSE, na.action = na.pass)$acf[2, 1, 1],
    rf = acf(resid_rf, plot = FALSE, na.action = na.pass)$acf[2, 1, 1],
    rf_w = acf(resid_rf_w, plot = FALSE, na.action = na.pass)$acf[2, 1, 1],
    rf_w2 = acf(resid_rf_w2, plot = FALSE, na.action = na.pass)$acf[2, 1, 1]
  ) 
ar1_est_bin <- resid_ts_bin %>% 
  bind_rows() %>% 
  filter(vemco_code %in% kept_tags_bin) %>% 
  group_by(vemco_code) %>% 
  summarize(
    ar1 = acf(resid_gam_bin, plot = FALSE, na.action = na.pass)$acf[2, 1, 1]
  ) %>% 
  mutate(model = "gam_bin")

ar1_est2 <- ar1_est %>% 
  pivot_longer(cols = c(gam, rf, rf_w, rf_w2),
               names_to = "model", values_to = "ar1") %>% 
  rbind(., ar1_est_bin) %>% 
  mutate(
    model_family = ifelse(grepl("gam", model), "gam", "random_forest")
  ) %>% 
  filter(
    !is.na(ar1)
  )

ar1_box <- ggplot(ar1_est2, aes(x = model, y = ar1, fill = model_family)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  labs(x = "Model", y = "Tag-Specific AR-1 Estimate")

png(here::here("figs", "model_comp", "ar1_plot.png"), units = "in",
    height = 3.5, width = 4.5, res = 250)
ar1_box
dev.off()


ar1_est2 %>% 
  group_by(model) %>% 
  summarize(
    med_ar1 = median(ar1),
    mean_ar1 = mean(ar1),
    sd_ar1 = sd(ar1)
  )



# spatial autocorrelation estimates
resid_list <- list(resid_dat$resid_gam, resid_dat$resid_rf, 
                   resid_dat$resid_rf_w, resid_dat$resid_rf_w2)
moran_full <- purrr::map(
  resid_list, ~ moranfast::moranfast(.x, resid_dat$utm_x, resid_dat$utm_y)
)
moran_bin <- moranfast::moranfast(train_bin$resid_gam_bin, train_bin$utm_x, 
                                  train_bin$utm_y)


## OUT OF SAMPLE PREDS ---------------------------------------------------------

rmse_dat <- ar1_est2 %>% 
  select(model, model_family) %>% 
  distinct() %>% 
  arrange(model_family)


gam_full_22 <- predict(
  fit_full, type = "response", newdata = test_depth_22,
  #estimate population level, excluding RIs for tag
  exclude = "s(vemco_code)"
) %>% 
  as.numeric()
gam_bin_22 <- predict(
  fit_sub, type = "response", newdata = test_depth_22,
  #estimate population level, excluding RIs for tag
  exclude = "s(vemco_code)"
) %>% 
  as.numeric()

rf_22 <- predict(
  ranger_rf, data = test_depth_22
)
rf_w_22 <- predict(
  ranger_rf_w, data = test_depth_22
)
rf_w2_22 <- predict(
  ranger_rf_w2, data = test_depth_22
)

pred_list <- list(gam_full_22, gam_bin_22, rf_22$predictions, 
                  rf_w_22$predictions, rf_w2_22$predictions)

rmse_dat$rmse <- purrr::map(
  pred_list, ~ Metrics::rmse(test_depth_22$rel_depth, .x)
) %>% 
  unlist()
rmse_dat$bias <- purrr::map(
  pred_list, ~ Metrics::bias(test_depth_22$rel_depth, .x)
) %>% 
  unlist()

rmse_bar <- ggplot(
  rmse_dat, aes(x = model, y = rmse, fill = model_family)
) +
  geom_bar(stat = "identity") +
  ggsidekick::theme_sleek() +
  labs(x = "Model", y = "RMSE") 

bias_bar <- ggplot(
  rmse_dat, aes(x = model, y = bias, fill = model_family)
) +
  geom_bar(stat = "identity") +
  ggsidekick::theme_sleek() +
  labs(x = "Model", y = "Bias") 

png(here::here("figs", "model_comp", "rmse_bar.png"), units = "in",
    height = 3.5, width = 4.5, res = 250)
rmse_bar
dev.off()

png(here::here("figs", "model_comp", "bias_bar.png"), units = "in",
    height = 3.5, width = 4.5, res = 250)
bias_bar
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


# biological data
bio_dat <- depth_dat_raw %>%
  filter(med_stage %in% c(0, 1)) %>% 
  select(vemco_code, fl, lipid, med_stage) %>%
  distinct()

# stratify predictions by non-spatial covariates
pred_dat1 <- bath_grid %>%
  mutate(
    fl = mean(bio_dat$fl),
    lipid = mean(bio_dat$lipid),
    day_night_dummy = 0.5,
    med_stage = 0.5,
    vemco_code = unique(depth_dat_raw$vemco_code)[1]
  )

stage_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, lipid, med_stage) %>% 
  filter(med_stage %in% c(0, 1)) %>% 
  distinct() %>% 
  group_by(med_stage) %>% 
  dplyr::summarize(
    fl = mean(fl),
    lipid = mean(lipid)
  ) %>% 
  mutate(stage = ifelse(med_stage == "1", "mature", "immature")) %>% 
  rbind(., 
        data.frame(
          med_stage = 0.5,
          fl = mean(bio_dat$fl),
          lipid = mean(bio_dat$lipid),
          stage = "average"
        ))

new_dat <- data.frame(
  utm_x = median(train_depth$utm_x),
  utm_y = median(train_depth$utm_y),
  fl = median(train_depth$fl),
  lipid = median(train_depth$lipid),
  mean_bathy = median(train_depth$mean_bathy),
  local_day = median(train_depth$local_day),
  mean_slope = median(train_depth$mean_slope),
  shore_dist = median(train_depth$shore_dist),
  u = median(train_depth$u),
  v = median(train_depth$v),
  w = median(train_depth$w),
  zoo = median(train_depth$zoo),
  oxygen = median(train_depth$oxygen),
  thermo_depth = median(train_depth$thermo_depth),
  roms_temp = median(train_depth$roms_temp),
  moon_illuminated = median(train_depth$moon_illuminated),
  day_night_dummy = 0.5,
  med_stage = 0.5,
  vemco_code = unique(depth_dat_raw$vemco_code)[1]
) %>% 
  mutate(
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
  ) %>% 
  #duplicate 100 times
  as_tibble() %>% 
  slice(rep(1:n(), each = 100))


# make tibble for different counterfacs
gen_pred_dat <- function(var_in) {
  # necessary to deal with string input
  varname <- ensym(var_in)
  if (grepl("utm", var_in)) {
    group_vals <- bath_grid %>% 
      dplyr::summarize(min_v = min(!!varname),
                       max_v = max(!!varname))
  } else {
    group_vals <- train_depth %>% 
      dplyr::summarize(min_v = min(!!varname),
                       max_v = max(!!varname))
  }
  # change to dataframe
  var_seq <- NULL
  for (i in 1:nrow(group_vals)) {
    var_seq <- c(var_seq, 
                 seq(group_vals$min_v[i], group_vals$max_v[i], 
                     length.out = 100))
  }
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
pred_foo <- function(fit, preds_in, ...) {
  preds1 <- predict(
    object = fit, 
    newdata = preds_in,
    type = "response",
    se.fit = TRUE,
    exclude = "s(vemco_code)",
    ...
  )
  
  preds_in %>%
    mutate(
      mean = as.numeric(preds1$fit),
      se = as.numeric(preds1$se.fit),
      lo = -1 * (mean + (qnorm(0.025) * se)),
      up = -1 * (mean + (qnorm(0.975) * se)),
      mean = -1 * mean
    ) 
}


pred_foo_ranger <- function(fit, preds_in, ...) {
  preds1 <- predict(
    fit, 
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
  var_in = c("mean_bathy", "utm_x", "utm_y", "local_day", "med_stage"
  )) %>%
  mutate(
    pred_dat_in = purrr::map(var_in,
                             gen_pred_dat),
    preds_gam_full = furrr::future_map(pred_dat_in, 
                                ~ pred_foo(fit = fit_full, preds_in = .x)),
    preds_gam_sub = furrr::future_map(pred_dat_in, 
                                ~ pred_foo(fit = fit_sub, preds_in = .x)))
counterfac_tbl$preds_rf <- purrr::map(
  counterfac_tbl$pred_dat_in, 
  ~ pred_foo_ranger(fit = ranger_rf, preds_in = .x, type = "se", 
                    se.method = "infjack")
)
counterfac_tbl$preds_rf_w <- purrr::map(
  counterfac_tbl$pred_dat_in, 
  ~ pred_foo_ranger(fit = ranger_rf_w, preds_in = .x, type = "se", 
                    se.method = "infjack")
)
counterfac_tbl$preds_rf_w2 <- purrr::map(
  counterfac_tbl$pred_dat_in, 
  ~ pred_foo_ranger(fit = ranger_rf_w2, preds_in = .x, type = "se", 
                    se.method = "infjack")
)

saveRDS(counterfac_tbl,
        here::here("data", "model_fits", "supplementary_counterfac_preds.rds"))

  
bind_foo <- function (cov_in) {
  dum <- counterfac_tbl %>% filter(var_in == {{ cov_in }})
  gam_full_dat <- dum$preds_gam_full[[1]] %>% mutate(model = "gam")
  gam_bin_dat <- dum$preds_gam_sub[[1]] %>% mutate(model = "gam_bin")
  rf_dat <- dum$preds_rf[[1]] %>% mutate(model = "rf")
  rf_w_dat <- dum$preds_rf_w[[1]] %>% mutate(model = "rf_w")
  rf_w2_dat <- dum$preds_rf_w2[[1]] %>% mutate(model = "rf_w2")
  
  out1 <- rbind(gam_full_dat, gam_bin_dat) %>% 
    select(model, {{ cov_in }}, mean, lo, up) %>% 
    mutate(model_family = "gam")
  out2 <- list(rf_dat, rf_w_dat, rf_w2_dat) %>%
    bind_rows() %>% 
    select(model, {{ cov_in }}, mean, lo, up) %>% 
    mutate(model_family = "random_forest")
  
  rbind(out1, out2)
}

plot_foo <- function (data, ...) {
  ggplot(data, mapping = aes(!!!ensyms(...))) +
    geom_line(aes(y = mean, colour = model_family)) +
    geom_ribbon(aes(ymin = lo, ymax = up, fill = model_family), alpha = 0.4) +
    ggsidekick::theme_sleek() +
    scale_x_continuous(expand = c(0, 0))+
    theme(
      axis.title.y = element_blank()
    ) +
    coord_cartesian(ylim = c(-1, -0))
}


png(here::here("figs", "model_comp", "bathy_cf_pred.png"), units = "in",
    height = 2.5, width = 6, res = 250)
bind_foo("mean_bathy") %>% 
  plot_foo(data = ., x = "mean_bathy") +  
  labs(x = "Mean Bottom Depth") +
  facet_wrap(~ model, nrow = 1) +
  theme(legend.position = "top")
dev.off()

png(here::here("figs", "model_comp", "day_cf_pred.png"), units = "in",
    height = 2.5, width = 6, res = 250)
bind_foo("local_day") %>% 
  plot_foo(data = ., x = "local_day") +  
  labs(x = "Day of Year") +
  facet_wrap(~ model, nrow = 1) +
  theme(legend.position = "top")
dev.off()

png(here::here("figs", "model_comp", "utmx_cf_pred.png"), units = "in",
    height = 2.5, width = 6, res = 250)
bind_foo("utm_x") %>% 
  plot_foo(data = ., x = "utm_x") +  
  labs(x = "Easting") +
  facet_wrap(~ model, nrow = 1) +
  theme(legend.position = "top")
dev.off()

png(here::here("figs", "model_comp", "utmy_cf_pred.png"), units = "in",
    height = 2.5, width = 6, res = 250)
bind_foo("utm_y") %>% 
  plot_foo(data = ., x = "utm_y") +  
  labs(x = "Northing") +
  facet_wrap(~ model, nrow = 1) +
  theme(legend.position = "top")
dev.off()

png(here::here("figs", "model_comp", "stage_cf_pred.png"), units = "in",
    height = 2.5, width = 6, res = 250)
bind_foo("med_stage") %>% 
  plot_foo(data = ., x = "med_stage") +  
  labs(x = "Maturity Stage") +
  facet_wrap(~ model, nrow = 1) +
  theme(legend.position = "top")
dev.off()
