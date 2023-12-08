## Fit hierarchical GAM

# Fit approximately equivalent model to ML alternative based on feedback from
# reviewers
# Fit model, estimate conditional effects and evaluate hold out performance on
# equivalent data
# POTENTIAL ISSUES ENCOUNTERED EARLIER:
# 1) autocorrelated residuals
# 2) slow convergence 

library(tidyverse)
library(mgcv)

depth_dat_raw1 <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  # approximately 7k detections have no available ROMS data; exclude 
  filter(!is.na(roms_temp)) %>% 
  # define time in hours for binning (all roms variables constant)
  mutate(
    date_time_hour = format(date_time_utm, format = "%Y-%m-%d %H"),
    day_night_dummy = ifelse(day_night == "night", 1, 0),
    stage_dummy = ifelse(stage == "mature", 1, 0)
  ) 


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

# with AR
depth_dat_ar <- depth_dat_raw %>% 
  group_by(vemco_code) %>% 
  mutate(
    start_time = min(date_time_utm),
    timestamp = difftime(start_time, date_time_utm, units = "mins"),
    timestamp = -1 * round(as.numeric(timestamp)),
    timestamp_f = cut_width(timestamp, width = 60, boundary = -0.1)
  ) %>% 
  group_by(vemco_code, 
           timestamp_f,
           fl, lipid, stage_dummy, utm_x, utm_y,
           day_night_dummy, local_day, mean_bathy, mean_slope, shore_dist,
           u, v, w, roms_temp, zoo, oxygen, thermo_depth) %>% 
  dplyr::mutate(
    moon_illuminated = mean(moon_illuminated),
    timestamp_n = mean(timestamp) + rnorm(1, 0, 0.01),
    date_time_utm = mean(date_time_utm),
    mean_rel_depth = mean(rel_depth),
    bin_id = row_number()#,
    # .groups = "drop"
  ) %>%
  ungroup() %>%
  left_join(., ind_folds, by = "vemco_code") %>% 
  # remove redundant observations
  filter(bin_id == "1")


# split by individual blocking (identical to random forest)
train_depth <- depth_dat_ar %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat_ar %>% filter(ind_block == "5") %>% droplevels()
test_depth_22 <- depth_dat_ar %>% 
  group_by(vemco_code, timestamp_n, fl, lipid, stage_dummy, utm_x, utm_y, 
           day_night_dummy, local_day, mean_bathy, mean_slope, shore_dist,
           u, v, w, roms_temp, zoo, oxygen, thermo_depth) %>% 
  summarize(
    moon_illuminated = mean(moon_illuminated),
    mean_rel_depth = mean(rel_depth),
    sd_rel_depth = sd(rel_depth),
    .groups = "drop"
  ) %>% 
  ungroup()


## SUBHOURLY VARIATION ---------------------------------------------------------

hist(depth_dat$sd_rel_depth)

# how balanced?
det_rates_binned <- depth_dat_ar %>% 
  group_by(vemco_code) %>% 
  summarize(
    n_dets = n(),
    n_days = length(unique(local_day))
  ) %>% 
  mutate(
    n_det_bin = cut(n_dets, breaks = c(0, 10, 100, 1000, 10000)),
    n_day_bin = cut(n_days, breaks = c(0, 10, 30, 100, 500))
  )

det_rates %>% 
  group_by(n_det_bin) %>% 
  tally()


## FIT -------------------------------------------------------------------------


fit1ar <- gamm(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k=c(10, 10)) +
    s(mean_bathy, k = 3) + 
    s(local_day, bs = "cc", k = 5) + 
    day_night_dummy + stage_dummy +# fl + lipid +
    s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit"),
  method = "REML"
)

gamm(
  pos_depth ~ region_f + #s(hour_c, bs = "cc", m = 2) +
    s(hour_c, by = region_f, bs = "cc") +
    s(max_bathy_c) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)

fit2 <- gam(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k = c(10, 10)) +
    s(mean_bathy, k = 3) + s(mean_slope, k = 3) + 
    s(local_day, bs = "cc", k = 5) + 
    day_night_dummy + stage_dummy +
    fl + lipid + zoo + shore_dist + moon_illuminated + 
    s(vemco_code, bs = "re"),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

fit3 <- gam(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k = c(10, 10)) +
    s(mean_bathy, k = 5) + s(mean_slope, k = 5) + 
    s(local_day, bs = "cc", k = 5) + 
    day_night_dummy + stage_dummy + thermo_depth +
    fl + lipid + s(zoo, k = 5) + shore_dist + moon_illuminated + u + v + w + 
    roms_temp +
    s(vemco_code, bs = "re"),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)


concurvity(fit1)
concurvity(fit2)

fit_list <- list(fit1, fit2, fit3)
saveRDS(fit_list, here::here("data", "model_fits", "gam_fits.rds"))


## CHECK RESIDUALS -------------------------------------------------------------

resids <- resid(fit2)

# normal distributed
hist(resids)

# qqplots
gam.check(fit2)

# temporal autocorrelation
acf(resids)

train_depth$resid <- resids

# look at 9 tags w/ most detections and plot temporal patterns in residuals
tt <- train_depth %>% 
  group_by(vemco_code) %>% 
  mutate(nn = n()) %>% 
  arrange(-nn) %>% 
  pull(vemco_code) %>% 
  unique()

ggplot(train_depth %>% 
         filter(vemco_code %in% tt[1:9])) + 
  geom_point(aes(x = local_day, y = resid)) +
  facet_wrap(~vemco_code, scales = "free") +
  ggsidekick::theme_sleek()


#
sims_list <- purrr::map(
  fit_list, 
  ~ simulate(.x, newdata = train_depth, nsim = 50)
)
fixed_pred_list <- purrr::map(
  fit_list, 
  ~ predict(.x, type = "response")
)
qq_list <- purrr::pmap(
  list(sims_list, fixed_pred_list), 
  function(y, z) {
    dharma_res <- DHARMa::createDHARMa(
      simulatedResponse = y,
      observedResponse = train_depth$mean_rel_depth,
      fittedPredictedResponse = z
    )
    plot(dharma_res)
  }
)
# looks good


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
    day_night_dummy = 0.5,
    stage_dummy = 0.5,
    vemco_code = unique(depth_dat_raw$vemco_code)[1]
  ) 

pred_fit <- predict(fit3, type = "response", se.fit = TRUE, newdata = pred_dat1,
                    #estimate population level, excluding RIs for tag
                    exclude = "s(vemco_code)")

pred_dat2 <- pred_dat1 %>%
  mutate(
    rel_pred_med = pred_fit$fit,
    rel_pred_se = pred_fit$se.fit,
    rel_pred_lo = rel_pred_med + (qnorm(0.025) * pred_fit$se.fit),
    rel_pred_up = rel_pred_med + (qnorm(0.975) * pred_fit$se.fit),
    pred_med = rel_pred_med * max_bathy,
    pred_lo = rel_pred_lo * max_bathy,
    pred_up = rel_pred_up * max_bathy,
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
              aes(x = utm_x_m, y = utm_y_m, fill = rel_pred_se)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Std. Error of Predictions", 
                       option = "C",
                       direction = -1)  +
  theme(legend.position = "top",
        axis.text = element_blank())

avg_depth1 <- cowplot::plot_grid(plotlist = list(rel_depth, mean_depth),
                                 ncol = 2)


## LATENT SPATIAL PROCESS ------------------------------------------------------

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

pred_latent_gam <- predict(
  fit3, type = "response", se.fit = TRUE, newdata = pred_latent,
  #estimate population level, excluding RIs for tag
  exclude = "s(vemco_code)"
)

pred_latent2 <- pred_latent %>%
  mutate(
    rel_pred_med = pred_latent_gam$fit %>% as.numeric(),
    rel_pred_se = pred_latent_gam$se.fit %>% as.numeric(),
    rel_pred_lo = rel_pred_med + (qnorm(0.025) * rel_pred_se),
    rel_pred_up = rel_pred_med + (qnorm(0.975) * rel_pred_se),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
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


## CONDITIONAL PREDICTIONS -----------------------------------------------------

stage_dat <- depth_dat_raw %>% 
  select(vemco_code, fl, lipid, stage) %>% 
  distinct() %>% 
  group_by(stage) %>% 
  dplyr::summarize(
    fl = mean(fl),
    lipid = mean(lipid)
  ) %>% 
  mutate(stage_dummy = ifelse(stage == "mature", 1, 0)) %>% 
  rbind(., 
        data.frame(
          stage = "average",
          fl = mean(depth_dat_raw$fl),
          lipid = mean(depth_dat_raw$lipid),
          stage_dummy = 0.5
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
  stage_dummy = 0.5,
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
  group_vals <- depth_dat %>% 
    dplyr::summarize(min_v = min(!!varname),
                     max_v = max(!!varname))
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
    dplyr::rename(!!varname := dum) 
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

counterfac_tbl <- tibble(
  var_in = c("mean_bathy", "fl", "utm_x", "utm_y",
             "local_day", "lipid", "thermo_depth",
             "roms_temp", "mean_slope",
             "day_night_dummy", "stage_dummy"
  )) %>%
  mutate(
    pred_dat_in = purrr::map(var_in,
                             gen_pred_dat),
    preds_ci1 = purrr::map(pred_dat_in, ~ pred_foo(fit = fit1, preds_in = .x)),
    preds_ci2 = purrr::map(pred_dat_in, ~ pred_foo(fit = fit2, preds_in = .x))#,
    # preds_ci3 = purrr::map(pred_dat_in, ~ pred_foo(fit = fit3, preds_in = .x))
  )


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
    coord_cartesian(ylim = c(-0.65, -0.1))
}



bathy_cond <- counterfac_tbl %>% 
  filter(var_in == "mean_bathy") %>% 
  pull(preds_ci2) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "mean_bathy") +  
  labs(x = "Mean Bottom\nDepth") +
  geom_text(x = -Inf, y = Inf, label = "a)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)

yday_cond <- counterfac_tbl %>% 
  filter(var_in == "local_day") %>% 
  pull(preds_ci2) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "local_day") +  
  labs(x = "Year\nDay") +
  geom_text(x = -Inf, y = Inf, label = "c)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)


slope_cond <- counterfac_tbl %>% 
  filter(var_in == "mean_slope") %>% 
  pull(preds_ci2) %>% 
  as.data.frame() %>% 
  plot_foo(data = .,
           x = "mean_slope") +  
  labs(x = "Mean\nSlope") +
  geom_text(x = -Inf, y = Inf, label = "b)", hjust = -0.5, vjust = 1.5,
            check_overlap = TRUE)


# maturity predictions
mat_pred_in <- rbind(
  new_dat %>% 
    mutate(
      stage_dummy = 0
    ),
  new_dat %>% 
    mutate(
      stage_dummy = 1
    )
) %>%
  # replace size and lipid with stage-specific averages
  select(-fl, -lipid) %>% 
  left_join(., stage_dat, by = "stage_dummy") %>% 
  distinct()

mat_preds <- pred_foo(fit = fit2, preds_in = mat_pred_in)

mat_cond <- mat_preds %>% 
  mutate(
    stage_f = factor(stage_dummy, levels = c(0, 1), 
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
                     limits = c(-0.65, -0.15)
  )


# day-night predictions
dn_pred_in <- rbind(
  new_dat %>% 
    mutate(
      day_night_dummy = 0
    ),
  new_dat %>% 
    mutate(
      day_night_dummy = 1
    )
) %>% 
  distinct()
dn_preds <- pred_foo(fit = fit2, preds_in = dn_pred_in)

dn_cond <- dn_preds %>% 
  mutate(
    dn_f = factor(day_night_dummy, levels = c(0, 1), 
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
                     limits = c(-0.65, -0.15)
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