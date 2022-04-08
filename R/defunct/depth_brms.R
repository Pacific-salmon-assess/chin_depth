## Based on model_depth.Rmd fit variety of brms models to identify drivers of 
# depth variability
# Ultimately include primary covariates, location scale structure and 
# autocorrelation parameters
# Dec. 20, 2021

library(tidyverse)
library(brms)
library(bayesplot)


# use only 2019/20 receivers for now
rec <- readRDS(here::here("data", "tagging_data", "for_glatos",
                          "receivers_2019_2020.RDS"))$rec_all 


# tibble made in genGlatosDetections.R (used for stage and GSI data)
chin_dat <- readRDS(here::here("data", "generated_data", 
                               "movement_tbl.RDS")) %>% 
  dplyr::select(-(dets:move_dat)) %>% 
  unnest(cols = c(chin_dat)) %>% 
  mutate(stage2 = case_when(
    stage == "unknown" ~ "mature",
    TRUE ~ stage
  )) 


# Moderately cleaned detections data (includes depth/temperature sensors)
depth_raw <- readRDS(here::here("data", "tagging_data", "for_glatos",
                                "detections_all.RDS")) %>% 
  filter(!grepl("2021", vemco_code))
  


# function to clean and bin depth data
depth1 <- depth_raw %>%
  filter(!is.na(depth))  %>% 
  #calculate timestep (width below in minutes) relative to first detection 
  group_by(vemco_code) %>% 
  mutate(start_time = min(date_time),
         timestamp = difftime(start_time, date_time, units = "mins"),
         timestamp = -1 * round(as.numeric(timestamp)),
         timestamp_f = cut_width(timestamp, width = 30, boundary = -0.1)) %>%
  # bin depth data by tag, receiver, and timestamp within a day
  group_by(vemco_code, timestamp_f, receiver, latitude, longitude,
           station_name) %>%
  summarize(
    #add small error to remove duplicates
    timestamp_n = mean(timestamp) + rnorm(1, 0, 0.01), 
    date_time = mean(date_time),
    depth = mean(depth),
    .groups = "drop") %>%
  ungroup() %>%
  left_join(., 
            rec %>% 
              dplyr::select(receiver = receiver_name, mean_bathy, max_bathy, 
                            sd_bathy, region),
            by = "receiver") %>%
  mutate(
    region_f = as.factor(region),
    region_f = fct_relevel(region_f, "swvi",  "nwwa", "jdf", 
                           "swwa", "sog", "puget", "columbia", "fraser"),
    depth = -1 * depth,
    # corrections for when depth is above the surface or deeper than max depth
    # in detection radius
    depth = ifelse(depth >= -0.1, -0.1, depth),
    depth_diff = depth - max_bathy,
    # adjust modest errors in depth relative to bathy
    depth = ifelse(depth_diff > -10 & depth_diff < 0, max_bathy - 0.5, depth),
    rel_depth = depth / max_bathy,
    # misc timestamp cleaning
    date_time_local = lubridate::with_tz(date_time, 
                                         tzone = "America/Los_Angeles"),
    hour = lubridate::hour(date_time_local),
    det_day = lubridate::yday(date_time_local),
    pos_mean_bathy = -1 * mean_bathy,
    pos_max_bathy = -1 * max_bathy) %>%
  left_join(., chin_dat, by = "vemco_code") %>%  
  dplyr::select(-c(depth_diff, region, date_time, injury:fin_dam, stage)) %>% 
  # remove large errors in depth relative to bottom bathymetry
  filter(!depth < max_bathy,
         !is.na(stage2)) %>% 
  mutate(vemco_code = as.factor(vemco_code),
         pos_depth = -1 * depth,
         # calculate expected fw entry date (Oct 18, i.e. max obs day for mat 
         # fish + 1 day)
         exit_day = ifelse(stage2 == "mature", 292, 292 + 365),
         days_to_exit = exit_day - det_day
  ) 

# subset for testing
set.seed(123)
unique_codes <- depth1 %>% 
  filter(stage2 == "mature") %>% 
  pull(vemco_code) %>% 
  sample(25)

# add sunrise/sunset data
sun_data <- data.frame(date = as.Date(depth1$date_time_local),
                       lat = depth1$latitude, 
                       lon = depth1$longitude)
temp <- suncalc::getSunlightTimes(data = sun_data,
                                  keep = c("sunrise", "sunset"),
                                  tz = "America/Los_Angeles")
depth30 <- cbind(depth1, temp %>% select(sunrise, sunset)) %>% 
  mutate(
    day_night = ifelse(date_time_local > sunrise & date_time_local < sunset, 
                       "day", "night")
  ) %>% 
  filter(vemco_code %in% unique_codes)


# exploratory plots
ggplot(depth30) +
  geom_point(aes(x = pos_max_bathy, y = pos_depth))
ggplot(depth30) +
  geom_point(aes(x = sqrt(pos_max_bathy), y = pos_depth))
ggplot(depth30) +
  geom_point(aes(x = det_day, y = depth)) +
  facet_grid(stage2~region_f)
ggplot(depth30) +
  geom_point(aes(x = hour, y = depth, fill = day_night), shape = 21) +
  facet_grid(stage2~region_f)




## FIT MODEL -------------------------------------------------------------------

# center variables for model fitting
depth30 <- depth30 %>% 
  mutate(
    sqrt_bathy_c = sqrt(pos_mean_bathy) - mean(sqrt(pos_mean_bathy)),
    hour_c = hour - mean(hour),
    yday_c = year_day - mean(year_day)
  )


# skeleton of model structure, still requires explanatory variables to be 
# finalized
brm_priors1 <- brm(
  pos_depth ~ s(sqrt_bathy_c, k = 4),
  data = depth30, family = Gamma(link = "log"), cores = 4, seed = 17,
  iter = 3500, warmup = 1000, thin = 10, refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
brm_priors2 <- brm(
  pos_depth ~ s(sqrt_bathy_c, k = 4),
  data = depth30, family = Gamma(link = "log"), cores = 4, seed = 17,
  iter = 3500, warmup = 1000, thin = 10, refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  prior=c(prior(normal(log(10), 2), class="Intercept"),
          prior(normal(0, 2), class="b"),
          prior(gamma(0.01,0.01),class="shape"))
)

saveRDS(brm1, 
        here::here("data", "generated_data", "depth_fits", "trim_brm_30.RDS"))

# rhat and neff look good
summary(brm1)

# check autocorrelation - seems reasonable
brm1_resid <- resid(brm1)[, "Estimate"]
plot(acf(brm1_resid, lag = 40))

# marginal smooths
msms <- conditional_smooths(brm1)
plot(msms)

# posterior predictive checks
pp_check(brm1)
pp_check(brm1, type='stat', stat='mean')
pp_check(brm1, type='error_scatter_avg')
pp_check(brm1, type='intervals')
pp_check(brm1, x = 'pos_max_bathy', type='error_scatter_avg_vs_x')


depth30 %>%
  droplevels() %>% 
  tidybayes::add_predicted_draws(brm1) %>%
  summarise(
    p_residual = mean(.prediction < y_star),
    z_residual = qnorm(p_residual),
    .groups = "drop_last"
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()


# shinystan for mixing/residuals 
launch_shinystan(brm1)






# more complex model may not be necessary with gamma distribution
brm1 <- brm(
  bf(pos_depth ~ s(sqrt(pos_max_bathy), k = 4) + 
       s(year_day, k = 4) + 
       s(hour, bs = "cc", m = 2) + s(hour, by = region_f, bs = "cc", m = 1) +
       (1 | vemco_code) + region_f + 
       ar(time = timestamp_n, gr = vemco_code, p = 1, cov = TRUE),
     shape ~ s(sqrt(pos_max_bathy), k = 4)),
  data = depth30, family = Gamma(link = "log"), cores = 4, seed = 17,
  iter = 3500, warmup = 1000, thin = 10, refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)