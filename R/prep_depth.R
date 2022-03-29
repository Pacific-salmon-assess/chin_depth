### Prep detections data for depth-specific analyses
## Jan. 13, 2022

library(tidyverse)

# receiver data (includes bathymetric data generated in chin_tagging repo)
rec <- readRDS(here::here("data", "receivers_all.RDS"))$rec_all 

# life stage estimates
stage_dat <- readRDS(here::here("data", "lifestage_df.RDS")) %>% 
   select(vemco_code, stage_predicted)

# hourly ROMS outputs matched to receiver stations (marine only and some 
# missing), restricted to 25m strata
roms_dat <- readRDS(here::here("data", "roms_25m_depth.RDS")) %>% 
  # exclude rho (highly correlated with temperature)
  select(-rho)


# moderately cleaned detections data (includes depth/temperature sensors)
depth_raw <- readRDS(here::here("data", "detections_all.RDS")) %>%
  filter(!flag == "yes",
         !is.na(depth),
         depth > 0) %>%
  left_join(., 
            rec %>% 
              dplyr::select(receiver = receiver_name, mean_bathy = mean_depth,
                            max_bathy = max_depth, mean_slope:shore_dist,
                            marine) %>% 
              distinct(),
            by = "receiver") %>% 
  filter(marine == "yes") %>% 
  left_join(., stage_dat %>% dplyr::rename(stage = stage_predicted), 
            by = "vemco_code") %>% 
  select(-temp, -flag, -station_name, -marine)


# B. Hendricks data
# NOTE: differences in receiver formatting may lead to issues...
depth_h <- readRDS(here::here("data", "hendricks_depth_dets.RDS")) %>% 
  mutate(receiver_sn = as.character(receiver_sn),
         stage = "mature") %>% 
  dplyr::rename(receiver = receiver_name)



# function to make hours continuous
time_foo <- function(x) {
  lubridate::hour(x) + (lubridate::minute(x) / 60) + 
    (lubridate::second(x) / 3600) 
}


## drop binning for now (see below)
depth_dat <- rbind(depth_raw, depth_h) %>%
  mutate(
    hour_int = lubridate::hour(date_time) + 1,
    day = lubridate::day(date_time),
    month = lubridate::month(date_time),
    year = lubridate::year(date_time),
    region_f = fct_relevel(as.factor(region), "swvi",  "nwwa", "jdf", 
                           "swwa", "sog", "puget", "columbia", "fraser"),
    # corrections for when depth is above the surface or deeper than max depth
    # in detection radius
    depth = ifelse(depth <= 0.05, 0.05, depth),
    depth_diff = max_bathy - depth,
    # adjust modest errors in depth relative to bathy
    depth = ifelse(depth_diff > -10 & depth_diff < 0, max_bathy - 0.5, depth),
    rel_depth = depth / max_bathy,
    # misc timestamp cleaning
    date_time_local = lubridate::with_tz(date_time, 
                                         tzone = "America/Los_Angeles"),
    hour = time_foo(date_time_local),
    det_day = lubridate::yday(date_time_local),
    vemco_code = as.factor(vemco_code)
  ) %>%
  left_join(
    ., roms_dat, 
    by = c("latitude", "longitude", "day", "month", "year", "hour_int")
  ) %>% 
  filter(
    # remove large errors in depth relative to bottom bathymetry
    depth < max_bathy
  ) %>% 
  droplevels()

# add sunrise/sunset data
sun_data <- data.frame(date = as.Date(depth_dat$date_time_local),
                       lat = depth_dat$latitude, 
                       lon = depth_dat$longitude)
temp <- suncalc::getSunlightTimes(data = sun_data,
                                  keep = c("sunrise", "sunset"),
                                  tz = "America/Los_Angeles")

# combine and subset
dat_out <- cbind(depth_dat, temp %>% dplyr::select(sunrise, sunset)) %>% 
  mutate(
    day_night = ifelse(date_time_local > sunrise & date_time_local < sunset, 
                       "day", "night")
  ) %>% 
  dplyr::select(
    vemco_code, stage, receiver:longitude, mean_bathy:shore_dist, 
    u, v, w, roms_temp, zoo,
    region_f, date_time_local, hour, day_night,
    det_day, year, pos_depth = depth, rel_depth
  ) 


saveRDS(dat_out,
        here::here("data", "depth_dat_nobin.RDS"))
