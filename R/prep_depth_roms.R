### ROMS Data Clean
## Export stations for query, then import and evaluate 
## Query conducted by Doug Jackson (doug@qedaconsulting.com)
## Feb 14, 2022
## Updated May 23, 2023

library(tidyverse)


## FUNCTIONS -------------------------------------------------------------------

# add coordinates in UTM space
lonlat_to_utm <- function(x, y, zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  sp::coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- sp::spTransform(
    xy, 
    sp::CRS(
      paste("+proj=utm +zone=",zone," +units=m +datum=WGS84 ellps=WGS84",
            sep = '')
      
    ))
  return(as.data.frame(res))
}


# function to make hours continuous
time_foo <- function(x) {
  lubridate::hour(x) + (lubridate::minute(x) / 60) + 
    (lubridate::second(x) / 3600) 
}


## PREP DEPTH DATA -------------------------------------------------------------

# receiver data (includes bathymetric data generated in chin_tagging repo)
rec_full <- readRDS(here::here("data", "receivers_all.RDS"))
rec <- rec_full$rec_all 


# initial life stage estimates (includes unknowns and fish IDd via PIT so 
# n_stagedat > n_chinin)
stage_dat <- readRDS(here::here("data", "agg_lifestage_df.RDS")) %>% 
  dplyr::select(fish, stage, agg) 
chin_in <- readRDS(here::here("data", "acousticOnly_GSI.RDS")) %>% 
  left_join(., stage_dat, by = "fish") %>% 
  filter(!year == "2023")


## predict life stage based on fitted model 
stage_mod <- readRDS(here::here("data", "stage_fl_hierA.RDS"))
# include RIs when stock ID is known
pred_dat <- tidybayes::add_epred_draws(
  object = stage_mod,
  newdata = chin_in %>% 
    filter(is.na(stage) | stage == "unknown",
           !is.na(agg)),
  re_formula = NULL, #scale = "response", 
  ndraws = 100
)

# otherwise marginalize RIs
pred_dat_na <- tidybayes::add_epred_draws(
  object = stage_mod,
  newdata = chin_in %>% 
    filter(is.na(stage) | stage == "unknown",
           is.na(agg)),
  re_formula = NA, #scale = "response", 
  ndraws = 100
)

pred_stage <- rbind(pred_dat, pred_dat_na) %>% 
  group_by(fish) %>% 
  summarize(med = median(.epred), 
            lo = quantile(.epred, prob = 0.05),
            up = quantile(.epred, prob = 0.95),
            .groups = "drop") 

chin_stage <- left_join(chin_in, pred_stage, by = "fish") %>% 
  mutate(
    med_stage = case_when(
      stage == "mature" ~ 1,
      stage == "immature" ~ 0,
      #use median predicted posterior value when stage unknown
      TRUE ~ med
    ),
    #correct missing agg ID
    agg = ifelse(cu_name == "South Puget Sound", "Puget Sound", agg),
    agg = agg_name
  ) %>% 
  filter(!is.na(acoustic_year)) %>%
  select(vemco_code = acoustic_year, tagging_year = year, acoustic_type, fl, 
         med_stage, lipid, cu_name, agg) %>% 
  filter(grepl("V13P", acoustic_type))


# moderately cleaned detections data (includes depth/temperature sensors)
depth_raw <- readRDS(here::here("data", "detections_all.RDS")) %>%
  left_join(., chin_stage, by = "vemco_code") %>% 
  left_join(., 
            rec %>% 
              select(receiver = receiver_name, trim_sn, mean_depth:shore_dist),
            by = "receiver") %>%
  filter(
    # remove stations that are freshwater or terminal and non-depth tags
    !station_name == "LakeWashington",
    grepl("V13P", acoustic_type)
  ) %>%
  select(vemco_code, tagging_year, med_stage, fl, lipid, cu_name, agg, 
         date_time, latitude, 
         longitude, receiver_name = receiver, trim_sn, region, 
         mean_bathy = mean_depth, max_bathy = max_depth, mean_slope, mean_aspect,
         shore_dist, depth)


# B. Hendricks data with energy meter rescaled to lipid
calc_lipid <- function(mean_energy) {
  ifelse(
    mean_energy < 1.059, 
    -1.973 + (6.758 * mean_energy),
    -1.973 + 6.758 * 1.059 + ((6.758 - 6.202) * (mean_energy - 1.059))
  )
}

# Import UBC data
hendricks_bio <- readRDS(here::here("data", "hendricks_chin_dat.RDS")) %>% 
  mutate(
    lipid = calc_lipid(exp(mean_log_e)),
    tagging_year = lubridate::year(date)
  )

depth_h <- readRDS(here::here("data", "hendricks_depth_dets.RDS")) %>% 
  mutate(
    trim_sn = as.character(receiver_sn),
    week = lubridate::week(date_time)
  ) %>% 
  left_join(., 
            hendricks_bio %>% 
              select(vemco_code, tagging_year, fl, lipid, cu_name, agg), 
            by = "vemco_code") %>%
  select(colnames(depth_raw)) 


# combine
depth_dets1 <- rbind(depth_raw, depth_h) %>% 
  mutate(
    hour = lubridate::hour(date_time) + 1,
    day = lubridate::day(date_time),
    month = lubridate::month(date_time),
    year = lubridate::year(date_time),
    cu_name = ifelse(is.na(cu_name), "unknown", cu_name),
    agg = ifelse(is.na(agg), "unknown", agg)
  ) %>% 
  filter(!region == "fraser",
         !is.na(depth)) 


# ~2% of tags lack energy density estimates, impute based on size, stock ID and 
# year
full_bio_dat <- depth_dets1 %>% 
  select(vemco_code, med_stage, fl, lipid, cu_name, agg, tagging_year) %>%
  distinct() 
interp_lipid <- full_bio_dat %>%
  select(-vemco_code) %>% 
  VIM::kNN(., k = 5) %>% 
  select(-ends_with("imp")) 
full_bio_dat$lipid <- interp_lipid$lipid

full_locs_dat <- depth_dets1 %>% 
  select(latitude, longitude, mean_bathy, max_bathy, mean_slope, mean_aspect, 
         shore_dist, receiver_name) %>% 
  distinct()
interp_bathy <- full_locs_dat %>% 
  select(-receiver_name) %>% 
  VIM::kNN(., k = 5) %>% 
  select(-ends_with("imp")) 
interp_bathy$receiver_name <- full_locs_dat$receiver_name

depth_dets2 <- depth_dets1 %>% 
  select(-c(lipid, latitude, longitude, mean_bathy, max_bathy, mean_slope, 
            mean_aspect, shore_dist)) %>% 
  left_join(., full_bio_dat %>% select(vemco_code, lipid), by = "vemco_code") %>% 
  left_join(., 
            interp_bathy, 
            by = "receiver_name")


## CHECK FOR FAULTY TAGS -------------------------------------------------------

deaths <- depth_dets2 %>%
  group_by(vemco_code, receiver_name) %>%
  mutate(nn = n()) %>%
  filter(
    nn > 300
  )

ggplot() +
  geom_point(data = depth_dets2 %>% filter(vemco_code %in% deaths$vemco_code),
             aes(x = date_time, y = depth, fill = region),
             shape = 21) +
  facet_wrap(~vemco_code, scales = "free_x") +
  ggsidekick::theme_sleek()
# no evidence that fish died within detection distance of receiver

# faulty tags
sensors_fail <- depth_dets2 %>% 
  mutate(time_stamp = paste(hour, day, month, year, sep = "_")) %>% 
  group_by(time_stamp, vemco_code) %>% 
  mutate(
    sd_depth = sd(depth),
    n = length(unique(date_time))
  ) %>% 
  filter(sd_depth < 0.1 & n > 5) %>% 
  droplevels()

ggplot(depth_dets2 %>% 
         filter(vemco_code %in% unique(sensors_fail$vemco_code)[1:12])) +
  geom_point(aes(x = date_time, y = depth, fill = region),
             shape = 21) +
  facet_wrap(~vemco_code, scales = "free_x") +
  ggsidekick::theme_sleek()
# no evidence that tags pressure sensors failed permanently


# check above surface tags
surface_tags <- depth_dets2 %>% 
  group_by(vemco_code) %>% 
  mutate(max_depth = min(depth)) %>% 
  filter(max_depth < 0) %>% 
  pull(vemco_code) %>% 
  unique()

ggplot(depth_dets2 %>% 
         filter(vemco_code %in% surface_tags)) +
  geom_point(aes(x = date_time, y = depth, fill = region),
             shape = 21) +
  facet_wrap(~vemco_code, scales = "free_x") +
  ggsidekick::theme_sleek()
# no evidence that above surface detections were consistent



## EXPORT STATIONS -------------------------------------------------------------

## Three different ROMS datasets
# 1) surface conditions for each variable-detection
# 2) temperature and oxygen profiles for each detection 
# 3) predictive grids for two time periods (July 1 and Feb 1)

# expand depth dets by different covariates and depths
trim_dets <- depth_dets2 %>% 
  select(year, month, day, hour, lat = latitude, lon = longitude) %>% 
  distinct()
var_list <- c("u", "v", "w", "temp", "zooplankton", "oxygen")
# focus on surface depth given detections in shallows will be excluded otherwise
# depth_list <- c(#5, 
#   , 25, 50
#   )

# surface conditions
surface_out <- expand.grid(
  variable = var_list,
  depth = 1#depth_list
) %>% 
  mutate(
    fac = paste(variable, depth, sep = "_")
  ) %>% 
  distinct() %>% 
  split(., .$fac) %>% 
  purrr::map(., function (x) {
    trim_dets %>% 
      mutate(variable = x$variable,
             depth = x$depth)
  }) %>% 
  bind_rows %>% 
  mutate(depthFrac = -999) %>% 
  dplyr::select(year, month, day, hour, variable, lat,
                lon, depth, depthFrac) 


# profiles
profiles_out <- expand.grid(
  variable = c("oxygen", "temp"),
  depth = seq(0, 350, by = 5)
) %>% 
  mutate(
    fac = paste(variable, depth, sep = "_")
  ) %>% 
  distinct() %>% 
  split(., .$fac) %>% 
  purrr::map(., function (x) {
    trim_dets %>% 
      mutate(variable = x$variable,
             depth = x$depth)
  }) %>% 
  bind_rows %>% 
  mutate(depthFrac = -999) %>% 
  dplyr::select(year, month, day, hour, variable, lat,
                lon, depth, depthFrac)


# prediction grid
pred_grid <- readRDS(here::here("data", "pred_bathy_grid_utm_no_bark.RDS"))

pred_grid_ll <- sf::st_transform(
  pred_grid, 
  crs = sp::CRS("+proj=longlat +datum=WGS84")
) 
pred_out <- expand.grid(
  variable = var_list,
  time_stamp = c("2020-07-15", "2021-02-15")
) %>% 
  mutate(
    time_stamp = as.POSIXct(
      time_stamp,
      format = "%Y-%m-%d"
    ),
    fac = paste(variable, as.factor(time_stamp), sep = "_")
    ) %>% 
  distinct() %>% 
  split(., .$fac) %>% 
  purrr::map(., function (x) {
    data.frame(
      sf::st_coordinates(pred_grid_ll[ , 1]),
      variable = x$variable,
      time_stamp = x$time_stamp,
      depth = 1,
      mean_bathy_depth = pred_grid_ll$depth,
      max_bathy_depth = pred_grid_ll$max_depth,
      slope = pred_grid_ll$slope,
      shore_dist = pred_grid_ll$shore_dist
      ) %>% 
      mutate(
        year = lubridate::year(time_stamp), 
        month = lubridate::month(time_stamp), 
        day = lubridate::day(time_stamp), 
        hour = 12
      ) 
  }) %>% 
  bind_rows %>% 
  mutate(depthFrac = -999) %>%
  rename(lat = Y, lon = X)

# export 
write.csv(surface_out, here::here("data", "stations_roms_no_infill.csv"),
          row.names = FALSE)
# write.csv(profiles_out, 
#           here::here("data", "vert_profiles_stations_roms_no_infill.csv"),
#           row.names = FALSE)
write.csv(pred_out, 
          here::here("data", "pred_stations_roms_no_infill.csv"),
          row.names = FALSE)


## IMPORT PROFILE DATA ---------------------------------------------------------

profile_dat <- read.csv(
  here::here("data", "raw_roms", 
             "stations_roms_no_infill_profiles_10may23_all.csv")) %>%
  left_join(., 
            rec %>% 
              select(lat = station_latitude, lon = station_longitude,
                     region) %>% 
              distinct(),
            by = c("lat", "lon"))  %>% 
  filter(!region == "fraser") %>%
  mutate(
    unique_id = paste(year, month, day, hour, lat, lon) %>%
      as.factor() %>%
      as.numeric(),
    date = paste(day, month, year, sep = "-"),
    time = paste(hour, "00", sep = ":"),
    time1 = as.POSIXct("01-07-2019 12:00", form = "%d-%m-%Y %H:%M"),
    timestamp = paste(date, time, sep = " ") %>% 
      as.POSIXct(., form = "%d-%m-%Y %H:%M") %>%
      difftime(time1, timestamp) %>% 
      as.numeric()
  ) %>% 
  filter(!depth < 0) 

# calculate thermocline depth at each station
thermo_dat<- profile_dat %>% 
  filter(variable == "temp") %>% 
  group_by(unique_id, lat, lon, timestamp, region) %>% 
  summarize(
    thermo_depth = rLakeAnalyzer::thermo.depth(wtr = value, depths = depth),
    .groups = "drop"
  ) %>% 
  group_by(region) %>% 
  mutate(median_thermo = median(thermo_depth, na.rm = TRUE)) %>% 
  ungroup()

pdf(here::here("figs", "thermocline_depth.pdf"))
ggplot(thermo_dat %>% filter(!is.na(thermo_depth))) +
  geom_histogram(aes(x = thermo_depth)) +
  geom_vline(aes(xintercept = median_thermo), colour = "red") +
  facet_wrap(~region) +
  ggsidekick::theme_sleek()
dev.off()


## IMPORT ROMS PULL ------------------------------------------------------------

# query generated by Doug; note some stations have missing values due to outside
# survey domain
roms_in <- read.csv(
  here::here("data", "raw_roms", "stations_roms_no_infill_10may23_all.csv"))

roms_dat <- roms_in %>%
  left_join(.,
            rec %>%
              select(lat = station_latitude, lon = station_longitude,
                     region, mean_depth, mean_slope,
                     shore_dist ) %>%
              distinct(),
            by = c("lat", "lon")
            ) %>%
  filter(!region == "fraser") %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  # add time step variable used during imputation and for ID 
  mutate(
    date = paste(day, month, year, sep = "-"),
    time = paste(hour, "00", sep = ":"),
    time1 = as.POSIXct("01-07-2019 12:00", form = "%d-%m-%Y %H:%M"),
    timestamp = paste(date, time, sep = " ") %>% 
      as.POSIXct(., form = "%d-%m-%Y %H:%M") %>%
      difftime(time1, timestamp) %>% 
      as.numeric()
  ) %>% 
  left_join(., 
            thermo_dat %>% select(-c(region, median_thermo)),
            by = c("lat", "lon", "timestamp")) %>% 
  select(-depth, -depthFrac, -c(date:time1)) %>% 
  # remove duplicates introduced by receivers being deployed at identical 
  # locations
  distinct() %>% 
  rename(zoo = zooplankton, latitude = lat, longitude = lon,
         roms_temp = temp) 

depth_utm <- lonlat_to_utm(roms_dat$longitude, roms_dat$latitude, 
                           zone = 10) 
roms_dat$utm_x <- depth_utm$X / 1000 
roms_dat$utm_y <- depth_utm$Y / 1000

# replace unknown values with NA
roms_dat[roms_dat == "-999"] <- NA



# join with receiver data (includes bathymetric variables) and impute using
# k nearest neighbors
# NOTE restricted to stations where at least some roms data were produced
roms_trim <- roms_dat %>% 
  filter(!if_all(c(oxygen:roms_temp, thermo_depth), ~ is.na(.))) 

# identify proportion of NA by each ROMS variable
roms_trim %>% 
  pivot_longer(., 
               cols = c(roms_temp, w, v, oxygen, u, zoo, thermo_depth), 
               names_to = "roms_var") %>% 
  filter(is.na(value)) %>% 
  group_by(roms_var) %>% 
  tally() %>% 
  mutate(n / nrow(roms_trim)) 


roms_interp <- VIM::kNN(roms_trim,
                        variable = c("roms_temp", "w", "v", "oxygen", "u", "zoo", 
                                     "thermo_depth"),
                        dist_var = c("timestamp", "utm_x", "utm_y",
                                     "mean_depth", "mean_slope", "shore_dist"),
                        k = 5)

# remove flags and non-ROMS data
roms_interp_trim <- roms_interp %>% 
  select(year:longitude, roms_temp, w, v, oxygen, u, zoo, thermo_depth) 


# as above but imputing ALL variables
roms_interp2 <- VIM::kNN(roms_dat,
                        variable = c("roms_temp", "w", "v", "oxygen", "u", "zoo", 
                                     "thermo_depth"),
                        dist_var = c("timestamp", "utm_x", "utm_y",
                                     "mean_depth", "mean_slope", "shore_dist"),
                        k = 5)

# remove flags and non-ROMS data
roms_interp2_trim <- roms_interp2 %>% 
  select(year:longitude, roms_temp, w, v, oxygen, u, zoo, thermo_depth) 



saveRDS(roms_interp_trim, 
        here::here("data", "interp_roms_surf_depth.RDS"))
saveRDS(roms_interp2_trim, 
        here::here("data", "interp_roms_surf_depth_FULL.RDS"))


## CLEAN DEPTH DATA ------------------------------------------------------------

roms_interp_trim <- readRDS(here::here("data", "interp_roms_surf_depth.RDS"))

depth_utm <- lonlat_to_utm(depth_dets2$longitude, depth_dets2$latitude, 
                           zone = 10) 
depth_dets2$utm_x <- depth_utm$X / 1000 
depth_dets2$utm_y <- depth_utm$Y / 1000


# add roms data (approximately 6.5k detections have no ROMS data)
depth_dat <- depth_dets2 %>% 
  left_join(
    ., 
    roms_interp_trim, 
    by = c("latitude", "longitude", "day", "month", "year", "hour"),
    
  ) %>% 
  ## some tags detected simultaneously at multiple receivers at the same time; 
  # to resolve take average of spatial variables, including roms data  
  group_by(vemco_code, date_time) %>% 
  mutate(
    nn = n(),
    rep_number = row_number(),
    # depth = mean(depth),
    mean_bathy = ifelse(nn > 1, mean(mean_bathy), mean_bathy),
    max_bathy = ifelse(nn > 1, mean(max_bathy), max_bathy),
    mean_slope = ifelse(nn > 1, mean(mean_slope), mean_slope),
    mean_aspect = ifelse(nn > 1, mean(mean_aspect), mean_aspect),
    shore_dist = ifelse(nn > 1, mean(shore_dist), shore_dist),
    depth = ifelse(nn > 1, mean(depth), depth),
    utm_x = ifelse(nn > 1, mean(utm_x), utm_x),
    utm_y = ifelse(nn > 1, mean(utm_y), utm_y),
    roms_temp = ifelse(nn > 1, mean(roms_temp), roms_temp),
    w = ifelse(nn > 1, mean(w), w),
    v = ifelse(nn > 1, mean(v), v),
    oxygen = ifelse(nn > 1, mean(oxygen), oxygen),
    u = ifelse(nn > 1, mean(u), u),
    zoo = ifelse(nn > 1, mean(zoo), zoo),
    thermo_depth = ifelse(nn > 1, mean(thermo_depth), thermo_depth)
  ) %>% 
  filter(
    rep_number == "1"
  ) %>% 
  ungroup() %>% 
  mutate(
    start_time = min(date_time),
    timestamp = difftime(start_time, date_time, units = "mins"),
    timestamp_n = -1 * round(as.numeric(timestamp)),
    region_f = fct_relevel(as.factor(region), "swvi",  "nwwa", "jdf", 
                           "swwa", "sog", "puget", "columbia"),
    # corrections for when depth deeper than max bathy depth
    depth_diff = max_bathy - depth,
    depth = ifelse(depth_diff > -2.5 & depth_diff < 0, max_bathy - 0.5, depth),
    # corrections for when depth slightly above surface
    depth = ifelse(depth < 0 & depth > -1.5, 0.5, depth),
    # misc timestamp cleaning
    date_time_local = lubridate::with_tz(date_time, 
                                         tzone = "America/Los_Angeles"),
    local_hour = time_foo(date_time_local),
    local_day = lubridate::yday(date_time_local),
    vemco_code = as.factor(vemco_code),
    tag_year = str_split(vemco_code, "_") %>% 
      purrr::map(., function (x) x[2]) %>% 
      as.numeric(),
    # redefine stage so that immature fish become mature May 1 of following
    # year
    med_stage = ifelse(
      med_stage < 0.1 & year > tag_year & local_day > 120,
      1,
      med_stage
    )
  ) %>% 
  filter(
    # remove large errors in depth relative to bottom bathymetry or above surface
    depth < max_bathy,
    depth > 0
  ) %>% 
  droplevels()


# add moon data
moon_data <- oce::moonAngle(depth_dat$date_time, 
                            depth_dat$longitude, 
                            depth_dat$latitude)$illuminatedFraction


# add sunrise/sunset data
sun_data <- data.frame(date = as.Date(depth_dat$date_time_local),
                       lat = depth_dat$latitude, 
                       lon = depth_dat$longitude)
temp <- suncalc::getSunlightTimes(data = sun_data,
                                  keep = c("sunrise", "sunset"),
                                  tz = "America/Los_Angeles")

depth_dat2 <- cbind(depth_dat, temp %>% dplyr::select(sunrise, sunset)) %>%
  mutate(
    rel_depth = depth / max_bathy,
    logit_rel_depth = boot::logit(rel_depth),
    day_night = ifelse(date_time_local > sunrise & date_time_local < sunset,
                       "day", "night"),
    moon_illuminated = moon_data,
    #create cyclical time steps representing year day
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365)
  ) %>%
  dplyr::select(
    vemco_code, cu_name, agg, fl, lipid, med_stage, trim_sn, 
    receiver_name, latitude, longitude, utm_y, utm_x, max_bathy,
    mean_bathy, mean_slope, 
    shore_dist, u, v, w, roms_temp, zoo, oxygen, thermo_depth,
    region_f, date_time_utm = date_time, date_time_local, timestamp_n, 
    local_day, det_dayx, det_dayy, year, 
    local_hour, day_night, moon_illuminated, pos_depth = depth,
    rel_depth, logit_rel_depth) 

saveRDS(depth_dat2, here::here("data", "depth_dat_nobin.RDS"))


## CLEAN PREDICTIVE DATA -------------------------------------------------------

# dataframe passed to Doug Jackson lacked bathy data and due to conversion to 
# lat/lon is no longer an evenly spaced grid; to correct merge ROMS data
# w/ original prediction grid after converting both to UTM

# roms data in lat/lon
pred_roms_in <- read.csv(
  here::here(
    "data", "raw_roms", "pred_stations_roms_no_infill_10may23_all.csv"
  )) %>% 
  mutate(
    # add NAs
    value = na_if(value, -999)
  )

# number of NAs by variable
pred_roms_in %>% 
  filter(is.na(value)) %>% 
  group_by(variable, month) %>% 
  tally()

pred_roms_in %>% 
  filter(variable == "v" & is.na(value)) %>% 
  glimpse()
pred_roms_in %>% 
  filter(variable == "u" & is.na(value)) %>% 
  glimpse()

pred_utm <- lonlat_to_utm(pred_roms_in$lon, pred_roms_in$lat, zone = 10) 

# combine and export
# NOTE: excludes biological data and thermocline depth (add means in modeling
# script)
pred_roms <- pred_roms_in %>% 
  mutate(
    x = pred_utm$X,
    y = pred_utm$Y,
    local_day = ifelse(month == "7", 211, 46),
    day_night_night = 0.5,
    moon_illuminated = 0.5,
    det_dayx = sin(2 * pi * local_day / 365),
    det_dayy = cos(2 * pi * local_day / 365),
    season = fct_recode(as.factor(month), "winter" = "2", "summer" = "7")
  ) %>% 
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  select(
    lat, lon, x, y, local_day:oxygen, roms_temp = temp, u, v, w, 
    zoo = zooplankton
  ) 

pred_roms_sf <- sf::st_as_sf(pred_roms, coords = c("x", "y"),
                                 crs = sp::CRS("+proj=utm +zone=10 +units=m"))

pred_grid_sf <- readRDS(here::here("data", "pred_bathy_grid_utm_no_bark.RDS"))
pred_grid_sf$utm_x <- sf::st_coordinates(pred_grid_sf[ , 1])[, "X"] / 1000
pred_grid_sf$utm_y <- sf::st_coordinates(pred_grid_sf[ , 1])[, "Y"] / 1000

dum <- sf::st_join(pred_roms_sf, pred_grid_sf, join = nngeo::st_nn, 
                   maxdist = 1000, k = 1, progress = FALSE)

pred_roms_out <- sf::st_drop_geometry(dum) %>% 
  rename(mean_bathy = depth, max_bathy = max_depth, mean_slope = slope)

ggplot(pred_roms_out %>% filter(season == "summer", !is.na(u))) +
  geom_raster(aes(x = utm_x, y = utm_y, fill = roms_temp))

# export removing cells that lack data for any ROMS variable
saveRDS(pred_roms_out %>% 
          filter(!is.na(u), !is.na(v)), 
        here::here("data", "pred_bathy_grid_roms.RDS"))


## BIOLOGICAL TRAITS FIGURE ----------------------------------------------------

bio_dat <- rbind(
  chin_stage %>% 
    select(vemco_code, lipid, fl, med_stage, cu_name, agg),
  hendricks_bio %>% 
    select(vemco_code, lipid, fl, med_stage, cu_name, agg)
) %>% 
  mutate(
    agg = case_when(
      grepl("East Vancouver Island", cu_name) ~ "ECVI",
      grepl("Willamette", cu_name) ~ "LowCol",
      grepl("Okanagan", cu_name) ~ "Up Col.",
      grepl("Upper Columbia", cu_name) ~ "Up Col.",
      grepl("_1.", cu_name) ~ "FraserYear",
      grepl("Vancouver Island", cu_name) ~ "WCVI",
      grepl("_0.", cu_name) | agg == "Fraser" ~ "FraserSub",
      cu_name == "Lower Columbia River" ~ "LowCol",
      TRUE ~ agg
    ),
    agg = factor(
      agg,
      levels = c("Cali", "WA_OR", "Up Col.", "LowCol", "Puget Sound", "FraserYear",
                 "FraserSub", "ECVI", "WCVI"),
      labels = c("Cali.", "WA/OR.", "Upriver\nCol.", "Lower\nCol.",
                 "Puget\nSound", "Fraser\nYear.", "Fraser\nSub.", "ECVI",
                 "WCVI")),
    stage = ifelse(med_stage > 0.5, "mature", "immature")
  )

n_tags_fl <- bio_dat %>% 
  filter(!is.na(fl)) %>% 
  group_by(stage, agg) %>% 
  tally()
n_tags_lipid <- bio_dat %>% 
  filter(!is.na(lipid)) %>% 
  group_by(stage, agg) %>% 
  tally()

# summary stats
bio_dat %>% 
  group_by(stage) %>% 
  summarize(
    mean_fl = mean(fl, na.rm = T),
    mean_lipid = mean(lipid, na.rm = T),
    sd_fl = sd(fl, na.rm = T),
    sd_lipid = sd(lipid, na.rm = T),
    n = n()
  ) 
  

lipid_box <- ggplot(
  bio_dat %>% 
    filter(!is.na(lipid)), 
  aes(x = agg, fill = stage)
) +
  geom_boxplot(aes(y = lipid)) +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  ggsidekick::theme_sleek() +
  geom_text(data = n_tags_lipid, aes(y = -Inf, label  = n),
            position = position_dodge(width = 0.75),
            hjust = 0.5, vjust = -0.5) +
  labs(y = "Lipid Content (% Wet Weight)") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

fl_box <- ggplot(
  bio_dat %>% 
    filter(!is.na(fl)), 
  aes(x = agg, fill = stage)
) +
  geom_boxplot(aes(y = fl)) +
  scale_fill_brewer(type = "seq", palette = "Purples") +
  ggsidekick::theme_sleek() +
  geom_text(data = n_tags_fl, aes(y = -Inf, label  = n),
            position = position_dodge(width = 0.75),
            hjust = 0.5, vjust = -0.5) +
  ylim(c(56, 102)) +
  labs(y = "Fork Length (cm)") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )


png(here::here("figs", "ms_figs_rel", "bio_boxplot.png"),
    height = 170, width = 170, units = "mm", res = 300)
cowplot::plot_grid(
  fl_box, lipid_box, ncol = 1
)
dev.off()


stock_tbl <- depth_dets1 %>% 
  select(cu_name, agg) %>% 
  distinct() %>% 
  arrange(agg)
write.csv(stock_tbl, here::here("data", "stock_definitions.csv"))


## DEPTH SCATTER FIGURE --------------------------------------------------------

depth_dat2 <- readRDS(here::here("data", "depth_dat_nobin.RDS"))

bottom_dot <- ggplot(depth_dat2 %>% 
                       filter(
                         !is.na(roms_temp)
                       )) +
  geom_point(aes(x = max_bathy, y = -1 * pos_depth, colour = med_stage), 
             alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = -1)) +
  labs(y = "Observed Depth (m)", 
       x = "Maximum Bottom Depth\nWithin Detection Radius (m)") +
  ggsidekick::theme_sleek() +
  scale_y_continuous(
    breaks = c(0, -100, -200, -300), 
    labels = seq(0, 300, by = 100)
  ) +
  scale_colour_gradient2(midpoint = 0.5, name = "Maturation Stage",
                         n.breaks = 3) +
  theme(legend.position = "top")

png(here::here("figs", "ms_figs_rel", "depth_vs_bathy.png"),
    height = 85, width = 85, units = "mm", res = 300)
bottom_dot
dev.off()


## MISCELLANEOUS DATA  ---------------------------------------------------------

# number of immature and mature fish tagged
bio_dat %>% 
  group_by(stage) %>%
  tally()

depth_dat_raw1 %>% 
  select(vemco_code, stage) %>% 
  distinct() %>% 
  group_by(stage) %>% 
  tally()

# number of detections per tag
n_det_dat <- depth_dat2 %>% 
  group_by(vemco_code) %>% 
  tally() 

# calculate timespan overwhich detections provided (merge with depth_raw to 
# include release data)
deployment_dates <- rbind(
  chin_in %>% 
    select(vemco_code = acoustic_year, release_date = deployment_time),
  hendricks_bio %>% 
    select(vemco_code, release_date = date)
)
n_day_dat <- depth_dat2 %>% 
  select(vemco_code, date_time_local) %>% 
  group_by(vemco_code) %>% 
  summarize(
    min_time = min(date_time_local),
    max_time = max(date_time_local)
  ) %>% 
  left_join(., 
            deployment_dates,
            by = "vemco_code") %>%
  mutate(
    days_before_det = difftime(release_date, min_time, units = "days"), 
    days_at_large = difftime(release_date, max_time, units = "days")
  )

det_hist <- n_det_dat %>% 
  filter(n < 600) %>% #excludes n=17 tags
  ggplot(., aes(x = n)) +
  geom_histogram(binwidth = 10) + 
  ggsidekick::theme_sleek() +
  labs(x = "Detections") +
  theme(axis.title.y = element_blank())

day_first_det_hist <- n_day_dat %>% 
  filter(days_before_det > -100) %>% #excludes n=16
  ggplot(., 
         aes(x = -1 * as.numeric(days_before_det))) +
  geom_histogram(binwidth = 5) + 
  ggsidekick::theme_sleek() +
  labs(x = "Days Between Release and First Detection") +
  theme(axis.title.y = element_blank())

day_at_large_hist <- ggplot(n_day_dat, 
                            aes(x = -1 * as.numeric(days_at_large))) +
  geom_histogram(binwidth = 10) + 
  ggsidekick::theme_sleek() +
  labs(x = "Days Between Release and Last Detection") +
  theme(axis.title.y = element_blank())


two_panel <- cowplot::plot_grid(
  det_hist, day_first_det_hist, day_at_large_hist, ncol = 1
)

png(here::here("figs", "ms_figs_rel", "detection_histograms.png"),
    height = 7, width = 4.5, res = 250, units = "in")
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    two_panel, 
    left = grid::textGrob("Number of Tags", rot = 90, 
                          gp = grid::gpar(fontface = "bold"))
  )
)
dev.off()



## DETECTION/BATHY MAPS --------------------------------------------------------

library(maptools)
library(rmapshaper)
library(mapdata)

# map data
w_can <- map_data("worldHires", region = c("usa", "canada")) %>%
  filter(long < -110) %>%
  fortify(.)
coast_plotting <- readRDS(here::here("data",
                                     "coast_major_river_sf_plotting.RDS"))
sf::st_crs(coast_plotting) <- 4326


## Detection Maps
depth_dat2 <- readRDS(here::here("data", "depth_dat_nobin.RDS")) %>% 
  filter(!is.na(roms_temp))

# base receiver plot
base_rec_plot <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  #removes extra border
  scale_x_continuous(
    breaks = c(-127, -125, -123), 
    expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_sf(xlim = c(-128, -122.15), ylim = c(46, 51.25)) +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey")) +
  facet_wrap(~year) +
  guides(fill = guide_legend(), 
         size = guide_legend()) 


# identify receivers w/ no detections
blank_rec <- rec_full$dup_rec_all %>%
  mutate(year = as.numeric(field_season)) %>% 
  filter(
    !(station_latitude %in% depth_dat2$latitude & 
        station_longitude %in% depth_dat2$longitude),
    !project_name == "LakeWashington",
    marine == "yes",
    #exclude regions not incldued in depth analysis
    region %in% depth_dets1$region,
    # exclude SOBaD redeployments that extended to 2023 to avoid dups
    !recover_year == "2023"
  ) 

# calculate number of detections by receiver
rec_dets <- depth_dat2 %>% 
  group_by(latitude, longitude, trim_sn, year) %>% 
  summarize(
    n_dets = n(),
    .groups = "drop"
  ) 

rec_plot_det <- base_rec_plot +
  geom_point(data = blank_rec, 
             aes(x = station_longitude, y = station_latitude),
             shape = 3, 
             inherit.aes = FALSE) +
  geom_point(data = rec_dets,
             aes(x = longitude, y = latitude, size = n_dets),
             shape = 21, fill = "red", alpha = 0.4,
             inherit.aes = FALSE) +
  scale_size_continuous(
    name = "Number of \nDetections",
    breaks = c(5, 50, 150, 250, 500),
    labels = c("5", "50", "150", "250", "500")
  ) +
  theme(legend.position = "top")


# calculate mean depth by by receiver
rec_depth <- depth_dat2 %>% 
  group_by(latitude, longitude, trim_sn, year) %>% 
  summarize(
    med_depth = median(pos_depth),
    .groups = "drop"
  ) 

rec_plot_depth <- base_rec_plot + 
  geom_point(data = blank_rec, 
             aes(x = station_longitude, y = station_latitude),
             shape = 3, 
             inherit.aes = FALSE) +
  geom_point(data = rec_depth,
             aes(x = longitude, y = latitude, size = med_depth),
             shape = 21, fill = "blue", alpha = 0.2,
             inherit.aes = FALSE) +
  scale_size_continuous(
    name = "Median Depth\nof Detections (m)",
    breaks = c(5, 50, 100, 175, 250),
    labels = c("5", "50", "100", "175", "250")
    ) +
  theme(legend.position = "top")


png(here::here("figs", "ms_figs_rel", "rec_locations_det.png"),
    height = 8.5, width = 8, res = 250, units = "in")
rec_plot_det
dev.off()

png(here::here("figs", "ms_figs_rel", "rec_locations_depth.png"),
    height = 8.5, width = 8, res = 250, units = "in")
rec_plot_depth
dev.off()


## Bathymetry Maps
coast2 <- rbind(rnaturalearth::ne_states( "United States of America",
                                          returnclass = "sf"),
                rnaturalearth::ne_states( "Canada", returnclass = "sf"))
coast_utm_bathy <- coast2 %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))
coast_utm_rec <- coast2 %>% 
  sf::st_crop(., 
              xmin = -128.5, ymin = 46, xmax = -122, ymax = 51.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))

bath_grid <- readRDS(here::here("data", "pred_bathy_grid_utm.RDS")) %>% 
  mutate(id = row_number(),
         shore_dist_km = shore_dist / 1000) %>% 
  filter(depth < 400,
         Y > (5095 * 1000))

# location labels
labs_dat <- data.frame(
  X = c(439000, 431000, 510800, 300000, 240000),
  Y = c(5303000, 5465000, 5187000, 5250000, 5400000),
  lab = c("Juan de Fuca\nStrait", "Strait of\n Georgia", "Puget\nSound", 
          "Washington\nCoast", "La Perouse\nBank")
)

# sampling locations for Port Renfrew and Ukee
loc_dat <- data.frame(
  X = c(313709, 385246),
  Y = c(5421491, 5378367),
  site = c("U", "PR")
)


# bounding box representing zone of predictions
domain_coords <- sf::st_coordinates(coast_utm_rec) %>% as.data.frame()
bb_coords <- sf::st_coordinates(coast_utm_bathy) %>% as.data.frame()

# receiver location data 
rec_locs <- rec %>%
  filter(marine == "yes",
         !station_latitude < 46) %>%
  select(station_latitude, station_longitude) %>%
  distinct()
rec_locs <- sdmTMB::add_utm_columns(
  rec_locs, 
  ll_names = c("station_longitude", "station_latitude"),
  units = "m"
)


# release location data 
dfo_releases <- chin_in %>% 
  filter(acoustic_type == "V13P") %>% 
  select(vemco_code = acoustic_year, lat, lon, year)
ubc_releases <- hendricks_bio %>% 
  filter(
    # remove missing latitude and one individual w/ incorrect location
    !is.na(lat),
    !(lat > 48.6 & lon > -124.6)
  ) %>% 
  mutate(
    year = str_split(vemco_code, "_") %>% 
      purrr::map(., ~ .x[2]) %>% 
      unlist(),
    # shift latitude to account for differences in projection (moves off land)
    lat = lat - 0.01
  ) %>% 
  select(vemco_code, lat, lon, year)

releases <- rbind(dfo_releases, ubc_releases)
releases <- sdmTMB::add_utm_columns(
  releases, 
  ll_names = c("lon", "lat"),
  units = "m"
)


# blank base map
blank_p <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        legend.position = "top") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


# inset map
w_can <- map_data("world", region = c("usa", "canada")) %>%
  fortify(.) 
inset_map <- ggplot() +
  geom_polygon(data = w_can, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "grey") + 
  labs(x = "", y = "") +
  geom_rect(aes(xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5),
            fill = NA, lty = 2, col = "red") +
  coord_map(ylim = c(28, 56), xlim = c(-132, -118)) +
  ggsidekick::theme_sleek() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) 

# make plots
rec_plot <-  blank_p +
  geom_sf(data = coast_utm_rec, fill = "darkgrey") +
  geom_rect(aes(xmax = max(bb_coords$X) - 25000,
                xmin = min(bb_coords$X),
                ymax = max(bb_coords$Y) - 45000,
                ymin = min(bb_coords$Y) + 5000),
            color = "black", fill = "transparent", lty = 2) + 
  # rectangle for release locations
  geom_rect(aes(xmax = max(releases$X) + 2000,
                xmin = min(releases$X) - 2000,
                ymax = max(releases$Y) + 2000,
                ymin = min(releases$Y) - 2000),
            color = "black", fill = "transparent", lty = 3) +
  geom_label(data = labs_dat, aes(x = X, y = Y, label = lab), size = 3) +
  geom_point(data = rec_locs, aes(x = X, y = Y), 
             fill = "red", shape = 21, alpha = 0.5, size = 0.9) +
  geom_point(data = loc_dat, aes(x = X, y = Y, shape = site), fill = "blue",
             size = 1.5) +
  scale_shape_manual(values = c(23, 24), guide = NULL) +
  coord_sf(ylim = c(min(domain_coords$Y), max(domain_coords$Y) - 15000),
           xlim = c(min(domain_coords$X), max(domain_coords$X) - 15000))

rec_plot2 <- cowplot::ggdraw(rec_plot) +
  cowplot::draw_plot(inset_map, x = 0.07, y = 0.1, #vjust = -0.2,
                     hjust = 0.1,
                     width = 0.37, height = 0.37)

depth_plot <- blank_p +
  geom_raster(data = bath_grid,
              aes(x = X, y = Y, fill = depth)) +
  geom_sf(data = coast_utm_bathy, fill = "darkgrey") +
  scale_shape_manual(values = c(21, 24), guide = NULL) +
  scale_fill_viridis_c(name = "Bottom\nDepth (m)", direction = -1)  +
  theme(legend.position = "right",
        axis.text = element_blank()) +
  coord_sf(ylim = c(min(bb_coords$Y) + 10000, 5480000),
           xlim = c(min(bb_coords$X), max(bb_coords$X) - 10000))
slope_plot <- blank_p +
  geom_raster(data = bath_grid, 
              aes(x = X, y = Y, fill = slope)) +
  geom_sf(data = coast_utm_bathy, fill = "darkgrey") +
  scale_fill_viridis_c(name = "Slope\n(degrees)", direction = -1) +
  theme(legend.position = "right",
        axis.text = element_blank()) +
  coord_sf(ylim = c(min(bb_coords$Y) + 10000, 5480000),
           xlim = c(min(bb_coords$X), max(bb_coords$X) - 10000))
dist_plot <- blank_p +
  geom_raster(data = bath_grid, 
              aes(x = X, y = Y, fill = shore_dist_km)) +
  geom_sf(data = coast_utm_bathy, fill = "darkgrey") +
  scale_fill_viridis_c(name = "Distance\nto Coast\n(km)", direction = -1) +
  theme(legend.position = "right",
        axis.text = element_blank()) +
  coord_sf(ylim = c(min(bb_coords$Y) + 10000, 5480000),
           xlim = c(min(bb_coords$X), max(bb_coords$X) - 10000))


bathy_vars1 <- cowplot::plot_grid(
  plotlist = list(depth_plot, slope_plot, dist_plot),
  ncol = 1
)
bathy_vars <- cowplot::plot_grid(
  plotlist = list(rec_plot2, bathy_vars1),
  ncol = 2,
  rel_widths = c(1.5, 1)
)

png(here::here("figs", "ms_figs_rel", "study_area.png"),
    height = 145, width = 170, units = "mm", res = 300)
bathy_vars
dev.off()


## release location maps
rel_locs <- blank_p +
  geom_sf(data = coast_utm_rec, fill = "darkgrey") +
  geom_point(data = releases, aes(x = X, y = Y, fill = year), 
             shape = 21, alpha = 0.7) +
  scale_fill_viridis_d(option = "A") +
  coord_sf(ylim = c(min(releases$Y) - 2000, max(releases$Y) + 2000),
           xlim = c(min(releases$X) - 2000, max(releases$X) + 2000))

png(here::here("figs", "ms_figs_rel", "release_locations.png"),
    height = 3, width = 6.5, res = 250, units = "in")
rel_locs
dev.off()

