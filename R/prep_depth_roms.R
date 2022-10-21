### ROMS Data Clean
## Export stations for query, then import and evaluate 
## Query conducted by Doug Jackson (doug@qedaconsulting.com)
## Feb 14, 2022

library(tidyverse)


## PREP DEPTH DATA -------------------------------------------------------------

# receiver data (includes bathymetric data generated in chin_tagging repo)
rec <- readRDS(here::here("data", "receivers_all.RDS"))$rec_all 

# life stage estimates
dum <- readRDS(here::here("data", "acousticOnly_GSI.RDS"))
stage_dat <- readRDS(here::here("data", "lifestage_df.RDS")) %>% 
  left_join(., 
            dum %>% dplyr::select(vemco_code = acoustic_year, mean_log_e, 
                                  cu_name), 
            by = "vemco_code") %>% 
  select(vemco_code, fl, stage_predicted, mean_log_e, cu_name, agg)  

# ~2% of tags lack energy density estimates, interpolate
interp_stage_dat <- VIM::kNN(stage_dat, k = 5) %>% 
  select(-ends_with("imp")) 



# moderately cleaned detections data (includes depth/temperature sensors)
depth_raw <- readRDS(here::here("data", "detections_all.RDS")) %>%
  filter(
    !is.na(depth),
    depth > 0,
    # remove stations that are freshwater or terminal
    !station_name == "LakeWashington",
    !region == "Fraser"
  ) %>% 
  left_join(., interp_stage_dat %>% dplyr::rename(stage = stage_predicted), 
            by = "vemco_code") %>% 
  left_join(., 
            rec %>% 
              select(receiver = receiver_name, trim_sn, mean_depth:shore_dist),
            by = "receiver") %>% 
  select(vemco_code, stage, fl, mean_log_e, cu_name, agg, date_time, latitude, 
         longitude, receiver_name = receiver, trim_sn, region, 
         mean_bathy = mean_depth, max_bathy = max_depth, mean_slope, shore_dist,
         depth)


# B. Hendricks data
hendricks_bio <- readRDS(here::here("data", "hendricks_chin_dat.RDS"))

depth_h <- readRDS(here::here("data", "hendricks_depth_dets.RDS")) %>% 
  mutate(trim_sn = as.character(receiver_sn),
         week = lubridate::week(date_time)) %>%
  group_by(vemco_code) %>%
  mutate(
    max_week = max(week),
    stage = case_when(
      max_week > 40 ~ "immature",
      TRUE ~ "mature")
  ) %>% 
  ungroup() %>% 
  left_join(., 
            hendricks_bio %>% 
              select(vemco_code, fl, mean_log_e, cu_name, agg = agg_name), 
            by = "vemco_code") %>%
  select(colnames(depth_raw)) 

# combine
depth_dets1 <- rbind(depth_raw, depth_h) %>% 
  mutate(
    hour = lubridate::hour(date_time) + 1,
    day = lubridate::day(date_time),
    month = lubridate::month(date_time),
    year = lubridate::year(date_time),
    agg = case_when(
      grepl("East Vancouver Island", cu_name) ~ "ECVI",
      grepl("_1.", cu_name) ~ "FraserYear",
      grepl("_0.", cu_name) ~ "FraserSub",
      cu_name == "Lower Columbia River" ~ "LowCol",
      TRUE ~ agg
    )
  ) %>% 
  filter(!region == "fraser")


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


## EXPORT STATIONS -------------------------------------------------------------

# expand depth dets by different covariates and depths
trim_dets <- depth_dets1 %>% 
  select(year, month, day, hour, lat = latitude, lon = longitude) %>% 
  distinct()
var_list <- c("u", "v", "w", "temp", "zooplankton", "phytoplankton", "oxygen")
# focus on one depth given strong correlations among them
depth_list <- c(#5, 
  25#, 50
  )

roms_dat_out <- expand.grid(
  variable = var_list,
  depth = depth_list
) %>% 
  mutate(
    # depth = ifelse(variable == "w", -999, depth),
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


# roms_dat_out %>% filter(
#   hour == "19", month == "9", day == "15", year == "2020",
#   lat == missing_roms$latitude[1], lon == missing_roms$longitude[1] 
# )
# roms_dat %>% filter(
#   hour %in% missing_roms$hour_int, 
#   month %in% lubridate::month(missing_roms$date_time), 
#   day %in% lubridate::day(missing_roms$date_time), 
#   year %in% lubridate::year(missing_roms$date_time),
#   lat %in% missing_roms$latitude, lon %in% missing_roms$longitude 
# )


# export 
write.csv(roms_dat_out, here::here("data", "stations_roms_no_infill.csv"),
          row.names = FALSE)


## IMPORT TEMPERATURE DATA -----------------------------------------------------

temp_dat <- read.csv(
  here::here("data", "raw_roms", 
             "stations_roms_no_infill_tempProfile_05oct22_all.csv")) %>%
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
  ) 

# calculate thermocline depth at each station
thermo_dat<- temp_dat %>% 
  group_by(unique_id, lat, lon, timestamp) %>% 
  summarize(
    thermo_depth = rLakeAnalyzer::thermo.depth(wtr = value, depths = depth),
    .groups = "drop"
  )

# visualize
# missed_dat <- thermo_dat %>% filter(is.na(thermo_depth))
# ids <- sample(unique(thermo_dat$unique_id), 5)
# ggplot(temp_dat %>% 
#          filter(unique_id %in% c("207", "213", "241", "478"))) +
#   geom_point(aes(x = value, y = (-1 * depth), colour = as.factor(date))) #+
  # geom_hline(aes(yintercept = -1 * thermo_depth, colour = as.factor(date)))


## IMPORT ROMS PULL ------------------------------------------------------------

# query generated by Doug; note some stations have missing values due to outside
# survey domain
# Feb pull has multiple depths, but latter do not because highly correlated 
# with one another
roms_dat <- read.csv(
  # here::here("data", "raw_roms", "stations_roms_no_infill_25mar22_all.csv")) %>%
  
  here::here("data", "raw_roms", "stations_roms_no_infill_05oct22_all.csv")) %>%
  distinct()


# histogram of values by depth and variable
# ggplot(roms_dat) +
#   geom_histogram(aes(x = value)) +
#   facet_wrap(~variable, scales = "free") +
#   ggsidekick::theme_sleek()
# 
# 
# # evaluate correlations between depth strata for each variable
# wide_roms <- roms_dat %>%
#   pivot_wider(names_from = "depth", names_prefix = "depth_",
#               values_from = "value") %>%
#   filter(!is.na(depth_25))
# 
# corr_list <- wide_roms %>%
#   split(., .$variable) %>%
#   purrr::map2(., names(.), function (x, y) {
#     corr <- cor(x %>% select(depth_5, depth_25, depth_50))
#     ggcorrplot::ggcorrplot(corr) +
#       labs(title = y)
#   })
# # u and v highly correlated through depth range; w not correlated at all;
# # zooplankton shows correlations at 25 m; use 25 for now


roms_25 <- roms_dat %>%
  left_join(., 
            rec %>% 
              select(lat = station_latitude, lon = station_longitude,
                     region, mean_depth, mean_slope, 
                     shore_dist ) %>% 
              distinct(),
            by = c("lat", "lon")) %>% 
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
  left_join(., thermo_dat, by = c("lat", "lon", "timestamp")) %>% 
  select(-depth, -depthFrac, -c(date:time1)) %>% 
  # remove duplicates introduced by receivers being deployed at identical 
  # locations
  distinct() %>% 
  rename(zoo = zooplankton, latitude = lat, longitude = lon,
         roms_temp = temp) %>% 
  # remove phyto because highly correlated with zoo
  select(-phytoplankton)

depth_utm <- lonlat_to_utm(roms_25$longitude, roms_25$latitude, 
                           zone = 10) 
roms_25$utm_x <- depth_utm$X / 1000 
roms_25$utm_y <- depth_utm$Y / 1000

# replace unknown values with NA
roms_25[roms_25 == "-999"] <- NA

# identify proportion fo NA by each ROMS variable
roms_25 %>% 
  pivot_longer(., 
               cols = c(roms_temp, w, v, oxygen, u, zoo, thermo_depth), 
               names_to = "roms_var") %>% 
  filter(is.na(value)) %>% 
  group_by(roms_var) %>% 
  tally() %>% 
  mutate(n / nrow(roms_25)) 


# plot locations of missing stations
# coast <- readRDS(here::here("data",
#                             "coast_major_river_sf_plotting.RDS")) %>% 
#   sf::st_crop(., 
#               xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5)
# 
# plot_missing_foo <- function(var = "roms_temp") {
#   dum <- roms_25 %>% 
#     filter(is.na(UQ(sym(var)))) %>% 
#     select(month, latitude, longitude, {{ var }}) %>% 
#     distinct() 
# 
#   ggplot() +
#     geom_sf(data = coast) +
#     geom_point(data = dum, aes(x = longitude, y = latitude), shape = 21, 
#                colour = "red") +
#     ggsidekick::theme_sleek()
# }
# 
# plot_missing_foo(var = "u")
# plot_missing_foo(var = "thermo_depth") +
#   facet_wrap(~month)
# 
# # check correlations among roms vars
# corr <- cor(roms_25 %>% 
#               select(roms_temp, w, v, oxygen, u, zoo, thermo_depth),
#             use = "pairwise.complete.obs")
# ggcorrplot::ggcorrplot(corr) 


# missing values due to a) unavailable data and b) delays between ROMS pulls

# join with receiver data (includes bathymetric variables) and impute using
# k nearest neighbors
# NOTE not all current ROMS pulls have associated receiver data but this will be 
# corrected by next extraction
roms_interp <- VIM::kNN(roms_25,
                        variable = c("roms_temp", "w", "v", "oxygen", "u", "zoo", 
                                     "thermo_depth"),
                        dist_var = c("timestamp", "utm_x", "utm_y",
                                     "mean_depth", "mean_slope", "shore_dist"),
                        k = 5)

# remove flags and non-ROMS data
roms_interp_trim <- roms_interp %>% 
  select(year:longitude, roms_temp, w, v, oxygen, u, zoo, thermo_depth) 


saveRDS(roms_interp_trim, 
        here::here("data", "interp_roms_25m_depth.RDS"))


## CLEAN DEPTH DATA ------------------------------------------------------------

roms_interp_trim <- readRDS(here::here("data", "interp_roms_25m_depth.RDS"))

depth_utm <- lonlat_to_utm(depth_dets1$longitude, depth_dets1$latitude, 
                           zone = 10) 
depth_dets1$utm_x <- depth_utm$X / 1000 
depth_dets1$utm_y <- depth_utm$Y / 1000


# add roms data (approximately 4.6k detections have no ROMS data)
depth_dat <- depth_dets1 %>% 
  left_join(
    ., roms_interp_trim, 
    by = c("latitude", "longitude", "day", "month", "year", "hour")
  ) %>% 
  # glimpse()
  mutate(
    start_time = min(date_time),
    timestamp = difftime(start_time, date_time, units = "mins"),
    timestamp_n = -1 * round(as.numeric(timestamp)),
    region_f = fct_relevel(as.factor(region), "swvi",  "nwwa", "jdf", 
                           "swwa", "sog", "puget", "columbia"),
    # corrections for when depth is above the surface or deeper than max depth
    # in detection radius
    depth = ifelse(depth <= 0.05, 0.05, depth),
    depth_diff = max_bathy - depth,
    # adjust modest errors in depth relative to bathy
    depth = ifelse(depth_diff > -10 & depth_diff < 0, max_bathy - 0.5, depth),
    # rel_depth = depth / max_bathy,
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
    stage = ifelse(
      stage == "immature" & year > tag_year & local_day > 120,
      "mature",
      stage
    )
  )  %>% 
  filter(
    # remove large errors in depth relative to bottom bathymetry
    depth < max_bathy
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
    day_night = ifelse(date_time_local > sunrise & date_time_local < sunset,
                       "day", "night"),
    moon_illuminated = moon_data
  ) %>%
  dplyr::select(
    vemco_code, cu_name, agg, fl, mean_log_e, stage, trim_sn, 
    receiver_name, latitude, longitude, utm_y, utm_x, mean_bathy, mean_slope, 
    shore_dist, u, v, w, roms_temp, zoo, oxygen, thermo_depth,
    region_f, date_time_local, timestamp_n, local_day, year, 
    local_hour, day_night, moon_illuminated, pos_depth = depth) 



# depth_dat %>%
#     select_if(function(x) any(is.na(x))) %>%
#     summarise_each(funs(sum(is.na(.))))

saveRDS(depth_dat2, here::here("data", "depth_dat_nobin.RDS"))



## CHECK DETS ------------------------------------------------------------------

# one receiver had > 1000 detections 
depth_dat2 %>% 
  filter(trim_sn == "108654", year == "2020") %>%
  group_by(vemco_code) %>% 
  tally()
# check tags w/ large number dets

depth_dat2 %>%
  filter(vemco_code %in% c("10123_2020", "10121_2020")) %>%
  ggplot(., aes(x = date_time_local, y = -1 * pos_depth, fill = region_f)) +
  geom_point(shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "") +
  labs(x = "Timestamp", y = "Depth (m)") +
  ggsidekick::theme_sleek() +
  facet_wrap(~vemco_code) +
  theme(legend.position = "top")
# seems reasonable...


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

# detections cleaned above
depth_dat2 <- readRDS(here::here("data", "depth_dat_nobin.RDS"))


# base receiver plot
base_rec_plot <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  #removes extra border
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_sf(xlim = c(-128, -122.15), ylim = c(46, 51.25)) +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey"),
        legend.position = "top") +
  facet_wrap(~year) +
  guides(fill = guide_legend(), 
         size = guide_legend()) 


# calculate number of detections by receiver
rec_dets <- depth_dat2 %>% 
  group_by(latitude, longitude, trim_sn, year) %>% 
  summarize(
    n_dets = n(),
    .groups = "drop"
  ) %>% 
  mutate(
    det_bin = cut(n_dets, breaks=c(0, 25, 100, 500))
  )

rec_plot <- base_rec_plot +
  geom_point(data = rec_dets,
             aes(x = longitude, y = latitude, size = n_dets),
             shape = 21, fill = "red", alpha = 0.4,
             inherit.aes = FALSE) +
  scale_size_continuous(name = "Number of \nDetections",
                        breaks = c(0, 25, 100, 500, 1000)) 

png(here::here("figs", "ms_figs", "rec_locations.png"),
    height = 5, width = 7, res = 250, units = "in")
rec_plot
dev.off()


# calculate mean depth by by receiver
rec_depth <- depth_dat2 %>% 
  group_by(latitude, longitude, trim_sn, year) %>% 
  summarize(
    mean_depth = mean(pos_depth),
    .groups = "drop"
  ) %>% 
  mutate(
    depth_bin = cut(mean_depth, breaks=c(0, 25, 50, 75, 100, 150, 250))
  )

rec_plot_depth <- base_rec_plot + 
  geom_point(data = rec_depth,
             aes(x = longitude, y = latitude, size = mean_depth),
             shape = 21, fill = "blue", alpha = 0.2,
             inherit.aes = FALSE) +
  scale_size_continuous(name = "Mean Depth of \nDetections",
                        breaks = c(0, 25, 50, 75, 100, 150, 250)) 


png(here::here("figs", "ms_figs", "rec_locations_det.png"),
    height = 5, width = 7, res = 250, units = "in")
rec_plot
dev.off()

png(here::here("figs", "ms_figs", "rec_locations_depth.png"),
    height = 5, width = 7, res = 250, units = "in")
rec_plot_depth
dev.off()


# utm grid for bathymetry maps
coast_utm <- coast_plotting %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


bath_grid <- readRDS(here::here("data", "pred_bathy_grid_utm.RDS")) %>% 
  mutate(id = row_number(),
         shore_dist_km = shore_dist / 1000) %>% 
  filter(depth < 400)

blank_p <- ggplot() + 
  geom_sf(data = coast_utm) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        legend.position = "top") 

depth_plot <- blank_p +
  geom_raster(data = bath_grid, 
              aes(x = X, y = Y, fill = depth)) +
  scale_fill_viridis_c(name = "Depth (m)") 
slope_plot <- blank_p +
  geom_raster(data = bath_grid, 
              aes(x = X, y = Y, fill = slope)) +
  scale_fill_viridis_c(name = "Slope\n(degrees)") 
dist_plot <- blank_p +
  geom_raster(data = bath_grid, 
              aes(x = X, y = Y, fill = shore_dist_km)) +
  scale_fill_viridis_c(name = "Distance\nto Coast (km)") 

bathy_vars <- cowplot::plot_grid(
  plotlist = list(depth_plot, slope_plot, dist_plot),
  nrow = 1
)


png(here::here("figs", "ms_figs", "bathy_preds.png"),
    height = 5, width = 8, res = 250, units = "in")
bathy_vars
dev.off()


## DEPTH PROFILES --------------------------------------------------------------

# individual depth distributions by time and terminal location
depth_dat2 <- depth_raw %>%
  mutate(date_time_local = lubridate::with_tz(date_time, 
                                              tzone = "America/Los_Angeles"),
         year = lubridate::year(date_time_local),
         n_det = length(unique(date_time)),
         fl_code = as.factor(paste(fl, vemco_code, sep = "_")),
         plot_group = case_when(
           stage == "immature" ~ "immature",
           TRUE ~ paste(agg, year, sep = "_")
         )) %>%
  filter(!n_det < 10,
         !is.na(agg)) %>%
  ungroup()
depth_list <- split(depth_dat2, depth_dat2$plot_group)


route_pal <- RColorBrewer::brewer.pal(length(unique(depth_dat2$region)),
                                      "Spectral")
names(route_pal) <- levels(fct_rev(depth_dat2$region))


trim_depth <- depth_dat2 %>%
  filter(vemco_code %in% c("7703_2019", "7707_2019", "7708_2019", "7696_2019",
                           "9969_2020", "10017_2020", "7700_2019", "9986_2020",
                           "7921_2019")) %>%
  ggplot(., aes(x = date_time_local, y = -1 * depth, fill = region)) +
  geom_point(shape = 21, alpha = 0.4) +
  scale_fill_manual(values = route_pal, name = "") +
  labs(x = "Timestamp", y = "Depth (m)") +
  ggsidekick::theme_sleek() +
  facet_wrap(~fct_reorder(vemco_code, desc(as.numeric(as.factor(stage))))) +
  theme(legend.position = "top")

png(here::here("figs", "trim_ind_profiles.png"), width = 8, height = 5,
    res = 200, units = "in")
trim_depth
dev.off()


# full plot list of absolute depth
depth_plots <- map2(depth_list, names(depth_list), .f = function(x, tit) {
  ggplot(x, aes(x = date_time_local, y = -1 * depth, fill = region)) +
    geom_point(shape = 21, alpha = 0.4) +
    scale_fill_manual(values = route_pal, name = "") +
    labs(title = tit, x = "Timestamp", y = "Depth (m)") +
    lims(y = c(-1 * max(depth_dat2$depth), 0)) +
    ggsidekick::theme_sleek() +
    facet_wrap(~fct_reorder(fl_code, desc(fl)))
})

pdf(here::here("figs", "ind_profiles.pdf"))
depth_plots
dev.off()

