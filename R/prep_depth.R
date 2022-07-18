### Prep detections data for depth-specific analyses
## Jan. 13, 2022

library(tidyverse)


## helper functions 
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
              dplyr::select(receiver = receiver_name, trim_sn,
                            mean_bathy = mean_depth,
                            max_bathy = max_depth, mean_slope:shore_dist,
                            marine) %>% 
              distinct(),
            by = "receiver") %>% 
  filter(marine == "yes") %>% 
  left_join(., stage_dat %>% dplyr::rename(stage = stage_predicted), 
            by = "vemco_code") %>% 
  select(-temp, -flag, -station_name, -marine)


# Add B. Hendricks data

# sampling data 
hendricks_bio <- readRDS(here::here("data", "hendricks_chin_dat.RDS"))
  
# NOTE: differences in receiver formatting may lead to issues...
# Only correct for maturity based on date, not survival (no terminal detections
# available)
depth_combined <- readRDS(here::here("data", "hendricks_depth_dets.RDS")) %>% 
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
  dplyr::rename(receiver = receiver_name) %>% 
  select(colnames(depth_raw)) %>% 
  rbind(depth_raw, .) %>% 
  #distinguish between aggregates based on CU 
  mutate(
    agg = case_when(
      grepl("East Vancouver Island", cu_name) ~ "ECVI",
      grepl("_1.", cu_name) ~ "FraserYear",
      grepl("_0.", cu_name) ~ "FraserSub",
      cu_name == "Lower Columbia River" ~ "LowCol",
      TRUE ~ agg
    )
  )

# depth_combined %>%
#   ggplot(.) +
#   geom_jitter(aes(x = a, y = depth, colour = stage))

# NOTE CHECK PRED GRID HAS BEEN CHANGED TO ZONE 10
depth_utm <- lonlat_to_utm(depth_combined$longitude, depth_combined$latitude, 
                          zone = 10) 
depth_combined$utm_x <- depth_utm$X / 1000 
depth_combined$utm_y <- depth_utm$Y / 1000


# function to make depth_data at different bin sizes 
depth_foo <- function(bin_size = 30) {
  if (!is.null(bin_size)) {
    depth_dat <- depth_combined %>%
      #calculate timestep (width below in minutes) relative to first detection
      mutate(
        start_time = min(date_time),
        timestamp = difftime(start_time, date_time, units = "mins"),
        timestamp = -1 * round(as.numeric(timestamp)),
        timestamp_f = cut_width(timestamp, width = bin_size, boundary = -0.1)
      ) %>%
      # bin depth data by tag, receiver, and timestep
      group_by(vemco_code, timestamp_f, receiver, latitude, longitude,
               region) %>%
      # dplyr::summarize(
      dplyr::mutate(
        timestamp_n = mean(timestamp) + rnorm(1, 0, 0.01),
        date_time = mean(date_time),
        depth = mean(depth),
        bin_id = row_number()#,
        # .groups = "drop"
      ) %>%
      ungroup() %>%
      # remove redundant observations
      filter(bin_id == "1")
  } else {
    depth_dat <- depth_combined %>%
      #calculate timestep (width below in minutes) relative to first detection
      # group_by(vemco_code) %>%
      mutate(
        start_time = min(date_time),
        timestamp = difftime(start_time, date_time, units = "mins"),
        timestamp_n = -1 * round(as.numeric(timestamp))
        )
  }
   depth_dat2 <- depth_dat %>%
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
   sun_data <- data.frame(date = as.Date(depth_dat2$date_time_local),
                          lat = depth_dat2$latitude, 
                          lon = depth_dat2$longitude)
   temp <- suncalc::getSunlightTimes(data = sun_data,
                                     keep = c("sunrise", "sunset"),
                                     tz = "America/Los_Angeles")
   
   cbind(depth_dat2, temp %>% dplyr::select(sunrise, sunset)) %>%
     mutate(
       day_night = ifelse(date_time_local > sunrise & date_time_local < sunset,
                          "day", "night")
     ) %>%
     dplyr::select(
       vemco_code, cu_name, agg, fl, mean_log_e, stage, trim_sn, 
       receiver:longitude, utm_y, utm_x, mean_bathy:shore_dist,
       u, v, w, roms_temp, zoo,
       region_f, date_time_local, timestamp_n, hour, day_night,
       det_day, year, pos_depth = depth, rel_depth
     ) 
  }


depth_dat_null <- depth_foo(bin_size = NULL)
depth_dat_15 <- depth_foo(bin_size = 15)
depth_dat_60 <- depth_foo(bin_size = 60)
 
 
# export
saveRDS(depth_dat_15,
        here::here("data", "depth_dat_15min.RDS"))
saveRDS(depth_dat_60,
        here::here("data", "depth_dat_60min.RDS"))
saveRDS(depth_dat_null,
        here::here("data", "depth_dat_nobin.RDS"))


## CHECK DETS ------------------------------------------------------------------

# one receiver had > 1000 detections 
depth_dat_null %>% 
  filter(trim_sn == "108654", year == "2020") %>%
  group_by(vemco_code) %>% 
  tally()
# check tags w/ large number dets

depth_dat_null %>%
  filter(vemco_code %in% c("10123_2020", "10121_2020")) %>%
  ggplot(., aes(x = date_time_local, y = -1 * pos_depth, fill = region_f)) +
  geom_point(shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "") +
  labs(x = "Timestamp", y = "Depth (m)") +
  ggsidekick::theme_sleek() +
  facet_wrap(~vemco_code) +
  theme(legend.position = "top")
# seems reasonable...


## DETECTION MAPS --------------------------------------------------------------

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

# calculate number of detections by receiver
rec_dets <- depth_dat_null %>% 
  group_by(latitude, longitude, trim_sn, year) %>% 
  summarize(
    n_dets = n(),
    .groups = "drop"
    ) %>% 
  mutate(
    det_bin = cut(n_dets, breaks=c(0, 25, 100, 500))
  )
  
  
rec_plot <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  geom_point(data = rec_dets,
             aes(x = longitude, y = latitude, size = n_dets),
             shape = 21, fill = "red", alpha = 0.4,
             inherit.aes = FALSE) +
  #removes extra border
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_sf(xlim = c(-128, -122.15), ylim = c(46, 51.25)) +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey")) +
  facet_wrap(~year) +
  guides(fill = guide_legend(), 
         size = guide_legend()) +
  scale_size_continuous(name = "Number of \nDetections",
                        breaks = c(0, 25, 100, 500, 1000)) 
  


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