## Temperature Selection
# Oct 19, 2022
# Identify temperature at depth of detections data based on ROMS outputs from DJ


library(tidyverse)


temp_dat <- read.csv(
  here::here("data", "raw_roms", 
             "stations_roms_no_infill_tempProfile_05oct22_all.csv")) %>% 
  # rename(temp = value, latitude = lat, longitude = lon, pos_depth = depth) %>% 
  select(year, month, day, hour, latitude = lat, longitude = lon, temp = value,
         depth)# %>%
  # mutate(depth_int = round(depth, digits = -1)) %>% 
  # select(-depth)



# moderately cleaned detections data (includes depth/temperature sensors)
# receiver data (includes bathymetric data generated in chin_tagging repo)
# depth_raw <- readRDS(here::here("data", "detections_all.RDS")) %>%
#   filter(
#     !is.na(depth),
#     depth > 0,
#     # remove stations that are freshwater or terminal
#     !station_name == "LakeWashington",
#     !region == "Fraser"
#   ) %>%
#   select(vemco_code, date_time, latitude, longitude, depth)
# 
# # B. Hendricks data
# depth_h <- readRDS(here::here("data", "hendricks_depth_dets.RDS")) %>%
#   select(vemco_code, date_time, latitude, longitude, depth)

# depth_dat_raw <- rbind(depth_raw, depth_h) %>% 
#   mutate(hour = lubridate::hour(date_time) + 1,
#          day = lubridate::day(date_time),
#          month = lubridate::month(date_time),
#          year = lubridate::year(date_time)) %>% 
#   glimpse()
# 
# depth_dat_raw %>% 
#   filter(year %in% temp_dat$year & month %in% temp_dat$month & 
#            day %in% temp_dat$day & 
#          latitude %in% temp_dat$latitude & 
#          longitude %in% temp_dat$longitude) %>% 
#   glimpse()

depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  mutate(stage = as.factor(stage),
         month = lubridate::month(date_time),
         day = lubridate::day(date_time),
         hour = lubridate::hour(date_time) + 1) %>%
  select(vemco_code:stage, year, month, day, hour, date_time:utm_x, region_f,
         day_night) 

d_trim <- depth_dat_raw %>% 
  select(vemco_code, year:hour, lat = latitude, lon = longitude, region_f,
         pos_depth)
  

depth_dat <- left_join(
  temp_dat, d_trim, 
  by = c("year", "month", "day", "hour", "lat", "lon")
) %>% 
  mutate(
    depth_diff =  abs(depth - pos_depth)
  ) %>% 
  group_by(year, month, day, hour, lat, lon, pos_depth) %>% 
  mutate(
    min_depth_diff = min(depth_diff)
  ) %>% 
  ungroup() %>% 
  filter(depth_diff == min_depth_diff) %>% 
  mutate(
    depth = pmin(-1 * depth, 0),
    season = case_when(
      month %in% c("12", "1", "2") ~ "winter",
      month %in% c("3", "4", "5") ~ "spring",
      month %in% c("6", "7", "8") ~ "summer",
      month %in% c("9", "10", "11") ~ "fall"
    ),
    region_f = case_when(
      region_f %in% c("sog", "disc_isl") ~ "sog",
      region_f %in% c("brooks", "swvi") ~ "wcvi",
      grepl("wwa", region_f) ~ "wa",
      TRUE ~ as.character(region_f)
    )
  ) %>% 
  select(year:lon, temperature = value, region_f, season, depth)

depth_dat %>% 
  ggplot(.) +
  geom_point(aes(x = temperature, y = depth), alpha = 0.3, fill = "red",
             shape = 21) +
  facet_grid(region_f~season)

# as above but with bins for temp and depth
depth_seq <- seq(-340, 0, by = 20)
temp_seq <- seq(6, 24, by = 0.5)


depth_bin <- depth_dat %>% 
  filter(!region_f == "columbia") %>% 
  mutate(
    depth_bin = cut(
      depth, 
      breaks = depth_seq) %>% 
      as.factor(),
    temp_bin = cut(
      temperature, 
      breaks = temp_seq) %>% 
      as.factor()
  ) %>% 
  group_by(depth_bin, temp_bin, region_f, season) %>% 
  tally()

ggplot(depth_bin) +
  geom_raster(aes(x = temp_bin, y = depth_bin, fill = n)) +
  scale_fill_viridis_c(
    trans = "sqrt"
  ) +
  facet_grid(region_f~season) 









# tt <- depth_dat_raw[900:950, ]
# depth_dat_raw$temp <- NaN
# 
# for (i in 1:nrow(depth_dat_raw)) {
#   dum <- depth_dat_raw[i, ]
#   xx <- temp_dat %>% 
#     filter(year == dum$year, month == dum$month, day == dum$day, 
#            hour == dum$month, latitude == dum$latitude, 
#            longitude == dum$longitude)
#   # whats the nearest depth
#   temp_index <- which(
#     abs(xx$depth - dum$depth) == min(abs(xx$depth - dum$depth))
#     )
#   
#   if (length(temp_index) > 1) stop("Multiple strata selected")
#   
#   if (nrow(xx) > 0) depth_dat_raw[i , "temp"] <- xx$temp[temp_index]
# }
