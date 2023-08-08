## Temperature Selection
# Oct 19, 2022
# Identify temperature at depth of detections data based on ROMS outputs from DJ


library(tidyverse)


profile_dat <- read.csv(
  here::here("data", "raw_roms", 
             "stations_roms_no_infill_profiles_10may23_all.csv")) %>%
  select(year, month, day, hour, latitude = lat, longitude = lon, variable, 
         value, depth)


depth_dat_raw <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  mutate(stage = as.factor(stage),
         month = lubridate::month(date_time_utm),
         day = lubridate::day(date_time_utm),
         hour = lubridate::hour(date_time_utm) + 1) 

d_trim <- depth_dat_raw %>% 
  select(vemco_code, year:hour, latitude, longitude, region_f, pos_depth)
  

depth_dat <- left_join(
  profile_dat, d_trim, 
  by = c("year", "month", "day", "hour", "latitude", "longitude"),
  relationship = "many-to-many"
) %>%
  mutate(
    depth_diff =  abs(depth - pos_depth)
  ) %>% 
  group_by(year, month, day, hour, latitude, longitude, pos_depth) %>% 
  mutate(
    min_depth_diff = min(depth_diff)
  ) %>% 
  ungroup() %>% 
  filter(depth_diff == min_depth_diff) %>% 
  pivot_wider(names_from = "variable", values_from = "value") %>% 
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
    ),
    oxy_mg_l = respR::convert_DO((oxygen / 1000), from = "mmol/L", to = "mg/L")
  ) %>% 
  select(year, month, day, hour, latitude, longitude, 
         temp, oxy_mg_l, region_f, season, depth)

depth_dat %>% 
  ggplot(.) +
  geom_point(aes(x = temp, y = depth), alpha = 0.3, fill = "red",
             shape = 21) +
  facet_grid(region_f~season)

depth_dat %>% 
  ggplot(.) +
  geom_point(aes(x = oxy_mg_l, y = depth), alpha = 0.3, fill = "blue",
             shape = 21) +
  facet_grid(region_f~season)


# as above but with bins for temp and depth
depth_seq <- seq(-340, 0, by = 20)
depth_labs <- as.character(depth_seq[-1])
temp_seq <- seq(6, 24, by = 0.5)
temp_labs <- as.character(temp_seq[-length(temp_seq)])
oxy_seq <- seq(1, 12, by = 0.5)
oxy_labs <- as.character(oxy_seq[-length(oxy_seq)])


depth_bin <- depth_dat %>% 
  filter(!region_f %in% c("columbia", "or_ca")) %>% 
  mutate(
    depth_bin = cut(
      depth, 
      breaks = depth_seq,
      labels = depth_labs) %>% 
      as.factor(),
    temp_bin = cut(
      temp, 
      breaks = temp_seq,
      labels = temp_labs) %>% 
      as.factor(),
    oxy_bin = cut(
      oxy_mg_l, 
      breaks = oxy_seq,
      labels = oxy_labs) %>% 
      as.factor(),
    season = fct_relevel(as.factor(season), 
                         "spring", "summer", "fall", "winter"),
    region_f = factor(region_f, 
                         levels = c("wcvi", "wa", "jdf", "sog", "puget"),
                      labels = c("WCVI", "WA", "JdF", "SoG", "PS"))
  ) 

temp_bin_plot <- depth_bin %>% 
  group_by(depth_bin, temp_bin, region_f, season) %>% 
  tally() %>% 
  ggplot(.) +
  geom_raster(aes(x = temp_bin, y = depth_bin, fill = n)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Number of\nDetections",
    breaks = c(50, 250, 450, 650, 850)
  ) +
  facet_grid(region_f~season) +
  ggsidekick::theme_sleek() +
  scale_y_discrete(name = "Depth (m)", 
                   breaks = seq(-300, 0, by = 25)) +
  scale_x_discrete(name = "Temperature (C)", 
                   breaks = seq(6, 24, by = 2)) +
  theme(axis.text.x = element_text(angle = 45))

png(here::here("figs", "ms_figs_rel", "temp_at_depth.png"),
    units = "in", res = 200, 
    width = 9, height = 6)
temp_bin_plot
dev.off()


oxy_bin_plot <- depth_bin %>% 
  group_by(depth_bin, oxy_bin, region_f, season) %>% 
  tally() %>% 
  ggplot(.) +
  geom_raster(aes(x = oxy_bin, y = depth_bin, fill = n)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Number of\nDetections",
    option = "A",
    breaks = c(50, 250, 450, 650, 850)
  ) +
  facet_grid(region_f~season) +
  ggsidekick::theme_sleek() +
  scale_y_discrete(name = "Depth (m)", 
                   breaks = seq(-300, 0, by = 25)) +
  scale_x_discrete(name = "Oxygen (mg/l)", 
                   breaks = seq(1, 12, by = 1)) +
  theme(axis.text.x = element_text(angle = 45))

png(here::here("figs", "ms_figs_rel", "oxy_at_depth.png"),
    units = "in", res = 200, 
    width = 9, height = 6)
oxy_bin_plot
dev.off()
