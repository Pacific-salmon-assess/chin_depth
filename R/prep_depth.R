### Prep detections data for depth-specific analyses
## Jan. 13, 2022

library(tidyverse)

# receiver data (includes bathymetric data generated in chin_tagging repo)
rec <- readRDS(here::here("data", "receivers_all.RDS"))$rec_all 

# life stage estimates
stage_dat <- readRDS(here::here("data", "lifestage_df.RDS"))

# biological data
chin_dat <- readRDS(
  here::here("data", "acousticOnly_GSI.RDS")) %>% 
  filter(!grepl("MAGNET", comment)) %>%
  dplyr::rename(vemco_code = acoustic_year) %>% 
  left_join(
    ., 
    stage_dat %>% dplyr::select(vemco_code, stage), by = "vemco_code"
  ) %>% 
  mutate(
    stage = ifelse(stage == "unknown", "mature", stage)
  )


# moderately cleaned detections data (includes depth/temperature sensors)
depth_raw <- readRDS(here::here("data", "detections_all.RDS")) %>%
  filter(!is.na(depth),
         depth > 0)


# hourly ROMS outputs matched to receiver stations (marine only and some 
# missing), restricted to 25m strata
roms_dat <- readRDS(here::here("data", "roms_25m_depth.RDS"))



# function to make hours continuous
time_foo <- function(x) {
  lubridate::hour(x) + (lubridate::minute(x) / 60) + 
    (lubridate::second(x) / 3600) 
}


# function to make depth_data at different bin sizes 
depth_foo <- function(bin_size = 30) {
  if (!is.null(bin_size)) {
    depth_dat <- depth_raw %>% 
      #calculate timestep (width below in minutes) relative to first detection 
      # group_by(vemco_code) %>% 
      mutate(
        start_time = min(date_time),
        timestamp = difftime(start_time, date_time, units = "mins"),
        timestamp = -1 * round(as.numeric(timestamp)),
        timestamp_f = cut_width(timestamp, width = bin_size, boundary = -0.1)
      ) %>%
      # bin depth data by tag, receiver, and hour within a day
      group_by(vemco_code, timestamp_f, receiver, latitude, longitude, 
               station_name, region) %>%
      dplyr::summarize(
        timestamp_n = mean(timestamp) + rnorm(1, 0, 0.01), 
        date_time = mean(date_time),
        depth = mean(depth),
        .groups = "drop"
      ) %>%
      ungroup()
  } else {
    depth_dat <- depth_raw %>% 
      #calculate timestep (width below in minutes) relative to first detection 
      # group_by(vemco_code) %>% 
      mutate(
        start_time = min(date_time),
        timestamp = difftime(start_time, date_time, units = "mins"),
        timestamp_n = -1 * round(as.numeric(timestamp))
        )
  }
   depth_dat2 <- depth_dat %>%
    left_join(., 
              rec %>% 
                dplyr::select(receiver = receiver_name, mean_bathy = mean_depth,
                              max_bathy = max_depth, mean_slope:shore_dist,
                              marine) %>% 
                distinct(),
              by = "receiver") %>%
    mutate(
      hour_int = lubridate::hour(date_time) + 1,
      day = lubridate::day(date_time),
      month = lubridate::month(date_time),
      year = lubridate::year(date_time),
      region_f = as.factor(region),
      region_f = fct_relevel(region_f, "swvi",  "nwwa", "jdf", 
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
     left_join(., chin_dat %>% select(-year), by = "vemco_code") %>%  
     filter(
       # remove large errors in depth relative to bottom bathymetry
       depth < max_bathy,
       # remove regions with suspect depth estimates
       marine == "yes"
     ) %>% 
     droplevels()
  
  # %>% 
    # mutate(# calculate expected fw entry date (Oct 18, i.e. max obs day for mat 
    #        # fish + 1 day)
    #        # exit_day = ifelse(stage == "mature", 292, 292 + 365),
    #        # days_to_exit = exit_day - det_day
    # ) #%>% 
  # add terminal distance
  # left_join(.,
  #           term_dist %>%
  #             dplyr::select(receiver = rec_id, total_dist, stock, cu,
  #                           agg = agg_name),
  #           by = c("receiver", "stock", "cu", "agg")) 
  
  
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
       vemco_code, stage, receiver:longitude, mean_bathy:shore_dist, u, v, w, 
       zoo, region_f, date_time_local, timestamp_n, hour, day_night,
       det_day, year, pos_depth = depth, rel_depth
     ) 
  }

depth_dat_null <- depth_foo(bin_size = NULL)
depth_dat_15 <- depth_foo(bin_size = 30)
depth_dat_60 <- depth_foo(bin_size = 60)


# export
saveRDS(depth_dat_null,
        here::here("data", "depth_dat_nobin.RDS"))
saveRDS(depth_dat_15,
        here::here("data", "depth_dat_15min.RDS"))
saveRDS(depth_dat_60,
        here::here("data", "depth_dat_60min.RDS"))



## EXPLORATORY FIGS ------------------------------------------------------------

# depth_dat <- det_dat %>% 
#   filter(grepl("V13P", acoustic_type),
#          !is.na(depth))  %>% 
#   left_join(., 
#             rec %>% 
#               dplyr::select(receiver = receiver_name, mean_bathy, max_bathy, sd_bathy),
#             by = "receiver") %>% 
#   mutate(
#     region_f = as.factor(region),
#     region_f = fct_relevel(region_f, "swvi",  "nwwa", "jdf", "swwa", "sog", 
#                            "puget", "fraser"),
#     depth = -1 * depth,
#     depth = ifelse(depth > 0, runif(1, min = -2, max = -0.01), depth),
#     depth_diff = depth - max_bathy,
#     depth = ifelse(depth_diff > -10 & depth_diff < 0, max_bathy, depth),
#     rel_depth = depth / max_bathy,
#     date_time_local = lubridate::with_tz(date_time, 
#                                          tzone = "America/Los_Angeles"),
#     hour = lubridate::hour(date_time_local),
#     day_night = ifelse(hour < 6 | hour > 19, "night", "day"),
#     det_day = lubridate::yday(date_time_local)
#   ) %>% 
#   filter(!depth < max_bathy)
# 
# 
# # day night contrasts by different areas
# ggplot(depth_dat) +
#   geom_boxplot(aes(x = day_night, y = depth, fill = stage2)) +
#   facet_wrap(~agg) +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~station_name, scales = "free")
# 
# 
# #calculate mean depth for each individual
# mean_depth <- depth_dat %>% 
#   group_by(vemco_code, fl, mean_log_e, clip, year_day, cu, agg, det_day, 
#            stage) %>% 
#   summarize(
#     n_depth_events = length(date_time_local), 
#     mean_depth = mean(depth),
#     se_depth = sd(depth) / n_depth_events,
#     ci_low = mean_depth + (qnorm(0.025) * se_depth),
#     ci_up = mean_depth + (qnorm(0.975) * se_depth),
#     .groups = "drop"
#   ) 
# 
# ggplot(mean_depth, aes(x = det_day, y = mean_depth, 
#                        fill = n_depth_events)) + 
#   geom_point(shape = 21, size = 2) +
#   geom_pointrange(aes(ymin = ci_low, ymax = ci_up), shape = 21) +
#   labs(y = "Mean Depth (m)",
#        x = "Detection Date (year-day)") +
#   scale_fill_viridis_c() +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~stage)
# 
# 
# # individual depth distributions by time and terminal location
# depth_dat2 <- depth_dat %>% 
#   mutate(n_det = length(unique(date_time_local)),
#          fl_code = as.factor(paste(fl, vemco_code, sep = "_")),
#          plot_group = case_when(
#            stage == "immature" ~ "immature",
#            TRUE ~ paste(agg, year, sep = "_")
#          )) %>% 
#   filter(!n_det < 10) %>% 
#   ungroup()
# 
# route_pal <- RColorBrewer::brewer.pal(length(levels(depth_dat2$region_f)),
#                                       "Spectral")
# names(route_pal) <- levels(fct_rev(depth_dat2$region_f))
# 
# depth_list <- split(depth_dat2, depth_dat2$plot_group)
# 
# # absolute depth
# depth_plots <- map2(depth_list, names(depth_list), .f = function(x, tit) {
#   ggplot(x, aes(x = date_time_local, y = depth, fill = region_f)) +
#     geom_point(shape = 21) +
#     scale_fill_manual(values = route_pal, name = "") +
#     labs(title = tit, x = "Timestamp", y = "Depth (m)") +
#     lims(y = c(min(depth_dat2$depth), 0)) +
#     ggsidekick::theme_sleek() +
#     facet_wrap(~fct_reorder(fl_code, desc(fl)))
# })
# 
# pdf(here::here("figs", "depth", "ind_profiles.pdf"))
# depth_plots
# dev.off()
# 
# # relative depth
# rel_depth_plots <- map2(depth_list, names(depth_list), .f = function(x, tit) {
#   ggplot(x, aes(x = date_time_local, y = rel_depth, fill = region_f)) +
#     geom_point(shape = 21) +
#     scale_fill_manual(values = route_pal, name = "") +
#     labs(title = tit, x = "Timestamp", y = "Depth (m)") +
#     lims(y = c(0, 1)) +
#     ggsidekick::theme_sleek() +
#     facet_wrap(~fct_reorder(fl_code, desc(fl)))
# })
# 
# pdf(here::here("figs", "depth", "ind_profiles_rel_depth.pdf"))
# rel_depth_plots
# dev.off()
# 
# 
# ## trimmed subset contrasting mature vs immature
# trim_depth <- depth_dat2 %>% 
#   filter(vemco_code %in% c("7703_2019", "7707_2019", "7708_2019", "7696_2019",
#                            "9969_2020", "10017_2020", "7700_2019", "9986_2020", 
#                            "7921_2019")) %>% 
#   ggplot(., aes(x = date_time_local, y = depth, fill = region_f)) +
#   geom_point(shape = 21) +
#   scale_fill_manual(values = route_pal, name = "") +
#   labs(x = "Timestamp", y = "Depth (m)") +
#   lims(y = c(min(depth_dat2$depth), 0)) +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~fct_reorder(vemco_code, desc(fl))) +
#   theme(legend.position = "top")
# 
# 
# pdf(here::here("figs", "depth", "ind_profiles_trim.pdf"))
# trim_depth
# dev.off()
# 
# 
# # bathymetry effects
# ggplot(depth_dat %>% filter(!is.na(mean_bathy), !stage2 == "immature"), 
#        aes(x = mean_bathy, y = depth, fill = region_f)) +
#   geom_point(shape = 21) +
#   scale_fill_manual(values = route_pal, name = "") +
#   labs(x = "Mean Bottom Depth", y = "Fish Depth (m)") +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~region_f) +
#   theme(legend.position = "top")
# 
# ggplot(depth_dat %>% filter(!is.na(mean_bathy), !stage2 == "immature"), 
#        aes(x = sd_bathy, y = depth, fill = region_f)) +
#   geom_point(shape = 21) +
#   scale_fill_manual(values = route_pal, name = "") +
#   labs(x = "SD Bottom Depth", y = "Fish Depth (m)") +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~region_f) +
#   theme(legend.position = "top")
