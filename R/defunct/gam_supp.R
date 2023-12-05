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


depth_dat <- depth_dat_raw %>% 
  group_by(vemco_code, fl, lipid, stage_dummy, utm_x, utm_y, day_night_dummy,
           local_day, mean_bathy, mean_slope, shore_dist,
           u, v, w, roms_temp, zoo, oxygen, thermo_depth) %>% 
  summarize(
    moon_illuminated = mean(moon_illuminated),
    mean_rel_depth = mean(rel_depth),
    sd_rel_depth = sd(rel_depth),
    .groups = "drop"
  ) %>% 
  left_join(., ind_folds, by = "vemco_code")


# split by individual blocking (identical to random forest)
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()
test_depth_22 <- depth_dat_raw1 %>% 
  group_by(vemco_code, fl, lipid, stage_dummy, utm_x, utm_y, day_night_dummy,
           local_day, mean_bathy, mean_slope, shore_dist,
           u, v, w, roms_temp, zoo, oxygen, thermo_depth) %>% 
  summarize(
    moon_illuminated = mean(moon_illuminated),
    mean_rel_depth = mean(rel_depth),
    sd_rel_depth = sd(rel_depth),
    .groups = "drop"
  ) 


## FIT -------------------------------------------------------------------------


fit1 <- gam(
  mean_rel_depth ~ te(utm_x, utm_y, bs=c("tp", "tp"), k=c(10, 10)) +
    s(mean_bathy, k = 3) + 
    s(local_day, bs = "cc", k = 5) + 
    day_night_dummy + stage_dummy +# fl + lipid +
    s(vemco_code, bs = "re"),
  data = train_depth,
  knots = list(local_day = c(0, 365)),
  family = betar(link = "logit")
)

plot(fit1)



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

pred_fit <- predict(fit1, type = "response", se.fit = TRUE, newdata = pred_dat1)
colnames(pred_fit$predictions) <- c("lo", "med", "up")

pred_dat2 <- pred_dat1 %>%
  mutate(
    rel_pred_med = pred_fit$fit,
    # rel_pred_lo = pred_fit$predictions[, "lo"],
    # rel_pred_up = pred_fit$predictions[, "up"],
    # pred_med = pred_fit$predictions[, "med"] * max_bathy,
    # pred_lo = pred_fit$predictions[, "lo"] * max_bathy,
    # pred_up = pred_fit$predictions[, "up"] * max_bathy,
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000#,
    # pred_int_width = pred_up - pred_lo
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
              aes(x = utm_x_m, y = utm_y_m, fill = pred_int_width)) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "80% Prediction\nInt. Width (m)", 
                       option = "C",
                       direction = -1)  +
  theme(legend.position = "top",
        axis.text = element_blank())

avg_depth1 <- cowplot::plot_grid(plotlist = list(rel_depth, mean_depth),
                                 ncol = 2)