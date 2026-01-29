## Depth Predictions for SRKW-Chinook-Rec Fishery CSAS
# Use fitted model (rel_depth_ranger_rf) to generate predicted depths for study
# domains in agent-based model
# Jan 29, 2026


library(tidyverse)


ranger_rf <- readRDS(here::here("data", "model_fits", "relative_rf_ranger.rds"))



## PREP NEW DATA ---------------------------------------------------------------

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m"))


# calculate thermocline depth for different seasons (using monthly averages)
mean_thermocline <- readRDS(here::here("data", "depth_dat_nobin.RDS")) %>%
  mutate(month = lubridate::month(date_time_local)) %>%
  filter(month %in% c(7)) %>%
  group_by(month) %>%
  dplyr::summarize(
    thermo_depth = mean(thermo_depth, na.rm = T)
  )


# input grid includes ROMS data for entire study area for specified dates 
# yday = 46 or 211
bath_grid_in <- readRDS(here::here("data", "pred_bathy_grid_roms.RDS")) 
bath_grid <- bath_grid_in %>% 
  filter(!mean_bathy > 400,
         !max_bathy > 500,
         utm_y > 5100, 
         season == "summer")


pts <- data.frame(
  subarea	= c("21A", "20D", "29DE"),
  lon = c(-124.928, -123.753, -123.36),
  lat = c(48.66694, 48.33989, 49.11563)
)

dist_m <- function(lon, lat, lon0, lat0) {
  R <- 6371000
  to_rad <- pi/180
  x <- (lon - lon0) * to_rad * cos((lat + lat0) * to_rad / 2)
  y <- (lat - lat0) * to_rad
  R * sqrt(x^2 + y^2)
}

nearest_idx <- sapply(seq_len(nrow(pts)), function(k) {
  d <- dist_m(bath_grid$lon, bath_grid$lat, pts$lon[k], pts$lat[k])
  which.min(d)
})

nearest_cells <- cbind(
  transform(pts, pt_id = seq_len(nrow(pts))),
  bath_grid[nearest_idx, setdiff(names(bath_grid), c("lat","lon")), drop = FALSE],
  lat_nn = bath_grid$lat[nearest_idx],
  lon_nn = bath_grid$lon[nearest_idx]
) %>% 
  mutate(
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000,
    thermo_depth = mean_thermocline$thermo_depth[1]
  )

nearest_cells


base_plot <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  # set limitsto avoid boundary effects from UTM projection
  scale_x_continuous(
    limits = c(210000, 560000), 
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(5100000, 5470000), expand = c(0, 0))



base_plot +
  geom_point(data = nearest_cells,
              aes(x = utm_x * 1000, y = utm_y * 1000), color = "red") +
  geom_sf(data = coast_utm)
  

## first set are average spatial predictions (i.e. mean biological and temporal 
## attributes)

# biological data
bio_dat <- depth_dat_raw %>% 
  # filter(med_stage %in% c("0", "1")) %>%
  select(vemco_code, fl, lipid) %>% 
  distinct()

# stratify predictions by non-spatial covariates
pred_dat1 <- nearest_cells %>% 
  mutate(
    fl = mean(bio_dat$fl),
    lipid = mean(bio_dat$lipid),
    med_stage = 1
  ) 

pred_rf1 <- predict(ranger_rf,
                    type = "quantiles",
                    quantiles = c(0.025, 0.5, 0.975),
                    data = pred_dat1)
colnames(pred_rf1$predictions) <- c("lo", "med", "up")



# Function to fit Beta distribution to quantiles
fit_beta_to_quantiles <- function(quantile_values, quantile_probs) {
  
  # Objective: minimize sum of squared errors
  objective <- function(params) {
    shape1 <- exp(params[1])  # α (ensure positive)
    shape2 <- exp(params[2])  # β (ensure positive)
    
    predicted_quantiles <- qbeta(quantile_probs, shape1, shape2)
    sum((quantile_values - predicted_quantiles)^2)
  }
  
  # Optimize
  result <- optim(c(0, 0), objective, method = "BFGS")
  
  list(
    shape1 = exp(result$par[1]),
    shape2 = exp(result$par[2]),
    convergence = result$convergence,
    value = result$value  # goodness of fit
  )
}

# Iterate over rows (observations)
quantile_probs <- c(0.025, 0.5, 0.975)

beta_fits <- purrr::map_dfr(1:nrow(pred_rf1$predictions), function(i) {
  
  q_vals <- pred_rf1$predictions[i, ]
  
  # Fit Beta distribution
  fit <- fit_beta_to_quantiles(q_vals, quantile_probs)
  
  tibble(
    obs = i,
    shape1 = fit$shape1,
    shape2 = fit$shape2,
    convergence = fit$convergence,
    fit_error = fit$value,
    # Original quantiles for reference
    q025 = q_vals[1],
    q50 = q_vals[2],
    q975 = q_vals[3]
  )
})



## generate list of vectors to calculate depth in real space relative to max
## bathymetry
n_draws <- 10000

# Using purrr (returns a list of vectors)
proportion_draws <- purrr::map(1:nrow(beta_fits), function(i) {
  rbeta(n_draws, shape1 = beta_fits$shape1[i], shape2 = beta_fits$shape2[i])
})

# multiply by max bathy
bathy_vec <- pred_dat1$max_bathy
name_vec <- pred_dat1$subarea

depth_dat <- purrr::pmap(
  list(proportion_draws, bathy_vec, name_vec),
  function (x, y, z) {
    depth_dist <- x * y
    data.frame(
      subarea = z,
      ppn_depth = x,
      real_depth = depth_dist,
      max_bathy = y
    )
  }
) %>% 
  bind_rows()

ggplot(depth_dat) +
  geom_histogram(aes(x = real_depth)) +
  geom_vline(aes(xintercept = max_bathy), color = "red") +
  facet_wrap(~subarea) +
  ggsidekick::theme_sleek()

