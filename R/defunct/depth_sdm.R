## Fit depth SDMs
#Oct 11, 2021


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(raster)
library(sf)


source(here::here("R", "functions", "add_utm.R"))


# import depth data 
depth <- readRDS(here::here("data", "generated_data", "trim_depth_dat.RDS")) %>% 
  add_utm(., lat_name = "latitude", long_name = "longitude") %>% 
  # downscale UTM data
  mutate(pos_depth = depth * -1,
         yUTM_ds = yUTM / 1000,
         xUTM_ds = xUTM / 1000,
         year = as.factor(year),
         sqrt_pos_mean_bathy = sqrt(pos_mean_bathy)) %>% 
  # filter(vemco_code == "7707_2019")
  filter(mat_stage != "immature") %>% 
  slice_sample(., n = 500)



## BUILD SPATIAL PREDICTION GRID -----------------------------------------------

# Select boundary box based on survey area
min_utmy <- min(floor(depth$yUTM))
max_utmy <- max(floor(depth$yUTM))
min_utmx <- min(floor(depth$xUTM))
max_utmx <- max(floor(depth$xUTM))
min_lon <- min(depth$longitude) - 0.1
max_lon <- max(depth$longitude) + 0.1
min_lat <- min(depth$latitude) - 0.1
max_lat <- max(depth$latitude) + 0.1


# First generate grid for entire area that just excludes landmass
projCRS <- "+proj=utm +zone=10 +datum=WGS84"
grid_res <- 2500
coast <- rbind(rnaturalearth::ne_states( "United States of America",
                                         returnclass = "sf"),
               rnaturalearth::ne_states( "Canada", returnclass = "sf"))
coastUTM <- st_transform(coast, crs = projCRS)
cropR <- raster::raster(extent(min_utmx, max_utmx, min_utmy, max_utmy),
                        crs = projCRS, res = grid_res)
g <- fasterize::fasterize(coastUTM, cropR)

## fast conversion pixel to polygons
p <- spex::polygonize(!is.na(g))

## layer is whether we are on land or not
# plot(subset(p, !layer)$geometry)
# points(depth$xUTM, depth$yUTM, col = "red")
# plot(ch, alpha = 0.2, add = TRUE)

# Extract coords and pass to grid
subgrid_coords <- subset(p, !layer)$geometry %>%
  st_sf(ID = seq(1, length(.), by = 1)) %>%
  st_coordinates(.)
pred_grid <- data.frame(X = subgrid_coords[, "X"],
                        Y = subgrid_coords[, "Y"]) %>%
  mutate(grid_id = 1:nrow(.))

pred_sf <- st_as_sf(pred_grid, coords = c("X", "Y"), crs = projCRS) %>%
  st_buffer(., dist = grid_res / 2)

# plot(subset(p, !layer)$geometry)
# plot(pred_sf$geometry, col = alpha("red", 0.2), add = TRUE)
# plot(pred_bathy$geometry, col = alpha("blue", ))

# Import bathymetry dataset
bathy_sf <- readRDS(here::here("data", "bathymetry", "bathy_lowres.RDS")) %>%
  # readRDS(here::here("data", "bathymetry", "big_bathymetry",
  #                    "bathy_hires.RDS")) %>%
  filter(lon > min_lon & lon < max_lon,
         lat > min_lat & lat < max_lat) %>%
  st_as_sf(., coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = projCRS)

# find intersection between prediction grid and bathymetry data
pred_bathy <- st_intersection(pred_sf, bathy_sf) %>%
  group_by(grid_id) %>%
  summarize(mean_depth = mean(depth, na.rm = T),
            sd_depth = sd(depth, na.rm = T)) %>%
  ungroup()

depth_pred_grid <- left_join(pred_grid,
                             pred_bathy %>%
                               as.data.frame() %>%
                               dplyr::select(-geometry),
                             by = "grid_id")
 
saveRDS(depth_pred_grid,
        here::here("data", "bathymetry", "depth_pred_grid_2500m.rds"))



## BUILD MESHES ----------------------------------------------------------------

depth_pred_grid <- readRDS(
  here::here("data", "bathymetry", "depth_pred_grid_2500m.rds"))


## Default spde mesh
# sp <- make_mesh(depth, c("xUTM_ds", "yUTM_ds"), n_knots = 250,
#                 type = "cutoff_search")
# sp2 <- make_mesh(depth, c("xUTM_ds", "yUTM_ds"), n_knots = 250,
#                  type = "kmeans")
# sp3 <- make_mesh(depth, c("xUTM_ds", "yUTM_ds"), cutoff = 4, 
#                  type = "cutoff")
# sp4 <- make_mesh(depth, c("xUTM_ds", "yUTM_ds"), cutoff = 5, 
#                  type = "cutoff")
# plot(sp)
# plot(sp2)
# plot(sp3)
# plot(sp4)


## Alternative mesh based on predictive grid (accounts for landmasses)

#downscaled predictive grid
ds_pred_grid <- depth_pred_grid %>% 
  mutate(X = X / 1000,
         Y = Y / 1000)

# make INLA mesh
bnd <- INLA::inla.nonconvex.hull(as.matrix(ds_pred_grid), convex = -0.03)
plot(bnd$loc)
mesh.loc <- SpatialPoints(as.matrix(cbind(depth$xUTM_ds, 
                                          depth$yUTM_ds)))
mesh <- INLA::inla.mesh.2d(loc=mesh.loc,
                           boundary=list(
                             bnd,
                             NULL),
                           # max.edge=c(2.5, 4.5),
                           min.angle=c(30, 21),
                           max.n=c(48000, 16000), ## Safeguard against large meshes.
                           max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                           cutoff = 5,#1.51, ## Filter away adjacent points.
                           offset=c(5, 10)
                           ) ## Offset for extra boundaries, if needed.
plot(mesh)

sp_inla <- make_mesh(depth, c("xUTM_ds", "yUTM_ds"), mesh = mesh)
plot(sp_inla)



## FIT MODELS ------------------------------------------------------------------

mod1 <- sdmTMB(
  pos_depth ~ sqrt_pos_mean_bathy, #year_f,
  data = depth,
  spde = sp_inla,
  silent = FALSE,
  anisotropy = TRUE,
  # include_spatial = TRUE,
  family = Gamma(link = "log"),
  time = "timestamp_n"
)
mod2 <- sdmTMB(
  pos_depth ~ sqrt_pos_mean_bathy, #year_f,
  data = depth,
  spde = sp_inla,
  silent = FALSE,
  anisotropy = TRUE,
  ar1_fields = TRUE,
  # include_spatial = TRUE,
  family = Gamma(link = "log"),
  time = "timestamp_n"
)

qqnorm(residuals(mod1)); abline(a = 0, b = 1)



## GENERATE PREDS --------------------------------------------------------------

pred_grid_out <- depth_pred_grid %>% 
  mutate(
    xUTM_ds = X / 1000,
    yUTM_ds = Y / 1000,
    mean_depth_z = (mean_depth - mu_obs_depth) / sig_obs_depth,
    sd_depth_z = (sd_depth - mu_sd_obs_depth) / sig_sd_obs_depth,
  )