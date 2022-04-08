### Depth models
## Trimmed version of depth_gam_gamma.R to share with S. Anderson


library(tidyverse)
library(mgcv)
library(gratia)
library(DHARMa)
library(mgcViz) #necessary for GAM/DHARMa compatibility


# Most columns in DF are self-explanatory, but depth is tag depth and bathy 
# (bathymetry) represents the mean, max, or SD within 800 meters of the moored
# receiver (approximate detection radius)
# Lat/longs refer to location ofreceiver mooring where detection occurred
# Maturation stage is inferred based on location/timing of final detection or 
# body size and is relevant because immature fish are at large longer and aren't
# undergoing a directed migration
# Depth estimates are means within a 60 minute timebin; this is partially to 
# make the number of datapoints tractable and to help deal with autocorrelation
# which is significant even with this pooling; bins were assigned relative to
# the first detection in the first year
# Region is a fairly arbitrary assignment, hence why I'd like to ultimately move
# to a spatially-explicit model


depth_dat <- readRDS(
  here::here("data", "generated_data", "depth_dat_60min.RDS")) 



# quick plots related to hypotheses
# 1) depth varies with bottom bathymetry
# 2) depth distributions vary seasonally particularly in immature fish
# 2) depth distributions vary spatially 
# 3) diurnal patterns to depth distribution but only in certain regions (perhaps
# related to foraging vs. traveling habitats)
ggplot(depth_dat) +
  geom_point(aes(x = max_bathy_c, y = pos_depth))
ggplot(depth_dat) +
  geom_point(aes(x = day_c, y = pos_depth), alpha = 0.4) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_point(aes(x = hour_c, y = pos_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_boxplot(aes(x = day_night, y = pos_depth)) +
  facet_grid(stage~region_f)


# subset dataset to make convergence times reasonable
# focus on immature fish because stronger signals in data than mature
imm_depth <- depth_dat %>% filter(stage == "immature")

set.seed(345)
unique_codes_imm <- imm_depth %>% 
  pull(vemco_code) %>% 
  unique() %>% 
  sample(20)

imm_depth_trim <- depth_dat %>% 
  filter(vemco_code %in% unique_codes_imm) %>% 
  droplevels()


# temporal structure
ggplot(imm_depth_trim) +
  geom_point(aes(x = vemco_code, y = timestamp_n))



# full model structure (this can take some time to run so open to 
# recommendations on control settings)
mod_i_60 <- gamm(
  pos_depth ~ 0 + region_f + s(hour_c, by = region_f, bs = "cc") + 
    s(day_c, k = 4) + 
    s(max_bathy_c, k = 4) + s(vemco_code, bs = "re"),
  # random= list(vemco_code = ~1),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = imm_depth_trim,
  family = Gamma(link = "log"),
  method = "REML",
  control = list(niterEM = 0, msVerbose = FALSE, opt = 'optim'),# optimMethod="L-BFGS-B"),
  niterPQL = 200
)

# large estimate for phi
summary(mod_i_60$lme)

# qq plot is terrible, definitely need to adjust something around the 
# distribution
appraise(mod_i_60$gam)

#check for autocorrelation in normalized residuals
res60 <- resid(mod_i_60$lme, type = "normalized")
acf(res60, lag.max = 36)
pacf(res60, lag.max = 36)
# autocorrelation is actually reasonable here, but I think that's a combination 
# of the random sample of tags in the subsetted data and using 60 minute bins; 
# normally it's still around ~0.3-0.4


## DHARMa residuals 
# I had trouble using the DHARMa package with gamm and didn't have time to pass
# the results as an unsupported model, so as a quick sanity check I just fit a
# standard gam to the full immature dataset (i.e. no AR term) so there are 
# caveats but still largely consistent with above


mod_m_60b <- gam(
  pos_depth ~ 0 + region_f + s(hour_c, by = region_f, bs = "cc") + 
    s(day_c, k = 4) + 
    s(max_bathy_c, k = 4) + s(vemco_code, bs = "re"),
  data = depth_dat %>% filter(stage == "mature"),
  family = Gamma(link = "log"),
  method = "REML"
)
mod_i_60b <- gam(
  pos_depth ~ 0 + region_f + s(hour_c, by = region_f, bs = "cc") + 
    s(day_c, k = 4) + 
    s(max_bathy_c, k = 4) + s(vemco_code, bs = "re"),
  data = imm_depth,
  family = Gamma(link = "log"),
  method = "REML"
)


dharma_resid <- simulateResiduals(fittedModel = mod_i_60b, plot = F)

plot(dharma_resid)
# QQ plot looks like there may be underdispersion; residual plot looks like 
# model underestimates large depth values

# predictor specific have outliers, but no strong pattern (except for region)
plotResiduals(dharma_resid, form = imm_depth$hour_c)
plotResiduals(dharma_resid, form = imm_depth$day_c)
plotResiduals(dharma_resid, form = imm_depth$max_bathy_c)
plotResiduals(dharma_resid, form = imm_depth$region_f)

# test dispersion (not sure if this is appropriate given RE structure, may need
# to resimulate residuals but not sure how with GAMs)
testDispersion(dharma_resid)


