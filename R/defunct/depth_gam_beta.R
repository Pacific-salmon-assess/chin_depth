### Relative depth models
# Use detections data to model variability in depth using gamma distribution
# Key hypotheses:
# Coarse scale -- Depth distribution varies among regions/with terminal distance
# and year day
# Fine scale -- Depth distribution varies with time of day/day-night and current
# strength/time relative to slack

## DEFUNCT SINCE GAMMs CANT HANDLE BETA ##

library(tidyverse)
library(mgcv)
library(gratia)


depth_dat <- readRDS(here::here("data", "generated_data", "depth_dat_30min.RDS"))



# quick plots
ggplot(depth_dat) +
  geom_point(aes(x = max_bathy_c, y = rel_depth))
ggplot(depth_dat) +
  geom_point(aes(x = sqrt(pos_max_bathy), y = rel_depth))
ggplot(depth_dat) +
  geom_point(aes(x = day_c, y = rel_depth), alpha = 0.4) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_point(aes(x = hour_c, y = rel_depth), alpha = 0.5) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_boxplot(aes(x = day_night, y = rel_depth)) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_boxplot(aes(x = day_night, y = pos_depth)) +
  facet_grid(stage~region_f)


## PRELIMINARY MODELS ----------------------------------------------------------


## These models are used to structure subsequent analyses

# subset for testing 
unique_codes <- depth_dat %>% 
  filter(stage == "mature") %>% 
  pull(vemco_code) %>% 
  unique() %>% 
  sample(30)

trim_depth <- depth_dat %>% 
  filter(vemco_code %in% unique_codes) %>% 
  droplevels()


# control options
ctrl <- list(niterEM = 0, msVerbose = FALSE, opt = 'optim')#optimMethod="L-BFGS-B")

mod1b <- gam(rel_depth ~ region_f + 
               s(det_day, bs = "tp", m = 2) +
               s(det_day, by = region_f, bs = "tp", m = 1) +
               s(hour, bs = "cc", m = 2) +
               s(hour, by = region_f, bs = "cc", m = 1) +
               s(vemco_code, bs = "re"),
             data = mat_depth,
             family = betar(link = "logit"),
             method = "REML")


## Is it necessary to treat time of day as a smooth or are categorical effects 
# adequate?
mod_hour <- gamm(
  rel_depth ~ region_f + s(hour_c, by = region_f, bs = "cc") +
    s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = betar(link = "logit"),
  method = "REML",
  control = ctrl
)
mod_dn <- gamm(
  rel_depth ~ region_f + day_night_region + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = betar(link = "logit"),
  method = "REML",
  control = ctrl
)

AIC(mod_hour$lme, mod_dn$lme)