### Absolute depth models
# Use detections data to model variability in depth using gamma distribution
# Key hypotheses:
# Coarse scale -- Depth distribution varies among regions/with terminal distance
# and year day
# Fine scale -- Depth distribution varies with time of day/day-night and current
# strength/time relative to slack
# NOTE: Fine scale hypotheses excluded for now based on preliminary analyses 
# in prep_currents.R

library(tidyverse)
library(mgcv)
library(gratia)
library(DHARMa)
library(mgcViz)


# Uses 30 minute bins for timestamps; preliminary models suggest equivalent/
# superior performance to 60 minute bins
depth_dat <- readRDS(here::here("data", "generated_data", "depth_dat_60min.RDS"))


# quick plots
ggplot(depth_dat) +
  geom_point(aes(x = max_bathy_c, y = pos_depth, alpha = 0.4))
ggplot(depth_dat) +
  geom_point(aes(x = log(pos_max_bathy), y = pos_depth, alpha = 0.4))
ggplot(depth_dat) +
  geom_point(aes(x = day_c, y = pos_depth), alpha = 0.3) +
  facet_grid(stage~region_f)
ggplot(depth_dat) +
  geom_point(aes(x = hour_c, y = pos_depth), alpha = 0.3) +
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


## Better to use mean or max bathymetry?
mean_bathy <- gamm(
  pos_depth ~ s(mean_bathy_c) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
max_bathy <- gamm(
  pos_depth ~ s(max_bathy_c) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
AIC(mean_bathy$lme, max_bathy$lme)
# substantially better fit from max_bathy





## Is it necessary to treat time of day as a smooth or are categorical effects 
# adequate?
mod_hour <- gamm(
  pos_depth ~ region_f + #s(hour_c, bs = "cc", m = 2) +
    s(hour_c, by = region_f, bs = "cc") +
    s(max_bathy_c) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
mod_dn <- gamm(
  pos_depth ~ day_night_region + s(max_bathy_c) + 
    s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl,
  niterPQL = 30
)
AIC(mod_hour$lme, mod_dn$lme)

res <- resid(mod_hour$lme, type = "normalized")
acf(res, lag.max = 36)
pacf(res, lag.max = 36)

intervals(mod_hour$lme, which = "var-cov")

appraise(mod_hour$gam)
draw(mod_hour$gam, residuals = TRUE)


## Detection day vs. exit day
mod_day <- gamm(
  pos_depth ~ region_f + s(day_c, bs = "tp") +
    s(max_bathy_c, k = 3) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl,
  niterPQL = 150
)
mod_exit <- gamm(
  pos_depth ~ region_f + s(exit_day_c, bs = "tp", k = 5) +
    s(max_bathy_c, k = 3) + s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = trim_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl,
  niterPQL = 150
)

AIC(mod_day$lme, mod_exit$lme)
## basically identical

res <- resid(mod_day$lme, type = "normalized")
acf(res, lag.max = 36)
pacf(res, lag.max = 36)

appraise(mod_day$gam)
draw(mod_day$gam, residuals = TRUE)



## FULL MODELS -----------------------------------------------------------------

ctrl <- list(niterEM = 0, msVerbose = FALSE, opt = 'optim',
             maxIter = 100)
             # optimMethod="L-BFGS-B")

# mature
mat_depth <- depth_dat %>% 
  filter(stage == "mature") %>%
  droplevels()
imm_depth <- depth_dat %>% 
  filter(stage == "immature") %>%
  droplevels()

# saturated model
mod_m_30 <- gamm(
  pos_depth ~ 0 + region_f + s(hour_c, by = region_f, bs = "cc") + 
    s(day_c) + s(max_bathy_c, k = 4) + s(vemco_code, bs = "re"),
  # random= list(vemco_code = ~1),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = imm_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl,
  niterPQL = 125
)
saveRDS(mod_m_30, 
        here::here("data", "generated_data", "depth_fits", "mature_gamm_30.RDS"))

# check diagnostics
appraise(mod_m_30$gam)


#check for autocorrelation in normalized residuals
res30 <- resid(mod_m_30$lme, type = "normalized")
acf(res30, lag.max = 36)
pacf(res30, lag.max = 36)
# still relatively large amount of autocorrelation

draw(mod_m_30$gam, residuals = TRUE)


# DHARMa residuals
dum <- gam(
  pos_depth ~ 0 + region_f + s(day_c) + s(hour_c, by = region_f, bs = "cc") + 
    s(max_bathy_c, k = 4) + s(vemco_code, bs = "re"),
  # correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = imm_depth,
  family = Gamma(link = "log"),
  method = "REML"
)

sim_resids <- simulateResiduals(fittedModel = dum, plot = F)

plot(sim_resids)
# predictor specific have outliers, but no strong pattern (except for region)
plotResiduals(sim_resids, form = imm_depth$hour_c)
plotResiduals(sim_resids, form = imm_depth$day_c)
plotResiduals(sim_resids, form = imm_depth$max_bathy_c)
plotResiduals(sim_resids, form = imm_depth$region_f)

# test dispersion (not sure if this is appropriate given RE structure, may need
# to resimulate residuals but not sure how with GAMs)
testDispersion(sim_resids)

# test temporal autocorrelation
agg_sim_resids <- recalculateResiduals(sim_resids, 
                                       group = imm_depth$timestamp_n)
testTemporalAutocorrelation(agg_sim_resids, time = imm_depth$timestamp_n)



## immature model
mod_i_30 <- gamm(
  pos_depth ~ 0 + region_f + s(hour_c, by = region_f, bs = "cc") + 
    s(max_bathy_c, k = 4) + s(day_c, k = 3) + 
    s(vemco_code, bs = "re"),
  correlation = corCAR1(form = ~ timestamp_n | vemco_code),
  data = imm_depth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl,
  niterPQL = 100
)
saveRDS(mod_i_30, 
        here::here("data", "generated_data", "depth_fits", "immature_gamm_30.RDS"))
mod_i <- readRDS(
  here::here("data", "generated_data", "depth_fits", "immature_gamm.RDS")
)

# check diagnostics
appraise(mod_i$gam)
appraise(mod_i_30$gam)

#check for autocorrelation in normalized residuals
res <- resid(mod_i$lme, type = "normalized")
acf(res, lag.max = 36)
pacf(res, lag.max = 36)

res30 <- resid(mod_i_30$lme, type = "normalized")
acf(res30, lag.max = 36)
pacf(res30, lag.max = 36)

draw(mod_i$gam, residuals = TRUE)
draw(mod_i_30$gam, residuals = TRUE)


## GENERATE PREDICTIONS --------------------------------------------------------

### Fixed effects

## generate counterfactual predictions
fe_pred_dat_hour <- expand.grid(
  region_f = unique(mat_depth$region_f)[1],
  hour_c = seq(-11.25, 11.75, by = 0.25),
  max_bathy_c = 0,
  day_c = 0,
  vemco_code = unique(mat_depth$vemco_code)[1]
)

# exclude coefficient estimates associated with random effects
exclude_par_locs <- grepl("vemco_code", names(coef(mod1$gam)))
exclude_pars <- names(coef(mod1))[exclude_par_locs]
pred_fe <- predict(mod1$gam, newdata = fe_pred_dat_hour, se.fit = T, 
                   exclude = exclude_pars)

mean_preds <- as.data.frame(pred_fe$fit) 
se_preds <- as.data.frame(pred_fe$se.fit) 
colnames(mean_preds) <- colnames(se_preds) <- c("mu", "power", "scale")

# function to rearrange predictions by parameter
pred_dat2 <- cbind(fe_pred_dat %>% select(-c(region_f, vemco_code)), 
                   mean_preds) %>% 
  pivot_longer(., cols = c(mu:scale), values_to = "mean", 
             names_to = "parameter")
se_preds2 <-  pivot_longer(se_preds, cols = c(mu:scale), values_to = "se", 
                          names_to = "parameter") 
pred_dat2$se <- se_preds2$se

plot_pred <- pred_dat2 %>% 
  mutate(
    upper = ifelse(parameter == "mu", exp(mean + (1.96 * se)), mean + 1.96 * se),
    lower = ifelse(parameter == "mu", exp(mean - (1.96 * se)), mean - 1.96 * se),
    fit = ifelse(parameter == "mu", exp(mean), mean)
  ) 

# visualize predictions
plot_foo <- function(dat, pred) {
  ggplot(dat, aes_string(x = pred)) +
    geom_line(aes(y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
    ggsidekick::theme_sleek()
}

mu_preds <- plot_pred %>% filter(parameter == "mu")

yday_counterfac <- mu_preds %>% 
  filter(pos_mean_bathy == "161",
         hour == "12") %>% 
  plot_foo(dat = ., pred = "det_day")
bathy_counterfac <- mu_preds %>% 
  filter(det_day == "251",
         hour == "12") %>% 
  plot_foo(dat = ., pred = "pos_mean_bathy")
hour_counterfac <- mu_preds %>% 
  filter(pos_mean_bathy == "161",
         det_day == "251") %>% 
  plot_foo(dat = ., pred = "hour")


png(here::here("figs", "depth", "depth_gams", "fe_counterfacs.png"), 
    height = 3, width = 7, res = 300, units = "in")
cowplot::plot_grid(
  yday_counterfac, bathy_counterfac, hour_counterfac,
  nrow = 1)
dev.off()

#counterfacs for regional estimates 
fe_pred_dat2 <- expand.grid(
  region_f = unique(mat_depth$region_f),
  det_day = mean(mat_depth$det_day),
  hour = 12,
  vemco_code = unique(mat_depth$vemco_code)[1]
) %>% 
  left_join(., 
            mat_depth %>% 
              group_by(region_f) %>% 
              summarize(pos_mean_bathy = mean(pos_mean_bathy)),
            by = "region_f")

exclude_par_locs2 <- grepl("vemco_code", names(coef(mod1)))
exclude_pars2 <- names(coef(mod1))[exclude_par_locs2]
pred_fe2 <- predict(mod1, newdata = fe_pred_dat2, se.fit = T, 
                   exclude = exclude_pars2)

mean_preds2 <- as.data.frame(pred_fe2$fit) 
se_preds2 <- as.data.frame(pred_fe2$se.fit) 
colnames(mean_preds2) <- colnames(se_preds2) <- c("mu", "power", "scale")

pred_dat2 <- cbind(fe_pred_dat2, mean_preds2) %>% 
  pivot_longer(., cols = c(mu:scale), values_to = "mean", 
               names_to = "parameter") 
se_preds2 <-  pivot_longer(se_preds2, cols = c(mu:scale), values_to = "se", 
                           names_to = "parameter") 
pred_dat2$se <- se_preds2$se

plot_pred2<- pred_dat2 %>% 
  mutate(
    upper = ifelse(parameter == "mu", exp(mean + (1.96 * se)), mean + 1.96 * se),
    lower = ifelse(parameter == "mu", exp(mean - (1.96 * se)), mean - 1.96 * se),
    fit = ifelse(parameter == "mu", exp(mean), mean)
  ) 


ggplot(plot_pred2 %>% filter(parameter == "mu")) +
  geom_pointrange(aes(x = region_f, y = fit, ymin = lower, ymax = upper)) +
  ggsidekick::theme_sleek()





### Random effects
## generate counterfactual predictions
re_pred_dat <- expand.grid(
  region_f = unique(mat_depth$region_f),
  det_day = seq(min(mat_depth$det_day), max(mat_depth$det_day), 2),
  hour = seq(0, 23, by = 1),
  pos_mean_bathy = seq(1, max(mat_depth$pos_mean_bathy), 5),
  vemco_code = unique(mat_depth$vemco_code)[1]
)

# exclude coefficient estimates associated with random effects
exclude_par_locs <- grepl("vemco_code", names(coef(mod1)))
exclude_pars <- names(coef(mod1))[exclude_par_locs]
pred_re <- predict(mod1, newdata = re_pred_dat, se.fit = T, 
                 exclude = exclude_pars)

mean_preds_re <- as.data.frame(pred_re$fit) 
se_preds_re <- as.data.frame(pred_re$se.fit) 
colnames(mean_preds_re) <- colnames(se_preds_re) <- c("mu", "power", "scale")


# function to rearrange predictions by parameter
pred_dat_re <- cbind(re_pred_dat %>% select(-c(vemco_code)), 
                     mean_preds_re) %>% 
  pivot_longer(., cols = c(mu:scale), values_to = "mean", 
               names_to = "parameter")
se_preds_re <-  pivot_longer(se_preds_re, cols = c(mu:scale), values_to = "se", 
                           names_to = "parameter") 
pred_dat_re$se <- se_preds_re$se

plot_pred_re <- pred_dat_re %>% 
  mutate(
    upper = ifelse(parameter == "mu", exp(mean + (1.96 * se)), mean + 1.96 * se),
    lower = ifelse(parameter == "mu", exp(mean - (1.96 * se)), mean - 1.96 * se),
    fit = ifelse(parameter == "mu", exp(mean), mean)
  ) 


mu_preds_re <- plot_pred_re %>% filter(parameter == "mu")

yday_counterfac_re <- mu_preds_re %>% 
  filter(pos_mean_bathy == "161",
         hour == "12") %>% 
  plot_foo(dat = ., pred = "det_day") +
  facet_wrap(~region_f, ncol = 1)
bathy_counterfac_re <- mu_preds_re %>% 
  filter(det_day == "251",
         hour == "12") %>% 
  plot_foo(dat = ., pred = "pos_mean_bathy") +
  facet_wrap(~region_f, ncol = 1)
hour_counterfac_re <- mu_preds_re %>% 
  filter(pos_mean_bathy == "161",
         det_day == "251") %>% 
  plot_foo(dat = ., pred = "hour") +
  facet_wrap(~region_f, ncol = 1)

png(here::here("figs", "depth", "depth_gams", "re_counterfacs.png"), 
    height = 10, width = 7, res = 300, units = "in")
cowplot::plot_grid(
  yday_counterfac_re, bathy_counterfac_re, hour_counterfac_re,
  nrow = 1)
dev.off()


## CURRENTS MODELS -------------------------------------------------------------

# test currents to evaluate whether reasonable to include
det_rec <- readRDS(here::here("data", "generated_data", "det_currents.RDS"))

mod_ln1 <- gamm(pos_depth ~ 1,
             data = det_rec,
             random = list(vemco_code =  ~ 1),
             family = tw, 
             method = "REML")
mod_gam1 <- gam(pos_depth ~ s(eastward, bs = "tp", k = 4, m = 2) +
               s(vemco_code, bs = "re"),
             data = det_rec, family = tw, method = "REML")
mod1 <- gamm(pos_depth ~ 1,
             data = det_rec,
             random = list(vemco_code =  ~ 1),
             family = tw, 
             method = "REML",
             correlation = corAR1(form = ~ date_time | vemco_code)
             # correlation = corARMA(form = ~ date_time | vemco_code, 
             #                       p = 1, q = 1)
)
acf(resid(mod0$lme, type = "normalized"))
acf(resid(mod1$lme, type = "normalized"))

mod1 <- gam(pos_depth ~ s(eastward, bs = "tp", k = 4, m = 2) +
              s(tag_id, bs = "re"),
            data = det_rec, family = tw, method = "REML")
# or with group level-smoothers
mod2 <- gam(pos_depth ~ s(eastward, bs = "tp", k = 4, m = 2) + 
              s(eastward, tag_id, k = 4, bs = "fs", m = 2),
            data = det_rec, family = tw, method = "REML")
AIC(mod0, mod1, mod2)


#generate predictions
pred_dat <- expand.grid(eastward = seq(-1, 1, length = 100),
                        tag_id = levels(det_rec$tag_id))
pred1 <-  predict(mod1, newdata = pred_dat, se.fit = T, exclude = "s(tag_id)")
pred2 <-  predict(mod2, newdata = pred_dat, se.fit = T, 
                  exclude = "s(eastward,tag_id)")
curr <- tibble(names = c("rand_int", "rand_spline"),
               gams = list(mod1, mod2),
               fe_preds = list(pred1, pred2),
               preds = map(gams, function (x) {
                 predict(x, newdata=pred_dat, se.fit = T) 
               })
)
curr$newdat_fe = map(curr$fe_preds, function (x) {
  pred_dat %>% 
    mutate(
      raw_fit = as.numeric(x$fit),
      se_raw_fit = as.numeric(x$se.fit),
      upper = -1 * exp(raw_fit + (1.96 * se_raw_fit)),
      lower = -1 * exp(raw_fit - (1.96 * se_raw_fit)),
      fit = -1 * exp(raw_fit)
    ) %>% 
    select(-tag_id) %>% 
    distinct()
})
curr$newdat = map(curr$preds, function (x) {
  pred_dat %>% 
    mutate(
      raw_fit = as.numeric(x$fit),
      se_raw_fit = as.numeric(x$se.fit),
      upper = -1 * exp(raw_fit + (1.96 * se_raw_fit)),
      lower = -1 * exp(raw_fit - (1.96 * se_raw_fit)),
      fit = -1 * exp(raw_fit)
    ) 
})

curr %>% 
  select(model = names, newdat_fe) %>% 
  unnest(newdat_fe) %>% 
  ggplot(.) +
  geom_ribbon(aes(x = eastward, ymin = lower, ymax = upper),
              alpha = 0.2) +
  geom_line(aes(x = eastward, y = fit)) + 
  facet_wrap(~model)

curr %>% 
  select(model = names, newdat) %>% 
  unnest(newdat) %>% 
  ggplot(.) +
  # geom_ribbon(aes(x = eastward, ymin = lower, ymax = upper),
  #             alpha = 0.2) +
  geom_line(aes(x = eastward, y = fit, colour = tag_id)) + 
  facet_wrap(~model, scales = "free_y") +
  theme(legend.position = "none")
# unrealistically deep predictions from random spline model


rand_int_pred <- curr %>% 
  select(model = names, newdat) %>% 
  unnest(newdat) %>% 
  filter(model == "rand_int")


ggplot() +
  geom_ribbon(data = rand_int_pred,
              aes(x = eastward, ymin = lower, ymax = upper, fill = tag_id),
              alpha = 0.2) +
  geom_line(data = rand_int_pred,
            aes(x = eastward, y = fit, colour = tag_id)) +
  geom_point(data = det_rec, aes(x = eastward, y = depth, colour = tag_id)) +
  facet_wrap(~tag_id) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none") 
