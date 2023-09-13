## Fit Regression Random Forest w/ Percentile Intervals

## Equivalent to rel_depth_ranger_rf.R except adds population identity
# as covariates

library(tidyverse)
library(caret)
library(recipes)
library(DALEX)
library(DALEXtra)
library(randomForest)


depth_dat_raw1 <- readRDS(
  here::here("data", "depth_dat_nobin.RDS")) %>%
  # approximately 7k detections have no available ROMS data; exclude 
  filter(!is.na(roms_temp))

# remove 2022 tag releases (~6k dets) for training model
depth_dat_raw <- depth_dat_raw1 %>% 
  filter(!grepl("2022", vemco_code))


## FIT -------------------------------------------------------------------------

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
  left_join(., ind_folds, by = "vemco_code") %>% 
  dplyr::select(
    depth = rel_depth, fl, lipid, stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    agg,
    ind_block
  )  

# split by individual blocking
train_depth <- depth_dat %>% filter(!ind_block == "5") %>% droplevels()
test_depth <- depth_dat %>% filter(ind_block == "5") %>% droplevels()
test_depth_22 <- depth_dat_raw1 %>% 
  filter(grepl("2022", vemco_code)) %>% 
  dplyr::select(
    depth = rel_depth, fl, lipid, stage, utm_x, utm_y, day_night,
    det_dayx, det_dayy,
    max_bathy, mean_bathy, mean_slope, shore_dist,
    u, v, w, roms_temp, zoo, oxygen, thermo_depth, moon_illuminated,
    agg)


depth_recipe <- recipe(depth ~ ., 
                       data = train_depth %>% 
                         dplyr::select(-ind_block, -max_bathy)) %>% 
  step_dummy(all_predictors(), -all_numeric())


train_folds <- caret::groupKFold(train_depth$ind_block,
                                 k = length(unique(train_depth$ind_block)))
depth_ctrl <-   caret::trainControl(
  method="repeatedcv",
  index = train_folds
)

# apply recipe to dataframe to make dummy variables and 
train_depth_baked <- prep(depth_recipe) %>%
  bake(., 
       new_data = train_depth %>% 
         dplyr::select(-ind_block, -max_bathy))

# pull model attributes from top ranger
# rf_list <- readRDS(here::here("data", "model_fits", "rf_model_comparison.rds"))
# top_mod <- rf_list[[2]]$top_model

# ranger_rf_gsi <- ranger::ranger(
#   depth ~ .,
#   data = train_depth_baked,
#   #hyperpars based on values from top model which is not saved on all locals
#   num.trees = 2500,
#   mtry = 13,
#   keep.inbag = TRUE,
#   quantreg = TRUE,
#   importance = "permutation",
#   splitrule = "extratrees"
# )
# 
# saveRDS(ranger_rf_gsi,
#         here::here("data", "model_fits", "relative_rf_ranger_gsi.rds"))
ranger_rf_gsi <- readRDS(here::here("data", "model_fits", "relative_rf_ranger_gsi.rds"))


# VARIABLE IMPORTANCE ----------------------------------------------------------

imp_vals <- ranger::importance(ranger_rf_gsi, type = "permutation", scale = F) 

# key for axis labels
var_name_key <- data.frame(
  var = names(imp_vals),
  var_f = c("Fork Length", "Lipid Content", "UTM X", "UTM Y", "Year Day 2", 
            "Year Day 1", "Bottom Depth", "Bottom Slope", "Shore Distance", 
            "Hor. Current 1", "Hor. Current 2", "Vert. Current", "Temperature", 
            "Zooplankton", "Oxygen", "Thermocline Depth", "Lunar Cycle", "Maturity",
            "Day/Night",  "Up. Col.", "ECVI", "Fraser Sub.", "Fraser Year.",
            "Low. Col", "Puget Sound", "WA/OR", "WCVI")
)

imp_dat <- data.frame(
  var = names(imp_vals),
  val = imp_vals
) %>% 
  left_join(., var_name_key, by = "var") %>% 
  arrange(-val) %>% 
  mutate(
    var = fct_reorder(as.factor(var), -val),
    category = case_when(
      var %in% c("mean_bathy", "shore_dist", "utm_x", "utm_y",
                 "mean_slope") ~ "spatial",
      var %in% c("det_day", "det_dayx", "det_dayy",
                 "day_night_night", "moon_illuminated") ~ "temporal",
      var %in% c("stage_mature", "fl", "lipid") ~ "biological",
      grepl("agg", var) ~ "stock",
      TRUE ~ "dynamic"
    )
  ) 


imp_plot <- ggplot(imp_dat, aes(x = fct_reorder(var_f, -val), y = val)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Covariate", y = "Relative Importance") +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png(here::here("figs", "ms_figs_rel", "importance_quantreg_gsi.png"),
    height = 4, width = 6, units = "in", res = 250)
imp_plot
dev.off()

