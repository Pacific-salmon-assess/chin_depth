# Analysis of Chinook depth distribution data

Data and code associated with analysis of Chinook depth distributions from acoustic telemetry study in southern BC (2019-21). 

Associated with: Chinook salmon depth distributions on the continental shelf are shaped by interactions between location, season, and individual condition (in press at Movement Ecology).

Scripts:
1. prep_depth_roms.R: data cleaning (biological, detections, covariates) and figures related to raw data
2. depth_caret_comparison.R: includes supplementary analysis to identify optimal machine learning model structure (response variable and hyperparameter tuning)
3. rel_depth_ranger_rf: primary analysis fitting random forest with top-ranked hyperparamters to data and associated figures
4. rel_depth_ranger_rf_gsi: as above but including stock identity to generate supplementary figure
5. depth_temperature.R: vertical profiles of oxygen and temperature experience by tagged fish based on ROMS model estimates
6. rel_depth_caret_comparison_weighted.R: fit machine learning models with weighted data for supplementary analysis
7. model_comparison_supp.R: compare hierarchical GAMs and weighted machine learning models to standard machine learning model for supplementary analysis

