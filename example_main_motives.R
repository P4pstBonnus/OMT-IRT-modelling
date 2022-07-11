# import helper functions and libraries from utils
source("glmm_model_fitting_utils.R")

# import your transformed data (transformation function can be found in the python script)
empirical_df <- read.csv("example_omt_data.csv")

# optional: removal of comparisions for underrepresented motives
# named list contains underrepresented motives for each picture
unwanted_picture_motive_mapping <- list("6" = list("DAFF", "DPOW"), "8" = list("DAFF"), "10" = list("DAFF"))
empirical_df_reduced <- remove_underrepresented_comparisions(empirical_df, unwanted_picture_motive_mapping)

# add picture effects
empirical_df_reduced$PICAFF <- as.factor(ifelse(empirical_df_reduced$DAFF!=0, empirical_df_reduced$Picture, 1))
empirical_df_reduced$PICACH <- as.factor(ifelse(empirical_df_reduced$DACH!=0, empirical_df_reduced$Picture, 1))
empirical_df_reduced$PICPOW <- as.factor(ifelse(empirical_df_reduced$DPOW!=0, empirical_df_reduced$Picture, 1))

# fit empirical model (without gender, could be added to formula via can_gender:DAFF)
empirical_model <- glmmTMB(cbind(Response1, Response2) ~ 0+(0+DAFF+DACH+DPOW|Person)
                               + PICAFF:DAFF + PICACH:DACH + PICPOW:DPOW
                               + DAFF:SDAFF + DACH:SDACH + DPOW:SDPOW,
                               family=binomial("probit"),data=empirical_df_reduced,
                               control=glmmTMBControl(optimizer=optim,
                                                      optArgs = list(method="BFGS"),
                                                      parallel=8))

# optional (recommended): check model summary for model quality
summary(empirical_model)

# fit n models for bias identification and n models for bias correction
correct_bias_and_calc_reliability_for_model(
  glmm_model = empirical_model, 
  number_of_sim_runs = 10, # choose your n
  log_file_name = "Bias_correction_results.txt", # choose your logfile name
  number_of_motivelevels = 3,
  number_of_picture_effects = 3,
  number_of_sustained_effects = 3,
  first_relevant_aff_picture = 1,
  first_relevant_ach_picture = 1,
  first_relevant_pow_picture = 1,
  number_of_cores = 8,
  unwanted_picture_motive_mapping = unwanted_picture_motive_mapping,
  include_gender = FALSE # only works for can_gender:DAFF effect
)
