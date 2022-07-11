# import helper functions and libraries from utils
source("glmm_model_fitting_utils.R")

# import your transformed data (transformation function can be found in the python script)
empirical_df <- read.csv("example_omt_data_motivelevel.csv")

# optional: removal of comparisions for underrepresented motives
# named list contains underrepresented motives for each picture
# unwanted_picture_motive_mapping <- list("6" = list("DAFF1", "DAFF2", "DAFF3", ... "DPOW5"), "8" = list(...), "10" = list(...))
unwanted_picture_motive_mapping <- list()
empirical_df_reduced <- remove_underrepresented_comparisions(empirical_df, unwanted_picture_motive_mapping)

# add picture effects
# ... for motives
empirical_df_reduced$PICAFF <- as.factor(ifelse(empirical_df_reduced$DAFF!=0, empirical_df_reduced$Picture, 1))
empirical_df_reduced$PICACH <- as.factor(ifelse(empirical_df_reduced$DACH!=0, empirical_df_reduced$Picture, 1))
empirical_df_reduced$PICPOW <- as.factor(ifelse(empirical_df_reduced$DPOW!=0, empirical_df_reduced$Picture, 1))

# ... for motivelevels
empirical_df_reduced$PICAFF <- as.factor(ifelse(empirical_df_reduced$DAFF1!=0,empirical_df_reduced$Picture,
                               ifelse(empirical_df_reduced$DAFF2!=0,empirical_df_reduced$Picture,
                                      ifelse(empirical_df_reduced$DAFF3!=0,empirical_df_reduced$Picture,
                                             ifelse(empirical_df_reduced$DAFF4!=0,empirical_df_reduced$Picture,
                                                    ifelse(empirical_df_reduced$DAFF5!=0,empirical_df_reduced$Picture, 1))))))
empirical_df_reduced$PICACH <- as.factor(ifelse(empirical_df_reduced$DACH1!=0,empirical_df_reduced$Picture,
                               ifelse(empirical_df_reduced$DACH2!=0,empirical_df_reduced$Picture,
                                      ifelse(empirical_df_reduced$DACH3!=0,empirical_df_reduced$Picture,
                                             ifelse(empirical_df_reduced$DACH4!=0,empirical_df_reduced$Picture,
                                                    ifelse(empirical_df_reduced$DACH5!=0,empirical_df_reduced$Picture, 1))))))
empirical_df_reduced$PICPOW <- as.factor(ifelse(empirical_df_reduced$DPOW1!=0,empirical_df_reduced$Picture,
                               ifelse(empirical_df_reduced$DPOW2!=0,empirical_df_reduced$Picture,
                                      ifelse(empirical_df_reduced$DPOW3!=0,empirical_df_reduced$Picture,
                                             ifelse(empirical_df_reduced$DPOW4!=0,empirical_df_reduced$Picture,
                                                    ifelse(empirical_df_reduced$DPOW5!=0,empirical_df_reduced$Picture, 1))))))

# fit empirical model (option 1: with all picture effects on motive levels; option 2: only aggregated picture effects)
# option 1:
empirical_ml_model <- glmmTMB(cbind(Response1,Response2) ~ 0+(0+DAFF1+DAFF2+DAFF3+DAFF4+DAFF5+DACH1+DACH2+DACH3+DACH4+DACH5+DPOW1+DPOW2+DPOW3+DPOW4+DPOW5|Person)
                        + PICAFF:DAFF1 + PICAFF:DAFF2 + PICAFF:DAFF3 + PICAFF:DAFF4 + PICAFF:DAFF5
                        + PICACH:DACH1 + PICACH:DACH2 + PICACH:DACH3 + PICACH:DACH4 + PICACH:DACH5 
                        + PICPOW:DPOW1 + PICPOW:DPOW2 + PICPOW:DPOW3 + PICPOW:DPOW4 + PICPOW:DPOW5
                        + DAFF1:SDAFF1 + DAFF2:SDAFF2 + DAFF3:SDAFF3 + DAFF4:SDAFF4 + DAFF5:SDAFF5 
                        + DACH1:SDACH1 + DACH2:SDACH2 + DACH3:SDACH3 + DACH4:SDACH4 + DACH5:SDACH5 
                        + DPOW1:SDPOW1 + DPOW2:SDPOW2 + DPOW3:SDPOW3 + DPOW4:SDPOW4 + DPOW5:SDPOW5,
                       family=binomial("probit"),data=empirical_df_reduced,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs = list(method="BFGS"),
                                              parallel=8))

# option 2:
empirical_ml_model <- glmmTMB(cbind(Response1,Response2) ~ 0+(0+DAFF1+DAFF2+DAFF3+DAFF4+DAFF5+DACH1+DACH2+DACH3+DACH4+DACH5+DPOW1+DPOW2+DPOW3+DPOW4+DPOW5|Person)
                        + PICAFF:DAFF + PICACH:DACH + PICPOW:DPOW
                        + DAFF1:SDAFF1 + DAFF2:SDAFF2 + DAFF3:SDAFF3 + DAFF4:SDAFF4 + DAFF5:SDAFF5 
                        + DACH1:SDACH1 + DACH2:SDACH2 + DACH3:SDACH3 + DACH4:SDACH4 + DACH5:SDACH5 
                        + DPOW1:SDPOW1 + DPOW2:SDPOW2 + DPOW3:SDPOW3 + DPOW4:SDPOW4 + DPOW5:SDPOW5,
                       family=binomial("probit"),data=empirical_df_reduced,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs = list(method="BFGS"),
                                              parallel=8))

# optional (recommended): check model summary for model quality
summary(empirical_ml_model)

# fit n models for bias identification and n models for bias correction
correct_bias_and_calc_reliability_for_model(
  glmm_model = empirical_ml_model, 
  number_of_sim_runs = 10, # choose your n
  log_file_name = "Bias_correction_results.txt", # choose your logfile name
  number_of_motivelevels = 15,
  number_of_picture_effects = 15, # for first option or else 3 for second option,
  number_of_sustained_effects = 15,
  first_relevant_aff_picture = 1,
  first_relevant_ach_picture = 1,
  first_relevant_pow_picture = 1,
  number_of_cores = 8, # choose your number of cpu cores
  unwanted_picture_motive_mapping = unwanted_picture_motive_mapping
)
