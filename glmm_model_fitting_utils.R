library(mvtnorm)
library(glmmTMB)
library(dplyr)

return_probit_vectors <- function(motivelvls) {
  base.names <- motivelvls
  base <- setNames(rep(0, length(base.names)), base.names)
  safe <- list()
  base_combinations <- combn(base.names, 2)
  for (x in 1:(length(base_combinations) / 2)) {
    line <- base
    line[toString(base_combinations[1,x])] <- 1
    line[toString(base_combinations[2,x])] <- -1
    safe <- c(safe, line)
  }
  return(matrix(safe, nrow=length(base_combinations) / 2, 
                ncol=length(base), byrow=TRUE))
}

get_number_of_excluded_picture_motive_effects <- function(unwanted_picture_motive_mapping) {
  result = 0
  for (picture in names(unwanted_picture_motive_mapping)) {
    for (motive in unname(unlist(unwanted_picture_motive_mapping[picture]))) {
      result = result + 1
    }
  }
  return(result)
}

get_empty_spot_number <- function(motive) {
  return(ifelse(motive %in% c("DAFF", "DAFF1"), 0,
             ifelse(motive %in% c("DACH", "DAFF2"), 1,
                    ifelse(motive %in% c("DPOW", "DAFF3"), 2,
                           ifelse(motive == "DAFF4", 3,
                                  ifelse(motive == "DAFF5", 4,
                                         ifelse(motive == "DACH1", 5,
                                                ifelse(motive == "DACH2", 6,
                                                       ifelse(motive == "DACH3", 7,
                                                              ifelse(motive == "DACH4", 8,
                                                                     ifelse(motive == "DACH5", 9,
                                                                            ifelse(motive == "DPOW1", 10,
                                                                                   ifelse(motive == "DPOW2", 11,
                                                                                          ifelse(motive == "DPOW3", 12,
                                                                                                 ifelse(motive == "DPOW45", 13, 14)))))))))))))))
}

add_dummies_to_complete_picture_effect_list <- function(picture_effects, unwanted_picture_motive_mapping, number_of_motive_categories) {
  empty_piceff_spots = c()
  for (picture in names(unwanted_picture_motive_mapping)) {
    for (motive in unname(unlist(unwanted_picture_motive_mapping[picture]))) {
      empty_piceff_spots = c(empty_piceff_spots, 
                             (((as.integer(picture) - 1) * number_of_motive_categories) + 
                              get_empty_spot_number(motive)))
    }
  }
  for (empty_spot_num in empty_piceff_spots) {
    picture_effects <- append(picture_effects, 0, empty_spot_num)
  }
  return(picture_effects)
}

add_base_motives <- function(data, motive) {
  data[motive] <- as.integer(rowSums(
      select(data, matches(paste("^", motive, ".*", sep="")))
  ))
  return(data)
}

generate_omt_data <- function(model, n_motivelevels, n_motivelevels_picture, 
                              n_dynamic_effects, picture_bias=NULL, 
                              sustained_effect_bias=NULL, include_gender=FALSE, unwanted_picture_motive_mapping=NULL) {
  # This function assumes that the order of the fixed effects in the model is
  # picture effects first, followed by dynamic effects
  # n_motivelevels means number of all motive levels except the 0-motive
  
  # Declare some variables
  number_of_omt_pictures <- 15
  list_of_motives <- list("AFF", "ACH", "POW")
  do_lvls_exist <- ifelse(n_motivelevels / length(list_of_motives) == 1, FALSE, TRUE)
  n_subjects <- unname(summary(model)$ngrps$cond)
  number_of_motive_lvl_comparisions <- ((n_motivelevels + 1)*n_motivelevels)/2
  number_of_missing_piceffs <- get_number_of_excluded_picture_motive_effects(
      unwanted_picture_motive_mapping=unwanted_picture_motive_mapping
  )
    
  # extract and reshape picture effects from model AND correct if picture bias is given
  picture_effects=unname(model$fit$par[1:(n_motivelevels_picture * number_of_omt_pictures - number_of_missing_piceffs)])
  if (!is.null(picture_bias)) {
    print("ATTENTION: Including picture effect bias:")
    print(paste(picture_bias[1], "...", picture_bias[length(picture_bias)]))
    picture_effects <- picture_effects + picture_bias
    print("Inclusion finished:")
    print(paste(picture_effects[1], "...", picture_effects[length(picture_effects)]))  
  }
  picture_effects <- add_dummies_to_complete_picture_effect_list(
    picture_effects=picture_effects,
    unwanted_picture_motive_mapping=unwanted_picture_motive_mapping, 
    number_of_motive_categories=n_motivelevels
  )
  pes <- matrix(picture_effects, nrow=number_of_omt_pictures, ncol=n_motivelevels_picture)
    
  dynamic_effects <- unname(model$fit$par[
    (n_motivelevels_picture * number_of_omt_pictures - number_of_missing_piceffs + 1):
      (n_motivelevels_picture * number_of_omt_pictures - number_of_missing_piceffs + n_dynamic_effects)
  ])
  if (!is.null(sustained_effect_bias)) {
    print("ATTENTION: Including sustained effect bias!")
    dynamic_effects <- dynamic_effects + sustained_effect_bias
  }

  # create random gender variable (0 = woman, 1 = man)
  gender <- rnorm(n_subjects, mean=.666477, sd=.471505)
  gender <- ifelse(gender < .5, 0, 1)
  # extract gender effect (only relevant for aff motive)
  if (include_gender) {
    gender_effects <- c(
      rep(unname(model$fit$par[
      (n_motivelevels_picture * number_of_omt_pictures - 
       number_of_missing_piceffs + n_dynamic_effects + include_gender)
      ]), (n_motivelevels/length(list_of_motives))),
      rep(0, (n_motivelevels/length(list_of_motives))), 
      rep(0, (n_motivelevels/length(list_of_motives)))
    )
  } else {
    gender_effects <- c(
      rep(0, (n_motivelevels/length(list_of_motives))),
      rep(0, (n_motivelevels/length(list_of_motives))), 
      rep(0, (n_motivelevels/length(list_of_motives)))
    )
  }
  print(sprintf("GENDER EFFECTS: %s", gender_effects))
  
  possible_motives <- list()
  possible_truemotives <- list()
  for (motive in list_of_motives) {
    for (motivelvl in 1:(n_motivelevels / 3)) {
      possible_motives <- c(possible_motives, 
                            paste("D", motive, ifelse(do_lvls_exist, motivelvl, ""), 
                                  sep=""))
      possible_truemotives <- c(possible_truemotives, paste("true", tolower(motive), ifelse(do_lvls_exist, motivelvl, ""), sep=""))
    }
  }
  possible_motives <- c(possible_motives, "DZERO")
  
  possible_sustainedeffects <- list()
  for (motive in list_of_motives) {
    for (motivelvl in 1:(n_dynamic_effects / 3)) {
      possible_sustainedeffects <- c(possible_sustainedeffects, paste("SD", motive, ifelse(do_lvls_exist, motivelvl, ""), sep=""))
    }
  }
  
  base_combinations <- combn(possible_motives, 2)
  probit_matrix <- return_probit_vectors(possible_motives)
  motive_number_map <- setNames(1:length(possible_motives), possible_motives)
  
  # Extract random effect SD and correlation matrix from model
  ranefsd <- unname(attributes(summary(model)$varcor$cond$Person)$stddev)
  ranefcor <- unname(attributes(summary(model)$varcor$cond$Person)$correlation)
  
  # create random effects for n motive(-levels) for simulation via random numbers
  # drawn from normal distribution for mean 0 and the given random effect covariance matrix sigma
  res <- rmvnorm(n_subjects, mean=rep(0, n_motivelevels), sigma=ranefcor*ranefsd%*%t(ranefsd))
  
  # Create shell/skeleton for simulation
  num_of_cols_till_response <- 1 + n_motivelevels + 1 # Person, raneff, picture
  num_of_needed_cols <- num_of_cols_till_response + 2 + (n_motivelevels + 1) + n_dynamic_effects + 1 # response1/2, probit-vector, dynamic effs, gender
  
  idata <- matrix(c(rep(1:n_subjects, each=number_of_omt_pictures*number_of_motive_lvl_comparisions),
                    rep(as.vector(res), each=number_of_omt_pictures*number_of_motive_lvl_comparisions),
                    rep(rep(1:number_of_omt_pictures, each=number_of_motive_lvl_comparisions), n_subjects),
                    rep(NA, n_subjects*number_of_omt_pictures*number_of_motive_lvl_comparisions*(num_of_needed_cols - num_of_cols_till_response))),
                  ncol=num_of_needed_cols,
                  nrow=n_subjects*number_of_omt_pictures*number_of_motive_lvl_comparisions)
  
  for (i in 1:n_subjects) {
    for (j in 1:number_of_omt_pictures) {
      if (j==1) { tonext <- rep(0,n_motivelevels) } # sustained dynamic effects
      #if (j==1) { tonext2 <- rep(0,n_motivelevels) } # temporary dynamic effects
      uvnr <- c(res[i,] + pes[j,] + tonext * dynamic_effects + gender[i] * gender_effects, 0) # calculate vector of activation probits
      index <- (i*number_of_omt_pictures*number_of_motive_lvl_comparisions-(number_of_omt_pictures*number_of_motive_lvl_comparisions-1))+(j-1)*number_of_motive_lvl_comparisions
      a <- ifelse(is.na(uvnr), -1e5, uvnr)
      # generate a vector of n normally distributed random numbers with mean and sd
      ra <- rnorm(n=n_motivelevels + 1, mean=a, sd=1)
      motivelvl_values_for_picture <- list()
      for (motive in possible_motives) {
        motivelvl_values_for_picture <- c(motivelvl_values_for_picture, setNames(as.integer(ifelse(ra[motive_number_map[motive]]==max(ra), 1, 0)), motive))
      }
      motivelvl_values_for_picture <- c(motivelvl_values_for_picture, "DZERO" = ifelse(ra[n_motivelevels + 1]==max(ra), 1, 0))
      # 
      for (x in 1:(length(base_combinations) / 2)) {
        idata[index + (x - 1), (num_of_cols_till_response+1):num_of_needed_cols] <- unlist(c(
          ifelse(is.na(uvnr[as.integer(motive_number_map[toString(base_combinations[1,x])])]), NA, motivelvl_values_for_picture[toString(base_combinations[1,x])]),
          ifelse(is.na(uvnr[as.integer(motive_number_map[toString(base_combinations[2,x])])]), NA, motivelvl_values_for_picture[toString(base_combinations[2,x])]), 
          probit_matrix[x,], tonext, gender[i]
        ))
      }
      # adapt the sustained dynamic effects by the current result
      for (l in 1:n_motivelevels) {
        tonext[l] <- ifelse(motivelvl_values_for_picture[l] > 0, tonext[l]+1, tonext[l])
      }
    }
  }
  # exclude entries where both responses = 0 and exclude NA's
  idata <- data.frame(idata[idata[,num_of_cols_till_response+1]!=0|idata[,num_of_cols_till_response+2]!=0,])
  idata <- na.exclude(idata)
  names(idata) <- c("Person", possible_truemotives, "Picture", "Response1",
                    "Response2", possible_motives, possible_sustainedeffects, "can_gender")
  # cast categorical vars to int for better glmmTMB performance  
  idata[, c("Person", "Picture", "Response1", "Response2", "can_gender")] <- sapply(idata[, c("Person", "Picture", "Response1", "Response2", "can_gender")], as.integer)
  idata[, unlist(possible_motives)] <- sapply(idata[, unlist(possible_motives)], as.integer)
  idata[, unlist(possible_sustainedeffects)] <- sapply(idata[, unlist(possible_sustainedeffects)], as.integer)
  
  for (main_motive in list_of_motives) {
    main_motive <- paste("D", main_motive, sep="")
    if(!main_motive %in% colnames(idata)) {
      print(paste("Adding", main_motive, "column to dataframe"))
      idata <- add_base_motives(data=idata, motive=main_motive)
    }
  }
  return(idata)
}

remove_underrepresented_comparisions <- function(dataframe, unwanted_picture_motive_mapping) {
  for (picture in names(unwanted_picture_motive_mapping)) {
    for (motive in unname(unlist(unwanted_picture_motive_mapping[picture]))) {
      dataframe <- dataframe[which(!((dataframe[,motive] %in% c(-1, 1)) & (dataframe$Picture == as.integer(picture)))),]
    }
  }
  rownames(dataframe) <- NULL
  return(dataframe)
}

calculate_mean_for_list_of_elements <- function(elements) {
  result = 0
  for (element in elements) {
    result <- result + element
  }
  return(result / length(elements))
}


correct_bias_and_calc_reliability_for_model <- function(glmm_model, number_of_sim_runs, log_file_name, number_of_motivelevels, 
                                                        number_of_picture_effects, number_of_sustained_effects, first_relevant_aff_picture,
                                                        first_relevant_ach_picture, first_relevant_pow_picture, number_of_cores, unwanted_picture_motive_mapping=NULL, 
                                                        include_gender=FALSE) {
  # Write output to log AND screen
  sink(log_file_name, append=TRUE, split=TRUE)
  print(Sys.time())

  number_of_pictures <- 15
  fixed_coefs <- list()
  reliability_scores <- list()
  ranef_std_devs <- list()
  ranef_correlations <- list()
    
  # Step 1: calculate the average over number_of_sim_runs fixed effect estimates
  print(sprintf("ESTIMATING %s MODELS FROM SIMULATED DATASETS", number_of_sim_runs))
  for (i in 1:number_of_sim_runs) {
    print("############")
    print(sprintf("Calculating test model %s ...", i))
    print(Sys.time())
    simulated_testdata_from_model <- generate_omt_data(
      model=glmm_model,
      n_motivelevels=number_of_motivelevels,
      n_motivelevels_picture=number_of_picture_effects,
      n_dynamic_effects=number_of_sustained_effects,
      include_gender=include_gender,
      unwanted_picture_motive_mapping=unwanted_picture_motive_mapping
    )
      
    if (!is.null(unwanted_picture_motive_mapping)) {
      print("ATTENTION: Removing underrepresented comparisions!")
      print(sprintf("Previous length: %s", nrow(simulated_testdata_from_model)))
      simulated_testdata_from_model <- remove_underrepresented_comparisions(
        dataframe=simulated_testdata_from_model,
        unwanted_picture_motive_mapping=unwanted_picture_motive_mapping
      )
      print(sprintf("Shortened length: %s", nrow(simulated_testdata_from_model)))
    }
      
    simulated_testdata_from_model$PICAFF <- as.factor(ifelse(simulated_testdata_from_model$DAFF != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_aff_picture))
    simulated_testdata_from_model$PICACH <- as.factor(ifelse(simulated_testdata_from_model$DACH != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_ach_picture))
    simulated_testdata_from_model$PICPOW <- as.factor(ifelse(simulated_testdata_from_model$DPOW != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_pow_picture))
    print(colnames(simulated_testdata_from_model))
    print(Sys.time())
    model_from_sim_data <- glmmTMB(formula(glmm_model),
                                   family=binomial("probit"), data=select(simulated_testdata_from_model, !matches("^true.*")),
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="BFGS"),
                                                          parallel=number_of_cores))
    print(Sys.time())
    print(sprintf("Model results for test model %s:", i))
    print(summary(model_from_sim_data))
    print(sprintf("Reliability: %s", diag(cor(ranef(model_from_sim_data)$cond$Person, unique(simulated_testdata_from_model[,2:(number_of_motivelevels + 1)]))^2)))
    print(sprintf("R-squared: %s", r2_nakagawa(model_from_sim_data)))
    fixed_coefs <- c(fixed_coefs, list(coef(summary(model_from_sim_data))$cond))
    reliability_scores <- c(reliability_scores, list(diag(cor(ranef(model_from_sim_data)$cond$Person, unique(simulated_testdata_from_model[,2:(number_of_motivelevels + 1)]))^2)))
    ranef_std_devs <- c(ranef_std_devs, list(unname(attributes(summary(model_from_sim_data)$varcor$cond$Person)$stddev)))
    ranef_correlations <- c(ranef_correlations, list(unname(attributes(summary(model_from_sim_data)$varcor$cond$Person)$correlation)))
  }
  
  # Step 2: Calculate average and the difference to the empirical model as bias
  print("CALCULATING THE AVERAGE OF THE RUNS AS WELL AS THE BIAS")
  coef_mean <- calculate_mean_for_list_of_elements(elements = fixed_coefs)
  coef_bias <- coef(summary(glmm_model))$cond - coef_mean
  picture_bias <- coef_bias[1:(number_of_picture_effects * number_of_pictures - 
                               get_number_of_excluded_picture_motive_effects(unwanted_picture_motive_mapping))]
  
  print("SUMMARY OF MODEL AVERAGES:")
  print("Average coefs:")
  print(coef_mean)
  print("Average reliability:")
  print(calculate_mean_for_list_of_elements(elements = reliability_scores))
  print("Average ranef sd:")
  print(calculate_mean_for_list_of_elements(elements = ranef_std_devs))
  print("Average ranef corr:")
  print(calculate_mean_for_list_of_elements(elements = ranef_correlations))
  
  # Step 3: Correct the simulation for the bias
  print(sprintf("ESTIMATING %s MODELS FROM SIMULATED AND BIAS CORRECTED DATASETS", number_of_sim_runs))
  fixed_coefs <- list()
  reliability_scores <- list()
  ranef_std_devs <- list()
  ranef_correlations <- list()
  for (i in 1:number_of_sim_runs) {
    print("############")
    print(sprintf("Calculating test model (WITH BIAS CORRECTION) %s ...", i))
    print(Sys.time())
    simulated_testdata_from_model <- generate_omt_data(
      model = glmm_model,
      n_motivelevels = number_of_motivelevels,
      n_motivelevels_picture = number_of_picture_effects,
      n_dynamic_effects = number_of_sustained_effects,
      picture_bias = picture_bias,
      include_gender=include_gender,
      unwanted_picture_motive_mapping=unwanted_picture_motive_mapping
    )
      
    if (!is.null(unwanted_picture_motive_mapping)) {
      print("ATTENTION: Removing underrepresented comparisions!")
      print(sprintf("Previous length: %s", nrow(simulated_testdata_from_model)))
      simulated_testdata_from_model <- remove_underrepresented_comparisions(
        dataframe=simulated_testdata_from_model,
        unwanted_picture_motive_mapping=unwanted_picture_motive_mapping
      )
      print(sprintf("Shortened length: %s", nrow(simulated_testdata_from_model)))
    }
      
    simulated_testdata_from_model$PICAFF <- as.factor(ifelse(simulated_testdata_from_model$DAFF != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_aff_picture))
    simulated_testdata_from_model$PICACH <- as.factor(ifelse(simulated_testdata_from_model$DACH != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_ach_picture))
    simulated_testdata_from_model$PICPOW <- as.factor(ifelse(simulated_testdata_from_model$DPOW != 0,
                                                             simulated_testdata_from_model$Picture, first_relevant_pow_picture))
    print(Sys.time())
    model_from_sim_data <- glmmTMB(formula(glmm_model),
                                   family=binomial("probit"),data=select(simulated_testdata_from_model, !matches("^true.*")),
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs = list(method="BFGS"),
                                                          parallel=number_of_cores))
    print(Sys.time())
    print(sprintf("Model results for test model %s:", i))
    print(summary(model_from_sim_data))
    print(sprintf("Reliability: %s", diag(cor(ranef(model_from_sim_data)$cond$Person, unique(simulated_testdata_from_model[,2:(number_of_motivelevels + 1)]))^2)))
    print(sprintf("R-squared: %s", r2_nakagawa(model_from_sim_data)))
    fixed_coefs <- c(fixed_coefs, list(coef(summary(model_from_sim_data))$cond))
    reliability_scores <- c(reliability_scores, list(diag(cor(ranef(model_from_sim_data)$cond$Person, unique(simulated_testdata_from_model[,2:(number_of_motivelevels + 1)]))^2)))
    ranef_std_devs <- c(ranef_std_devs, list(unname(attributes(summary(model_from_sim_data)$varcor$cond$Person)$stddev)))
    ranef_correlations <- c(ranef_correlations, list(unname(attributes(summary(model_from_sim_data)$varcor$cond$Person)$correlation)))
  }
  
  # Step 4: Summarize bias correction
  print("SUMMARIZE BIAS CORRECTION:")
  print("Average coefs:")
  print(calculate_mean_for_list_of_elements(elements = fixed_coefs))
  print("Average reliability:")
  print(calculate_mean_for_list_of_elements(elements = reliability_scores))
  print("Average ranef sd:")
  print(calculate_mean_for_list_of_elements(elements = ranef_std_devs))
  print("Average ranef corr:")
  print(calculate_mean_for_list_of_elements(elements = ranef_correlations))
  
  print(Sys.time())
  sink()
}
