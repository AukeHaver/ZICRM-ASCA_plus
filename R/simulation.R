### ZICRM-ASCA+ simulation
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA

# FUNCTION FOR SIMULATION INFO 
# Generates list with design description
# Used in later functions
gen_sim_info <- function(n_treatments,
                         n_timepoints,
                         timepoints,
                         treatments,
                         n_PC,
                         n_variables,
                         n_subjects,
                         n_dummy_replicates,
                         family,
                         prob_zero,
                         sigma,
                         jackknife_folds){
  output <- list()
  output$n_timepoints       <- n_timepoints # Number of levels in time
  output$timepoints         <- timepoints   # Levels of time effect
  output$n_treatments       <- n_treatments # Number of levels in treatment
  output$treatments         <- treatments   # Levels of treatment effect
  output$n_PC               <- n_PC         # Number of Principal components
  output$n_variables        <- n_variables  # Number of variables
  output$variables          <- LETTERS[1:n_variables] # Levels of Variables
  output$n_subjects         <- n_subjects   # Number of subjects
  output$n_dummy_replicates <- n_dummy_replicates # Dummy replicates
  # For unbalanced design
  output$n_observations     <- prod(n_timepoints,
                                    n_treatments,
                                    n_dummy_replicates)
  output$family             <- family
  output$prob_zero          <- prob_zero
  output$sigma               <- sigma
  output$jackknife_folds    <- jackknife_folds
  return(output)
}

# MAKE DESIGN DATAFRAME FOR SIMULATION 
# Generates balanced design for simulation
# Used sum-coding for contrasts
# Columns are ordered timepoint -> treatment 
gen_sim_design <- function(design_info_list){
  # Create output dataframe
  output_df <- tibble(
    # Column of subjects
    subject = 
      seq(design_info_list$n_treatments*
            design_info_list$n_dummy_replicates) %>% 
      rep(design_info_list$n_timepoints) %>% 
      factor(),
    # Column of timepoints
    timepoint = rep(
      x = seq(1:design_info_list$n_timepoints) -1,
      each = design_info_list$n_treatments * design_info_list$n_dummy_replicates
    ),
    # Column of treatments
    treatment = rep(
      x = rep(
        x = seq(1:design_info_list$n_treatments),
        each=design_info_list$n_dummy_replicates
      ),
      times = design_info_list$n_timepoints)
  )
  # Set contrasts
  return(output_df)}


# FUNCTION FOR MATRIX COMPOSITION 
compose_matrices <- function(design_info_list,
                             function_seed = 1,
                             mu = 5,
                             effects = c("time","interaction"),
                             scores_time = NULL,
                             loadings_time = NULL,
                             scores_interaction = NULL,
                             loadings_interaction = NULL,
                             include_subject = FALSE,
                             subject_mu = 0,
                             subject_Sigma = .05
                             
){
  output<-list()
  set.seed(function_seed)
  # Matrix of mean values
  output$M0 <- rep(mu, 
                   prod(design_info_list$n_observations,
                        design_info_list$n_variables)) %>%
    matrix(ncol=design_info_list$n_variables)
  output$logYhat <- output$M0
  # Time effect matrix
  if("time" %in% effects){
    output$Ma <- scores_time %*% t(loadings_time)
    output$logYhat <- output$logYhat + output$Ma}
  # Interaction effect matrix
  if("interaction" %in% effects){
    output$Mab <- scores_interaction %*% t(loadings_interaction)
    output$logYhat <- output$logYhat + output$Mab}
  # Subject Effect
  if(include_subject){
    temp_list <- list()
    output$Ms <- vector(mode="numeric")
    for(i in 1:design_info_list$n_variables){
      output$Ms <- c(output$Ms,
                     MASS::mvrnorm(prod(
                       design_info_list$n_dummy_replicates,
                       design_info_list$n_treatments),
                       empirical = TRUE,
                       Sigma=subject_Sigma,
                       mu=subject_mu))}
    output$Ms <- matrix(output$Ms,ncol=design_info_list$n_variables)
    for(i in 1:design_info_list$n_timepoint){
      temp_list[[i]] <- output$Ms
    }
    output$Ms <- do.call("rbind", temp_list)
    output$logYhat <- output$logYhat + output$Ms
  }
  return(output)
}

# RESPONSE MATRIX FUNCTION FOR DIFFERENT DISTRIBUTIONS 
# Function for generation of a response matrix with a certain distribution
# current distributions: Poisson, Zero-inflated Poisson, Negative Binomial, Zero-inflated Negative Binomial
# FUNCTION FOR SAMPLING FROM EXPONENTIAL FAMILY
generate_response_matrix <- function(input_matrix,
                                     type = "pois",
                                     sigma = NULL,
                                     prob_str_0 = NULL){
  if(type == "pois"){ # Poisson Distribution
    return(input_matrix %>%
             exp() %>%
             matrix(ncol=ncol(input_matrix)) %>%
             sapply(function(x){stats::rpois(n=1,
                                             lambda=x)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
  else if(type == "zipois"){ # Zero-inflated Poisson Distribution
    return(input_matrix %>%
             exp() %>%
             matrix(ncol=ncol(input_matrix)) %>%
             sapply(function(x){VGAM::rzipois(n=1,
                                              lambda=x,
                                              pstr0=prob_str_0)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
  else if(type == "nbinom"){ # Not currently working
    return(input_matrix %>%
             exp() %>%
             matrix(ncol=ncol(input_matrix)) %>%
             sapply(function(x){stats::rnbinom(
               n=1,
               mu=x,
               size=sigma)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
  else if(type == "zinbinom"){ # Appears to be working
    return(input_matrix %>%
             exp() %>%
             matrix(ncol=ncol(input_matrix)) %>%
             sapply(function(x){VGAM::rzinegbin(
               n=1,
               size = sigma,
               pstr0 = prob_str_0,
               munb=x)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
}

# GENERATE SIMULATION SUBSET MATRICES 
# Function For Unbalancing datasets
gen_sim_sub_matrices <- function(design_info_list,
                                 matrix_of_design,
                                 function_seed = 1,
                                 response_matrix,
                                 simulation_matrices,
                                 fixed_model_parameter_matrix){
  # Set Seed
  set.seed(function_seed)
  # Make list to store output
  output <- list()
  # Subsample population size from all dummy replicates, store ID
  subset_subjects <- sample(prod(design_info_list$n_dummy_replicates,
                                 design_info_list$n_treatments),
                            design_info_list$n_subjects,
                            replace=FALSE) %>% 
    sort()
  # Make object to store observation for each ID
  subset_observations <- vector(mode="integer")
  for(i in 1:design_info_list$n_timepoints){
    subset_observations <-
      c(subset_observations, 
        subset_subjects + (i-1)*prod(design_info_list$n_dummy_replicates,
                                     design_info_list$n_treatments))}
  # Unbalance simulation matrices (this could be a seperate function)
  output$response <- response_matrix[subset_observations,]
  output$M0 <- simulation_matrices$M0[subset_observations,]
  output$Ma <- simulation_matrices$Ma[subset_observations,]
  output$Mab<- simulation_matrices$Mab[subset_observations,]
  output$Ms <- simulation_matrices$Ms[subset_observations,]
  output$X  <- fixed_model_parameter_matrix[subset_observations,]
  # Modify Design File
  output$design <- matrix_of_design %>% filter(subject %in% subset_subjects)
  # Calculate replicates per treatment group
  return(output)
}


# GENERATE LONG RESPONSE DATAFRAME 
# Function for generating long-pivoted response dataframe from matrix list
gen_long_response_df <- function(matrix_list){
  colnames(matrix_list$response) <- paste0("zOTU_", seq(1:ncol(matrix_list$response)))
  response_tibble <- as_tibble(
    cbind(
      matrix_list$design %>%
        select(subject, timepoint, treatment),
      matrix_list$response
    )
  ) %>%
    pivot_longer(cols=!c("subject","timepoint","treatment"),
                 values_to = "count",
                 names_to = "zOTU",
                 names_prefix = "zOTU_")
  return(response_tibble)
}
