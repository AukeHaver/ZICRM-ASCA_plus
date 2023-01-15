### ZICRM-ASCA+ simulation
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA






# FUNCTION FOR MATRIX COMPOSITION 
compose_matrices <- function(design_info_list,
                             model_matrix,
                             function_seed = 1,
                             mu = 5,
                             scores_time = NULL,
                             loadings_time = NULL,
                             scores_interaction = NULL,
                             loadings_interaction = NULL,
                             subject_mu = 0,
                             subject_Sigma = .05
){
  n_treatments <- length(design_info_list$treatments)
  n_timepoints <- length(design_info_list$timepoints)
  n_observations <- n_treatments * n_timepoints * design_info_list$n_dummy_replicates
  output<-list()
  set.seed(function_seed)
  # Matrix of mean values
  output$M0 <- rep(mu, 
                   prod(n_observations,
                        length(design_info_list$variables))) %>%
    matrix(ncol=length(design_info_list$variables))
  output$logYhat <- output$M0
  # Time effect matrix
  output$Ma <- as.matrix(design$dummy_design_matrix[,5:6]) %*% scores_time %*% t(loadings_time)
  output$logYhat <- output$logYhat + output$Ma
  # Interaction effect matrix
  output$Mab <- as.matrix(design$dummy_design_matrix[,7:14]) %*% scores_interaction %*% t(loadings_interaction)
  output$logYhat <- output$logYhat + output$Mab
  # Subject Effect
  temp_list <- list()
  output$Ms <- vector(mode="numeric")
  for(i in 1:length(design_info_list$variables)){
    output$Ms <- c(output$Ms,
                   MASS::mvrnorm(prod(
                     design_info_list$n_dummy_replicates,
                     n_treatments),
                     empirical = TRUE,
                     Sigma=subject_Sigma,
                     mu=subject_mu))}
  output$Ms <- matrix(output$Ms,ncol=length(design_info_list$variables))
  for(i in 1:n_timepoints){
    temp_list[[i]] <- output$Ms
  }
  output$Ms <- do.call("rbind", temp_list)
  output$logYhat <- output$logYhat + output$Ms
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
                                 function_seed = 1,
                                 response_matrix,
                                 simulation_matrices
                                 ){
  # Set Seed
  set.seed(function_seed)
  # Make list to store output
  output <- list()
  # Subsample population size from all dummy replicates, store ID
  subset_subjects <- sample(prod(design_info_list$n_dummy_replicates,
                                 length(design_info_list$treatments)),
                            design_info_list$n_subjects,
                            replace=FALSE) %>% 
    sort()
  # Make object to store observation for each ID
  subset_observations <- vector(mode="integer")
  for(i in 1:length(design_info_list$timepoints)){
    subset_observations <-
      c(subset_observations, 
        subset_subjects + (i-1)*prod(design_info_list$n_dummy_replicates,
                                     length(design_info_list$treatments)))}
  # Unbalance simulation matrices (this could be a separate function)
  output$response <- response_matrix[subset_observations,]
  output$M0 <- simulation_matrices$M0[subset_observations,]
  output$Ma <- simulation_matrices$Ma[subset_observations,]
  output$Mab<- simulation_matrices$Mab[subset_observations,]
  output$Ms <- simulation_matrices$Ms[subset_observations,]
  output$X  <- design_info_list$X_dummy[subset_observations,]
  # Modify Design File
  output$design <- design_info_list$dummy_design_matrix %>% filter(subject %in% subset_subjects)
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
