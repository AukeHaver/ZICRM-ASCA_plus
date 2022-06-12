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
    timepoint = seq(1:design_info_list$n_timepoints)%>%
      rep(each=design_info_list$n_treatments*
            design_info_list$n_dummy_replicates) %>% 
      factor(),
    # Column of treatments
    treatment = seq(1:design_info_list$n_treatments) %>%
      rep(each=design_info_list$n_dummy_replicates) %>%
      rep(design_info_list$n_timepoints) %>%
      factor())
  # Set contrasts
  contrasts(output_df$timepoint) <- contr.sum(design_info_list$n_timepoints)
  contrasts(output_df$treatment) <- contr.sum(design_info_list$n_treatments)
  return(output_df)}

# SIMULATE FIXED PARAMETER MATRIX 
# Simulate a fixed parameter matrix with 2 effects and their interactions
simulate_fixed_param_matrix <- function(
  levels = c(3,4), # The levels for each effects (e.g. 3 timepoints and 4 treatments )
  effects = c("A","B","AB"), # The effects in the output matrix
  replicates = 1){ # The number of replicates per level combination
  n_interactions = prod(levels) # Calculate possible interactions between each level of each effect
  A <- rep(contr.sum(levels[1]),each=levels[2]) # A sum-coded matrix for effect A
  B <- apply(contr.sum(levels[2]), # A sum-coded matrix for effect B
             FUN=function(x)rep(x,levels[1]),MARGIN=2) 
  AB <- matrix(rep(0,levels[1]*(levels[2]-1)*levels[1]*levels[2]), # A zero matrix for interaction effects
               ncol=levels[1]*(levels[2]-1)) 
  # Consider x levels for effect A and y levels for effect B.
  # An interaction matrix could be viewed as a singular matrix of size x*x.
  # Each zero in this matrix is instead a matrix of zeroes of size (y*(y-1)).
  # Each one in this matrix is instead a treatment contrast matrix for y treatments.
  # This replacement is performed in the three lines below.
  for(i in 1:levels[1]){ 
    AB[(levels[2]*(i-1)+1):(levels[2]*i),
       ((levels[2]-1)*(i-1)+1):((levels[2]-1)*i)]<-contr.sum(levels[2])}
  out <- rep(1,prod(levels)) # output matrix, initially formatted as a list with ones for the mean.
  if("A" %in% effects){ # Only add effect A to the final matrix if it should be included
    out<-c(out,A)}
  if("B" %in% effects){ # Only add effect B to the final matrix if it should be included
    out<-c(out,B)}
  if("AB" %in% effects){ # Only add effect AB to the final matrix if it should be included
    out<-c(out,AB)}
  return(matrix(rep(out,each=replicates),nrow=prod(levels)*replicates))}

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
  output$replicates <- output$design %>% 
    select(-timepoint) %>%
    distinct()%>%
    group_by(treatment)%>%
    summarize(replicates=n(),.groups="drop") 
  return(output)
}




# GENERATE LONG RESPONSE DATAFRAME 
# Function for generating long-pivoted response dataframe from matrix list
gen_long_response_df <- function(matrix_list){
  return(
    cbind(
      matrix_list$design,
      matrix_list$response %>%
        as_tibble(.name_repair) %>%
        set_names(seq(1:ncol(.)))) %>%
      pivot_longer(cols=!c("subject","timepoint","treatment"),
                   values_to = "count",
                   names_to = "zOTU") %>%
      mutate(zOTU = factor(zOTU, 
                           levels = seq(1:ncol(matrix_list$response)))))
}
