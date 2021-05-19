### PEM-SCA supplementary Functions
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA
### Date 05-05-2021

# LOAD PACKAGES (modified 04-05-2021)
library("lme4") # GLMM Fitting
library("tidyverse") # For 
library("kableExtra") # For tables
library("ggpubr") # Extension on ggplot2
library("vegan") # Orthogonal procrustus rotation
library("glmmTMB") # zero-inflated GLMM fitting
library("DHARMa") # Residual analysis 
library("furrr") # Multi-core processing

# Virids color objects (modified 05-05-2021)
dark_purple <- rgb(72/255,40/255,120/255)
teal <- rgb(38/255,130/255,142/255)


# FUNCTION FOR SIMULATION INFO (modified 04-05-2021)
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
  output$jackknife_folds    <- jackknife_folds
  return(output)
}

# MAKE DESIGN DATAFRAME FOR SIMULATION (modified 03-05-2021)
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

# SIMULATE FIXED PARAMETER MATRIX (modified 03-05-2021)
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

# FUNCTION FOR MATRIX COMPOSITION (modified 03-05-2021)
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
    output$Ms <- vector(mode="numeric")
    for(i in 1:prod(design_info_list$n_variables,
                    design_info_list$n_timepoints,
                    design_info_list$n_treatments)){
      output$Ms <- c(output$Ms,
                     MASS::mvrnorm(design_info_list$n_dummy_replicates,
                                   empirical = TRUE,
                                   Sigma=subject_Sigma,
                                   mu=subject_mu))}
    output$Ms <- matrix(output$Ms,ncol=design_info_list$n_variables)
    output$logYhat <- output$logYhat + output$Ms
  }
  return(output)
}

# RESPONSE MATRIX FUNCTION FOR DIFFERENT DISTRIBUTIONS  (modified 03-05-2021)
# Function for generation of a response matrix with a certain distribution
# current distributions: Poisson, Zero-inflated Poisson, Negative Binomial, Zero-inflated Negative Binomial
# FUNCTION FOR SAMPLING FROM EXPONENTIAL FAMILY
generate_response_matrix <- function(input_matrix,
                                     type = "pois",
                                     disp = NULL,
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
             sapply(function(x){MASS::rnegbin(n=1,
                                        mu=x,
                                        theta=disp)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
  else if(type == "zinbinom"){ # Not currently working
    return(input_matrix %>%
             exp() %>%
             matrix(ncol=ncol(input_matrix)) %>%
             sapply(function(x){VGAM::rzinegbin(n=1,
                                      size = disp,
                                      prob = NULL,
                                      pstr0 = prob_str_0,
                                      munb=x)}) %>% 
             matrix(ncol=ncol(input_matrix)))
  }
}

# GENERATE SIMULATION SUBSET MATRICES (modified 04-05-2021)
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

# GLMM FITTING FUNCTION (modified 03-05-2021)
# Function for fitting (zero-inflated) Poisson or NB GLMMS 
fit_glmm <- function(long_response_df,
                     dist_family = "nbinom2",
                     zeroinfl = TRUE){
  if(zeroinfl){
    return(long_response_df %>%
             group_by(zOTU)%>%
             group_nest() %>%
             mutate(model = map(data,function(x)
               glmmTMB(count ~ timepoint + timepoint:treatment+(1|subject),
                       family = dist_family,
                       data=x,
                       ziformula = ~1,
                       REML = TRUE))))
  } else {
    return(long_response_df %>%
             group_by(zOTU)%>%
             group_nest() %>%
             mutate(model = map(data,function(x)
               glmmTMB(count ~ timepoint + timepoint:treatment+(1|subject),
                       family = dist_family,
                       ziformula = ~1,
                       data=x,
                       REML = TRUE))))
  }
}

# MUTLICORE GLMM FITTING FUNCTION (modified 14-05-2021)
# Function for fitting (zero-inflated) Poisson or NB GLMMS 
# Includes skip for failed models
fit_glmm2 <- function(long_response_df,
                      dist_family = "nbinom2",
                      zeroinfl = TRUE,
                      cores = 1){
  plan(multisession, workers = cores)
  
  if(zeroinfl){
    return(long_response_df %>%
             group_by(zOTU)%>%
             group_nest() %>%
             mutate(model = future_map(data,possibly(function(x)
               glmmTMB(count ~ timepoint + timepoint:treatment+(1|subject),
                       family = dist_family,
                       data=x,
                       ziformula = ~1,
                       REML = TRUE),otherwise=NA))))
  } else {
    return(long_response_df %>%
             group_by(zOTU)%>%
             group_nest() %>%
             mutate(model = future_map(data,possibly(function(x)
               glmmTMB(count ~ timepoint + timepoint:treatment+(1|subject),
                       family = dist_family,
                       ziformula = ~1,
                       data=x,
                       REML = TRUE),otherwise=NA))))
  }
}


# GENERATE LONG RESPONSE DATAFRAME (modified 03-05-2021)
# Function for generating long-pivoted response dataframe from matrix list
gen_long_response_df <- function(matrix_list){
  return(cbind(matrix_list$design,
               matrix_list$response %>%
                 as_tibble(.name_repair) %>%
                 set_names(seq(1:ncol(.)))) %>%
           pivot_longer(cols=!c("subject","timepoint","treatment"),
                        values_to = "count",
                        names_to = "zOTU") %>%
           mutate(zOTU = factor(zOTU, levels = seq(1:ncol(matrix_list$response)))))
}

# EXTRACT GLMM COEFFICIENTS FROM glmmTMB (modified 04-05-2021)
# Currently only working for zero_inflated negative binomial
extr_glmm_coef <- function(model_df,
                           zeroinfl = TRUE,
                           family = "nbinom2"){
  if(family =="nbinom2"){
    return( model_df %>%
              mutate(
                # Conditional Fixed Effects
                fixed_cond = map(model,possibly(function(x) fixef(x)$cond,
                                                otherwise=NA)),
                # Zero-inflated Fixed Effects 
                fixed_zi = map(model,possibly(function(x) fixef(x)$zi,
                               otherwise=NA)),
                # Random effects
                random = map(model,possibly(function(x) ranef(x)$cond$subject[,1],
                             otherwise=NA)),
                # Dispersion Parameter
                disp = map(model,possibly(function(x) summary(x)$sigma,
                                          otherwise = NA)) %>% unlist(),
                # Posterior zero-inflation parameter
                # Summary value requires reverse logit link exp(x)/(1+exp(x))
                zi_prob = map(model,possibly(function(x) exp(x$fit$par[1])/(1+exp(x$fit$par[1])),
                                             otherwise=NA))))}
}

# EXTRACT GLMM COEF MULTICORE FROM glmmTMB (modified 14-05-2021)
# Currently only working for zero_inflated negative binomial
extr_glmm_coef2 <- function(model_df,
                            zeroinfl = TRUE,
                            family = "nbinom2"){
  if(family =="nbinom2"){
    return( model_df %>%
              mutate(
                # Conditional Fixed Effects
                fixed_cond = map(model,
                                        possibly(function(x) fixef(x)$cond,
                                                 otherwise = NA)),
                # Zero-inflated Fixed Effects 
                fixed_zi = map(model,
                                      possibly(function(x) fixef(x)$zi,
                                               otherwise=NA)),
                # Random effects
                random = map(model,
                                    possibly(function(x) ranef(x)$cond$subject[,1],
                                             otherwise=NA)),
                # Dispersion Parameter
                disp = map(model,
                                  possibly(function(x) summary(x)$sigma,
                                           otherwise=NA)) %>% unlist(),
                # Posterior zero-inflation parameter
                # Summary value requires reverse logit link exp(x)/(1+exp(x))
                zi_prob = map(model,
                                     possibly(function(x) exp(x$fit$par[1])/(1+exp(x$fit$par[1])),
                                              otherwise = NA))))}
}

# COMPOSE EFFECT MATRICES FROM COEFFICIENTS (Modified 04-05-2021)
# Currently very-much hard-coded
# Integration of effects in design info list could allow for broader analysis
comp_mat_from_coef <- function(glmm_results_df,  # List with GLMM fit results
                               design_info_list, # List with information about simulation
                               fixef_mod_mat){   # Fixed effect model matrix
  
  # Initialize list for storing results
  output <- list()
  # Create coefficient matrix for fixed effects
  output$fixed_coef_mat <- do.call("cbind",
                                   glmm_results_df$fixed_cond) 
  row.names(output$fixed_coef_mat)<-NULL # Required, otherwise rownames error
  # Recompose effect matrices
  # Mean (might need different term for matrix)
  output$M0 <- matrix(fixef_mod_mat)[,1] %*% 
    t(output$fixed_coef_mat[1,])
  colnames(output$M0)<-design_info_list$variables
  # Time effect
  output$Ma <- fixef_mod_mat[,2:3] %*% output$fixed_coef_mat[2:3,]
  colnames(output$Ma)<-design_info_list$variables
  # Interaction effect
  output$Mab <- fixef_mod_mat[,4:15] %*% output$fixed_coef_mat[4:15,]
  colnames(output$Mab)<-design_info_list$variables
  # Subject effect
  output$Ms <- do.call("cbind",glmm_results_df$random)%>%
    apply(FUN=function(x)rep(x,design_info_list$n_timepoints),MARGIN=2) %>%
    matrix(ncol=design_info_list$n_variables)
  colnames(output$Mab)<-design_info_list$variables
  # Return results
  return(output)
}

# PERFORM PCA (modified 04-05-2021)
# Return scores, loadings and singular values of singular value decomposition (SVD)
# For a Principal Component Analysis (PCA).
perform_pca <- function(input_matrix){
  USV <- svd(input_matrix)
  output <- list()
  output$scores <- USV$u %*% diag(USV$d)
  output$loadings <- USV$v
  output$singular_values = USV$d
  return(output)
}

# PERFORM PCA ON EFFECT MATRICES (modified 04-05-2021)
# performs an apply of pca function on all elements of effect matrix list
pca_effect_matrices <- function(glmm_efmat_list){
  return(lapply(glmm_efmat_list,
                FUN="perform_pca"))}


# CALCULATE EXPLAINED VARIANCE FROM SINGULAR VALUES (modified 05-05-2021)
expl_var_from_svd <- function(singular_values,
                              decimals = 1){
  return(round(singular_values^2/sum(singular_values^2)*100,decimals))
}

# RETURN SUM OF SQUARE MATRIX
ssq <- function(number_matrix){
  return(sum(number_matrix^2))
}

# Turn a matrix into a long tibble
param_matrix_to_long_tibble <- function(matrix){
  return(matrix %>% 
           as_tibble() %>%
           pivot_longer(cols = colnames(matrix),
                        names_to = "variable",
                        values_to = "param_value") %>%
           group_by(variable) %>%
           mutate(parameter = row_number(),
                  variable = gsub("variable_","",variable)%>% factor()) %>%
           ungroup())
}
# Column-center a matrix
col_center <- function(data_matrix){
  return(sweep(x      = data_matrix,
               MARGIN = 2,
               STATS  = colMeans(data_matrix),
               FUN    = "-"))}

# Make vector unit vector
unit_vector <- function(input_vector){
  return(input_vector/ sqrt(ssq(input_vector)))
}




# Extract unique scores from matrix (rounded to 10 decimals)
unique_loadings<- function(PCA_res){
  return(
    unique(
      round(PCA_res$loadings,10)))
}




# Function to generate resampling groups
# Population_size = input list of treatment group sizes
# Folds = input vector of number of resampling folds
generate_resampling_groups<- function(population_size,folds){
  group_size = ceiling(population_size/folds)
  group_labels = vector("numeric")
  for(i in 1:length(group_size)){
    group_labels <- rep(seq(1:folds),group_size[i])[1:population_size[i]] %>% sample()%>%c(group_labels,.)
  }
  return(group_labels)
}


# Function for generation of a sum-coded interaction matrix
generate_interaction_matrix <- function(n_levels_A,n_levels_B){
  sum_interactions <- n_levels_A * n_levels_B
  coded_interactions <- (n_levels_B-1)*n_levels_A
  interaction_matrix <- matrix(data=rep(0,sum_interactions*coded_interactions),
                               ncol = coded_interactions)
  for(i in 1:n_levels_A){
    row_start <- n_levels_B*(i-1)+1
    row_stop  <- n_levels_B*i
    col_start <- (n_levels_B-1)*(i-1)+1
    col_stop  <- (n_levels_B-1)*i
    interaction_matrix[row_start:row_stop,col_start:col_stop] <- contr.sum(n_levels_B)
  }
  return(interaction_matrix)}
  


# Orthogonal procrustes rotation
orthprocr_pca <- function(target,query){
  sol <- query
  sol$loadings <- cds::orthprocr(Z = target$loadings,
                                 X = query$loadings)$XQ
  sol$scores <- sol$scores %*% (cds::orthprocr(Z = target$loadings,
                                             X = query$loadings)$Q %>% solve())
  #sol$scores <- sol$scores * solve(sol$loadings)
  return(sol)
}

# PCA PLOT TIBBLE FUNCTION - edited 14-05-2021
## Function for making long tibbles for PCA plots
pca_df_for_plot <- function(pca_object, # a list containing scores and loadings matrices and singular values
                            design_df, # a dataframe with design specifics
                            variable_list # a list of variables corresponding to the loadings
){
  df_list <- list()
  # SCORES
  df_list$scores <-
    # Make dataframe with design and specified components
    cbind(design_df,
          pca_object$scores) %>%
    # Pivot into long table
    pivot_longer(cols=!colnames(design_df),
                 names_to = "PC",
                 values_to = "score") 
  # Loadings
  df_list$loadings <-
    # Make dataframe
    tibble(zOTU = variable_list) %>%
    cbind(pca_object$loadings) %>%
    # Pivot into long tibble
    pivot_longer(cols=!zOTU,
                 names_to = "PC",
                 values_to = "Loadings")
  # Scree plot
  df_list$screeplot <- 
    pca_object$singular_values %>%
    # Calculate explained variance per component
    expl_var_from_svd() %>%
    # Use values as a column in a new tibble
    tibble(percentage = .)%>%
    # Filter for components explaining more than 1% of variation
    filter(percentage > 0)%>%
    # Add Components numbers
    mutate(PC = row_number())
  return(df_list)
}

# PERFORM GLMM JACKKNIFE (modified 15-05-2021)
# long function with hard-coded options
# could require improvements
fit_jackknife_models <- function(fitted_models, # Fitted models on full dataset
                                 design_info_list){   # design info
  
  # Recreate data from fitted models
  original_data <- do.call("rbind",fitted_models$data) %>%
    mutate(zOTU = rep(fitted_models$zOTU, each = nrow(fitted_models$data[[1]]))) %>%
    arrange(timepoint,treatment,subject)
  contrasts(original_data$timepoint) <- contr.sum(n_distinct(original_data$timepoint))
  contrasts(original_data$treatment) <- contr.sum(n_distinct(original_data$treatment))
  
  # Make a dataframe for assigning group labels to each subject
  subject_group_df <- 
    # subjects have the same ID as before
    tibble(subject = original_data$subject %>% unique(),
           subgroup = generate_resampling_groups(original_data %>%
                                                   select(subject,treatment) %>%
                                                   distinct()%>%
                                                   group_by(treatment)%>%
                                                   summarize(replicates = n()) %>%
                                                   pull("replicates"),
                                                 folds = design_info_list$jackknife_folds))
  # Create list to store models
  jackknife_models <- list()
  jackknife_zOTUs  <- list()
  jackknife_PCA <- list()
  jackknife_res <- list()
  # Add subgroup info to original data
  data_with_subgroups <- original_data %>%
    left_join(subject_group_df,
              by="subject") 
  # Iterate over jackknife folds
  for(i in 1:design_info_list$jackknife_folds
  ){
    print(paste0("starting fold ",i))
    fold_start_time <- Sys.time()
    jackknife_models[[i]]<-
      data_with_subgroups %>% 
      filter(!subgroup==i) %>%
      fit_glmm2(long_response_df = .,
                dist_family = "nbinom2",
                zeroinfl = TRUE,
                cores=8)%>%
      extr_glmm_coef(model_df = .,
                     zeroinfl = TRUE,
                     family = "nbinom2") %>% 
      filter(!is.na(model)) %>%
      select(-c(model,data))
    jackknife_zOTUs[[i]]<- jackknife_models[[i]]$zOTU
    print(paste0("fold ",i," finished in ",format(Sys.time()-fold_start_time, digits=2)))
  }
  # Make a list of zOTUs which were succesfully fitted in each fold
  shared_fitted_zOTUs <- Reduce(intersect,jackknife_zOTUs)
  
  output <- list()
  output$zOTUs <- shared_fitted_zOTUs
  output$models<- jackknife_models
  output$full_model <- fitted_models
  output$subject_groups <- subject_group_df
  return(output)}


process_jackknife_models <- function(jackknife_models_list,
                                     design_info_list){
  shared_fitted_zOTUs <- jackknife_models_list$zOTUs
  jackknife_models <- jackknife_models_list$models
  fitted_models <- jackknife_models_list$full_model 
  subject_group_df <- jackknife_models_list$subject_groups
  jackknife_res <- list()
  
  # Filter the results of the overal fit by these zOTUs
  reference_fit <- fitted_models %>% 
    filter(zOTU %in% shared_fitted_zOTUs)
  # Compose parameter matrices and perform pca on effect matrices
  reference_pca <-comp_mat_from_coef2(reference_fit,
                                      design_info_list$X) %>%
    pca_effect_matrices()
  # Iterate over models, compose parameter matrices, perform pca
  # Finally perform procrustus rotation
  for(i in 1:length(jackknife_models)){
    jackknife_fold_design <- 
      design_info_list$design %>%
      left_join(subject_group_df,
                by="subject") %>%
      filter(subgroup!=i)%>%
      select(-subgroup)
    # Modify fixed model design matrix
    fold_X <- design_info_list$design %>%
      cbind(design_info_list$X) %>%
      left_join(subject_group_df,
                by="subject") %>%
      filter(!subgroup==i)%>%
      select(-subgroup) %>%
      select(-colnames(design_info_list$design))%>%
      as.matrix()
    
    
    # Transform models into pca matrices per effect
    fold_pca_res <- comp_mat_from_coef2(filter(jackknife_models[[i]], 
                                               zOTU %in% shared_fitted_zOTUs),
                                        fold_X) %>%
      pca_effect_matrices()
    
    # Perform Procrustes rotation on time and interaction matrices
    fold_pca_res$Ma <- orthprocr_pca(reference_pca$Ma,
                                     fold_pca_res$Ma)
    fold_pca_res$Mab <- orthprocr_pca(reference_pca$Mab,
                                      fold_pca_res$Mab)
    # Calculate explained variance
    
    fold_expl_a <- tibble(fold = i,
                          PC = seq(1:length(fold_pca_res$Ma$singular_values)),
                          perc_expl = fold_pca_res$Ma$singular_values %>% expl_var_from_svd())
    fold_expl_ab <- tibble(fold = i,
                           PC = seq(1:length(fold_pca_res$Mab$singular_values)),
                           perc_expl = fold_pca_res$Mab$singular_values %>% expl_var_from_svd())
    
    # Turning scores and loadings into long tibble
    fold_Ta <- fold_pca_res$Ma$scores %>%
      as_tibble(.name_repair)%>%
      cbind(jackknife_fold_design) %>%
      pivot_longer(cols=-colnames(jackknife_fold_design),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V")
    fold_Tab <- fold_pca_res$Mab$scores %>%
      as_tibble(.name_repair)%>%
      cbind(jackknife_fold_design) %>%
      pivot_longer(cols=-colnames(jackknife_fold_design),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V")
    fold_Pa <- fold_pca_res$Ma$loadings %>%
      as_tibble(.name_repair)%>%
      mutate(var = shared_fitted_zOTUs,
             fold = i)%>%
      pivot_longer(col=-var,
                   names_to = "PC",
                   values_to = "loading",
                   names_prefix = "V")
    fold_Pab <- fold_pca_res$Mab$loadings %>%
      as_tibble(.name_repair)%>%
      mutate(var = shared_fitted_zOTUs,
             fold = i)%>%
      pivot_longer(col=-var,
                   names_to = "PC",
                   values_to = "loading",
                   names_prefix = "V")
    # Merge long tibble list
    if(i == 1){
      jackknife_res$Ta <- fold_Ta
      jackknife_res$Pa <- fold_Pa
      jackknife_res$Tab <- fold_Tab
      jackknife_res$Pab <- fold_Pab
      jackknife_res$perc_a <- fold_expl_a
      jackknife_res$perc_ab<- fold_expl_ab
    }else{
      jackknife_res$Ta <- rbind(jackknife_res$Ta,fold_Ta)
      jackknife_res$Pa <- rbind(jackknife_res$Pa,fold_Pa)
      jackknife_res$Tab <- rbind(jackknife_res$Tab,fold_Tab)
      jackknife_res$Pab <- rbind(jackknife_res$Pab,fold_Pab)
      jackknife_res$perc_a <- rbind(jackknife_res$perc_a,fold_expl_a)
      jackknife_res$perc_ab<- rbind(jackknife_res$perc_ab,fold_expl_a)
    }
  }
  ####################################
  # Calculate 2.5th and 97.5 quantile and median over all folds
  ## Scores
  ### Time effect
  jackknife_res$Ta<-jackknife_res$Ta %>% 
    mutate(score = round(score,
                         digits=10))%>%
    select(-c(subject,treatment)) %>%
    distinct() %>%
    group_by(timepoint,PC)%>%
    summarize(quantile_0.025 = quantile(score,.025),
              median = quantile(score,.5),
              quantile_0.975 = quantile(score,.975),
              .groups="drop_last")
  ### Interaction effect
  jackknife_res$Tab<-jackknife_res$Tab %>% 
    mutate(score = round(score,
                         digits=10))%>%
    select(-subject) %>%
    distinct() %>%
    group_by(timepoint,treatment,PC)%>%
    summarize(quantile_0.025 = quantile(score,.025),
              median = quantile(score,.5),
              quantile_0.975 = quantile(score,.975),
              .groups="drop_last")
  ## Loadings
  ### Time effect
  jackknife_res$Pa <- 
    jackknife_res$Pa %>%
    group_by(var,PC) %>%
    summarize(quantile_0.025 = quantile(loading,.025),
              median = quantile(loading,.5),
              quantile_0.975 = quantile(loading,.975),
              .groups="drop_last")
  ### Interaction effect
  jackknife_res$Pab <-
    jackknife_res$Pab %>%
    group_by(var,PC) %>%
    summarize(quantile_0.025 = quantile(loading,.025),
              median = quantile(loading,.5),
              quantile_0.975 = quantile(loading,.975),
              .groups="drop_last")
  ## Percentage Explained
  ### Time effect
  jackknife_res$perc_a <-
    jackknife_res$perc_a %>%
    group_by(PC) %>%
    summarize(quantile_0.025 = quantile(perc_expl,.025),
              median =         quantile(perc_expl,.5),
              quantile_0.975 = quantile(perc_expl,.975),
              .groups="drop_last")
  ### Time effect
  jackknife_res$perc_ab <-
    jackknife_res$perc_ab %>%
    group_by(PC) %>%
    summarize(quantile_0.025 = quantile(perc_expl,.025),
              median =         quantile(perc_expl,.5),
              quantile_0.975 = quantile(perc_expl,.975),
              .groups="drop_last")
  return(jackknife_res)
}

# COMPOSE EFFECT MATRICES FROM COEFFICIENTS (Modified 15-05-2021)
# Currently very-much hard-coded
# Integration of effects in design info list could allow for broader analysis
comp_mat_from_coef2 <- function(glmm_results_df,  # List with GLMM fit results
                                fixef_mod_mat){   # Fixed effect model matrix
  
  # Initialize list for storing results
  output <- list()
  # Create coefficient matrix for fixed effects
  output$fixed_coef_mat <- do.call("cbind",
                                   glmm_results_df$fixed_cond) 
  row.names(output$fixed_coef_mat)<-NULL # Required, otherwise rownames error
  # Recompose effect matrices
  # Mean (might need different term for matrix)
  output$M0 <- matrix(fixef_mod_mat)[,1] %*% 
    t(output$fixed_coef_mat[1,])
  colnames(output$M0)<-glmm_results_df$zOTU
  # Time effect
  output$Ma <- fixef_mod_mat[,2:3] %*% output$fixed_coef_mat[2:3,]
  colnames(output$Ma)<-glmm_results_df$zOTU
  # Interaction effect
  output$Mab <- fixef_mod_mat[,4:15] %*% output$fixed_coef_mat[4:15,]
  colnames(output$Mab)<-glmm_results_df$zOTU
  # Subject effect
  output$Ms <- do.call("cbind",glmm_results_df$random)%>%
    apply(FUN=function(x)rep(x,3),MARGIN=2) %>%
    matrix(ncol=nrow(glmm_results_df))
  colnames(output$Mab)<-glmm_results_df$zOTU
  # Return results
  return(output)
}