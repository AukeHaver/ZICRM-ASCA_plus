### ZICRM-ASCA+ supplementary Functions
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA
### Date 
# MUTLICORE GLMM FITTING FUNCTION 
# Function for fitting (zero-inflated) Poisson or NB GLMMS 
# Includes skip for failed models
fit_glmm <- function(long_response_df,
                     dist_family,
                     cores = 1){
  if(dist_family %in% c("zipois", "zinbinom")){
    if(dist_family == "zipois"){
      family <- "poisson"
    } else {
      family <- "nbinom2"
    }
    return(
      long_response_df %>%
        group_by(zOTU)%>%
        group_nest() %>%
        mutate(
          model = mclapply(
            X = data,
            FUN = possibly(
              function(x) glmmTMB(
                count ~ timepoint + timepoint:treatment+(1|subject),
                family = family,
                data=x,
                ziformula = ~1,
                REML = TRUE
              ),
              otherwise=NA
            ),
            mc.cores = cores
          )
        )
    )
  } else if(dist_family %in% c("pois", "nbinom")) {
    if(dist_family == "zipois"){
      family <- "poisson"
    } else {
      family <- "nbinom2"
    }
    return(
      long_response_df %>%
        group_by(zOTU)%>%
        group_nest() %>%
        mutate(
          model = mclapply(
            X = data,
            FUN = possibly(
              function(x) glmmTMB(
                count ~ timepoint + timepoint:treatment+(1|subject),
                family = family,
                data=x,
                REML = TRUE
              ),
              otherwise=NA
            ),
            mc.cores = cores
          )
        )
    )
  }
}


# COMPOSE EFFECT MATRICES FROM COEFFICIENTS 
# Integration of effects in design info list could allow for broader analysis
comp_mat_from_coef <- function(glmm_results_df,  # List with GLMM fit results
                               fixef_mod_mat,   # Fixed effect model matrix
                               include_residuals=FALSE){
  
  # Initialize list for storing results
  output <- list()
  # Create coefficient matrix for fixed effects
  output$fixed_coef_mat <- do.call("cbind",
                                   glmm_results_df$fixed_cond) 
  row.names(output$fixed_coef_mat)<-NULL # Required, otherwise rownames error
  # Recompose effect matrices
  # Mean (might need different term for matrix)
  output$M0 <- as.matrix(fixef_mod_mat)[,1] %*% 
    t(as.matrix(output$fixed_coef_mat[1,]))
  colnames(output$M0)<-glmm_results_df$zOTU
  # Time effect
  output$Ma <- fixef_mod_mat[,2:3] %*% output$fixed_coef_mat[2:3,]
  colnames(output$Ma)<-glmm_results_df$zOTU
  # Interaction effect
  output$Mab <- fixef_mod_mat[,4:15] %*% output$fixed_coef_mat[4:15,]
  colnames(output$Mab)<-glmm_results_df$zOTU
  # Combined Time and Interaction effect
  output$Ma_ab <- output$Ma + output$Mab
  # Subject effect
  output$Ms <- do.call("cbind",glmm_results_df$random)%>%
    apply(FUN=function(x)rep(x,3),MARGIN=2) %>%
    matrix(ncol=nrow(glmm_results_df))
  colnames(output$Ms)<-glmm_results_df$zOTU
  # Simulated residuals
  if(include_residuals){
    output$Mres <- do.call("cbind",glmm_results_df$resids)
    colnames(output$Ms)<-glmm_results_df$zOTU}
  # Return results
  return(output)
}

# EXTRACT GLMM COEF FROM glmmTMB 
# Currently only working for zero_inflated negative binomial
extr_glmm_coef <- function(model_df,
                           family,
                           cores = 1){
  output<- model_df %>%
    mutate(
      # Conditional Fixed Effects
      fixed_cond = mclapply(
        X = model,
        FUN = possibly(
          function(x) fixef(x)$cond,
          otherwise = NA),
        mc.cores = cores
      ),
      # Zero-inflated Fixed Effects 
      fixed_zi = mclapply(
        X = model,
        FUN = possibly(
          function(x) fixef(x)$zi,
          otherwise=NA
        ),
        mc.cores = cores
      ),
      # Random effects
      random = mclapply(
        X = model,
        FUN = possibly(
          function(x) ranef(x)$cond$subject[,1],
          otherwise=NA
        ),
        mc.cores = cores
      ),
      # Dispersion Parameter
      sigma = mclapply(
        X = model,
        FUN = possibly(
          function(x) summary(x)$sigma,
          otherwise=NA
        ),
        mc.cores = cores
      ) %>% unlist(),
      # Posterior zero-inflation parameter
      # Summary value requires reverse logit link exp(x)/(1+exp(x))
      zi_prob = mclapply(
        X = model,
        FUN = possibly(
          function(x) exp(x$fit$par[1])/(1+exp(x$fit$par[1])),
          otherwise = NA
        ),
        mc.cores = cores
      )
    )
  return(output)
}


# PERFORM PCA 
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

# PERFORM PCA ON EFFECT MATRICES 
# performs an apply of pca function on all elements of effect matrix list
pca_effect_matrices <- function(glmm_efmat_list){
  return(
    lapply(
      glmm_efmat_list,
      FUN="col_center"
    )%>%
      lapply(
        FUN="perform_pca"
      )
  )
}


# CALCULATE EXPLAINED VARIANCE FROM SINGULAR VALUES 
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


#Orthogonal procrustes rotation with VEGAN (RETURNS POOR ROTATION)
# orthprocr_pca <- function(target,query, n_comp){
#   sol <- query
#   procrustes_rotation <- vegan::procrustes(
#     X = target$loadings[,1:n_comp],
#     Y = query$loadings[,1:n_comp],
#     scale = FALSE
#   )
#   sol$loadings <- procrustes_rotation$Yrot
#   sol$scores <- sol$scores[,1:n_comp] %*% solve(procrustes_rotation$rotation)
#   return(sol)
# }

#Orthogonal procrustes rotation with CDS (WORKS)
orthprocr_pca <- function(target,query, n_comp){
  sol <- query
  sol$loadings <- cds::orthprocr(Z = target$loadings[,1:n_comp],
                                 X = query$loadings[,1:n_comp])$XQ
  sol$scores <- sol$scores %*% (cds::orthprocr(Z = target$loadings[,1:n_comp],
                                               X = query$loadings[,1:n_comp])$Q %>% solve())
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


# PERFORM GLMM JACKKNIFE 
# long function with hard-coded options
# could require improvements
fit_jackknife_models <- function(fitted_models, # Fitted models on full dataset
                                 design_info_list,
                                 cores = 1){   # design info
  
  # Recreate data from fitted models
  original_data <- do.call("rbind",fitted_models$data) %>%
    mutate(zOTU = rep(fitted_models$zOTU, each = nrow(fitted_models$data[[1]]))) %>%
    arrange(timepoint,treatment,subject)
  contrasts(original_data$timepoint) <- contr.sum(n_distinct(original_data$timepoint))
  contrasts(original_data$treatment) <- contr.sum(n_distinct(original_data$treatment))
  
  # Make a dataframe for assigning group labels to each subject
  subject_group_df <- 
    # subjects have the same ID as before
    tibble(
      subject = original_data$subject %>% unique(),
      subgroup = generate_resampling_groups(
        original_data %>%
          select(subject,treatment) %>%
          distinct()%>%
          group_by(treatment)%>%
          summarize(replicates = n()) %>%
          pull("replicates"),
        folds = design_info_list$jackknife_folds
      )
    )
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
      fit_glmm(
        long_response_df = .,
        design_info_list$family,
        cores=cores
      )%>%
      extr_glmm_coef(
        model_df = .,
        design_info_list$family,
        cores = cores) %>% 
      filter(!is.na(model)) %>%
      select(-c(model,data))
    jackknife_zOTUs[[i]]<- jackknife_models[[i]]$zOTU
    print(paste0("fold ",i," finished in ",Sys.time()-fold_start_time))
  }
  # Make a list of zOTUs which were succesfully fitted in each fold
  shared_fitted_zOTUs <- Reduce(intersect,jackknife_zOTUs)
  
  output <- list()
  output$zOTUs <- shared_fitted_zOTUs
  output$models<- jackknife_models
  output$full_model <- fitted_models
  output$subject_groups <- subject_group_df
  return(output)}


# Transform a matrix of scores to a long tibble
score_matrix_to_long_tibble <- function(score_matrix, design, fold){
  score_matrix %>% 
    as_tibble(.name_repair) %>%
    cbind(design) %>%
    pivot_longer(
      cols = -colnames(design),
      names_to = "PC",
      values_to = "score",
      names_prefix = "V"
    ) %>%
    mutate(
      fold = as.character(fold)
    ) %>%
    return()
}
# Transform a matrix of loadings to a long tibble
loading_matrix_to_long_tibble <- function(loading_matrix, varnames, fold){
  loading_matrix %>%
    as_tibble(.name_repair)%>%
    mutate(
      var = varnames,
      fold = fold
    )%>%
    pivot_longer(
      col=-c(var, fold),
      names_to = "PC",
      values_to = "loading",
      names_prefix = "V"
    ) %>%
    mutate(fold = as.character(fold))
}

# Generate tibble with explained variance per PC for a specific jackknife fold
fold_expl_var_tibble_from_singval <- function(singular_values, fold){
  tibble(
    fold = as.character(fold),
    PC = seq(1:length(singular_values)),
    perc_expl = singular_values %>% expl_var_from_svd()
  )
}

# PROCESS GLMM JACKKNIFE 
process_jackknife_models <- function(jackknife_models_list,
                                     design_info_list){
  jackknife_models <- jackknife_models_list$models
  subject_group_df <- jackknife_models_list$subject_groups
  jackknife_res <- list()
  jackknife_res$zOTUs <- jackknife_models_list$zOTUs
  # Filter the results of the overal fit by these zOTUs
  # Compose parameter matrices and perform pca on effect matrices
  jackknife_res$reference_pca  <- jackknife_models_list$full_model %>% 
    filter(zOTU %in% jackknife_res$zOTUs) %>%
    comp_mat_from_coef(design_info_list$X) %>%
    pca_effect_matrices()
  # Iterate over models, compose parameter matrices, perform pca
  # Finally perform procrustus rotation
  # Start of fold loop
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
      select(
        -subgroup,
        -colnames(design_info_list$design)) %>%
      as.matrix()
    # Transform models into pca matrices per effect
    fold_pca_res <- jackknife_models[[i]] %>% 
      filter(zOTU %in% jackknife_res$zOTUs) %>%
      comp_mat_from_coef(fold_X) %>%
      pca_effect_matrices()
    # Calculate the minimum number of components to include
    # If the reference matrix and fold matrix have different number of components. The rotation matrix will not be square and therefore no inverse can be calculated.
    n_components = min(ncol(jackknife_res$reference_pca$Ma$loadings), ncol(fold_pca_res$Ma$loadings))
    # Perform Procrustes rotation on time and interaction matrices
    fold_pca_res$Ma <- orthprocr_pca(jackknife_res$reference_pca$Ma, fold_pca_res$Ma, n_comp = n_components)
    fold_pca_res$Mab <- orthprocr_pca(jackknife_res$reference_pca$Mab, fold_pca_res$Mab, n_comp = n_components)
    fold_pca_res$Ma_ab <- orthprocr_pca(jackknife_res$reference_pca$Ma_ab, fold_pca_res$Ma_ab, n_comp = n_components)
    # Calculate explained variance and turn scores and loadings into long tibble
   for(effect in c("a","ab","a_ab")){
      jackknife_res[[paste0("perc_", effect)]][[i]] <- fold_expl_var_tibble_from_singval(fold_pca_res[[paste0("M", effect)]]$singular_values, i)
      jackknife_res[[paste0("T", effect)]][[i]] <- score_matrix_to_long_tibble(fold_pca_res[[paste0("M", effect)]]$scores, jackknife_fold_design, i)
      jackknife_res[[paste0("P", effect)]][[i]] <- loading_matrix_to_long_tibble(fold_pca_res[[paste0("M", effect)]]$loadings, jackknife_res$zOTUs, i)
    }
  }
  # End of fold loop
  for(effect in c("a","ab","a_ab")){
    jackknife_res[[paste0("perc_", effect)]][["ref"]] <- fold_expl_var_tibble_from_singval(jackknife_res$reference_pca[[paste0("M", effect)]]$singular_values, "ref")
    jackknife_res[[paste0("T", effect)]][["ref"]] <- score_matrix_to_long_tibble(jackknife_res$reference_pca[[paste0("M", effect)]]$scores, design_info_list$design, "ref")
    jackknife_res[[paste0("P", effect)]][["ref"]] <- loading_matrix_to_long_tibble(jackknife_res$reference_pca[[paste0("M", effect)]]$loadings, jackknife_res$zOTUs, "ref")
  }
  for(value in c("Ta", "Tab", "Ta_ab","Pa", "Pab", "Pa_ab","perc_a", "perc_ab", "perc_a_ab")){
    jackknife_res[[value]] <- do.call("rbind", jackknife_res[[value]])
  }
  return(jackknife_res)}



