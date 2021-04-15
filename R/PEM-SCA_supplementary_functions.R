# Return the sum of a squared matrix
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

# Return the scores, loadings and singular values of a singular value decomposition (SVD) for a Principal Component Analysis (PCA) analysis.
perform_pca <- function(input_matrix){
  USV <- svd(input_matrix)
  output <- list()
  output$scores <- USV$u %*% diag(USV$d)
  output$loadings <- USV$v
  output$singular_values = USV$d
  return(output)
}


# Extract unique scores from matrix (rounded to 10 decimals)
unique_loadings<- function(PCA_res){
  return(
    unique(
      round(PCA_res$loadings,10)))
}
