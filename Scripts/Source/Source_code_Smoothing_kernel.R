
## Function: noyau ---------------------------------------------------------
#' 
#' This function computes kernel values for a given distance vector,
#' according to a specified kernel type.
#'
#' @param u Numeric vector. Standardized distances.
#' @param kernel_type Character. Type of kernel function. Supported values:
#'  - "Gaussian"
#'  - "Linear"
#'  - "Quadratic"
#'  - "Power3"
#'  - "Power4"
#'  - "ExpPower"
#' @param b Numeric. Shape parameter for the "ExpPower" kernel, required only for this kernel.
#'
#' @details
#' The function implements several kernel functions commonly used in spatial smoothing.
#' For compact-support kernels, values outside [0, 1] are set to 0.
#'
#' @return Numeric vector. Kernel weights corresponding to input distances.
#'


noyau <- function(u, kernel_type, b = NULL) { 
  if(kernel_type=="Gaussian"){
    out <- dnorm(u)
  }
  if(kernel_type=="Linear"){
    out <- (1 - u)*  (u >= 0 & u <= 1)  
  }
  if(kernel_type=="Quadratic"){
    out <- (1 - u^2)*  (u >= 0 & u <= 1)  
  }
  if(kernel_type=="Power3"){
    out <- (1 - u^3)*  (u >= 0 & u <= 1)  
  }
  if(kernel_type=="Power4"){
    out <- (1 - u^4)*  (u >= 0 & u <= 1)  
  }
  if(kernel_type == "ExpPower"){
    if(is.null(b)) stop("Parameter b must be provided for ExpPower kernel")
    out <- (b / (2 * pi * gamma(2 / b))) * exp(-abs(u)^b)
  }
  return(out)
}

## Function: rmse ----------------------------------------------------------
#' 
#' This function computes the Root Mean Squared Error (RMSE) between two vectors.
#'
#' @param y1 Numeric vector. First set of values.
#' @param y2 Numeric vector. Second set of values.
#'
#' @details
#' RMSE is defined as sqrt(mean((y1 - y2)^2))
#'
#' @return Numeric value. The RMSE.
#'


rmse <- function(y1,y2){
  
  return(sqrt(mean((y1-y2)^2)))
  
}



## Function: kernel_smoothing ----------------------------------------------
#' 
#' This function performs spatial kernel smoothing to interpolate values
#' from source locations to new locations.
#'
#' @param coords_source Matrix/data.frame. Coordinates of source points (n × 2).
#' @param values_source Matrix/data.frame. Values at source points (n × p).
#' @param coords_new Matrix/data.frame. Coordinates of target points (m × 2).
#' @param bw Numeric. Bandwidth parameter controlling the smoothing scale.
#' @param kernel_type Character. Type of kernel used (see function "noyau").
#' @param dist_seuil Numeric. Maximum distance threshold beyond which predictions
#'  are set to NA. Default is Inf.
#' @param b Numeric. Shape parameter for the "ExpPower" kernel (optional).
#'
#' @details
#' The function:
#'
#' 1. Compute pairwise distances between new and source coordinates
#'     and convert them into weights using the selected kernel
#' 2. Compute weighted averages for each variable, normalized by row-wise sum of weights
#' 3. Apply a distance threshold to filter unreliable predictions
#'
#' If no weights are available (sum = 0), the result is set to NA.
#'
#' @return A matrix of size (m × p) containing smoothed values at new locations.
#'

kernel_smoothing <- function(coords_source , values_source, coords_new,
                             bw , kernel_type, dist_seuil = Inf, b = NULL){
  
  # Compute weights and distances
  res_weights <- get_weights(coords_source, coords_new, bw, kernel_type, b)
  distmin <- apply(res_weights$dist, 1, min)
  weights <- res_weights$weights 
  
  # Initialize output matrix
  values_new <- matrix(0 , nrow = nrow(coords_new), ncol = ncol(values_source))
  
  # Weighted averaging
  for (i in 1:ncol(values_source)) { 
    den <- rowSums(weights)
    den[den == 0] <- NA # Avoid div by 0
    values_new[,i] <- as.numeric(weights %*% values_source[,i])/den
  }
  
  # Apply distance threshold
  values_new[distmin > dist_seuil, 1:ncol(values_new)] <- NA
  
  return(values_new)
}


## Function: get_weights ---------------------------------------------------
#' 
#' This function computes pairwise distances and associated kernel weights
#' between source and target coordinates.
#'
#' @param coords_source Matrix/data.frame. Source coordinates (n × 2).
#' @param coords_new Matrix/data.frame. Target coordinates (m × 2).
#' @param bw Numeric. Bandwidth parameter.
#' @param kernel_type Character. Type of kernel function (see "noyau").
#' @param b Numeric. Shape parameter for the "ExpPower" kernel (optional).
#'
#' @details
#' The function:
#'
#' 1. Computes Euclidean distances between each pair of points
#' 2. Standardizes distances by the bandwidth
#' 3. Applies the kernel function to obtain weights
#'
#' @return A list with:
#'  - dist: Matrix (m × n) of pairwise distances
#'  - weights: Matrix (m × n) of kernel weights
#'



get_weights <- function(coords_source, coords_new, bw, kernel_type, b = NULL){ 
  
  dist <- sqrt(outer(coords_new[,1], coords_source[,1], "-")^2 + 
                 outer(coords_new[,2], coords_source[,2], "-")^2)
  
  return(list(dist = dist,
              weights = noyau( u = dist/bw, 
                               kernel_type = kernel_type,
                               b = b)))
  
} 



## Function: Optimize_bandwidth_CV -----------------------------------------
#' 
#' This function selects the optimal bandwidth parameter for spatial kernel
#' smoothing using cross-validation.
#'
#' The optimization is based on minimizing RMSE computed
#' across predefined groups (e.g., spatial folds or clusters).
#'
#' @param IDs Vector. Identifiers for observations (not directly used in computation,
#'  but included for consistency checks).
#' @param values Matrix/data.frame. Response variables (n × p).
#' @param coords Matrix/data.frame. Coordinates of observations (n × 2).
#' @param groups Vector. Group labels used for cross-validation (length n).
#' @param interval Numeric vector of length 2. Search interval for bandwidth.
#' @param kernel_type Character. Type of kernel function (default: "Gaussian").
#' @param b Numeric. Shape parameter for the "ExpPower" kernel (optional).
#'
#' @details
#' The procedure:
#'
#' 1. Perform group-based cross-validation:
#'    - Each group is used as a test set
#'    - Remaining data form the training set
#'
#' 2. For each candidate bandwidth:
#'    - Apply kernel smoothing from training to test locations
#'    - Compute RMSE for each variable
#'    - Average RMSE across variables and groups
#'
#' 3. Select the bandwidth minimizing the total RMSE
#'
#' @return A list containing:
#'  - minimum: Optimal bandwidth
#'  - objective: Corresponding RMSE value
#'


Optimize_bandwidth_CV <- function(IDs, values, coords, groups, 
                                  interval, 
                                  kernel_type = "Gaussian", 
                                  b = NULL) {
  
  # Quick input check
  if (length(IDs) != nrow(values) | length(IDs) != nrow(coords) | length(IDs) != length(groups)) {
    stop("IDs, values, coords, and groups must have same length")
  }
  
  unique_groups <- unique(groups)
  n_groups <- length(unique_groups)
  
  # Objective function to minimize
  f <- function(bw) {
    
    rmse_total <- 0
    
    for (g in unique_groups) {
      
      # Train/test split
      idx_test  <- which(groups == g)
      idx_train <- which(groups != g)
      
      # Skip case
      if (length(idx_train) == 0 | length(idx_test) == 0) next
      
      # Predict test from smoothing train
      pred <- kernel_smoothing(coords_source = coords[idx_train, , drop = FALSE],
                               values_source = values[idx_train, , drop = FALSE],
                               coords_new    = coords[idx_test, , drop = FALSE],
                               bw = bw,
                               kernel_type = kernel_type,
                               dist_seuil = Inf,
                               b = b)
      
      # RMSE
      rmse_vec <- sapply(1:ncol(values), function(j)
        rmse(pred[, j], values[idx_test, j])
      )
      
      # Aggregate across variables
      rmse_group <- mean(rmse_vec, na.rm = TRUE)
      
      # Sum (works cause same number of elements for every distance tested)
      rmse_total <- rmse_total + rmse_group
    }
    
    return(rmse_total)
  }
  
  # Optimization
  res <- optimize(f, interval)
  
  return(res)
}


## Function: Optimize_bandwidth_ExpPower -----------------------------------
#' 
#' This function jointly selects the optimal bandwidth and shape parameter (b)
#' for the exponential power kernel using cross-validation.
#'
#' It performs a grid search over the shape parameter (b), and for each value,
#' optimizes the bandwidth using "Optimize_bandwidth_CV".
#'
#' @param IDs Vector. Identifiers for observations.
#' @param values Matrix/data.frame. Response variables (n × p).
#' @param coords Matrix/data.frame. Coordinates of observations (n × 2).
#' @param groups Vector. Group labels used for cross-validation.
#' @param interval_bw Numeric vector of length 2. Search interval for bandwidth.
#' @param b_grid Numeric vector. Grid of candidate values for the shape parameter b.
#'
#' @details
#' The procedure:
#'
#' 1. Loop over candidate values of b
#' 2. For each b:
#'    - Optimize bandwidth using cross-validation
#' 3. Select the combination (b, bandwidth) minimizing RMSE
#'
#' @return A data.frame row corresponding to the best parameter combination:
#'  - b: Optimal shape parameter
#'  - bw: Optimal bandwidth
#'  - score: Corresponding RMSE
#'

Optimize_bandwidth_ExpPower <- function(IDs, values, coords, groups,
                                        interval_bw,  b_grid) {
  
  # Grid search over shape parameter b
  res <- lapply(b_grid, function(b) {
    
    opt <- Optimize_bandwidth_CV(IDs = IDs,
                                 values = values,
                                 coords = coords,
                                 groups = groups,
                                 interval = interval_bw,
                                 kernel_type = "ExpPower",
                                 b = b )
    
    list(b = b,
         bw = opt$minimum,
         score = opt$objective)
  })
  
  # Aggregate results
  res_df <- do.call(rbind, lapply(res, as.data.frame))
  
  # Return best combination
  res_df[which.min(res_df$score), ]
}
