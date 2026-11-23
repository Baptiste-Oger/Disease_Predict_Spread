
## Function: dpc -----------------------------------------------------------
#' 
#' This function computes a Disease Progress Curve (DPC) using a logistic model.
#'
#' The curve represents disease severity over time, bounded between 0 and 100.
#'
#' @param time Numeric vector. Time points.
#' @param par Numeric vector of length 2. Model parameters:
#'  - par[1]: Growth rate
#'  - par[2]: Inflection point (time of maximum growth rate)
#'
#' @details
#' The model is defined as a logistic function:
#'
#' To avoid numerical overflow, values exceeding a large exponential threshold
#' are capped at 100.
#'
#' @return Numeric vector. Predicted disease severity (0–100).
#'

dpc <- function(time, par){
  
  # Compute logistic disease progress curve
  res <- 100*exp(par[1]*(time - par[2]))/(1 + exp(par[1]*(time - par[2])))
  
  # Cap very large values
  res[exp(par[1]*(time - par[2])) > 10^9] <- 100
  
  return(res)
}




## Function: fit_dpc -------------------------------------------------------
#' 
#' This function fits a Disease Progress Curve (DPC) model to observed data
#' using numerical optimization.
#'
#' @param time Numeric vector. Observed time points.
#' @param y Numeric vector. Observed disease severity values.
#' @param time_pred Optional numeric vector. Time points for prediction.
#'  If NULL, predictions are made at observed time points.
#'
#' @details
#' The fitting procedure:
#'
#' 1. Initializes parameters:
#'  - Growth rate initialized at 0.1
#'  - Inflection point initialized at mean(time)
#'
#' 2. Optimizes parameters based on RMSE as bjective function
#'
#' 3. Generates predictions using the optimized parameters
#'
#' @return A list with:
#'  - optim: Output from "optim" (estimated parameters and diagnostics)
#'  - time: Time points used for prediction
#'  - pred: Predicted disease severity values
#'

fit_dpc <- function(time, y, time_pred = NULL){
  
  # Initial parameter values: growth rate and inflection point 
  par0 <- c(0.1, time[which.min(abs(y - 50))]) 
  
  # Minimize RMSE between observed and predicted values 
  oo <- optim(par0, 
              function(par) rmse(dpc(time, par), y), 
              method = "L-BFGS-B", 
              lower = c(0.001, min(time) - 365), 
              upper = c(0.2, max(time) + 365))
  
  # Use observed time if no prediction time is provided
  if(is.null(time_pred)){
    time_pred <- time
  }
  
  # Predictions 
  pred <- dpc(time_pred, oo$par)
  
  # Output
  return(list(optim = oo,
              time = time_pred,
              pred = pred))
}


## Function: rmse ----------------------------------------------------------
#' 
#' This function computes the Root Mean Squared Error (RMSE) between two vectors.
#'
#' @param y1 Numeric vector. First set of values.
#' @param y2 Numeric vector. Second set of values.
#'
#' @details
#' Computed RMSE, defined as sqrt(mean((y1 - y2)^2))
#'
#' @return Numeric value. The RMSE.
#'

rmse <- function(y1, y2){
  
  return(sqrt(mean((y1 - y2)^2)))
}




## Function: Apply_dpc_by_id ----------------------------------------------
#'
#' This function fits Disease Progress Curves (DPC) for multiple spatial units
#' (identified by IDs) and returns predicted trajectories over a common
#' time sequence.
#'
#'
#' @param df Data frame. Input dataset containing at least:
#'  - id: Identifier for each spatial unit
#'  - date: Time variable
#'  - y: Disease severity variable used to fit the DPC
#'       (scaled between 0 and 1)
#' @param time_seq Numeric vector. Time sequence on which predictions
#'  are computed.
#' @param min_obs_dpc Integer. Minimum number of unique observation dates
#'  required to fit a DPC model.
#'
#' @details
#' The procedure:
#'
#' 1. Loop over each spatial unit (id)
#' 2. Check if sufficient observations are available
#' 3. Extract time and severity values
#' 4. Fit a DPC model using `fit_dpc()`
#' 5. Predict the fitted curve over a common time sequence
#'
#' Disease severity values (0-1) are internally rescaled to percentages
#' (×100) before fitting.
#'
#' @return A data.frame with:
#'  - uuid: Spatial unit identifier
#'  - time: Prediction time points
#'  - pred: Predicted DPC values
#'

Apply_dpc_by_id <- function(df,
                            time_seq,
                            min_obs_dpc) {
  
  # Storage
  df_dpc_list <- list()
  
  ids <- unique(df$id)
  
  # Loop over ids
  for(i in ids){
    
    # Subset
    df_i <- df[df$id == i, ]
    
    # If enough data point
    if(length(unique(df_i$date)) >= min_obs_dpc){
      
      # Time
      x <- as.numeric(df_i$date)
      
      # Severity (0-100)
      y <- df_i$y * 100
      
      # Fit DPC from y and x
      dpc_mod <- fit_dpc(time = x,
                         y = y,
                         time_pred = as.numeric(time_seq))
      
      # Storage
      df_tmp <- data.frame(id = i,
                           time = time_seq,
                           pred = dpc_mod$pred)
      
      df_dpc_list[[length(df_dpc_list) + 1]] <- df_tmp
    }
  }
  
  # Concatenate
  df_dpc <- do.call(rbind, df_dpc_list)
  
  return(df_dpc)
}


## Function: Apply_dpc_by_id_low_obs ---------------------------------------
#'
#' This function fits Disease Progress Curves (DPC) for spatial units (IDs)
#' with insufficient observations by combining local observations with
#' neighborhood-smoothed DPC information.
#'
#' A single DPC model is fitted for each ID using:
#'  - observed values (y) repeated nrepeat times
#'  - neighborhood-smoothed DPC values
#'
#' @param df Data frame containing:
#'  - id: spatial unit identifier
#'  - date: time variable
#'  - y: disease severity (0–1 scale)
#' @param time_seq Numeric vector of prediction time points.
#' @param min_obs_dpc Integer. Minimum number of observation dates required.
#' @param nrepeat Integer. Number of repetitions of observed values
#'  to increase their influence in the fitting (default = length(time_seq)).
#' @param local_smoothed_dpc Matrix/data.frame. Smoothed DPC values
#'  (rows = IDs, columns = time_seq).
#'
#' @details
#' Procedure:
#' 1. Select IDs with low number of observations (< min_obs_dpc)
#' 2. Extract observed values
#' 3. Combine observations + smoothed DPC
#' 4. Fit DPC model
#' 5. Predict over time_seq
#'
#' Values are internally rescaled to percentages (×100).
#'
#' @return data.frame with:
#'  - id: Spatial identifier
#'  - time: Prediction time points
#'  - pred: Predicted DPC based on observations + smoothing
#'

Apply_dpc_by_id_low_obs <- function(df,
                                    time_seq,
                                    min_obs_dpc,
                                    nrepeat = length(time_seq),
                                    local_smoothed_dpc) {
  
  df_dpc_list <- list()
  
  ids <- unique(df$id)
  
  for(i in ids){
    
    df_i <- df[df$id == i, ]
    
    # Only low observation cases
    if(length(unique(df_i$date)) > 0 &&
       length(unique(df_i$date)) < min_obs_dpc){
      
      x <- as.numeric(df_i$date)
      y <- df_i$y * 100
      
      # Retrieve smoothed DPC
      y_smoothed <- local_smoothed_dpc[as.character(i), ]
      
      # Skip if no smoothing available
      if(is.na(y_smoothed[1])){
        
        pred <- rep(NA, length(time_seq))
        
      } else {
        
        dpc_mod <- fit_dpc(
          time = as.numeric(c(rep(x, nrepeat), time_seq)),
          y = c(rep(y, nrepeat), y_smoothed),
          time_pred = as.numeric(time_seq)
        )
        
        pred <- dpc_mod$pred
      }
      
      df_dpc_list[[length(df_dpc_list) + 1]] <- data.frame(
        id = i,
        time = time_seq,
        pred = pred
      )
    }
  }
  
  do.call(rbind, df_dpc_list)
}


## Function: Spatio_temporal_spread ----------------------------------------
#'
#' This function reconstructs and spatially interpolates Disease Progress
#' Curves (DPC) over a prediction grid from sparse field observations or
#' model-derived values.
#'
#' The workflow combines:
#'  - temporal DPC fitting
#'  - neighborhood-based smoothing
#'  - spatial interpolation
#'
#' @param y_col Numeric vector. Disease severity values (0–1 scale).
#' @param id_col Vector. Spatial unit identifiers.
#' @param date_col Date vector. Observation dates.
#' @param geom_col sf geometry column.
#' @param grp_col Vector. Grouping variable used during bandwidth optimization
#' @param time_step_dpc Integer. Temporal resolution used for DPC prediction.
#' @param min_obs_dpc Integer. Minimum number of observations required
#'  to directly fit a DPC.
#' @param distance_interval Numeric vector of length 2.
#'  Search interval for spatial bandwidth optimization.
#' @param grid Data frame or matrix containing prediction coordinates (x y).
#' @param nrepeat Integer. Number of repetitions used to weight observed
#'  values in low-observation DPC fitting.
#'  Default = NULL, replaced internally by length(time_seq).
#'
#' @details
#' Workflow:
#'
#' 1. Build standardized spatial-temporal dataset
#' 2. Fit DPCs for IDs with sufficient observations
#' 3. Spatially smooth fitted DPCs
#' 4. Reconstruct DPCs for low-observation IDs using neighborhood information
#' 5. Refit complete DPC matrix
#' 6. Optimize spatial smoothing bandwidth
#' 7. Interpolate DPCs over the prediction grid
#'
#' Disease severity values are internally converted from 0–1 to percentages.
#'
#' @return A list containing:
#'  - res: Spatially interpolated DPC matrix on the prediction grid
#'  - time_seq: Prediction time sequence
#'  - bw_dpc: Optimal bandwidth used for neighborhood DPC smoothing
#'  - bw_grid: Optimal bandwidth used for final grid interpolation
#'  - mat_dpc_complete: Final DPC matrix
#'

Spatio_temporal_spread <- function(y_col,
                                   id_col,
                                   date_col,
                                   geom_col,
                                   grp_col,
                                   time_step_dpc,
                                   min_obs_dpc,
                                   distance_interval,
                                   grid,
                                   nrepeat = NULL) {
  
  ## Check inputs
  
  input_lengths <- c(length(id_col),
                     length(date_col),
                     length(y_col),
                     length(geom_col),
                     length(grp_col))
  
  if(length(unique(input_lengths)) != 1){
    stop("All input columns must have the same length.")
  }
  
  ## Prepare inputs
  # Standardize input data 
  df_dpc_input <- st_as_sf(data.frame(id = id_col,
                                      date = date_col,
                                      y = y_col,
                                      geometry = geom_col,
                                      group = grp_col))
  
  
  # Coordinates by id
  df_coords <- df_dpc_input %>%
    group_by(id) %>%
    summarise(geometry = st_centroid(st_union(geometry)),
      group = first(group))
  
  
  # Define time sequence for prediction
  time_range <- range(df_dpc_input$date)
  time_seq <- seq(from = time_range[1],
                  to = time_range[2],
                  by = time_step_dpc)
  
  if(is.null(nrepeat)){
    nrepeat <- length(time_seq)
  }
  
  
  ## DPC fitting - sufficient observations
  
  df_dpc <- Apply_dpc_by_id(df = df_dpc_input,
                            time_seq = time_seq,
                            min_obs_dpc = min_obs_dpc)
  
  
  # convert long format to matrix 
  mat_dpc_df <- reshape(df_dpc[, c("id", "time", "pred")],
                        idvar = "id",
                        timevar = "time",
                        direction = "wide")
  
  rownames(mat_dpc_df) <- mat_dpc_df$id
  mat_dpc <- mat_dpc_df[, -1]
  
  ## DPC fitting - low observation ids
  
  ids <- rownames(mat_dpc)
  
  df_coords_sub <- df_coords[df_coords$id %in% ids, ]
  df_coords_sub <- df_coords_sub[match(ids, df_coords_sub$id),]
  
  # Bandwidth optimization (DPC smoothing)
  bw_dpc <- Optimize_bandwidth_CV(IDs = df_coords_sub$id,
                                  values = mat_dpc,
                                  coords = st_coordinates(df_coords_sub$geometry),
                                  groups = df_coords_sub$group,
                                  interval = distance_interval)
  
  # Neighborhood smoothing 
  local_smoothed_dpc <- kernel_smoothing(coords_source = st_coordinates(df_coords_sub$geometry),
                                         values_source = mat_dpc,
                                         coords_new = st_coordinates(df_coords$geometry),
                                         bw = bw_dpc$minimum,
                                         kernel_type = "Gaussian")
  
  rownames(local_smoothed_dpc) <- df_coords$id
  
  # DPC fitting based on smoothing and available values
  df_dpc_low <- Apply_dpc_by_id_low_obs(df = df_dpc_input,
    time_seq = time_seq,
    min_obs_dpc = min_obs_dpc,
    nrepeat = nrepeat,
    local_smoothed_dpc = local_smoothed_dpc)
  
  # Gather all DPC
  df_dpc_complete <- rbind(df_dpc, df_dpc_low)
  
  mat_dpc_df_complete <- reshape(df_dpc_complete[, c("id", "time", "pred")],
                                 idvar = "id",
                                 timevar = "time",
                                 direction = "wide")
  
  rownames(mat_dpc_df_complete) <- mat_dpc_df_complete$id
  mat_dpc_complete <- mat_dpc_df_complete[, -1]
  
  ## Final spatial interpolation 
  
  ids_complete <- rownames(mat_dpc_df_complete)
  df_coords_complete <- df_coords[df_coords$id %in% ids_complete,]
  df_coords_complete <- df_coords_complete[match(ids_complete, df_coords_complete$id),]
  
  # Bandwidth optimization
  bw_grid <- Optimize_bandwidth_CV(IDs = df_coords_complete$id,
                                   values = mat_dpc_complete,
                                   coords = st_coordinates(df_coords_complete$geometry),
                                   groups = df_coords_complete$group,
                                   interval = distance_interval)
  
  # Spatial interpolation to grid 
  res <- kernel_smoothing(coords_source = st_coordinates(df_coords_complete$geometry),
                          values_source = mat_dpc_complete,
                          coords_new = grid,
                          bw = bw_grid$minimum,
                          kernel_type = "Gaussian") / 100
  
  colnames(res) <- as.character(time_seq)
  
  ## Output
  
  return(list(res = res,
              time_seq = time_seq,
              bw_dpc = bw_dpc,
              bw_grid = bw_grid,
              mat_dpc_complete = mat_dpc_complete))
}