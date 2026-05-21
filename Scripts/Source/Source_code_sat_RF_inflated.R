library(ranger)
library(caret)
library(pROC)




## Function: fit_RF_zero_inflated ------------------------------------------
#' 
#' This function fits a hurdle-type Random Forest model. It combines:
#'  - A classification model to predict the probability of presence (non-zero values)
#'  - A regression model fitted only on positive observations
#'
#' Hyperparameters of both models are tuned via cross-validation using
#' the function "tune_RF_hyperparams", see bellow.
#'
#' @param training.severity Numeric vector. Response variable, potentially zero-inflated.
#' @param training.covariates Data frame or matrix. Predictor variables.
#' @param mtry_grid Integer vector. Grid of "mtry" values to test.
#' @param ntree_grid Integer vector. Grid of number of trees to test.
#' @param n_folds Integer. Number of folds for cross-validation. Default is 10.
#'
#' @details
#' The modeling procedure follows two steps to handle zero-inflated response variable
#'
#' 1. Classification step
#'  - Converts the response into a binary variable (0 vs >0)
#'  - Tunes and fits a Random Forest classifier using package "ranger"
#'  - Outputs class probabilities
#'  - Also find optimal threshold  using "roc" from the "pROC" package
#'
#' 2. Regression step
#'  - Subsets data to strictly positive observations
#'  - Tunes and fits a Random Forest regression model
#'
#' @return A list with the following elements:
#'  - rf_class: Fitted Random Forest classification model ("ranger" object)
#'  - rf_reg: Fitted Random Forest regression model ("ranger" object)
#'  - threshold: Optimal probability threshold for classification
#'  - tuning: List containing tuning results for classification and regression
#'  - info: Metadata about the model fitting procedure
#'


fit_RF_zero_inflated <- function(training.severity,
                                 training.covariates,
                                 mtry_grid = 1:10,
                                 ntree_grid = seq(500, 1500, 100),
                                 n_folds = 10) {
  
  ## 1. Classification
  
  # Convert response to binary (presence/absence)
  training.severity.bin <- as.factor(ifelse(training.severity > 0, 1, 0))
  df_class <- data.frame(training.severity.bin, training.covariates)
  
  # Hyperparameter tuning for classification
  tune_class <- tune_RF_hyperparams(Eq = training.severity.bin ~ .,
                                    data = df_class,
                                    type = "classification",
                                    mtry_grid = mtry_grid,
                                    ntree_grid = ntree_grid,
                                    n_folds = n_folds)
  
  best_class <- tune_class$best_params
  
  # Fit classification model
  rf_class <- ranger(training.severity.bin ~ .,
                     data = df_class,
                     classification = TRUE,
                     num.trees = best_class$num.trees,
                     mtry = best_class$mtry,
                     probability = TRUE,
                     importance = "impurity")
  
  # Find optimal custom threshold
  pred_class <- predict(rf_class, data = df_class)
  prob_train <- pred_class$predictions[, "1"]
  roc_obj <- roc(response = training.severity.bin, predictor = prob_train, quiet = TRUE)
  threshold_opt <- coords(roc_obj, x = "best", best.method = "youden", ret = "threshold")
  
  ## 2. Regression
  
  # Select positive observations
  select_positive <- training.severity > 0
  df_reg <- data.frame(training.severity = training.severity[select_positive],
                       training.covariates[select_positive, , drop = FALSE])
  
  # Hyperparameter tuning for regression
  tune_reg <- tune_RF_hyperparams(Eq = training.severity ~ .,
                                  data = df_reg,
                                  type = "regression",
                                  mtry_grid = mtry_grid,
                                  ntree_grid = ntree_grid,
                                  n_folds = n_folds)
  
  best_reg <- tune_reg$best_params
  
  # Fit regression model
  rf_reg <- ranger(training.severity ~ .,
                   data = df_reg,
                   num.trees = best_reg$num.trees,
                   mtry = best_reg$mtry,
                   importance = "impurity")
  
  # Output
  output <- list(rf_class = rf_class,
                 rf_reg = rf_reg,
                 threshold = threshold_opt,
                 tuning = list(class = tune_class,
                               reg = tune_reg),
                 info = list(type = "fit_RF_zero_inflated output",
                             n_folds = n_folds))
  
  return(output)
}




## Function: predict_RF_zero_inflated --------------------------------------
#' 
#' This function generates predictions from a fitted hurdle-type Random Forest model
#' created with "fit_RF_zero_inflated".
#'
#' @param covariates Data frame or matrix. Predictor variables for prediction.
#' @param RF_zero_inflated List. Output object from "fit_RF_zero_inflated".
#'
#' @details
#' The prediction procedure follows three steps:
#'
#' 1. Classification prediction
#'  - Predicts probability of non-zero values using the classification model
#'  - Extracts probability of class "1" (non-zero)
#'  - Converts probabilities into binary predictions using the stored threshold
#'
#' 2. Regression prediction
#'  - Predicts continuous severity using the regression model
#'
#' 3. Final prediction
#'  - Expected prediction: probability × regression
#'  - Hard prediction: regression if presence is predicted, 0 otherwise
#'
#' @return A data.frame with:
#'  - prob_nonzero: Probability of being non-zero
#'  - pred_nonzero: Binary prediction (0/1)
#'  - pred_reg: Regression prediction (continuous)
#'  - pred: Expected prediction (probability × regression)
#'  - pred_hard: Hard prediction using threshold
#'


predict_RF_zero_inflated <- function(covariates,
                                     RF_zero_inflated) {
  # Input check
  if (is.null(RF_zero_inflated$info$type) ||
      RF_zero_inflated$info$type != "fit_RF_zero_inflated output") {
    stop("The provided model object is not a valid 'fit_RF_zero_inflated' output.")
  }
  
  # Extract models
  rf_class <- RF_zero_inflated$rf_class
  rf_reg <- RF_zero_inflated$rf_reg
  threshold <- RF_zero_inflated$threshold # Or can be set at 0.5
  
  # Predict probability of being non-zero 
  pred_class <- predict(rf_class, data = covariates)$predictions
  
  # If the classification RF outputs probabilities, keep the prob for class '1'
  if (is.matrix(pred_class)) {
    prob_nonzero <- pred_class[, "1"]
  } else {
    # In some ranger versions, predictions might be directly 0/1 labels
    prob_nonzero <- as.numeric(pred_class)
  }
  
  # Convert probabilities into binary decisions using the threshold
  pred_nonzero <- ifelse(prob_nonzero >= threshold, 1, 0)
  
  # Predict continuous severity for all observations ---
  pred_reg <- predict(rf_reg, data = covariates)$predictions
  
  # Combine both predictions:
  pred_final <- prob_nonzero * pred_reg # Expectation
  pred_hard <- ifelse(pred_nonzero == 1, pred_reg, 0)  # Hard prediction

  # Output
  return(data.frame(prob_nonzero = prob_nonzero,
                    pred_nonzero = pred_nonzero,
                    pred_reg = pred_reg,
                    pred = pred_final,
                    pred_hard = pred_hard))
}



## Function: tune_RF_hyperparams -------------------------------------------
#' 
#' This function performs hyperparameter tuning for a Random Forest model
#' using cross-validation via the "caret" package. It evaluates combinations of:
#'  - mtry (number of variables randomly sampled at each split)
#'  - num.trees (number of trees)
#'
#' The function supports both classification and regression tasks.
#'
#' @param Eq Formula. Model formula (e.g., y ~ .).
#' @param data Data frame. Dataset containing response and predictors.
#' @param type Character. Either "classification" or "regression".
#' @param mtry_grid Integer vector. Grid of "mtry" values to test.
#' @param ntree_grid Integer vector. Grid of number of trees to test.
#' @param n_folds Integer. Number of folds for cross-validation. Default is 5.
#'
#' @details
#' The tuning procedure:
#'
#' 1. Grid search
#'  - Generates all combinations of "mtry" and "num.trees"
#'  - k-fold cross-validation via "caret::train"
#'  - Model evaluation
#'      Classification: Accuracy is maximized
#'      Regression: RMSE is minimized 
#'
#' The best parameter set is selected based on the highest score.
#'
#' @return A list with:
#'  - best_params: Data frame with optimal "mtry" and "num.trees"
#'  - tuning_table: Data frame with all tested combinations and their scores
#'


tune_RF_hyperparams <- function(Eq, 
                                data, 
                                type = c("classification", "regression"),
                                mtry_grid = 1:10,
                                ntree_grid = seq(500, 1500, 100),
                                n_folds = 5) {
  type <- match.arg(type)
  
  # Define grid
  param_grid <- expand.grid(mtry = mtry_grid, num.trees = ntree_grid)
  param_grid$score <- NA
  
  # Cross-validation
  ctrl <- trainControl(method = "cv", number = n_folds)
  
  # Grid search
  for (i in seq_len(nrow(param_grid))) {
    
    mtry_i <- param_grid$mtry[i]
    ntree_i <- param_grid$num.trees[i]
    
    if (type == "classification") {
      
      model_cv <- train(Eq, 
                        data = data,
                        method = "ranger",
                        trControl = ctrl,
                        tuneGrid = data.frame(mtry = mtry_i, splitrule = "gini", min.node.size = 1),
                        num.trees = ntree_i,
                        importance = "none")
      
      # Maximize accuracy
      param_grid$score[i] <- model_cv$results$Accuracy
      
    } else {
      model_cv <- train(Eq, 
                        data = data,
                        method = "ranger",
                        trControl = ctrl,
                        tuneGrid = data.frame(mtry = mtry_i, splitrule = "variance", min.node.size = 1),
                        num.trees = ntree_i,
                        importance = "none")
      
      # Minimize RMSE
      param_grid$score[i] <- - min(model_cv$results$RMSE)
    }
  }
   
  # Get best configuration
  best_idx <- which.max(param_grid$score)
  
  best_params <- param_grid[best_idx, c("mtry", "num.trees")]
  
  # Output
  return(list(best_params = best_params,
              tuning_table = param_grid))
}




## Function: Select_training_data ------------------------------------------
#' 
#' This function creates a logical index to select training data for a given
#' prediction year. Training data include
#'  - All observations from years different than the prediction year
#'  - A random subset of observations from the prediction year, sampled by group
#'
#' @param vect_id Vector. Identifier for sampling units (e.g. field).
#' @param vect_year Vector. Year associated with each observation.
#' @param sampling_group Vector. Grouping variable used for stratified sampling.
#' @param predict_year Numeric. Year to predict.
#' @param prop Numeric. Proportion of unique IDs to keep within each group
#' for the prediction year (between 0 and 1).
#' @param randomness_seed Integer or NA. Optional seed for reproducibility.
#'
#' @details
#' The selection procedure:
#'
#' 1. All observations from years different than the prediction year are included
#'
#' 2. For the prediction year:
#'  - Observations are grouped using "sampling_group"
#'  - A proportion ("prop") of unique IDs is randomly selected within each group
#'  - All observations corresponding to selected IDs are included
#'
#' @return A logical vector indicating which observations are selected for training
#'

Select_training_data <- function(vect_id,
                                 vect_year,
                                 sampling_group,
                                 predict_year,
                                 prop,
                                 randomness_seed = NA) {
  
  # Set seed for reproducibility
  if (!is.na(randomness_seed)) {set.seed(randomness_seed)}
  
  # Initialization
  train_index <- rep(FALSE, length(vect_id))
  
  # Include data from all other years
  train_index[vect_year != predict_year] <- TRUE
  
  # Sampling within prediction year
  select_same_year <- vect_year == predict_year
  
  #For each group
  for (grp in unique(sampling_group)) {
    
    grp_id <- unique(vect_id[select_same_year & sampling_group == grp])
    
    if (length(grp_id) > 0) {
      
      # Keep the same proportion
      n_keep <- ceiling(prop * length(grp_id))
      keep_ids <- sample(grp_id, n_keep)
      
      # TRUE for selected ids
      train_index[vect_id %in% keep_ids] <- TRUE
    }
  }
  
  return(train_index)
}



## Function: compute_spatial_neighbourhood ---------------------------------
#'
#' This function computes spatial neighbourhood meata-features based on predictions.
#' For each observation, it calculates the mean prediction within a given
#' spatial radius, restricted to observations from the same year.
#'
#' @param geometry Optional sfc object (sf geometry column).
#' @param vect_pred Numeric vector. Predictions for each observation.
#' @param vect_year Vector. Year associated with each observation.
#' @param thresholds Numeric vector. Distance thresholds (km).
#'
#' @details
#' The function:
#'  - Computes pairwise distances between observations
#'  - Defines neighbours within each distance threshold
#'  - Restricts neighbours to the same year
#'  - Computes mean prediction within neighbourhood
#'
#' @return A matrix with one column per threshold (predN.*)
#'


compute_spatial_neighbourhood <- function(geometry,
                                          vect_pred,
                                          vect_year,
                                          dist_thresholds = c(0, 10, 25, 50, 100, 200, 1000)) {
  
  # Spatial distances (in km)
  dist_S <- matrix(round(as.numeric(st_distance(geometry)) / 1000), 
                   nrow = length(geometry))
  
  # Temporal distances (0 or 1)
  dist_T <- outer(vect_year, vect_year, "!=") + 0
  
  # Initialize matrix to store smoothed predictions
  pred_neighbourhood <- NULL
  
  # Compute average prediction of neighbors within distance thresholds
  for (threshold in dist_thresholds) {
    neighbours <- sapply(1:length(vect_pred),
                         function(i) (1:length(vect_pred))[dist_S[i,] <= threshold & dist_T[i,]==0])
    
    pred_neighbourhood <- cbind(pred_neighbourhood, 
                                as.numeric(lapply(neighbours, function(u) mean(vect_pred[u], na.rm = TRUE))))
  }
  
  colnames(pred_neighbourhood) <- paste("predN.", dist_thresholds, sep = "")
  
  return(pred_neighbourhood)
}



## Function: get_hurdle_metrics --------------------------------------------
#'
#' This function computes performance metrics for a hurdle / zero-inflated model
#' based on predictions for total, training, and validation datasets.
#'
#' It calculates:
#' - Root Mean Squared Error (RMSE)
#' - Coefficient of determination (R-squared)
#'
#' @param pred_total Data frame or list. Predictions for the full dataset.
#' Must contain a numeric vector pred.
#' @param pred_train Data frame or list. Predictions for the training subset.
#' Must contain a numeric vector pred.
#' @param pred_valid Data frame or list. Predictions for the validation subset.
#' Must contain a numeric vector pred.
#' @param Y_total Numeric vector. Observed values for the full dataset.
#' @param Y_train Numeric vector. Observed values for the training subset.
#' @param Y_valid Numeric vector. Observed values for the validation subset.
#'
#' @details
#' The function:
#'
#' 1. Computes RMSE for each dataset (total, train, validation):
#' - Measures average prediction error magnitude
#'
#' 2. Computes R-squared using a linear model:
#' - Regression of predicted values on observed values
#' - Indicates proportion of explained variance
#'
#' Predictions must be aligned with corresponding observed values
#'
#' @return A list with two named vectors:
#' - RMSE: RMSE for total, train, and validation datasets
#' - R2: R-squared for total, train, and validation datasets
#'

get_hurdle_metrics <- function(pred_total, pred_train, pred_valid,
                               Y_total, Y_train, Y_valid) {
  
  # Quick check
  stopifnot(length(Y_total) == nrow(pred_total))
  stopifnot(length(Y_train) == nrow(pred_train))
  stopifnot(length(Y_valid) == nrow(pred_valid))
  
  # RMSE 
  RMSE_total <- sqrt(mean((Y_total - pred_total$pred)^2, na.rm = TRUE))
  RMSE_train <- sqrt(mean((Y_train - pred_train$pred)^2, na.rm = TRUE))
  RMSE_valid <- sqrt(mean((Y_valid - pred_valid$pred)^2, na.rm = TRUE))
  
  # R2
  R2_total <- summary(lm(pred_total$pred ~ Y_total))$r.squared
  R2_train <- summary(lm(pred_train$pred ~ Y_train))$r.squared
  R2_valid <- summary(lm(pred_valid$pred ~ Y_valid))$r.squared
  
  # Output 
  return(list(
    RMSE = list(total = RMSE_total,
                train = RMSE_train,
                valid = RMSE_valid),
    
    R2 = list(total = R2_total,
              train = R2_train,
              valid = R2_valid)
  ))
}




##  Function: stacked_hurdle_model ------------------------------------------
#'
#' This function fits a stacked hurdle model.
#' It first fits a baseline RF model, then augments covariates with spatial
#' neighbourhood meta-features and fits a second RF model.
#'
#' @param train_select Logical vector. Indicates training observations.
#' @param X_total Data frame or matrix. Covariates for the full dataset.
#' @param Y_total Numeric vector. Response variable for the full dataset.
#' @param vect_year Vector. Year associated with each observation.
#' @param geometry sfc object. Spatial geometries (sf column).
#' @param dist_thresholds Numeric vector. Distance thresholds (km) used to compute
#' spatial neighbourhood features.
#'
#' @details
#' The function:
#'
#' 1. Splits data into training and validation sets using train_select
#'
#' 2. Fits an initial hurdle model (Model I) using original covariates
#'
#' 3. Computes predictions and performance metrics (initial model)
#'
#' 4. Builds mean predictions within spatial neighbourhoods 
#'    restricted to same-year observations as spatial meta-features
#'
#' 5. Fits a second hurdle model (Model F) 
#'
#' 6. Computes predictions and performance metrics for the final model
#'
#' @return A list containing:
#' - metrics_I, model_I: metrics and model for initial model
#' - metrics_F, model_F: metrics and model for final stacked model
#' - predictions (total, train, valid) for both models
#' - Y_train, Y_valid, Y_total
#'

stacked_hurdle_model <- function(train_select,
                                 X_total,
                                 Y_total,
                                 vect_year,
                                 geometry,
                                 dist_thresholds = c(0, 10, 25, 50, 100, 200, 1000)) {
  
  # Training and validation split
  valid_select <- !train_select
  
  # Response variables 
  Y_train <- Y_total[train_select]
  Y_valid <- Y_total[valid_select]
  
  # Covariates 
  X_train <- X_total[train_select, ]
  
  # Fit Model I 
  hurdle_mod_I <- fit_RF_zero_inflated(training.severity = Y_train,
                                       training.covariates = X_train,
                                       mtry_grid = 1:10,
                                       ntree_grid = seq(500, 1500, 200),
                                       n_folds = 5)
  
  # Predictions Model I
  pred_total_I <- predict_RF_zero_inflated(X_total, hurdle_mod_I)
  pred_train_I <- pred_total_I[train_select, ]
  pred_valid_I <- pred_total_I[valid_select, ]
  
  # Metrics Model I
  hurdle_metrics_I <- get_hurdle_metrics(pred_total_I, pred_train_I, pred_valid_I,
                                         Y_total, Y_train, Y_valid)
  
  # Spatial meta-features
  pred_neighbourhood <- compute_spatial_neighbourhood(geometry = geometry,
                                                      vect_pred = pred_total_I$pred,
                                                      vect_year = vect_year,
                                                      dist_thresholds = dist_thresholds)
  
  X_total_F <- cbind(X_total, pred_neighbourhood)
  X_train_F <- X_total_F[train_select, ]
  
  # Fit Model F (Final)
  hurdle_mod_F <- fit_RF_zero_inflated(training.severity = Y_train,
                                       training.covariates = X_train_F,
                                       mtry_grid = 1:10,
                                       ntree_grid = seq(500, 1500, 200),
                                       n_folds = 5)
  
  # Predictions Model F
  pred_total_F <- predict_RF_zero_inflated(X_total_F, hurdle_mod_F)
  pred_train_F <- pred_total_F[train_select, ]
  pred_valid_F <- pred_total_F[valid_select, ]
  
  # Metrics Model F
  hurdle_metrics_F <- get_hurdle_metrics(pred_total_F, pred_train_F, pred_valid_F,
                                         Y_total, Y_train, Y_valid)
  
  # Outputs
  return(list(metrics_I     = hurdle_metrics_I,
              model_I       = hurdle_mod_I,
              
              metrics_F     = hurdle_metrics_F,
              model_F       = hurdle_mod_F,
              
              pred_total_I  = pred_total_I,
              pred_train_I  = pred_train_I,
              pred_valid_I  = pred_valid_I,
              
              pred_total_F  = pred_total_F,
              pred_train_F  = pred_train_F,
              pred_valid_F  = pred_valid_F,
              
              Y_train       = Y_train,
              Y_valid       = Y_valid,
              Y_total       = Y_total))
}