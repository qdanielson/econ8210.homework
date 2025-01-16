# Install and load the xgboost package
# install.packages("xgboost")
# install.packages("caret")

# Other useful libraries
library(xgboost)  # Machine Learning
library(dplyr)    # Data manipulation
library(caret)    # For creating a train-test split
library(pROC)
library(ggplot2)
library(ParBayesianOptimization)


# # set working directory
# setwd("C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/replication")
# 
# # Set seed for reproducibility
# set.seed(123)
# 
# panel <- readRDS("panel.RDS")
# 
# ## arrange for grouping
# panel <- panel %>% arrange(gender)
# 
# # Add unique identifier to the dataset
# panel$id <- 1:nrow(panel)
# 
# # tabulate count by gender
# table(panel$gender)
# 
# # Result:
# # 0     1 
# # 46731 76269
# 
# f_indices <- 1:46731
# m_indices <- 46732:123000
# 
# panel_f <- panel[f_indices,]
# panel_m <- panel[m_indices,]
# 
# # Pull
# y <- panel$y
# 
# # Train-test split
# train_index <- createDataPartition(y, p = 0.8, list = FALSE)
# 
# # Exclude the "id" and "y" column for training and testing
# X_train <- panel[train_index, !(names(panel) %in% c("id","y"))]
# y_train <- y[train_index]
# X_test <- panel[-train_index, !(names(panel) %in% c("id","y"))]
# y_test <- y[-train_index]
# 
# # Convert to xgb.DMatrix format
# dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
# dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
# 
# # Female-specific train-test split
# train_index_f <- train_index[train_index %in% f_indices]
# test_index_f <- f_indices[!(f_indices %in% train_index_f)]
# 
# # Exclude "id" and "y" columns for female data
# X_train_f <- panel[train_index_f, !(names(panel) %in% c("id", "y"))]
# y_train_f <- y[train_index_f]
# X_test_f <- panel[test_index_f, !(names(panel) %in% c("id", "y"))]
# y_test_f <- y[test_index_f]
# 
# # Convert female data to xgb.DMatrix format
# dtrain_f <- xgb.DMatrix(data = as.matrix(X_train_f), label = y_train_f)
# dtest_f <- xgb.DMatrix(data = as.matrix(X_test_f), label = y_test_f)
# 
# # Male-specific train-test split
# train_index_m <- train_index[train_index %in% m_indices]
# test_index_m <- m_indices[!(m_indices %in% train_index_m)]
# 
# # Exclude "id" and "y" columns for male data
# X_train_m <- panel[train_index_m, !(names(panel) %in% c("id", "y"))]
# y_train_m <- y[train_index_m]
# X_test_m <- panel[test_index_m, !(names(panel) %in% c("id", "y"))]
# y_test_m <- y[test_index_m]
# 
# # Convert male data to xgb.DMatrix format
# dtrain_m <- xgb.DMatrix(data = as.matrix(X_train_m), label = y_train_m)
# dtest_m <- xgb.DMatrix(data = as.matrix(X_test_m), label = y_test_m)
# 
# # Define the search space
# search_space <- list(
#   max_depth = c(1, 100),
#   eta = c(exp(-9), exp(0)),  # Log-uniform
#   subsample = c(0.5, 1),
#   colsample_bytree = c(0.5, 1),
#   min_child_weight = c(exp(-2), exp(3))
# )
# 
# # Define the objective function
# objective_function <- function(max_depth, eta, subsample, colsample_bytree, min_child_weight) {
#   # Transform log-scale parameters back to normal scale
#   eta <- exp(eta)
#   min_child_weight <- exp(min_child_weight)
#   
#   # Define XGBoost parameters
#   params <- list(
#     objective = "binary:logistic",
#     eval_metric = "logloss",
#     max_depth = as.integer(max_depth),
#     eta = eta,
#     subsample = subsample,
#     colsample_bytree = colsample_bytree,
#     min_child_weight = min_child_weight
#   )
#   
#   # Perform cross-validation
#   cv_results <- xgb.cv(
#     params = params,
#     data = dtrain,
#     nrounds = 100,
#     nfold = 5,
#     verbose = FALSE,
#     early_stopping_rounds = 10
#   )
#   
#   # Return the logloss
#   list(Score = min(cv_results$evaluation_log$test_logloss_mean),
#        Pred = min(cv_results$evaluation_log$test_logloss_mean))
# }
# 
# # Run Bayesian Optimization
# opt_results <- bayesOpt(
#   FUN = objective_function,
#   bounds = search_space,
#   initPoints = 10,  # Initial random points
#   iters.n = 500,     # Iterations
#   acq = "ucb"       # Acquisition function
# )
# 
# # Pull optimal parameters from Bayesian Optimization results
# final_params <- function(best_params) {
#   list(
#     objective = "binary:logistic",          # Binary classification
#     eval_metric = "logloss",               # Use logloss as evaluation metric
#     max_depth = as.integer(best_params$max_depth),
#     eta = best_params$eta,
#     colsample_bytree = best_params$colsample_bytree,
#     subsample = best_params$subsample,
#     min_child_weight = best_params$min_child_weight
#   )
# }
# 
# # Use the best parameters for pooled data
# best_params_pooled <- getBestPars(opt_results)  # Best parameters for pooled model
# params_pooled <- final_params(best_params_pooled)
# 
# # Train the final XGBoost model for pooled data
# final_model_pooled <- xgb.train(
#   params = params_pooled,
#   data = dtrain,
#   nrounds = 100,  # Number of boosting rounds
#   watchlist = list(train = dtrain, test = dtest),
#   early_stopping_rounds = 10,
#   verbose = 1
# )
# 
# # ----------------------
# # Repeat Bayesian Optimization and Model Training for Male and Female Subsamples
# # ----------------------
# 
# # Female-Specific Optimization and Training
# objective_function_female <- function(max_depth, eta, subsample, colsample_bytree, min_child_weight) {
#   eta <- exp(eta)
#   min_child_weight <- exp(min_child_weight)
#   
#   params <- list(
#     objective = "binary:logistic",
#     eval_metric = "logloss",
#     max_depth = as.integer(max_depth),
#     eta = eta,
#     subsample = subsample,
#     colsample_bytree = colsample_bytree,
#     min_child_weight = min_child_weight
#   )
#   
#   cv_results <- xgb.cv(
#     params = params,
#     data = dtrain_f,
#     nrounds = 100,
#     nfold = 5,
#     verbose = FALSE,
#     early_stopping_rounds = 10
#   )
#   
#   list(Score = min(cv_results$evaluation_log$test_logloss_mean),
#        Pred = min(cv_results$evaluation_log$test_logloss_mean))
# }
# 
# opt_results_female <- bayesOpt(
#   FUN = objective_function_female,
#   bounds = search_space,
#   initPoints = 10,
#   iters.n = 500,
#   acq = "ucb"
# )
# 
# best_params_female <- getBestPars(opt_results_female)  # Best parameters for female model
# params_female <- final_params(best_params_female)
# 
# # Train final model for female data
# final_model_female <- xgb.train(
#   params = params_female,
#   data = dtrain_f,
#   nrounds = 100,
#   watchlist = list(train = dtrain_f, test = dtest_f),
#   early_stopping_rounds = 10,
#   verbose = 1
# )
# 
# # Male-Specific Optimization and Training
# objective_function_male <- function(max_depth, eta, subsample, colsample_bytree, min_child_weight) {
#   eta <- exp(eta)
#   min_child_weight <- exp(min_child_weight)
#   
#   params <- list(
#     objective = "binary:logistic",
#     eval_metric = "logloss",
#     max_depth = as.integer(max_depth),
#     eta = eta,
#     subsample = subsample,
#     colsample_bytree = colsample_bytree,
#     min_child_weight = min_child_weight
#   )
#   
#   cv_results <- xgb.cv(
#     params = params,
#     data = dtrain_m,
#     nrounds = 100,
#     nfold = 5,
#     verbose = FALSE,
#     early_stopping_rounds = 10
#   )
#   
#   list(Score = min(cv_results$evaluation_log$test_logloss_mean),
#        Pred = min(cv_results$evaluation_log$test_logloss_mean))
# }
# 
# opt_results_male <- bayesOpt(
#   FUN = objective_function_male,
#   bounds = search_space,
#   initPoints = 10,
#   iters.n = 500,
#   acq = "ucb"
# )
# 
# best_params_male <- getBestPars(opt_results_male)  # Best parameters for male model
# params_male <- final_params(best_params_male)
# 
# # Train final model for male data
# final_model_male <- xgb.train(
#   params = params_male,
#   data = dtrain_m,
#   nrounds = 100,
#   watchlist = list(train = dtrain_m, test = dtest_m),
#   early_stopping_rounds = 10,
#   verbose = 1
# )

# Save best parameters for pooled, female, and male models as .RDS files
saveRDS(best_params_pooled, "best_params_pooled.RDS")
saveRDS(best_params_female, "best_params_female.RDS")
saveRDS(best_params_male, "best_params_male.RDS")


