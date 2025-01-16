# Install and load the xgboost package
# install.packages("xgboost")
# install.packages("caret")

# Other useful libraries
library(xgboost)  # Machine Learning
library(dplyr)    # Data manipulation
library(caret)    # For creating a train-test split
library(pROC)
library(ggplot2)

# set working directory
setwd("C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/replication")

# Set seed for reproducibility
set.seed(123)

# Load tuned params:
params_pooled <- readRDS("best_params_pooled.RDS")
params_female <- readRDS("best_params_female.RDS")
params_male <- readRDS("best_params_male.RDS")

params_pooled$max_depth <- round(params_pooled$max_depth,0)
params_female$max_depth <- round(params_female$max_depth,0)
params_male$max_depth <- round(params_male$max_depth,0)

# Load data
panel <- readRDS("panel.RDS")

## arrange for grouping
panel <- panel %>% arrange(gender)

# Add unique identifier to the dataset
panel$id <- 1:nrow(panel)

# tabulate count by gender
table(panel$gender)

# Result:
# 0     1 
# 46731 76269

f_indices <- 1:46731
m_indices <- 46732:123000

panel_f <- panel[f_indices,]
panel_m <- panel[m_indices,]

# Pull
y <- panel$y

# Train-test split
train_index <- createDataPartition(y, p = 0.8, list = FALSE)

# Exclude the "id" and "y" column for training and testing
X_train <- panel[train_index, !(names(panel) %in% c("id","y"))]
y_train <- y[train_index]
X_test <- panel[-train_index, !(names(panel) %in% c("id","y"))]
y_test <- y[-train_index]

# Convert to xgb.DMatrix format
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

# Female-specific train-test split
train_index_f <- train_index[train_index %in% f_indices]
test_index_f <- f_indices[!(f_indices %in% train_index_f)]

# Exclude "id" and "y" columns for female data
X_train_f <- panel[train_index_f, !(names(panel) %in% c("id", "y"))]
y_train_f <- y[train_index_f]
X_test_f <- panel[test_index_f, !(names(panel) %in% c("id", "y"))]
y_test_f <- y[test_index_f]

# Convert female data to xgb.DMatrix format
dtrain_f <- xgb.DMatrix(data = as.matrix(X_train_f), label = y_train_f)
dtest_f <- xgb.DMatrix(data = as.matrix(X_test_f), label = y_test_f)

# Male-specific train-test split
train_index_m <- train_index[train_index %in% m_indices]
test_index_m <- m_indices[!(m_indices %in% train_index_m)]

# Exclude "id" and "y" columns for male data
X_train_m <- panel[train_index_m, !(names(panel) %in% c("id", "y"))]
y_train_m <- y[train_index_m]
X_test_m <- panel[test_index_m, !(names(panel) %in% c("id", "y"))]
y_test_m <- y[test_index_m]

# Convert male data to xgb.DMatrix format
dtrain_m <- xgb.DMatrix(data = as.matrix(X_train_m), label = y_train_m)
dtest_m <- xgb.DMatrix(data = as.matrix(X_test_m), label = y_test_m)

################################################################


# Train the pooled model
watchlist <- list(train = dtrain, test = dtest)
bst_pooled <- xgb.train(
  params = params_pooled,
  data = dtrain,
  nrounds = 500,      # Number of boosting rounds
  early_stopping_rounds = 10,
  watchlist = watchlist,
  verbose = 1
)

# Train the female model
watchlist_f <- list(train = dtrain_f, test = dtest_f)
bst_female <- xgb.train(
  params = params_female,
  data = dtrain_f,
  nrounds = 500,      # Number of boosting rounds
  early_stopping_rounds = 10,
  watchlist = watchlist_f,
  verbose = 1
)

# Train the male model
watchlist_m <- list(train = dtrain_m, test = dtest_m)
bst_male <- xgb.train(
  params = params_male,
  data = dtrain_m,
  nrounds = 500,      # Number of boosting rounds
  early_stopping_rounds = 10,
  watchlist = watchlist_m,
  verbose = 1
)

# Predict probabilities on the test data
pred_probs_pooled <- predict(bst_pooled, dtest)
pred_probs_female <- predict(bst_female, dtest_f)
pred_probs_male <- predict(bst_male, dtest_m)

# Convert probabilities to binary predictions (threshold = 0.5)
pred_labels_pooled <- ifelse(pred_probs_pooled > 0.2, 1, 0)
pred_labels_female <- ifelse(pred_probs_female > 0.2, 1, 0)
pred_labels_male <- ifelse(pred_probs_male > 0.2, 1, 0)

# Compute ROC curves for pooled, female, and male models
roc_curve_pooled <- roc(y_test, pred_probs_pooled)
roc_curve_female <- roc(y_test_f, pred_probs_female)
roc_curve_male <- roc(y_test_m, pred_probs_male)

# Plot the ROC curve for the pooled model
plot(roc_curve_pooled, col = "blue", lwd = 2, main = "ROC Curves for Pooled, Female, and Male Models", 
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line for reference

# Add the female model ROC curve
plot(roc_curve_female, col = "red", lwd = 2, add = TRUE)

# Add the male model ROC curve
plot(roc_curve_male, col = "green", lwd = 2, add = TRUE)

# Add a legend to distinguish between models
legend("bottomright", 
       legend = c("Pooled Model", "Female Model", "Male Model"), 
       col = c("blue", "red", "green"), 
       lwd = 2)

## Comparison plots

# Extract test IDs and genders for the pooled test set
test_ids <- panel[-train_index, "id"]
test_genders <- panel[-train_index, "gender"]

# Initialize the dataframe with test IDs and genders
results <- data.frame(
  id = test_ids,
  gender = test_genders,
  pooled_prob = pred_probs_pooled,  # Add probabilities from the pooled model
  gender_specific_prob = NA         # Initialize gender-specific probabilities as NA
)

# Populate gender-specific probabilities based on gender
results$gender_specific_prob <- ifelse(
  results$gender == 0, 
  pred_probs_female[match(results$id, panel[test_index_f, "id"])],  # Female-specific probs
  pred_probs_male[match(results$id, panel[test_index_m, "id"])]     # Male-specific probs
)

# Create the scatter plot with female points on top
ggplot() +
  # Male points (plotted first)
  geom_point(data = results[results$gender == 1, ],
             aes(x = pooled_prob, y = gender_specific_prob),
             color = "blue", alpha = 0.25) +  # Male points (blue)
  # Female points (plotted second)
  geom_point(data = results[results$gender == 0, ],
             aes(x = pooled_prob, y = gender_specific_prob),
             color = "red", alpha = 0.25) +  # Female points (red)
  # Dashed lines
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
  # Labels and legend
  labs(
    title = "Scatter Plot of Pooled vs Gender-Specific Probabilities",
    x = "Pooled Probability",
    y = "Gender-Specific Probability"
  ) +
  theme_minimal()  # Minimal theme
