##########################################
### INITIAL DATA CLEANING/VAR CREATION ###
##########################################

df <- read.csv('cleaned_transfusion_data.csv')

library(glmnet)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pROC)
library(caret)
library(rpart)
library(rpart.plot)
library(car)
set.seed(123)

# create df_bin, bin stands for binary classification
df <- df %>% mutate("RBC_tfsd" = case_when(`Total.24hr.RBC` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("PLT_tfsd" = case_when(`Total.24hr.Plt` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("FFP_tfsd" = case_when(`Total.24hr.FFP` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("Cryo_tfsd" = case_when(`Total.24hr.Cryo` == 0 ~ 0, TRUE ~ 1))
  
df_bin <- df[,c(3:26,28,29,63:66)]
df_bin[c(3:6,17:24)] <- scale(df_bin[c(3:6,17:24)])

####################################################
### Meeting Assumptions for Lasso Classification ###
####################################################

# 1) Outcome is binary: Yes
# 2) Independence of Datapoints: Yes via study ID
# 3) Multicollinearity:
vif(lm(RBC_tfsd ~ ., data = df_bin))
# Remove Height and Weight due to high ViF, keep BMI
df_bin <- df_bin[,c(-3,-4)]
# Check ViF again
vif(lm(RBC_tfsd ~ ., data = df_bin))
# Remove Pre_PT bc Pre_INR measures it better and Hb because Hct is percentage of Hb in blood
df_bin <- df_bin[,c(-16,-19)]
# Check ViF again
vif(lm(RBC_tfsd ~ ., data = df_bin))
# 4) Scaled variables
df_bin[c(1,2,5:14,21:26)] <- lapply(df_bin[c(1,2,5:14,21:26)], factor)
# 5) Continuous Predictors are normally distributed
# 6) Continuous Predictors have homoscedasticity
# 7) Removal of outliers
Q1 <- quantile(df_bin$LAS.score, 0.25)
Q3 <- quantile(df_bin$LAS.score, 0.75)
IQR <- Q3 - Q1
LB <- Q1 - 1.5 * IQR
UB <- Q3 + 1.5 * IQR

# Remove outliers
df_bin <- df_bin %>% filter(LAS.score >= LB & LAS.score <= UB)

Q1 <- quantile(df_bin$Pre_INR, 0.25)
Q3 <- quantile(df_bin$Pre_INR, 0.75)
IQR <- Q3 - Q1
LB <- Q1 - 1.5 * IQR
UB <- Q3 + 1.5 * IQR

# Remove outliers
df_bin <- df_bin %>% filter(Pre_INR >= LB & Pre_INR <= UB)

Q1 <- quantile(df_bin$Pre_PTT, 0.25)
Q3 <- quantile(df_bin$Pre_PTT, 0.75)
IQR <- Q3 - Q1
LB <- Q1 - 1.5 * IQR
UB <- Q3 + 1.5 * IQR

# Remove outliers
df_bin <- df_bin %>% filter(Pre_PTT >= LB & Pre_PTT <= UB)


################################################################
### LASSO MODEL FOR PREDICTING IF RBC TRANSFUSION WILL OCCUR ###
################################################################
df_bin1 <- df_bin[,c(-24,-25,-26)]
set.seed(123) # for reproducibility
n_iterations <- 100
auc_values <- numeric(n_iterations)

for(i in 1:n_iterations) {
  # Sample the data
  train.set <- sample(nrow(df_bin1), round(nrow(df_bin1) * 0.7), replace = TRUE)
  
  # Train predictor
  x.train <- model.matrix(RBC_tfsd ~., df_bin1)[train.set, -1]
  # Train response
  y.train <- df_bin1$RBC_tfsd[train.set]
  
  # Perform cross-validation for lambda selection
  cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial", type.measure = "auc")
  
  lambda_optimal <- cv.lasso$lambda.min
  lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1, lambda = lambda_optimal)
  
  # Predict
  pred.lasso <- as.numeric(predict(lasso.model, 
                                   newx = model.matrix(RBC_tfsd ~., df_bin1)[-train.set, -1], 
                                   s = cv.lasso$lambda.min, 
                                   type = "response"))
  
  # Create ROC curve
  myroc <- roc(response = df_bin1$RBC_tfsd[-train.set], predictor = pred.lasso)
  auc_values[i] <- myroc$auc
}

# Calculate the average AUC
average_auc <- mean(auc_values)
se_auc <- sd(auc_values) / sqrt(n_iterations)
list(average_auc = average_auc, SE = se_auc)

################################################################
### LASSO MODEL FOR PREDICTING IF PLT TRANSFUSION WILL OCCUR ###
################################################################
df_bin2 <- df_bin[,c(-23,-25,-26)]
set.seed(123) # for reproducibility
n_iterations <- 100
auc_values <- numeric(n_iterations)

for(i in 1:n_iterations) {
  # Sample the data
  train.set <- sample(nrow(df_bin2), round(nrow(df_bin2) * 0.7), replace = TRUE)
  
  # Train predictor
  x.train <- model.matrix(PLT_tfsd ~., df_bin2)[train.set, -1]
  # Train response
  y.train <- df_bin2$PLT_tfsd[train.set]
  
  # Perform cross-validation for lambda selection
  cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial", type.measure = "auc")
  
  lambda_optimal <- cv.lasso$lambda.min
  lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1, lambda = lambda_optimal)
  
  # Predict
  pred.lasso <- as.numeric(predict(lasso.model, 
                                   newx = model.matrix(PLT_tfsd ~., df_bin2)[-train.set, -1], 
                                   s = cv.lasso$lambda.min, 
                                   type = "response"))
  
  # Create ROC curve
  myroc <- roc(response = df_bin2$PLT_tfsd[-train.set], predictor = pred.lasso)
  auc_values[i] <- myroc$auc
}

average_auc <- mean(auc_values)
se_auc <- sd(auc_values) / sqrt(n_iterations)
list(average_auc = average_auc, SE = se_auc)
################################################################
### LASSO MODEL FOR PREDICTING IF FFP TRANSFUSION WILL OCCUR ###
################################################################
df_bin3 <- df_bin[,c(-23,-24,-26)]
set.seed(123) # for reproducibility
n_iterations <- 100
auc_values <- numeric(n_iterations)

for(i in 1:n_iterations) {
  # Sample the data
  train.set <- sample(nrow(df_bin3), round(nrow(df_bin3) * 0.7), replace = TRUE)
  
  # Train predictor
  x.train <- model.matrix(FFP_tfsd ~., df_bin3)[train.set, -1]
  # Train response
  y.train <- df_bin3$FFP_tfsd[train.set]
  
  # Perform cross-validation for lambda selection
  cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial", type.measure = "auc")
  
  lambda_optimal <- cv.lasso$lambda.min
  lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1, lambda = lambda_optimal)
  
  # Predict
  pred.lasso <- as.numeric(predict(lasso.model, 
                                   newx = model.matrix(FFP_tfsd ~., df_bin3)[-train.set, -1], 
                                   s = cv.lasso$lambda.min, 
                                   type = "response"))
  
  # Create ROC curve
  myroc <- roc(response = df_bin3$FFP_tfsd[-train.set], predictor = pred.lasso)
  auc_values[i] <- myroc$auc
}

average_auc <- mean(auc_values[1:50])
se_auc <- sd(auc_values) / sqrt(n_iterations)
list(average_auc = average_auc, SE = se_auc)
#################################################################
### LASSO MODEL FOR PREDICTING IF CRYO TRANSFUSION WILL OCCUR ###
#################################################################
df_bin4 <- df_bin[,c(-23,-24,-25)]
set.seed(123) # for reproducibility
n_iterations <- 100
auc_values <- numeric(n_iterations)

for(i in 1:n_iterations) {
  # Sample the data
  train.set <- sample(nrow(df_bin4), round(nrow(df_bin4) * 0.7), replace = TRUE)
  
  # Train predictor
  x.train <- model.matrix(Cryo_tfsd ~., df_bin4)[train.set, -1]
  # Train response
  y.train <- df_bin4$Cryo_tfsd[train.set]
  
  # Perform cross-validation for lambda selection
  cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial", type.measure = "auc")
  
  lambda_optimal <- cv.lasso$lambda.min
  lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1, lambda = lambda_optimal)
  
  # Predict
  pred.lasso <- as.numeric(predict(lasso.model, 
                                   newx = model.matrix(Cryo_tfsd ~., df_bin4)[-train.set, -1], 
                                   s = cv.lasso$lambda.min, 
                                   type = "response"))
  
  # Create ROC curve
  myroc <- roc(response = df_bin4$Cryo_tfsd[-train.set], predictor = pred.lasso)
  auc_values[i] <- myroc$auc
}

average_auc <- mean(auc_values[1:50])
se_auc <- sd(auc_values) / sqrt(n_iterations)
list(average_auc = average_auc, SE = se_auc)

#########################################################
### WHY USE A MULTIVARIATE MODEL OVER INDIVIDUAL ONES ###
#########################################################
# 1) Correlation matrix 
library(ggplot2)
df_mvar_outcomes <- df_mvar[,27:30]
cor_matrix <- cor(df_mvar_outcomes)
plot(cor_matrix)

# Convert the correlation matrix to a long-format data frame
library(reshape2)
cor_matrix_long <- melt(cor_matrix)

# Create a ggplot heatmap with values inside each square
ggplot(data = cor_matrix_long, aes(Var1, Var2, fill = value, label = round(value, 2))) +
  geom_tile() +
  geom_text(aes(label = ifelse(value != 1, as.character(round(value, 2)), "")), vjust = 1) + # Display values (exclude diagonal)
  scale_fill_gradient(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################################
### MULTIVARIATE LASSO MODEL TO PREDICT AMOUNT OF TRANSFUSION ###
#################################################################
set.seed(123)
# create df_mvar, nvar stands for dataset for multivariate regression
df_mvar <- df[,c(3:26,28,29,57,60:62)]
# remove multicollinearity vars
df_mvar <- df_mvar[,c(-3,-4, -18, -21)]
# check vif
vif(lm(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~ ., data = df_mvar))
# good

df_mvar[c(1, 2, 5:14, 21:22)] <- lapply(df_mvar[c(1, 2, 5:14, 21:22)], factor)

# Function to remove outliers from a numeric vector
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  LB <- Q1 - 1.5 * IQR
  UB <- Q3 + 1.5 * IQR
  x_filtered <- ifelse(x < LB | x > UB, NA, x)
  return(x_filtered)
}

excluded_columns <- c(1,2,5:14,21:26)

# Apply the function to numeric columns excluding the specified ones
df_mvar <- df_mvar %>%
  mutate(across(
    .cols = -all_of(excluded_columns),
    .fns = list(remove_outliers)
  ))

# Remove rows with NA values (outliers)
df_mvar <- df_mvar %>%
  na.omit()

# START OF MODEL #
set.seed(123)
train.set <- sample(nrow(df_mvar),round(nrow(df_mvar)*.7))

# train predictor 
x.train <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~.,df_mvar)[train.set,-1]
# train response
y.train <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                        df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[train.set, -1]

# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, nfolds = 3, alpha = 1, family = "mgaussian",
                      type.measure = "mse")
lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "mgaussian", alpha = 1,
                      lambda = lambda_optimal)


plot(cv.lasso)
title("MSE values at different values of log(lambda)", line = 3, cex.main = 0.9)

coefficients(lasso.model)

# PREDICT AND CALCULATE MSE
x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[-train.set, -1]

pred.lasso <- predict(lasso.model, newx = x.test, s = cv.lasso$lambda.min, type = "response")
pred.lasso <- drop(pred.lasso)

# Actual values for the test set
y.test <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                         df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[-train.set, -1]

# Calculate MSE
mse <- colMeans((pred.lasso - y.test)^2)
mse

################################
### RANDOM FOREST REGRESSION ###
################################



###########################################################
### DECISION REGRESSION TREE FOR RBC TRANSFUSION AMOUNT ###
###########################################################
df_tree1 <- df[,c(3:26,28,29,57,60:62)]
df_tree1[c(1,2,7:16,25:26)] <- lapply(df_tree1[c(1,2,7:16,25:26)], factor)

# Splitting the data into training and test sets
set.seed(123) # for reproducibility
num_iterations <- 100
mse_values <- numeric(num_iterations) # To store MSE values

for(i in 1:num_iterations) {
  # Bootstrapped sampling for training data
  boot_indexes <- sample(nrow(df_tree1), size = round(nrow(df_tree1) * 0.7), replace = TRUE)
  test_indexes <- setdiff(1:nrow(df_tree1), boot_indexes)
  
  train_set <- df_tree1[boot_indexes, ]
  test_set <- df_tree1[test_indexes, ]
  
  fit <- rpart(Total.24hr.RBC  ~ . - Total.24hr.Plt - Total.24hr.FFP - Total.24hr.Cryo,  
               method = "anova", data = df_tree1)
  
  predictions <- predict(fit, test_set)
  
  # Actual values from the test set
  actual_values <- test_set$Total.24hr.RBC
  
  # Calculate MSE on the test set
  mse_values[i] <- postResample(pred = predictions, obs = actual_values)['RMSE']^2
}

mean_mse <- mean(mse_values)
se_mse <- sd(mse_values) / sqrt(length(mse_values))

list(mean_mse = mean_mse, se_mse = se_mse)

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.030878)

# Plot the pruned tree
rpart.plot(pruned_fit)

# Make predictions with the pruned tree
pruned_predictions <- predict(pruned_fit, test_set)

# Calculate MSE with the pruned tree predictions
pruned_mse <- postResample(pred = pruned_predictions, obs = actual_values)
pruned_mse <- unname(pruned_mse['RMSE'])
pruned_mse
# PRUNING THE TREE REDUCED MSE BY 0.002 NICE #

###########################################################
### DECISION REGRESSION TREE FOR PLT TRANSFUSION AMOUNT ###
###########################################################
library(caret)
library(rpart)
library(rpart.plot)

set.seed(123) # for reproducibility
num_iterations <- 100
mse_values <- numeric(num_iterations) # To store MSE values

for(i in 1:num_iterations) {
  # Bootstrapped sampling for training data
  boot_indexes <- sample(nrow(df_tree1), size = round(nrow(df_tree1) * 0.7), replace = TRUE)
  test_indexes <- setdiff(1:nrow(df_tree1), boot_indexes)
  
  train_set <- df_tree1[boot_indexes, ]
  test_set <- df_tree1[test_indexes, ]
  
  fit <- rpart(Total.24hr.Plt  ~ . - Total.24hr.RBC -
                  Total.24hr.FFP - Total.24hr.Cryo,  
                method = "anova", data = train_set)
  
  predictions <- predict(fit, test_set)
  
  # Actual values from the test set
  actual_values <- test_set$Total.24hr.Plt
  
  # Calculate MSE on the test set
  mse_values[i] <- postResample(pred = predictions, obs = actual_values)['RMSE']^2
}

mean_mse <- mean(mse_values)
se_mse <- sd(mse_values) / sqrt(length(mse_values))

list(mean_mse = mean_mse, se_mse = se_mse)

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.018)

# Plot the pruned tree
rpart.plot(pruned_fit)

# Make predictions with the pruned tree
pruned_predictions <- predict(pruned_fit, test_set)

# Calculate MSE with the pruned tree predictions
pruned_mse <- postResample(pred = pruned_predictions, obs = actual_values)
pruned_mse <- unname(pruned_mse['RMSE'])
pruned_mse

###########################################################
### DECISION REGRESSION TREE FOR FFP TRANSFUSION AMOUNT ###
###########################################################
set.seed(123) # for reproducibility
num_iterations <- 100
mse_values <- numeric(num_iterations) # To store MSE values

for(i in 1:num_iterations) {
  # Bootstrapped sampling for training data
  boot_indexes <- sample(nrow(df_tree1), size = round(nrow(df_tree1) * 0.7), replace = TRUE)
  test_indexes <- setdiff(1:nrow(df_tree1), boot_indexes)
  
  train_set <- df_tree1[boot_indexes, ]
  test_set <- df_tree1[test_indexes, ]
  
  fit <- rpart(Total.24hr.FFP  ~ . - Total.24hr.RBC -
                  Total.24hr.Plt - Total.24hr.Cryo,  
                method = "anova", data = train_set) 
  
  predictions <- predict(fit, test_set)
  
  # Actual values from the test set
  actual_values <- test_set$Total.24hr.FFP
  
  # Calculate MSE on the test set
  mse_values[i] <- postResample(pred = predictions, obs = actual_values)['RMSE']^2
}

mean_mse <- mean(mse_values)
se_mse <- sd(mse_values) / sqrt(length(mse_values))

list(mean_mse = mean_mse, se_mse = se_mse)


# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.010000)

# Plot the pruned tree
rpart.plot(pruned_fit)

# Make predictions with the pruned tree
pruned_predictions <- predict(pruned_fit, test_set)

# Calculate MSE with the pruned tree predictions
pruned_mse <- postResample(pred = pruned_predictions, obs = actual_values)
pruned_mse <- unname(pruned_mse['RMSE'])
pruned_mse

############################################################
### DECISION REGRESSION TREE FOR CRYO TRANSFUSION AMOUNT ###
############################################################
set.seed(123) # for reproducibility
num_iterations <- 100
mse_values <- numeric(num_iterations) # To store MSE values

for(i in 1:num_iterations) {
  # Bootstrapped sampling for training data
  boot_indexes <- sample(nrow(df_tree1), size = round(nrow(df_tree1) * 0.7), replace = TRUE)
  test_indexes <- setdiff(1:nrow(df_tree1), boot_indexes)
  
  train_set <- df_tree1[boot_indexes, ]
  test_set <- df_tree1[test_indexes, ]
  
  fit <- rpart(Total.24hr.Cryo  ~ . - Total.24hr.RBC -
                 Total.24hr.Plt - Total.24hr.FFP,  
               method = "anova", data = train_set)  
  
  predictions <- predict(fit, test_set)
  
  # Actual values from the test set
  actual_values <- test_set$Total.24hr.Cryo
  
  # Calculate MSE on the test set
  mse_values[i] <- postResample(pred = predictions, obs = actual_values)['RMSE']^2
}

mean_mse <- mean(mse_values)
se_mse <- sd(mse_values) / sqrt(length(mse_values))

list(mean_mse = mean_mse, se_mse = se_mse)



# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.017310)

# Plot the pruned tree
rpart.plot(pruned_fit)

# Make predictions with the pruned tree
pruned_predictions <- predict(pruned_fit, test_set)

# Calculate MSE with the pruned tree predictions
pruned_mse <- postResample(pred = pruned_predictions, obs = actual_values)
pruned_mse <- unname(pruned_mse['RMSE'])
pruned_mse

####################################################
### ELASTIC NET MODEL FOR RBC TRANSFUSION AMOUNT ###
####################################################

elastic_net_alpha = 0.5  # For example, 0.5 can be used for a balanced Elastic Net

# Perform cross-validation for lambda selection with Elastic Net
cv.elastic_net <- cv.glmnet(x.train, y.train, alpha = elastic_net_alpha, family = "mgaussian",
                            type.measure = "mse")

lambda_optimal <- cv.elastic_net$lambda.min
elastic_net_model <- glmnet(x.train, y.train, family = "mgaussian", alpha = elastic_net_alpha,
                            lambda = lambda_optimal)

# Plotting the CV results
plot(cv.elastic_net)
title("MSE values at different values of log(lambda)", line = 3, cex.main = 0.9)

# Coefficients of the Elastic Net Model
coefficients(elastic_net_model)

# Predict and Calculate MSE
x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[-train.set, -1]
pred.elastic_net <- predict(elastic_net_model, newx = x.test, s = lambda_optimal, type = "response")

# Actual values for the test set
y.test <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                         df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[-train.set, -1]

# Ensure dimensions match
pred.elastic_net <- drop(pred.elastic_net)

# Calculate MSE
mse <- colMeans((pred.elastic_net - y.test)^2)

# Output the MSE
mse

########################################################
### RIDGE REGRESSION MODEL FOR AMOUNT OF TRANSFUSION ###
########################################################

# Assuming you have already set the seed and prepared df_mvar, train.set, x.train, and y.train

# Ridge Regression parameter: alpha set to 0
ridge_alpha = 0  # For Ridge Regression

# Perform cross-validation for lambda selection with Ridge Regression
cv.ridge <- cv.glmnet(x.train, y.train, alpha = ridge_alpha, family = "mgaussian",
                      type.measure = "mse")

lambda_optimal <- cv.ridge$lambda.min
ridge_model <- glmnet(x.train, y.train, family = "mgaussian", alpha = ridge_alpha,
                      lambda = lambda_optimal)

# Plotting the CV results
plot(cv.ridge)
title("MSE values at different values of log(lambda)", line = 3, cex.main = 0.9)

# Coefficients of the Ridge Model
coefficients(ridge_model)

# Predict and Calculate MSE
x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[-train.set, -1]
pred.ridge <- predict(ridge_model, newx = x.test, s = lambda_optimal, type = "response")

# Actual values for the test set
y.test <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                         df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[-train.set, -1]

# Ensure dimensions match
pred.ridge <- drop(pred.ridge)

# Calculate MSE
mse <- colMeans((pred.ridge - y.test)^2)

# Output the MSE
mse

###########################################################
### XGBOOST GRADIENT BOOSTING FOR AMOUNT OF TRANSFUSION ###
###########################################################
# Load necessary libraries
library(xgboost)
library(readr)
library(caret) # for data splitting
library(Metrics) # for calculating MSE

# Read the dataset
data <- read_csv("df_mvar.csv")

# Function to prepare data, train the model, and evaluate MSE
train_and_evaluate_xgboost <- function(data, target_variable) {
  # Excluding the target variable and other outcome variables from predictors
  predictors <- setdiff(names(data), c(target_variable, "Total.24hr.RBC", "Total.24hr.Plt", "Total.24hr.FFP", "Total.24hr.Cryo"))
  
  # Split data into train and test sets
  set.seed(123) # for reproducibility
  trainIndex <- createDataPartition(data[[target_variable]], p = .8, list = FALSE)
  data_train <- data[trainIndex, ]
  data_test <- data[-trainIndex, ]
  
  # Prepare data for xgboost
  dtrain <- xgb.DMatrix(data = as.matrix(data_train[predictors]), label = data_train[[target_variable]])
  dtest <- xgb.DMatrix(data = as.matrix(data_test[predictors]), label = data_test[[target_variable]])
  
  # Parameters for the xgboost model
  params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    eta = 0.3,
    max_depth = 4,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1
  )
  
  # Training the model
  xgb_model <- xgb.train(params, dtrain, nrounds = 100)
  
  # Predicting and calculating MSE
  predictions <- predict(xgb_model, dtest)
  mse_value <- mse(data_test[[target_variable]], predictions)
  
  return(list(model = xgb_model, mse = mse_value))
}

# Train models for each outcome and evaluate MSE
results_rbc <- train_and_evaluate_xgboost(data, "Total.24hr.RBC")
results_plt <- train_and_evaluate_xgboost(data, "Total.24hr.Plt")
results_ffp <- train_and_evaluate_xgboost(data, "Total.24hr.FFP")
results_cryo <- train_and_evaluate_xgboost(data, "Total.24hr.Cryo")

results_rbc
results_plt
results_ffp
results_cryo


############################################
### BOOTSTRAPPING LASSO REGRESSION MODEL ###
############################################

set.seed(123) # for reproducibility
num_iterations <- 100
num_outcomes <- 4 # Assuming there are 4 outcomes
mse_values <- array(dim = c(num_iterations, num_outcomes)) # To store MSE values

for(i in 1:num_iterations) {
  boot_indexes <- sample(nrow(df_mvar), size = round(nrow(df_mvar) * 0.7), replace = TRUE)
  test_indexes <- setdiff(1:nrow(df_mvar), boot_indexes)
  
  x.train <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[boot_indexes, -1]
  y.train <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP + df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[boot_indexes, -1]
  
  cv.lasso <- cv.glmnet(x.train, y.train, nfolds = 3, alpha = 1, family = "mgaussian", type.measure = "mse")
  lasso.model <- glmnet(x.train, y.train, family = "mgaussian", alpha = 1, lambda = cv.lasso$lambda.min)
  
  x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[test_indexes, -1]
  y.test <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP + df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[test_indexes, -1]
  
  pred.lasso <- predict(lasso.model, newx = x.test, s = cv.lasso$lambda.min)
  pred.lasso <- matrix(pred.lasso, ncol = ncol(y.test))  # Reshape if necessary
  
  if(is.matrix(y.test) && is.matrix(pred.lasso)) {
    # Calculate MSE for each column (variable) and take the mean
    mse_values[i] <- mean(colMeans((pred.lasso - y.test)^2))
  } else {
    mse_values[i] <- mean((pred.lasso - y.test)^2)
  }
  
  for(j in 1:num_outcomes) {
    mse_values[i, j] <- mean((pred.lasso[, j] - y.test[, j])^2)
  }
}

# Calculating the mean and standard error of MSE for each outcome
mean_mse <- colMeans(mse_values)
se_mse <- apply(mse_values, 2, function(x) sd(x) / sqrt(length(x)))

list(mean_mse = mean_mse, se_mse = se_mse)





