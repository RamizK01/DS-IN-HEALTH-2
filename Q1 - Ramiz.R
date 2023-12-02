##########################################
### INITIAL DATA CLEANING/VAR CREATION ###
##########################################

df <- read.csv('cleaned_transfusion_data.csv')

library(glmnet)
set.seed(123)

# create df_bin, bin stands for binary classification
df <- df %>% mutate("RBC_tfsd" = case_when(`Total.24hr.RBC` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("PLT_tfsd" = case_when(`Total.24hr.Plt` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("FFP_tfsd" = case_when(`Total.24hr.FFP` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("Cryo_tfsd" = case_when(`Total.24hr.Cryo` == 0 ~ 0, TRUE ~ 1))
  
### THIS DF WILL ONLY WORK FOR RBC, NEED TO ADD NEW COLUMNS 64:66 FOR PLT, FFP, CRYO
df_bin <- df[,c(3:26,28,29,63)]

correlation_matrix <- cor(df_bin)
find_high_correlations <- function(cor_matrix, threshold) {
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
  high_correlations <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  return(data.frame(Variable1 = rownames(cor_matrix)[high_correlations[, 1]],
                    Variable2 = colnames(cor_matrix)[high_correlations[, 2]],
                    Correlation = cor_matrix[high_correlations]))
}

# Set your correlation threshold
threshold <- 0.7

# Find and display correlated variable pairs
high_corr_pairs <- find_high_correlations(correlation_matrix, threshold)
print(high_corr_pairs)

library(ggplot2)

# heatmap
ggplot(data = melt(correlation_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +  # Choose your color palette
  theme_minimal() +
  labs(title = "Correlation Matrix Heatmap")

df_bin[c(1, 2, 7:16, 25:26)] <- lapply(df_bin[c(1, 2, 7:16, 25:26)], factor)



################################################################
### LASSO MODEL FOR PREDICTING IF RBC TRANSFUSION WILL OCCUR ###
################################################################
library(pROC)
library(glmnet)
train.set <- sample(nrow(df_bin),round(nrow(df_bin)*.7))

# train predictor 
x.train <- model.matrix(RBC_tfsd ~.,df_bin)[train.set,-1]
# train response
y.train <- df_bin$RBC_tfsd[train.set]
# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial",
                      type.measure = "auc")

lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1,
                      lambda = lambda_optimal)

# predict
pred.lasso <- as.numeric(predict(lasso.model, 
                                 newx = model.matrix(RBC_tfsd ~.,df_bin)
                                 [-train.set,-1], s=cv.lasso$lambda.min,
                                 type = "response"))

# create ROC curve
myroc <- roc(RBC_tfsd ~ pred.lasso, data=df_bin[-train.set,])
auc.lasso <- myroc$auc
plot(myroc)

### WE SHOULD REALLY BOOTSTRAP THE LASSO ###
library(SparseLearner)
bstrap.lasso.model <- BRLasso(x, y, B = 5, Boots = 100, kfold = 10, seed = 0123)


################################################################
### LASSO MODEL FOR PREDICTING IF PLT TRANSFUSION WILL OCCUR ###
################################################################
library(pROC)
library(glmnet)
train.set <- sample(nrow(df_bin),round(nrow(df_bin)*.7))

# train predictor 
x.train <- model.matrix(PLT_tfsd ~.,df_bin)[train.set,-1]
# train response
y.train <- df_bin$PLT_tfsd[train.set]
# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial",
                      type.measure = "auc")

lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1,
                      lambda = lambda_optimal)

# predict
pred.lasso <- as.numeric(predict(lasso.model, 
                                 newx = model.matrix(PLT_tfsd ~.,df_bin)
                                 [-train.set,-1], s=cv.lasso$lambda.min,
                                 type = "response"))

# create ROC curve
myroc <- roc(PLT_tfsd ~ pred.lasso, data=df_bin[-train.set,])
auc.lasso <- myroc$auc
plot(myroc)

################################################################
### LASSO MODEL FOR PREDICTING IF FFP TRANSFUSION WILL OCCUR ###
################################################################
library(pROC)
library(glmnet)
train.set <- sample(nrow(df_bin),round(nrow(df_bin)*.7))

# train predictor 
x.train <- model.matrix(FFP_tfsd ~.,df_bin)[train.set,-1]
# train response
y.train <- df_bin$PLT_tfsd[train.set]
# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial",
                      type.measure = "auc")

lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1,
                      lambda = lambda_optimal)

# predict
pred.lasso <- as.numeric(predict(lasso.model, 
                                 newx = model.matrix(FFP_tfsd ~.,df_bin)
                                 [-train.set,-1], s=cv.lasso$lambda.min,
                                 type = "response"))

# create ROC curve
myroc <- roc(FFP_tfsd ~ pred.lasso, data=df_bin[-train.set,])
auc.lasso <- myroc$auc
plot(myroc)

#################################################################
### LASSO MODEL FOR PREDICTING IF CRYO TRANSFUSION WILL OCCUR ###
#################################################################
library(pROC)
library(glmnet)
train.set <- sample(nrow(df_bin),round(nrow(df_bin)*.7))

# train predictor 
x.train <- model.matrix(Cryo_tfsd ~.,df_bin)[train.set,-1]
# train response
y.train <- df_bin$PLT_tfsd[train.set]
# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial",
                      type.measure = "auc")

lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1,
                      lambda = lambda_optimal)

# predict
pred.lasso <- as.numeric(predict(lasso.model, 
                                 newx = model.matrix(FFP_tfsd ~.,df_bin)
                                 [-train.set,-1], s=cv.lasso$lambda.min,
                                 type = "response"))

# create ROC curve
myroc <- roc(FFP_tfsd ~ pred.lasso, data=df_bin[-train.set,])
auc.lasso <- myroc$auc
plot(myroc)


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
df_mvar[c(1, 2, 7:16, 25:26)] <- lapply(df_mvar[c(1, 2, 7:16, 25:26)], factor)

train.set <- sample(nrow(df_mvar),round(nrow(df_mvar)*.7))

# train predictor 
x.train <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~.,df_mvar)[train.set,-1]
# train response
y.train <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                        df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[train.set, -1]

# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "mgaussian",
                      type.measure = "mse")
lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "mgaussian", alpha = 1,
                      lambda = lambda_optimal + 1)


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


################################
### RANDOM FOREST REGRESSION ###
################################
library(randomForest)
library(ggplot2)

rf.fit1 <- randomForest(Total.24hr.RBC + Total.24hr.Plt + 
                       Total.24hr.FFP + Total.24hr.Cryo ~ ., data=df_mvar, ntree=1000,
                       keep.forest=FALSE, importance=TRUE)

ImpData <- as.data.frame(importance(rf.fit1))
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# remove variables to decrease MSE 
df_rfr <- df_mvar[,c(-5,-8,-10,-24,-12,-1)]

rf.fit2 <- randomForest(Total.24hr.RBC + Total.24hr.Plt + 
                         Total.24hr.FFP + Total.24hr.Cryo ~ ., data=df_rfr, ntree=1000,
                       keep.forest=TRUE, importance=TRUE)

single_tree <- getTree(rf.fit2, k = 1, labelVar=TRUE)

ggplot(single_tree) +
  geom_point(aes(x = `split var`, y = `split point`))


###########################################################
### DECISION REGRESSION TREE FOR RBC TRANSFUSION AMOUNT ###
###########################################################
library(caret)
library(rpart)
library(rpart.plot)

# Splitting the data into training and test sets
set.seed(123) # for reproducibility
index <- createDataPartition(df_mvar$Total.24hr.RBC, p = 0.8, list = FALSE)
train_set <- df_mvar[index, ]
test_set <- df_mvar[-index, ]

fit <- rpart(Total.24hr.RBC  ~ . - Total.24hr.Plt -
               Total.24hr.FFP - Total.24hr.Cryo,  
             method = "anova", data = train_set) 

predictions <- predict(fit, test_set)

# Actual values from the test set
actual_values <- test_set$Total.24hr.RBC

# Calculate MSE on the test set
mse <- postResample(pred = predictions, obs = actual_values)
mse <- unname(mse['RMSE'])

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.032572)

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

# Splitting the data into training and test sets
set.seed(123) # for reproducibility
index <- createDataPartition(df_mvar$Total.24hr.Plt, p = 0.8, list = FALSE)
train_set <- df_mvar[index, ]
test_set <- df_mvar[-index, ]

fit <- rpart(Total.24hr.Plt  ~ . - Total.24hr.RBC -
               Total.24hr.FFP - Total.24hr.Cryo,  
             method = "anova", data = train_set) 

predictions <- predict(fit, test_set)

# Actual values from the test set
actual_values <- test_set$Total.24hr.Plt

# Calculate MSE on the test set
mse <- postResample(pred = predictions, obs = actual_values)
mse <- unname(mse['RMSE'])

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.045798)

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
library(caret)
library(rpart)
library(rpart.plot)

# Splitting the data into training and test sets
set.seed(123) # for reproducibility
index <- createDataPartition(df_mvar$Total.24hr.Plt, p = 0.8, list = FALSE)
train_set <- df_mvar[index, ]
test_set <- df_mvar[-index, ]

fit <- rpart(Total.24hr.FFP  ~ . - Total.24hr.RBC -
               Total.24hr.Plt - Total.24hr.Cryo,  
             method = "anova", data = train_set) 

predictions <- predict(fit, test_set)

# Actual values from the test set
actual_values <- test_set$Total.24hr.FFP

# Calculate MSE on the test set
mse <- postResample(pred = predictions, obs = actual_values)
mse <- unname(mse['RMSE'])

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.076547)

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
library(caret)
library(rpart)
library(rpart.plot)

# Splitting the data into training and test sets
set.seed(123) # for reproducibility
index <- createDataPartition(df_mvar$Total.24hr.Cryo, p = 0.8, list = FALSE)
train_set <- df_mvar[index, ]
test_set <- df_mvar[-index, ]

fit <- rpart(Total.24hr.Cryo  ~ . - Total.24hr.RBC -
               Total.24hr.Plt - Total.24hr.FFP,  
             method = "anova", data = train_set) 

predictions <- predict(fit, test_set)

# Actual values from the test set
actual_values <- test_set$Total.24hr.Cryo

# Calculate MSE on the test set
mse <- postResample(pred = predictions, obs = actual_values)
mse <- unname(mse['RMSE'])

# Plot 
plot(fit, uniform = TRUE)
text(fit, use.n = TRUE, cex = .7) 

printcp(fit)

# Plot the CP table to help in visualizing the right CP to choose
plotcp(fit)

# Choose a CP value (for example, based on the CP table)
optimal_cp <- fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"]

# Prune the tree
pruned_fit <- prune(fit, cp = 0.053521)

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
