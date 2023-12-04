df <- readxl::read_excel("transfusion data.xlsx")

# index columns of interest
df <- df[,c(1,3,4:15,25,26,28:41,46,47,49,50,57,60:65,77:92,116,117)]

### Rename and convert variables to binary ###
library(dplyr)
df <- rename(df, "Gender" = "Gender (male)")
df[,28] <- is.character(df[,28])

# Columns to convert T/F to 1 and 0
cols_convert <- c(4,9:18,28,29,30)
for (col_index in cols_convert) {
  df[, col_index] <- ifelse(df[, col_index] == "TRUE", 1, 0)
}

# rename variables 
df <- df %>% mutate(Type = case_when(Type == "Bilateral" ~ 0,
                                     Type == "Single Right Lung" | Type == "Single Left Lung" ~ 1))

# CREATE VARIABLE FOR SURVIVAL TIME IN DAYS
df$`OR Date` <- as.Date(df$`OR Date`, format = "%Y-%m-%d")
df$DEATH_DATE <- as.Date(strptime(df$DEATH_DATE, format = "%d-%b-%Y"))

df$Survival_Time <- as.numeric(difftime(df$DEATH_DATE, df$`OR Date`, units = "days"))
df$Survival_Time[is.na(df$DEATH_DATE)] <- NA
# categorize survival time?
# we can create a censored group (> 669 amount of days) then categorize all groups less than 669
# 669 being the max survival time reported 
max(df$Survival_Time, na.rm=T)

create_report(df)

### Dealing with missing data ###

# All these columns missing data should be 0
df[, 42:56] <- replace(df[, 42:56], is.na(df[, 42:56]), 0)
# Now evaluate missing data

colSums(is.na(df))
# Fibrinogen: 97.4% missing, remove
df <- df[,-26]
# Death Date and Survival Time: 83.33% missing, no worries this is systematic
# LAS score is 6.25% missing, imputation?
imp_df <- df[,c(3:14,19)]
library(missForest)
imp_df <- as.data.frame(imp_df)
imp_df <- missForest(imp_df)$ximp

ggplot() +
  geom_point(data = imp_df, aes(x = seq_along(`LAS score`), y = `LAS score`), color = "red") +
  geom_point(data = df, aes(x = seq_along(`LAS score`), y = `LAS score`), color = "black") +
  ggtitle("Original vs Imputed Data") +
  xlab("Index") +
  ylab("LAS score") +
  theme_minimal()

df$`LAS score` <- imp_df$`LAS score`

# Pre_PTT: 1 missing, mean impute
mean_pre_ptt <- mean(df$Pre_PTT, na.rm = TRUE)
df$Pre_PTT <- ifelse(is.na(df$Pre_PTT), mean_pre_ptt, df$Pre_PTT)

### Creating 24hr blood product variables ###
df <- df %>% mutate('Total 24hr Plt' = Intra_Platelets + `Plt 0-24hrs`) %>%
  mutate('Total 24hr FFP' = `Intra_Fresh Frozen Plasma` + `FFP 0-24hrs`) %>%
  mutate('Total 24hr Cryo' = `Intra_Cryoprecipitate` + `Cryo 0-24hrs`)

### K-MEANS CLUSTERING ###
library(tidyverse)
library(scales)
library(cluster)
library(factoextra)

# Assuming your data is in a DataFrame named df
# Selecting the specified columns
selected_columns <- c(2:39, 56:61)
selected_data <- df[, selected_columns]

# Factorize binary columns and scale numerical columns
selected_data <- selected_data %>%
  mutate_if(~n_distinct(.) == 2 & all(. %in% c(0, 1)), as.factor) %>%
  mutate_if(is.numeric, scale)

# Remove columns with any NaN values
cleaned_data <- selected_data %>% select_if(~sum(is.na(.)) == 0)

# Determine the optimal number of clusters using the Elbow Method
fviz_nbclust(cleaned_data, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow Method")

# Applying K-means clustering with the chosen number of clusters (e.g., 4)
set.seed(42) # for reproducibility
kmeans_result <- kmeans(cleaned_data, centers = 4, nstart = 25)

# Adding cluster labels to the data
cleaned_data$Cluster <- as.factor(kmeans_result$cluster)

# Performing PCA for visualization
pca_result <- prcomp(cleaned_data, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data$Cluster <- cleaned_data$Cluster

# Visualizing the clusters based on PCA
fviz_cluster(list(data = pca_data, cluster = cleaned_data$Cluster))

### EXPORTING FINAL CLEANED DATA ###
write.csv(df, file = "cleaned_transfusion_data.csv", row.names = FALSE)

#############################################
### Q1 INITIAL DATA CLEANING/VAR CREATION ###
#############################################
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
n_iterations <- 50
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


#######################################################
# Load necessary libraries
library(e1071)
library(caret)
library(pROC)
cleaned_transfusion_data<- read.csv("cleaned_transfusion_data.csv")

# Splitting the data into training and testing sets
set.seed(123)  

#take out some variables 
SVM_df <- cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62) )]

index <- createDataPartition(SVM_df$Massive.Transfusion, p = 0.8, list = FALSE)

trainData <- SVM_df[index, ]
testData <- SVM_df[-index, ]

# Factor and Continuous Variables Handling
cols_to_factor <- c(1:2, 7:16, 25:27)
continuous_vars <- c("Height", "Weight", "Age", "BMI", "LAS.score", "Pre_Hb", "Pre_Hct", 
                     "Pre_Platelets", "Pre_PT", "Pre_INR", "Pre_PTT", "Pre_Creatinine")

# Adjust factor levels and scale continuous variables
for (var in cols_to_factor) {
  trainData[[var]] <- factor(trainData[[var]])
  testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
}
train_scaled <- trainData
test_scaled <- testData
train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])

# Hyperparameter Tuning and Cross-Validation
tune_result <- tune(svm, `Massive.Transfusion` ~ ., data = train_scaled, 
                    kernel = "radial", 
                    ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))

best_model <- tune_result$best.model

# Train SVM model with best parameters
svm_model <- svm(`Massive.Transfusion` ~ ., data = train_scaled, type = 'C-classification', 
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)

# Predict on test data
svm_pred <- predict(svm_model, test_scaled, probability = T)

# Evaluate model
confusionMatrix(svm_pred, test_scaled$Massive.Transfusion)


# Extract the probabilities for the positive class
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]

# Calculate ROC
roc_response <- roc(test_scaled$Massive.Transfusion, svm_pred_prob)
plot(roc_response, main = "Massive Transfusion ROC Curve")

#AUC
auc_value <- roc_response$auc


################## SVM for 24 hour values ######################## 
#Binary variables
cleaned_transfusion_data$binary.Total.24hr.RBC<- as.integer(cleaned_transfusion_data$Total.24hr.RBC > 0, 1, 0)
cleaned_transfusion_data$binary.Total.24hr.Plt<- as.integer(cleaned_transfusion_data$Total.24hr.Plt > 0, 1, 0)
cleaned_transfusion_data$binary.Total.24hr.FFP<- as.integer(cleaned_transfusion_data$Total.24hr.FFP > 0, 1, 0)
cleaned_transfusion_data$binary.Total.24hr.Cryo<- as.integer(cleaned_transfusion_data$Total.24hr.Cryo > 0, 1, 0)


##### Binary Total 24hr RBC ########

#take out some variables 
SVM_df_2<-cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62), -(64:66))]

index <- createDataPartition(SVM_df_2$binary.Total.24hr.RBC, p = 0.8, list = FALSE)

trainData <- SVM_df_2[index, ]
testData <- SVM_df_2[-index, ]

# Factor and Continuous Variables Handling
cols_to_factor <- c(1:2, 7:16, 25:26, 28)
continuous_vars <- c("Height", "Weight", "Age", "BMI", "LAS.score", "Pre_Hb", "Pre_Hct", 
                     "Pre_Platelets", "Pre_PT", "Pre_INR", "Pre_PTT", "Pre_Creatinine")

# Adjust factor levels and scale continuous variables
for (var in cols_to_factor) {
  trainData[[var]] <- factor(trainData[[var]])
  testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
}
train_scaled <- trainData
test_scaled <- testData
train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])

# Hyperparameter Tuning and Cross-Validation
tune_result <- tune(svm, binary.Total.24hr.RBC ~ ., data = train_scaled, 
                    kernel = "radial", 
                    ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))

best_model <- tune_result$best.model

# Train SVM model with best parameters
svm_model <- svm(binary.Total.24hr.RBC ~ ., data = train_scaled, type = 'C-classification', 
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)

# Predict on test data
svm_pred <- predict(svm_model, test_scaled, probability = T)

# Evaluate model
confusionMatrix(svm_pred, test_scaled$binary.Total.24hr.RBC)

# Extract the probabilities for the positive class
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]

# Calculate ROC
roc_response <- roc(test_scaled$binary.Total.24hr.RBC, svm_pred_prob)
plot(roc_response, main = "Total 24hr RBC ROC Curve")

#AUC
auc_value <- auc(roc_response)
print(auc_value)

##### Binary Total 24hr Plt ########


#take out some variables 
SVM_df_3 <- cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62), -63, -(65:66))]

index <- createDataPartition(SVM_df_3$binary.Total.24hr.Plt, p = 0.8, list = FALSE)

trainData <- SVM_df_3[index, ]
testData <- SVM_df_3[-index, ]

# Factor and Continuous Variables Handling
cols_to_factor <- c(1:2, 7:16, 25:26, 28)
continuous_vars <- c("Height", "Weight", "Age", "BMI", "LAS.score", "Pre_Hb", "Pre_Hct", 
                     "Pre_Platelets", "Pre_PT", "Pre_INR", "Pre_PTT", "Pre_Creatinine")

# Adjust factor levels and scale continuous variables
for (var in cols_to_factor) {
  trainData[[var]] <- factor(trainData[[var]])
  testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
}
train_scaled <- trainData
test_scaled <- testData
train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])

# Hyperparameter Tuning and Cross-Validation
tune_result <- tune(svm, binary.Total.24hr.Plt ~ ., data = train_scaled, 
                    kernel = "radial", 
                    ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))

best_model <- tune_result$best.model

# Train SVM model with best parameters
svm_model <- svm(binary.Total.24hr.Plt ~ ., data = train_scaled, type = 'C-classification', 
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)

# Predict on test data
svm_pred <- predict(svm_model, test_scaled, probability = T)

# Evaluate model
confusionMatrix(svm_pred, test_scaled$binary.Total.24hr.Plt)

# Extract the probabilities for the positive class
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]

# Calculate ROC
roc_response <- roc(test_scaled$binary.Total.24hr.Plt, svm_pred_prob)
plot(roc_response, main = "Total 24hr Plt ROC Curve")

#AUC
auc_value <- auc(roc_response)
print(auc_value)

##### binary.Total.24hr.FFP ########

#take out some variables 
SVM_df_4 <- cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62), -(63:64), -66)]

index <- createDataPartition(SVM_df_4$binary.Total.24hr.FFP, p = 0.8, list = FALSE)

trainData <- SVM_df_4[index, ]
testData <- SVM_df_4[-index, ]

# Factor and Continuous Variables Handling
cols_to_factor <- c(1:2, 7:16, 25:26, 28)
continuous_vars <- c("Height", "Weight", "Age", "BMI", "LAS.score", "Pre_Hb", "Pre_Hct", 
                     "Pre_Platelets", "Pre_PT", "Pre_INR", "Pre_PTT", "Pre_Creatinine")

# Adjust factor levels and scale continuous variables
for (var in cols_to_factor) {
  trainData[[var]] <- factor(trainData[[var]])
  testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
}
train_scaled <- trainData
test_scaled <- testData
train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])

# Hyperparameter Tuning and Cross-Validation
tune_result <- tune(svm, binary.Total.24hr.FFP ~ ., data = train_scaled, 
                    kernel = "radial", 
                    ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))

best_model <- tune_result$best.model

# Train SVM model with best parameters
svm_model <- svm(binary.Total.24hr.FFP ~ ., data = train_scaled, type = 'C-classification', 
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)

# Predict on test data
svm_pred <- predict(svm_model, test_scaled, probability = T)

# Evaluate model
confusionMatrix(svm_pred, test_scaled$binary.Total.24hr.FFP)

# Extract the probabilities for the positive class
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]

# Calculate ROC
roc_response <- roc(test_scaled$binary.Total.24hr.FFP, svm_pred_prob)
plot(roc_response, main = "Total 24hr FFP ROC Curve")

#AUC
auc_value <- auc(roc_response)
print(auc_value)

##### binary.Total.24hr.Cryo ########

#take out some variables 
SVM_df_5 <- cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62), -(63:65))]

index <- createDataPartition(SVM_df_5$binary.Total.24hr.Cryo, p = 0.8, list = FALSE)

trainData <- SVM_df_5[index, ]
testData <- SVM_df_5[-index, ]

# Factor and Continuous Variables Handling
cols_to_factor <- c(1:2, 7:16, 25:26, 28)
continuous_vars <- c("Height", "Weight", "Age", "BMI", "LAS.score", "Pre_Hb", "Pre_Hct", 
                     "Pre_Platelets", "Pre_PT", "Pre_INR", "Pre_PTT", "Pre_Creatinine")

# Adjust factor levels and scale continuous variables
for (var in cols_to_factor) {
  trainData[[var]] <- factor(trainData[[var]])
  testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
}
train_scaled <- trainData
test_scaled <- testData
train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])

# Hyperparameter Tuning and Cross-Validation
tune_result <- tune(svm, binary.Total.24hr.Cryo ~ ., data = train_scaled, 
                    kernel = "radial", 
                    ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))

best_model <- tune_result$best.model

# Train SVM model with best parameters
svm_model <- svm(binary.Total.24hr.Cryo ~ ., data = train_scaled, type = 'C-classification', 
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma)

# Predict on test data
predictions <- predict(svm_model, test_scaled)

# Evaluate model
confusionMatrix(predictions, test_scaled$binary.Total.24hr.Cryo)

# Extract the probabilities for the positive class
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]

# Calculate ROC
roc_response <- roc(test_scaled$binary.Total.24hr.Cryo, svm_pred_prob)
plot(roc_response, main = "Total 24hr Cryo ROC Curve")

#AUC
auc_value <- auc(roc_response)
print(auc_value)

###########################
# BOOTSTRAP RBC MODEL df2 #

library(caret)
library(e1071)
library(pROC)

set.seed(123)  # Set a random seed for reproducibility

n_iterations <- 100
auc_values <- numeric(n_iterations)
valid_iterations <- 0

for (i in 1:n_iterations) {
  # Bootstrap resampling
  resampled_index <- sample(nrow(SVM_df_2), size = round(0.8 * nrow(SVM_df_2)), replace = TRUE)
  trainData <- SVM_df_2[resampled_index, ]
  testData <- SVM_df_2[-resampled_index, ]
  
  # Check if all factor variables have at least two levels in the training data
  sufficient_data <- all(sapply(cols_to_factor, function(var) {
    length(levels(factor(trainData[[var]]))) >= 2
  }))
  
  if (!sufficient_data) {
    next  # Skip this iteration if not all factor variables have sufficient data
  }
  
  # Adjust factor levels and scale continuous variables
  for (var in cols_to_factor) {
    trainData[[var]] <- factor(trainData[[var]])
    testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
  }
  train_scaled <- trainData
  test_scaled <- testData
  train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
  test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])
  
  # Hyperparameter Tuning and Cross-Validation
  tune_result <- tune(svm, binary.Total.24hr.RBC ~ ., data = train_scaled, 
                      kernel = "radial", 
                      ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))
  best_model <- tune_result$best.model
  
  # Train SVM model with best parameters
  svm_model <- svm(binary.Total.24hr.RBC ~ ., data = train_scaled, type = 'C-classification', 
                   kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)
  
  # Predict on test data
  svm_pred <- predict(svm_model, test_scaled, probability = T)
  
  # Extract the probabilities for the positive class
  svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
  
  # Calculate ROC
  roc_response <- roc(test_scaled$binary.Total.24hr.RBC, svm_pred_prob)
  
  # Store AUC
  auc_values[i] <- roc_response$auc
  valid_iterations <- valid_iterations + 1
}

# Calculate the mean AUC and its standard error
mean_auc <- mean(auc_values[1:valid_iterations])
se_auc <- sd(auc_values[1:valid_iterations]) / sqrt(valid_iterations)

list(mean_auc = mean_auc, se_auc = se_auc)


##########################
# BOOTSTRAP Plt MODEL df3#

set.seed(123)  # Set a random seed for reproducibility

n_iterations <- 100
auc_values <- numeric(n_iterations)
valid_iterations <- 0

for (i in 1:n_iterations) {
  # Bootstrap resampling
  resampled_index <- sample(nrow(SVM_df_3), size = round(0.8 * nrow(SVM_df_3)), replace = TRUE)
  trainData <- SVM_df_3[resampled_index, ]
  testData <- SVM_df_3[-resampled_index, ]
  
  # Check if all factor variables have at least two levels in the training data
  sufficient_data <- all(sapply(cols_to_factor, function(var) {
    length(levels(factor(trainData[[var]]))) >= 2
  }))
  
  if (!sufficient_data) {
    next  # Skip this iteration if not all factor variables have sufficient data
  }
  
  # Adjust factor levels and scale continuous variables
  for (var in cols_to_factor) {
    trainData[[var]] <- factor(trainData[[var]])
    testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
  }
  train_scaled <- trainData
  test_scaled <- testData
  train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
  test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])
  
  # Hyperparameter Tuning and Cross-Validation
  tune_result <- tune(svm, binary.Total.24hr.Plt ~ ., data = train_scaled, 
                      kernel = "radial", 
                      ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))
  best_model <- tune_result$best.model
  
  # Train SVM model with best parameters
  svm_model <- svm(binary.Total.24hr.Plt ~ ., data = train_scaled, type = 'C-classification', 
                   kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)
  
  # Predict on test data
  svm_pred <- predict(svm_model, test_scaled, probability = T)
  
  # Extract the probabilities for the positive class
  svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
  
  # Calculate ROC
  roc_response <- roc(test_scaled$binary.Total.24hr.Plt, svm_pred_prob)
  
  # Store AUC
  auc_values[i] <- roc_response$auc
  valid_iterations <- valid_iterations + 1
}

# Calculate the mean AUC and its standard error
mean_auc <- mean(auc_values[1:valid_iterations])
se_auc <- sd(auc_values[1:valid_iterations]) / sqrt(valid_iterations)

list(mean_auc = mean_auc, se_auc = se_auc)


###########################
# BOOTSTRAP FFP MODEL df4 #

set.seed(123)  # Set a random seed for reproducibility

n_iterations <- 100
auc_values <- numeric(n_iterations)
valid_iterations <- 0

for (i in 1:n_iterations) {
  # Bootstrap resampling
  resampled_index <- sample(nrow(SVM_df_4), size = round(0.8 * nrow(SVM_df_4)), replace = TRUE)
  trainData <- SVM_df_4[resampled_index, ]
  testData <- SVM_df_4[-resampled_index, ]
  
  # Check if all factor variables have at least two levels in the training data
  sufficient_data <- all(sapply(cols_to_factor, function(var) {
    length(levels(factor(trainData[[var]]))) >= 2
  }))
  
  if (!sufficient_data) {
    next  # Skip this iteration if not all factor variables have sufficient data
  }
  
  # Adjust factor levels and scale continuous variables
  for (var in cols_to_factor) {
    trainData[[var]] <- factor(trainData[[var]])
    testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
  }
  train_scaled <- trainData
  test_scaled <- testData
  train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
  test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])
  
  # Hyperparameter Tuning and Cross-Validation
  tune_result <- tune(svm, binary.Total.24hr.FFP ~ ., data = train_scaled, 
                      kernel = "radial", 
                      ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))
  best_model <- tune_result$best.model
  
  # Train SVM model with best parameters
  svm_model <- svm(binary.Total.24hr.FFP ~ ., data = train_scaled, type = 'C-classification', 
                   kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)
  
  # Predict on test data
  svm_pred <- predict(svm_model, test_scaled, probability = T)
  
  # Extract the probabilities for the positive class
  svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
  
  # Calculate ROC
  roc_response <- roc(test_scaled$binary.Total.24hr.FFP, svm_pred_prob)
  
  # Store AUC
  auc_values[i] <- roc_response$auc
  valid_iterations <- valid_iterations + 1
}

# Calculate the mean AUC and its standard error
mean_auc <- mean(auc_values[1:valid_iterations])
se_auc <- sd(auc_values[1:valid_iterations]) / sqrt(valid_iterations)

list(mean_auc = mean_auc, se_auc = se_auc)

############################
# BOOTSTRAP Cryo MODEL df5 #

set.seed(123)  # Set a random seed for reproducibility

n_iterations <- 100
auc_values <- numeric(n_iterations)
valid_iterations <- 0

for (i in 1:n_iterations) {
  # Bootstrap resampling
  resampled_index <- sample(nrow(SVM_df_5), size = round(0.8 * nrow(SVM_df_5)), replace = TRUE)
  trainData <- SVM_df_5[resampled_index, ]
  testData <- SVM_df_5[-resampled_index, ]
  
  # Check if all factor variables have at least two levels in the training data
  sufficient_data <- all(sapply(cols_to_factor, function(var) {
    length(levels(factor(trainData[[var]]))) >= 2
  }))
  
  if (!sufficient_data) {
    next  # Skip this iteration if not all factor variables have sufficient data
  }
  
  # Adjust factor levels and scale continuous variables
  for (var in cols_to_factor) {
    trainData[[var]] <- factor(trainData[[var]])
    testData[[var]] <- factor(testData[[var]], levels = levels(trainData[[var]]))
  }
  train_scaled <- trainData
  test_scaled <- testData
  train_scaled[continuous_vars] <- scale(train_scaled[continuous_vars])
  test_scaled[continuous_vars] <- scale(test_scaled[continuous_vars])
  
  # Hyperparameter Tuning and Cross-Validation
  tune_result <- tune(svm, binary.Total.24hr.Cryo ~ ., data = train_scaled, 
                      kernel = "radial", 
                      ranges = list(cost = 10^(-1:2), gamma = 10^(-2:1)))
  best_model <- tune_result$best.model
  
  # Train SVM model with best parameters
  svm_model <- svm(binary.Total.24hr.Cryo ~ ., data = train_scaled, type = 'C-classification', 
                   kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, probability = T)
  
  # Predict on test data
  svm_pred <- predict(svm_model, test_scaled, probability = T)
  
  # Extract the probabilities for the positive class
  svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
  
  # Calculate ROC
  roc_response <- roc(test_scaled$binary.Total.24hr.Cryo, svm_pred_prob)
  
  # Store AUC
  auc_values[i] <- roc_response$auc
  valid_iterations <- valid_iterations + 1
}

# Calculate the mean AUC and its standard error
mean_auc <- mean(auc_values[1:valid_iterations])
se_auc <- sd(auc_values[1:valid_iterations]) / sqrt(valid_iterations)

list(mean_auc = mean_auc, se_auc = se_auc)

####################################
### QUESTION 2 SURVIVAL ANALYSIS ###
####################################

# Convert status variables to binary
data1$ALIVE_30DAYS_YN <- as.factor(data1$ALIVE_30DAYS_YN)
data1$ALIVE_30DAYS_YN <- ifelse(data1$ALIVE_30DAYS_YN == "Y", 1, 0)
data1$ALIVE_90DAYS_YN <- as.factor(data1$ALIVE_90DAYS_YN)
data1$ALIVE_90DAYS_YN <- ifelse(data1$ALIVE_90DAYS_YN == "Y", 1, 0)
df$ALIVE_12MTHS_YN <- as.factor(df$ALIVE_12MTHS_YN)
df$ALIVE_12MTHS_YN <- ifelse(df$ALIVE_12MTHS_YN == "Y", 1, 0)

# Create Survival Object
df$DEATH_DATE <- as.numeric(df$DEATH_DATE)
df$DEATH_DATE <- ifelse(is.na(df$DEATH_DATE), 0, 1)

# Factorize RBC transfusion variable
df$RBC_tfsd <- as.factor(df$RBC_tfsd)

# Survival Curve - Mortality and RBC_tfsd

sf1 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data = df)
print(sf1)
plot(sf1, col = 1:2, ylim = c(0.8, 1), xlab = "Days from Surgery", ylab = "Survival", main = "Survival Probability based on RBC Transfusion")
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df), fun = "cloglog")

# Perform log-rank test
survdiff(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df)

# Cox PH model
coxmod1 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data = df)
summary(coxmod1)

cox.zph(coxmod1)

# Survival Curve 2 - Duration of ICU stay based on RBC transfusion. 
# Can make a composite variable, choosing a cutoff for duration of ICU stay in days for status. 
# ICU stays longer than 4 are typically considered long

df$Long_ICU_Stay <- ifelse(df$Duration.of.ICU.Stay..days. > 7, 1, 0)

sf2 <- survfit(Surv(Duration.of.ICU.Stay..days., Long_ICU_Stay == 0) ~ RBC_tfsd, data = df)
print(sf2)
plot(sf2, col = 1:2, xlim = c(0,20), main = "Duration of ICU Stay based on RBC Transfusion", xlab = "ICU Time in Days", ylab = "Probability of ICU Discharge")
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(Duration.of.ICU.Stay..days., Long_ICU_Stay == 0) ~ RBC_tfsd, data=df), fun = "S")
# Doesnt meet PH assumptions 

# plot a cloglog plot against log(t)
plot(survfit(Surv(Duration.of.ICU.Stay..days., Long_ICU_Stay == 0) ~ RBC_tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod2 <- coxph(Surv(Duration.of.ICU.Stay..days., Long_ICU_Stay == 0) ~ RBC_tfsd, data = df)
summary(coxmod2)

cox.zph(coxmod2)


# Survival Curve 3 -  Mortality and FFP_tfsd
sf3 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data = df)
print(sf3)

plot(sf3, col = 1:2)
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FPP_tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data=df), fun = "cloglog")

# Perform log-rank test
survdiff(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df)

# Cox PH model
coxmod3 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data = df)
summary(coxmod3)

cox.zph(coxmod3)

# Survival curve 4 - Length of Hospital stay
# Create composite variable for long hospital stay
# typical for lung transplant patients is 3-4 weeks, so ill choose 30 days as a cutoff

df$Long_Hospital_Stay <- ifelse(df$HOSPITAL_LOS > 30, 1, 0)

sf4 <- survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data = df)
print(sf4)
plot(sf4, col = 1:2, xlim = c(8,35), main = "Duration of Hospital Stay based on RBC Transfusion", xlab = "Length of Hospital Stay (Days)", ylab = "Probability of Hospital Discharge")
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data=df), fun = "S")
# Doesnt meet PH assumptions 
# plot a cloglog plot against log(t)
plot(survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data=df), fun = "cloglog")

# Perform log-rank test
survdiff(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df)

# Cox PH model
coxmod4 <- coxph(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data = df)
summary(coxmod4)

cox.zph(coxmod4)

# Survival Curve 5 -  Mortality and Plt_tfsd
sf5 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Plt_Tfsd, data = df)
print(sf5)

plot(sf5, col = 1:2)
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 
# bad curve

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Plt_Tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Plt_Tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod5 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Plt_Tfsd, data = df)
summary(coxmod3)

cox.zph(coxmod5)


# Survival Curve 6 -  Mortality and Cryo_tfsd

sf6 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Cryo_tfsd, data = df)
print(sf6)

plot(sf6, col = 1:2)
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 
# bad curve

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Cryo_tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Cryo_tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod6 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Cryo_tfsd, data = df)
summary(coxmod6)

cox.zph(coxmod6)

# Survival Curve 7 - Mortality and any one of FFP, Cryo, Plt
# create composite variable
df$Other_transfusion <- ifelse(df$Plt_Tfsd == 1 | df$Cryo_tfsd == 1 | df$FFP_tfsd == 1, 1, 0)
sf7 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Other_transfusion, data = df)
print(sf7)

plot(sf7, col = 1:2, ylim = c(0.6,1), xlab = "Days from Surgery", ylab = "Survival", main = "Survival Probability based on any Transfusion")
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 
# bad curve

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Other_transfusion, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Other_transfusion, data=df), fun = "cloglog")

# Cox PH model
coxmod7 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ Other_transfusion, data = df)
summary(coxmod7)

cox.zph(coxmod7)
# not significant 