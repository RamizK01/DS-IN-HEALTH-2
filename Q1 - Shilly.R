############## Massive Transfusions ###################
# Load necessary libraries
library(e1071)
library(caret)
cleaned_transfusion_data<- read.csv("cleaned_transfusion_data.csv")

# Check for missing values
if (sum(is.na(cleaned_transfusion_data)) > 0) {
}

# Splitting the data into training and testing sets
set.seed(123)  

#take out some variables 
SVM_df <- cleaned_transfusion_data[, c(-1, -2, -27, -(30:57), -(59:62) )]

# scale numerical variables
SVM_df[,c(3:6,17:24)] <- scale(SVM_df[,c(3:6,17:24)])

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
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma)

# Predict on test data
predictions <- predict(svm_model, test_scaled)

# Evaluate model
confusionMatrix(predictions, test_scaled$Massive.Transfusion)

# CALCULATE AUC

library(pROC)

# Retrain SVM model with probability predictions
svm_model_prob <- svm(Massive.Transfusion ~ ., data = train_scaled, type = 'C-classification', 
                      kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, 
                      probability = TRUE)

# Predict probabilities on test data
prob_predictions <- predict(svm_model_prob, test_scaled, probability = TRUE)

# Extract the probabilities for the positive class
probabilities <- attr(prob_predictions, "probabilities")[,2]

# Calculate the AUC
roc_curve <- roc(test_scaled$Massive.Transfusion, probabilities)
auc_value <- auc(roc_curve)

# Print the AUC value
print(auc_value)



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
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma)

# Predict on test data
predictions <- predict(svm_model, test_scaled)

# Evaluate model
confusionMatrix(predictions, test_scaled$binary.Total.24hr.RBC)

# AUC CALCULATION
library(pROC)

# Extract decision values
decision_values <- attr(prob_predictions, "decision.values")[, 1]

# Calculate ROC and AUC
roc_response <- roc(response = test_scaled$binary.Total.24hr.RBC, predictor = decision_values)
auc(roc_response)


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
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma)

# Predict on test data
predictions <- predict(svm_model, test_scaled)

# Evaluate model
confusionMatrix(predictions, test_scaled$binary.Total.24hr.Plt)

# CALCULATE AUC
library(pROC)

# Retrain SVM model with probability predictions
svm_model_prob <- svm(binary.Total.24hr.Plt ~ ., data = train_scaled, type = 'C-classification', 
                      kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, 
                      probability = TRUE)

# Predict probabilities on test data
prob_predictions <- predict(svm_model_prob, test_scaled, probability = TRUE)

# Extract the probabilities for the positive class
probabilities <- attr(prob_predictions, "probabilities")[,2]

# Calculate the AUC
roc_curve <- roc(test_scaled$binary.Total.24hr.Plt, probabilities)
auc_value <- auc(roc_curve)

# Print the AUC value
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
                 kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma)

# Predict on test data
predictions <- predict(svm_model, test_scaled)

# Evaluate model
confusionMatrix(predictions, test_scaled$binary.Total.24hr.FFP)

# CALCULATE AUC 
library(pROC)

# Retrain SVM model with probability predictions
svm_model_prob <- svm(binary.Total.24hr.FFP ~ ., data = train_scaled, type = 'C-classification', 
                      kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, 
                      probability = TRUE)

# Predict probabilities on test data
prob_predictions <- predict(svm_model_prob, test_scaled, probability = TRUE)

# Extract the probabilities for the positive class
probabilities <- attr(prob_predictions, "probabilities")[,2]

# Calculate the AUC
roc_curve <- roc(test_scaled$binary.Total.24hr.FFP, probabilities)
auc_value <- auc(roc_curve)

# Print the AUC value
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

# auc time!!!
library(pROC)

# Retrain SVM model with probability predictions
svm_model_prob <- svm(binary.Total.24hr.Cryo ~ ., data = train_scaled, type = 'C-classification', 
                      kernel = 'radial', cost = best_model$cost, gamma = best_model$gamma, 
                      probability = TRUE)

# Predict probabilities on test data
prob_predictions <- predict(svm_model_prob, test_scaled, probability = TRUE)

# Extract the probabilities for the positive class
probabilities <- attr(prob_predictions, "probabilities")[,2]

# Calculate the AUC
roc_curve <- roc(test_scaled$binary.Total.24hr.Cryo, probabilities)
auc_value <- auc(roc_curve)

# Print the AUC value
print(auc_value)