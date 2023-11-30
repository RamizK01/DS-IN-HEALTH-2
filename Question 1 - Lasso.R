##########################################
### INITIAL DATA CLEANING/VAR CREATION ###
##########################################

df <- read.csv('cleaned_transfusion_data.csv')

library(glmnet)
set.seed(123)

# create df_las 
df <- df %>% mutate("RBC_tfsd" = case_when(`Total.24hr.RBC` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("PLT_tfsd" = case_when(`Total.24hr.Plt` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("FFP_tfsd" = case_when(`Total.24hr.FFP` == 0 ~ 0, TRUE ~ 1)) %>%
             mutate("Cryo_tfsd" = case_when(`Total.24hr.Cryo` == 0 ~ 0, TRUE ~ 1))
  
### THIS DF WILL ONLY WORK FOR RBC, NEED TO ADD NEW COLUMNS 64:66 FOR PLT, FFP, CRYO
df_las <- df[,c(3:26,28,29,63)]
df_las[c(1, 2, 7:16, 25:26)] <- lapply(df_las[c(1, 2, 7:16, 25:26)], factor)


################################################################
### LASSO MODEL FOR PREDICTING IF RBC TRANSFUSION WILL OCCUR ###
################################################################
library(pROC)

train.set <- sample(nrow(df_las),round(nrow(df_las)*.7))

# train predictor 
x.train <- model.matrix(RBC_tfsd ~.,df_las)[train.set,-1]
# train response
y.train <- df_las$RBC_tfsd[train.set]
# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial",
                      type.measure = "auc")
lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "binomial", alpha = 1,
                      lambda = lambda_optimal)

# predict
pred.lasso <- as.numeric(predict(lasso.model, 
                                 newx = model.matrix(RBC_tfsd ~.,df_las)
                                 [-train.set,-1], s=cv.lasso$lambda.min,
                                 type = "response"))

# create ROC curve
myroc <- roc(RBC_tfsd ~ pred.lasso, data=df_las[-train.set,])
auc.lasso <- myroc$auc
plot(myroc)

### WE SHOULD REALLY BOOTSTRAP THE LASSO ###
library(SparseLearner)
bstrap.lasso.model <- BRLasso(x, y, B = 5, Boots = 100, kfold = 10, seed = 0123)


################################################################
### LASSO MODEL FOR PREDICTING IF PLT TRANSFUSION WILL OCCUR ###
################################################################



################################################################
### LASSO MODEL FOR PREDICTING IF FFP TRANSFUSION WILL OCCUR ###
################################################################



#################################################################
### LASSO MODEL FOR PREDICTING IF CRYO TRANSFUSION WILL OCCUR ###
#################################################################



#################################################################
### MULTIVARIATE LASSO MODEL TO PREDICT AMOUNT OF TRANSFUSION ###
#################################################################

# create df_mlas
df_mlas <- df[,c(3:26,28,29,57,60:62)]
df_mlas[c(1, 2, 7:16, 25:26)] <- lapply(df_mlas[c(1, 2, 7:16, 25:26)], factor)

train.set <- sample(nrow(df_mlas),round(nrow(df_mlas)*.7))

# train predictor 
x.train <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~.,df_mlas)[train.set,-1]
# train response
y.train <- model.matrix(~ df_mlas$Total.24hr.RBC + df_mlas$Total.24hr.FFP +
                        df_mlas$Total.24hr.Plt + df_mlas$Total.24hr.Cryo)[train.set, -1]

# Perform cross-validation for lambda selection
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "mgaussian",
                      type.measure = "auc")
lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "mgaussian", alpha = 1,
                      lambda = lambda_optimal)

coefficients(lasso.model)[1]

# PREDICT AND CALCULATE MSE
x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mlas)[-train.set, -1]
pred.lasso <- as.numeric(predict(lasso.model, newx = x.test, s = cv.lasso$lambda.min, type = "response"))

# Actual values for the test set
y.test <- model.matrix(~ df_mlas$Total.24hr.RBC + df_mlas$Total.24hr.FFP +
                         df_mlas$Total.24hr.Plt + df_mlas$Total.24hr.Cryo)[-train.set, -1]

# Calculate MSE
mse <- mean((pred.lasso - y.test)^2)


