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



################################################################
### LASSO MODEL FOR PREDICTING IF FFP TRANSFUSION WILL OCCUR ###
################################################################



#################################################################
### LASSO MODEL FOR PREDICTING IF CRYO TRANSFUSION WILL OCCUR ###
#################################################################



#################################################################
### MULTIVARIATE LASSO MODEL TO PREDICT AMOUNT OF TRANSFUSION ###
#################################################################

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
                      type.measure = "auc")
lambda_optimal <- cv.lasso$lambda.min
lasso.model <- glmnet(x.train, y.train, family = "mgaussian", alpha = 1,
                      lambda = lambda_optimal)

coefficients(lasso.model)

# PREDICT AND CALCULATE MSE
x.test <- model.matrix(Total.24hr.RBC + Total.24hr.Plt + Total.24hr.FFP + Total.24hr.Cryo ~., df_mvar)[-train.set, -1]
pred.lasso <- as.numeric(predict(lasso.model, newx = x.test, s = cv.lasso$lambda.min, type = "response"))

# Actual values for the test set
y.test <- model.matrix(~ df_mvar$Total.24hr.RBC + df_mvar$Total.24hr.FFP +
                         df_mvar$Total.24hr.Plt + df_mvar$Total.24hr.Cryo)[-train.set, -1]

# Calculate MSE
mse <- mean((pred.lasso - y.test)^2)


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

