# read df
df <- readxl::read_excel("transfusion data.xlsx")

# index columns of interest
df <- df[,c(1,2,3,4:15,25,26,28:41,46,47,49,50,57,60:65,77:92,116,117)]

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


### Dealing with missing data ###

# All these columns missing data should be 0
df[, 42:56] <- replace(df[, 42:56], is.na(df[, 42:56]), 0)
# Now evaluate missing data

colSums(is.na(df))
# Fibrinogen: 97.4% missing, can we impute using domain knowledge?
# Death Date and Survival Time: 83.33% missing, no worries this is systematic
# LAS score is 6.25% missing, imputation?
# Duration of ICU stays: 1 missing, we can just remove pt
# Pre_PTT: 1 missing, we can impute i guess using domain knowledge

### EDA REPORT EXPORTED ###
library(DataExplorer)
create_report(df)
data()

### EXPORTING FINALCLEANED DATA ###
# Assuming your cleaned dataframe is named df


df <- read_csv("cleaned_transfusion_data.csv")
create_report(df)

as.factor(df$Gender)
# Support Vector Machines
# load packages
library(ISLR)
library(e1071)
library(cluster)
library(epiR)
library(pROC)

# Split data into training and test groups
as.factor(df$Massive.Transfusion)

set.seed(123)

train.index <- train.I <- sample(nrow(df), round(nrow(df)/2))

train <- df[train.index, c(28,29,58)]
test <- df[-train.index, c(28,29,58)]

svmfit <- svm(Massive.Transfusion ~ ECLS_ECMO + ECLS_CPB, data = train, kernel = "linear", cost = 1, scale=FALSE, decision.values=T, probability = T)

head(cbind(Actual = test$Massive.Transfusion, Predicted = predictions))
predictions <- predict(svmfit, newdata = test)

confusion_matrix <- table(Actual = test$Massive.Transfusion, Predicted = predictions)

coefficients <- coef(svmfit)

x <- plot(svmfit, df[, 28:29])



dat <- data.frame(x=x, y= train$Massive.Transfusion) # Create data frame



plot(svmfit, iris, Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = mean(iris$Sepal.Width), Sepal.Length = mean(iris$Sepal.Length)))


set.seed(1)
x <- matrix(rnorm(20*2), ncol=2)
y <- c(rep(-1,10), rep(1,10))
x[y==1,] <- x[y==1,] + 1 # Add 1 to separate classes (Class 1 & -1)

plot(x, col=(3-y)) 

dat <- data.frame(x=x, y=as.factor(y)) # Create data frame

# Run the SVM model
svmfit <- svm(y ~ ., data=dat , kernel = "linear", cost = 2, scale=FALSE, decision.values=T, probability = T)

plot(svmfit, dat)





# predictions with decision values and probabilities
# Predictions on training data
pr <- predict(svmfit, dat, decision.values = T, probability = T)

# Confusion Matrix to evaluate the model
t <- table(pr, dat$y)
t

# Calculate accuracy
accuracy <- sum(diag(t))/sum(t)
accuracy

colSums(is.na(df))
