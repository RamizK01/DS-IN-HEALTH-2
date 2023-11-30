#### RANDOM FORESTS: ####

# Loading necessary libraries
library(randomForest)
library(dplyr)
library(readxl)

# Reading the data set
df <- readxl::read_excel("cleaned_transfusion_data.xlsx")

# Data preparation 
# Changing the column names to better suit random forests 
names(df) <- gsub(" ", "_", names(df))
names(df) <- gsub("-", "_", names(df))

# Renaming columns to remove spaces and special characters
names(df) <- gsub("[[:punct:]]", "_", names(df))
names(df) <- gsub(" ", "_", names(df))

# Random Forest for Total 24hr Plt
rf_Plt <- randomForest(`Total_24hr_Plt` ~ ., data = df, na.action = na.omit)

# Random Forest for Total 24hr FFP
rf_FFP <- randomForest(`Total_24hr_FFP` ~ ., data = df, na.action = na.omit)

# Random Forest for Total 24hr Cryo
rf_Cryo <- randomForest(`Total_24hr_Cryo` ~ ., data = df, na.action = na.omit)

# Random Forest for Total 24hr RBC
rf_RBC <- randomForest(`Total_24hr_RBC` ~ ., data = df, na.action = na.omit)


# Creating Random Forest models
rf_Plt <- randomForest(Total_24hr_Plt ~ ., data = df, na.action = na.omit)
rf_FFP <- randomForest(Total_24hr_FFP ~ ., data = df, na.action = na.omit)
rf_Cryo <- randomForest(Total_24hr_Cryo ~ ., data = df, na.action = na.omit)
rf_RBC <- randomForest(Total_24hr_RBC ~ ., data = df, na.action = na.omit)

# Extracting variable importance
importance_Plt <- importance(rf_Plt)
importance_FFP <- importance(rf_FFP)
importance_Cryo <- importance(rf_Cryo)
importance_RBC <- importance(rf_RBC)

# Plotting variable importance
varImpPlot(rf_Plt, main = "Variable Importance for Total 24hr Plt")
varImpPlot(rf_FFP, main = "Variable Importance for Total 24hr FFP")
varImpPlot(rf_Cryo, main = "Variable Importance for Total 24hr Cryo")
varImpPlot(rf_RBC, main = "Variable Importance for Total 24hr RBC")

