# QUESTION 1: 

# Decision Trees: 

# Installing packages and loading into library: 
install.packages("rpart")
install.packages("rpart.plot")
install.packages("tree")
install.packages("dplyr")
library(rpart)
library(rpart.plot)
library(tree)
library(dplyr)


# Converting the outcome variables into factors: 
df <- df %>%
  mutate(`Total 24hr Plt` = as.factor(Intra_Platelets + `Plt 0-24hrs`),
         `Total 24hr FFP` = as.factor(`Intra_Fresh Frozen Plasma` + `FFP 0-24hrs`),
         `Total 24hr Cryo` = as.factor(`Intra_Cryoprecipitate` + `Cryo 0-24hrs`))
as.factor(df$`Total 24hr RBC`)

df <- df %>%
  mutate(`Total 24hr Plt` = as.factor(Intra_Platelets + `Plt 0-24hrs`),
         `Total 24hr FFP` = as.factor(`Intra_Fresh Frozen Plasma` + `FFP 0-24hrs`),
         `Total 24hr Cryo` = as.factor(`Intra_Cryoprecipitate` + `Cryo 0-24hrs`),
         `Total 24hr RBC` = as.factor(`Total 24hr RBC`))


# Set a seed for reproducibility
set.seed(123)

# Sample for training data
train.1 <- sample(nrow(df), round(nrow(df)/2))

# Tree 1: for Total 24hr RBC
tree_RBC <- tree(`Total 24hr RBC` ~ ., data = df, subset = train.1)
plot(tree_RBC)
text(tree_RBC, pretty = 0, cex = 0.5)

# Tree 2: for Total 24hr Plt
tree_Plt <- tree(`Total 24hr Plt` ~ ., data = df, subset = train.1)
plot(tree_Plt)
text(tree_Plt, pretty = 0, cex = 0.5)

# Tree 3: for Total 24hr FFP
tree_FFP <- tree(`Total 24hr FFP` ~ ., data = df, subset = train.1)
plot(tree_FFP)
text(tree_FFP, pretty = 0, cex = 0.5)

# Tree 4: for Total 24hr Cryo
tree_Cryo <- tree(`Total 24hr Cryo` ~ ., data = df, subset = train.1)
plot(tree_Cryo)
text(tree_Cryo, pretty = 0, cex = 0.5)

## OPTIONAL: 
# Cross-validation for pruning: 
cv_Plt <- cv.tree(tree_Plt, FUN = prune.tree)
best_size_Plt <- cv_Plt$size[which.min(cv_Plt$dev)]
pruned_Plt <- prune.tree(tree_Plt, best = best_size_Plt)

# Plot the pruned tree for Plt
plot(pruned_Plt)
text(pruned_Plt, pretty = 0, cex = 0.9)
