# read df
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

# Pre_PTT: 1 missing, mean impute
mean_pre_ptt <- mean(df$Pre_PTT, na.rm = TRUE)
df$Pre_PTT <- ifelse(is.na(df$Pre_PTT), mean_pre_ptt, df$Pre_PTT)

### Creating 24hr blood product variables ###
df <- df %>% mutate('Total 24hr Plt' = Intra_Platelets + `Plt 0-24hrs`) %>%
             mutate('Total 24hr FFP' = `Intra_Fresh Frozen Plasma` + `FFP 0-24hrs`) %>%
             mutate('Total 24hr Cryo' = `Intra_Cryoprecipitate` + `Cryo 0-24hrs`)

### EXPORTING FINAL  CLEANED DATA ###
write.csv(df, file = "cleaned_transfusion_data.csv", row.names = FALSE)

### EDA REPORT EXPORTED ###
library(DataExplorer)
create_report(df)
data()

# QUESTION 2: 

data1$ALIVE_30DAYS_YN <- as.factor(data1$ALIVE_30DAYS_YN)
data1$ALIVE_30DAYS_YN <- ifelse(data1$ALIVE_30DAYS_YN == "Y", 1, 0)
data1$ALIVE_90DAYS_YN <- as.factor(data1$ALIVE_90DAYS_YN)
data1$ALIVE_90DAYS_YN <- ifelse(data1$ALIVE_90DAYS_YN == "Y", 1, 0)
df$ALIVE_12MTHS_YN <- as.factor(df$ALIVE_12MTHS_YN)
df$ALIVE_12MTHS_YN <- ifelse(df$ALIVE_12MTHS_YN == "Y", 1, 0)

# Create Survival Object
df$DEATH_DATE <- as.numeric(df$DEATH_DATE)
df$DEATH_DATE <- ifelse(is.na(df$DEATH_DATE), 0, 1)


surv1 <- Surv(df$Survival_Time, !is.na(df$DEATH_DATE))

# Example with only ECMO and CPB as predictors
sf1 <- survfit(surv1 ~ ECLS_ECMO + ECLS_CPB, data = df)
print(sf1)

plot(sf1, col = 1:2)
legend("bottomright", legend = c("ECMO", "CPB"), lty = 1, col = 1:2) 


# typical stratified KM curve
plot(survfit(surv1 ~ ECLS_ECMO + ECLS_CPB, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(surv1 ~ ECLS_ECMO + ECLS_CPB, data=df), fun = "cloglog")

survdiff(Surv(time, status==1) ~ sex,data=melanoma)

# Cox PH model
coxmod1 <- coxph(surv1 ~ ECLS_ECMO + ECLS_CPB, data = df)
summary(coxmod1)

cox.zph(coxmod1)



# Have to find probability of survival for each patient at each time point (30 days, 90 days, and 1 year)
# 160 censored observations. Known that 32 patients died. 
# Convert variables into factors



# Use the characteristics of the patients found in the first question to build cox model

data1$`OR Date` <- as.numeric(data1$`OR Date`)

sf1 <- survfit(Surv(`OR Date`, ALIVE_30DAYS_YN == 0) ~ 1, data=data1)
print(sf)

sf2 <- survfit(Surv(`OR Date`, ALIVE_90DAYS_YN == 0) ~ 1, data=data1)

sf3 <- survfit(Surv(`OR Date`, ALIVE_12MTHS_YN == 0) ~ 1, data=data1)

plot(sf2,xscale = 365.25, xlab = "years from surgery", ylab="Survival", col=1:2) 
legend("bottomright", legend = c("female", "male"), lty = 1, col = 1:2) 

# Then, do coxph for each three models as well. 

# Visualize using kaplan meier curves
