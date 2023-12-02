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

df$Survival_Time <- ifelse(is.na(df$Survival_Time), 365, df$Survival_Time )

data1$ALIVE_30DAYS_YN <- as.factor(data1$ALIVE_30DAYS_YN)
data1$ALIVE_30DAYS_YN <- ifelse(data1$ALIVE_30DAYS_YN == "Y", 1, 0)
data1$ALIVE_90DAYS_YN <- as.factor(data1$ALIVE_90DAYS_YN)
data1$ALIVE_90DAYS_YN <- ifelse(data1$ALIVE_90DAYS_YN == "Y", 1, 0)
df$ALIVE_12MTHS_YN <- as.factor(df$ALIVE_12MTHS_YN)
df$ALIVE_12MTHS_YN <- ifelse(df$ALIVE_12MTHS_YN == "Y", 1, 0)

# Create Survival Object
df$DEATH_DATE <- as.numeric(df$DEATH_DATE)
df$DEATH_DATE <- ifelse(is.na(df$DEATH_DATE), 0, 1)


df$RBC_tfsd <- as.factor(df$RBC_tfsd)

# Survival Curve - Mortality and RBC_tfsd

sf1 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data = df)
print(sf1)

plot(sf1, col = 1:2)
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod1 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ RBC_tfsd, data = df)
summary(coxmod1)

cox.zph(coxmod1)

# Survival Curve 2 - Duration of ICU stay based on RBC transfusion. 
# Can make a composite variable, choosing a cutoff for duration of ICU stay in days for status. 
# ICU stays longer than 4 are typically considered long

df$Long_ICU_Stay <- ifelse(df$Duration.of.ICU.Stay..days. > 4, 1, 0)

sf2 <- survfit(Surv(Duration.of.ICU.Stay..days., Long_ICU_Stay == 0) ~ RBC_tfsd, data = df)
print(sf2)
plot(sf2, col = 1:2, xlim = c(0,20))

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

# Survival curve 4 - Length of Hospital stay
# Create composite variable for long hospital stay
# typical for lung transplant patients is 3-4 weeks, so ill choose 30 days as a cutoff

df$Long_Hospital_Stay <- ifelse(df$HOSPITAL_LOS > 30, 1, 0)

sf4 <- survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data = df)
print(sf4)
plot(sf4, col = 1:2, xlim = c(0,35))

legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data=df), fun = "S")
# Doesnt meet PH assumptions 
# plot a cloglog plot against log(t)
plot(survfit(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod4 <- coxph(Surv(HOSPITAL_LOS, Long_Hospital_Stay == 0) ~ RBC_tfsd, data = df)
summary(coxmod4)

cox.zph(coxmod4)

# Survival Curve 3 -  Mortality and FFP_tfsd
sf3 <- survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data = df)
print(sf3)

plot(sf3, col = 1:2)
legend("bottomright", legend = c("No Transfusion", "Transfusion"), lty = 1, col = 1:2) 

# typical stratified KM curve
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FPP_tfsd, data=df), fun = "S")

# plot a cloglog plot against log(t)
plot(survfit(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data=df), fun = "cloglog")

# Cox PH model
coxmod3 <- coxph(Surv(Survival_Time, ALIVE_12MTHS_YN == 0) ~ FFP_tfsd, data = df)
summary(coxmod3)

cox.zph(coxmod3)

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

plot(sf7, col = 1:2)
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