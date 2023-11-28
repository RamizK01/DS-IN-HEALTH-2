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
library(mice)
df_PMM <- rename(df, "LAS_score" = "LAS score")

df_PMM <- df[, c(3:25)]

df_PMM <- rename(df_PMM, "alpha1_Antitrypsin_Deficiency" = "alpha1-Antitrypsin Deficiency")
df_PMM <- rename(df_PMM, "Cystic_Fibrosis" = "Cystic Fibrosis")
df_PMM <- rename(df_PMM, "Idiopathic_Pulmonary_Hypertension" = "Idiopathic Pulmonary Hypertension")
df_PMM <- rename(df_PMM, "Interstitial_Lung_Disease" = "Interstitial Lung Disease")
df_PMM <- rename(df_PMM, "First_Lung_Transplant" = "First Lung Transplant")
df_PMM <- rename(df_PMM, "Redo_Lung_Transplant" = "Redo Lung Transplant")
df_PMM <- rename(df_PMM, "ExVIVO_Lung_Perfusion" = "ExVIVO Lung Perfusion")
df_PMM <- rename(df_PMM, "Preoperative_ECLS" = "Preoperative ECLS")
df_PMM <- rename(df_PMM, "LAS_score" = "LAS score")


# Setting up the method for PMM specifically for LAS_score
method <- make.method(df_PMM)
method["LAS_score"] <- "pmm"

# Performing the imputation
imp <- mice(df_PMM, method = method, m = 20, seed = 123, print = FALSE)

fit <- with(imp, lm(LAS_score ~ Type + Gender + Height + Weight + Age + BMI + COPD + alpha1_Antitrypsin_Deficiency
                    + Cystic_Fibrosis + Idiopathic_Pulmonary_Hypertension + Interstitial_Lung_Disease + Pulm_Other 
                    + First_Lung_Transplant + Redo_Lung_Transplant + ExVIVO_Lung_Perfusion + Preoperative_ECLS + Pre_Hb
                    + Pre_Hct + Pre_Platelets + Pre_PT + Pre_INR + Pre_PTT))
summary(imp)
summary(pool(fit))
df_complete <- complete(imp, 1)

new_LAS_scores <- df_complete[,17]

# Insert new LAS scores
df$`LAS score` <- new_LAS_scores
# Duration of ICU stays: 1 missing, we can just remove pt

# Pre_PTT: 1 missing, we can impute i guess using domain knowledge
mean_pre_ptt <- mean(df$Pre_PTT, na.rm = TRUE)
df$Pre_PTT <- ifelse(is.na(df$Pre_PTT), mean_pre_ptt, df$Pre_PTT)

### Creating 24hr blood product variables ###
df <- df %>% mutate('Total 24hr Plt' = Intra_Platelets + `Plt 0-24hrs`) %>%
             mutate('Total 24hr FFP' = `Intra_Fresh Frozen Plasma` + `FFP 0-24hrs`) %>%
             mutate('Total 24hr Cryo' = `Intra_Cryoprecipitate` + `Cryo 0-24hrs`)

### EDA REPORT EXPORTED ###
library(DataExplorer)
create_report(df)
data()

### EXPORTING FINAL  CLEANED DATA ###
write.csv(df, file = "cleaned_transfusion_data.csv", row.names = FALSE)



