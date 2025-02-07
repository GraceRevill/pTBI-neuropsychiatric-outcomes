
################################################################
#
# 
# Supplementary analysis (comparing SMD for 1:1 vs 2:1 matching) of psychiatric symptoms, 
# diagnoses and service use using MatchIt for propensity score matching in the study:
#
# Only Anxiety Remains Reliably Associated with Paediatric 
# Mild TBI at Two Years Follow-up after Adjusting for Pre-Existing Mental Health Problems 
#
# Revill et al
#
#
################################################################

#
# This script loads ABCD 4.0 data from the data_dir directory
#
# It writes imputed datasets and publication ready tables to the output_dir directory
#
library(dplyr)
library(lme4)
library(jtools)
library(sjPlot)
library(MatchIt)
library(visdat)
library(gtsummary)
library(cobalt)

# Remove everything from memory if needed to start afresh
# rm(list = ls())

# Set data and output directories
data_dir = "S:/ABCDPsychosis/ABCD4.0/"
output_dir = "S:/ABCDPsychosis/Grace/Output/"

setwd(data_dir)


###########################################################
#
# Load CBCL and KSADS data 
#
###########################################################

# Load imputed dataset for baseline 
imputed_data_baseline <- read.csv(paste(output_dir, "imputed_baseline_KSAD_CBCL.csv", sep = ""))

# Load imputed dataset for 2-year follow-up 
imputed_data_t2 <- read.csv(paste(output_dir, "imputed_t2_KSAD_CBCL.csv", sep = ""))

###########################################################
#
# Load psychiatric service use data 
#
###########################################################

# Load psyc service imputed dataset for baseline 
imputed_psyc_data_baseline <- read.csv(paste(output_dir, "imputed_baseline_psyc_services.csv", sep = ""))

# Load psyc service imputed dataset for 2-year follow-up 
imputed_psyc_data_t2 <- read.csv(paste(output_dir, "imputed_t2_psyc_services.csv", sep = ""))

###########################################################
#
# Rename and recode CBCL and KSADS variables
#
###########################################################

# Rename CBCL and KSADS imputed baseline variables 
baseline_vars <- imputed_data_baseline %>%
  select(subject_key, baseline_cbcl_tot = cbcl_tot_problems,
         baseline_anxiety_disorders = anxiety_disorders,
         baseline_behavioural_disorders = behavioural_disorders, 
         baseline_cbcl_int = cbcl_int_problems,
         baseline_cbcl_ext = cbcl_ext_problems,
         baseline_cbcl_anx = cbcl_anx_score,
         baseline_cbcl_dep = cbcl_dep_score,
         baseline_cbcl_adhd = cbcl_adhd_score,
         baseline_cbcl_odd = cbcl_odd_score,
         baseline_cbcl_conduct = cbcl_conduct_score)


# Add baseline variables into follow-up data
imputed_data_t2 <- imputed_data_t2 %>%
  left_join(baseline_vars, by = "subject_key")


# Dummy variable generation
ethnicity_dummy_variables <- model.matrix(~ race_ethnicity - 1, data = imputed_data_t2)
imputed_data_t2 <- cbind(imputed_data_t2, ethnicity_dummy_variables)

# Recode variables to 0 and 1
imputed_data_t2 <- imputed_data_t2 %>%
  mutate(male = ifelse(sex == "Male", 1, 0))

imputed_data_t2 <- imputed_data_t2 %>%
  mutate(advlife_cat = ifelse(advlife_cat == "Yes", 1, 0))

# Make sure no missing values present 
imputed_data_t2 <- na.omit(imputed_data_t2)

# Recode treatment variable (TBI vs no injury) to 0 and 1 (doesn't like labels in treatment)
imputed_data_t2$tbi_present <- paste(imputed_data_t2$injury_type)

imputed_data_t2$tbi_present[imputed_data_t2$tbi_present == 'None'] <- 0
imputed_data_t2$tbi_present[imputed_data_t2$tbi_present == 'TBI'] <- 1
imputed_data_t2$tbi_present[imputed_data_t2$tbi_present == 'Ortho'] <- NA

# Remove NA values and create new data frame with TBI and no injury cases 
tbi_data <- filter(imputed_data_t2, tbi_present != "NA")

tbi_data$tbi_present <- as.numeric(tbi_data$tbi_present)

# Recode treatment (OI vs no injury) variable to 0 and 1 (doesn't like labels in treatment)
imputed_data_t2$oi_present <- paste(imputed_data_t2$injury_type)

imputed_data_t2$oi_present[imputed_data_t2$oi_present == 'Ortho'] <- 1
imputed_data_t2$oi_present[imputed_data_t2$oi_present == 'TBI'] <- NA
imputed_data_t2$oi_present[imputed_data_t2$oi_present == 'None'] <- 0

# Remove NA values and create new data frame with OI and no injury cases 
oi_data <- filter(imputed_data_t2, oi_present != "NA")

oi_data$oi_present <- as.numeric(oi_data$oi_present)


# Supplementary material: Recode treatment variable (TBI vs OI) to 0 and 1 (doesn't like labels in treatment)
imputed_data_t2$tbi_oi <- paste(imputed_data_t2$injury_type)

imputed_data_t2$tbi_oi[imputed_data_t2$tbi_oi == 'None'] <- NA
imputed_data_t2$tbi_oi[imputed_data_t2$tbi_oi == 'TBI'] <- 1
imputed_data_t2$tbi_oi[imputed_data_t2$tbi_oi == 'Ortho'] <- 0

# Remove NA values and create new data frame with TBI and no injury cases 
tbi_oi_data <- filter(imputed_data_t2, tbi_oi != "NA")

tbi_oi_data$tbi_oi <- as.numeric(tbi_oi_data$tbi_oi)

# Visualise missing data
vis_miss(tbi_data)
vis_miss(oi_data)
vis_miss(tbi_oi_data)

# Set seed to make results reproducible
set.seed(7550822)

###########################################################
#
# Rename and recode psyc service variables
#
###########################################################

# Rename psyc service imputed baseline variables 
baseline_psyc_vars <- imputed_psyc_data_baseline %>%
  select(subject_key, baseline_any_support = any_support,
         baseline_outpatient = outpatient,
         baseline_inpatient = inpatient, 
         baseline_psychotherapy = psychotherapy,
         baseline_medication = medication,
         baseline_psychotherapy_medication = psychotherapy_medication)


# Add baseline variables into follow-up data
imputed_psyc_data_t2 <- imputed_psyc_data_t2 %>%
  left_join(baseline_psyc_vars, by = "subject_key")


# Dummy variable generation
ethnicity_dummy_variables <- model.matrix(~ race_ethnicity - 1, data = imputed_psyc_data_t2)
imputed_psyc_data_t2 <- cbind(imputed_psyc_data_t2, ethnicity_dummy_variables)

# Recode variables to 0 and 1
imputed_psyc_data_t2 <- imputed_psyc_data_t2 %>%
  mutate(male = ifelse(sex == "Male", 1, 0))

imputed_psyc_data_t2 <- imputed_psyc_data_t2 %>%
  mutate(advlife_cat = ifelse(advlife_cat == "Yes", 1, 0))

# Make sure no missing values present 
imputed_psyc_data_t2 <- na.omit(imputed_psyc_data_t2)

# Recode treatment variable (TBI vs no injury) to 0 and 1 (doesn't like labels in treatment)
imputed_psyc_data_t2$tbi_present <- paste(imputed_psyc_data_t2$injury_type)

imputed_psyc_data_t2$tbi_present[imputed_psyc_data_t2$tbi_present == 'None'] <- 0
imputed_psyc_data_t2$tbi_present[imputed_psyc_data_t2$tbi_present == 'TBI'] <- 1
imputed_psyc_data_t2$tbi_present[imputed_psyc_data_t2$tbi_present == 'Ortho'] <- NA

# Remove NA values and create new data frame with TBI and no injury cases 
tbi_psyc_data <- filter(imputed_psyc_data_t2, tbi_present != "NA")

tbi_psyc_data$tbi_present <- as.numeric(tbi_psyc_data$tbi_present)

# Recode treatment (Ortho vs no injury) variable to 0 and 1 (doesn't like labels in treatment)
imputed_psyc_data_t2$oi_present <- paste(imputed_psyc_data_t2$injury_type)

imputed_psyc_data_t2$oi_present[imputed_psyc_data_t2$oi_present == 'Ortho'] <- 1
imputed_psyc_data_t2$oi_present[imputed_psyc_data_t2$oi_present == 'TBI'] <- NA
imputed_psyc_data_t2$oi_present[imputed_psyc_data_t2$oi_present == 'None'] <- 0

# Remove NA values and create new data frame with OI and no injury cases 
oi_psyc_data <- filter(imputed_psyc_data_t2, oi_present != "NA")

oi_psyc_data$oi_present <- as.numeric(oi_psyc_data$oi_present)

# Supplementary material: Recode treatment variable (TBI vs OI) to 0 and 1 (doesn't like labels in treatment)
imputed_psyc_data_t2$tbi_oi <- paste(imputed_psyc_data_t2$injury_type)

imputed_psyc_data_t2$tbi_oi[imputed_psyc_data_t2$tbi_oi == 'None'] <- NA
imputed_psyc_data_t2$tbi_oi[imputed_psyc_data_t2$tbi_oi == 'TBI'] <- 1
imputed_psyc_data_t2$tbi_oi[imputed_psyc_data_t2$tbi_oi == 'Ortho'] <- 0

# Remove NA values and create new data frame with TBI and no injury cases 
tbi_oi_psyc_data <- filter(imputed_psyc_data_t2, tbi_oi != "NA")

tbi_oi_psyc_data$tbi_oi <- as.numeric(tbi_oi_psyc_data$tbi_oi)

# Visualise missing data
vis_miss(tbi_psyc_data)
vis_miss(oi_psyc_data)
vis_miss(tbi_oi_psyc_data)

# Set seed to make results reproducible
set.seed(7550822)


###########################################################
#
# CBCL and KSADS propensity score matching using MatchIt
#
###########################################################

# Make sure TBI vs no injury variables are factors 
as.factor(tbi_data$baseline_anxiety_disorders)
as.factor(tbi_data$tbi_present)
as.factor(tbi_data$anxiety_disorders)

# Make sure TBI vs Ortho variables are factors 
as.factor(oi_data$baseline_anxiety_disorders)
as.factor(oi_data$oi_present)
as.factor(oi_data$anxiety_disorders)

## KSADS Anxiety (TBI vs no injury) 

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


#  2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

##  KSADS Anxiety (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## KSADS Anxiety (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## KSADS Behavioural Disorders (TBI vs no injury)

# 1:1 matching  
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## KSADS Behavioural Disorders (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = oi_data,
                 method = "optimal", ratio = 2,  distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## KSADS Behavioural Disorders (TBI vs OI) 

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2,distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## CBCL Total Problems (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL Total Problems (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL Total Problems (TBI vs OI) 

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


### CBCL subscales 

## CBCL internalising symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = tbi_data,
                 method = "optimal", ratio = 2,  distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL internalising symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## CBCL internalising symptoms (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL externalising symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL externalising symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL externalising symptoms (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL anxiety symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

#2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL anxiety symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL anxiety symptoms (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL depression symptoms (TBI vs no injury)

# 1:1 matching
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL depression symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL depression symptoms (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## CBCL ADHD symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL ADHD symptoms (OI vs no injury)

# 1:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL ADHD symptoms (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL ODD symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## CBCL ODD symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL ODD symptoms (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL conduct symptoms (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = tbi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## CBCL conduct symptoms (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## CBCL conduct symptoms (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = tbi_oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = tbi_oi_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")



###########################################################
#
# Psychiatric service use propensity score matching using MatchIt
#
###########################################################

## Any mental health support (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Any mental health support (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Any mental health support (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Outpatient support (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Outpatient support (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Outpatient support (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Inpatient support (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

## Inpatient support (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Inpatient support (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Psychotherapy (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")


## Psychotherapy (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Psychotherapy (TBI vs OI)

# 1:1 matching
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Medication (TBI vs no injury)

# 1:1 matching
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Medication (OI vs no injury)

# 1:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Medication (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)


# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Psychotherapy and medication (TBI vs no injury)

# 1:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = tbi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Psychotherapy and medication (OI vs no injury)

# 1:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

## Psychotherapy and medication (TBI vs OI)

# 1:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = tbi_oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

# 2:1 matching 
model <- matchit(tbi_oi ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = tbi_oi_psyc_data,
                 method = "optimal", ratio = 2, distance = "glm")

summary(model, standardize = TRUE)

# Check SMD before and after matching 
balance_table <- bal.tab(model, stars = TRUE)
print(balance_table)
love.plot(model, stats = "mean.diffs")

