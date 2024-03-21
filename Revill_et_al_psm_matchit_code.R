
################################################################
#
# 
# Analysis of psychiatric symptoms, diagnoses and service use using MatchIt 
# for propensity score matching in the study:
#
# Mild TBI but not Orthopaedic Injury Associated with Subsequent
# Psychiatric Problems and Increased Mental Health Service Use in
# 9-12-year-old Children
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




# Create column merging anxiety at baseline and 2 year follow-up 
imputed_data_t2$anxiety_both_years <- paste(imputed_data_t2$baseline_anxiety_disorders, imputed_data_t2$anxiety_disorders)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'No No'] <- 0
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'No NA'] <- 0
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'NA No'] <- 0
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'Yes NA'] <- 0
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'NA Yes'] <- 1
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'Yes No'] <- 0
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'No Yes'] <- 1
imputed_data_t2$anxiety_both_years[imputed_data_t2$anxiety_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_data_t2$anxiety_both_years <- factor(imputed_data_t2$anxiety_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ anxiety_both_years , data = imputed_data_t2)



# Check if variabe is coded correctly 
xtabs(~ anxiety_both_years + injury_type , data = imputed_data_t2)


# Check if variabe is coded correctly 
xtabs(~ anxiety_disorders + injury_type , data = imputed_data_t2)





# Create column merging behavioural disorders at baseline and 2 year follow-up 
imputed_data_t2$behavioural_both_years <- paste(imputed_data_t2$baseline_behavioural_disorders, imputed_data_t2$behavioural_disorders)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'No No'] <- 0
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'No NA'] <- 0
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'NA No'] <- 0
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'Yes NA'] <- 0
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'NA Yes'] <- 1
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'Yes No'] <- 0
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'No Yes'] <- 1
imputed_data_t2$behavioural_both_years[imputed_data_t2$behavioural_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_data_t2$behavioural_both_years <- factor(imputed_data_t2$behavioural_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ behavioural_both_years + injury_type, data = imputed_data_t2)


# Any present anxiety disorder at 2-year follow-up 
#

# Change reference categories 
imputed_data_t2$anxiety_both_years <- relevel(imputed_data_t2$anxiety_both_years, ref = "No")
imputed_data_t2$injury_type <- relevel(imputed_data_t2$injury_type, ref = "None")

# Run unadjusted model 
anxiety_both_years_model_unadjusted_t2 <- glmer(anxiety_both_years ~ injury_type +
                                                 (1|site_id_l),
                                               data = imputed_data_t2,
                                               family = binomial(link = "logit"),
                                               control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
anxiety_both_years_model_adjusted1_t2 <- glmer(anxiety_both_years ~ injury_type +
                                                interview_age +
                                                sex +
                                                race_ethnicity +
                                                household_income +
                                                (1|site_id_l),
                                              data = imputed_data_t2,
                                              family = binomial(link = "logit"),
                                              control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
anxiety_both_years_model_adjusted2_t2 <- glmer(anxiety_both_years ~ injury_type +
                                                interview_age +
                                                sex +
                                                race_ethnicity +
                                                household_income +
                                                advlife_cat +
                                                safety +
                                                family_conflict_scale +
                                                parent_mh_score +
                                                (1|site_id_l),
                                              data = imputed_data_t2,
                                              family = binomial(link = "logit"),
                                              control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)





# Any present behavioural disorder at 2-year follow-up
#

# Change reference categories 
imputed_data_t2$behavioural_both_years <- relevel(imputed_data_t2$behavioural_both_years, ref = "No")

# Run unadjusted model 
behavioural_both_years_model_unadjusted_t2 <- glmer(behavioural_both_years ~ injury_type +
                                                     (1|site_id_l),
                                                   data = imputed_data_t2, family = binomial(link = "logit"),
                                                   control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
behavioural_both_years_model_adjusted1_t2 <- glmer(behavioural_both_years ~ injury_type +
                                                    interview_age +
                                                    sex +
                                                    race_ethnicity +
                                                    household_income +
                                                    (1|site_id_l),
                                                  data = imputed_data_t2,
                                                  family = binomial(link = "logit"),
                                                  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
behavioural_both_years_model_adjusted2_t2 <- glmer(behavioural_both_years ~ injury_type +
                                                    interview_age +
                                                    sex +
                                                    race_ethnicity +
                                                    household_income +
                                                    advlife_cat +
                                                    safety +
                                                    family_conflict_scale +
                                                    parent_mh_score + 
                                                    (1|site_id_l),
                                                  data = imputed_data_t2,
                                                  family = binomial(link = "logit"),
                                                  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)





# Anxiety disorders at 2-year follow-up 
#

# Create unadjusted model table
anxiety_unadjusted_table_t2 <- tbl_regression(anxiety_both_years_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Any present anxiety disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
anxiety_adjusted_table_model1_t2 <- tbl_regression(anxiety_both_years_model_adjusted1_t2, exponentiate = TRUE, 
                                                   pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                   include = injury_type, 
                                                   label = injury_type ~ "Any present anxiety disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
anxiety_adjusted_table_model2_t2 <- tbl_regression(anxiety_both_years_model_adjusted2_t2, exponentiate = TRUE, 
                                                   pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                   include = injury_type, 
                                                   label = injury_type ~ "Any present anxiety disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_anxiety_t2 <- tbl_merge(tbls = list(anxiety_unadjusted_table_t2, anxiety_adjusted_table_model1_t2, anxiety_adjusted_table_model2_t2),
                                    tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Any behavioural disorder at 2-year follow-up 
#

# Create unadjusted model table
behavioural_unadjusted_table_t2 <- tbl_regression(behavioural_both_years_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Any present behavioural disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
behavioural_adjusted_table_model1_t2 <- tbl_regression(behavioural_both_years_model_adjusted1_t2, exponentiate = TRUE, 
                                                       pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                       include = injury_type, 
                                                       label = injury_type ~ "Any present behavioural disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
behavioural_adjusted_table_model2_t2 <- tbl_regression(behavioural_both_years_model_adjusted2_t2, exponentiate = TRUE, 
                                                       pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                       include = injury_type, 
                                                       label = injury_type ~ "Any present behavioural disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_behavioural_t2 <- tbl_merge(tbls = list(behavioural_unadjusted_table_t2, behavioural_adjusted_table_model1_t2, behavioural_adjusted_table_model2_t2),
                                        tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))




# Create table for any present anxiety and behavioural disorder and total problems score at 2-year follow-up 
#

# Merge anxiety and behavioural disorder tables 
ksad_table_2year <- tbl_stack(tbls = list(table_merge_anxiety_t2, table_merge_behavioural_t2)) %>%
  modify_caption("**Table 3. Association between mental health at age 11-12 and new mTBI / orthopaedic injury in the previous 24 months**")

# Save table 
table_2_filename = paste(output_dir, "new_imputed_ksad_disorders_t2.html", sep="")
gt::gtsave(as_gt(ksad_table_2year), file = table_2_filename)
head(table_2_filename)






























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

# Visualise missing data
vis_miss(tbi_data)


# Set seed to make results reproducible
set.seed(7550822)

# Visualise missing data
vis_miss(oi_data)

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



# Create column merging behavioural disorders at baseline and 2 year follow-up 
imputed_psyc_data_t2$any_support_both_years <- paste(imputed_psyc_data_t2$baseline_any_support, imputed_psyc_data_t2$any_support)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'No No'] <- 0
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'No NA'] <- 0
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'NA No'] <- 0
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'Yes NA'] <- 0
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'NA Yes'] <- 1
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'Yes No'] <- 0
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'No Yes'] <- 1
imputed_psyc_data_t2$any_support_both_years[imputed_psyc_data_t2$any_support_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_psyc_data_t2$any_support_both_years <- factor(imputed_psyc_data_t2$any_support_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ any_support_both_years + injury_type, data = imputed_psyc_data_t2)







# Create column merging behavioural disorders at baseline and 2 year follow-up 
imputed_psyc_data_t2$medication_both_years <- paste(imputed_psyc_data_t2$baseline_medication, imputed_psyc_data_t2$medication)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'No No'] <- 0
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'No NA'] <- 0
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'NA No'] <- 0
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'Yes NA'] <- 0
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'NA Yes'] <- 1
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'Yes No'] <- 0
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'No Yes'] <- 1
imputed_psyc_data_t2$medication_both_years[imputed_psyc_data_t2$medication_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_psyc_data_t2$medication_both_years <- factor(imputed_psyc_data_t2$medication_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ medication_both_years + injury_type, data = imputed_psyc_data_t2)





# Create column merging behavioural disorders at baseline and 2 year follow-up 
imputed_psyc_data_t2$psychotherapy_both_years <- paste(imputed_psyc_data_t2$baseline_psychotherapy, imputed_psyc_data_t2$psychotherapy)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'No No'] <- 0
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'No NA'] <- 0
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'NA No'] <- 0
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'Yes NA'] <- 0
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'NA Yes'] <- 1
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'Yes No'] <- 0
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'No Yes'] <- 1
imputed_psyc_data_t2$psychotherapy_both_years[imputed_psyc_data_t2$psychotherapy_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_psyc_data_t2$psychotherapy_both_years <- factor(imputed_psyc_data_t2$psychotherapy_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ psychotherapy_both_years + injury_type, data = imputed_psyc_data_t2)





# Create column merging behavioural disorders at baseline and 2 year follow-up 
imputed_psyc_data_t2$psychotherapy_medication_both_years <- paste(imputed_psyc_data_t2$baseline_psychotherapy_medication, imputed_psyc_data_t2$psychotherapy_medication)

# Recode column so that TBI at baseline is removed (only new cases included)
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'No No'] <- 0
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'No NA'] <- 0
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'NA No'] <- 0
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'Yes NA'] <- 0
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'NA Yes'] <- 1
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'Yes No'] <- 0
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'No Yes'] <- 1
imputed_psyc_data_t2$psychotherapy_medication_both_years[imputed_psyc_data_t2$psychotherapy_medication_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
imputed_psyc_data_t2$psychotherapy_medication_both_years <- factor(imputed_psyc_data_t2$psychotherapy_medication_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ psychotherapy_medication_both_years + injury_type, data = imputed_psyc_data_t2)


# Ensure variables are stored as factors 
imputed_psyc_data_t2$any_support <- as.factor(imputed_psyc_data_t2$any_support)
imputed_psyc_data_t2$outpatient <- as.factor(imputed_psyc_data_t2$outpatient)
imputed_psyc_data_t2$inpatient <- as.factor(imputed_psyc_data_t2$inpatient)
imputed_psyc_data_t2$psychotherapy <- as.factor(imputed_psyc_data_t2$psychotherapy)
imputed_psyc_data_t2$medication <- as.factor(imputed_psyc_data_t2$medication)
imputed_psyc_data_t2$injury_type <- as.factor(imputed_psyc_data_t2$injury_type)
imputed_psyc_data_t2$psychotherapy_medication <- as.factor(imputed_psyc_data_t2$psychotherapy_medication)

# Any mental health support measured at 2-year follow-up 
#

# Change reference category 
imputed_psyc_data_t2$any_support_both_years <- relevel(imputed_psyc_data_t2$any_support_both_years, ref = "No")
imputed_psyc_data_t2$injury_type <- relevel(imputed_psyc_data_t2$injury_type, ref = "None")

# Run unadjusted model 
any_support_model_unadjusted_t2 <- glmer(any_support_both_years ~ injury_type +
                                           (1|site_id_l),
                                         data = imputed_psyc_data_t2,
                                         family = binomial(),
                                         control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
any_support_model_adjusted1_t2 <- glmer(any_support_both_years ~ injury_type +
                                          interview_age +
                                          sex +
                                          race_ethnicity +
                                          household_income +
                                          (1|site_id_l),
                                        data = imputed_psyc_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
any_support_model_adjusted2_t2 <- glmer(any_support_both_years ~ injury_type +
                                          advlife_cat +
                                          safety +
                                          family_conflict_scale +
                                          parent_mh_score +
                                          (1|site_id_l),
                                        data = imputed_psyc_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Psychotherapy measured at 2-year follow-up 
#

# Change reference category 
imputed_psyc_data_t2$psychotherapy_both_years <- relevel(imputed_psyc_data_t2$psychotherapy_both_years, ref = "No")

# Run unadjusted model 
psychotherapy_model_unadjusted_t2 <- glmer(psychotherapy_both_years ~ injury_type + 
                                             (1|site_id_l),
                                           data = imputed_psyc_data_t2,
                                           family = binomial(link = "logit"),
                                           control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_model_adjusted1_t2 <- glmer(psychotherapy_both_years ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income + 
                                            (1|site_id_l),
                                          data = imputed_psyc_data_t2,
                                          family = binomial(link = "logit"),
                                          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
psychotherapy_model_adjusted2_t2 <- glmer(psychotherapy_both_years ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            advlife_cat +
                                            safety +
                                            family_conflict_scale +
                                            parent_mh_score +
                                            (1|site_id_l),
                                          data = imputed_psyc_data_t2,
                                          family = binomial(link = "logit"),
                                          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Medication for mental health measured at 2-year follow-up 
#

# Change reference category 
imputed_psyc_data_t2$medication_both_years <- relevel(imputed_psyc_data_t2$medication_both_years, ref = "No")

# Run unadjusted model 
medication_model_unadjusted_t2 <- glmer(medication_both_years ~ injury_type + 
                                          (1|site_id_l),
                                        data = imputed_psyc_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
medication_model_adjusted1_t2 <- glmer(medication_both_years ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income + 
                                         (1|site_id_l),
                                       data = imputed_psyc_data_t2,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
medication_model_adjusted2_t2 <- glmer(medication_both_years ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income +
                                         advlife_cat +
                                         safety +
                                         family_conflict_scale +
                                         parent_mh_score +
                                         (1|site_id_l),
                                       data = imputed_psyc_data_t2,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Psychotherapy and medication for mental health measured at 2-year follow-up 
#

# Change reference category 
imputed_psyc_data_t2$psychotherapy_medication_both_years <- relevel(imputed_psyc_data_t2$psychotherapy_medication_both_years, ref = "No")

# Run unadjusted model 
psychotherapy_medication_model_unadjusted_t2 <- glmer(psychotherapy_medication_both_years ~ injury_type + 
                                                        (1|site_id_l),
                                                      data = imputed_psyc_data_t2,
                                                      family = binomial(link = "logit"),
                                                      control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_medication_model_adjusted1_t2 <- glmer(psychotherapy_medication_both_years ~ injury_type +
                                                       interview_age +
                                                       sex +
                                                       race_ethnicity +
                                                       household_income + 
                                                       (1|site_id_l),
                                                     data = imputed_psyc_data_t2,
                                                     family = binomial(link = "logit"),
                                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
psychotherapy_medication_model_adjusted2_t2 <- glmer(psychotherapy_medication_both_years ~ injury_type +
                                                       interview_age +
                                                       sex +
                                                       race_ethnicity +
                                                       household_income +
                                                       advlife_cat +
                                                       safety +
                                                       family_conflict_scale +
                                                       parent_mh_score +
                                                       (1|site_id_l),
                                                     data = imputed_psyc_data_t2,
                                                     family = binomial(link = "logit"),
                                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)





# Any mental health support measured at 2-year follow-up 
#

# Create unadjusted model table
any_support_unadjusted_table_t2 <- tbl_regression(any_support_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Any mental health service use") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
any_support_adjusted_table_short_model1_t2 <- tbl_regression(any_support_model_adjusted1_t2, exponentiate = TRUE, 
                                                             pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Any mental health service use") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
any_support_adjusted_table_short_model2_t2 <- tbl_regression(any_support_model_adjusted2_t2, exponentiate = TRUE, 
                                                             pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Any mental health service use") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_any_support_short_t2 <- tbl_merge(tbls = list(any_support_unadjusted_table_t2, any_support_adjusted_table_short_model1_t2, any_support_adjusted_table_short_model2_t2),
                                              tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Psychotherapy measured at 2-year follow-up
#

# Create unadjusted model table
psychotherapy_unadjusted_table_t2 <- tbl_regression(psychotherapy_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Psychotherapy") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
psychotherapy_adjusted_table_short_model1_t2 <- tbl_regression(psychotherapy_model_adjusted1_t2, exponentiate = TRUE, 
                                                               pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Psychotherapy") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
psychotherapy_adjusted_table_short_model2_t2 <- tbl_regression(psychotherapy_model_adjusted2_t2, exponentiate = TRUE, 
                                                               pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Psychotherapy") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))


# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_psychotherapy_short_t2 <- tbl_merge(tbls = list(psychotherapy_unadjusted_table_t2, psychotherapy_adjusted_table_short_model1_t2, psychotherapy_adjusted_table_short_model2_t2),
                                                tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Medication measured at 2-year follow-up
#

# Create unadjusted model table
medication_unadjusted_table_t2 <- tbl_regression(medication_model_unadjusted_t2, exponentiate = TRUE, 
                                                 label = injury_type ~ "Medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
medication_adjusted_table_short_model1_t2 <- tbl_regression(medication_model_adjusted1_t2, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type,
                                                            label = injury_type ~ "Medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
medication_adjusted_table_short_model2_t2 <- tbl_regression(medication_model_adjusted2_t2, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type,
                                                            label = injury_type ~ "Medication for mental health") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_medication_short_t2 <- tbl_merge(tbls = list(medication_unadjusted_table_t2, medication_adjusted_table_short_model1_t2, medication_adjusted_table_short_model2_t2),
                                             tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Psychotherapy and medication measured at 2-year follow-up
#

# Create unadjusted model table
psychotherapy_medication_unadjusted_table_t2 <- tbl_regression(psychotherapy_medication_model_unadjusted_t2, exponentiate = TRUE, 
                                                               label = injury_type ~ "Psychotherapy and/or medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
psychotherapy_medication_adjusted_table_short_model1_t2 <- tbl_regression(psychotherapy_medication_model_adjusted1_t2, exponentiate = TRUE, 
                                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                          include = injury_type,
                                                                          label = injury_type ~ "Psychotherapy and/or medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
psychotherapy_medication_adjusted_table_short_model2_t2 <- tbl_regression(psychotherapy_medication_model_adjusted2_t2, exponentiate = TRUE, 
                                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                          include = injury_type,
                                                                          label = injury_type ~ "Psychotherapy and/or medication for mental health") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_psychotherapy_medication_short_t2 <- tbl_merge(tbls = list(psychotherapy_medication_unadjusted_table_t2, psychotherapy_medication_adjusted_table_short_model1_t2, 
                                                                       psychotherapy_medication_adjusted_table_short_model2_t2),
                                                           tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# # Create table for psychiatric service use at 2-year follow-up
#

# Merge  psychotherapy, medication and any mental health support tables 
psyc_service_t2_table <- tbl_stack(tbls = list(table_merge_psychotherapy_short_t2, table_merge_medication_short_t2, table_merge_psychotherapy_medication_short_t2, table_merge_any_support_short_t2)) %>%
  modify_caption("**Table 5. Association between new mTBI / orthopaedic injury and mental health service use in the 12-24 months following baseline (age 11-12)**")

# Save table 
table_2_filename = paste(output_dir, "new_imputed_psyc_service_t2.html", sep="")
gt::gtsave(as_gt(psyc_service_t2_table), file = table_2_filename)
head(table_2_filename)




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


# Visualise missing data
vis_miss(tbi_psyc_data)
vis_miss(oi_psyc_data)

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

# KSADS Anxiety (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$anxiety_disorders <- as.factor(matched_data$anxiety_disorders)

# Run model 
tbi_anxiety_disorders_model_adjusted1 <- glm(anxiety_disorders ~ tbi_present + 
                                               (1|site_id_l),
                                         data = matched_data,
                                         family = binomial(link = "logit"))
tab_model(tbi_anxiety_disorders_model_adjusted1)


# KSADS Anxiety (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income +
                   advlife_cat + safety + family_conflict_scale +
                   baseline_anxiety_disorders,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$anxiety_disorders <- as.factor(matched_data$anxiety_disorders)

# Run model 
oi_anxiety_disorders_model_adjusted1 <- glm(anxiety_disorders ~ oi_present,
                                         data = matched_data,
                                         family = binomial(link = "logit"))
tab_model(oi_anxiety_disorders_model_adjusted1)



# KSADS Behavioural Disorders (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$behavioural_disorders <- as.factor(matched_data$behavioural_disorders)

# Run model 
tbi_behavioural_disorders_model_adjusted1 <- glm(behavioural_disorders ~ tbi_present, 
                                             data = matched_data,
                                             family = binomial(link = "logit"))
tab_model(tbi_behavioural_disorders_model_adjusted1)


tbi_behavioural_disorders_model_adjusted1 <- glmer(behavioural_disorders ~ tbi_present + 
                                                 (1|site_id_l),
                                               data = matched_data,
                                               family = binomial(link = "logit"))
tab_model(tbi_behavioural_disorders_model_adjusted1)



# KSADS Behavioural Disorders (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_behavioural_disorders,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$behavioural_disorders <- as.factor(matched_data$behavioural_disorders)

# Run model 
oi_behavioural_disorders_model_adjusted1 <- glm(behavioural_disorders ~ oi_present, 
                                             data = matched_data,
                                             family = binomial(link = "logit"))
tab_model(oi_behavioural_disorders_model_adjusted1)


oi_behavioural_disorders_model_adjusted1 <- glmer(behavioural_disorders ~ oi_present + 
                                                 (1|site_id_l),
                                               data = matched_data,
                                               family = binomial(link = "logit"))
tab_model(oi_behavioural_disorders_model_adjusted1)



# CBCL Total Problems (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Run model
tbi_cbcl_tot_symptoms_model_adjusted1 <- lm(cbcl_tot_problems ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_tot_symptoms_model_adjusted1)



# CBCL Total Problems (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_tot,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Run model
oi_cbcl_tot_symptoms_model_adjusted1 <- lm(cbcl_tot_problems ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_tot_symptoms_model_adjusted1)


# CBCL subscales 

# CBCL internalising symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Run model
tbi_cbcl_int_symptoms_model_adjusted1 <- lm(cbcl_int_problems ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_int_symptoms_model_adjusted1)



# CBCL internalising symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_int,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Run model
oi_cbcl_int_symptoms_model_adjusted1 <- lm(cbcl_int_problems ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_int_symptoms_model_adjusted1)


# CBCL externalising symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_ext_symptoms_model_adjusted1 <- lm(cbcl_ext_problems ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_ext_symptoms_model_adjusted1)


# CBCL externalising symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_ext,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_ext_symptoms_model_adjusted1 <- lm(cbcl_ext_problems ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_ext_symptoms_model_adjusted1)


# CBCL anxiety symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_anx_symptoms_model_adjusted1 <- lm(cbcl_anx_score ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_anx_symptoms_model_adjusted1)


# CBCL anxiety symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_anx,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_anx_symptoms_model_adjusted1 <- lm(cbcl_anx_score ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_anx_symptoms_model_adjusted1)



# CBCL depression symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_dep_symptoms_model_adjusted1 <- lm(cbcl_dep_score ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_dep_symptoms_model_adjusted1)



# CBCL depression symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_dep,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_dep_symptoms_model_adjusted1 <- lm(cbcl_dep_score ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_dep_symptoms_model_adjusted1)


# CBCL ADHD symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_adhd_symptoms_model_adjusted1 <- lm(cbcl_adhd_score ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_adhd_symptoms_model_adjusted1)




# CBCL ADHD symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_adhd,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_adhd_symptoms_model_adjusted1 <- lm(cbcl_adhd_score ~ oi_present,
                                         data = matched_data)
tab_model(oi_cbcl_adhd_symptoms_model_adjusted1)


# CBCL ODD symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_odd_symptoms_model_adjusted1 <- lm(cbcl_odd_score ~ tbi_present,
                                         data = matched_data)
tab_model(tbi_cbcl_odd_symptoms_model_adjusted1)



# CBCL ODD symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_odd,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_odd_symptoms_model_adjusted1 <- lm(cbcl_odd_score ~ oi_present,
                                        data = matched_data)
tab_model(oi_cbcl_odd_symptoms_model_adjusted1)



# CBCL conduct symptoms (TBI vs no injury)
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = tbi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
tbi_cbcl_conduct_symptoms_model_adjusted1 <- lm(cbcl_conduct_score ~ tbi_present,
                                        data = matched_data)
tab_model(tbi_cbcl_conduct_symptoms_model_adjusted1)




# CBCL conduct symptoms (OI vs no injury)
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_cbcl_conduct,
                 data = oi_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)


# Run model
oi_cbcl_conduct_symptoms_model_adjusted1 <- lm(cbcl_conduct_score ~ oi_present,
                                            data = matched_data)
tab_model(oi_cbcl_conduct_symptoms_model_adjusted1)

###########################################################
#
# Psychiatric service use propensity score matching using MatchIt
#
###########################################################

# Any mental health support 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$any_support <- as.factor(matched_data$any_support)

# Run model 
tbi_any_support_model_adjusted1 <- glm(any_support ~ tbi_present, 
                                             data = matched_data,
                                             family = binomial(link = "logit"))
tab_model(tbi_any_support_model_adjusted1)



# Any mental health support 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_any_support,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$any_support <- as.factor(matched_data$any_support)

# Run model 
oi_any_support_model_adjusted1 <- glm(any_support ~ oi_present, 
                                       data = matched_data,
                                       family = binomial(link = "logit"))
tab_model(oi_any_support_model_adjusted1)



# Outpatient support 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$outpatient <- as.factor(matched_data$outpatient)

# Run model 
tbi_outpatient_model_adjusted1 <- glm(outpatient ~ tbi_present, 
                                   data = matched_data,
                                   family = binomial(link = "logit"))
tab_model(tbi_outpatient_model_adjusted1)




# Outpatient support 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_outpatient,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$outpatient <- as.factor(matched_data$outpatient)

# Run model 
oi_outpatient_model_adjusted1 <- glm(outpatient ~ oi_present, 
                                      data = matched_data,
                                      family = binomial(link = "logit"))
tab_model(oi_outpatient_model_adjusted1)







# Inpatient support 
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$inpatient <- as.factor(matched_data$inpatient)

# Run model 
tbi_inpatient_model_adjusted1 <- glm(inpatient ~ tbi_present, 
                                  data = matched_data,
                                  family = binomial(link = "logit"))
tab_model(tbi_inpatient_model_adjusted1)


#inpatient_model_adjusted1 <- glmer(inpatient ~ tbi_present + 
#                                      (1|site_id_l),
#                                    data = matched_data,
#                                    family = binomial(link = "logit"))
#tab_model(inpatient_model_adjusted1)


# Inpatient support 
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_inpatient,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$inpatient <- as.factor(matched_data$inpatient)

# Run model 
oi_inpatient_model_adjusted1 <- glm(inpatient ~ oi_present, 
                                     data = matched_data,
                                     family = binomial(link = "logit"))
tab_model(oi_inpatient_model_adjusted1)

# Psychotherapy
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$psychotherapy <- as.factor(matched_data$psychotherapy)

# Run model 
tbi_psychotherapy_model_adjusted1 <- glm(psychotherapy ~ tbi_present, 
                                 data = matched_data,
                                 family = binomial(link = "logit"))
tab_model(tbi_psychotherapy_model_adjusted1)


#psychotherapy_model_adjusted1 <- glmer(psychotherapy ~ tbi_present + 
#                                     (1|site_id_l),
#                                   data = matched_data,
#                                   family = binomial(link = "logit"))
#tab_model(psychotherapy_model_adjusted1)



# Psychotherapy
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$psychotherapy <- as.factor(matched_data$psychotherapy)

# Run model 
oi_psychotherapy_model_adjusted1 <- glm(psychotherapy ~ oi_present, 
                                         data = matched_data,
                                         family = binomial(link = "logit"))
tab_model(oi_psychotherapy_model_adjusted1)



# Medication
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$medication <- as.factor(matched_data$medication)

# Run model 
tbi_medication_model_adjusted1 <- glm(medication ~ tbi_present, 
                                     data = matched_data,
                                     family = binomial(link = "logit"))
tab_model(tbi_medication_model_adjusted1)


#medication_model_adjusted1 <- glmer(medication ~ tbi_present + 
#                                         (1|site_id_l),
#                                      data = matched_data,
#                                      family = binomial(link = "logit"))
#tab_model(medication_model_adjusted1)





# Medication
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_medication,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$medication <- as.factor(matched_data$medication)

# Run model 
oi_medication_model_adjusted1 <- glm(medication ~ oi_present, 
                                      data = matched_data,
                                      family = binomial(link = "logit"))
tab_model(oi_medication_model_adjusted1)


# Psychotherapy and medication
model <- matchit(tbi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = tbi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$psychotherapy_medication <- as.factor(matched_data$psychotherapy_medication)

# Run model 
tbi_psychotherapy_medication_model_adjusted1 <- glm(psychotherapy_medication ~ tbi_present, 
                                  data = matched_data,
                                  family = binomial(link = "logit"))
tab_model(tbi_psychotherapy_medication_model_adjusted1)


#psychotherapy_medication_model_adjusted1 <- glmer(psychotherapy_medication ~ tbi_present + 
#                                      (1|site_id_l),
#                                    data = matched_data,
#                                    family = binomial(link = "logit"))
#tab_model(psychotherapy_medication_model_adjusted1)



# Psychotherapy and medication
model <- matchit(oi_present ~
                   interview_age + male + race_ethnicity + household_income + 
                   advlife_cat + safety + family_conflict_scale +
                   baseline_psychotherapy_medication,
                 data = oi_psyc_data,
                 method = "optimal", distance = "glm")

summary(model, standardize = TRUE)

matched_data <- match.data(model)

# Make sure variable is a factor 
matched_data$psychotherapy_medication <- as.factor(matched_data$psychotherapy_medication)

# Run model 
oi_psychotherapy_medication_model_adjusted1 <- glm(psychotherapy_medication ~ oi_present, 
                                                    data = matched_data,
                                                    family = binomial(link = "logit"))
tab_model(oi_psychotherapy_medication_model_adjusted1)


##############################################################################
#
# Create tables 
#
##############################################################################

#
# Tables for KSADS-5 analyses 
tbi_anxiety_disorder <- tbl_regression(tbi_anxiety_disorders_model_adjusted1, exponentiate = TRUE,  
                                   pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                   label = tbi_present ~ "Any present anxiety disorder")

tbi_behavioural_disorder <- tbl_regression(tbi_behavioural_disorders_model_adjusted1, exponentiate = TRUE, 
                                       pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                       label = tbi_present ~ "Any present behavioural disorder")

oi_anxiety_disorder <- tbl_regression(oi_anxiety_disorders_model_adjusted1, exponentiate = TRUE,  
                                       pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                       label = oi_present ~ "Any present anxiety disorder")


oi_behavioural_disorder <- tbl_regression(oi_behavioural_disorders_model_adjusted1, exponentiate = TRUE, 
                                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                           label = oi_present ~ "Any present behavioural disorder")


# Merge anxiety and behavioural disorder tables 
tbi_ksad_table <- tbl_stack(tbls = list(tbi_anxiety_disorder, tbi_behavioural_disorder)) %>%
  modify_caption("**Association between mental health and new mTBI using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "tbi_ksad_psm.html", sep="")
gt::gtsave(as_gt(tbi_ksad_table), file = table_2_filename)
head(table_2_filename)




# Merge anxiety and behavioural disorder tables 
oi_ksad_table <- tbl_stack(tbls = list(oi_anxiety_disorder, oi_behavioural_disorder)) %>%
  modify_caption("**Association between mental health and new Ortho using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "oi_ksad_psm.html", sep="")
gt::gtsave(as_gt(oi_ksad_table), file = table_2_filename)
head(table_2_filename)


#
# Tables for CBCL analyses 
tbi_cbcl_tot <- tbl_regression(tbi_cbcl_tot_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Total problems score")

tbi_cbcl_int <- tbl_regression(tbi_cbcl_int_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Internalising symptoms")


tbi_cbcl_ext <- tbl_regression(tbi_cbcl_ext_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Externalising symptoms")

tbi_cbcl_anx <- tbl_regression(tbi_cbcl_anx_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Anxiety symptoms")

tbi_cbcl_dep <- tbl_regression(tbi_cbcl_dep_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Depression symptoms")

tbi_cbcl_adhd <- tbl_regression(tbi_cbcl_adhd_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "ADHD symptoms")

tbi_cbcl_conduct <- tbl_regression(tbi_cbcl_conduct_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "Conduct symptoms")

tbi_cbcl_odd <- tbl_regression(tbi_cbcl_odd_symptoms_model_adjusted1, 
                           pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                           label = tbi_present ~ "ODD symptoms")



#
# Tables for CBCL analyses 
oi_cbcl_tot <- tbl_regression(oi_cbcl_tot_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "Total problems score")

oi_cbcl_int <- tbl_regression(oi_cbcl_int_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "Internalising symptoms")


oi_cbcl_ext <- tbl_regression(oi_cbcl_ext_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "Externalising symptoms")

oi_cbcl_anx <- tbl_regression(oi_cbcl_anx_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "Anxiety symptoms")

oi_cbcl_dep <- tbl_regression(oi_cbcl_dep_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "Depression symptoms")

oi_cbcl_adhd <- tbl_regression(oi_cbcl_adhd_symptoms_model_adjusted1, 
                                pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                label = oi_present ~ "ADHD symptoms")

oi_cbcl_conduct <- tbl_regression(oi_cbcl_conduct_symptoms_model_adjusted1, 
                                   pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                   label = oi_present ~ "Conduct symptoms")

oi_cbcl_odd <- tbl_regression(oi_cbcl_odd_symptoms_model_adjusted1, 
                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                               label = oi_present ~ "ODD symptoms")

# Merge CBCL score tables 
tbi_cbcl_table <- tbl_stack(tbls = list(tbi_cbcl_tot, tbi_cbcl_int, tbi_cbcl_ext, tbi_cbcl_anx, tbi_cbcl_dep, tbi_cbcl_adhd, 
                                    tbi_cbcl_conduct, tbi_cbcl_odd)) %>%
  modify_caption("**Association between mental health and new mTBI using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "tbi_cbcl_psm.html", sep="")
gt::gtsave(as_gt(tbi_cbcl_table), file = table_2_filename)
head(table_2_filename)

# Merge CBCL score tables 
oi_cbcl_table <- tbl_stack(tbls = list(oi_cbcl_tot, oi_cbcl_int, oi_cbcl_ext, oi_cbcl_anx, oi_cbcl_dep, oi_cbcl_adhd, 
                                       oi_cbcl_conduct, oi_cbcl_odd)) %>%
  modify_caption("**Association between mental health and new Ortho using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "oi_cbcl_psm.html", sep="")
gt::gtsave(as_gt(oi_cbcl_table), file = table_2_filename)
head(table_2_filename)


#
# Tables for psychiatric service use analyses 
tbi_psychotherapy <- tbl_regression(tbi_psychotherapy_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Psychotherapy")

tbi_medication  <- tbl_regression(tbi_medication_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Medication for mental health")

tbi_psychotherapy_medication <- tbl_regression(tbi_psychotherapy_medication_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Psychotherapy and/or medication")


tbi_outpatient <- tbl_regression(tbi_outpatient_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Outpatient")


tbi_inpatient <- tbl_regression(tbi_inpatient_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Inpatient")


tbi_any_support <- tbl_regression(tbi_any_support_model_adjusted1, exponentiate = TRUE,  
                              pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                              label = tbi_present ~ "Any mental health support")


# Merge psyc service tables 
tbi_psyc_service_table <- tbl_stack(tbls = list(tbi_psychotherapy, tbi_medication, tbi_psychotherapy_medication, tbi_outpatient, 
                                            tbi_inpatient, tbi_any_support)) %>%
  modify_caption("**Association between mental health and new mTBI using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "tbi_psyc_service_psm.html", sep="")
gt::gtsave(as_gt(tbi_psyc_service_table), file = table_2_filename)
head(table_2_filename)



#
# Tables for psychiatric service use analyses 
oi_psychotherapy <- tbl_regression(oi_psychotherapy_model_adjusted1, exponentiate = TRUE,  
                                    pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                    label = oi_present ~ "Psychotherapy")

oi_medication  <- tbl_regression(oi_medication_model_adjusted1, exponentiate = TRUE,  
                                  pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                  label = oi_present ~ "Medication for mental health")

oi_psychotherapy_medication <- tbl_regression(oi_psychotherapy_medication_model_adjusted1, exponentiate = TRUE,  
                                               pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                               label = oi_present ~ "Psychotherapy and/or medication")


oi_outpatient <- tbl_regression(oi_outpatient_model_adjusted1, exponentiate = TRUE,  
                                 pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                 label = oi_present ~ "Outpatient")


oi_inpatient <- tbl_regression(oi_inpatient_model_adjusted1, exponentiate = TRUE,  
                                pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                label = oi_present ~ "Inpatient")


oi_any_support <- tbl_regression(oi_any_support_model_adjusted1, exponentiate = TRUE,  
                                  pvalue_fun = function(x) style_pvalue(x, digits = 2), 
                                  label = oi_present ~ "Any mental health support")


# Merge psyc service tables 
oi_psyc_service_table <- tbl_stack(tbls = list(oi_psychotherapy, oi_medication, oi_psychotherapy_medication, oi_outpatient, 
                                                oi_inpatient, oi_any_support)) %>%
  modify_caption("**Association between mental health and new Ortho using propensity score matching**")

# Save table 
table_2_filename = paste(output_dir, "oi_psyc_service_psm.html", sep="")
gt::gtsave(as_gt(oi_psyc_service_table), file = table_2_filename)
head(table_2_filename)


