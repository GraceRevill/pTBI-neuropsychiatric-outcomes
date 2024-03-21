################################################################
#
# Analysis of psychiatric service use from the study:
#
# Mild TBI but not Orthopaedic Injury Associated with Subsequent
# Psychiatric Problems and Increased Mental Health Service Use in
# 9-12-year-old Children
#
# Revill et al
#
################################################################

#
# This script loads ABCD 4.0 data from the data_dir directory
#
# It writes imputed datasets and publication ready tables to the output_dir directory
#

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(survey)
library(lme4)
library(Hmisc)
library(table1)
library(epiR)
library(jtools)
library(gtsummary)
library(gt)
library(afex)
library(missForest)
library(doParallel)
library(visdat)
library(doRNG)
library(ggplot2)

# Remove everything from memory if needed to start afresh
# rm(list = ls())

# Set data and output directories
data_dir = "S:/ABCDPsychosis/ABCD4.0/"
output_dir = "S:/ABCDPsychosis/Grace/Output/"

setwd(data_dir)

# The following function uses readr to read in an ABCD dataframe to address some formatting issues - 
# the file was tab delimited/separated plus the first row was the description of variable which
# set the variable type wrong. The code below reads in the column names, then reads in the datafile
# ignoring the first two rows. It then uses column names which are already read in, meaning 
# that read_delim can identify the data type when it reads in the dataframe

read_abcd_file <- function(abcd_filename) {
  col_names <- names(read_delim(abcd_filename, delim="\t", n_max = 0))
  abcd_df <- read_delim(abcd_filename, delim="\t", na = "", col_names = col_names, skip = 3, guess_max = 2000)
  return(abcd_df)
}

###########################################################
#
# Load data
#
###########################################################

# TBI
tbi_baseline <- read_abcd_file("abcd_tbi01.txt")
tbi_long <- read_abcd_file("abcd_lsstbi01.txt")

# Orthopaedic injury
oi_baseline <- read_abcd_file("abcd_mx01.txt")
oi_long <- read_abcd_file("abcd_lssmh01.txt")

# Psychiatric service use 
psyc_services_baseline <- read_abcd_file("dibf01.txt")
psyc_services_long <- read_abcd_file("abcd_lpksad01.txt")

# Covariates
ethnicity <- read_abcd_file("acspsw03.txt")
parent_demographics <- read_abcd_file("pdem02.txt") 
advlife <- read_abcd_file("abcd_ptsd01.txt")
parenting <- read_abcd_file("crpbi01.txt")
site <- read_abcd_file("abcd_lt01.txt")
neighbourhood_safety <- read_abcd_file("abcd_nsc01.txt")
family_conflict <- read_abcd_file("abcd_fes01.txt")
parent_mh <- read_abcd_file("abcd_asrs01.txt")

###########################################################
#
# Rename and recode variables
#
###########################################################

# Rename and select TBI severity at baseline
tbi_baseline <- tbi_baseline %>%
  rename(tbi_severity = tbi_ss_worst_overall) %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(subjectkey, interview_age, sex, tbi_severity)  

# Recode TBI severity variable at baseline (include only possible mild and mild TBI cases)
tbi_baseline$tbi_severity[tbi_baseline$tbi_severity == 1] <- 0 #None/improbable 
tbi_baseline$tbi_severity[tbi_baseline$tbi_severity == 2] <- 1 #Possible 
tbi_baseline$tbi_severity[tbi_baseline$tbi_severity == 3] <- 1 #Mild 
tbi_baseline$tbi_severity[tbi_baseline$tbi_severity == 4] <- NA #Moderate 
tbi_baseline$tbi_severity[tbi_baseline$tbi_severity == 5] <- NA #Severe 

# Set data type and level labels
tbi_baseline$tbi_severity <- factor(tbi_baseline$tbi_severity, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select TBI severity at 1-year follow-up
tbi_1_year <- tbi_long %>%
  rename(tbi_severity_1_year = tbi_ss_worst_overall_l) %>%
  filter(eventname == '1_year_follow_up_y_arm_1') %>%
  select(subjectkey, tbi_severity_1_year)  

# Recode TBI severity variable at 1-year follow-up (include only possible mild and mild TBI cases)
tbi_1_year$tbi_severity_1_year[tbi_1_year$tbi_severity_1_year == 1] <- 0 #None/improbable 
tbi_1_year$tbi_severity_1_year[tbi_1_year$tbi_severity_1_year == 2] <- 1 #Possible 
tbi_1_year$tbi_severity_1_year[tbi_1_year$tbi_severity_1_year == 3] <- 1 #Mild 
tbi_1_year$tbi_severity_1_year[tbi_1_year$tbi_severity_1_year == 4] <- NA #Moderate 
tbi_1_year$tbi_severity_1_year[tbi_1_year$tbi_severity_1_year == 5] <- NA #Severe 

# Set data type and level labels
tbi_1_year$tbi_severity_1_year <- factor(tbi_1_year$tbi_severity_1_year, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select TBI severity at 2-year follow-up
tbi_2_year <- tbi_long %>%
  rename(tbi_severity_2_year = tbi_ss_worst_overall_l) %>%
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  select(subjectkey, interview_age, sex, tbi_severity_2_year)  

# Recode TBI severity variable at 2-year follow-up (include only possible mild and mild TBI cases)
tbi_2_year$tbi_severity_2_year[tbi_2_year$tbi_severity_2_year == 1] <- 0 #None/improbable 
tbi_2_year$tbi_severity_2_year[tbi_2_year$tbi_severity_2_year == 2] <- 1 #Possible 
tbi_2_year$tbi_severity_2_year[tbi_2_year$tbi_severity_2_year == 3] <- 1 #Mild 
tbi_2_year$tbi_severity_2_year[tbi_2_year$tbi_severity_2_year == 4] <- NA #Moderate 
tbi_2_year$tbi_severity_2_year[tbi_2_year$tbi_severity_2_year == 5] <- NA #Severe 

# Set data type and level labels
tbi_2_year$tbi_severity_2_year <- factor(tbi_2_year$tbi_severity_2_year, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select orthopaedic injury at baseline
oi_baseline <- oi_baseline %>%
  rename(broken_bones = medhx_6a) %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(subjectkey, broken_bones)

# Recode orthopaedic injury variable at baseline - yes (1) or no (0) 
oi_baseline$broken_bones[oi_baseline$broken_bones == 0] <- 0
oi_baseline$broken_bones[oi_baseline$broken_bones == 1] <- 1

# Set data type and level labels
oi_baseline$broken_bones <- factor(oi_baseline$broken_bones, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select orthopaedic injury at 1-year follow-up 
oi_1_year <- oi_long %>%
  rename(broken_bones_1_year = medhx_ss_6a_times_p_l) %>%
  filter(eventname == '1_year_follow_up_y_arm_1') %>%
  select(subjectkey, broken_bones_1_year)

# Recode orthopaedic injury variable at 1-year follow-up - yes (1) or no (0)
oi_1_year$broken_bones_1_year[oi_1_year$broken_bones_1_year == 0] <- 0
oi_1_year$broken_bones_1_year[oi_1_year$broken_bones_1_year == 1] <- 1

# Set data type and level labels
oi_1_year$broken_bones_1_year <- factor(oi_1_year$broken_bones_1_year, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select orthopaedic injury at 2-year follow-up 
oi_2_year <- oi_long %>%
  rename(broken_bones_2_year = medhx_ss_6a_times_p_l) %>%
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  select(subjectkey, broken_bones_2_year)

# Recode orthopaedic injury variable at 2-year follow-up - yes (1) or no (0)
oi_2_year$broken_bones_2_year[oi_2_year$broken_bones_2_year == 0] <- 0
oi_2_year$broken_bones_2_year[oi_2_year$broken_bones_2_year == 1] <- 1

# Set data type and level labels
oi_2_year$broken_bones_2_year <- factor(oi_2_year$broken_bones_2_year, levels = c(1,0), labels = c("Yes", "No"))

# Rename and select psychiatric service use variables at baseline
psyc_services_baseline <- psyc_services_baseline %>%
  rename(any_support = kbi_p_c_mh_sa) %>% 
  rename(outpatient = kbi_p_c_scheck1) %>%
  rename(inpatient = kbi_p_c_scheck3) %>%
  rename(psychotherapy = kbi_p_c_scheck7) %>%
  rename(medication = kbi_p_c_scheck8) %>%
  select(subjectkey, any_support, outpatient, inpatient, 
         psychotherapy, medication)  

# Recode any mental health support at baseline as yes (1) or no (0)
psyc_services_baseline$any_support[psyc_services_baseline$any_support == 0] <- 0
psyc_services_baseline$any_support[psyc_services_baseline$any_support == 1] <- 1
psyc_services_baseline$any_support[psyc_services_baseline$any_support == 3] <- 0

# Set data type and level labels
psyc_services_baseline$any_support <- factor(psyc_services_baseline$any_support, levels = c(1,0), labels = c("Yes", "No"))

# Recode outpatient support at baseline as yes (1) or no (0)
psyc_services_baseline$outpatient[psyc_services_baseline$outpatient == 0] <- 0
psyc_services_baseline$outpatient[psyc_services_baseline$outpatient == 1] <- 1
psyc_services_baseline$outpatient[psyc_services_baseline$outpatient == 3] <- 0

# Set data type and level labels
psyc_services_baseline$outpatient <- factor(psyc_services_baseline$outpatient, levels = c(1,0), labels = c("Yes", "No"))

# Recode inpatient support at baseline as yes (1) or no (0)
psyc_services_baseline$inpatient[psyc_services_baseline$inpatient == 0] <- 0
psyc_services_baseline$inpatient[psyc_services_baseline$inpatient == 1] <- 1
psyc_services_baseline$inpatient[psyc_services_baseline$inpatient == 3] <- 0

# Set data type and level labels
psyc_services_baseline$inpatient <- factor(psyc_services_baseline$inpatient, levels = c(1,0), labels = c("Yes", "No"))

# Recode psychotherapy at baseline as yes (1) or no (0)
psyc_services_baseline$psychotherapy[psyc_services_baseline$psychotherapy == 0] <- 0
psyc_services_baseline$psychotherapy[psyc_services_baseline$psychotherapy == 1] <- 1
psyc_services_baseline$psychotherapy[psyc_services_baseline$psychotherapy == 3] <- 0

# Set data type and level labels
psyc_services_baseline$psychotherapy <- factor(psyc_services_baseline$psychotherapy, levels = c(1,0), labels = c("Yes", "No"))

# Recode medication for mental health at baseline as yes (1) or no (0)
psyc_services_baseline$medication[psyc_services_baseline$medication == 0] <- 0
psyc_services_baseline$medication[psyc_services_baseline$medication == 1] <- 1
psyc_services_baseline$medication[psyc_services_baseline$medication == 3] <- 0

# Set data type and level labels
psyc_services_baseline$medication <- factor(psyc_services_baseline$medication, levels = c(1,0), labels = c("Yes", "No"))


# Create overall psychotherapy and medication variable at baseline using psychotherapy and medication columns 
psyc_services_baseline$psychotherapy_medication <- paste(psyc_services_baseline$psychotherapy, psyc_services_baseline$medication)

# Recode injury type column 
psyc_services_baseline$psychotherapy_medication [psyc_services_baseline$psychotherapy_medication  == 'No No'] <- 0
psyc_services_baseline$psychotherapy_medication [psyc_services_baseline$psychotherapy_medication  == 'Yes No'] <- 1
psyc_services_baseline$psychotherapy_medication [psyc_services_baseline$psychotherapy_medication  == 'No Yes'] <- 1
psyc_services_baseline$psychotherapy_medication [psyc_services_baseline$psychotherapy_medication  == 'Yes Yes'] <- 1

# Set data type and level labels
psyc_services_baseline$psychotherapy_medication <- factor(psyc_services_baseline$psychotherapy_medication, levels = c(1,0), labels = c("Yes", "No"))


# Create injury type variable at baseline using TBI severity and ortho injury columns 
psyc_data$injury_type <- paste(psyc_data$tbi_severity, psyc_data$broken_bones)

# Recode injury type column 
psyc_data$injury_type[psyc_data$injury_type == 'No No'] <- 0
psyc_data$injury_type[psyc_data$injury_type == 'Yes No'] <- 2
psyc_data$injury_type[psyc_data$injury_type == 'No Yes'] <- 1
psyc_data$injury_type[psyc_data$injury_type == 'Yes Yes'] <- 2

# Set data type and level labels
psyc_data$injury_type <- factor(psyc_data$injury_type, levels = c(2,1,0), labels = c("TBI", "Ortho", "None" ))

# Rename and select psychiatric service use variables at 2-year follow-up 
psyc_services_2_year <- psyc_services_long %>%
  rename(any_support = kbi_p_c_mh_sa_l) %>%
  rename(outpatient = kbipcserviceschecklistl1) %>%
  rename(inpatient = kbipcserviceschecklistl3) %>%
  rename(psychotherapy = kbipcserviceschecklistl7) %>%
  rename(medication = kbipcserviceschecklistl8) %>%
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  select(subjectkey, any_support,outpatient, inpatient, 
         psychotherapy, medication)

# Recode any mental health support at follow-up as yes (1) or no (0)
psyc_services_2_year$any_support[psyc_services_2_year$any_support == 2] <- 0
psyc_services_2_year$any_support[psyc_services_2_year$any_support == 1] <- 1
psyc_services_2_year$any_support[psyc_services_2_year$any_support == 3] <- 0

# Set data type and level labels
psyc_services_2_year$any_support <- factor(psyc_services_2_year$any_support, levels = c(1,0), labels = c("Yes", "No"))

# Recode outpatient support at follow-up as yes (1) or no (0)
psyc_services_2_year$outpatient[psyc_services_2_year$outpatient == 0] <- 0
psyc_services_2_year$outpatient[psyc_services_2_year$outpatient == 1] <- 1
psyc_services_2_year$outpatient[psyc_services_2_year$outpatient == 3] <- 0

# Set data type and level labels
psyc_services_2_year$outpatient <- factor(psyc_services_2_year$outpatient, levels = c(1,0), labels = c("Yes", "No"))



# Recode inpatient support at follow-up as yes (1) or no (0)
psyc_services_2_year$inpatient[psyc_services_2_year$inpatient == 0] <- 0
psyc_services_2_year$inpatient[psyc_services_2_year$inpatient == 1] <- 1
psyc_services_2_year$inpatient[psyc_services_2_year$inpatient == 3] <- 0

# Set data type and level labels
psyc_services_2_year$inpatient <- factor(psyc_services_2_year$inpatient, levels = c(1,0), labels = c("Yes", "No"))


# Recode psychotherapy at follow-up as yes (1) or no (0)
psyc_services_2_year$psychotherapy[psyc_services_2_year$psychotherapy == 0] <- 0
psyc_services_2_year$psychotherapy[psyc_services_2_year$psychotherapy == 1] <- 1
psyc_services_2_year$psychotherapy[psyc_services_2_year$psychotherapy == 3] <- 0

# Set data type and level labels
psyc_services_2_year$psychotherapy <- factor(psyc_services_2_year$psychotherapy, levels = c(1,0), labels = c("Yes", "No"))


# Recode medication for mental health at follow-up as yes (1) or no (0)
psyc_services_2_year$medication[psyc_services_2_year$medication == 0] <- 0
psyc_services_2_year$medication[psyc_services_2_year$medication == 1] <- 1
psyc_services_2_year$medication[psyc_services_2_year$medication == 3] <- 0

# Set data type and level labels
psyc_services_2_year$medication <- factor(psyc_services_2_year$medication, levels = c(1,0), labels = c("Yes", "No"))

# Create overall psychotherapy and medication variable at follow-up using psychotherapy and medication columns 
psyc_services_2_year$psychotherapy_medication <- paste(psyc_services_2_year$psychotherapy, psyc_services_2_year$medication)

# Recode injury type column 
psyc_services_2_year$psychotherapy_medication [psyc_services_2_year$psychotherapy_medication  == 'No No'] <- 0
psyc_services_2_year$psychotherapy_medication [psyc_services_2_year$psychotherapy_medication  == 'Yes No'] <- 1
psyc_services_2_year$psychotherapy_medication [psyc_services_2_year$psychotherapy_medication  == 'No Yes'] <- 1
psyc_services_2_year$psychotherapy_medication [psyc_services_2_year$psychotherapy_medication  == 'Yes Yes'] <- 1

# Set data type and level labels
psyc_services_2_year$psychotherapy_medication <- factor(psyc_services_2_year$psychotherapy_medication, levels = c(1,0), labels = c("Yes", "No"))


# Rename and recode ethnicity 
ethnicity <- ethnicity %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, race_ethnicity, rel_family_id) %>%
  mutate(race_ethnicity_cat = ifelse(race_ethnicity == "1", 1, 2)) 

# Set data type and level labels
ethnicity$race_ethnicity <- factor(ethnicity$race_ethnicity,
                                   levels = c(1, 2, 3, 4, 5),
                                   labels = c("White", "Black", "Hispanic", "Asian", "Other"))

# Rename and recode household income
parent_demographics <- parent_demographics %>%
  rename(household_income = demo_comb_income_v2) %>%
  mutate(household_income = recode(household_income,
                                   "1" = "1",  # [<50K]
                                   "2" = "1",  # [<50K]
                                   "3" = "1",  # [<50K]
                                   "4" = "1",  # [<50K]
                                   "5" = "1",  # [<50K]
                                   "6" = "1",  # [<50K]
                                   "7" = "2",  # [>=50K & <100K]
                                   "8" = "2",  # [>=50K & <100K]
                                   "9" = "3",  # [>=100K]
                                   "10" ="3",  # [>=100K]
                                   "777" = "777",  # To be NA (replaced below) 
                                   "999" = "999")) # To be NA (replaced below) 

parent_demographics <- parent_demographics %>%
  mutate(household_income = na_if(household_income, "777")) %>%
  mutate(household_income = na_if(household_income, "999"))

parent_demographics <- parent_demographics %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, household_income)

# Set data type and level labels
parent_demographics$household_income <- factor(parent_demographics$household_income,
                                               levels= 1:3,
                                               labels = c("[<50K]", "[>=50K & <100K]", "[>=100K]"))

# Recode interview age to be measured in years not months 
# Interview age at baseline 
tbi_baseline <- tbi_baseline %>% 
  mutate(interview_age = interview_age / 12)

# Interview age at 2-year follow-up 
tbi_2_year <- tbi_2_year %>% 
  mutate(interview_age = interview_age / 12)

# Recode sex 
# Set data type and level labels for baseline data
tbi_baseline$sex <- factor(tbi_baseline$sex,
                           levels = c("F", "M"), 
                           labels = c("Female", "Male"))

# Set data type and level labels for 2-year follow-up data
tbi_2_year$sex <- factor(tbi_2_year$sex,
                         levels = c("F", "M"), 
                         labels = c("Female", "Male"))

# Rename parent mental health 
parent_mh <- parent_mh %>%
  rename(parent_depression = asr_scr_depress_r) %>%
  rename(parent_anxiety = asr_scr_anxdisord_r) %>%  
  rename(parent_adhd = asr_scr_adhd_r) %>%
  rename(parent_antisocial = asr_scr_antisocial_r) 

# Sum parent mental health into 'number of parent mental health disorders' 
parent_mh <- parent_mh %>%
  mutate(parent_mh_score =
           parent_depression + 
           parent_anxiety + 
           parent_adhd + 
           parent_antisocial) 

# Select parent mental health score at baseline
parent_mh <- parent_mh %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(subjectkey, parent_mh_score)  

# Select study site ID
site <- site %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, site_id_l)

# Create adverse life event variable (0 = absence; 1 = presence)
advlife <- advlife %>%
  mutate(advlife_total = 
           ksads_ptsd_raw_754_p + 
           ksads_ptsd_raw_755_p + 
           ksads_ptsd_raw_756_p +
           ksads_ptsd_raw_757_p +
           ksads_ptsd_raw_758_p +
           ksads_ptsd_raw_759_p +
           ksads_ptsd_raw_760_p +
           ksads_ptsd_raw_761_p +
           ksads_ptsd_raw_762_p +
           ksads_ptsd_raw_763_p +
           ksads_ptsd_raw_764_p +
           ksads_ptsd_raw_765_p +
           ksads_ptsd_raw_766_p +
           ksads_ptsd_raw_767_p +
           ksads_ptsd_raw_768_p +
           ksads_ptsd_raw_769_p +
           ksads_ptsd_raw_770_p) %>%
  mutate(advlife_cat = ifelse(advlife_total >=1, 1, 0)) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, advlife_cat)

# Set data type and level labels
advlife$advlife_cat <- factor(advlife$advlife_cat,
                              levels = c(0,1),
                              labels = c("No", "Yes"))

# Create neighbourhood safety score 
neighbourhood_safety_baseline <- neighbourhood_safety %>%
  rename(safety = neighborhood_crime_y) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, safety)

# Set data type and level labels
neighbourhood_safety_baseline$safety <- factor(neighbourhood_safety_baseline$safety,
                                               levels = c(1,2,3,4,5),
                                               labels = c("Strongly Disagree", "Disagree", "Neutral", "Agree", "Strongly Agree"))

# Create family conflict scale
family_conflict_baseline <- family_conflict %>%
  mutate(family_conflict_scale = 
           fes_youth_q1 + 
           fes_youth_q3 + 
           fes_youth_q6) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, family_conflict_scale)


###########################################################
#
# Merge datasets to create baseline and 2-year follow-up data
#
###########################################################

# Merge baseline files
psyc_data <- tbi_baseline %>%
  left_join(oi_baseline, by = "subjectkey") %>%
  left_join(psyc_services_baseline, by = "subjectkey") %>%
  left_join(ethnicity, by = "subjectkey") %>%
  left_join(parent_demographics, by = "subjectkey") %>%
  left_join(site, by = "subjectkey") %>%
  left_join(advlife, by = "subjectkey") %>%
  left_join(neighbourhood_safety_baseline, by = "subjectkey") %>%
  left_join(family_conflict_baseline, by = "subjectkey") %>%
  left_join(parent_mh, by = "subjectkey") %>%
  unique()

# Merge 2-year follow-up files
psyc_data_t2 <- tbi_1_year %>%
  left_join(tbi_baseline, by = "subjectkey") %>%
  left_join(tbi_2_year, by = "subjectkey") %>%
  left_join(oi_baseline, by = "subjectkey") %>%
  left_join(oi_1_year, by = "subjectkey") %>%
  left_join(oi_2_year, by = "subjectkey") %>%
  left_join(psyc_services_2_year, by = "subjectkey") %>%
  left_join(ethnicity, by = "subjectkey") %>%
  left_join(parent_mh, by = "subjectkey") %>% 
  left_join(parent_demographics, by = "subjectkey") %>%
  left_join(site, by = "subjectkey") %>%
  left_join(advlife, by = "subjectkey") %>%
  left_join(neighbourhood_safety_baseline, by = "subjectkey") %>%
  left_join(family_conflict_baseline, by = "subjectkey") %>%
  unique()


##############################################################################
#
# Remove baseline interview age and sex variables for 2-year follow-up data
#
#############################################################################

# Remove age and sex duplicates
psyc_data_t2 <- subset(psyc_data_t2, select = -c(interview_age.x, sex.y))

# Rename interview age and sex variables at 2-year follow-up 
psyc_data_t2 <- psyc_data_t2 %>%
  rename(interview_age = interview_age.y) %>%
  rename(sex = sex.x)  

##############################################################################
#
# Create injury type variable consisting of TBI, ortho injury and no injury cases
#
#############################################################################

# Create injury type variable at baseline using TBI severity and ortho injury columns 
psyc_data$injury_type <- paste(psyc_data$tbi_severity, psyc_data$broken_bones)

# Recode injury type column 
psyc_data$injury_type[psyc_data$injury_type == 'No No'] <- 0
psyc_data$injury_type[psyc_data$injury_type == 'Yes No'] <- 2
psyc_data$injury_type[psyc_data$injury_type == 'No Yes'] <- 1
psyc_data$injury_type[psyc_data$injury_type == 'Yes Yes'] <- 2

# Set data type and level labels
psyc_data$injury_type <- factor(psyc_data$injury_type, levels = c(2,1,0), labels = c("TBI", "Ortho", "None" ))

# Combine TBI cases from 1- and 2-year follow-up 
#
# Create column merging TBI at 1 and 2-year follow-up 
psyc_data_t2$tbi_both_years <- paste(psyc_data_t2$tbi_severity_1_year, psyc_data_t2$tbi_severity_2_year)

# Recode column 
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No No'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No NA'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes NA'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'NA Yes'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes No'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No Yes'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes Yes'] <- 1

# Set data type and level labels
psyc_data_t2$tbi_both_years <- factor(psyc_data_t2$tbi_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Create column merging TBI at baseline, 1 and 2 year follow-up 
psyc_data_t2$tbi_both_years <- paste(psyc_data_t2$tbi_severity, psyc_data_t2$tbi_both_years)

# Recode column so that TBI at baseline is removed (only new cases included)
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No No'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No NA'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'NA No'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes NA'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'NA Yes'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes No'] <- 0
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'No Yes'] <- 1
psyc_data_t2$tbi_both_years[psyc_data_t2$tbi_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
psyc_data_t2$tbi_both_years <- factor(psyc_data_t2$tbi_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Check if variabe is coded correctly 
xtabs(~ tbi_both_years, data = psyc_data_t2)


# Create column merging broken bones at 1 and 2-year follow-up 
psyc_data_t2$broken_bones_both_years <- paste(psyc_data_t2$broken_bones_1_year, psyc_data_t2$broken_bones_2_year)

# Recode column 
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No No'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'NA No'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No NA'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes NA'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'NA Yes'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes No'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No Yes'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes Yes'] <- 1

# Set data type and level labels
psyc_data_t2$broken_bones_both_years <- factor(psyc_data_t2$broken_bones_both_years, levels = c(1,0), labels = c("Yes", "No" ))


# Create column merging OI at baseline, 1-year and 2-year follow-up 
psyc_data_t2$broken_bones_both_years <- paste(psyc_data_t2$broken_bones, psyc_data_t2$broken_bones_both_years)

# Recode column so that OI at baseline is not included (only new cases included)
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No No'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'NA No'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No NA'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes NA'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'NA Yes'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes No'] <- 0
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'No Yes'] <- 1
psyc_data_t2$broken_bones_both_years[psyc_data_t2$broken_bones_both_years == 'Yes Yes'] <- NA

# Set data type and level labels
psyc_data_t2$broken_bones_both_years <- factor(psyc_data_t2$broken_bones_both_years, levels = c(1,0), labels = c("Yes", "No" ))

# Create injury type variable at 2-year follow-up using TBI severity and ortho injury columns
psyc_data_t2$injury_type <- paste(psyc_data_t2$tbi_both_years, psyc_data_t2$broken_bones_both_years)

# Recode injury type column 
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'No No'] <- 0
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'NA No'] <- 0
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'No NA'] <- 0
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'Yes No'] <- 2
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'No Yes'] <- 1
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'Yes Yes'] <- 2
psyc_data_t2$injury_type[psyc_data_t2$injury_type == 'Yes NA'] <- 2

# Set data type and level labels
psyc_data_t2$injury_type <- factor(psyc_data_t2$injury_type, levels = c(2, 1,0), labels = c("TBI", "Ortho", "None" ))


##############################################################################
#
# Rename variables for better readability 
# Create tables for demographic and descriptive statistics of sample 
#
#############################################################################

# Rename baseline variables for better readability
label(psyc_data$interview_age) <- "Age"
label(psyc_data$sex) <- "Sex"
label(psyc_data$tbi_severity) <- "TBI"
label(psyc_data$race_ethnicity) <- "Ethnicity"
label(psyc_data$household_income) <- "Combined annual household income"
label(psyc_data$safety) <- "Neighbourhood safety"
label(psyc_data$family_conflict_scale) <- "Family conflict scale"
label(psyc_data$advlife_cat) <- "Traumatic events"
label(psyc_data$broken_bones) <- "Orthopaedic injury"
label(psyc_data$parent_mh_score) <- "Parent mental health score"
label(psyc_data$any_support) <- "Any mental health/sustance abuse services"
label(psyc_data$outpatient) <- "Outpatient for mental health"
label(psyc_data$inpatient) <- "Inpatient for mental health"
label(psyc_data$psychotherapy) <- "Psychotherapy"
label(psyc_data$medication) <- "Medication management"
label(psyc_data$psychotherapy_medication) <- "Psychotherapy or medication management"

# Create demographic and decriptive table for sample 
descriptive_table <- psyc_data %>%
  select(injury_type, sex, interview_age, race_ethnicity, household_income, safety, family_conflict_scale, advlife_cat) %>%
  tbl_summary(by = injury_type, type = list(family_conflict_scale ~ 'categorical', advlife_cat ~ 'categorical')) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Injury type**") %>%
  modify_caption("**Table 1. Demographics and descriptive statistics for the sample**")

# Save table
table_1_filename = paste(output_dir, "descriptive_table.html", sep="")
gt::gtsave(as_gt(descriptive_table), file = table_1_filename)
display_html(file=table_1_filename)
head(table_1_filename)


# Create table showing percentages of service use at baseline  
psyc_service_percentage_table <- psyc_data %>%
  select(injury_type, any_support, outpatient, inpatient, psychotherapy, medication) %>%
  tbl_summary(by = injury_type) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Injury type**") %>%
  modify_caption("**Table 1. Proportions and percentages of psychiatric service use at baseline**")

# Save table
table_1_filename = paste(output_dir, "psyc_use_percentage_baseline.html", sep="")
gt::gtsave(as_gt(psyc_service_percentage_table), file = table_1_filename)
display_html(file=table_1_filename)
head(table_1_filename)


# Rename 2-year follow-up variables for better readability 
label(psyc_data_t2$interview_age) <- "Age"
label(psyc_data_t2$sex) <- "Sex"
label(psyc_data_t2$household_income) <- "Household Income"
label(psyc_data_t2$race_ethnicity) <- "Ethnicity"
label(psyc_data_t2$any_support) <- "Any mental health/sustance abuse services"
label(psyc_data_t2$outpatient) <- "Outpatient for mental health"
label(psyc_data_t2$inpatient) <- "Inpatient for mental helth"
label(psyc_data_t2$psychotherapy) <- "Psychotherapy"
label(psyc_data_t2$medication) <- "Medication management"
label(psyc_data_t2$advlife_cat) <- "Traumatic events"
label(psyc_data_t2$safety) <- "Neighbourhood safety"
label(psyc_data_t2$family_conflict_scale) <- "Family conflict scale"
label(psyc_data_t2$psychotherapy_medication) <- "Psychotherapy or medication management"



# Create table showing percentages of service use at 2-year follow-up  
psyc_service_percentage_table_t2 <- psyc_data_t2 %>%
  select(injury_type, any_support, outpatient, inpatient, psychotherapy, medication) %>%
  tbl_summary(by = injury_type) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Injury type**") %>%
  modify_caption("**Table 1. Proportions and percentages of psychiatric service use at two-year follow-up**")

# Save table
table_1_filename = paste(output_dir, "psyc_use_percentage_t2.html", sep="")
gt::gtsave(as_gt(psyc_service_percentage_table_t2), file = table_1_filename)
display_html(file=table_1_filename)
head(table_1_filename)

##############################################################################
#
# Multiple imputation using missForest package
#
#############################################################################

# Impute baseline dataset
#

# Register process with doParallel 
registerDoParallel()

# Set seed for reproducibility 
registerDoRNG(seed = 123)

# Ensure study site ID is stored as a factor  
class(psyc_data$site_id_l)
psyc_data$site_id_l <- as.factor(psyc_data$site_id_l)

# Select variables for imputed dataset 
imputed_data <- select(psyc_data, 
                       subjectkey,
                       interview_age, sex, race_ethnicity, parent_mh_score, household_income, advlife_cat, 
                       safety, family_conflict_scale, tbi_severity, broken_bones, site_id_l, injury_type,
                       any_support, outpatient, inpatient, medication, psychotherapy, psychotherapy_medication)

# Store subject key
id_sk <- select(imputed_data, subjectkey)

# Remove it from imputed_data to not impute it
imputed_data$subjectkey <- NULL 

# Create imputed dataframe 
imputed_data <- as.data.frame(imputed_data)

# Visualise missing data 
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 200)
vis_miss(imputed_data, sort_miss = TRUE, warn_large_data = FALSE)
ggsave(paste(output_dir, "missing_psych_service_1.png", sep = ""))


set.seed(123)

# Record start time 
start_time <- Sys.time() 

# Impute missing data
imputed_data_baseline <- missForest(imputed_data, parallelize = "forests")

# Record end time 
end_time <- Sys.time() 

# Report time duration of missForest imputation 
end_time - start_time 

# Check normalised root mean squared error 
imputed_data_baseline$OOBerror 

# Create data frame with imputed covariates and outcomes 
imputed_data_baseline <- imputed_data_baseline$ximp 

# Restore subject key
imputed_data_baseline$subject_key <- id_sk$subjectkey

# Save new file for analyses 
write.csv(imputed_data_baseline, paste(output_dir, "imputed_baseline_psyc_services.csv", sep = ""))



# Impute 2-year follow-up dataset
#

# Register process with doParallel 
registerDoParallel()

# Set seed for reproducibility 
registerDoRNG(seed = 123)

# Ensure study site ID is stored as a factor  
class(psyc_data_t2$site_id_l)
psyc_data_t2$site_id_l <- as.factor (psyc_data_t2$site_id_l)

# Select variables for imputed dataset 
imputed_data_2 <- select(psyc_data_t2, 
                       subjectkey,
                       interview_age, sex, race_ethnicity, parent_mh_score, household_income, advlife_cat, 
                       safety, family_conflict_scale, injury_type, site_id_l, tbi_both_years, 
                       broken_bones_both_years,
                       any_support, outpatient, inpatient, medication, psychotherapy, psychotherapy_medication)

# Store subject key
id_sk <- select(imputed_data_2, subjectkey)

# Make sure subject ID is not imputed 
imputed_data_2$subjectkey <- NULL 

# Create imputed dataframe 
imputed_data_2 <- as.data.frame(imputed_data_2)

# Visualise missing data 
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 200)
vis_miss(imputed_data_2, sort_miss = TRUE, warn_large_data = FALSE)
ggsave(paste(output_dir, "missing_psych_service_2.png", sep = ""))

set.seed(123)

# Record start time 
start_time <- Sys.time() 

# Impute missing data 
imputed_data_t2 <- missForest(imputed_data_2, parallelize = "forests")

# Record end time 
end_time <- Sys.time() 

# Report time duration of missForest imputation 
end_time - start_time 

# Check normalized root mean squared error 
imputed_data_t2$OOBerror 

# Create data frame with imputed covariates and outcomes  
imputed_data_t2 <- imputed_data_t2$ximp 

# Restore subject key
imputed_data_t2$subject_key <- id_sk$subjectkey

# Save new file for analyses 
write.csv(imputed_data_t2, paste(output_dir, "imputed_t2_psyc_services.csv", sep = ""))

##############################################################################
#
# Load imputed datasets to avoid running imputation again 
#
##############################################################################

# Load imputed dataset for baseline 
imputed_data_baseline <- read.csv(paste(output_dir, "imputed_baseline_psyc_services.csv", sep = ""))

# Load imputed dataset for 2-year follow-up 
imputed_data_t2 <- read.csv(paste(output_dir, "imputed_t2_psyc_services.csv", sep = ""))

##############################################################################
#
# Psychiatric service use 
#
# Adjusted model 1 covariates are age, sex, ethnicity and household income
# Adjusted model 2 covariates are adverse life event, neighbourhood safety,
#  parent mental health and family conflict 
# 
#############################################################################

# Baseline
# Ensure variables are stored as factors 
imputed_data_baseline$any_support <- as.factor(imputed_data_baseline$any_support)
imputed_data_baseline$outpatient <- as.factor(imputed_data_baseline$outpatient)
imputed_data_baseline$inpatient <- as.factor(imputed_data_baseline$inpatient)
imputed_data_baseline$psychotherapy <- as.factor(imputed_data_baseline$psychotherapy)
imputed_data_baseline$medication <- as.factor(imputed_data_baseline$medication)
imputed_data_baseline$injury_type <- as.factor(imputed_data_baseline$injury_type)
imputed_data_baseline$psychotherapy_medication <- as.factor(imputed_data_baseline$psychotherapy_medication)

# Any mental health support measured at baseline
#

# Change reference category 
imputed_data_baseline$any_support <- relevel(imputed_data_baseline$any_support, ref = "No")
imputed_data_baseline$injury_type <- relevel(imputed_data_baseline$injury_type, ref = "None")

# Run unadjusted model 
any_support_model_unadjusted <- glmer(any_support ~ injury_type +
                                        (1|site_id_l),
                                      data = imputed_data_baseline, family = binomial(),
                                      control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
any_support_model_adjusted1 <- glmer(any_support ~ injury_type +
                                       interview_age +
                                       sex +
                                       race_ethnicity +
                                       household_income +
                                       (1|site_id_l),
                                     data = imputed_data_baseline,
                                     family = binomial(link = "logit"),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)


# Run adjusted model 2
any_support_model_adjusted2 <- glmer(any_support ~ injury_type +
                                       interview_age +
                                       sex +
                                       race_ethnicity +
                                       household_income +
                                       advlife_cat +
                                       safety +
                                       family_conflict_scale +
                                       parent_mh_score +
                                       (1|site_id_l),
                                     data = imputed_data_baseline,
                                     family = binomial(link = "logit"),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Psychotherapy measured at baseline

# Change reference category
imputed_data_baseline$psychotherapy <- relevel(imputed_data_baseline$psychotherapy, ref = "No")

# Run unadjusted model 
imputed_data_baseline$psychotherapy <- as.factor(imputed_data_baseline$psychotherapy)


psychotherapy_model_unadjusted <- glmer(psychotherapy ~ injury_type + 
                                          (1|site_id_l),
                                        data = imputed_data_baseline,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_model_adjusted1 <- glmer(psychotherapy ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income + 
                                         (1|site_id_l),
                                       data = imputed_data_baseline,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)



# Run adjusted model 2
psychotherapy_model_adjusted2 <- glmer(psychotherapy ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income + 
                                         advlife_cat +
                                         safety +
                                         family_conflict_scale +
                                         parent_mh_score +
                                         (1|site_id_l),
                                       data = imputed_data_baseline,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Medication for mental health measured at baseline
#

# Change reference category
imputed_data_baseline$medication <- relevel(imputed_data_baseline$medication, ref = "No")

# Run unadjusted model 
medication_model_unadjusted <- glmer(medication ~ injury_type + 
                                       (1|site_id_l),
                                     data = imputed_data_baseline,
                                     family = binomial(link = "logit"),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
medication_model_adjusted1 <- glmer(medication ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income + 
                                      (1|site_id_l),
                                    data = imputed_data_baseline, family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
medication_model_adjusted2 <- glmer(medication ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income + 
                                      advlife_cat +
                                      safety +
                                      family_conflict_scale +
                                      parent_mh_score +
                                      (1|site_id_l),
                                    data = imputed_data_baseline,
                                    family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)


# Psychotherapy and mmedication for mental health measured at baseline
#

# Change reference category
imputed_data_baseline$psychotherapy_medication <- relevel(imputed_data_baseline$psychotherapy_medication, ref = "No")

# Run unadjusted model 
psychotherapy_medication_model_unadjusted <- glmer(psychotherapy_medication ~ injury_type + 
                                       (1|site_id_l),
                                     data = imputed_data_baseline,
                                     family = binomial(link = "logit"),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_medication_model_adjusted1 <- glmer(psychotherapy_medication ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income + 
                                      (1|site_id_l),
                                    data = imputed_data_baseline, family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
psychotherapy_medication_model_adjusted2 <- glmer(psychotherapy_medication ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income + 
                                      advlife_cat +
                                      safety +
                                      family_conflict_scale +
                                      parent_mh_score +
                                      (1|site_id_l),
                                    data = imputed_data_baseline,
                                    family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)


# Outpatient support measured at baseline
#

# Change reference category
imputed_data_baseline$outpatient <- relevel(imputed_data_baseline$outpatient, ref = "No")

# Run unadjusted model 
outpatient_model_unadjusted <- glmer(outpatient ~ injury_type +
                                        (1|site_id_l),
                                     data = imputed_data_baseline,
                                     family = binomial(),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_unadjusted)

# Run adjusted model 1
outpatient_model_adjusted1 <- glmer(outpatient ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income +
                                      (1|site_id_l),
                                     data = imputed_data_baseline,
                                    family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_adjusted1)

# Run adjusted model 2
outpatient_model_adjusted2 <- glmer(outpatient ~ injury_type +
                                      interview_age +
                                      sex +
                                      race_ethnicity +
                                      household_income +
                                      advlife_cat +
                                      safety +
                                      family_conflict_scale +
                                      parent_mh_score +
                                      (1|site_id_l),
                                    data = imputed_data_baseline, family = binomial(link = "logit"),
                                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_adjusted2)

# Inpatient support measured at baseline
#

# Change reference category
imputed_data_baseline$inpatient <- relevel(imputed_data_baseline$inpatient, ref = "No")

# Run unadjusted model 
inpatient_model_unadjusted <- glmer(inpatient ~ injury_type +
                                       (1|site_id_l),
                                     data = imputed_data_baseline, family = binomial(),
                                     control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_unadjusted)

# Run adjusted model 1
inpatient_model_adjusted1 <- glmer(inpatient ~ injury_type +
                                     interview_age +
                                     sex +
                                     race_ethnicity +
                                     household_income +
                                     (1|site_id_l),
                                   data = imputed_data_baseline, family = binomial(link = "logit"),
                                   control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_adjusted1)

# Run adjusted model 2
inpatient_model_adjusted2 <- glmer(inpatient ~ injury_type +
                                     interview_age +
                                     sex +
                                     race_ethnicity +
                                     household_income +
                                     advlife_cat +
                                     safety +
                                     family_conflict_scale +
                                     parent_mh_score +
                                     (1|site_id_l),
                                   data = imputed_data_baseline, family = binomial(link = "logit"),
                                   control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_adjusted2)

#
# 2-year follow-up
#

# Ensure variables are stored as factors 
imputed_data_t2$any_support <- as.factor(imputed_data_t2$any_support)
imputed_data_t2$outpatient <- as.factor(imputed_data_t2$outpatient)
imputed_data_t2$inpatient <- as.factor(imputed_data_t2$inpatient)
imputed_data_t2$psychotherapy <- as.factor(imputed_data_t2$psychotherapy)
imputed_data_t2$medication <- as.factor(imputed_data_t2$medication)
imputed_data_t2$injury_type <- as.factor(imputed_data_t2$injury_type)
imputed_data_t2$psychotherapy_medication <- as.factor(imputed_data_t2$psychotherapy_medication)

# Any mental health support measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$any_support <- relevel(imputed_data_t2$any_support, ref = "No")
imputed_data_t2$injury_type <- relevel(imputed_data_t2$injury_type, ref = "None")

# Run unadjusted model 
any_support_model_unadjusted_t2 <- glmer(any_support ~ injury_type +
                                           (1|site_id_l),
                                         data = imputed_data_t2,
                                         family = binomial(),
                                         control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
any_support_model_adjusted1_t2 <- glmer(any_support ~ injury_type +
                                          interview_age +
                                          sex +
                                          race_ethnicity +
                                          household_income +
                                          (1|site_id_l),
                                        data = imputed_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
any_support_model_adjusted2_t2 <- glmer(any_support ~ injury_type +
                                          advlife_cat +
                                          safety +
                                          family_conflict_scale +
                                          parent_mh_score +
                                          (1|site_id_l),
                                        data = imputed_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Psychotherapy measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$psychotherapy <- relevel(imputed_data_t2$psychotherapy, ref = "No")

# Run unadjusted model 
psychotherapy_model_unadjusted_t2 <- glmer(psychotherapy ~ injury_type + 
                                             (1|site_id_l),
                                           data = imputed_data_t2,
                                           family = binomial(link = "logit"),
                                           control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_model_adjusted1_t2 <- glmer(psychotherapy ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income + 
                                            (1|site_id_l),
                                          data = imputed_data_t2,
                                          family = binomial(link = "logit"),
                                          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
psychotherapy_model_adjusted2_t2 <- glmer(psychotherapy ~ injury_type +
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
                                          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Medication for mental health measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$medication <- relevel(imputed_data_t2$medication, ref = "No")

# Run unadjusted model 
medication_model_unadjusted_t2 <- glmer(medication ~ injury_type + 
                                          (1|site_id_l),
                                        data = imputed_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
medication_model_adjusted1_t2 <- glmer(medication ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income + 
                                         (1|site_id_l),
                                       data = imputed_data_t2,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
medication_model_adjusted2_t2 <- glmer(medication ~ injury_type +
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
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Psychotherapy and medication for mental health measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$psychotherapy_medication <- relevel(imputed_data_t2$psychotherapy_medication, ref = "No")

# Run unadjusted model 
psychotherapy_medication_model_unadjusted_t2 <- glmer(psychotherapy_medication ~ injury_type + 
                                          (1|site_id_l),
                                        data = imputed_data_t2,
                                        family = binomial(link = "logit"),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
psychotherapy_medication_model_adjusted1_t2 <- glmer(psychotherapy_medication ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income + 
                                         (1|site_id_l),
                                       data = imputed_data_t2,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 2
psychotherapy_medication_model_adjusted2_t2 <- glmer(psychotherapy_medication ~ injury_type +
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
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Outpatient support measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$outpatient <- relevel(imputed_data_t2$outpatient, ref = "No")

# Run unadjusted model 
outpatient_model_unadjusted_t2 <- glmer(outpatient ~ injury_type +
                                          (1|site_id_l),
                                        data = imputed_data_t2,
                                        family = binomial(),
                                        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_unadjusted_t2)

# Run adjusted model 1
outpatient_model_adjusted1_t2 <- glmer(outpatient ~ injury_type +
                                         interview_age +
                                         sex +
                                         race_ethnicity +
                                         household_income +
                                         (1|site_id_l),
                                       data = imputed_data_t2,
                                       family = binomial(link = "logit"),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_adjusted1_t2)

# Run adjusted model 2
outpatient_model_adjusted2_t2 <- glmer(outpatient ~ injury_type +
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
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(outpatient_model_adjusted2_t2)

# Inpatient support measured at 2-year follow-up 
#

# Change reference category 
imputed_data_t2$inpatient <- relevel(imputed_data_t2$inpatient, ref = "No")



# Run unadjusted model 
inpatient_model_unadjusted_t2 <- glmer(inpatient ~ injury_type +
                                         (1|site_id_l),
                                       data = imputed_data_t2,
                                       family = binomial(),
                                       control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_unadjusted_t2)

# Run adjusted model 1
inpatient_model_adjusted1_t2 <- glmer(inpatient ~ injury_type +
                                        interview_age +
                                        sex +
                                        race_ethnicity +
                                        household_income +
                                        (1|site_id_l),
                                      data = imputed_data_t2,
                                      family = binomial(link = "logit"),
                                      control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_adjusted1_t2)

# Run adjusted model 2
inpatient_model_adjusted2_t2 <- glmer(inpatient ~ injury_type +
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
                                      control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# View model values
summary(inpatient_model_adjusted2_t2)


##############################################################################
#
# Create tables for psychiatric service use measured at baseline 
#
#############################################################################

# Any mental health support measured at baseline
#

# Create unadjusted model table
any_support_unadjusted_table <- tbl_regression(any_support_model_unadjusted, exponentiate = TRUE, label = injury_type ~ "Any mental health service use") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
any_support_adjusted_table_short_model1 <- tbl_regression(any_support_model_adjusted1, exponentiate = TRUE, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Any mental health service use") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
any_support_adjusted_table_short_model2 <- tbl_regression(any_support_model_adjusted2, exponentiate = TRUE, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Any mental health service use") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_any_support_short <- tbl_merge(tbls = list(any_support_unadjusted_table, any_support_adjusted_table_short_model1, any_support_adjusted_table_short_model2),
                                           tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Psychotherapy measured at baseline
#

# Create unadjusted model table
psychotherapy_unadjusted_table <- tbl_regression(psychotherapy_model_unadjusted, exponentiate = TRUE, label = injury_type ~ "Psychotherapy") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
psychotherapy_adjusted_table_short_model1 <- tbl_regression(psychotherapy_model_adjusted1, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Psychotherapy") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
psychotherapy_adjusted_table_short_model2 <- tbl_regression(psychotherapy_model_adjusted2, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Psychotherapy") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_psychotherapy_short <- tbl_merge(tbls = list(psychotherapy_unadjusted_table, psychotherapy_adjusted_table_short_model1, psychotherapy_adjusted_table_short_model2),
                                             tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Medication for mental health measured at baseline
#

# Create unadjusted model table
medication_unadjusted_table <- tbl_regression(medication_model_unadjusted, exponentiate = TRUE, 
                                              label = injury_type ~ "Medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
medication_adjusted_table_short_model1 <- tbl_regression(medication_model_adjusted1, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
medication_adjusted_table_short_model2 <- tbl_regression(medication_model_adjusted2, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Medication for mental health") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_medication_short <- tbl_merge(tbls = list(medication_unadjusted_table, medication_adjusted_table_short_model1, medication_adjusted_table_short_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Psychotherapy and medication for mental health measured at baseline
#

# Create unadjusted model table
psychotherapy_medication_unadjusted_table <- tbl_regression(psychotherapy_medication_model_unadjusted, exponentiate = TRUE, 
                                              label = injury_type ~ "Psychotherapy and/or medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
psychotherapy_medication_adjusted_table_short_model1 <- tbl_regression(psychotherapy_medication_model_adjusted1, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Psychotherapy and/or medication for mental health") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
psychotherapy_medication_adjusted_table_short_model2 <- tbl_regression(psychotherapy_medication_model_adjusted2, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Psychotherapy and/or medication for mental health") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_psychotherapy_medication_short <- tbl_merge(tbls = list(psychotherapy_medication_unadjusted_table, 
                                                                    psychotherapy_medication_adjusted_table_short_model1, psychotherapy_medication_adjusted_table_short_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Create table for psychiatric service use at baseline
#

# Merge  psychotherapy, medication and any mental health support tables 
psyc_service_baseline_table <- tbl_stack(tbls = list(table_merge_psychotherapy_short, table_merge_medication_short, table_merge_psychotherapy_medication_short, table_merge_any_support_short)) %>%
  modify_caption("**Table 4. Association between history of mTBI / orthopaedic injury and mental health service use at baseline (age 9-10)**")

# Save table 
table_2_filename = paste(output_dir, "imputed_psyc_service_baseline.html", sep="")
gt::gtsave(as_gt(psyc_service_baseline_table), file = table_2_filename)
head(table_2_filename)


##############################################################################
#
# Create tables for psychiatric service use measured at 2-year follow-up 
#
#############################################################################

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
table_2_filename = paste(output_dir, "imputed_psyc_service_t2.html", sep="")
gt::gtsave(as_gt(psyc_service_t2_table), file = table_2_filename)
head(table_2_filename)


# Create outpatient and inpatient tables 

# Outpatient measured at baseline
#

# Create unadjusted model table
outpatient_unadjusted_table <- tbl_regression(outpatient_model_unadjusted, exponentiate = TRUE, 
                                              label = injury_type ~ "Outpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
outpatient_adjusted_table_short_model1 <- tbl_regression(outpatient_model_adjusted1, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Outpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
outpatient_adjusted_table_short_model2 <- tbl_regression(outpatient_model_adjusted2, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Outpatient") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_outpatient_short <- tbl_merge(tbls = list(outpatient_unadjusted_table, outpatient_adjusted_table_short_model1, outpatient_adjusted_table_short_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Create unadjusted model table
inpatient_unadjusted_table <- tbl_regression(inpatient_model_unadjusted, exponentiate = TRUE, 
                                              label = injury_type ~ "Inpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
inpatient_adjusted_table_short_model1 <- tbl_regression(inpatient_model_adjusted1, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Inpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
inpatient_adjusted_table_short_model2 <- tbl_regression(inpatient_model_adjusted2, exponentiate = TRUE, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type,
                                                         label = injury_type ~ "Inpatient") %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_inpatient_short <- tbl_merge(tbls = list(inpatient_unadjusted_table, inpatient_adjusted_table_short_model1, inpatient_adjusted_table_short_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# 2-year follow-up 

# Outpatient
# Create unadjusted model table
outpatient_unadjusted_table_t2 <- tbl_regression(outpatient_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Outpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
outpatient_adjusted_table_short_model1_t2 <- tbl_regression(outpatient_model_adjusted1_t2, exponentiate = TRUE, 
                                                             pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Outpatient") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
outpatient_adjusted_table_short_model2_t2 <- tbl_regression(outpatient_model_adjusted2_t2, exponentiate = TRUE, 
                                                             pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Outpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_outpatient_short_t2 <- tbl_merge(tbls = list(outpatient_unadjusted_table_t2, outpatient_adjusted_table_short_model1_t2, outpatient_adjusted_table_short_model2_t2),
                                              tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Inpatient
# Create unadjusted model table
inpatient_unadjusted_table_t2 <- tbl_regression(inpatient_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Inpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
inpatient_adjusted_table_short_model1_t2 <- tbl_regression(inpatient_model_adjusted1_t2, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Inpatient") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
inpatient_adjusted_table_short_model2_t2 <- tbl_regression(inpatient_model_adjusted2_t2, exponentiate = TRUE, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Inpatient") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_inpatient_short_t2 <- tbl_merge(tbls = list(inpatient_unadjusted_table_t2, inpatient_adjusted_table_short_model1_t2, inpatient_adjusted_table_short_model2_t2),
                                             tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Create table for outpatient and inpatient at baseline
#

# Merge  outpatient and inpatient tables 
outpatient_inpatient_baseline_table <- tbl_stack(tbls = list(table_merge_outpatient_short, table_merge_inpatient_short)) %>%
  modify_caption("**Table 4. Association between history of mTBI / orthopaedic injury and mental health service use at baseline (age 9-10)**")

# Save table 
table_2_filename = paste(output_dir, "oupatient_inpatient_baseline.html", sep="")
gt::gtsave(as_gt(outpatient_inpatient_baseline_table), file = table_2_filename)
head(table_2_filename)


# # Create table for outpatient and inpatient at 2-year follow-up
#

# Merge  outpatient and inpatient tables 
outpatient_inpatient_t2_table <- tbl_stack(tbls = list(table_merge_outpatient_short_t2, table_merge_inpatient_short_t2)) %>%
  modify_caption("**Table 5. Association between new mTBI / orthopaedic injury and mental health service use in the 12-24 months following baseline (age 11-12)**")

# Save table 
table_2_filename = paste(output_dir, "outpatient_inpatient_t2.html", sep="")
gt::gtsave(as_gt(outpatient_inpatient_t2_table), file = table_2_filename)
head(table_2_filename)
