################################################################
#
# Analysis of Child Behaviour Checklist (CBCL) and DSM diagnoses
# outcomes (from KSADS) from the study:
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
rm(list = ls())

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

# CBCL
parent_cbcl <- read_abcd_file("abcd_cbcls01.txt")

# KSAD 
ksad_parent <- read_abcd_file("abcd_ksad01.txt")

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

# Rename CBCL variables 
parent_cbcl <- parent_cbcl %>%
  rename(cbcl_tot_problems = cbcl_scr_syn_totprob_r) %>%
  rename(cbcl_int_problems = cbcl_scr_syn_internal_r) %>%
  rename(cbcl_ext_problems = cbcl_scr_syn_external_r) %>%
  rename(cbcl_dep_score = cbcl_scr_dsm5_depress_r) %>%
  rename(cbcl_anx_score = cbcl_scr_dsm5_anxdisord_r) %>%
  rename(cbcl_somatic_score = cbcl_scr_dsm5_somaticpr_r) %>%
  rename(cbcl_adhd_score = cbcl_scr_dsm5_adhd_r) %>%
  rename(cbcl_odd_score = cbcl_scr_dsm5_opposit_r) %>%
  rename(cbcl_conduct_score = cbcl_scr_dsm5_conduct_r) 

# Select CBCL variables at baseline
cbcl_t0 <- parent_cbcl %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(subjectkey, cbcl_tot_problems, cbcl_int_problems, cbcl_ext_problems, cbcl_dep_score, cbcl_anx_score, 
         cbcl_somatic_score, cbcl_adhd_score, cbcl_odd_score, cbcl_conduct_score)

# Select CBCL variables at 1-year follow-up
cbcl_t1 <- parent_cbcl %>%
  filter(eventname == "1_year_follow_up_y_arm_1") %>%
  select(subjectkey, cbcl_tot_problems, cbcl_int_problems, cbcl_ext_problems, cbcl_dep_score, cbcl_anx_score, 
         cbcl_somatic_score, cbcl_adhd_score, cbcl_odd_score, cbcl_conduct_score)

# Select CBCL variables at 2-year follow-up
cbcl_t2 <- parent_cbcl %>%
  filter(eventname == "2_year_follow_up_y_arm_1") %>%
  select(subjectkey, cbcl_tot_problems, cbcl_int_problems, cbcl_ext_problems, cbcl_dep_score, cbcl_anx_score, 
         cbcl_somatic_score, cbcl_adhd_score, cbcl_odd_score, cbcl_conduct_score)

# Rename KSADS variables
ksad_parent <- ksad_parent %>%
  rename(p_panic_present = ksads_5_857_p) %>%
  rename(p_sepanx_present = ksads_7_861_p) %>%
  rename(p_socanx_present = ksads_8_863_p) %>%
  rename(p_genanx_present = ksads_10_869_p) %>%
  rename(p_ocd_present = ksads_11_917_p) %>%
  rename(p_odd_present = ksads_15_901_p) %>%
  rename(p_conduct_present = ksads_16_897_p) 

# Sum KSADS anxiety disorders into 'number of anxiety disorders' 
ksad_parent <- ksad_parent %>%
  mutate(anxiety_disorders = 
           p_genanx_present + 
           p_panic_present + 
           p_sepanx_present + 
           p_socanx_present + 
           p_ocd_present) 

# Recode 'number of anxiety disorders' into 'any anxiety disorder present' yes (1) or no (0)
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 0] <- 0
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 1] <- 1
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 2] <- 1
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 3] <- 1
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 4] <- 1
ksad_parent$anxiety_disorders[ksad_parent$anxiety_disorders == 5] <- 1

# Set data type and level labels
ksad_parent$anxiety_disorders <- factor(ksad_parent$anxiety_disorders, levels = c( 1,0), labels = c("Yes", "No" ))

# Sum KSADS behavioural disorders into 'number of behavioural disorders' 
ksad_parent <- ksad_parent %>%
  mutate(behavioural_disorders = 
           p_odd_present + 
           p_conduct_present) 

# Recode 'number of behavioural disorders' into 'any behavioural disorder present' yes (1) or no (0)
ksad_parent$behavioural_disorders[ksad_parent$behavioural_disorders == 0] <- 0
ksad_parent$behavioural_disorders[ksad_parent$behavioural_disorders == 1] <- 1
ksad_parent$behavioural_disorders[ksad_parent$behavioural_disorders == 2] <- 1

# Set data type and level labels
ksad_parent$behavioural_disorders <- factor(ksad_parent$behavioural_disorders, levels = c( 1,0), labels = c("Yes", "No" ))

# Select KSAD variables at baseline
ksad_parent_baseline <- ksad_parent %>%  
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(subjectkey, anxiety_disorders, behavioural_disorders)   

# Select KSAD variables at 2-year follow-up 
ksad_parent_2_year <- ksad_parent %>%  
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  select(subjectkey, anxiety_disorders, behavioural_disorders)   

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
#
# Interview age at baseline 
tbi_baseline <- tbi_baseline %>% 
  mutate(interview_age = interview_age / 12)

# Interview age at 2-year follow-up 
tbi_2_year <- tbi_2_year %>% 
  mutate(interview_age = interview_age / 12)

# Recode sex 
#
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
  left_join(ksad_parent_baseline, by = "subjectkey") %>%
  left_join(cbcl_t0, by = "subjectkey") %>%
  left_join(ethnicity, by = "subjectkey") %>%
  left_join(parent_mh, by = "subjectkey") %>%
  left_join(parent_demographics, by = "subjectkey") %>%
  left_join(site, by = "subjectkey") %>%
  left_join(advlife, by = "subjectkey") %>%
  left_join(neighbourhood_safety_baseline, by = "subjectkey") %>%
  left_join(family_conflict_baseline, by = "subjectkey") %>%
  unique()


# Merge 2-year follow-up files
psyc_data_t2 <- tbi_1_year %>%
  left_join(tbi_2_year, by = "subjectkey") %>%
  left_join(oi_1_year, by = "subjectkey") %>%
  left_join(oi_2_year, by = "subjectkey") %>%
  left_join(ksad_parent_2_year, by = "subjectkey") %>%
  left_join(cbcl_t2, by = "subjectkey") %>%
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
# Create column
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

# Combine ortho injury cases from 1- and 2-year follow-up 
#
# Create column
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
label(psyc_data$broken_bones) <- "Orthopaedic injury"
label(psyc_data$anxiety_disorders) <- "Present anxiety disorders"
label(psyc_data$behavioural_disorders) <- "Present behavioural disorders"
label(psyc_data$parent_mh_score) <- "Parent mental health score"
label(psyc_data$cbcl_anx_score) <- "Anxiety symptoms scale"
label(psyc_data$cbcl_dep_score) <- "Depression symptoms scale"
label(psyc_data$cbcl_adhd_score) <- "ADHD symptoms"
label(psyc_data$cbcl_odd_score) <- "Oppositional defiant problems"
label(psyc_data$cbcl_conduct_score) <- "Conduct problems"
label(psyc_data$cbcl_int_problems) <- "Internalising symptoms scale"
label(psyc_data$cbcl_ext_problems) <- "Externalising symptoms scale"
label(psyc_data$cbcl_tot_problems) <- "Total problems scale"

# Create demographic and descriptive table for sample 
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
head(table_1_filename)

# Rename 2-year follow-up variables for better readability 
label(psyc_data_t2$household_income) <- "Household Income "
label(psyc_data_t2$sex) <- "Sex"
label(psyc_data_t2$race_ethnicity) <- "Ethnicity"
label(psyc_data_t2$interview_age) <- "Age at baseline"
label(psyc_data_t2$advlife_cat) <- "Traumatic events"
label(psyc_data_t2$safety) <- "Neighbourhood safety"
label(psyc_data_t2$family_conflict_scale) <- "Family conflict scale"
label(psyc_data_t2$parent_mh_score) <- "Parent mental health score"
label(psyc_data_t2$tbi_both_years) <- "TBI"
label(psyc_data_t2$anxiety_disorders) <- "Any present anxiety disorder"
label(psyc_data_t2$behavioural_disorders) <- "Any present behavioural disorder"
label(psyc_data_t2$cbcl_anx_score) <- "Anxiety symptoms scale"
label(psyc_data_t2$cbcl_dep_score) <- "Depression symptoms scale"
label(psyc_data_t2$cbcl_adhd_score) <- "ADHD symptoms"
label(psyc_data_t2$cbcl_odd_score) <- "Oppositional defiant problems"
label(psyc_data_t2$cbcl_conduct_score) <- "Conduct problems"
label(psyc_data_t2$cbcl_int_problems) <- "Internalising symptoms scale"
label(psyc_data_t2$cbcl_ext_problems) <- "Externalising symptoms scale"
label(psyc_data_t2$cbcl_tot_problems) <- "Total problems scale"



##############################################################################
#
# Multiple imputation for missing covariates using missForest package
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
                       anxiety_disorders, behavioural_disorders, cbcl_anx_score, cbcl_dep_score,
                       cbcl_adhd_score, cbcl_odd_score, cbcl_conduct_score, cbcl_int_problems,
                       cbcl_ext_problems, cbcl_tot_problems)

# Make sure subject ID is not imputed 
imputed_data$subjectkey <- NULL 

# Create imputed dataframe 
imputed_data <- as.data.frame(imputed_data)

# Visualise missing data 
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 200)
vis_miss(imputed_data, sort_miss = TRUE, warn_large_data = FALSE)
ggsave(paste(output_dir, "missing_psych_symptoms_1.png", sep = ""))

         
set.seed(123)

# Record start time 
start_time <- Sys.time() 

# Impute covariates but not outcomes  
imputed_data_baseline <- missForest(imputed_data[c(1:12)])

# Record end time 
end_time <- Sys.time() 

# Report time duration of missForest imputation 
end_time - start_time 

# Check normalised root mean squared error 
imputed_data_baseline$OOBerror 

# Create data frame with imputed covariates and outcomes 
imputed_data_baseline <- imputed_data_baseline$ximp 
imputed_data_baseline <- data.frame(imputed_data[13:22], imputed_data_baseline)

# Save new file for analyses 
write.csv(imputed_data_baseline, paste(output_dir, "imputed_baseline_KSAD_CBCL.csv", sep = ""))

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
                         safety, family_conflict_scale, tbi_both_years, broken_bones_both_years, injury_type, site_id_l, 
                         anxiety_disorders, behavioural_disorders, cbcl_anx_score, cbcl_dep_score,
                         cbcl_adhd_score, cbcl_odd_score, cbcl_conduct_score, cbcl_int_problems,
                         cbcl_ext_problems, cbcl_tot_problems)

# Make sure subject ID is not imputed 
imputed_data_2$subjectkey <- NULL 

# Create imputed dataframe 
imputed_data_2 <- as.data.frame(imputed_data_2)

# Visualise missing data 
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 200)
vis_miss(imputed_data_2, sort_miss = TRUE, warn_large_data = FALSE)
ggsave(paste(output_dir, "missing_psych_symptoms_2.png", sep = ""))

set.seed(123)

# Record start time 
start_time <- Sys.time() 

# Impute covariates but not outcomes  
imputed_data_t2 <- missForest(imputed_data_2[c(1:12)])

# Record end time 
end_time <- Sys.time() 

# Report time duration of missForest imputation 
end_time - start_time 

# Check normalized root mean squared error 
imputed_data_t2$OOBerror 

# Create data frame with imputed covariates and outcomes  
imputed_data_t2 <- imputed_data_t2$ximp 
imputed_data_t2 <- data.frame(imputed_data_2[13:22], imputed_data_t2)

# Save new file for analyses 
write.csv(imputed_data_t2, paste(output_dir, "imputed_t2_KSAD_CBCL.csv", sep = ""))

##############################################################################
#
# Load imputed datasets to avoid running imputation again 
#
#############################################################################

# Load imputed dataset for baseline 
imputed_data_baseline <- read.csv(paste(output_dir, "imputed_baseline_KSAD_CBCL.csv", sep = ""))

# Load imputed dataset for 2-year follow-up 
imputed_data_t2 <- read.csv(paste(output_dir, "imputed_t2_KSAD_CBCL.csv", sep = ""))

##############################################################################
#
# Complete statistical analyses
#
# Mixed-effects models with study site as a random effect 
# Adjusted model 1 covariates are age, sex, ethnicity and household income
# Adjusted model 2 covariates are adverse life event, neighbourhood safety, 
#  parent mental health and family conflict 
#
#############################################################################

#
# KSADS-5 analyses 
#

# Baseline
#

# Ensure baseline variables are stored as factors
imputed_data_baseline$sex <- as.factor(imputed_data_baseline$sex)
imputed_data_baseline$injury_type <- as.factor(imputed_data_baseline$injury_type)
imputed_data_baseline$anxiety_disorders <- as.factor(imputed_data_baseline$anxiety_disorders)
imputed_data_baseline$behavioural_disorders <- as.factor(imputed_data_baseline$behavioural_disorders)

# Any present anxiety disorder at baseline
#

# Change reference categories 
imputed_data_baseline$injury_type <- factor(imputed_data_baseline$injury_type, 
                                            levels = c("TBI", "Ortho", "None"))

imputed_data_baseline$anxiety_disorders <- relevel(imputed_data_baseline$anxiety_disorders, ref = "No")
imputed_data_baseline$injury_type <- relevel(imputed_data_baseline$injury_type, ref = "None")

# Run unadjusted model 
anxiety_disorders_model_unadjusted <- glmer(anxiety_disorders ~ injury_type +
                                              (1|site_id_l),
                                            data = imputed_data_baseline,
                                            family = binomial(link = "logit"),
                                            control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
anxiety_disorders_model_adjusted1 <- glmer(anxiety_disorders ~ injury_type +
                                             interview_age +
                                             sex +
                                             race_ethnicity +
                                             household_income +
                                             (1|site_id_l),
                                           data = imputed_data_baseline,
                                           family = binomial(link = "logit"),
                                           control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
anxiety_disorders_model_adjusted2 <- glmer(anxiety_disorders ~ injury_type +
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
                                           control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)


# Any present behavioural disorder at baseline
#

# Change reference categories 
imputed_data_baseline$behavioural_disorders <- relevel(imputed_data_baseline$behavioural_disorders, ref = "No")
imputed_data_baseline$injury_type <- relevel(imputed_data_baseline$injury_type, ref = "None")

# Run unadjusted model
behavioural_disorders_model_unadjusted <- glmer(behavioural_disorders ~ injury_type +
                                                  (1|site_id_l),
                                                data = imputed_data_baseline,
                                                family = binomial(link = "logit"),
                                                control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)



# Run adjusted model 1
behavioural_disorders_model_adjusted1 <- glmer(behavioural_disorders ~ injury_type +
                                                 interview_age +
                                                 sex +
                                                 race_ethnicity +
                                                 household_income +
                                                 (1|site_id_l),
                                               data = imputed_data_baseline,
                                               family = binomial(link = "logit"),
                                               control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
behavioural_disorders_model_adjusted2 <- glmer(behavioural_disorders ~ injury_type +
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
                                               control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

#
# KSADS-5 2-year follow-up
#

# Check variables are stored as factors 
imputed_data_t2$sex <- as.factor(imputed_data_t2$sex)
imputed_data_t2$injury_type <- as.factor(imputed_data_t2$injury_type)
imputed_data_t2$anxiety_disorders <- as.factor(imputed_data_t2$anxiety_disorders)
imputed_data_t2$behavioural_disorders <- as.factor(imputed_data_t2$behavioural_disorders)
imputed_data_t2$injury_type <- factor(imputed_data_t2$injury_type, 
                                            levels = c("TBI", "Ortho", "None"))

# Any present anxiety disorder at 2-year follow-up 
#

# Change reference categories 
imputed_data_t2$anxiety_disorders <- relevel(imputed_data_t2$anxiety_disorders, ref = "No")
imputed_data_t2$injury_type <- relevel(imputed_data_t2$injury_type, ref = "None")

# Run unadjusted model 
anxiety_disorders_model_unadjusted_t2 <- glmer(anxiety_disorders ~ injury_type +
                                              (1|site_id_l),
                                            data = imputed_data_t2,
                                            family = binomial(link = "logit"),
                                            control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
anxiety_disorders_model_adjusted1_t2 <- glmer(anxiety_disorders ~ injury_type +
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
anxiety_disorders_model_adjusted2_t2 <- glmer(anxiety_disorders ~ injury_type +
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
imputed_data_t2$behavioural_disorders <- relevel(imputed_data_t2$behavioural_disorders, ref = "No")

# Run unadjusted model 
behavioural_disorders_model_unadjusted_t2 <- glmer(behavioural_disorders ~ injury_type +
                                                  (1|site_id_l),
                                                data = imputed_data_t2, family = binomial(link = "logit"),
                                                control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
behavioural_disorders_model_adjusted1_t2 <- glmer(behavioural_disorders ~ injury_type +
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
behavioural_disorders_model_adjusted2_t2 <- glmer(behavioural_disorders ~ injury_type +
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

##############################################################################
#
# CBCL internalising, externalising and total symptoms score analyses 
#
#############################################################################

# Internalising symptoms score at baseline
#

# Run unadjusted model 
cbcl_int_symptoms_model_unadjusted <- lmer(cbcl_int_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_baseline
)

# Run adjusted model 1
cbcl_int_symptoms_model_adjusted1 <- lmer(cbcl_int_problems ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_baseline, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
cbcl_int_symptoms_model_adjusted2 <- lmer(cbcl_int_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Externalising symptoms score at baseline
#

# Run unadjusted model 
cbcl_ext_symptoms_model_unadjusted <- lmer(cbcl_ext_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_baseline
)


# Run adjusted model 1
cbcl_ext_symptoms_model_adjusted1 <- lmer(cbcl_ext_problems ~ injury_type + 
                                            interview_age + 
                                            sex + 
                                            race_ethnicity + 
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_baseline, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)



# Run adjusted model 2
cbcl_ext_symptoms_model_adjusted2 <- lmer(cbcl_ext_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Total problems score at baseline 
#

# Run unadjusted model 
cbcl_tot_symptoms_model_unadjusted <- lmer(cbcl_tot_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_baseline
                                           
)

# Run adjusted model 1
cbcl_tot_symptoms_model_adjusted1 <- lmer(cbcl_tot_problems ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_baseline, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
cbcl_tot_symptoms_model_adjusted2 <- lmer(cbcl_tot_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)


# Internalising symptoms score at 2-year follow-up
#

# Run unadjusted model 
cbcl_int_symptoms_model_unadjusted_t2 <- lmer(cbcl_int_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_t2
)

# Run adjusted model 1
cbcl_int_symptoms_model_adjusted1_t2 <- lmer(cbcl_int_problems ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_t2, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
cbcl_int_symptoms_model_adjusted2_t2 <- lmer(cbcl_int_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Externalising symptoms score at 2-year follow-up
#

# Run unadjusted model 
cbcl_ext_symptoms_model_unadjusted_t2 <- lmer(cbcl_ext_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_t2
)

# Run adjusted model 1
cbcl_ext_symptoms_model_adjusted1_t2 <- lmer(cbcl_ext_problems ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_t2, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
cbcl_ext_symptoms_model_adjusted2_t2 <- lmer(cbcl_ext_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)


# Total problems score at 2-year follow-up
#

# Run unadjusted model 
cbcl_tot_symptoms_model_unadjusted_t2 <- lmer(cbcl_tot_problems ~ injury_type +
                                             (1|site_id_l),
                                           data = imputed_data_t2
)

# Run adjusted model 1
cbcl_tot_symptoms_model_adjusted1_t2 <- lmer(cbcl_tot_problems ~ injury_type +
                                            interview_age +
                                            sex +
                                            race_ethnicity +
                                            household_income +
                                            (1|site_id_l),
                                          data = imputed_data_t2, 
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
cbcl_tot_symptoms_model_adjusted2_t2 <- lmer(cbcl_tot_problems ~ injury_type +
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
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)


##############################################################################
#
# CBCL symptom subscale analyses 
#
#############################################################################

# Anxiety symptoms at baseline
#

# Run unadjusted model 
anxiety_symptoms_model_unadjusted <- lmer(cbcl_anx_score ~ injury_type +
                                            (1|site_id_l),
                                          data = imputed_data_baseline,
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
anxiety_symptoms_model_adjusted1 <- lmer(cbcl_anx_score ~ injury_type +
                                           interview_age +
                                           sex +
                                           race_ethnicity +
                                           household_income +
                                           (1|site_id_l),
                                         data = imputed_data_baseline, 
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
anxiety_symptoms_model_adjusted2 <- lmer(cbcl_anx_score ~ injury_type +
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
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Depression symptoms at baseline
#

# Run unadjusted model 
depression_symptoms_model_unadjusted <- lmer(cbcl_dep_score ~ injury_type +
                                               (1|site_id_l),
                                             data = imputed_data_baseline,
                                             control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
depression_symptoms_model_adjusted1 <- lmer(cbcl_dep_score ~ injury_type +
                                              interview_age +
                                              sex +
                                              race_ethnicity +
                                              household_income +
                                              (1|site_id_l),
                                            data = imputed_data_baseline, 
                                            control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
depression_symptoms_model_adjusted2 <- lmer(cbcl_dep_score ~ injury_type +
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
                                            control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# ADHD symptoms at baseline
#

# Run unadjusted model 
adhd_symptoms_model_unadjusted <- lmer(cbcl_adhd_score ~ injury_type +
                                         (1|site_id_l),
                                       data = imputed_data_baseline,
                                       control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
adhd_symptoms_model_adjusted1 <- lmer(cbcl_adhd_score ~ injury_type +
                                        interview_age +
                                        sex +
                                        race_ethnicity +
                                        household_income +
                                        (1|site_id_l),
                                      data = imputed_data_baseline, 
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
adhd_symptoms_model_adjusted2 <- lmer(cbcl_adhd_score ~ injury_type +
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
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# ODD symptoms at baseline
#

# Run unadjusted model 
odd_symptoms_model_unadjusted <- lmer(cbcl_odd_score ~ injury_type +
                                        (1|site_id_l),
                                      data = imputed_data_baseline,
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
odd_symptoms_model_adjusted1 <- lmer(cbcl_odd_score ~ injury_type +
                                       interview_age +
                                       sex +
                                       race_ethnicity +
                                       household_income +
                                       (1|site_id_l),
                                     data = imputed_data_baseline, 
                                     control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
odd_symptoms_model_adjusted2 <- lmer(cbcl_odd_score ~ injury_type +
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
                                     control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Conduct symptoms at baseline
#

# Run unadjusted model 
conduct_symptoms_model_unadjusted <- lmer(cbcl_conduct_score ~ injury_type +
                                            (1|site_id_l),
                                          data = imputed_data_baseline
)

# Run adjusted model 1
conduct_symptoms_model_adjusted1 <- lmer(cbcl_conduct_score ~ injury_type +
                                           interview_age +
                                           sex +
                                           race_ethnicity +
                                           household_income +
                                           (1|site_id_l),
                                         data = imputed_data_baseline, 
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
conduct_symptoms_model_adjusted2 <- lmer(cbcl_conduct_score ~ injury_type +
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
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

#
# CBCL symptom subscales at 2-year follow-up
#

# Anxiety symptoms at 2-year follow-up
#
# Run unadjusted model 
anxiety_symptoms_model_unadjusted_t2 <- lmer(cbcl_anx_score ~ injury_type +
                                            (1|site_id_l),
                                          data = imputed_data_t2,
                                          control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
anxiety_symptoms_model_adjusted1_t2 <- lmer(cbcl_anx_score ~ injury_type +
                                           interview_age +
                                           sex +
                                           race_ethnicity +
                                           household_income +
                                           (1|site_id_l),
                                         data = imputed_data_t2, 
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
anxiety_symptoms_model_adjusted2_t2 <- lmer(cbcl_anx_score ~ injury_type +
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
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Depression symptoms at 2-year follow-up
#

# Run unadjusted model 
depression_symptoms_model_unadjusted_t2 <- lmer(cbcl_dep_score ~ injury_type +
                                               (1|site_id_l),
                                             data = imputed_data_t2,
                                             control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
depression_symptoms_model_adjusted1_t2 <- lmer(cbcl_dep_score ~ injury_type +
                                              interview_age +
                                              sex +
                                              race_ethnicity +
                                              household_income +
                                              (1|site_id_l),
                                            data = imputed_data_t2, 
                                            control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
depression_symptoms_model_adjusted2_t2 <- lmer(cbcl_dep_score ~ injury_type +
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
                                            control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# ADHD symptoms at 2-year follow-up
#

# Run unadjusted model 
adhd_symptoms_model_unadjusted_t2 <- lmer(cbcl_adhd_score ~ injury_type +
                                         (1|site_id_l),
                                       data = imputed_data_t2,
                                       control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
adhd_symptoms_model_adjusted1_t2 <- lmer(cbcl_adhd_score ~ injury_type +
                                        interview_age +
                                        sex +
                                        race_ethnicity +
                                        household_income +
                                        (1|site_id_l),
                                      data = imputed_data_t2, 
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
adhd_symptoms_model_adjusted2_t2 <- lmer(cbcl_adhd_score ~ injury_type +
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
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# ODD symptoms at 2-year follow-up
#

# Run unadjusted model 
odd_symptoms_model_unadjusted_t2 <- lmer(cbcl_odd_score ~ injury_type +
                                        (1|site_id_l),
                                      data = imputed_data_t2,
                                      control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

# Run adjusted model 1
odd_symptoms_model_adjusted1_t2 <- lmer(cbcl_odd_score ~ injury_type +
                                       interview_age +
                                       sex +
                                       race_ethnicity +
                                       household_income +
                                       (1|site_id_l),
                                     data = imputed_data_t2, 
                                     control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
odd_symptoms_model_adjusted2_t2 <- lmer(cbcl_odd_score ~ injury_type +
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
                                     control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Conduct symptoms at 2-year follow-up
#

# Run unadjusted model 
conduct_symptoms_model_unadjusted_t2 <- lmer(cbcl_conduct_score ~ injury_type +
                                            (1|site_id_l),
                                          data = imputed_data_t2
)

# Run adjusted model 1
conduct_symptoms_model_adjusted1_t2 <- lmer(cbcl_conduct_score ~ injury_type +
                                           interview_age +
                                           sex +
                                           race_ethnicity +
                                           household_income +
                                           (1|site_id_l),
                                         data = imputed_data_t2, 
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)

# Run adjusted model 2
conduct_symptoms_model_adjusted2_t2 <- lmer(cbcl_conduct_score ~ injury_type +
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
                                         control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
)


##############################################################################
#
# Create tables 
#
##############################################################################

#
# Tables for KSADS-5 analyses 
#

# Any anxiety disorder at baseline
#

# Create unadjusted model table
anxiety_unadjusted_table <- tbl_regression(anxiety_disorders_model_unadjusted, exponentiate = TRUE, label = injury_type ~ "Any present anxiety disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
anxiety_adjusted_table_model1 <- tbl_regression(anxiety_disorders_model_adjusted1, exponentiate = TRUE, 
                                                      pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                      include = injury_type, 
                                                      label = injury_type ~ "Any present anxiety disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
anxiety_adjusted_table_model2 <- tbl_regression(anxiety_disorders_model_adjusted2, exponentiate = TRUE, 
                                                      pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                      include = injury_type, 
                                                      label = injury_type ~ "Any present anxiety disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_anxiety <- tbl_merge(tbls = list(anxiety_unadjusted_table, anxiety_adjusted_table_model1, anxiety_adjusted_table_model2),
                                       tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Any behavioural disorder at baseline
#

# Create unadjusted model table
behavioural_unadjusted_table <- tbl_regression(behavioural_disorders_model_unadjusted, exponentiate = TRUE, label = injury_type ~ "Any present behavioural disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
behavioural_adjusted_table_model1 <- tbl_regression(behavioural_disorders_model_adjusted1, exponentiate = TRUE, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Any present behavioural disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
behavioural_adjusted_table_model2 <- tbl_regression(behavioural_disorders_model_adjusted2, exponentiate = TRUE, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Any present behavioural disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_behavioural <- tbl_merge(tbls = list(behavioural_unadjusted_table, behavioural_adjusted_table_model1, behavioural_adjusted_table_model2),
                                           tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))



# CBCL total problems score at baseline
#

# Create unadjusted model table
cbcl_tot_symptoms_unadjusted_table <- tbl_regression(cbcl_tot_symptoms_model_unadjusted, label = injury_type ~ "Total problems scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
cbcl_tot_symptoms_adjusted_table_model1 <- tbl_regression(cbcl_tot_symptoms_model_adjusted1, 
                                                          estimate_fun = function(x) style_number(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Total problems scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_tot_symptoms_adjusted_table_model2 <- tbl_regression(cbcl_tot_symptoms_model_adjusted2, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Total problems scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_tot_symptoms <- tbl_merge(tbls = list(cbcl_tot_symptoms_unadjusted_table, cbcl_tot_symptoms_adjusted_table_model1, cbcl_tot_symptoms_adjusted_table_model2),
                                           tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))



# Create table for any present anxiety and behavioural disorder and total problems at baseline
#

# Merge anxiety and behavioural disorder tables 
ksad_table_baseline <- tbl_stack(tbls = list(table_merge_anxiety, table_merge_behavioural, table_merge_cbcl_tot_symptoms)) %>%
  modify_caption("**Table 2. Association between mental health at age 9-10 and lifetime mTBI / orthopaedic injury**")

# Save table 
table_2_filename = paste(output_dir, "imputed_ksad_disorders_baseline.html", sep="")
gt::gtsave(as_gt(ksad_table_baseline), file = table_2_filename)
head(table_2_filename)



# Anxiety disorders at 2-year follow-up 
#

# Create unadjusted model table
anxiety_unadjusted_table_t2 <- tbl_regression(anxiety_disorders_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Any present anxiety disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
anxiety_adjusted_table_model1_t2 <- tbl_regression(anxiety_disorders_model_adjusted1_t2, exponentiate = TRUE, 
                                                      pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                      include = injury_type, 
                                                      label = injury_type ~ "Any present anxiety disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
anxiety_adjusted_table_model2_t2 <- tbl_regression(anxiety_disorders_model_adjusted2_t2, exponentiate = TRUE, 
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
behavioural_unadjusted_table_t2 <- tbl_regression(behavioural_disorders_model_unadjusted_t2, exponentiate = TRUE, label = injury_type ~ "Any present behavioural disorder") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
behavioural_adjusted_table_model1_t2 <- tbl_regression(behavioural_disorders_model_adjusted1_t2, exponentiate = TRUE, 
                                                          pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                          include = injury_type, 
                                                          label = injury_type ~ "Any present behavioural disorder") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
behavioural_adjusted_table_model2_t2 <- tbl_regression(behavioural_disorders_model_adjusted2_t2, exponentiate = TRUE, 
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


# CBCL total problems score at 2-year follow-up 
#

# Create unadjusted model table
cbcl_tot_symptoms_unadjusted_table_t2 <- tbl_regression(cbcl_tot_symptoms_model_unadjusted_t2, label = injury_type ~ "Total problems scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))


# Create adjusted model 1 table
cbcl_tot_symptoms_adjusted_table_model1_t2 <- tbl_regression(cbcl_tot_symptoms_model_adjusted1_t2, 
                                                             estimate_fun = function(x) style_number(x, digits = 3),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Total problems scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_tot_symptoms_adjusted_table_model2_t2 <- tbl_regression(cbcl_tot_symptoms_model_adjusted2_t2, 
                                                             pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                             include = injury_type, 
                                                             label = injury_type ~ "Total problems scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_tot_symptoms_t2 <- tbl_merge(tbls = list(cbcl_tot_symptoms_unadjusted_table_t2, cbcl_tot_symptoms_adjusted_table_model1_t2, cbcl_tot_symptoms_adjusted_table_model2_t2),
                                              tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))



# Create table for any present anxiety and behavioural disorder and total problems score at 2-year follow-up 
#

# Merge anxiety and behavioural disorder tables 
ksad_table_2year <- tbl_stack(tbls = list(table_merge_anxiety_t2, table_merge_behavioural_t2, table_merge_cbcl_tot_symptoms_t2)) %>%
  modify_caption("**Table 3. Association between mental health at age 11-12 and new mTBI / orthopaedic injury in the previous 24 months**")

# Save table 
table_2_filename = paste(output_dir, "imputed_ksad_disorders_t2.html", sep="")
gt::gtsave(as_gt(ksad_table_2year), file = table_2_filename)
head(table_2_filename)

#
#
# Create tables for CBCL symptom subscale analyses 
#

# Internalising symptoms at baseline
#

# Create unadjusted model table
cbcl_int_symptoms_unadjusted_table <- tbl_regression(cbcl_int_symptoms_model_unadjusted, label = injury_type ~ "Internalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
cbcl_int_symptoms_adjusted_table_model1 <- tbl_regression(cbcl_int_symptoms_model_adjusted1, 
                                                                estimate_fun = function(x) style_number(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Internalising symptoms scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_int_symptoms_adjusted_table_model2 <- tbl_regression(cbcl_int_symptoms_model_adjusted2, 
                                                                pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Internalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_int_symptoms <- tbl_merge(tbls = list(cbcl_int_symptoms_unadjusted_table, cbcl_int_symptoms_adjusted_table_model1, cbcl_int_symptoms_adjusted_table_model2),
                                                 tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Externalising symptoms at baseline
#

# Create unadjusted model table
cbcl_ext_symptoms_unadjusted_table <- tbl_regression(cbcl_ext_symptoms_model_unadjusted, label = injury_type ~ "Externalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
cbcl_ext_symptoms_adjusted_table_model1 <- tbl_regression(cbcl_ext_symptoms_model_adjusted1, 
                                                                estimate_fun = function(x) style_number(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Externalising symptoms scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_ext_symptoms_adjusted_table_model2 <- tbl_regression(cbcl_ext_symptoms_model_adjusted2, 
                                                                pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Externalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_ext_symptoms <- tbl_merge(tbls = list(cbcl_ext_symptoms_unadjusted_table, cbcl_ext_symptoms_adjusted_table_model1, cbcl_ext_symptoms_adjusted_table_model2),
                                                 tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))



# Anxiety symptoms at baseline 
#

# Create unadjusted model table
anxiety_symptoms_unadjusted_table <- tbl_regression(anxiety_symptoms_model_unadjusted, label = injury_type ~ "Anxiety symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
anxiety_symptoms_adjusted_table_model1 <- tbl_regression(anxiety_symptoms_model_adjusted1,  
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type, 
                                                         label = injury_type ~ "Anxiety symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
anxiety_symptoms_adjusted_table_model2 <- tbl_regression(anxiety_symptoms_model_adjusted2,  
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type, 
                                                         label = injury_type ~ "Anxiety symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_anxiety_symptoms <- tbl_merge(tbls = list(anxiety_symptoms_unadjusted_table, anxiety_symptoms_adjusted_table_model1, anxiety_symptoms_adjusted_table_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Depression symptoms at baseline 
#

# Create unadjusted model table
depression_symptoms_unadjusted_table <- tbl_regression(depression_symptoms_model_unadjusted, label = injury_type ~ "Depression symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
depression_symptoms_adjusted_table_model1 <- tbl_regression(depression_symptoms_model_adjusted1, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Depression symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
depression_symptoms_adjusted_table_model2 <- tbl_regression(depression_symptoms_model_adjusted2, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "Depression symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_depression_symptoms <- tbl_merge(tbls = list(depression_symptoms_unadjusted_table, depression_symptoms_adjusted_table_model1, depression_symptoms_adjusted_table_model2),
                                             tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# ADHD symptoms at baseline 
#
#

# Create unadjusted model table
adhd_symptoms_unadjusted_table <- tbl_regression(adhd_symptoms_model_unadjusted, label = injury_type ~ "ADHD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
adhd_symptoms_adjusted_table_model1 <- tbl_regression(adhd_symptoms_model_adjusted1,
                                                      pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                      include = injury_type, 
                                                      label = injury_type ~ "ADHD symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
adhd_symptoms_adjusted_table_model2 <- tbl_regression(adhd_symptoms_model_adjusted2, 
                                                      pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                      include = injury_type, 
                                                      label = injury_type ~ "ADHD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))


# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_adhd_symptoms <- tbl_merge(tbls = list(adhd_symptoms_unadjusted_table, adhd_symptoms_adjusted_table_model1, adhd_symptoms_adjusted_table_model2),
                                       tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# ODD symptoms at baseline 
#

# Create unadjusted model table
odd_symptoms_unadjusted_table <- tbl_regression(odd_symptoms_model_unadjusted, label = injury_type ~ "ODD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
odd_symptoms_adjusted_table_model1 <- tbl_regression(odd_symptoms_model_adjusted1, 
                                                     pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                     include = injury_type, 
                                                     label = injury_type ~ "ODD symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
odd_symptoms_adjusted_table_model2 <- tbl_regression(odd_symptoms_model_adjusted2,  
                                                     pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                     include = injury_type, 
                                                     label = injury_type ~ "ODD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))


# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_odd_symptoms <- tbl_merge(tbls = list(odd_symptoms_unadjusted_table, odd_symptoms_adjusted_table_model1, odd_symptoms_adjusted_table_model2),
                                      tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Conduct symptoms at baseline 
#

# Create unadjusted model table
conduct_symptoms_unadjusted_table <- tbl_regression(conduct_symptoms_model_unadjusted, label = injury_type ~ "Conduct problems score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
conduct_symptoms_adjusted_table_model1 <- tbl_regression(conduct_symptoms_model_adjusted1, 
                                                         estimate_fun = function(x) style_number(x, digits = 3),
                                                         include = injury_type, 
                                                         label = injury_type ~ "Conduct problems score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
conduct_symptoms_adjusted_table_model2 <- tbl_regression(conduct_symptoms_model_adjusted2, 
                                                         pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                         include = injury_type, 
                                                         label = injury_type ~ "Conduct problems score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_conduct_symptoms <- tbl_merge(tbls = list(conduct_symptoms_unadjusted_table, conduct_symptoms_adjusted_table_model1, conduct_symptoms_adjusted_table_model2),
                                          tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Create table for CBCL symptom subscales at baseline  
#

# Merge internalising, externalising, anxiety, depression, ADHD, ODD and conduct symptom tables
cbcl_symptoms_baseline_table <- tbl_stack(tbls = list(table_merge_cbcl_int_symptoms, table_merge_cbcl_ext_symptoms, table_merge_anxiety_symptoms, table_merge_depression_symptoms, table_merge_adhd_symptoms, table_merge_conduct_symptoms, table_merge_odd_symptoms)) %>%
  modify_caption("**Supplementary Table 1. Association between mental health symptoms at age 9-10 and lifetime mTBI / orthopaedic injury**")

# Save table 
table_1_filename = paste(output_dir, "imputed_cbcl_subscale_symptoms_baseline.html", sep="")
gt::gtsave(as_gt(cbcl_symptoms_baseline_table), file = table_1_filename)
head(table_1_filename)


# Internalising symptoms at 2-year follow-up 
#

# Create unadjusted model table
cbcl_int_symptoms_unadjusted_table_t2 <- tbl_regression(cbcl_int_symptoms_model_unadjusted_t2, label = injury_type ~ "Internalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
cbcl_int_symptoms_adjusted_table_model1_t2 <- tbl_regression(cbcl_int_symptoms_model_adjusted1_t2, 
                                                                estimate_fun = function(x) style_number(x, digits = 3),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Internalising symptoms scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_int_symptoms_adjusted_table_model2_t2 <- tbl_regression(cbcl_int_symptoms_model_adjusted2_t2, 
                                                                pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Internalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_int_symptoms_t2 <- tbl_merge(tbls = list(cbcl_int_symptoms_unadjusted_table_t2, cbcl_int_symptoms_adjusted_table_model1_t2, cbcl_int_symptoms_adjusted_table_model2_t2),
                                                 tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Externalising symptoms at 2-year follow-up 
#

# Create unadjusted model table
cbcl_ext_symptoms_unadjusted_table_t2 <- tbl_regression(cbcl_ext_symptoms_model_unadjusted_t2, label = injury_type ~ "Externalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
cbcl_ext_symptoms_adjusted_table_model1_t2 <- tbl_regression(cbcl_ext_symptoms_model_adjusted1_t2, 
                                                                estimate_fun = function(x) style_number(x, digits = 3),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Externalising symptoms scale") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
cbcl_ext_symptoms_adjusted_table_model2_t2 <- tbl_regression(cbcl_ext_symptoms_model_adjusted2_t2, 
                                                                pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                include = injury_type, 
                                                                label = injury_type ~ "Externalising symptoms scale") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_cbcl_ext_symptoms_t2 <- tbl_merge(tbls = list(cbcl_ext_symptoms_unadjusted_table_t2, cbcl_ext_symptoms_adjusted_table_model1_t2, cbcl_ext_symptoms_adjusted_table_model2_t2),
                                                 tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))


# Anxiety symptoms at 2-year follow-up 
#

# Create unadjusted model table
anxiety_symptoms_unadjusted_table_t2 <- tbl_regression(anxiety_symptoms_model_unadjusted_t2, label = injury_type ~ "Anxiety symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
anxiety_symptoms_adjusted_table_model1_t2 <- tbl_regression(anxiety_symptoms_model_adjusted1_t2,  
                                                               pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Anxiety symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
anxiety_symptoms_adjusted_table_model2_t2 <- tbl_regression(anxiety_symptoms_model_adjusted2_t2,  
                                                               pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Anxiety symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_anxiety_symptoms_t2 <- tbl_merge(tbls = list(anxiety_symptoms_unadjusted_table_t2, anxiety_symptoms_adjusted_table_model1_t2, anxiety_symptoms_adjusted_table_model2_t2),
                                                tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Depression symptoms at 2-year follow-up 
#

# Create unadjusted model table
depression_symptoms_unadjusted_table_t2 <- tbl_regression(depression_symptoms_model_unadjusted_t2, label = injury_type ~ "Depression symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
depression_symptoms_adjusted_table_model1_t2 <- tbl_regression(depression_symptoms_model_adjusted1_t2, 
                                                                  pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                  include = injury_type, 
                                                                  label = injury_type ~ "Depression symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
depression_symptoms_adjusted_table_model2_t2 <- tbl_regression(depression_symptoms_model_adjusted2_t2, 
                                                                  pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                                  include = injury_type, 
                                                                  label = injury_type ~ "Depression symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_depression_symptoms_t2 <- tbl_merge(tbls = list(depression_symptoms_unadjusted_table_t2, depression_symptoms_adjusted_table_model1_t2, depression_symptoms_adjusted_table_model2_t2),
                                                   tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# ADHD symptoms at 2-year follow-up 

# Create unadjusted model table
adhd_symptoms_unadjusted_table_t2 <- tbl_regression(adhd_symptoms_model_unadjusted_t2, label = injury_type ~ "ADHD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
adhd_symptoms_adjusted_table_model1_t2 <- tbl_regression(adhd_symptoms_model_adjusted1_t2,
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "ADHD symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
adhd_symptoms_adjusted_table_model2_t2 <- tbl_regression(adhd_symptoms_model_adjusted2_t2, 
                                                            pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                            include = injury_type, 
                                                            label = injury_type ~ "ADHD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_adhd_symptoms_t2 <- tbl_merge(tbls = list(adhd_symptoms_unadjusted_table_t2, adhd_symptoms_adjusted_table_model1_t2, adhd_symptoms_adjusted_table_model2_t2),
                                             tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# ODD symptoms at 2-year follow-up 
#

# Create unadjusted model table
odd_symptoms_unadjusted_table_t2 <- tbl_regression(odd_symptoms_model_unadjusted_t2, label = injury_type ~ "ODD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
odd_symptoms_adjusted_table_model1_t2 <- tbl_regression(odd_symptoms_model_adjusted1_t2, 
                                                           pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                           include = injury_type, 
                                                           label = injury_type ~ "ODD symptoms score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
odd_symptoms_adjusted_table_model2_t2 <- tbl_regression(odd_symptoms_model_adjusted2_t2,  
                                                           pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                           include = injury_type, 
                                                           label = injury_type ~ "ODD symptoms score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_odd_symptoms_t2 <- tbl_merge(tbls = list(odd_symptoms_unadjusted_table_t2, odd_symptoms_adjusted_table_model1_t2, odd_symptoms_adjusted_table_model2_t2),
                                            tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Conduct symptoms at 2-year follow-up 
#

# Create unadjusted model table
conduct_symptoms_unadjusted_table_t2 <- tbl_regression(conduct_symptoms_model_unadjusted_t2, label = injury_type ~ "Conduct problems score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 1 table
conduct_symptoms_adjusted_table_model1_t2 <- tbl_regression(conduct_symptoms_model_adjusted1_t2, 
                                                               estimate_fun = function(x) style_number(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Conduct problems score") %>% 
  bold_labels() %>%
  italicize_levels() %>%
  remove_row_type(variable = injury_type, type = 'reference') %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Create adjusted model 2 table
conduct_symptoms_adjusted_table_model2_t2 <- tbl_regression(conduct_symptoms_model_adjusted2_t2, 
                                                               pvalue_fun = function(x) style_pvalue(x, digits = 2),
                                                               include = injury_type, 
                                                               label = injury_type ~ "Conduct problems score") %>% 
  remove_row_type(variable = injury_type, type = 'reference') %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_column_merge(pattern = "{estimate} ({ci}) [{p.value}]", rows = !is.na(estimate))

# Merge unadjusted, adjusted model 1 and adjusted model 2 tables
table_merge_conduct_symptoms_t2 <- tbl_merge(tbls = list(conduct_symptoms_unadjusted_table_t2, conduct_symptoms_adjusted_table_model1_t2, conduct_symptoms_adjusted_table_model2_t2),
                                                tab_spanner = c("**Unadjusted**","**Adjusted, model 1**", "**Adjusted, model 2**"))

# Create table for CBCL symptom subscales at 2-year follow-up 
#

# Merge internalising, externalising, anxiety, depression, ADHD, ODD and conduct symptom tables
cbcl_symptoms_t2_table <- tbl_stack(tbls = list(table_merge_cbcl_int_symptoms_t2, table_merge_cbcl_ext_symptoms_t2, table_merge_anxiety_symptoms_t2, table_merge_depression_symptoms_t2, table_merge_adhd_symptoms_t2, table_merge_conduct_symptoms_t2, table_merge_odd_symptoms_t2)) %>%
  modify_caption("**Supplementary Table 2. Association between mental health symptoms at age 11-12 and new mTBI / orthopaedic injury in the previous 24 months**")

# Save table 
table_1_filename = paste(output_dir, "imputed_cbcl_subscale_symptoms_t2.html", sep="")
gt::gtsave(as_gt(cbcl_symptoms_t2_table), file = table_1_filename)
head(table_1_filename)
