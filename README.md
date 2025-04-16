Analysis code and results for the study

#  Only Anxiety Remains Reliably Associated with Paediatric Mild Traumatic Brain Injury at Two Years Follow-up After Adjusting for Pre-Existing Mental Health


<p align="center">
	<a href="https://en.wikipedia.org/wiki/R_(programming_language)"><img
		alt="R Programming Language"
		src="https://img.shields.io/badge/Language-R-%232268BB.svg"></a>
	<a href="https://opensource.org/licenses/MIT"><img
		alt="MIT License"
		src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

Pre-printed as: 

> Revill, G., Poole, N., Carlisi, C., David, A. S., Bell, V. (2024). Association Between Paediatric Mild Traumatic Brain Injury and Two-Year Psychiatric Outcomes Largely Explained by Pre-Existing Mental Health Problems. medRxiv, [pre-print doi: 10.1101/2024.03.22.24304723](https://doi.org/10.1101/2024.03.22.24304723)

This archive contains the [R code](https://en.wikipedia.org/wiki/R_(programming_language)) for the analysis reported in the above study. 

This repository contains:

1.  [Revill_et_al_psychiatric_symptoms_diagnoses.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_psychiatric_symptoms_diagnoses.R) - R script to complete the psychiatric diagnosis and scales outcomes analysis
2.  [Revill_et_al_psychiatric_service_use.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_psychiatric_service_use.R) - R script to complete the psychiatric services use analysis
3.  [Revill_et_al_2_1_psm_matchit_code.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_2_1_psm_matchit_code.R) - R script to complete propensity score matching analysis
4.  [Revill_et_al_senspoweranalysis_mTBI_v_ortho.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_senspoweranalysis_mTBI_v_ortho.R) - R script to calculate sensitivity power analysis for mTBI vs orthopaedic controls
5.  [Revill_et_al_senspoweranalysis_mTBI_v_uninjured.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_senspoweranalysis_mTBI_v_uninjured.R) - R script to calculate sensitivity power analysis for mTBI vs uninjured controls
6.  [Revill_et_al_supple_psm_matchit_code.R](https://github.com/GraceRevill/pTBI-neuropsychiatric-outcomes/blob/main/Revill_et_al_supple_psm_matchit_code.R) - Supplementary analysis comparing 1:1 vs 2:1 propensity score matching


### Dataset

This study used data from version 4.0 of the [Adolescent Brain Cognitive Development Study](https://en.wikipedia.org/wiki/ABCD_Study) (ABCD).

The ABCD Study is a "prospective longitudinal study starting at the ages of 9-10 and following participants for 10 years. The study includes a diverse sample of nearly 12,000 youth enrolled at 21 research sites across the USA. It measures brain development (via structural, task functional, and resting state functional imaging), social, emotional, and cognitive development, mental health, substance use and attitudes, gender identity and sexual health, biospecimens, as well as a variety of physical health, and environmental factors."

The data can be downloaded from the NIMH [ABCD Data Repository](https://nda.nih.gov/abcd).

### Platform and package versions

R language version, and package versions used to generate the results are:

R Version 4.2.1<br>
Platform x86_64 Windows<br>
<br>
Package version for readr is 2.1.4<br>
Package version for dplyr is 1.1.2<br>
Package version for tidyr is 1.3.0<br>
Package version for survey is 4.0<br>
Package version for lme4 is 1.1.32<br>
Package version for Hmisc is 5.0.1<br>
Package version for table1 is 1.4.3<br>
Package version for epiR is 2.0.61<br>
Package version for jtools is 2.2.1<br>
Package version for gtsummary is 1.7.1<br>
Package version for gt is 0.9.0<br>
Package version for afex is 1.3.0<br>
Package version for missForest is 1.5<br>
Package version for doParallel is 1.0.17<br>
Package version for visdat is 0.6.0<br>
Package version for doRNG is 1.8.6<br>
Package version for ggplot2 is 3.3.5<br>
Package version for MatchIt is 4.5.5<br>
Package version for sjPlot is 2.8.14<br>
