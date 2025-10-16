# Import library
library(tidyverse)


#################################################################### I ### Load data ----
roswell_data <-
  readxl::read_xlsx(paste0(here::here(),
                           "/De-Identified Limited Body Comp Common Data Set for Moffitt Collab 01252023.xlsx"))


#################################################################### II ### harmonizing data with moffitt's----
roswell_data <- roswell_data %>% 
  rename(id = Study_ID, age_at_diagnosis = ageDX, 
         tnm_cs_mixed_group_stage = stage, 
         treatment_type = first_line_tx,
         bmi = BMIstartchemo, 
         sat_area_cm2 = PRE_SAT, imat_area_cm2 = PRE_IMAT, 
         vat_area_cm2 = PRE_VAT, total_fat_area = PRE_TAT,
         os_event = survstat, os_time = survtime, 
         pfs_time = rectime
         ) %>%
  mutate(tnm_cs_mixed_group_stage = case_when(
    tnm_cs_mixed_group_stage == 9                       ~ NA_character_,
    tnm_cs_mixed_group_stage == 1 |
      tnm_cs_mixed_group_stage == 2                     ~ "I-II",
    tnm_cs_mixed_group_stage == 3                       ~ "III",
    tnm_cs_mixed_group_stage == 4                       ~ "IV",
    TRUE                                                ~ as.character(tnm_cs_mixed_group_stage)
  )) %>% 
  mutate(treatment_type = case_when(
    treatment_type == "Neo"                             ~ "Upfront Neoadjuvant",
    treatment_type == "Adjuvant"                        ~ "Upfront Surgery",
    TRUE                                                ~ as.character(treatment_type)
  )) %>% 
  mutate(bmi_cat = case_when(
    bmi < 25                                            ~ "Underweight and normal weight",
    bmi >= 25 &
      bmi < 30                                          ~ "Overweight",
    bmi >= 30                                           ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight and normal weight", "Overweight", "Obese")))  %>% 
  mutate(pfs_event = case_when(
    recstat == 0                                        ~ 0,
    recstat == 1                                        ~ 1,
    recstat == 2                                        ~ 1
  ))

# write_rds(roswell_data, "roswell_data.rds")




