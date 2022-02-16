# Import library
library(tidyverse)
library(lubridate)


#################################################################### I ### Load data ----

###_______________________ custom function _______________________###
fct_name_repair <- function(colnms) {
  tolower(gsub("[, ()=/]", "_", colnms))
}
###_______________________________________________________________###

path <- fs::path("","Volumes","Peres_Research", "Ovarian - Radiomics", "CTbased adiposity")
radiomics <-
  readxl::read_xlsx(paste0(path,"/dataset/CT data 09272021.xlsx"),
                    .name_repair = fct_name_repair) %>% 
  `colnames<-`(str_replace(colnames(.), "__", "_"))
clinical_data <-
  read_csv(paste0(path,"/dataset/ForDan_clinical.csv")) %>% 
  select(-"...1")


#################################################################### II ### Data cleaning ----
adipose_data <- left_join(radiomics, clinical_data, 
                          by=c("mrn", "date_of_diagnosis", "baseline_ct_scan_date")) %>% # keep patient with ct image
  purrr::keep(~!all(is.na(.))) %>%
  # Eight patients were excluded due to coverage artifacts
  filter(is.na(should_exclude)) %>% 
  # Eight patients were excluded due to incomplete or missing CT images 
  filter(!is.na(muscle_area_cm2)) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  mutate(weight_kg = case_when(
    weight_unit == "kg" |
      weight_unit == "Kg"           ~ weight,
    weight_unit == "lbs"            ~ weight / 2.205
  )) %>% 
  # Create variable
  mutate(bmi = weight_kg / (height_m_ * height_m_)) %>% 
  mutate(bmi_cat = case_when(
    bmi < 25                    ~ "Underweight and normal weight",
    bmi >= 25 &
      bmi < 30                  ~ "Overweight",
    bmi >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight and normal weight", "Overweight", "Obese")))  %>% 
  mutate(bmi_cat2 = case_when(
    bmi < 25                                         ~ "<25",
    bmi >= 25 &
      bmi < 30                                       ~ "25-29",
    bmi >= 30 &
      bmi < 35                                       ~ "30-34",
    bmi >= 35                                        ~ "≥35"
  )) %>% 
  mutate(bmi_cat2 = factor(bmi_cat2, levels = c("<25", "25-29", "30-34", "≥35")))  %>% 
  mutate(SMI = muscle_area_cm2 / (height_m_ * height_m_)) %>% 
  mutate(sarcopenia = case_when(
    SMI <= 38.73                                     ~ "Yes",
    TRUE                                             ~ "No"
  )) %>% 
  mutate(martin = case_when(
    SMI < 41                                         ~ "Sarcopenia",
    TRUE                                             ~ "No sarcopenia"
  ))%>% 
  mutate(tnm_cs_mixed_group_stage = case_when(
    tnm_cs_mixed_group_stage == 1 |
      tnm_cs_mixed_group_stage == 2           ~ "I-II",
    tnm_cs_mixed_group_stage == 3             ~ "III",
    tnm_cs_mixed_group_stage == 4             ~ "IV"
  )) %>% 
  mutate(tnm_cs_mixed_group_stage = factor(tnm_cs_mixed_group_stage, levels = c("IV", "III", "I-II"))) %>% 
  mutate(weight_date = as.Date(weight_date, format = "%m/%d/%Y")) %>% 
  mutate(raceeth1 = case_when(
    raceeth == "White Non-Hispanic"        ~ "NHWhite",
    TRUE                                   ~ "Others"
  )) %>% 
  mutate(raceeth1 = factor(raceeth1, levels = c("NHWhite", "Others"))) %>% 
  mutate(ascites = case_when(
    is.na(ascites)                         ~ "absence",
    str_detect(ascites, "tumor")           ~ "absence",
    str_detect(ascites, "moderate")        ~ "moderate",
    str_detect(ascites, "severe")          ~ "severe",
    str_detect(ascites, "mild")            ~ "mild"
  )) %>% 
  mutate(ascites2 = case_when(
    str_detect(ascites, "absence")                   ~ "absence",
    str_detect(ascites, "moderate|severe|mild")      ~ "presence"
  )) %>% 
  mutate(debulking_status = as.factor(debulking_status))

adipose_data %>% nrow()

# Create and merge indexed data
adipose_data <- adipose_data %>% 
  select(mrn, muscle_area_cm2, imat_area_cm2, vat_area_cm2, sat_area_cm2, total_fat_area) %>% 
  mutate(across(where(is.numeric), ~ (. / (adipose_data$height_m_ ^ 2)))) %>% 
  `colnames<-`(paste(colnames(.), "ind", sep = "_")) %>% 
  full_join(adipose_data, ., by = c("mrn" = "mrn_ind"))

adipose_data <- adipose_data %>% 
  select(mrn, muscle_area_cm2, imat_area_cm2, vat_area_cm2, sat_area_cm2, total_fat_area,
         ends_with("_ind"),
         muscle_mean_hu, imat_mean_hu, muscle_hu_sd, v_s_ratio) %>% 
  mutate(across(where(is.numeric), ~ sd(.))) %>% 
  `colnames<-`(paste(colnames(.), "SD", sep = "_")) %>% 
  full_join(adipose_data, ., by = c("mrn" = "mrn_SD")) %>% 
  mutate(muscle_area_cm2_HR = muscle_area_cm2 / muscle_area_cm2_SD,
         imat_area_cm2_HR = imat_area_cm2 / imat_area_cm2_SD,
         vat_area_cm2_HR = vat_area_cm2 / vat_area_cm2_SD,
         sat_area_cm2_HR = sat_area_cm2 / sat_area_cm2_SD,
         total_fat_area_HR = total_fat_area / total_fat_area_SD,
         muscle_area_cm2_ind_HR = muscle_area_cm2_ind / muscle_area_cm2_ind_SD,
         imat_area_cm2_ind_HR = imat_area_cm2_ind / imat_area_cm2_ind_SD,
         vat_area_cm2_ind_HR = vat_area_cm2_ind / vat_area_cm2_ind_SD,
         sat_area_cm2_ind_HR = sat_area_cm2_ind / sat_area_cm2_ind_SD,
         total_fat_area_ind_HR = total_fat_area_ind / total_fat_area_ind_SD,
         muscle_mean_hu_HR = muscle_mean_hu / muscle_mean_hu_SD,
         imat_mean_hu_HR = imat_mean_hu / imat_mean_hu_SD,
         v_s_ratio_HR = v_s_ratio / v_s_ratio_SD)


write_rds(adipose_data, "adipose_data.rds")

# End cleaning
