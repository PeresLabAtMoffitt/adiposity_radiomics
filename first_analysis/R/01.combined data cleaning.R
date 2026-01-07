# Import library
library(tidyverse)
library(lubridate)


#################################################################### I ### Load data ----

###_______________________ custom function _______________________###
fct_name_repair <- function(colnms) {
  tolower(gsub("[, ()=/]", "_", colnms))
}
###_______________________________________________________________###

path <- fs::path("","Volumes","Peres_Research", "Ovarian - Radiomics")
body_comp <-
  readxl::read_xlsx(paste0(path, "/CTbased adiposity",
                           "/dataset/CT data 09272021.xlsx"),
                    .name_repair = fct_name_repair) %>% 
  `colnames<-`(str_replace(colnames(.), "__", "_"))
# clinical_data <-
#   read_csv(paste0(path, "/CTbased adiposity",
#                   "/dataset/ForDan_clinical_modif_06092022.csv"))
updated_survival <- 
  read_rds(paste0(path, 
                  "/Grant applications/Miles for Moffitt/data/",
                  "processed data/Ovarian_Normalized_Radiomics_Features_30Apr2025.rds")) %>% 
  `colnames<-`(str_replace(colnames(.), "__", "_"))

roswell_data <-
  readxl::read_xlsx(paste0(path, "/Roswell CTs (May 2024)/Dataset",
                           "/De-ID Data Set for Moffitt 04032024.xlsx"))
# roswell_data_old <-
#   readxl::read_xlsx(paste0(here::here(), 
#                            "/De-Identified Limited Body Comp Common Data Set for Moffitt Collab 01252023.xlsx"))

#################################################################### II ### Cleaning data ----
body_comp <- body_comp %>% 
  purrr::keep(~!all(is.na(.))) %>%
  mutate(mrn = as.character(mrn)) %>% 
  mutate(baseline_ct_scan_date = as.Date(baseline_ct_scan_date))

updated_survival1 <- updated_survival %>% 
  mutate(mrn = as.character(mrn)) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  # Outcome
  mutate(os_time_from_dx = round(interval(start = date_of_diagnosis, end = fwdate_most_recent)/
                                   duration(n=1, units = "months"), 2)) %>%
  mutate(pfs_enddate = case_when(
    rec_event == 0 & 
      os_event == 1               ~ fwdate_most_recent, #only use date of death when death occurred but no recurrence
    TRUE                          ~ rec_event_date, # use recurrence date for everything else
  )) %>% 
  mutate(pfs_time_from_dx = round(interval(start = date_of_diagnosis, end = pfs_enddate)/
                                   duration(n=1, units = "months"), 2)) %>%
  # Date for merging
  mutate(baseline_ct_scan_date = as.Date(baseline_ct_scan_date, format = "%m/%d/%Y")) %>% 
  # clinical
  mutate(debulking_status = case_when(
    debulking_status == "incomplete records"             ~ NA_character_,
    str_detect(debulking_status, "^optimal \\(")         ~ "Optimal",
    str_detect(debulking_status, "^suboptimal \\(")      ~ "Suboptimal",
    str_detect(debulking_status, "^complete \\(")        ~ "Complete"
  )) %>% 
  mutate(race = str_to_sentence(race_cancer_registry)) %>% 
  # Remove unnecessary variables
  select(-c(starts_with("nor_"), matches("^f[1-9]")))
  # select(mrn, vital_status_date : rec_time_from_dx) %>% 
  # `colnames<-`(c("mrn", paste0(colnames(.)[2:ncol(.)], "_Jan2025")
  #                ))

# clinical_data1 <- full_join(updated_survival1, 
#                  clinical_data #%>% 
#                    # select(mrn, date_of_diagnosis, age_at_diagnosis,
#                    #        baseline_ct_scan_date, 
#                    #        date_of_first_recurrence, date_of_last_followup,
#                    #        fwdate_most_recent,
#                    #        tnm_cs_mixed_group_stage, debulking_status,
#                    #        age_at_first_recurrence : ncol(.)
#                           ) %>% 
                   # Need to recode recurrence date and time to use date_of_last_followup and not fwdate_most_recent
                 #   mutate(rec_event_date = coalesce(date_of_first_recurrence, date_of_last_followup)) %>% 
                 #   mutate(recurrence_time = round(interval(start = first_treatment_date, end = rec_event_date)/
                 #                                    duration(n=1, units = "months"), 2)) %>% 
                 #   select(-date_of_first_recurrence, -date_of_last_followup),
                 # by = "mrn") #%>% 
  # Add old survival data to new for the missing patients in the new data
  # mutate(vital_new = coalesce(vital_status_Nov2024, vital_new)) %>% 
  # mutate(vital_date_new = coalesce(vital_status_date_Nov2024, vital_date_new)) %>% 
  # 
  # mutate(os_event_Jan2025 = case_when(
  #   vital_status_Jan2025 == "ALIVE"      ~ 0,
  #   vital_status_Jan2025 == "DEAD"      ~ 1
  # ))
  # mutate(os_time = coalesce(os_time_Nov2024, os_time)) %>% 
  # mutate(rec_event = coalesce(rec_event_Nov2024, rec_event)) %>% 
  # mutate(recurrence_time = coalesce(recurrence_time_Nov2024, recurrence_time)) %>% 
  # mutate(rec_event_date = coalesce(rec_event_date_Nov2024, rec_event_date))


moffitt_data <- left_join(body_comp, updated_survival1, 
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
  mutate(SMI = muscle_area_cm2 / (height_m_ * height_m_)) %>% 
  mutate(sarcopenia = case_when(
    SMI <= 38.73                                     ~ "Yes",
    TRUE                                             ~ "No"
  )) %>% 
  mutate(martin = case_when(
    SMI < 41                                         ~ "Sarcopenia",
    TRUE                                             ~ "No sarcopenia"
  )) %>% 
  # mutate(weight_date = as.Date(weight_date, format = "%m/%d/%Y")) %>% 
  
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
  ))

# adipose_data %>% nrow()
write_rds(moffitt_data, 
          paste0("data/", "updated moffitt body comp data_", today(), ".rds"))
write_csv(moffitt_data, 
          paste0("data/", "updated moffitt body comp data_", today(), ".csv"))

write_rds(moffitt_data, 
          paste0(path, "/CTbased adiposity",
                 "/dataset/processed data/",
                 "updated moffitt body comp data_", today(), ".rds"))
write_csv(moffitt_data, 
          paste0(path, "/CTbased adiposity",
                 "/dataset/processed data/",
                 "updated moffitt body comp data_", today(), ".csv"))



# Create and merge indexed data
# adipose_data <- adipose_data %>% 
#   select(mrn, muscle_area_cm2, imat_area_cm2, vat_area_cm2, sat_area_cm2, total_fat_area) %>% 
#   mutate(across(where(is.numeric), ~ (. / (adipose_data$height_m_ ^ 2)))) %>% 
#   `colnames<-`(paste(colnames(.), "ind", sep = "_")) %>% 
#   full_join(adipose_data, ., by = c("mrn" = "mrn_ind"))

# Clean Roswell data - get same names etc as Moffit
roswell_data <- roswell_data %>% 
  mutate(recstat = case_when(
    recstat == 9             ~ NA_real_,
    TRUE                     ~ recstat
  )) %>% 
  rename(# id = studyID, 
         age_at_diagnosis = ageDX, 
         tnm_cs_mixed_group_stage = combined_stage, 
         histology = HistologyDescription,
         # treatment_type = first_line_tx,
         debulking_status = DEBULK,
         weight_kg = chemo_start_wgt,
         height_m_ = ht_first_met,
         bmi = BMIstartchemo_new,
         sat_area_cm2 = PRE_SUBQ_1, imat_area_cm2 = PRE_IMA_1, 
         vat_area_cm2 = PRE_VAT_1, total_fat_area = PRE_TAT_1,
         muscle_area_cm2 = PRE_SMA_1,
         os_event = survstat, os_time_from_dx = survtime,
         pfs_event = recstat, pfs_time_from_dx = rectime
  ) %>% 
  mutate(ascites2 = case_when(
    Ascites_1_Yes_No == 0                            ~ "absence",
    Ascites_1_Yes_No == 1                            ~ "presence"
  )) %>% 
  mutate(debulking_status = case_when(
    debulking_status == "SUB"           ~ "Suboptimal",
    debulking_status == "UNKNOWN"       ~ NA_character_,
    TRUE                                ~ debulking_status
  ))

comb_dat <- bind_rows(moffitt_data, 
                      roswell_data,
                      .id = "site") %>% 
  mutate(site = case_when(
    site == 1              ~ "Moffitt",
    site == 2              ~ "Roswell"
  ))

comb_dat <- comb_dat %>% 
  mutate(age_at_diagnosis = round(age_at_diagnosis, 0)) %>% 
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
  mutate(BCI = 
           (1.4079 + 
              (
                (0.0039 * sat_area_cm2 - 0.0036 * vat_area_cm2 - 0.0309 * muscle_area_cm2) / 
                  (height_m_ * height_m_)
                
              )) * 
           (muscle_area_cm2 / imat_area_cm2)
  ) %>% 
  mutate(BCI_roswell_cat = case_when(
    BCI < 3.48                             ~ "Low",
    BCI >= 3.48                            ~ "High"
  ), BCI_roswell_cat = factor(BCI_roswell_cat, levels = c("High", "Low"))) %>% 
  # mutate(BCI_moffitt_cat = case_when(
  #   BCI < 3.44                             ~ "Low",
  #   BCI >= 3.44                            ~ "High"
  # ), BCI_moffitt_cat = factor(BCI_moffitt_cat, levels = c("High", "Low"))) %>% 
  # mutate(BCI_cat = case_when(
  #   BCI < 3.54                             ~ "Low",
  #   BCI >= 3.54                            ~ "High"
  # ), BCI_cat = factor(BCI_cat, levels = c("High", "Low"))) %>% 
  # mutate(raceeth1 = case_when(
  #   raceeth == "White Non-Hispanic"        ~ "NHWhite",
  #   TRUE                                   ~ "Others"
  # )) %>% 
  # mutate(raceeth1 = factor(raceeth1, levels = c("NHWhite", "Others"))) %>% 
  mutate(debulking_status = str_to_sentence(debulking_status),
         debulking_status = case_when(
           debulking_status == "Complete"            ~ "Optimal",
           is.na(debulking_status)                   ~ "Unknown",
           TRUE                                      ~ debulking_status
         ), debulking_status = factor(debulking_status, 
                                      levels = c("Suboptimal", "Optimal", "Unknown"))
  ) %>% 
  mutate(tnm_cs_mixed_group_stage = case_when(
    tnm_cs_mixed_group_stage == 1 |
      tnm_cs_mixed_group_stage == 2           ~ "I-II",
    tnm_cs_mixed_group_stage == 3 |#             ~ "III",
      tnm_cs_mixed_group_stage == 4           ~ "III-IV",
    tnm_cs_mixed_group_stage == 9             ~ "Unknown",
    is.na(tnm_cs_mixed_group_stage)           ~ "Unknown",
    TRUE                                      ~ as.character(tnm_cs_mixed_group_stage)
  )) %>% 
  mutate(tnm_cs_mixed_group_stage = factor(tnm_cs_mixed_group_stage, levels = c("I-II", 
                                                                                "III-IV",# "III", 
                                                                                "Unknown")))

write_rds(comb_dat, 
          paste0("data/", "combined Moffitt-Roswell body comp data_", today(), ".rds"))
write_csv(comb_dat, 
          paste0("data/", "combined Moffitt-Roswell body comp data_", today(), ".csv"))

write_rds(comb_dat, 
          paste0(path, "/CTbased adiposity",
                 "/dataset/processed data/",
                 "combined Moffitt-Roswell body comp data_", today(), ".rds"))
write_csv(comb_dat, 
          paste0(path, "/CTbased adiposity",
                 "/dataset/processed data/",
                 "combined Moffitt-Roswell body comp data_", today(), ".csv"))





  