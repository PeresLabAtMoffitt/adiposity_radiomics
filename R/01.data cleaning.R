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
  # Create variable
  mutate(bmi = weight / (height_m_ * height_m_)) %>% 
  mutate(bmi_cat = case_when(
    bmi < 25                    ~ "Underweight and normal weight",
    bmi >= 25 &
      bmi < 30                  ~ "Overweight",
    bmi >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight and normal weight", "Overweight", "Obese"))) %>% 
  mutate(tnm = factor(tnm, levels = c("Underweight and normal weight", "Overweight", "Obese"))) %>% 
  mutate(weight_date = as.Date(weight_date, format = "%m/%d/%Y"))

adipose_data %>% nrow()

# Create and merge indexed data
adipose_data <- adipose_data %>% 
  select(mrn, muscle_area_cm2, imat_area_cm2, vat_area_cm2, sat_area_cm2, total_fat_area) %>% 
  mutate(across(where(is.numeric), ~ (. / (adipose_data$height_m_ ^ 2)))) %>% 
  `colnames<-`(paste(colnames(.), "ind", sep = "_")) %>% 
  full_join(adipose_data, ., by = c("mrn" = "mrn_ind"))















write_rds(adipose_data, "adipose_data.rds")
