# Import library
library(tidyverse)


#################################################################### I ### Load data ----
path <- fs::path("","Volumes","Peres_Research", "Ovarian - Radiomics", "CTbased adiposity")

radiomics <-
  read_csv(paste0(path,"/ForDan.csv"))
clinical_data <-
  read_csv(paste0(path,"/ForDan_clinical.csv"))


#################################################################### II ### Data cleaning ----
adipose_data <- left_join(radiomics, clinical_data, by="mrn") %>% # keep patient with ct image
  # Eight patients were excluded due to coverage artifacts
  filter(is.na(exclude)) %>% 
  # Eight patients were excluded due to incomplete or missing CT images 
  filter(!is.na(muscle_area_cm2)) 

adipose_data %>% nrow()

write_rds(adipose_data, "adipose_data.rds")
