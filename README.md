# AIM
The project aim to evaluate pre-treatment adiposity features measured by computed tomography and clinical features in the survival of women with ovarian cancer.

# Cleaning
Then you can run the "01.data cleaning.R" script which does the clinical recoding/cleaning/create new variables.

# Statistical analyses
The Rmd files run the statistical analyses.
Investigate data.Rmd will run the whole analysis (upfront chemo and upfront surgery included) so you shouldn't need to run the stratified Rmd.

Note : In case we need to change the variables included in the imputation, the code is dynamic so you don't need to select the best "Number of multiple imputations" aka "m" argument in the mice function. So just changing the variables names that should be included is needed.
