## Prepare `DATASET` dataset

impact_missing <- read.csv(file = "impact_missing.csv", header = TRUE)
impact_imputed <- read.csv(file = "impact_imputed.csv", header = TRUE)
immune_original <- read.csv(file = "immune_original.csv", header = TRUE,
                            encoding = "UTF-8")

usethis::use_data(impact_missing, impact_imputed, immune_original,
                  overwrite = TRUE)
