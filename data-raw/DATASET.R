# How to prepare the codon_usage dataset
codon_usage <- read.table("data-raw/codon_usage.csv", sep = ",",
                          stringsAsFactors = FALSE, header = TRUE)
usethis::use_data(codon_usage, overwrite = TRUE)
