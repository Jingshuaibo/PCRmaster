## code to prepare `DATASET` dataset goes here
# read data
sample_list_test <- readxl::read_excel('data-raw/sample_list.xlsx')
primer_list_test <- readxl::read_excel('data-raw/primer_list.xlsx')

# save data as .rda
usethis::use_data(sample_list_test, overwrite = TRUE)
usethis::use_data(primer_list_test, overwrite = TRUE)
