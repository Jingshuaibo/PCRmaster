#' Sample list data
#'
#' An example sample data containing 54 samples and 17 blanks.
#'
#' @format A tibble data frame with 71 rows and 4 columns：
#' \describe{
#'   \item{Lable}{character, type of each sample ('Sample' or 'Blank')}
#'   \item{ID}{character, unique sample identifiers}
#'   \item{Date_1}{character, the 1st processing date in YYYYMMDD format}
#'   \item{Date_2}{character, the 2nd processing date in YYYYMMDD format}
#' }
#' @source Laboratory sample processing records
"sample_list_test"

#' Primer list data
#'
#' An example primer dataset containing primer ID, tag sequence, and tag ID information for 200 primers.
#'
#' @format A tibble data frame with 200 rows and 3 columns:
#' \describe{
#'   \item{PrimerID}{character, unique primer identifiers in "fwh1-XXX" format}
#'   \item{tag}{character, 8-base primer tag sequences}
#'   \item{tagID}{numeric, tag identifiers, ranging from 1 to 200}
#' }
#' @source Laboratory primer management system
"primer_list_test"
