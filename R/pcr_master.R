# Blank_index----
#' @title Generate standard sample-blank matrix
#' @description Function \code{Blank_index} generates standard sample-blank matrix for subsequent analysis in package *PCRmaster*.
#'              It is usually the first step for PCR design.
#'
#' @param date_list The original sample and blank list as data.frame, with each sample as one row.
#'                  It should include:
#'                  1) **Label** (column 1): The type of each sample (\code{Sample} or \code{Blank});
#'                  2) **Sample** (column 2): The name of each sample;
#'                  3) **Tag** (column 3 ~ ...): Tags for connecting samples and blanks. Samples and their corresponding blanks
#'                  should have the same tag in one column. More than one columns can be used.
#'
#' @returns A list contains:
#'          1) **Blank_matrix**: The standard sample-blank matrix as a 0/1 table
#'          2) **Blank_index**: A table contains all samples (column 1) and their corresponding blanks (column 2)
#'
#' @export
#' @import dplyr
#'
#' @examples
#' data(sample_list_test)
#' cat('Sample list example:')
#' head(sample_list_test)
#'
#' bnk_idx_res <- Blank_index(sample_list_test)
#' cat('Standard sample-blank matrix:')
#' head(bnk_idx_res[['Blank_matrix']])
#' cat('Sample-blank table:')
#' head(bnk_idx_res[['Blank_index']])

Blank_index <- function(date_list){
  nDate <- ncol(date_list)-2
  colnames(date_list) <- c('Lable', 'Sample', paste0('Date_',1:nDate))

  sample_table <- filter(date_list, Lable == 'Sample')
  sample_name <- sample_table$Sample
  nSample <- nrow(sample_table)

  blank_table <- filter(date_list, Lable != 'Sample')
  blank_name <- blank_table$Sample
  nBlank <- nrow(blank_table)

  blank_mt <- data.frame(matrix(nrow = nSample, ncol = nBlank+1))
  colnames(blank_mt) <- c('Sample', blank_name)
  blank_mt$Sample <- sample_name
  for (i in 1:nSample) {
    for (j in 1:nBlank) {
      i_sample <- blank_mt$Sample[i]
      i_blank <- colnames(blank_mt)[j+1]
      i_sam_date <- filter(sample_table, Sample==i_sample)[3:(nDate+2)] %>% paste()
      i_bnk_date <- filter(blank_table, Sample==i_blank)[3:(nDate+2)] %>% paste()
      blank_mt[i, j+1] <- ifelse(any(i_sam_date == i_bnk_date), 1, 0)
    }
  }

  res_tp <- data.frame(Sample=sample_name, Blank=NA)
  for (i in 1:nrow(res_tp)) {
    i_sample <- res_tp$Sample[i]
    res_tp$Blank[i] <- colnames(blank_mt)[which(blank_mt[i,]==1)] %>% paste(collapse = ',')
  }

  res_list <- list(Blank_matrix=blank_mt, Blank_index=res_tp)
  return(res_list)
}

# Blank_statistic----
#' @title Calculation of the basic library information based on standard sample-blank matrix
#' @description Function \code{Blank_statistic} records the start and ending sample of a library in a standard sample-blank matrix,
#'              and calculates the number of samples and blanks and the sample-blank ratio (SBR, sample/blank).
#'
#' @param blank_mt A standard sample-blank matrix.
#' @param start The row index of the start sample in \code{blank_mt}.
#' @param end The row index of the ending sample in \code{blank_mt}.
#'
#' @returns A vector contains:
#'          1) Number of samples (\code{nSample});
#'          2) Number of blanks (\code{nBlank});
#'          3) Sum number of samples and blanks (\code{nTotal});
#'          4) Ratio of \code{nSample} to \code{nBlank}.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' set.seed(2026)
#' blank_mt_test <- data.frame('Sample'=paste0('S', 1:32),
#'                             'B1'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B2'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B3'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B4'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)))
#' cat('Standard sample-blank matrix example:')
#' print(blank_mt_test)
#'
#' blank_sta_res <- Blank_statistic(blank_mt_test, start=1, end=16)
#' print(blank_sta_res)

Blank_statistic <- function(blank_mt, start=NA, end=NA){
  n_sample <- nrow(blank_mt)
  start_tp <- ifelse(is.na(start), 1, start)
  end_tp <- ifelse(is.na(end), n_sample, end)

  data_tp <- blank_mt[start_tp:end_tp, ]
  res_tp <- rep(NA, 4)
  names(res_tp) <- c('nSample', 'nBlank', 'nTotal', 'SBR')
  res_tp[1] <- nrow(data_tp)
  res_tp[2] <- sum(colSums(data_tp[-1]) > 0)
  res_tp[3] <- res_tp[1] + res_tp[2]
  res_tp[4] <- res_tp[1]/res_tp[2]

  return(res_tp)
}

# Lib_division----
#' @title Divide samples into libraries.
#' @description Function \code{Lib_division} divides samples into multiple libraries and export basic information of each library
#'              under each division case.
#'
#' @param blank_matrix A standard sample-blank matrix.
#' @param step The step size for library division.
#' @param max_lib The maximum library number.
#' @param min_lib The minimum library number.
#' @param num_lib Apart from the maximum and minimum library number, user can also set a fixed library number.
#' @param print_summary When this parameter is 'TRUE', the function will print the final statistic information
#'                      for each library design scheme.
#'
#' @returns A list containing:
#'          1) **Blank_matrix**: the sample-blank matrix as input.
#'          2) **Library**: A list containing the detailed information of each library scheme named as \code{Library_x}
#'          (x: the library number). Under each library number, a library scheme can contain multiple specific
#'          secondary schemes named as \code{Case_y} (y: ID of secondary scheme). The information of each secondary
#'          scheme including:
#'            a) **Start_index**: The start sample index of each library.
#'            b) **Cut_index**: The division index of each library, that is, **Start_index - 1**.
#'            c) **lib_data**: The information of each library including its start (\code{Start}) and end (\code{End}) index,
#'            numbers of samples (\code{nSample}), blanks (\code{nBlank}), their sums (\code{nTotal}) and ratios
#'            (\code{SBR} = \code{nSample}/\code{nBlank}).
#'          3) **Summary**: The summary of all library schemes, including the total sample numbers (including blanks)
#'          of each library and the maximum, minimum and mean value of the total sample number and SBR across these libraries.
#'
#' @export
#' @import dplyr
#'
#' @examples
#' set.seed(2026)
#' blank_mt_test <- data.frame('Sample'=paste0('S', 1:32),
#'                             'B1'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B2'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B3'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B4'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)))
#' cat('Standard sample-blank matrix example:')
#' head(blank_mt_test)
#'
#' lib_div_res <- Lib_division(blank_mt_test, step=8, max_lib=3, min_lib=1, print_summary=F)
#' cat('Standard sample-blank matrix:')
#' head(lib_div_res[['Blank_matrix']])
#' cat('Library division summary:')
#' print(lib_div_res[['Summary']])
#'
Lib_division <- function(blank_matrix, step=8, max_lib=2, min_lib=1, num_lib=NULL, print_summary=FALSE){
  bnk_mt <-blank_matrix
  n_sample <- nrow(bnk_mt)
  if(n_sample %% step != 0){
    cut_idx <- seq(step, n_sample, step) }
  else{
    cut_idx <- seq(step, n_sample-1, step) }

  if(!is.null(num_lib)){
    min_lib <- num_lib
    max_lib <- num_lib }
  else{}

  lib_list <- as.vector(rep(NA, max_lib-min_lib+1), mode = 'list')
  names(lib_list) <- paste0('Library_', min_lib:max_lib)
  nlib <- length(lib_list)

  #i: number of library
  #lib_index: ID of list
  lib_index <- 1
  for (i in min_lib:max_lib) {
    cut_set <- combn(cut_idx, i-1)
    ncut <- choose(length(cut_idx), i-1)

    case_list <- as.vector(rep(NA, ncut), mode = 'list')
    names(case_list) <- paste0('Case_', 1:ncut)
    #j:cut_index_group
    for (j in 1:ncut) {
      case_list[[j]] <- list('Cut_index'=NA, 'Start_index'=NA, 'Library_info'=NA)
      lib_data <- data.frame(matrix(nrow = i, ncol = 7))
      colnames(lib_data) <- c('Library', 'Start', 'End', 'nSample', 'nBlank', 'nTotal', 'SBR')
      lib_data$Library <- paste0('Lib_', 1:i)

      if(i == 1){
        case_list[[j]][['Cut_index']] <- 0
        case_list[[j]][['Start_index']] <- 1
        lib_data$Start[i] <- 1
        lib_data$End[i] <- n_sample
        lib_data[1, 4:7] <- Blank_statistic(bnk_mt, 1, n_sample)
      }
      else{
        case_list[[j]][['Cut_index']] <- cut_set[,j]
        case_list[[j]][['Start_index']] <- c(1, cut_set[,j] + 1)
        lib_data$Start <- c(1, cut_set[,j] + 1)
        lib_data$End <- c(cut_set[,j], n_sample)
        for (k in 1:nrow(lib_data)) {
          lib_data[k, 4:7] <- Blank_statistic(bnk_mt, lib_data$Start[k], lib_data$End[k])
        }}

      case_list[[j]][['Library_info']] <- lib_data
    }

    lib_list[[lib_index]] <- case_list
    lib_index <- lib_index + 1
  }

  #Summary
  smr_tp <- data.frame(matrix(nrow = sum(choose(length(cut_idx), (min_lib-1):(max_lib-1))), ncol = max_lib+8))
  colnames(smr_tp) <- c('Library','Case', paste0('Lib_', 1:max_lib),
                        'Max_sam', 'Min_sam', 'Mean_sam', 'Max_SBR', 'Min_SBR', 'Mean_SBR')
  smr_tp[3:(2+max_lib)] <- 0

  i_index <- 1
  for (i in 1:nlib) {
    nlib_tp <- length(lib_list[[i]])
    smr_tp$Library[i_index:(i_index + nlib_tp - 1)] <- rep(names(lib_list)[i], nlib_tp)
    smr_tp$Case[i_index:(i_index + nlib_tp - 1)] <- names(lib_list[[i]])
    i_index <- i_index + nlib_tp
  }

  for (i in 1:nrow(smr_tp)) {
    lib_tp <- smr_tp$Library[i]
    case_tp <- smr_tp$Case[i]
    data_tp <- lib_list[[lib_tp]][[case_tp]][["Library_info"]]
    smr_tp[i, 3:(3+nrow(data_tp)-1)] <- data_tp$nTotal
    smr_tp$Max_sam[i] <- max(data_tp$nTotal)
    smr_tp$Min_sam[i] <- min(data_tp$nTotal)
    smr_tp$Mean_sam[i] <- mean(data_tp$nTotal)
    smr_tp$Max_SBR[i] <- max(data_tp$SBR)
    smr_tp$Min_SBR[i] <- min(data_tp$SBR)
    smr_tp$Mean_SBR[i] <- mean(data_tp$SBR)
  }
  if(print_summary){
    print(smr_tp) }
  else {}
  res_tp <- list(Blank_matrix=bnk_mt, Library=lib_list, Summary=smr_tp)

  return(res_tp)
}

# Lib_design----
#' @title Generate library design schemes
#' @description Function \code{Lib_design} generate different library design schemes based on original sample list and
#'              calculate basic information for each library including the number of samples and the ratio
#'              of samples to blanks.
#'
#' @param sample_list The original sample and blank list as data.frame, with each sample as one row.
#'                    It should include:
#'                    1) **Label** (column 1): The type of each sample ('Sample' or 'Blank);
#'                    2) **Sample** (column 2): The name of each sample;
#'                    3) **Tag** (column 3 to ...): Tags for connecting samples and blanks. Samples and their
#'                    corresponding blanks should have the same tag in one column. More than one columns can be used.
#' @param div_step The step size for library division.
#' @param div_max_lib The maximum library number.
#' @param div_min_lib The minimum library number.
#' @param div_num_lib Apart from the maximum and minimum library number, user can also set a fixed library number.
#' @param div_print_summary When \code{div_print_summary} is \code{TRUE}, the function will print the statistic information
#'                          for each library design scheme.
#'
#' @returns A list contains:
#'          1) **Blank_matrix**: Standard sample-blank matrix.
#'          2) **Blank_index**: A table contains all samples (column 1) and their corresponding blanks (column 2).
#'          3) **Library**: A list containing the detailed information of each library scheme named as 'Library_x'
#'          (x: the library number). Under each library number, a library scheme can contain multiple specific
#'          secondary schemes named as 'Case_y' (y: ID of secondary scheme). The information of each secondary
#'          scheme including:
#'            a) **Start_index**: The start sample index of each library.
#'            b) **Cut_index**: The division index of each library, that is, **Start_index - 1**.
#'            c) **lib_data**: The information of each library including its start (\code{Start}) and end (\code{End}) index,
#'            numbers of samples (\code{nSample}), blanks (\code{nBlank}), their sums (\code{nTotal}) and ratios
#'            (\code{SBR} = \code{nSample}/\code{nBlank}).
#'          4) **Library_summary**: The summary of all library schemes, including the total sample numbers (including blanks)
#'          of each library and the maximum, minimum and mean value of the total sample number and SBR across these libraries.
#'
#' @export
#'
#' @examples
#' data_test <- data(sample_list_test)
#' cat('Sample list example:')
#' head(data_test)
#'
#' lib_dsg_res <- Lib_design(sample_list_test, div_step=8, div_max_lib=3, div_min_lib=1)
#' cat('Standard sample-blank matrix:')
#' head(lib_dsg_res[['Blank_matrix']])
#' cat('Sample-blank table:')
#' head(lib_dsg_res[['Blank_index']])
#' cat('Library division summary:')
#' print(lib_dsg_res[['Library_summary']])

Lib_design <- function(sample_list, div_step=8, div_max_lib=2, div_min_lib=1,
                           div_num_lib=NULL, div_print_summary=FALSE){
  #sample-blank table
  bnk_idx_res <- Blank_index(sample_list)
  bnk_matrix <- bnk_idx_res[["Blank_matrix"]]

  #library division
  lib_div_res <- Lib_division(bnk_matrix, step=div_step, print_summary=div_print_summary,
                              max_lib=div_max_lib, min_lib=div_min_lib, num_lib=div_num_lib)

  #result
  res_list <- list("Blank_matrix"=bnk_matrix,
                   "Blank_index"=bnk_idx_res[["Blank_index"]],
                   "Library"=lib_div_res[["Library"]],
                   "Library_summary"=lib_div_res[["Summary"]])

  return(res_list)
}

# PCR table----
#' @title Generate 96-well PCR plates based on PCR list.
#'
#' @param pcr_list A data frame containing PCR information. It should **at least** contains:
#'                 1) **Sample**: Sample ID used to distinguish different samples.
#'                 2) **PCR_ID**: PCR ID used to distinguish different PCR replications for the same sample
#'                 (such as "X-1", "X-2", "X-3" for the sample "X"). The name of this column must correspond
#'                 with the parameter \code{by.id}.
#' @param by.id Character, the name of column which is used to be shown in the final PCR plates. This columnn
#'              should distinguish different PCR replications.
#'
#' @return A list containing the information of all 96-well PCR plates named as "PLATE-x" (x: number of each plate).
#'         Each plate is a 12-row × 8-column data frame corresponding to the layout of a 96-well plate, filled with
#'         PCR ID by **by.id** column. Every 3 rows corresponds to 8 samples (i.e., each sample occupies 3 rows × 1 column).
#'         Specially, if \code{by.id} is \code{NULL}, returns \code{0}.
#'
#' @export
#' @import dplyr
#' @examples
#' data_test <- data.frame('Sample'=paste0('S', rep(1:32, each=3)),
#'                         'PCR_ID'=paste0('S', rep(1:32, each=3), '-', rep(1:3, each=32)))
#' cat("PCR list example:")
#' head(data_test)
#'
#' pcr_table_res <- PCR_table(data_test, by.id='PCR_ID')
#' cat("Plate template:")
#' print(pcr_table_res[[1]])
#'
PCR_table <- function(pcr_list, by.id=NULL){
  if(is.null(by.id)){
    print('Please enter ID column!')
    return(0)
  }
  else{}

  n_pcr <- nrow(pcr_list)
  sample <- unique(pcr_list$Sample)
  n_sample <- length(sample)

  tri_table <- data.frame(matrix(nrow = 3, ncol = n_sample))
  for (i in 1:n_sample) {
    tri_table[i] <- filter(pcr_list, Sample==sample[i])[[by.id]]
  }

  cut_tp <- seq(1, n_pcr, 96)
  n_plate <- length(cut_tp)

  plate_list <- as.vector(rep(NA, n_plate), mode = 'list')
  names(plate_list) <- paste0('PLATE-', 1:n_plate)
  for (i in 1:n_plate) {
    plate_list[[i]] <- data.frame(matrix(nrow = 12, ncol = 8))
  }

  for (i in 1:n_plate) {
    for (j in 1:4) {
      if((i-1)*32+j*8 < n_sample){
        data_tp <- tri_table[((i-1)*32+(j-1)*8+1) : ((i-1)*32+j*8)] }
      else{
        if((i-1)*32+(j-1)*8+1 <= n_sample){
          data_tp <- tri_table[((i-1)*32+(j-1)*8+1) : n_sample] }
        else{data_tp <- data.frame(matrix(nrow = 3, ncol = 8))}
      }
      plate_list[[i]][((j-1)*3+1):(j*3), 1:ncol(data_tp)] <- data_tp
    }}

  return(plate_list)
}

# PCR_list----
#' @title Generate PCR list and 96-well PCR plates using sample-blank matrix
#'
#' @param blank_mt Standard sample-blank matrix.
#' @param cut_index A vector containing the start index of each library in the sample-blank matrix.
#' @param tail If \code{tail} is not \code{NULL}, this function will add it after each \code{PCR_ID} in the PCR list.
#' @param include_PB If \code{include_PB} is \code{TRUE}, this function will add 3 PCR blanks (PB-1/2/3) in the PCR list.
#' @param primer_list The primer list (data.frame) for assigning primers to each PCR. It is suggested to contain the sequence
#'                    and ID of each primer.
#' @param primer_index The start primer index for assignment in the primer list. If \code{primer_index} is \code{NULL},
#'                     the primer assignment will start with the first primer.
#' @param byID Name of column (character) which is used to be shown in the final PCR plates. This columnn
#'             should distinguish different PCR replications.
#' @param print_list If \code{print_list} is \code{TRUE}, this function will export PCR list as an excel file (.xlsx) under the
#'                   working directory.
#' @param print_plate If \code{print_plate} is \code{TRUE}, this function will export PCR plates as an excel file (.xlsx) under the
#'                    working directory.
#' @param print_pcrnum If \code{print_pcrnum} is \code{TRUE}, this function will export PCR number information as an excel file
#'                     (.xlsx) under the working directory.
#'
#' @returns This function returns a list containing:
#'          1) **Library-X**: The sample-blank matrix, PCR list and PCR plates for library X (X: ID number for each library);
#'          2) **PCR_number**: The number information for each library, including the presence/absence of each sample
#'          and the sum sample number for each library.
#'
#' @import dplyr
#' @import writexl
#' @export
#'
#' @examples
#' set.seed(2026)
#' blank_mt_test <- data.frame('Sample'=paste0('S', 1:32),
#'                             'B1'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B2'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B3'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)),
#'                             'B4'=sample(0:1, 32, replace = T, prob=c(0.75,0.25)))
#' cat("Standard sample-blank matrix:")
#' print(blank_mt_test)
#'
#' primer_test <- data.frame('PrimerID'=paste0('Primer-',1:100), 'Sequence'=NA)
#' for (i in 1:100){
#'   primer_test$Sequence[i] <- paste0(sample(c('a','c','t','g'), 6, replace = T), collapse = '')
#' }
#' cat("Primer list:")
#' head(primer_test)
#'
#' pcr_list_res <- PCR_list(blank_mt_test, cut_index=c(1,17), include_PB = TRUE,
#'                          primer_list = primer_test, primer_index = c(1,20), byID='PCR_ID',
#'                          print_list = FALSE,  print_plate = FALSE, print_pcrnum = FALSE)
#' cat("PCR list:")
#' head(pcr_list_res[[1]][['PCR_list']])
#' cat("PCR plate:")
#' head(pcr_list_res[[1]][['Plate']])
#' cat("PCR number information:")
#' head(pcr_list_res[['PCR_number']])
#'
PCR_list <- function(blank_mt, cut_index=1, tail=NULL, include_PB=TRUE,
                     primer_list=NULL, primer_index=NULL, byID='PCR_ID',
                     print_list=FALSE, print_plate=FALSE, print_pcrnum=FALSE){
  n_sample <- nrow(blank_mt)
  n_blank <- ncol(blank_mt)-1
  nlib <- length(cut_index)

  res_list <- as.list(rep(NA, nlib))
  names(res_list) <- paste0('Library_', 1:nlib)
  for (i in 1:nlib) {
    if(i < nlib){
      lib_mt <- blank_mt[cut_index[i]:(cut_index[i+1]-1),] }
    else{
      lib_mt <- blank_mt[cut_index[i]:n_sample,] }
    if(include_PB){
      lib_mt[paste0('PB',i)] <- 1 }
    else{}
    res_list[[i]] <- list('Sample_matrix'=lib_mt, 'PCR_list'=NA, 'Plate'=NA)
  }

  for (i in 1:nlib) {
    mt_tp <- res_list[[i]][['Sample_matrix']]
    sample_id <- c(mt_tp$Sample, colnames(mt_tp)[(which(colSums(mt_tp[-1]) > 0))+1])
    n_total <- length(sample_id)
    n_pcr <- 3*n_total

    #PCR_list
    list_tp <- data.frame(matrix(ncol = 3, nrow = n_pcr))
    colnames(list_tp) <- c('Sample', 'PCR_ID', 'ShortID')
    for (j in 1:ceiling(n_total/8)) {
      sam_sta <- (j-1)*8 + 1
      pcr_sta <- (j-1)*24 + 1
      sam_end <- ifelse(j*8 < n_total, j*8, n_total)
      pcr_end <- ifelse(j*8 < n_total, j*24, n_pcr)

      list_tp$Sample[pcr_sta:pcr_end] <- rep(sample_id[sam_sta:sam_end], 3)
      list_tp$PCR_ID[pcr_sta:pcr_end] <- paste0(list_tp$Sample[pcr_sta:pcr_end], '-', rep(1:3, each=length(sam_sta:sam_end)))
      list_tp$ShortID[pcr_sta:pcr_end] <- paste0(sam_sta:sam_end, '-', rep(1:3, each=length(sam_sta:sam_end)))
    }

    if(!is.null(tail)){
      list_tp$Sample <- paste0(list_tp$Sample, tail)}
    else{}

    res_list[[i]][['PCR_list']] <- list_tp
  }

  #add primer information
  if(!is.null(primer_list)){
    for (i in 1:nlib) {
      k <- ifelse(i <= length(primer_index), i, length(primer_index))
      pmr_sta <- ifelse(is.null(primer_index), 1, primer_index[k])
      pmr_end <- pmr_sta + nrow(res_list[[i]][['PCR_list']]) - 1
      res_list[[i]][['PCR_list']] <- cbind(res_list[[i]][['PCR_list']], primer_list[pmr_sta:pmr_end,])
    }
  }
  else{}

  #plate_table
  for (i in 1:nlib) {
    res_list[[i]][['Plate']] <- PCR_table(res_list[[i]][['PCR_list']], by.id = byID)
  }

  #PCR_count
  pcr_num <- data.frame(matrix(nrow = n_sample+n_blank+2, ncol = 3 + nlib))
  colnames(pcr_num) <- c('Lable', 'Sample', names(res_list), 'Total')
  pcr_num$Lable <- c(rep('Sample', n_sample), rep('Blank', n_blank), 'nSample', 'nPCR')
  pcr_num$Sample <- c(blank_mt$Sample, colnames(blank_mt)[-1], 'nSample', 'nPCR')
  for (i in 1:nlib) {
    for (j in 1:(nrow(pcr_num)-2)) {
      sample_tp <- pcr_num$Sample[j]
      pcr_num[j, i+2] <- ifelse(sample_tp %in% res_list[[i]][['PCR_list']][["Sample"]], 1, 0)
    }
    pcr_num[nrow(pcr_num)-1, i+2] <- sum(pcr_num[1:(nrow(pcr_num)-2), i+2])
    pcr_num[nrow(pcr_num), i+2] <- 3*pcr_num[nrow(pcr_num)-1, i+2]
  }
  pcr_num$Total <- rowSums(pcr_num[3:(2+nlib)])
  res_list[['PCR_number']] <- pcr_num

  #print_results
  if(print_list){
    for (i in 1:nlib) {
      write_xlsx(res_list[[i]][['PCR_list']], paste0('PCR_list_lib_',i,'.xlsx'))
    }}

  if(print_plate){
    for (i in 1:nlib) {
      nplate <- length(res_list[[i]][['Plate']])
      for (j in 1:nplate) {
        write_xlsx(res_list[[i]][['Plate']][[j]], paste0('PCR_Plate_',i,'-',j,'.xlsx'), col_names=F)
      }}}

  if(print_pcrnum){
    write_xlsx(res_list[['PCR_number']], 'PCR_number.xlsx')
  }

  return(res_list)
}


