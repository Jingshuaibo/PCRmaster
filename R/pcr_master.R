usethis::use_package(package = 'dplyr', type = 'Imports')
usethis::use_package(package = 'readxl', type = 'Imports')
usethis::use_package(package = 'writexl', type = 'Imports')

#1 Blank_index
#' @title Generate standard sample-blank matrix
#' @description Function Blank_index generates standard sample-blank matrix for subsequent analysis in package PCRmaster.
#'     It is usually the first step for PCR design.
#' @param date_list The original sample and blank list as data.frame, with each sample as one row.
#'     It should include:
#'     1) Label (column 1): The type of each sample ('Sample' or 'Blank).
#'     2) Sample (column 2): The name of each sample.
#'     3) Tag (column 3~...): Tags for connecting samples and blanks. Samples and their corresponding blanks should have
#'     the same tag in one column. More than one columns can be used.
#'
#' @returns A list contains:
#'     1) Blank_matrix: The standard sample-blank matrix as a 0/1 table
#'     2) Blank_index: A table contains all samples (column 1) and their corresponding blanks (column 2)
#' @export
#' @importFrom dplyr filter
#'
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

#1.2 library_division----
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

Lib_division <- function(blank_idx_output, step=8, max_lib=2, min_lib=1, num_lib=NULL, print_summary=FALSE){
  bnk_mt <-blank_idx_output
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

#1.3 PCR_designer----
Library_design <- function(sample_list, div_step=8, div_max_lib=2, div_min_lib=1,
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

#1.4 PCR table----
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
  n_board <- length(cut_tp)

  board_list <- as.vector(rep(NA, n_board), mode = 'list')
  names(board_list) <- paste0('BOARD-', 1:n_board)
  for (i in 1:n_board) {
    board_list[[i]] <- data.frame(matrix(nrow = 12, ncol = 8))
  }

  for (i in 1:n_board) {
    for (j in 1:4) {
      if((i-1)*32+j*8 < n_sample){
        data_tp <- tri_table[((i-1)*32+(j-1)*8+1) : ((i-1)*32+j*8)] }
      else{
        if((i-1)*32+(j-1)*8+1 <= n_sample){
          data_tp <- tri_table[((i-1)*32+(j-1)*8+1) : n_sample] }
        else{data_tp <- data.frame(matrix(nrow = 3, ncol = 8))}
      }
      board_list[[i]][((j-1)*3+1):(j*3), 1:ncol(data_tp)] <- data_tp
    }}

  return(board_list)
}

#1.5 PCR_list----
PCR_list <- function(blank_mt, cut_index=1, tail=NULL, include_PB=TRUE,
                     primer_list=NULL, primer_index=NULL, byID=NULL,
                     print_list=FALSE, print_board=FALSE){
  n_sample <- nrow(blank_mt)
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
    res_list[[i]] <- list('Sample_matrix'=lib_mt, 'PCR_list'=NA, 'Board'=NA)
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

  #board_table
  for (i in 1:nlib) {
    res_list[[i]][['Board']] <- PCR_table(res_list[[i]][['PCR_list']], by.id = byID)
  }

  #print_results
  if(print_list){
    for (i in 1:nlib) {
      write_xlsx(res_list[[i]][['PCR_list']], paste0('PCR_list_lib_',i,'.xlsx'))
    }}
  if(print_board){
    for (i in 1:nlib) {
      nboard <- length(res_list[[i]][['Board']])
      for (j in 1:nboard) {
        write_xlsx(res_list[[i]][['Board']][[j]], paste0('PCR_Board_',i,'-',j,'.xlsx'), col_names=F)
      }}}

  return(res_list)
}

