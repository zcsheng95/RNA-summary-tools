# Helper function ----------------------------------------------------------------------

#' logCombine
#' It combines log columns of different samples
#' @import dplyr
#' @param df1 The star log file for sample 1
#' @param df2 The star log file sample 2
#' @return The joint df from df1 and df2

Combine <- function(df1, df2) {
  dplyr::full_join(df1, df2, by=names(df1)[1])
}


# Export function-----------------------------------------------------------

utils::globalVariables(c("stardir","."))
#' getSTARsum
#' Produce a summary table for STAR alignment process
#' @param files The list of files to read in
#' @param libnames Sample names, make sure they are paired with files 
#' @param verbose True or False. Show read in progress 
#' @return A joint table for STAR alignment staitistics
#' @export

getSTARsum <- function(files, libnames, verbose = TRUE){
  ### Predefine column types
  coltypes <- list(col_character(), col_character())
  mylogfiles <- tibble(myfname = files, mylibname = libnames)  
  out <- foreach(myfile = iter(mylogfiles, by ="row"), .combine = Combine)%do%{
    myfname <- myfile[["myfname"]]
    mylibname <-  myfile[["mylibname"]]
    if(verbose == TRUE){
    cat("Read in",mylibname,"\n")
    }
    suppressWarnings(readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>%
      dplyr::rename_at(vars(names(.)),~c("item", mylibname)))
    
  }
  return(out)
 }
  
  
  
  
  