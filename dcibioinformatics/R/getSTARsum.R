# Helper function ----------------------------------------------------------------------

#' logCombine
#' It combines log columns of different samples
#' @import dplyr
#' @param df1 The star log file for sample 1
#' @param df2 The star log file sample 2
#' @return The joint df from df1 and df2

logCombine <- function(df1, df2) {
  dplyr::full_join(df1, df2, by="item")
}




#' mystarLogfile
#' Generates the full path to a STAR log file.
#' It assume that the log file is of the form
#' rootdir/stardir/suffix
#' @param rootdir The root directory under which the log file is stored
#' @param stardir The directory holding the files generated by STAR
#' @param suffix Suffix used by STAR to identify log file
#' @return Full path to the STAR log file
mystarLogfile <- function(rootdir, stardir, suffix = "Log.final.out") {
  file.path(rootdir, stardir, paste(stardir, suffix, sep=""))
}


# Export function-----------------------------------------------------------

utils::globalVariables(c("stardir","."))
#' getSTARsum
#' Produce a summary table for STAR alignment process
#' @import foreach
#' @import readr
#' @import dplyr
#' @param rootdir The root directory under which the count file is stored
#' @return A joint table for STAR alignment staitistics
#' @export

getSTARsum <- function(rootdir){
  coltypes <- list(col_character(), col_character())
  stardirs <- list.files(rootdir,full.names = FALSE)
  out <- foreach(stardir = stardirs, .combine = logCombine)%do%{
    logfile <- mystarLogfile(rootdir, stardir)
    suppressWarnings(readr::read_tsv(logfile, col_names = FALSE, col_types = coltypes))%>%
      dplyr::rename_at(vars(names(.)),~c("item",stardir))
  }
  return(out)
 }
  
  
  
  
  