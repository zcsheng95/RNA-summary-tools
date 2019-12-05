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
#' @param rootdir The root directory under which the count file is stored
#' @param pipeline The rnaseq pipeline version
#' @param level indicate the structure of star output, level = 0 means basename contains libname, level = 1 and 2 indicates different nested dir structure.
#' @return A joint table for STAR alignment staitistics
#' @export

getSTARsum <- function(rootdir, level){
  ### Predefine column types
  coltypes <- list(col_character(), col_character())
  mylogfiles <- getSTARcntfnames(level, rootdir, "Log.final.out", recursive = TRUE)
  out <- foreach(myfile = iter(mylogfiles, by ="row"), .combine = Combine)%do%{
    myfname <- myfile[["myfname"]]
    mylibname <-  myfile[["mylibname"]]
    cat("Read in",mylibname,"\n")
    suppressWarnings(readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>%
      dplyr::rename_at(vars(names(.)),~c("item", mylibname)))
    
  }
  return(out)
 }
  
  
  
  
  