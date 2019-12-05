# Helper function ----------------------------------------------------------------------

#' This is a convenience function for extracting the STAR count files
#' It will return a tibble with two colums:
#' 1. Absolute file name for count file
#' 2. library name
#' Perhaps, it is best to extend it to contain  absolute file name
#' for the final log file
#' @param mylevel indicate the structure of star output, level = 0 means basename contains libname, level = 1 and 2 indicates different nested dir structure.
#' @param rootdir rootdir rootdir for STAR output
#' @param mypattern pattern for star out put file
#' @param recursive boolean value for locating files recursively or not
#' @param emptylab If no lib name is found, use this to replace
#' @param duplab functions to make unique names for each sample

getSTARcntfnames <- function(mylevel, rootdir, mypattern, recursive, emptylab = "emptylib", duplab = function(x){make.unique(x, sep = "-")}) {
  myfiles <- normalizePath(list.files(path = rootdir,
                                      pattern = mypattern,
                                      include.dirs = TRUE,
                                      full.names = TRUE,
                                      recursive = recursive))
  
  ## Use the file basename to get library name
  if (mylevel == 0) {
    ## Get basename after removing pattern
    mylibname <- str_remove(basename(myfiles), mypattern)
    ## If filename equals pattern (no prefix) set it to emptylab
    mylibname <- dplyr::case_when(mylibname == "" ~ emptylab, TRUE ~ mylibname)
  }
  ## Use the directory structure to get library names
  ## For this format foo/a/star/ReadsPerGene.out.tab , foo/b/star/ReadsPerGene.out.tab
  ## Use mylevel =2 to get a and b
  ## For this format foo/star/a/ReadsPerGene.out.tab , foo/star/b/ReadsPerGene.out.tab
  ## Use mylevel =1 to get a and b
  else {
    mylibname <- sapply(str_split(dirname(myfiles),"/"), dplyr::nth, n = -mylevel)
  }
  tibble::tibble(myfname = myfiles, mylibname = duplab(mylibname))
}


# Export function -------------------------------------------------------------------
utils::globalVariables(c("stardir","."))

#' getSTARcounts
#' Extract strand-specific STAR counts for each sample and then merge all counts into a single table
#' @import dplyr
#' @import foreach
#' @import readr
#' @import tibble
#' @importFrom utils read.table
#' @importFrom iterators iter
#' @param rootdir The root directory under which the count file is stored
#' @param strand The strand-specific protocol used in RNA sequencing. Choose one from \code{c("first", "second", "unstranded")}
#' @param level indicate the structure of star output, level = 0 means basename contains libname, level = 1 and 2 indicates different nested dir structure.
#' @return A strand-specific STAR count summary table
#' @export 

getSTARcounts <- function(rootdir,strand, level) {
  if(!strand %in% c("first","second","unstranded")) stop("Missing strand argument! You must specify strand-specific protocol from c('first','second','unstranded')")
  col = switch(strand, second = 4,first= 3 ,unstranded= 2)  
    ### Predefine column types
    coltypes <- list(col_character(), col_integer(), col_integer(), col_integer())
    stardirs <- list.files(rootdir, full.names = FALSE)
    mycntfiles <- getSTARcntfnames(level, rootdir, "ReadsPerGene.out.tab", recursive = TRUE)
    out <- foreach(myfile = iter(mycntfiles, by ="row"), .combine = Combine)%do%{
      myfname <- myfile[["myfname"]]
      mylibname <-  myfile[["mylibname"]]
      cat("Read in",mylibname,"\n")
      readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>%
        dplyr::select(1, col) %>%
        dplyr::rename_at(vars(names(.)),~c("gene", mylibname))
      
    }
    if(sum(is.na(out))!= 0)stop("NA exists!")
    out <- out %>%
        as.data.frame()%>%
        tibble::column_to_rownames(var = "gene")
    
    return(out)
}

