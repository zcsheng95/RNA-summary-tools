
# Export function -------------------------------------------------------------------
utils::globalVariables(c("stardir","."))

#' getSTARcounts
#' Extract strand-specific STAR counts for each sample and then merge all counts into a single table
#' Make sure the files and library names are paired 
#' @import dplyr
#' @import foreach
#' @import readr
#' @import tibble
#' @importFrom iterators iter
#' @param files The list of files to read in
#' @param libnames Sample names, make sure they are paired with files 
#' @param strand The strand-specific protocol used in RNA sequencing. Choose one from \code{c("first", "second", "unstranded")}
#' @param verbose True or False. Show read in progress 
#' @return A strand-specific STAR count summary table
#' @export 

getSTARcounts <- function(files, libnames, strand, verbose = TRUE) {
  if(!strand %in% c("first","second","unstranded")) stop("Missing strand argument! You must specify strand-specific protocol from c('first','second','unstranded')")
  col = switch(strand, second = 4,first= 3 ,unstranded= 2)  
  ### Predefine column types
  coltypes <- list(col_character(), col_integer(), col_integer(), col_integer())
  mycntfiles <- tibble(myfname = files, mylibname = libnames)
  out <- foreach(myfile = iter(mycntfiles, by ="row"), .combine = Combine)%do%{
    myfname <- myfile[["myfname"]]
    mylibname <-  myfile[["mylibname"]]
    if(verbose == TRUE){
    cat("Read in",mylibname,"\n")
    }
    readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>%
      dplyr::select(1, col) %>%
      dplyr::rename_at(vars(names(.)),~c("gene", mylibname))
    
  }
  if(sum(is.na(out))!= 0)stop("NA exists! Check if the same reference file is used!")
  out <- out %>%
    as.data.frame()%>%
    tibble::column_to_rownames(var = "gene")
  
  return(out)
}

