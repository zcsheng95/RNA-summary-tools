# Helper function -------------------------------------------------------------

#' picardCombine
#' Row bind of picard summary for df1 and df2
#' @param df1 picard summary for sample 1
#' @param df2 picard summary for sample 2
#' @return A joint summary dataframe from df1 and df2

picardCombine <- function(df1, df2){
  rbind(df1, df2)
}



#' myPidfile
#' get picard \code{.metrics} files from root directory
#' @param rootdir Directory that stores all picard \code{.metrics} file
#' @return A list contains all the \code{.metrics} files

myPidfile <- function(rootdir){
  list.files(rootdir,pattern = "*.metrics",full.names = TRUE, recursive = TRUE)
}


# Export function ------------------------------------------------------------------------
utils::globalVariables(c("pidfile","PF_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES","RIBOSOMAL_BASES",
                         "Coding+UTR","Intronic","Intergenic","Ribosomal"))

#' getPIcounts
#' Extract counting bases summary from PICARD output
#' @import dplyr
#' @import magrittr
#' @param rootdir Directory that stores all picard \code{.metrics} file
#' @return joint table for counting bases from PICARD output
#' @export


getPIcounts <- function(rootdir) {
  pidfiles <- myPidfile(rootdir)  
  out <- foreach(pidfile = pidfiles, .combine = picardCombine)%do%{
      gene <- gsub(".metrics","",basename(pidfile))
      readr::read_tsv(pidfile, col_names = TRUE, skip = 6, n_max = 1,col_types = cols()) %>%
        dplyr::transmute(`Total Bases` = format(PF_BASES,
                                         big.mark = ",",
                                         scientific = F),
                  `Coding+UTR` = CODING_BASES + UTR_BASES,
                  Intronic = INTRONIC_BASES,
                  Intergenic = INTERGENIC_BASES, 
                  Ribosomal = RIBOSOMAL_BASES,
                  Unaligned = PF_BASES - `Coding+UTR` - Intronic - Intergenic - Ribosomal)%>%
        magrittr::set_rownames(gene)
  }
  out[, 2:6] <- round(prop.table(out[, 2:6]) * 100, 2)
  return(out)
  }