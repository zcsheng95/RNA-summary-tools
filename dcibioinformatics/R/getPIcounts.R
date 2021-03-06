# Helper function -------------------------------------------------------------

#' picardCombine
#' Row bind of picard summary for df1 and df2
#' @param df1 picard summary for sample 1
#' @param df2 picard summary for sample 2
#' @return A joint summary dataframe from df1 and df2

picardCombine <- function(df1, df2){
  rbind(df1, df2)
}



# Export function ------------------------------------------------------------------------
utils::globalVariables(c("pidfile","PF_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES","RIBOSOMAL_BASES",
                         "Coding+UTR","Intronic","Intergenic","Ribosomal"))

#' getPIcounts
#' Extract counting bases summary from PICARD output
#' @import dplyr
#' @importFrom magrittr set_rownames
#' @param files all picard \code{.metrics} file
#' @param libnames Library Names
#' @return joint table for counting bases from PICARD output
#' @export


getPIcounts <- function(files, libnames) {
  mymetricfiles <- tibble(myfname = files, mylibname = libnames)
  out <- foreach(myfile = iter(mymetricfiles, by ="row"), .combine = picardCombine)%do%{
    pidfile <- myfile[["myfname"]]
    mylibname <-  myfile[["mylibname"]]  
    readr::read_tsv(pidfile, col_names = TRUE, skip = 6, n_max = 1,col_types = cols()) %>%
        dplyr::transmute(`Total Bases` = format(PF_BASES,
                                         big.mark = ",",
                                         scientific = F),
                  `Coding+UTR` = CODING_BASES + UTR_BASES,
                  Intronic = INTRONIC_BASES,
                  Intergenic = INTERGENIC_BASES, 
                  Ribosomal = RIBOSOMAL_BASES,
                  Unaligned = PF_BASES - `Coding+UTR` - Intronic - Intergenic - Ribosomal)%>%
        as.data.frame()%>%
        magrittr::set_rownames(mylibname)
  }
  out[, 2:6] <- round(prop.table(as.matrix(out[,2:6]),margin = 1) * 100, 2)
  return(out)
}
