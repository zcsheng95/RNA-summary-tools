#### Some low level functions for internal usage
#### Not export
### coolcat function is borrowed from S4Vectors: https://github.com/Bioconductor/S4Vectors/blob/master/R/show-utils.R
#' @importFrom Biobase selectSome
coolcat <- function(fmt, vals=character(), exdent=2, ...)
{
  vals <- ifelse(nzchar(vals), vals, "''")
  lbls <- paste(selectSome(vals), collapse=" ")
  txt <- sprintf(fmt, length(vals), lbls)
  cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}
