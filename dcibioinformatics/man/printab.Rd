% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printab.R
\name{printab}
\alias{printab}
\title{printab
Print summary table in LaTex format}
\usage{
printab(tab, cap = NULL, top = 20, scale = 0.7, align = NULL,
  grepcols = "pval|padj", ...)
}
\arguments{
\item{tab}{Summary table from getSTARcounts, getSTARsum or getPIcounts}

\item{cap}{Character vector of length 1 or 2 containing the table's caption or title. If length is 2, the second item is the "short caption" used when LaTeX generates a "List of Tables". Set to NULL to suppress the caption. Default value is NULL}

\item{top}{number of rows to display, usually samples}

\item{scale}{parameters used to control the scale of the table}

\item{align}{Character vector of length equal to the number of columns of the resulting table, indicating the alignment of the corresponding columns. Also, "|" may be used to produce vertical lines between columns in LaTeX tables, but these are effectively ignored when considering the required length of the supplied vector. If a character vector of length one is supplied, it is split as strsplit(align, "")[[1]] before processing. Since the row names are printed in the first column, the length of align is one greater than ncol(x) if x is a data.frame. Use "l", "r", and "c" to denote left, right, and center alignment, respectively. Use "p{3cm}" etc. for a LaTeX column of the specified width. For HTML output the "p" alignment is interpreted as "l", ignoring the width request. Default depends on the class of x.}

\item{grepcols}{The columns that needs to be presented in scientific mode}

\item{...}{Other visual options to adjust the latex output. See more information from package \code{xtable}}
}
\value{
output table in LaTex format
}
\description{
printab
Print summary table in LaTex format
}
