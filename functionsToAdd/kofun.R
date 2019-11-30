#' This is a convenience function for extracting the STAR count files
#' It will return a tibble with two colums:
#' 1. Absolute file name for count file
#' 2. library name
#' Perhaps, it is best to extend it to contain  absolute file name
#' for the final log file

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


library(foreach)
library(iterators)
library(tidyverse)





countCombine <- function(df1, df2) {
    dplyr::full_join(df1, df2, by="gene")
}  
coltypes <- list(col_character(), col_integer(), col_integer(), col_integer())

mycntfiles <- getSTARcntfnames(0, ".", "ReadsPerGene.out.tab", recursive = FALSE)

  
out <- foreach(myfile = iter(mycntfiles, by ="row"), .combine = countCombine)%do%{
  myfname <- myfile[["myfname"]]
  mylibname <-  myfile[["mylibname"]]
  print(mylibname)
  readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>%
    dplyr::select(1, 4) %>%
      dplyr::rename_at(vars(names(.)),~c("gene", mylibname))
  
}



#' Import GTF anotation
library(rtracklayer)
gtffile <- "/mnt/data1/Annotation/Gencode/08292019/Human/gencode.v31.primary_assembly.annotation.gtf"

#' Only keep "genes"
rtracklayer::import(gtffile) -> gtfobj
#' as.data.frame(gtfobj) -> gtfdf



quickDESeq2obj <- function(starcnt, rm4, gidx = NULL , lidx = NULL, gtfobj = NULL, geneid = NULL) {
  #' This is a a convenience function for creating a DESeq2 data object. By default it assumes
  #' assumes that the first column contains the gene names and the remaining are the library names
  #' It also provides the option for adding annotation data as a GRanges object.

  ## Should the first for rows be removed?
  if (rm4)
    starcnt <- starcnt[-(1:4),]
  ## By default assume that the gene id is in the first column
  if (is.null(gidx))
    gidx <- 1
  ## By default assume that the libraries are in columns 2 ... ncol(starcnt)
  if (is.null(lidx))
    lidx <- 2:ncol(starcnt)
  ## Create a data.frame of the counts using the gene ids as row names
  cntdat <- data.frame(starcnt[, lidx], row.names = starcnt[[gidx]])
  ## Create a dummy colData DataFrame
  coldat <- S4Vectors::DataFrame(idx = lidx, row.names = colnames(starcnt[,lidx]))
  # Create a barebones DESeq2 data object
  dedat <- DESeqDataSetFromMatrix(countData = cntdat, colData = coldat, design = ~1)
  ## To add annotation data first check that a Granges object has been provided
  ## and that a geneid column has been given
  if (class(gtfobj) == "GRanges" & !is.null(geneid)) {
    ## Confirm that the gene ids in the data object are identical to those in the
    ## gtf object (same length/order) and that the geneid slot is present
    if (identical(rownames(dedat), mcols(gtfobj)[[geneid]]) & geneid %in% colnames(mcols(gtfobj))) {
      rowRanges(dedat) <- gtfobj
    }
  }
  return(dedat)
}

if(FALSE) {
  deobj <- quickDESeq2obj(out, rm4 = TRUE)
  deobj <- quickDESeq2obj(out, rm4 = TRUE, gidx = 1, lidx = 2:5)
  deobj <- quickDESeq2obj(out, rm4 = TRUE, gidx = "gene" , lidx = colnames(out)[-1])
  deobj <- quickDESeq2obj(out, rm4 = TRUE, gtfobj = gtfobj[gtfobj$type == "gene", ], geneid="gene_id")
  
  deobj <-estimateSizeFactors(deobj)
  deobj <-estimateDispersions(deobj)
  deobj1 <- deobj[mcols(deobj)$gene_type %in% c("protein_coding", "lncRNA", "Mt_rRNA"),]

  mcols(deobj)[,c("gene_id", "gene_type"), drop = FALSE] %>%
    as_tibble %>%
    full_join(out, by = c("gene_id" = "gene")) %>%
    group_by(gene_type) %>%
    summarize_at(vars(2:6), sum) %>%
    arrange(desc(`pt-2-pre-DUKE`)) %>% data.frame
}

