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

# Quick DESeq2 data object

mycnts <- DESeqDataSetFromMatrix(data.frame(out[-seq_len(4),-1], row.names = out[["gene"]][-seq_len(4)]),
                                 colData= S4Vectors::DataFrame(idx = seq_len(ncol(out)-1),
                                                               row.names = colnames(out)[-1]),
                                 design = ~1)


#' Import GTF anotation
library(rtracklayer)

gtffile <- "/mnt/data1/Annotation/Gencode/08292019/Human/gencode.v31.primary_assembly.annotation.gtf"

#' Only keep "genes"
rtracklayer::import(gtffile) %>%
  dplyr::filter(type == "gene") ->
  gtfdf


quickDESeq2obj <- function(starcnt, rm4, gidx = NULL , lidx = NULL, gtfdf = NULL, gtfid = NULL) {
  ## Assumes that the first column contains the gene names
  ## and the remaining are the library names
  if (rm4)
    starcnt <- starcnt[-(1:4),]
  if (is.null(gidx))
    gidx <- 1
  if (is.null(lidx))
    lidx <- 2:ncol(starcnt)
  cntdat <- data.frame(starcnt[, lidx], row.names = starcnt[[gidx]])
  coldat <- S4Vectors::DataFrame(idx = lidx, row.names = colnames(starcnt)[lidx])
  out <- DESeqDataSetFromMatrix(cntdat, colData = coldat, design = ~1)
  ## Fix this! You have to feed this using granges
  if (!is.null(gtfdf) & !is.null(gtfid)) {
    if (identical(rownames(out), gtfdf[[gtfid]]) & gtfid %in% colnames(gtfdf)) {
      rownames(gtfdf) <- gtfdf[[gtfid]]
      mcols(out) <- gtfdf
    }
  }
  out
}

deobj <- quickDESeq2obj(out,rm4 = TRUE)
deobj <- quickDESeq2obj(out,rm4 = TRUE, gtfdf = gtfdf, gtfid = "gene_id")

deobj <-estimateSizeFactors(deobj)
deobj <-estimateDispersions(deobj) 

gtfanno  %>%
  dplyr::select(seqnames, eid, gene_name, gene_id, gene_type) %>%
  dplyr::full_join(out[-(1:4),], by = c("gene_id"="gene")) %>%
  dplyr::filter(seqnames %in% paste0("chr",c(1:22, "X", "Y", "M"))) %>%
  dplyr::filter(gene_type %in% c("protein_coding","lncRNA")) %>%
  dplyr::rename(chr = seqnames) ->
  out1
