getSTARcntfnames <- function(mylevel, rootdir, mypattern, recursive, emptylab = "emptylib", duplab = function(x){make.unique(x, sep = "-")}) {
  myfiles <- normalizePath(list.files(path = rootdir,
                                      pattern = mypattern,
                                      include.dirs = TRUE,
                                      full.names = TRUE,
                                      recursive = recursive))
  ### Use the file basename to get library name
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

library(rtracklayer)
library(tidyverse)

gtffile <- "/mnt/data1/Annotation/Gencode/08292019/Human/gencode.v31.primary_assembly.annotation.gtf"

rtracklayer::import(gtffile) %>%
  (tibble::as_tibble) %>%
  subset(type == "gene") %>%
  mutate(eid = str_remove(gene_id,"\\.[0-9]")) ->
  gtfanno



gtfanno  %>%
  dplyr::select(seqnames, eid, gene_name, gene_id, gene_type) %>%
  dplyr::full_join(out[-(1:4),], by = c("gene_id"="gene")) %>%
  dplyr::filter(seqnames %in% paste0("chr",c(1:22, "X", "Y", "M"))) %>%
  dplyr::filter(gene_type %in% c("protein_coding","lncRNA")) %>%
  dplyr::rename(chr = seqnames) ->
  out1


quickDESeq2obj <- function(starcnt, genelab, libcols) {
  cntdat <- data.frame(starcnt[, libcols], row.names = starcnt[[genelab]])
  coldat <- S4Vectors::DataFrame(idx = seq_len(length(libcols)), row.names = colnames(starcnt)[libcols])
  out <- DESeqDataSetFromMatrix(cntdat, colData = coldat, design = ~1)
  annodat <- S4Vectors::DataFrame(starcnt[,-libcols], row.names = starcnt[[genelab]])
# revise to use granges instead
  mcols(out) <- cbind(mcols(out), annodat)
  if(!identical(starcnt[[genelab]], rownames(mcols(out))))
    out <- NULL
  return(out)
}
deobj <- quickDESeq2obj(out1, "eid", 6:10)
deobj <-estimateSizeFactors(deobj)

deobj <-estimateDispersions(deobj) 
