getSTARcntfnames <- function(mylevel, rootdir, mypattern, recursive, emptylab = "emptylib", 
    duplab = function(x) {
        make.unique(x, sep = "-")
    }) {
    myfiles <- normalizePath(list.files(path = rootdir, pattern = mypattern, include.dirs = TRUE, 
        full.names = TRUE, recursive = recursive))
    if (mylevel == 0) {
        ## Use the file basename to get library name Get basename after removing pattern
        mylibname <- str_remove(basename(myfiles), mypattern)
        ## If filename equals pattern (no prefix) set it to emptylab
        mylibname <- dplyr::case_when(mylibname == "" ~ emptylab, TRUE ~ mylibname)
    } else {
        ## Use the directory structure to get library names For this format
        ## foo/a/star/ReadsPerGene.out.tab , foo/b/star/ReadsPerGene.out.tab Use mylevel
        ## =2 to get a and b For this format foo/star/a/ReadsPerGene.out.tab ,
        ## foo/star/b/ReadsPerGene.out.tab Use mylevel =1 to get a and b
        mylibname <- sapply(str_split(dirname(myfiles), "/"), dplyr::nth, n = -mylevel)
    }
    tibble::tibble(myfname = myfiles, mylibname = duplab(mylibname))
}


library(foreach)
library(iterators)
library(tidyverse)





countCombine <- function(df1, df2) {
    dplyr::full_join(df1, df2, by = "gene")
}
coltypes <- list(col_character(), col_integer(), col_integer(), col_integer())

mycntfiles <- getSTARcntfnames(0, ".", "ReadsPerGene.out.tab", recursive = FALSE)


out <- foreach(myfile = iter(mycntfiles, by = "row"), .combine = countCombine) %do% 
    {
        myfname <- myfile[["myfname"]]
        mylibname <- myfile[["mylibname"]]
        print(mylibname)
        readr::read_tsv(myfname, col_names = FALSE, col_types = coltypes) %>% dplyr::select(1, 
            4) %>% dplyr::rename_at(vars(names(.)), ~c("gene", mylibname))
        
    }

# Quick DESeq2 data object

mycnts <- DESeqDataSetFromMatrix(data.frame(out[-seq_len(4), -1], row.names = out[["gene"]][-seq_len(4)]), 
    colData = S4Vectors::DataFrame(idx = seq_len(ncol(out) - 1), row.names = colnames(out)[-1]), 
    design = ~1)


