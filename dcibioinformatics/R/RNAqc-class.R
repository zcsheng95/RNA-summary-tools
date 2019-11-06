# --RNAqc--
#' Create new RNAqc class
#' Inherit directly from DESeqDataSet class 
#' @export
#' @import methods
#' @import SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors setValidity2
#' @importClassesFrom DESeq2 DESeqDataSet
#' @rdname RNAqc-class
#' @name RNAqc-class

.RNAqc <- setClass("RNAqc",
                   slots  = representation(
                     picard = "DataFrame",
                     Nmap = "matrix"
                   ),
                   contains = "DESeqDataSet",
                   prototype = prototype(assays = Assays(SimpleList(counts=matrix(0))),
                                         colData = DataFrame(0)))
#' Constructor for RNAqc Class
#' 
#'@export
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom DESeq2 DESeqDataSet
#'@importFrom S4Vectors mcols
#'@param counts STAR counts matrix
#'@param colData manifest information
#'@param picard picard summary table
#'@param anno path to the annotation file 
#'@param ... for potential extra arguments in constructing SummarizedExperiment
#'@rdname RNAqc-class
#'@name RNAqc-constructor
#'@docType methods

RNAqc <- function(counts,colData,picard = DataFrame(), anno = NULL,...){
  if(is(colData,"data.frame")){
    colData = as(colData, "DataFrame")
  }
  if(is(picard,"data.frame")){
    picard = as(picard, "DataFrame")
  }
  if(!is(counts,"matrix")){
    counts = as(counts, "matrix")
  }
  Nmapp = counts[1:4,]
  rcounts = counts[-1:-4,]
  se <- SummarizedExperiment(assays = list(counts = rcounts),colData = colData,...)
  dds <- DESeqDataSet(se, design = ~1)
  if(!is.null(anno)){
    annot <- read.table(anno, sep = '\t', header = T, stringsAsFactors = F)
    annot %>% transmute(ens_id_ver = Geneid, ens_id = gsub("\\.\\d+", "", Geneid),
                        symbol = GeneSymbol, genetype = Class) -> annot
    rownames(annot) <- annot$ens_id_ver
    if(!all(row.names(annot) == row.names(rcounts))) annot <-annot[order(row.names(rcounts)),]
    mcols(dds) <- DataFrame(mcols(dds),annot[rownames(dds),])
  }
  .RNAqc(dds,picard = picard,Nmap = Nmapp)
}


.valid.RNAqc <- function(object){
  msg <- NULL
  
  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }
  
  if (min(assay(object)) < 0) {
    msg <- c(msg, "negative values in 'counts'")
  }
  
  if (!is(colData(object),"DataFrame")){
    msg <- c(msg, "colData must be a DataFrame")
  }
  if(!is(piData(object),"DataFrame")){
    msg <- c(msg, "picard must be a DataFrame")
  }
  if(min(Nmap(object) <0)){
    msg <- c(msg, "Mapping statistics can not be negative")
  }
  if(nrow(mcols(object)) != nrow(counts(object))){
    msg <- c(msg, "Number of genes in the annotation does not match number of genes in the count matrix")
  }
  if (is.null(msg)) {
    TRUE
  }else msg
}

setValidity2("RNAqc", .valid.RNAqc)

#-- Accessor Methods
#' Methods for RNAqc object
#' @export
#' @rdname RNAqc-methods
setGeneric("piData",function(object,...)standardGeneric("piData"))

#' @export
#' @param object An RNAqc instance 
#' @rdname RNAqc-methods
setMethod("piData",signature = "RNAqc", function(object){
  out <- object@picard
  return(out)
})

#' @export
#' @rdname RNAqc-methods
setGeneric("Nmap", function(object,...)standardGeneric("Nmap"))

#' @export
#' @rdname RNAqc-methods
setMethod("Nmap", signature = "RNAqc", function(object){
  out <- object@Nmap
  return(out)
})
#-- Modify show methods

#' @export
#' @importMethodsFrom SummarizedExperiment show
#' @rdname RNAqc-methods
#'
setMethod("show", signature = "RNAqc", function(object){
  callNextMethod()
  coolcat("picard names(%d): %s\n", rownames(piData(object)))
  coolcat("gene symbols(%d): %s\n", mcols(object)$symbol)
  
})


#-- Setter Methods

#' @export
#' @rdname RNAqc-methods
setGeneric("piData<-", function(object,...,value)standardGeneric("piData<-"))

#' @export
#' @rdname RNAqc-methods
setReplaceMethod("piData","RNAqc",function(object,value){
  object@picard <- value
  validObject(object)
  return(object)
})

#' @export
#' @rdname RNAqc-methods
setGeneric("Nmap<-", function(object,...,value)standardGeneric("Nmap<-"))

#' @export
#' @rdname RNAqc-methods
setReplaceMethod("Nmap","RNAqc", function(object,value){
  object@Nmap <- value
  validObject(object)
  return(object)
})

#-- Enable Subsetting operations

#' @export
#' @rdname RNAqc-methods
#' 
setMethod("[","RNAqc",function(x, i, j, drop = TRUE){
  picard <- piData(x)
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        j, colnames(x), fmt
      )
    }
    j <- as.vector(j)
    picard <- picard[j,]
    
  }
  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, picard = picard, check=FALSE)
})


#-- Subsetting Assignment
#' @export
#' @rdname RNAqc-methods
setReplaceMethod("[", c("RNAqc", "ANY", "ANY", "RNAqc"),function(x, i, j, ..., value) {
  picard <- piData(x)
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        j, colnames(x), fmt
      )
    }
    j <- as.vector(j)
    picard[j,] <- piData(value)
  }
  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, picard = picard, check=FALSE)
})



#-- Define Coerce
#' @rdname RNAqc-methods
#' @name coerce
#' @aliases coerce,DESeqDataSet,RNAqc-method
#' @exportMethod coerce


setAs("DESeqDataSet", "RNAqc", function(from) {
  new("RNAqc", from, 
      picard=DataFrame(row.names = colnames(from)), 
  )
})

# which works as expected:
# se <- SummarizedExperiment(matrix(rpois(100, lambda=1), ncol=5))
# as(se, "RNAqc")





#### Reference: http://127.0.0.1:27386/library/SummarizedExperiment/doc/Extensions.html
