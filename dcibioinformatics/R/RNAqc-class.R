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
                     Nmap = "matrix",
                     normcounts = "matrix"
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
#'@param anno path to the annotation file, first two columns must be gene id and gene symbols
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
    colnames(anno)[1:2] <- c("gene_id","gene_name")
    # anno %>% transmute( gene_id = Geneid, id_nover = gsub("\\.\\d+", "", Geneid),
    #                     gene_name = GeneSymbol, type = Class) -> anno
    rownames(anno) <- anno$gene_id
    if(!all(row.names(anno) == row.names(rcounts))) anno <-anno[order(row.names(rcounts)),]
    mcols(dds) <- DataFrame(mcols(dds),anno[rownames(dds),])
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

# Calculate normalized expression using vst() beforehand and store the results
#' @export
#' @rdname RNAqc-class
#' @param obj a RNAqc object

normMatrix <- function(obj){
  out <- vst(obj)
  normAssay(obj) <- assay(out)
  colData(obj) <- colData(out)
  return(obj)
}



#-- Add gtf Annotations
#' @export
#' @rdname RNAqc-class
#' @param obj a RNAqc object
#' @param gtfobj a GRange object, typically import from rtracklayer
#' @param id gene ids or transcript ids
#' @param type regions from gtf to match
addGTF <- function(obj, gtfobj,id,type){
  if (class(gtfobj) == "GRanges" & !is.null(id)) {
    gtfobj <- gtfobj[gtfobj$type == type, ]
    if(identical(sort(rownames(obj)), sort(unique(mcols(gtfobj)[[id]])))){
      rowRanges(obj) <- gtfobj
      rownames(obj) <- mcols(obj)$gene_id
    }
  }
  else{
    stop("gtf must be a GRanges object or gene id column is not given!")
  }
  return(obj)
}
#-- Accessor Methods
#' Methods for RNAqc object
#' @export
#' @rdname RNAqc-methods
setGeneric("normAssay", function(object,...)standardGeneric("normAssay"))


#' @export
#' @param object an RNAqc instance
#' @rdname RNAqc-methods
setMethod("normAssay",signature = "RNAqc",function(object){
  out <- object@normcounts
  return(out)
})


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
  coolcat("gene symbols(%d): %s\n", mcols(object)$gene_name)
  
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



#' @export
#' @rdname RNAqc-methods
setGeneric("normAssay<-", function(object,...,value)standardGeneric("normAssay<-"))

#' @export
#' @rdname RNAqc-methods
setReplaceMethod("normAssay","RNAqc",function(object,value){
  object@normcounts <- value
  validObject(object)
  return(object)
})
#-- Enable Subsetting operations

#' @export
#' @rdname RNAqc-methods
#' 
setMethod("[","RNAqc",function(x, i, j, drop = TRUE){
  picard <- piData(x)
  normcounts <- normAssay(x)
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        j, colnames(x), fmt
      )
    }
    j <- as.vector(j)
    picard <- picard[j,]
    normcounts <- normcounts[,j]
    
  }
  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, picard = picard, normcounts = normcounts,check=FALSE)
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
