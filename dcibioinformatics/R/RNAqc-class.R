# --RNAqc--
#' Create new RNAqc class
#' Inherit directly from DESeqDataSet class 
#' @export
#' @import methods
#' @importClassesFrom DESeq2 DESeqDataSet

.RNAqc <- setClass("RNAqc",
                   slots  = representation(
                     picard = "DataFrame",
                     Nmap = "matrix"
                   ),
                   contains = "DESeqDataSet",
                   prototype = prototype(assays = Assays(SimpleList(counts=matrix(0))),
                                         colData = DataFrame(0)))

#'@export
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom DESeq2 DESeqDataSet
#'

RNAqc <- function(counts,colData,picard = DataFrame(),...){
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
  if (is.null(msg)) {
    TRUE
  }else msg
}

setValidity2("RNAqc", .valid.RNAqc)

#-- Accessor Methods

#' @export
setGeneric("piData",function(object,...)standardGeneric("piData"))

#' @export
setMethod("piData",signature = "RNAqc", function(object){
  out <- object@picard
  return(out)
})

#' @export
setGeneric("Nmap", function(object,...)standardGeneric("Nmap"))

#' @export
setMethod("Nmap", signature = "RNAqc", function(object){
  out <- object@Nmap
  return(out)
})
#-- Modify show methods

#' @export
#' @importMethodsFrom SummarizedExperiment show
#'
#'
setMethod("show", signature = "RNAqc", function(object){
  callNextMethod()
  coolcat("picard names(%d): %s\n", rownames(piData(object)))
})


#-- Setter Methods

#' @export
#' 
setGeneric("piData<-", function(object,...,value)standardGeneric("piData<-"))

#' @export
#'
setReplaceMethod("piData","RNAqc",function(object,value){
  object@picard <- value
  validObject(object)
  return(object)
})

#' @export
#' 
setGeneric("Nmap<-", function(object,...,value)standardGeneric("Nmap<-"))

#' @export
#' 
setReplaceMethod("Nmap","RNAqc", function(object,value){
  object@Nmap <- value
  validObject(object)
  return(object)
})

#-- Enable Subsetting operations

#' @export
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
#' @exportMethod coerce

setAs("DESeqDataSet", "RNAqc", function(from) {
  new("RNAqc", from, 
      picard=DataFrame(row.names = colnames(from)), 
  )
})

# which works as expected:
# se <- SummarizedExperiment(matrix(rpois(100, lambda=1), ncol=5))
# as(se, "QNAqc")





#### Reference: http://127.0.0.1:27386/library/SummarizedExperiment/doc/Extensions.html
