# --RNAqc--
#' Create new RNAqc class
#' Inherit directly from SummarizedExperiment class 
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment

.RNAqc <- setClass("RNAqc",
         slots  = representation(
           picard = "DataFrame"
           ),
         contains = "SummarizedExperiment",
         prototype = prototype(assays = Assays(SimpleList(counts=matrix(0))),
                               colData = DataFrame(0)))

#'@export
#'@importFrom SummarizedExperiment SummarizedExperiment
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
  se <- SummarizedExperiment(assays = list(counts = counts),colData = colData,...)
  .RNAqc(se,picard = picard)
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


#-- Modify show methods

#' @export
#' @importMethodsFrom methods SummarizedExperiment
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
#' @exportMethods coerce

setAs("SummarizedExperiment", "RNAqc", function(from) {
  new("RNAqc", from, 
      picard=DataFrame(row.names = colnames(from)), 
      )
})

# which works as expected:
# se <- SummarizedExperiment(matrix(rpois(100, lambda=1), ncol=5))
# as(se, "QNAqc")





#### Reference: http://127.0.0.1:27386/library/SummarizedExperiment/doc/Extensions.html