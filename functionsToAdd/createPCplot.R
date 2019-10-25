#' get_legend
#' Helper function for getting legend
#' @param gplot ggplot type of plot
#' @return
#'

get_legend<-function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



#' starPCA
#' Function for calculating PCA 
#' @param object an RNAqc object
#' @param ntop use top genes to calculate PCA
#' @param vars grouping variable
#' @return
#'
starPCA <- function(object,ntop = 1000){
  if(!all(colnames(object) == colnames(counts(object)))){stop("Sample names not in same order!")}
  vsd <- vst(object = object)
  rv <- matrixStats::rowVars(assay(vsd))
  idx <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(vsd)[idx, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  dat_pca <- as.data.frame(colData(vsd)) %>%
    rownames_to_column(var = "label") %>%
    mutate(pc1 = pca$x[, 1],
           pc2 = pca$x[, 2],
           pc3 = pca$x[, 3])
  return(list(dat_pca,percentVar))
}


#' PCplot
#' Create PCA plot 
#' @param df input dataframe
#' @param vars grouping variables
#' @param labels sample id
#' @param pc principal components to compare
#' @return 
#' 
#' 

PCplot <- function(df, vars, labels, pc, colors.group, labels.group.pca, shapes.group){
  
  pc_1 = paste0("pc",pc)[1]
  pc_2 = paste0("pc",pc)[2]
  g <- ggplot(df, aes_string(x = pc_1, y = pc_2,
                            color=vars,
                            shape=vars,
                            label = "label")) +
    geom_point(size = 4, alpha=1) +
    # geom_text()+
    labs( x = xl, y = yl,
          shape = "", color="") +
    scale_color_manual(values=colors.group,
                       labels = labels.group.pca) +
    scale_shape_manual(values=shapes.group,
                       labels = labels.group.pca) +
    geom_vline(xintercept=0, linetype = "dashed") +
    geom_hline(yintercept=0, linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "bottom")
   if(labels == TRUE){
     g <- g + geom_text_repel(segment.size = 0.2, segment.color = "grey50", direction =
                                "both",force = 5)
       
   }
  return(g)
}

#' createPCplot
#' This function takes the meta data and star counts 
#' use DESeq2 for estimating size factors before considering any experiment designs.
#' It requires the column names of the input count data matrix and the row names of the metadata have the same order
#' @import ggrepel
#' @import RColorBrewer
#' @param cts star counts matrix
#' @param coldata meta information
#' @param vars grouping variables
#' @param labels sample ids
#' @return 
#' @export


createPCplot <- function(object, var, labels = FALSE){
  out <- starPCA(object)
  pcadf <- out[[1]]
  percentVar <- out[[2]]
  xl<-paste0("PC1", ": ", round(percentVar[1] * 100), "% Variance")
  yl<-paste0("PC2", ": ", round(percentVar[2] * 100), "% Variance")
  zl<-paste0("PC3", ": ", round(percentVar[3] * 100), "% Variance")
  
  myColor =  brewer.pal(n = length(unique(pcadf[[var]])), name = "Paired")
  
  group = factor(pcadf[[var]])
  colors.group <- myColor[factor(levels(group))]
  labels.group.pca = levels(group)
  shapes.group = factor(levels(group))
  
  
  g1 <- PCplot(df = pcadf, vars = var, labels = labels, pc = c(1,2),colors.group = colors.group, labels.group.pca = labels.group.pca, shapes.group = shapes.group)
  g2 <- PCplot(df = pcadf, vars = var, labels = labels, pc = c(1,3),colors.group = colors.group, labels.group.pca = labels.group.pca, shapes.group = shapes.group)
  g3 <- PCplot(df = pcadf, vars = var, labels = labels, pc = c(2,3),colors.group = colors.group, labels.group.pca = labels.group.pca, shapes.group = shapes.group)
  
  mylegend<-get_legend(g1)
  final <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                           g2 + theme(legend.position="none"),
                           g3 + theme(legend.position="none"),
                           nrow=1,top = "Principal Components Analysis"),
               mylegend, nrow=2,heights=c(10, 1))
  
  
  return(final)
  
}




