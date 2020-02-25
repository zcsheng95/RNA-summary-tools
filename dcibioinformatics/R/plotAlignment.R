# Summarize STAR output---------------------------------------------------------------------
#' summary.RNAqc
#' summary of an RNAqc object to generate counts proportion table
#' @param object An RNAqc instance
#' @return A dataframe ready for plotting
#' 
#' 

summary.RNAqc <- function(object){
  counts <- assay(object)
  sid <- as.character(colnames(counts))
  tmp <- data.frame(t(Nmap(object)[, sid]),
                    N_gene = colSums(counts[, sid]))[, 5:1]
  colnames(tmp) <- gsub("N_", "Prop_", colnames(tmp))
  cols2 <- data.frame(Depth = formatC(rowSums(tmp), 
                                      format = "d", 
                                      big.mark = ","))
  cols3 <- round(prop.table(as.matrix(tmp), margin = 1) * 100, 
                 digits = 2)
  summ <- as.data.frame(cbind(cols2, cols3)) 
  return(summ)
}

#' @export
#' @rdname summary.RNAqc
setMethod("summary","RNAqc",function(object){
  summary.RNAqc(object)
})

# export ----------------------------------------------------------------------------
### Define global variables to avoid warnings in R CMD Check
utils::globalVariables(c("object","variable",".","value","sid","Depth","Total.Bases","Prop_gene","Prop_ambiguous",
                         "Prop_noFeature"))

#' plotAlignment
#' Plot stacked barplot for bases counting proportion for STAR and PICARD
#' @rdname plotAll
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @param object A RNAqc instance
#' @param group grouping variable, default to NULL. If specified, will also return scatter plot for mapping rate.
#' @param textsize Annotation size on the plot
#' @param angle x axis angle
#' @param ... settings for geom_point
#' @return A list of ggplot.
#' 


plotSTAR <- function(object, textsize = 3,text = TRUE,angle = 60,...){
  summtable <- summary(object)
  plotdata <- summtable %>%
    rownames_to_column(var = "sid") %>%
    tidyr::gather(variable, value,-Depth,-sid) %>%
    mutate(variable = factor(variable, 
                             levels = unique(variable),
                             labels = c("Gene", 
                                        "Ambiguous", 
                                        "No Feature", 
                                        "Multimapping", 
                                        "Unmapped")),
           value = as.numeric(value),
           sid = sid)
  
  
  colors.alignment <- c("Gene" = "#0571B0", 
                        "Ambiguous" = "#92C5DE",
                        "No Feature" = "#F7F7F7", 
                        "Multimapping" = "#F4A582",
                        "Unmapped" = "#CA0020")
  g_star <- ggplot2::ggplot(plotdata, aes_string(x = "sid", y = "value", fill = "variable")) +
    scale_x_discrete(limits = unique(plotdata$sid)) + 
    geom_bar(stat="identity", color = "black", size =0.3) + 
    scale_fill_manual(values=colors.alignment) +
    labs(y = "Percentage", x = "Sample", fill = "Type") +
    theme_classic()+
    theme(axis.text.x=element_text(angle=angle, hjust=1))
  if(text == TRUE){
    g_star <- g_star + geom_text(aes(label = paste0(value,"%")),  
              position = position_stack(vjust = 0.5),size = textsize)
  }
  return(g_star)
}

plotPicard <- function(object, textsize = 3,text = TRUE,angle = 60,...){
   if(!identical(piData(object),DataFrame(0))){
    counts <- piData(object)
    plotdata <- counts %>%
      as.data.frame() %>%
      rownames_to_column(var = "sid") %>%
      tidyr::gather(variable,value, -sid, -Total.Bases) %>%
      mutate(value = as.numeric(value))
    
    colors.alignment <- c("Coding+UTR" = "#E1FFFF",
                          "Intronic" = "#FFF8DD",
                          "Intergenic" = "#FEE4E1", 
                          "Ribosomal" = "#9df19d",
                          "Unaligned" = "#E6E6F9")
    g_pi<-ggplot(plotdata, aes_string(x = "sid", y = "value", fill = "variable")) +
      scale_x_discrete(limits = unique(plotdata$sid)) + 
      geom_bar(stat="identity", color = "black", size =0.3) + 
      scale_fill_manual(values=colors.alignment) +
      labs(y = "Percentage", x = "Sample", fill = "Type") +
      theme_classic()+
      theme(axis.text.x=element_text(angle=angle, hjust=1))
    if(text == TRUE){
      g_pi <- g_pi + geom_text(aes(label = paste0(value,"%")),  
                                   position = position_stack(vjust = 0.5),size = textsize)
    }
    return(g_pi)
  }
  else(print("No picard outputs provided!"))
}

plotMapp <- function(object, group = NULL, textsize = 3,angle = 60,...){
  if(!is.null(group)){
    summtable <- summary(object)
    plot_dat <- summtable %>% 
      rownames_to_column(var = "sid") %>%
      dplyr::mutate(unimap = Prop_gene + Prop_ambiguous+ Prop_noFeature)
    plot_dat$sid <- factor(plot_dat$sid, levels = plot_dat$sid)
    g_map <- ggplot(plot_dat, aes_string(x = "sid", y = "unimap",color = as.factor(colData(object)[[group]])))+
      geom_point(...) +
      labs(x = "Samples", y = "Unique Mapping Rate",color = group) +
      scale_x_discrete(limits= plot_dat$sid) +
      theme_classic() +
      theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12, face="bold"),
          plot.title=element_text(face="italic"),
          axis.text.x=element_text(hjust=1, angle = angle))
  return(g_map)
  }
  else stop("Grouping variable must be provided to visualize mapping rate!")
}


# --- plotAll functions
#' @param object an RNAqc instance
#' @param type type of plot 
#' @param ... additional settings
#' @export

plotAll <- function(object,type,...){
  if(type == "STAR"){
    plotSTAR(object,...)
  }
  else if(type == "PICARD"){
    plotPicard(object,...)
  }
  else if(type == "MAPR"){
    plotMapp(object,...)
  }
  else if(type == "PCA"){
    createPCplot(object,...)
  }
  else if(type == "SizeFactor"){
    plotSizefactor(object,...)
  }
  else if(type == "Expression"){
    checkExpression(object,...)
  }
}

#' @export 
#' @rdname plotAll

setMethod("plot",
          signature = signature(x = "RNAqc",y = "missing"),
          definition =function(x,y,...)plotAll(x,...))
