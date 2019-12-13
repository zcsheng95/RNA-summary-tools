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
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @param object A RNAqc instance
#' @param group grouping variable, default to NULL. If specified, will also return scatter plot for mapping rate.
#' @param textsize Annotation size on the plot
#' @param ... Other visual options used for geom_text
#' @return A list of ggplot.
#' @export
#' 


plotAlignment <- function(object, group = NULL, textsize = 3,...){
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
    theme(axis.text.x=element_text(angle=60, hjust=1))+
    geom_bar(stat="identity", color = "black", size =0.3) + 
    scale_fill_manual(values=colors.alignment) +
    labs(y = "Percentage", x = "Sample", fill = "Type") +
    geom_text(aes(label = paste0(value,"%")),  
              position = position_stack(vjust = 0.5),size = textsize,...)+
    theme_classic()
  
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
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      geom_bar(stat="identity", color = "black", size =0.3) + 
      scale_fill_manual(values=colors.alignment) +
      labs(y = "Percentage", x = "Sample", fill = "Type") +
      geom_text(aes(label = paste0(value,"%")), 
                position = position_stack(vjust = 0.5),size=textsize,...)+
      theme_classic()
  }
  else(g_pi = "No picard outputs provided!")
  
  if(!is.null(group)){
  plot_dat <- summtable %>% 
    rownames_to_column(var = "sid") %>%
    dplyr::mutate(unimap = Prop_gene + Prop_ambiguous+ Prop_noFeature)
  
  g_map <- ggplot(plot_dat, aes_string(x = "sid", y = "unimap",color = as.factor(colData(object)[[group]])))+
    geom_point(size = 4) +
    labs(x = "Samples", y = "Unique Mapping Rate",color = group) +
    scale_x_discrete(limits= plot_dat$sid) +
    theme_classic() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=12, face="bold"),
          plot.title=element_text(face="italic"),
          axis.text.x=element_text(hjust=1))
  return(list(g_star,g_pi,g_map))
  }
  else return(list(g_star,g_pi))
}

#' @export 
#' @rdname plotAlignment

setMethod("plot",
          signature = signature(x = "RNAqc",y = "missing"),
          definition =function(x,y,...)plotAlignment(x,...))
