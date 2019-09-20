# Helper Function---------------------------------------------------------------------
#' sumHTSeq
#' Helper function that makes adjustments to output from getSTARcounts to be ready for plotting
#' @param counts Counts table from STAR alignment
#' @return A dataframe ready for plotting
#' 


sumHTSeq <- function(counts){
  sid <- as.character(colnames(counts))
  tmp <- data.frame(t(counts[1:4, sid]),
                    N_gene = colSums(counts[-1:-4, sid]))[, 5:1]
  colnames(tmp) <- gsub("N_", "Prop_", colnames(tmp))
  cols2 <- data.frame(Depth = formatC(rowSums(tmp), 
                                      format = "d", 
                                      big.mark = ","))
  cols3 <- round(prop.table(as.matrix(tmp), margin = 1) * 100, 
                 digits = 2)
  summ <- as.data.frame(cbind(cols2, cols3)) 
  return(summ)
}


# Export Function ----------------------------------------------------------------------------
#' plotAlignment
#' Plot stacked barplot for bases counting proportion for STAR and PICARD
#' @param counts Summary counts table from getSTARcounts or getPIcounts
#' @param type Counting reads or counting bases, select one from \code{c("HTSeq","Picard")}
#' @param textsize Annotation size on the plot
#' @return A ggplot type plot
#' @export
#' 


plotAlignment <- function(counts, type = "HTSeq",textsize = 2,...){
  if(type == "HTSeq"){
      summtable <- sumHTSeq(counts)
      plotdata <- summtable %>%
        rownames_to_column(var = "sid") %>%
        melt(measure.vars = grep("Prop_", colnames(.), value = T)) %>%
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
    
    ggplot(plotdata, aes_string(x = "sid", y = "value", fill = "variable")) +
      scale_x_discrete(limits = unique(plotdata$sid)) + 
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      geom_bar(stat="identity", color = "black", size =0.3) + 
      scale_fill_manual(values=colors.alignment) +
      labs(y = "Percentage", x = "Sample", fill = "Type") +
      geom_text(aes(label = paste0(value,"%")), 
                position = position_stack(vjust = 0.5), size = textsize)
  }
  else if(type == "Picard"){
   
     plotdata <- counts %>%
      rownames_to_column(var = "sid") %>%
      melt(id.vars = c("sid","Total Bases")) %>%
      mutate(value = as.numeric(value))
    
      colors.alignment <- c("Coding+UTR" = "#E1FFFF",
                          "Intronic" = "#FFF8DD",
                          "Intergenic" = "#FEE4E1", 
                          "Ribosomal" = "#9df19d",
                          "Unaligned" = "#E6E6F9")
    ggplot(plotdata, aes_string(x = "sid", y = "value", fill = "variable")) +
      scale_x_discrete(limits = unique(plotdata$sid)) + 
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      geom_bar(stat="identity", color = "black", size =0.3) + 
      scale_fill_manual(values=colors.alignment) +
      labs(y = "Percentage", x = "Sample", fill = "Type") +
      geom_text(aes(label = paste0(value,"%")), 
                position = position_stack(vjust = 0.5), size = textsize)
  }
}