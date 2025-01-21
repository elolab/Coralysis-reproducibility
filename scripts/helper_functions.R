#-------------------------------Helper script----------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Date: 20/01/2025
# Last update: 20/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
NormalizeData <- function(sce) {
  
  # 'NormalizeData()' function applies the basic Seurat normalization to 
  #a SingleCellExperiment object with a 'counts' assay. Normalized data 
  #is saved in the 'logcounts' assay.
  
  logcounts(sce) <- apply(counts(sce), 2, function(x) {
    log1p(x/sum(x)*10000)
  }) # log1p normalization w/ 10K scaling factor
  logcounts(sce) <- as(logcounts(sce), "sparseMatrix")
  return(sce)
}
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Customize original ILoReg 'GeneScatterPlot' function: https://github.com/elolab/ILoReg/blob/250a4facdaa66977a529fee02b4a8927364b43df/R/CoreMethods.R#L1006
CustomGeneScatterPlot <- function (object, genes = "", return.plot = FALSE, dim.reduction.type = "tsne", 
                                   point.size = 0.7, title = "", plot.expressing.cells.last = FALSE, 
                                   nrow = NULL, ncol = NULL) {
  if (dim.reduction.type == "umap") {
    two.dim.data <- reducedDim(object, "UMAP")
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  }
  else if (dim.reduction.type == "tsne") {
    two.dim.data <- reducedDim(object, "TSNE")
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  }
  else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }
  if (length(genes) == 1) {
    df <- as.data.frame(two.dim.data)
    if (!(genes %in% rownames(object))) {
      stop("invalid gene name")
    }
    color.by <- logcounts(object)[genes, ]
    df$group <- color.by
    colnames(df) <- c("dim1", "dim2", "group")
    if (title == "") {
      if (plot.expressing.cells.last) {
        df <- df[order(df$group, decreasing = FALSE), 
        ]
      }
      p <- ggplot(df, aes_string(x = "dim1", y = "dim2")) + 
        geom_point(size = point.size, aes_string(color = "group"), stroke = 0) + 
        scale_colour_gradient2(low = scales::muted("red"), mid = "lightgrey", 
                               high = "blue", name = genes) + 
        theme(axis.title = element_blank(), panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), panel.background = element_blank(), 
              axis.line = element_line(colour = "black"), 
              legend.title = element_blank(), 
              legend.key.size = unit(0.3, "cm"), 
              plot.margin = unit(c(0,0,0,0), "cm"), 
              legend.margin = margin(c(0,5,0,-5)), 
              plot.title = element_text(margin = margin(c(2.5,0,-1,0)), size = 10))
    }
    else {
      p <- ggplot(df, aes_string(x = "dim1", y = "dim2")) + 
        geom_point(size = point.size, aes_string(color = "group"), stroke = 0) + 
        scale_colour_gradient2(low = scales::muted("red"), mid = "lightgrey", 
                               high = "blue", name = genes) + 
        theme(axis.title = element_blank(), panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), panel.background = element_blank(), 
              axis.line = element_line(colour = "black")) + 
        ggtitle(paste(title, gene, sep = " : ")) + theme(legend.title = element_blank(), 
                                                         legend.key.size = unit(0.3, "cm"), 
                                                         plot.margin = unit(c(0,0,0,0), "cm"), 
                                                         legend.margin = margin(c(0,5,0,-5)), 
                                                         plot.title = element_text(hjust = 0.5, 
                                                                                   margin = margin(c(2.5,0,-1,0)), size = 10))
    }
    if (return.plot) {
      return(p)
    }
    else {
      print(p)
    }
  }
  else {
    plot_list <- list()
    for (gene in genes) {
      df <- as.data.frame(two.dim.data)
      if (!(gene %in% rownames(object))) {
        stop(paste0("invalid gene name: ", gene))
      }
      color.by <- logcounts(object)[gene, ]
      df$group <- color.by
      colnames(df) <- c("dim1", "dim2", "group")
      if (plot.expressing.cells.last) {
        df <- df[order(df$group, decreasing = FALSE), 
        ]
      }
      if (title == "") {
        p <- ggplot(df, aes_string(x = "dim1", y = "dim2")) + 
          geom_point(size = point.size, aes_string(color = "group"), stroke = 0) + 
          scale_colour_gradient2(low = scales::muted("red"), 
                                 mid = "lightgrey", high = "blue", name = gene) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black"), 
                axis.title=element_blank(), 
                legend.title = element_blank(), 
                legend.key.size = unit(0.3, "cm"), 
                plot.margin = unit(c(0,0,0,0), "cm"), 
                legend.margin = margin(c(0,5,0,-5)), 
                plot.title = element_text(margin = margin(c(2.5,0,-1,0)), size = 10))
      }
      else {
        p <- ggplot(df, aes_string(x = "dim1", y = "dim2")) + 
          geom_point(size = point.size, aes_string(color = "group"), stroke = 0) + 
          scale_colour_gradient2(low = scales::muted("red"), 
                                 mid = "lightgrey", high = "blue", name = gene) + 
          theme(panel.grid.major = element_blank(), 
                axis.title=element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black")) + 
          ggtitle(paste(title, gene, sep = " : ")) + theme(legend.title = element_blank(), 
                                                           legend.key.size = unit(0.3, "cm"), 
                                                           plot.margin = unit(c(0,0,0,0), "cm"), 
                                                           legend.margin = margin(c(0,5,0,-5)), 
                                                           plot.title = element_text(hjust = 0.5, margin = margin(c(2.5,0,-1,0)), size = 10)) 
      }
      plot_list[[gene]] <- p
    }
    p <- cowplot::plot_grid(plotlist = plot_list, align = "hv", nrow = nrow, ncol = ncol)
    if (return.plot) {
      return(p)
    }
    else {
      print(p)
    }
  }
}
#
#------------------------------------------------------------------------------#