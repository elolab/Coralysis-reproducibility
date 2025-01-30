#-------------------------------Helper script----------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Date: 20/01/2025
# Last update: 30/01/2025
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

#------------------------------------------------------------------------------#
#
plot_umap <- function(x, by = "batch", size = 0.01, use.color = NULL, legend.nrow = 8,
                      seed = 123, no_legend = TRUE, title = NULL, 
                      draw.line = FALSE) {
  if (by == "batch") group <- colnames(x)[3] else group <- colnames(x)[2]  
  if (is.null(use.color)) {
    color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
    ngroups <- nlevels(x[,group, drop=TRUE])
    set.seed(seed)
    use.color <- sample(color.palette, ngroups)
  }
  p <- ggplot(data = x, mapping = aes(x = UMAP1, y = UMAP2, color = .data[[group]])) + 
    scale_color_manual(values = use.color) + 
    theme_void() + 
    theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0), 
          legend.key.size = unit(0.1, 'mm'), legend.justification = "left", 
          legend.text = element_text(size = 7)) +
    guides(color = guide_legend(title="", nrow = legend.nrow, bycol = TRUE, 
                                override.aes = list(size=2.5)))
  if (nrow(x)<3e4) { # rasterize dots if >=30K cells
    p <- p + geom_point(size = size, stroke = 0) 
  } else {
    p <- p + ggrastr::rasterise(geom_point(size=size,  stroke = 0), dpi = 300)
  }
  if (no_legend) {
    p <- p + theme(legend.position = "none")
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title[1]) + 
      theme(plot.title = element_text(hjust = 0.5))
    if (length(title)>1) {
      p <- p + labs(subtitle = title[2])
    }
  }
  if (draw.line) {
    p <- p + annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                      x = Inf, xend = Inf, linewidth = 0.75) 
  }
  return(p + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
}
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
plot_umap_dataset <- function(emb.list, dataset, size = 0.01) {
  #size <- ifelse(dataset %in% c("Lung", "Immune (human)"), 0.1, 0.01)
  plts <- list()
  plts[[1]] <- plot_grid(plotlist = lapply(c("batch", "cell"), function(x) {
    if (x=="batch") {
      plot_umap(x = emb.list[[1]], by = x, size = size) + 
        theme(legend.position = "none", 
              plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + 
        labs(subtitle = "full(-)")
    } else {
      plot_umap(x = emb.list[[1]], by = x, size = size) + 
        theme(legend.position = "none")  + 
        annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                 x = Inf, xend = Inf, linewidth = 0.5)
    }
  }), ncol=2, align = "vh")
  plts[[2]] <- plot_grid(plotlist = lapply(c("batch", "cell"), function(x) {
    if (x=="batch") {
      plot_umap(x = emb.list[[2]], by = x, size = size) + 
        theme(legend.position = "none", 
              plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + 
        labs(subtitle = "full(+)")
    } else {
      plot_umap(x = emb.list[[2]], by = x, size = size) + 
        theme(legend.position = "none") + 
        annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                 x = Inf, xend = Inf, linewidth = 0.5) 
    }
  }), ncol=2, align = "vh")
  plts[[3]] <- plot_grid(plotlist = lapply(c("batch", "cell"), function(x) {
    if (x=="batch") {
      plot_umap(x = emb.list[[3]], by = x, size = size) + 
        labs(subtitle = "hvg(-)") + 
        theme(legend.key.spacing.x = unit(0.5, "mm"), 
              legend.key.spacing.y = unit(0.5, "mm"), 
              legend.text = element_text(margin = margin(l = unit(0.5, "mm"))),
              legend.justification = "left", 
              legend.position = c(0, -0.4), 
              plot.margin = unit(c(0, 0, 3, 0), "cm"))
    } else {
      plot_umap(x = emb.list[[3]], by = x, size = size) + 
        theme(legend.position = "none") + 
        annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                 x = Inf, xend = Inf, linewidth = 0.5)
    }
  }), ncol=2, align = "vh")
  plts[[4]] <- plot_grid(plotlist = lapply(c("batch", "cell"), function(x) {
    if (x=="batch") {
      plot_umap(x = emb.list[[4]], by = x, size = size) + 
        theme(legend.position = "none") + 
        labs(subtitle = "hvg(+)")
    } else {
      plot_umap(x = emb.list[[4]], by = x, size = size) + 
        annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                 x = Inf, xend = Inf, linewidth = 0.5) + 
        theme(legend.key.spacing.x = unit(0.5, "mm"), 
              legend.key.spacing.y = unit(0.5, "mm"),
              legend.text = element_text(margin = margin(l = unit(0.5, "mm"))),
              legend.justification = "left", 
              legend.position = c(-1, -0.4), 
              plot.margin = unit(c(0, 0, 3, 0), "cm"))
    }
  }), ncol=2, align = "vh")
  plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]], ncol = 2, labels = c(dataset), 
            rel_heights = c(0.375, 0.625), label_x = 0, label_y = 0.9,
            hjust = 0, vjust = 0)
}
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
plot_scib_metrics_heatmap <- function(scaled_data, meta_data, metrics_list, 
                                      legend_title = "", 
                                      show_annot_legend = TRUE, ...) {
  # Parse data
  row.colors <- c(RColorBrewer::brewer.pal(9, "Reds")[c(3, 5, 7, 9)], 
                  RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)])
  names(row.colors) <- levels(as.factor(meta_data[,"Output-Input",drop=T]))
  col.names <- unlist(metrics_list)
  names(col.names) <- gsub(pattern = paste(paste0(names(metrics_list), "\\."), collapse = "|"),
                           replacement = "", x = names(col.names))
  col.split <- unlist(lapply(X = names(metrics_list), FUN = function(x) rep(x, length(metrics_list[[x]]))))
  names(col.split) <- names(col.names)
  symbols.pos <- unlist(lapply(split(1:nrow(scaled_data), 
                                     factor(unlist(lapply(X = row.names(scaled_data), FUN = function(x) {
                                       y <- strsplit(x, split = "_")[[1]]
                                       y[length(y)-1]})), 
                                       levels = c("Coralysis", "scVI", "Scanorama", "fastMNN", "Harmony", 
                                                  "Seurat v4 CCA", "Seurat v4 RPCA", "Unintegrated"))), 
                               function(x) floor(mean(x))
  ))
  symbols.pos <- symbols.pos[!is.na(symbols.pos)]
  method.symbols <- rep(NA, nrow(scaled_data))
  symbols <- c(24, 0:1, 3:7)
  names(symbols) <- c("Coralysis", "scVI", "Scanorama", "fastMNN", "Harmony", 
                      "Seurat v4 CCA", "Seurat v4 RPCA", "Unintegrated")
  method.symbols[symbols.pos] <- symbols[names(symbols.pos)]
  # Plot heatmap
  heat.plt <- Heatmap(matrix = scaled_data, name = legend_title, 
                      cluster_rows = F, cluster_columns = F, 
                      split = factor(unlist(lapply(X = row.names(scaled_data), FUN = function(x) {
                        y <- strsplit(x, split = "_")[[1]]
                        y[length(y)-1]})), levels = c("Coralysis", "scVI", "Scanorama", "fastMNN", "Harmony", 
                                                      "Seurat v4 CCA", "Seurat v4 RPCA", "Unintegrated")),
                      row_labels = unlist(lapply(X = row.names(scaled_data), FUN = function(x) {
                        y <- strsplit(x, split = "_")[[1]]
                        y[length(y)]})),
                      row_title_gp = gpar(fontsize = 10),
                      left_annotation = HeatmapAnnotation(Methods = anno_simple(
                        unlist(lapply(X = row.names(scaled_data), FUN = function(x) {
                          y <- strsplit(x, split = "_")[[1]]
                          y[length(y)-1]})), 
                        pch = method.symbols, 
                        col = c("Coralysis" = "white", "scVI" = "white", "Scanorama" = "white", "fastMNN" = "white", "Harmony" = "white", 
                                "Seurat v4 CCA" = "white", "Seurat v4 RPCA" = "white", "Unintegrated" = "white")
                      ), 
                      which = "row", show_annotation_name = FALSE),
                      row_title = NULL,
                      top_annotation = HeatmapAnnotation(Metrics = col.split[colnames(scaled_data)], 
                                                         col = list(Metrics = c("Batch correction" = "#8a94c5ff", "Bio conservation" = "#f899b3ff")), 
                                                         show_annotation_name = F, show_legend = show_annot_legend), 
                      # HeatmapAnnotation(Metrics = anno_block(labels = c("Batch correction", "Bio conservation"), 
                      #                                                       labels_gp = gpar(col = "white"),
                      #                                                       gp = gpar(fill = c("#8a94c5ff", "#f899b3ff")))),
                      right_annotation = HeatmapAnnotation(Input_Output = unlist(lapply(X = row.names(scaled_data), FUN = function(x) {
                        y <- strsplit(x, split = "_")[[1]]
                        y[length(y)]})), which = "row", col = list(Input_Output = row.colors), show_annotation_name = F, 
                        show_legend = show_annot_legend), 
                      show_row_names = F, 
                      row_gap = unit(2.5, "mm"),
                      column_split = col.split[colnames(scaled_data)],
                      column_labels = col.names[colnames(scaled_data)], 
                      column_names_gp = gpar(fontsize = 5.5), column_names_rot = 60, 
                      column_gap = unit(2.5, "mm"),
                      col = circlize::colorRamp2(c(-4, -2, 0, 2, 4), c("#451077FF", "#721F81FF", "white", "#6DCD59FF", "#35B779FF")), # circlize::colorRamp2(c(-4, -2, 0, 2, 4), c("blue4", "blue1", "white", "red1", "red4")) 
                      ...)
  heat.plt <- grid.grabExpr(draw(heat.plt, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right"))
  return(heat.plt)
}
#
#------------------------------------------------------------------------------#