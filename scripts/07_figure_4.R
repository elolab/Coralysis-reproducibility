#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make Figure 4. 
# Date: 23/01/2025
# Last update: 23/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("ComplexHeatmap") # v.2.14.0
library("SingleCellExperiment") # v.1.20.1
#require("scran") # v. 1.26.2

# Folders to save results 
analysis <- "07_figure_4"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Process data 

# Import data
input <- file.path("data", analysis, "sce.rds")
sce <- readRDS(file = input)

# Params
batch.label <- "batch"
cell.label <- "cell_type"
nhvg <- 2000

# Split object into 'ref' & 'query'
sce[[batch.label]] <- factor(sce[[batch.label]])
batches <- levels(sce[[batch.label]])
ref.cells <- (sce[[batch.label]] %in% batches[c(1:4, 6:7)])
query.cells <- (!ref.cells)
ref <- sce[,ref.cells]
query <- sce[,query.cells]
rm(sce)
gc()
ref[[batch.label]] <- factor(x = ref[[batch.label]], 
                             levels = batches[c(1:4, 6:7)])
query <- list(
  "smartseq2" =  query[,query$batch==batches[8]],
  "indrop2" =  query[,query$batch==batches[5]]
)
query$smartseq2[[batch.label]] <- factor(x = query$smartseq2[[batch.label]], 
                                         levels = batches[8])
query$indrop2[[batch.label]] <- factor(x = query$indrop2[[batch.label]], 
                                       levels = batches[5])

# Colors to use 
color.palette <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                   "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                   "#8C564B", "#C49C94", "#E377C2") 
cell.cols <- color.palette[1:length(sort(unique(ref$cell_type)))]
names(cell.cols) <- sort(unique(ref$cell_type))
batch.cols <- scales::brewer_pal(palette = "Dark2")(8)
names(batch.cols) <- c(levels(ref$batch), names(query))
use.colors <- list("batch" = batch.cols, "cell_type" = cell.cols)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 4 A

## Dimensional reduction: before reference-mapping

## Reference
set.seed(1024)
m.hvg <- scran::modelGeneVar(x = ref)
top.hvg <- scran::getTopHVGs(stats = m.hvg, n = nhvg)
set.seed(123)
ref <- RunPCA(object = ref, assay.name = "logcounts", dimred.name = "unintPCA", 
              features = top.hvg)
set.seed(123)
ref <- RunUMAP(object = ref, dimred.type = "unintPCA", dimred.name = "unintUMAP", 
               umap.method = "umap")

# Plot 
vars2plt <- c(batch.label, cell.label)
names(vars2plt) <- vars2plt
unint.ref.umap.plts <- lapply(X = vars2plt, FUN = function(x) {
  PlotDimRed(object = ref, color.by = x, dimred = "unintUMAP", use.color = use.colors[[x]], 
             point.size = 0.5, point.stroke = 0) + 
    theme(legend.position = "none",  
          plot.margin =  margin(0.25, 0.1, 0.1, 0.1, "cm"), #unit(c(0, 0, 0, 0), "cm"), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    ggtitle(ifelse(x == batch.label, "reference - unintegrated", "")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "UMAP2", "")
    )
  
})
unint.ref.umap.plts <- cowplot::plot_grid(plotlist = unint.ref.umap.plts, 
                                          labels = c("A", ""), ncol = 2, align = "vh")
pdf(file = file.path(res.dir[1], "figure_4_A.pdf"), width = 5, height = 2.55)
print(unint.ref.umap.plts)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Select HVGs & train the reference (+ PCA/UMAP)

# HVG selection - reference object
ref[[batch.label]] <- factor(ref[[batch.label]])
m.hvg <- scran::modelGeneVar(ref, block=ref[[batch.label]])
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg.batch <- row.names(ref)[hvg.ordered[1:nhvg]]
rowData(ref) <- cbind("highly_variable" = row.names(ref) %in% top.hvg.batch, m.hvg)
ref <- ref[top.hvg.batch,]

# Multi-level integration - train reference
set.seed(1024)
ref <- RunParallelDivisiveICP(object = ref, batch.label = batch.label, threads = 4) # 17.28977 mins

# PCA
set.seed(123)
ref <- RunPCA(object = ref, assay.name = "joint.probability", dimred.name = "intPCA", 
              pca.method = "stats", return.model = TRUE)

# UMAP
set.seed(123)
ref <- RunUMAP(object = ref, dimred.type = "intPCA", dimred.name = "intUMAP", 
               umap.method = "umap", return.model = TRUE)

# Plot
int.ref.umap.plts <- lapply(X = vars2plt, FUN = function(x) {
  PlotDimRed(object = ref, color.by = x, dimred = "intUMAP", use.color = use.colors[[x]], 
             point.size = 0.5, point.stroke = 0) + 
    theme(legend.position = "bottom",  
          legend.key.size = unit(0.2, 'mm'), 
          legend.text = element_text(size = 10),
          legend.justification = "left", 
          plot.margin =  margin(0.25, 0.1, 0.1, 0.1, "cm"),
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3), nrow = 7, title = NULL)) + 
    ggtitle(ifelse(x == batch.label, "reference - integrated", "")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "UMAP2", "")
    )
})
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 4 C-D

## Reference mapping 
query.map <- lapply(X = query, FUN = function(x) {
  set.seed(1024)
  ReferenceMapping(ref = ref, query = x, ref.label = cell.label, 
                   project.umap = TRUE, dimred.name.prefix = "ref")
})

# Combine single SCE objects
ref.query <- cbind(query.map$smartseq2, query.map$indrop2) 
colData(ref)[,c("coral_labels", "coral_probability")] <- NA
reducedDimNames(ref)[3:4] <- gsub(pattern = "int", replacement = "ref", x = reducedDimNames(ref)[3:4])
unint.ref.dimreds <- reducedDims(ref)[1:2]
reducedDims(ref)[1:2] <- NULL
ref.query <- cbind(ref, ref.query[top.hvg.batch,])
ref.query$type <- ifelse(ref.query$batch %in% levels(ref$batch), "reference", 
                         ifelse(ref.query$batch == "smartseq2", "query - smartseq2", "query - indrop2"))
reducedDims(ref)[3:4] <- unint.ref.dimreds
reducedDimNames(ref)[3:4] <- names(unint.ref.dimreds)
reducedDims(ref) <- reducedDims(ref)[c(3, 4, 1, 2)]
reducedDimNames(ref)[3:4] <- gsub(pattern = "ref", replacement = "int", x = reducedDimNames(ref)[3:4])
rm(unint.ref.dimreds)

# Plot: fig. 3 C
use.colors[["type"]] <- c("reference" = "gray85", "query - smartseq2" = "#A6761D", "query - indrop2" = "#666666")
map.query.ref.umap.plts <- lapply(X = c("type", "cell_type"), FUN = function(x) {
  PlotDimRed(object = ref.query, color.by = x, dimred = "refUMAP", use.color = use.colors[[x]], 
             point.size = 0.5, point.stroke = 0) + 
    theme(legend.position = "none",  
          plot.margin =  margin(0.25, 0.1, 0.1, 0.1, "cm"), #unit(c(0, 0, 0, 0), "cm"), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    ggtitle(ifelse(x == cell.label, "ref + query - ground-truth", "")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "UMAP2", "")
    )
})
map.query.ref.umap.plts <- cowplot::plot_grid(plotlist = map.query.ref.umap.plts, 
                                              labels = c("C", ""), 
                                              ncol = 2, align = "vh")
pdf(file = file.path(res.dir[1], "figure_4_C.pdf"), width = 5, height = 2.55)
print(map.query.ref.umap.plts)
dev.off()

# Plot: fig. 3 D
use.colors[["coral_labels"]] <- use.colors$cell_type
map.query.umap.plts <- lapply(X = c("type", "coral_labels"), FUN = function(x) {
  PlotDimRed(object = ref.query[,ref.query$type!="reference"], color.by = x, dimred = "refUMAP", use.color = use.colors[[x]], 
             point.size = 0.5, point.stroke = 0) + 
    theme(legend.position = "none",  
          plot.margin =  margin(0.25, 0.1, 0.1, 0.1, "cm"), #unit(c(0, 0, 0, 0), "cm"), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    ggtitle(ifelse(x == "coral_labels", "query - predictions", "")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "UMAP2", "")
    )
})
map.query.umap.plts[[1]] <- map.query.umap.plts[[1]] + 
  theme(legend.position = "bottom",  
        legend.key.size = unit(0.2, 'mm'), 
        legend.text = element_text(size = 10),
        legend.justification = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3, title = NULL)) 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 4 E,G

# Quantify accuracy 
query.coldata <- as.data.frame(colData(ref.query)) %>% 
  filter(batch %in% batches[c(5, 8)]) %>% 
  mutate("Correct_classification" = (cell_type == as.character(coral_labels))) %>% 
  mutate("Correct_classification" = factor(ifelse(Correct_classification, "Correct", "Incorrect")))

# Accuracy and precision 
query.coldata$coral_labels <- as.character(query.coldata$coral_labels)
acc.prec.data <- split(x = query.coldata[,c("cell_type", "coral_labels")], 
                       f = as.character(query.coldata$batch))
conf.mtx <- lapply(X = acc.prec.data, FUN = function(x) {
  table(x$coral_labels, x$cell_type)
})
stopifnot(all(row.names(conf.mtx$indrop2) == colnames(conf.mtx$indrop2)))
stopifnot(all(row.names(conf.mtx$smartseq2) == colnames(conf.mtx$smartseq2)))
acc <- lapply(X = conf.mtx, FUN = function(x) (sum(diag(x)) / sum(x)) )

# Plot: fig. 4 E
map.indrop2.query.umap.plts <- lapply(X = c("cell_type", "coral_labels", "coral_probability"), FUN = function(x) {
  if (x %in% c("cell_type", "coral_labels")) {
    PlotDimRed(object = ref.query[,ref.query$type=="query - indrop2"], color.by = x, dimred = "refUMAP", use.color = use.colors[[x]], 
               point.size = 0.5, point.stroke = 0) + 
      theme(legend.position = "none",  
            plot.margin =  margin(0.75, 0.1, 0.1, 0.1, "cm"), 
            plot.title = ggtext::element_markdown(size = 7.5, hjust = 0.5)) + 
      ggtitle(ifelse(x == cell.label, "indrop2 - ground-truth", "indrop2 - predictions")) + 
      labs(
        x = ifelse(x == cell.label, "UMAP1", ""),
        y = ifelse(x == cell.label, "UMAP2", "")
      )
  } else {
    PlotExpression(object = ref.query[,ref.query$type=="query - indrop2"], color.by = x, dimred = "refUMAP", 
                   color.scale = "viridis", point.size = 0.5, point.stroke = 0) + 
      scale_colour_continuous(limits = c(0, 1), type = "viridis") + 
      theme(legend.position = "none",  
            plot.margin =  margin(0.75, 0.1, 0.1, 0.1, "cm"), 
            plot.title = ggtext::element_markdown(size = 7.5, hjust = 0.5)) + 
      ggtitle("indrop2 - confidence") + 
      labs(x = "", y = "")
  }
})
map.indrop2.query.umap.plts <- cowplot::plot_grid(
  NULL, cowplot::plot_grid(
    plotlist = map.indrop2.query.umap.plts, 
    labels = paste("Coralysis prediction accuracy for indrop2:", round(acc$indrop2, 3)), 
    label_size = 12, hjust = -0.15, label_fontface = "plain", align = "vh", ncol = 3), 
  labels = c("E", ""), rel_widths = c(0.025, 0.975), ncol = 2 
) 
pdf(file = file.path(res.dir[1], "figure_4_E.pdf"), width = 7.5, height = 2.75)
print(map.indrop2.query.umap.plts)
dev.off()

# Plot: fig. 4 G
map.smartseq2.query.umap.plts <- lapply(X = c("cell_type", "coral_labels", "coral_probability"), FUN = function(x) {
  if (x %in% c("cell_type", "coral_labels")) {
    PlotDimRed(object = ref.query[,ref.query$type=="query - smartseq2"], color.by = x, dimred = "refUMAP", use.color = use.colors[[x]], 
               point.size = 0.5, point.stroke = 0) + 
      theme(legend.position = "none",  
            plot.margin =  margin(0.75, 0.1, 0.1, 0.1, "cm"), 
            plot.title = ggtext::element_markdown(size = 7.5, hjust = 0.5)) + 
      ggtitle(ifelse(x == cell.label, "smartseq2 - ground-truth", "smartseq2 - predictions")) + 
      labs(
        x = ifelse(x == cell.label, "UMAP1", ""),
        y = ifelse(x == cell.label, "UMAP2", "")
      )
  } else {
    PlotExpression(object = ref.query[,ref.query$type=="query - smartseq2"], color.by = x, dimred = "refUMAP", 
                   color.scale = "viridis", point.size = 0.5, point.stroke = 0) + 
      scale_colour_continuous(limits = c(0, 1), type = "viridis") + 
      theme(legend.position = "bottom",  
            plot.margin =  margin(0.75, 0.1, 0.1, 0.1, "cm"), 
            plot.title = ggtext::element_markdown(size = 7.5, hjust = 0.5), 
            legend.title = element_blank(), 
            legend.text = element_text(size=8),  
            legend.margin = margin(0, 0, 0, 0, "cm"), 
            legend.direction = "horizontal", legend.key.size = unit(0.55, "cm"), 
            legend.key.height = unit(0.25, "cm")) + 
      ggtitle("smartseq2 - confidence") + 
      labs(x = "", y = "")
  }
})
map.smartseq2.query.umap.plts <- cowplot::plot_grid(
  NULL, cowplot::plot_grid(
    plotlist = map.smartseq2.query.umap.plts, 
    labels = paste("Coralysis prediction accuracy for smartseq2:", round(acc$smartseq2, 3)), 
    label_size = 12, hjust = -0.15, label_fontface = "plain", align = "vh", ncol = 3), 
  labels = c("G", ""), rel_widths = c(0.025, 0.975), ncol = 2 
) 
pdf(file = file.path(res.dir[1], "figure_4_G.pdf"), width = 7.5, height = 3.2)
print(map.smartseq2.query.umap.plts)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 4 F,H

# Plot: fig. 4 F
# Annotations
col.text <- matrix("black", ncol = 13, nrow = 13)
diag(col.text) <- "white"
col.text <- as.vector(col.text)
indrop2_bottom_ha = HeatmapAnnotation("Cells" = anno_barplot(log10(colSums(conf.mtx$indrop2))), 
                                      height = unit(0.5, "cm"))
indrop2_row_ha = rowAnnotation(
  "Prediction" = anno_text(
    paste(paste0(round((rowSums(conf.mtx$indrop2) - diag(conf.mtx$indrop2))/(rowSums(conf.mtx$indrop2))*100, 1), "%"), 
          paste0("(", (rowSums(conf.mtx$indrop2) - diag(conf.mtx$indrop2)), ")")),
    location = 0.5, just = "center", gp = gpar(fontsize = 9)
  )
)
# >colSums(conf.mtx$indrop2)
# acinar activated_stellate              alpha               beta              delta             ductal        endothelial 
# 115                 79                659                373                127                199                 23 
# epsilon              gamma         macrophage               mast quiescent_stellate            schwann 
# 2                 89                 17                 12                 23                  6 
# Heatmap
conf.mtx.indrop2.heat <- grid.grabExpr(draw(Heatmap(matrix = ((conf.mtx$indrop2 / rowSums(conf.mtx$indrop2))*100), 
                                                    name = "Cell %",
                                                    cluster_rows = F, cluster_columns = F,
                                                    show_column_names = F, row_names_side = "left",
                                                    col = circlize::colorRamp2(c(0, 10, 90, 95, 100), c("#F5F5F5", "gray50", "lightblue", "pink", "firebrick")), 
                                                    #circlize::colorRamp2(c(0, 50, 100), c("gray90", "red4", "red1")),
                                                    column_names_side = "top", column_names_rot = 45,
                                                    bottom_annotation = indrop2_bottom_ha,
                                                    top_annotation = HeatmapAnnotation(df=data.frame("GroundTruth" = colnames(conf.mtx$indrop2)),
                                                                                       col = list("GroundTruth" = use.colors$cell_type[colnames(conf.mtx$indrop2)]),
                                                                                       show_legend = FALSE, 
                                                                                       show_annotation_name = FALSE, gap = unit(0.1, "cm"), 
                                                                                       annotation_legend_param = list(title = "Ground-Truth")),
                                                    layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                                                      v = pindex(((conf.mtx$indrop2 / rowSums(conf.mtx$indrop2))*100), i, j)
                                                      grid.text(sprintf(as.character(as.numeric(round(v, 1)))), x, y, 
                                                                gp = gpar(fontsize = 7, col = col.text, fontface = "plain"))
                                                      if(slice_r != slice_c) {
                                                        grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
                                                      }
                                                    },
                                                    right_annotation = indrop2_row_ha, 
                                                    rect_gp = gpar(col = "white", lwd = 4), column_title = "Ground-truth", row_title = "Predicted", 
                                                    show_heatmap_legend = FALSE, row_names_gp = gpar(fontsize = 9), 
                                                    column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10)
)))

# Plot: fig. 4 H
# Annotations
smartseq2_bottom_ha = HeatmapAnnotation("Cells" = anno_barplot(log10(colSums(conf.mtx$smartseq2))), 
                                        height = unit(0.5, "cm"))
smartseq2_row_ha = rowAnnotation(
  "Prediction" = anno_text(
    paste(paste0(round((rowSums(conf.mtx$smartseq2) - diag(conf.mtx$smartseq2))/(rowSums(conf.mtx$smartseq2))*100, 1), "%"), 
          paste0("(", (rowSums(conf.mtx$smartseq2) - diag(conf.mtx$smartseq2)), ")")),
    location = 0.5, just = "center", gp = gpar(fontsize = 9)
  )
)
# > colSums(conf.mtx$smartseq2)
# acinar activated_stellate              alpha               beta              delta             ductal        endothelial 
# 188                 55               1008                308                127                444                 21 
# epsilon              gamma         macrophage               mast quiescent_stellate            schwann 
# 8                213                  7                  7                  6                  2 
conf.mtx.smartseq2.heat <- grid.grabExpr(draw(Heatmap(matrix = ((conf.mtx$smartseq2 / rowSums(conf.mtx$smartseq2))*100), 
                                   name = "Cell %",
                                   cluster_rows = F, cluster_columns = F,
                                   show_column_names = F, row_names_side = "left",
                                   col = circlize::colorRamp2(c(0, 10, 90, 95, 100), c("#F5F5F5", "gray50", "lightblue", "pink", "firebrick")),
                                   column_names_side = "top", column_names_rot = 45,
                                   bottom_annotation = smartseq2_bottom_ha,
                                   top_annotation = HeatmapAnnotation(df=data.frame("GroundTruth" = colnames(conf.mtx$smartseq2)),
                                                                      col = list("GroundTruth" = use.colors$cell_type[colnames(conf.mtx$smartseq2)]),
                                                                      show_annotation_name = FALSE, gap = unit(0.1, "cm"),
                                                                      show_legend = FALSE, 
                                                                      annotation_legend_param = list(title = "Ground-Truth")),
                                   layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                                     v = pindex(((conf.mtx$smartseq2 / rowSums(conf.mtx$smartseq2))*100), i, j)
                                     grid.text(sprintf(as.character(as.numeric(round(v, 1)))), x, y, gp = gpar(fontsize = 7, col = col.text))
                                     if(slice_r != slice_c) {
                                       grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
                                     }
                                   },
                                   right_annotation = smartseq2_row_ha, 
                                   rect_gp = gpar(col = "white", lwd = 4), column_title = "Ground-truth", row_title = "Predicted", 
                                   show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 9), 
                                   column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
                                   heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", 
                                                               title_gp = gpar(fontsize = 9), grid_height = unit(0.25, "cm"))
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))

# Export plots
conf.mtx.indrop2.heat <- cowplot::plot_grid(plotlist = list(conf.mtx.indrop2.heat), labels = "F")
pdf(file = file.path(res.dir[1], "figure_4_F.pdf"), width = 5.5, height = 2.75)
print(conf.mtx.indrop2.heat)
dev.off()
conf.mtx.smartseq2.heat <- cowplot::plot_grid(plotlist = list(conf.mtx.smartseq2.heat), labels = "H")
pdf(file = file.path(res.dir[1], "figure_4_H.pdf"), width = 5.5, height = 3.25)
print(conf.mtx.smartseq2.heat)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 4 A-H

# Plot altogether
plot_A_C <- cowplot::plot_grid(unint.ref.umap.plts, map.query.ref.umap.plts, ncol = 2)
plot_B_D <- cowplot::plot_grid(plotlist = c(int.ref.umap.plts, map.query.umap.plts), ncol = 4, align = "vh", 
                               labels = c("B", "", "D", ""))
plot_E_H <- cowplot::plot_grid(plotlist = list(map.indrop2.query.umap.plts, conf.mtx.indrop2.heat,
                                               map.smartseq2.query.umap.plts, conf.mtx.smartseq2.heat), 
                               nrow = 2, align = "vh", rel_heights = c(0.45, 0.55), rel_widths = c(0.55, 0.45))
pdf(file = file.path(res.dir[1], "figure_4_A_H.pdf"), width = 11, height = 12.5) #width = 9.25
cowplot::plot_grid(plot_A_C, plot_B_D, plot_E_H, ncol = 1, rel_heights = c(0.225, 0.38, 0.45))
dev.off()

# Plot individual plots missing
# B
int.ref.umap.plts <- cowplot::plot_grid(plotlist = int.ref.umap.plts,
                                        labels = c("B", ""), ncol = 2, align = "vh")
pdf(file = file.path(res.dir[1], "figure_4_B.pdf"), width = 5, height = 4.2)
print(int.ref.umap.plts)
dev.off()
# D
map.query.umap.plts <- cowplot::plot_grid(plotlist = map.query.umap.plts,
                                          labels = c("D", ""),
                                          ncol = 2, align = "vh")
pdf(file = file.path(res.dir[1], "figure_4_D.pdf"), width = 5, height = 3.2)
print(map.query.umap.plts)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object
saveRDS(object = ref, file = file.path(res.dir[3], "ref.rds"))
saveRDS(object = query, file = file.path(res.dir[3], "query.rds"))
saveRDS(object = query.map, file = file.path(res.dir[3], "query_map.rds"))
saveRDS(object = ref.query, file = file.path(res.dir[3], "ref_query.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Clean environment variables
rm(list=ls())
gc()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Session
sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
# [5] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
# [9] MatrixGenerics_1.10.0       matrixStats_1.1.0           ComplexHeatmap_2.14.0       Coralysis_1.0.0            
# [13] ggplot2_3.5.1               dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.7.1          colorspace_2.1-0          rjson_0.2.21              class_7.3-20              modeltools_0.2-23        
# [6] circlize_0.4.15           scuttle_1.8.4             markdown_1.5              bluster_1.8.0             XVector_0.38.0           
# [11] GlobalOptions_0.1.2       BiocNeighbors_1.16.0      gridtext_0.1.5            ggtext_0.1.2              clue_0.3-64              
# [16] rstudioapi_0.14           farver_2.1.1              RSpectra_0.16-1           fansi_1.0.4               xml2_1.3.3               
# [21] codetools_0.2-18          sparseMatrixStats_1.10.0  doParallel_1.0.17         jsonlite_1.8.4            Cairo_1.6-0              
# [26] umap_0.2.10.0             cluster_2.1.3             png_0.1-8                 pheatmap_1.0.12           compiler_4.2.1           
# [31] dqrng_0.3.0               Matrix_1.6-5              limma_3.54.2              aricode_1.0.2             cli_3.6.1                
# [36] BiocSingular_1.14.0       tools_4.2.1               rsvd_1.0.5                igraph_1.4.1              gtable_0.3.3             
# [41] glue_1.6.2                GenomeInfoDbData_1.2.9    RANN_2.6.1                reshape2_1.4.4            LiblineaR_2.10-22        
# [46] doRNG_1.8.6               Rcpp_1.0.10               vctrs_0.6.1               iterators_1.0.14          DelayedMatrixStats_1.20.0
# [51] xfun_0.37                 stringr_1.5.0             beachmat_2.14.0           lifecycle_1.0.3           irlba_2.3.5.1            
# [56] rngtools_1.5.2            statmod_1.5.0             edgeR_3.40.2              zlibbioc_1.44.0           scales_1.3.0             
# [61] doSNOW_1.0.20             parallel_4.2.1            SparseM_1.81              RColorBrewer_1.1-3        reticulate_1.34.0        
# [66] ggrastr_1.0.1             stringi_1.7.12            foreach_1.5.2             ScaledMatrix_1.6.0        scran_1.26.2             
# [71] BiocParallel_1.32.6       shape_1.4.6               commonmark_1.9.0          rlang_1.1.0               pkgconfig_2.0.3          
# [76] bitops_1.0-7              lattice_0.20-45           labeling_0.4.2            cowplot_1.1.1             tidyselect_1.2.0         
# [81] plyr_1.8.8                magrittr_2.0.3            R6_2.5.1                  snow_0.4-4                generics_0.1.3           
# [86] metapod_1.6.0             DelayedArray_0.24.0       pillar_1.9.0              withr_2.5.0               RCurl_1.98-1.10          
# [91] tibble_3.2.1              crayon_1.5.2              utf8_1.2.3                GetoptLong_1.0.5          locfit_1.5-9.7           
# [96] flexclust_1.4-1           digest_0.6.31             openssl_2.0.6             munsell_0.5.0             viridisLite_0.4.1        
# [101] beeswarm_0.4.0            vipor_0.4.5               askpass_1.1        
#
#------------------------------------------------------------------------------#