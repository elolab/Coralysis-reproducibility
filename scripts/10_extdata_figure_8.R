#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make extended data figure 8. 
# Date: 24/01/2025
# Last update: 24/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("ggalluvial") # v.0.12.5
library("ComplexHeatmap") # v.2.14.0
library("SingleCellExperiment") # v.1.20.1
source("scripts/helper_functions.R")
#require("scran") # v. 1.26.2
#require("cowplot") # v.1.1.1

# Folders to save results 
analysis <- "10_extdata_figure_8"
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
emb.dims <- 30
color.palette <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                   "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                   "#8C564B", "#C49C94", "#E377C2") 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 8 A

# HVG
set.seed(1024)
m.hvg <- scran::modelGeneVar(sce)
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg <- row.names(sce)[hvg.ordered[1:nhvg]]

# PCA
set.seed(123)
sce <- scran::fixedPCA(sce, subset.row=top.hvg, rank=emb.dims)

# Select colors for cell type
use.cols <- color.palette[1:length(unique(sce$cell_type))]
names(use.cols) <- levels(sce$cell_type)

# Calculate the distance of cell type centroids between query-reference
centroids <- reducedDim(sce, "PCA") %>% 
  as.data.frame(.) %>% cbind(., colData(sce)[,c(batch.label, cell.label)]) %>% 
  mutate("cell_id" = row.names(.)) %>%
  group_by(batch, cell_type) %>% 
  summarise_at(vars(PC1:PC30), mean) %>% 
  as.data.frame(.) %>%
  `rownames<-`(paste(.$batch, .$cell_type, sep="_"))
distance <- dist(centroids)
dist.conf.mtx <- as.matrix(distance)[1:13,14:26]
all(gsub("CTRL_", "", row.names(dist.conf.mtx)) == gsub("STIM_", "", colnames(dist.conf.mtx)))
diag.dist <- diag(dist.conf.mtx)
names(diag.dist) <- gsub("STIM_", "", colnames(dist.conf.mtx))
diag.dist <- sort(diag.dist)
dist.conf.mtx.trans <- t(apply(dist.conf.mtx, 1, function(x) (1 - (x / max(x))) ))
row.names(dist.conf.mtx.trans) <- colnames(dist.conf.mtx.trans) <- gsub("STIM_", "", colnames(dist.conf.mtx.trans)) 

# PCA plot
pca.data2plot <- reducedDim(sce, type="PCA")[,1:2] %>% 
  as.data.frame(.) %>% cbind(., colData(sce))
unint.pca.plot <- ggplot(data = pca.data2plot, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(fill = cell_type, shape = batch), size = 1.5, color = "white") + 
  scale_shape_manual(values = c(21, 24), breaks = c("CTRL", "STIM"), labels = c("reference", "query")) +
  scale_fill_manual(values = use.cols[levels(sce$cell_type)]) + 
  guides(fill = guide_legend(override.aes = list(size = 3.5, shape=21), nrow = 5, title = NULL),
         shape = guide_legend(override.aes = list(size = 3.5, alpha = 0.7, 
                                                  color = scales::hue_pal()(2), 
                                                  fill = scales::hue_pal()(2)), 
                              nrow = 4, title = NULL)) +
  theme_minimal() + 
  theme(plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"), 
        legend.position = "bottom", 
        legend.margin = margin(0, 0, 0, 0, "cm"), 
        legend.spacing.x = unit(5, "pt")) + 
  theme(legend.key.spacing.y = unit(0.1, "pt"), legend.key.spacing.x = unit(0.1, "pt"))
unint.pca.plot.xdens <- cowplot::axis_canvas(unint.pca.plot, axis = "x") +
  geom_density(data = pca.data2plot, aes(x = PC1, fill = cell_type, color = cell_type),
               alpha = 0.7, size = 0.2) +
  scale_fill_manual(values = use.cols[levels(sce$cell_type)]) +
  scale_color_manual(values = use.cols[levels(sce$cell_type)])
unint.pca.plot.ydens <-cowplot:: axis_canvas(unint.pca.plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = pca.data2plot, aes(x = PC2, fill = batch, color = batch),
               alpha = 0.7, size = 0.2) +
  geom_vline(data = centroids, mapping = aes(xintercept = PC2),
             color = use.cols[centroids$cell_type], 
             linetype = "dashed", size = 0.75) +
  coord_flip()
unint.pca.plot2 <- cowplot::insert_xaxis_grob(unint.pca.plot, unint.pca.plot.xdens, grid::unit(.2, "null"), position = "top")
unint.pca.plot2 <- cowplot::insert_yaxis_grob(unint.pca.plot2, unint.pca.plot.ydens, grid::unit(.2, "null"), position = "right")
panel_A <- cowplot::plot_grid(cowplot::ggdraw(unint.pca.plot2), labels = "A")
pdf(file = file.path(res.dir[1], "extdata_figure_8_A.pdf"), width = 4.5, height = 5.5)
print(panel_A)
dev.off()
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#
## Extended Data Figure 8 B

# Heatmap
col.heat.annot <- HeatmapAnnotation(df = data.frame(Cell_labels = names(diag.dist)), 
                                    col = list(Cell_labels = use.cols[names(diag.dist)]), 
                                    which = "col", show_annotation_name = FALSE, show_legend = FALSE)
row.annots <- data.frame(Dissimilarity = ((diag.dist -min(diag.dist)) / (max(diag.dist) - min(diag.dist))), 
                         Cell_labels = names(diag.dist))
row.heat.annot <- HeatmapAnnotation(df = row.annots, 
                                    which = "row", 
                                    col = list(Dissimilarity = circlize::colorRamp2(c(0, 0.25, 0.75, 1), 
                                                                                    RColorBrewer::brewer.pal(9,"Blues")[seq(2,9,2)]), 
                                               Cell_labels = use.cols[names(diag.dist)]), 
                                    show_legend = c(TRUE, FALSE), show_annotation_name = c(FALSE, FALSE), 
                                    annotation_legend_param = list(Dissimilarity = list(direction = "horizontal", title_position = "topcenter")))
col_fun <- circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8), 
                                RColorBrewer::brewer.pal(9,"Reds")[seq(1,9,2)])
heat.cell_type.centroids.dist <- Heatmap(matrix = dist.conf.mtx.trans[names(diag.dist), names(diag.dist)], 
                                         name = "Closeness", 
                                         cluster_columns = FALSE, cluster_rows = FALSE, 
                                         show_row_dend = FALSE, show_column_dend = FALSE, 
                                         rect_gp = gpar(col = "white", lwd = 2), 
                                         left_annotation = row.heat.annot, 
                                         top_annotation = col.heat.annot, 
                                         show_column_names = FALSE, show_row_names = FALSE, 
                                         col = col_fun, row_title = "reference", column_title = "query", 
                                         heatmap_legend_param = list(direction = "horizontal", 
                                                                     title_position = "topcenter", 
                                                                     legend_width = unit(2, "cm")))
pdf(file = file.path(res.dir[1], "extdata_figure_8_B.pdf"), width = 3.5, height = 3.5)
draw(heat.cell_type.centroids.dist, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
heat.cell_type.centroids.dist2 <- grid.grabExpr(draw(
  Heatmap(matrix = dist.conf.mtx.trans[names(diag.dist), names(diag.dist)], 
          name = "Closeness", 
          cluster_columns = FALSE, cluster_rows = FALSE, 
          show_row_dend = FALSE, show_column_dend = FALSE, 
          rect_gp = gpar(col = "white", lwd = 2), 
          left_annotation = row.heat.annot, 
          top_annotation = col.heat.annot, 
          show_column_names = FALSE, show_row_names = FALSE, 
          col = col_fun, row_title = "reference", column_title = "query", 
          column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), 
          heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", legend_width = unit(2, "cm"))), 
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom"
))
panel_A_B <- cowplot::plot_grid("", unint.pca.plot2, 
                                cowplot::plot_grid(heat.cell_type.centroids.dist2, "", nrow=2, rel_heights = c(0.7, 0.3)), 
                                "",
                                rel_widths = c(0.07, 0.35, 0.3, 0.05), ncol = 4, 
                                labels = c("", "A", "B", ""))  
pdf(file = file.path(res.dir[1], "extdata_figure_8_A_B.pdf"), width = 7, height = 4.5)
print(panel_A_B)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Reference mapping: STIM vs. CTRL (i.e., query vs reference)

# Split objects into ref-query
map <- list(
  "ref" = sce[,sce$batch=="CTRL"],
  "query" = sce[,sce$batch=="STIM"]
)

# HVG selection in ref 
set.seed(1024)
m.hvg <- scran::modelGeneVar(map$ref)
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg <- row.names(map$ref)[hvg.ordered[1:nhvg]]
map$ref <- map$ref[top.hvg,]

# Train the reference
set.seed(123)
map$ref <- RunParallelDivisiveICP(object = map$ref, threads = 4)
map$ref <- RunPCA(object = map$ref, return.model = TRUE)
map$ref <- RunUMAP(map$ref, return.model=TRUE)

# Map the query against the reference
map$map <- ReferenceMapping(ref = map$ref, query = map$query, ref.label = cell.label, project.umap = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 8 D-F

# Plot UMAP
obj <- c("ref", "map", "ref", "map", "map", "map")
label <- c("batch", "batch", "cell_type", "cell_type", "coral_labels", "coral_probability")
title <- c("reference", "query", "reference - ground-truth", "query - ground-truth", "query - predictions", "query - confidence")
out.plots <- list()
for (i in seq_along(obj)) {
  o <- obj[i]
  l <- label[i]
  if (l == "batch") {
    cols <- scales::pal_hue()(2)
    cols <- ifelse(o == "ref", cols[1], cols[2]) 
  } else {
    cols <- use.cols
  }
  if (l != "coral_probability") {
    p <- PlotDimRed(object = map[[o]], color.by = l, use.color = cols, 
                    point.size = 0.75, point.stroke = 0, 
                    plot.theme = theme_minimal()) + 
      theme(legend.position = "none") 
  } else {
    p <- PlotExpression(object = map[[o]], color.by = l, color.scale = "viridis", 
                        point.size = 0.5, point.stroke = 0, 
                        plot.theme = theme_minimal()) + 
      scale_colour_continuous(limits = c(0, 1), type = "viridis") + 
      theme(legend.title = element_blank(), 
            legend.text = element_text(size=8),  
            legend.margin = margin(0, 0, 0, 0, "cm"), 
            legend.direction = "horizontal", legend.key.size = unit(0.55, "cm"), 
            legend.key.height = unit(0.25, "cm"), 
            legend.position = c(0.75, 0.1)) 
  }
  if (!((o == "map") && (l == "batch"))) {
    p <- p + theme(axis.title = element_blank())
  }
  out.plots[[i]] <- p + ggtitle(title[i]) +
    theme(plot.title = element_text(size=9.5, hjust=0.5))
}
panel_D_F <- cowplot::plot_grid(plotlist = out.plots, byrow = FALSE, ncol = 3, 
                                labels = c("D", "E", "F", "", "", ""))
pdf(file = file.path(res.dir[1], "extdata_figure_8_D_F.pdf"), width = 8.5, height = 5.5)
print(panel_D_F)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 8 C

# Classifications
query.coldata <- as.data.frame(colData(map$map)) %>%
  filter(batch %in% "STIM") %>%
  mutate("Correct_classification" = (cell_type == as.character(coral_labels))) %>%
  mutate("Correct_classification" = factor(ifelse(Correct_classification, "Correct", "Incorrect")))

# Accuracy and precision
query.coldata$coral_labels <- as.character(query.coldata$coral_labels)
query.coldata$cell_type <- as.character(query.coldata$cell_type)
acc.prec.data <- split(x = query.coldata[,c("cell_type", "coral_labels")],
                       f = as.character(query.coldata$batch))
conf.mtx <- lapply(X = acc.prec.data, FUN = function(x) {
  table(x$coral_labels, x$cell_type)
})
stopifnot(all(row.names(conf.mtx$STIM) == colnames(conf.mtx$STIM)))
acc <- lapply(X = conf.mtx, FUN = function(x) (sum(diag(x)) / sum(x)) ) # 0.8856529

# See cell numbers 
cells.batch <- colData(sce)[,c(batch.label, cell.label)] %>% 
  as.data.frame(.) %>% 
  mutate("cell_id" = row.names(.)) %>% 
  group_by(batch, cell_type) %>% 
  tally(., name = "abs") %>% ungroup(.) %>% 
  group_by(batch) %>% 
  mutate("perc" = abs / sum(abs) * 100, 
         "cells" = 1:length(unique(cell_type))) %>% 
  as.data.frame(.)

# Barplot with correct/misclassifications
cells.batch.preds <- query.coldata[,c(batch.label, cell.label, "coral_labels")] %>% 
  mutate("cell_id" = row.names(.), 
         "coral_labels" = factor(coral_labels, levels = levels(sce$cell_type)), 
         "cell_type" = factor(cell_type, levels = levels(sce$cell_type))) %>% 
  group_by(batch, cell_type, coral_labels) %>% 
  tally(., name = "abs") %>% ungroup(.) %>% 
  group_by(batch) %>% 
  mutate("perc" = abs / sum(abs) * 100) %>% 
  as.data.frame(.) %>% 
  left_join(., y=cells.batch[cells.batch$batch=="STIM", c("cell_type", "cells"), drop=FALSE])
perc.cells.query.preds <- ggplot(cells.batch.preds, aes(axis1 = cell_type, axis2 = coral_labels, y = perc)) +
  geom_alluvium(aes(color = cell_type, fill = cell_type)) +
  geom_stratum(color="gray") +
  geom_text(stat = "stratum", aes(label = paste(after_stat(stratum))), size = 2) + 
  scale_fill_manual(values = use.cols[levels(cells.batch.preds$cell_type)]) + 
  scale_color_manual(values = use.cols[levels(cells.batch.preds$cell_type)]) + 
  theme_classic() + 
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.y = element_line(color="darkgray"),
        axis.ticks.y = element_line(color="darkgray"), 
        legend.key.size = unit(0.3, 'cm'), 
        legend.spacing.y = unit(0.2, 'cm'), title = element_text(size=9), 
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        legend.margin = margin(0, 0, 0, 0, "cm")) +  
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0, 0)) +
  guides(color = guide_legend(ncol = 1, title = NULL, byrow = TRUE), 
         fill = guide_legend(ncol = 1, title = NULL, byrow = TRUE)) +
  ylab("Percentage of cell labels:\nground-truth vs predictions") + 
  ggtitle(paste("Coralysis prediction accuracy:", round(acc$STIM, 2)))
panel_C <- cowplot::plot_grid(perc.cells.query.preds, labels = "C")
pdf(file = file.path(res.dir[1], "extdata_figure_8_C.pdf"), width = 5, height = 4.5)
print(panel_C)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 8

# Plot altogether
pdf(file = file.path(res.dir[1], "extdata_figure_8.pdf"), width = 9.5, height = 8.5)
cowplot::plot_grid(
  cowplot::plot_grid(panel_A_B, panel_C, ncol = 2, rel_widths  = c(0.58, 0.42)), 
  cowplot::plot_grid(NULL, panel_D_F, NULL, ncol = 3, rel_widths = c(0.06, 0.8, 0.1)), 
  nrow = 2, rel_heights = c(0.43, 0.57) 
)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" 
saveRDS(object = map, file = file.path(res.dir[3], "sce.rds"))

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
# [9] MatrixGenerics_1.10.0       matrixStats_1.1.0           ComplexHeatmap_2.14.0       ggalluvial_0.12.5          
# [13] Coralysis_1.0.0             ggplot2_3.5.1               dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.7.1          colorspace_2.1-0          rjson_0.2.21              class_7.3-20              modeltools_0.2-23        
# [6] circlize_0.4.15           scuttle_1.8.4             bluster_1.8.0             XVector_0.38.0            GlobalOptions_0.1.2      
# [11] BiocNeighbors_1.16.0      clue_0.3-64               rstudioapi_0.14           farver_2.1.1              RSpectra_0.16-1          
# [16] fansi_1.0.4               codetools_0.2-18          sparseMatrixStats_1.10.0  doParallel_1.0.17         jsonlite_1.8.4           
# [21] Cairo_1.6-0               umap_0.2.10.0             cluster_2.1.3             png_0.1-8                 pheatmap_1.0.12          
# [26] compiler_4.2.1            dqrng_0.3.0               Matrix_1.6-5              limma_3.54.2              aricode_1.0.2            
# [31] cli_3.6.1                 BiocSingular_1.14.0       tools_4.2.1               rsvd_1.0.5                igraph_1.4.1             
# [36] gtable_0.3.3              glue_1.6.2                GenomeInfoDbData_1.2.9    RANN_2.6.1                reshape2_1.4.4           
# [41] LiblineaR_2.10-22         doRNG_1.8.6               Rcpp_1.0.10               vctrs_0.6.1               iterators_1.0.14         
# [46] DelayedMatrixStats_1.20.0 stringr_1.5.0             beachmat_2.14.0           lifecycle_1.0.3           irlba_2.3.5.1            
# [51] rngtools_1.5.2            statmod_1.5.0             edgeR_3.40.2              zlibbioc_1.44.0           scales_1.3.0             
# [56] doSNOW_1.0.20             parallel_4.2.1            SparseM_1.81              RColorBrewer_1.1-3        reticulate_1.34.0        
# [61] ggrastr_1.0.1             stringi_1.7.12            foreach_1.5.2             ScaledMatrix_1.6.0        scran_1.26.2             
# [66] BiocParallel_1.32.6       shape_1.4.6               rlang_1.1.0               pkgconfig_2.0.3           bitops_1.0-7             
# [71] lattice_0.20-45           purrr_1.0.1               labeling_0.4.2            cowplot_1.1.1             tidyselect_1.2.0         
# [76] plyr_1.8.8                magrittr_2.0.3            R6_2.5.1                  snow_0.4-4                generics_0.1.3           
# [81] metapod_1.6.0             DelayedArray_0.24.0       pillar_1.9.0              withr_2.5.0               RCurl_1.98-1.10          
# [86] tibble_3.2.1              crayon_1.5.2              utf8_1.2.3                GetoptLong_1.0.5          locfit_1.5-9.7           
# [91] flexclust_1.4-1           digest_0.6.31             tidyr_1.3.0               openssl_2.0.6             munsell_0.5.0            
# [96] viridisLite_0.4.1         beeswarm_0.4.0            vipor_0.4.5               askpass_1.1               
#
#------------------------------------------------------------------------------#