#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make figure 6 and supplementary figure 14.
# Date: 25/01/2025
# Last update: 25/01/2025
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
#require("cowplot") # v.1.1.1
#require("zellkonverter") # v.1.8.0

# Folders to save results 
analysis <- "12_figure_6"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)

# Download the CD34+ HSCs data from Persad et al., 2023, Nat Biotechnol (PMID:36973557):
# https://zenodo.org/records/6383269
# options(timeout = 300)
# file.down <- "https://zenodo.org/records/6383269/files/cd34_multiome_rna.h5ad?download=1"
# download.file(url = file.down, destfile = file.path("data", analysis, "cd34_multiome_rna.h5ad"))
# Import the data
input <- file.path("data", analysis, "cd34_multiome_rna.h5ad")
sce <- zellkonverter::readH5AD(file = input) 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 A

## Multi-level integration

# Params
assayNames(sce)[1] <- "logcounts"
batch.label <- "sample"
cell.label <- "celltype" # "cell_type"

# HVG
set.seed(1204)
top.hvg <- scran::getTopHVGs(sce, n = 2000)
sce <- sce[top.hvg,]

# Multi-level integration
set.seed(1024)
sce <- RunParallelDivisiveICP(object = sce, threads = 4)

# PCA
set.seed(123)
sce <- RunPCA(object = sce)

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "PCA", umap.method = "uwot",
               n_neighbors = 100, min_dist = 0.5)

#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 A

# Plot 
celltype.persad.umap <- PlotDimRed(object = sce, color.by = cell.label, dimred = "UMAP", 
                                   point.size = 0.75, point.stroke = 0, legend.nrow = 2, 
                                   seed.color = 1983) + 
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.1, 'mm'), 
        legend.justification = "center", 
        plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
        legend.text = element_text(size = 10), 
        plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
  guides(color = guide_legend(title="", nrow = 2, bycol = TRUE, 
                              override.aes = list(size=4))) + 
  ggtitle("Original annotations<br>(as in Persad *et al*., 2023)") 
pdf(file = file.path(file.path(res.dir[1], "figure_6_A.pdf")), width = 3.25, height = 4.25)
print(cowplot::plot_grid(celltype.persad.umap, labels = "A"))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 B

## Cell cluster probability score - mean
sce <- SummariseCellClusterProbability(object = sce, icp.round = 4)

# Plot probability scores
prob.persad.umap <- PlotExpression(object = sce, color.by = "mean_probs", dimred = "UMAP", 
                                   point.size = 0.75, point.stroke = 0, color.scale = "viridis") + 
  scale_color_viridis_c(name = "Mean Cell Cluster Probability", breaks = c(0.4, 0.6, 0.8)) +
  labs(x = "", y = "") + 
  theme(legend.position = "bottom", 
        legend.justification = "center", 
        plot.margin =  unit(c(0, 0.2, 0, 0), "cm"), 
        legend.text = element_text(size = 10), 
        legend.title.align=0.5, 
        legend.title.position = "top", 
        plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
  guides(color = guide_colourbar(theme = theme(
    legend.key.width  = unit(10, "lines"), 
    legend.key.height  = unit(1, "lines")
  ))) + 
  ggtitle("Coralysis cell cluster probability<br>(mean across *L*=50 runs)") 
pdf(file = file.path(file.path(res.dir[1], "figure_6_B.pdf")), width = 3.25, height = 4.25)
print(cowplot::plot_grid(prob.persad.umap, labels = "B"))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 C

## Correlation between Palantir pseudotime vs mean cell cluster probability
color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
ngroups <- nlevels(sce[[cell.label]])
set.seed(1983)
use.color <- sample(color.palette, ngroups)
pseudo_vs_probs <- cor(x = sce$palantir_pseudotime, y = sce$mean_probs, method = "pearson")
corr.prob.pseudo.plt <- ggplot(data = as.data.frame(colData(sce)),
                               mapping = aes(x = scaled_mean_probs, y = palantir_pseudotime, color = celltype))  +
  geom_point(size = 0.75, stroke = 0) +
  theme_classic() + 
  scale_color_manual(values = use.color) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.1, 'mm'), 
        legend.justification = "center", 
        plot.margin =  unit(c(0, 0.2, 0, 0), "cm"), 
        legend.text = element_text(size = 10), 
        plot.title = element_text(size = 9.5, hjust = 0.5), 
        axis.title = ggtext::element_markdown()) + 
  guides(color = guide_legend(title="", nrow = 2, bycol = TRUE, 
                              override.aes = list(size=4))) + 
  labs(x = "Mean cell cluster probability<br>(scaled 0-1)", y = "Palantir pseudotime<br>(as in Persad *et al*., 2023)", 
       title = bquote(rho ~ "=" ~ .(round(pseudo_vs_probs, 2))))  + 
  xlim(0, 1)
pdf(file = file.path(file.path(res.dir[1], "figure_6_C.pdf")), width = 4, height = 4.25)
print(cowplot::plot_grid(corr.prob.pseudo.plt, labels = "C"))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 A-C

## Plot altogether
plts.cd34.all <- cowplot::plot_grid(celltype.persad.umap, prob.persad.umap, corr.prob.pseudo.plt, 
                                    ncol = 3, align = "vh", labels = LETTERS[1:3])
pdf(file = file.path(file.path(res.dir[1], "figure_6_A_C.pdf")), width = 12, height = 4.75)
print(plts.cd34.all)
dev.off()

# Export SCE object
saveRDS(object = sce, file = file.path(res.dir[3], "sce_persad.rds"))

# Clean environment variables
rm(list=ls())
gc()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import the next dataset

# Variables
analysis <- "12_figure_6"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))

## Embryonic development
# download embryoid body from Weiler et al., 2024, Nat Methods (PMID:38871986):
# https://figshare.com/ndownloader/files/41674599
# options(timeout = 300)
# file.down <- "https://figshare.com/ndownloader/files/41674599" 
# download.file(url = file.down, destfile = file.path("data", analysis, "embryoid_body.h5ad"))

# Import the data
input <- file.path("data", analysis, "embryoid_body.h5ad")
sce <- zellkonverter::readH5AD(file = input) 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 D-E

# Dimensional reduction
assayNames(sce)[1] <- "logcounts"
batch.label <- "stage"
cell.label <- "cell_type" 
nhvg <- 6000
sce[[batch.label]] <- factor(sce[[batch.label]])
m.hvg <- scran::modelGeneVar(sce, block=sce[[batch.label]])
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg <- row.names(sce)[hvg.ordered[1:nhvg]]
rowData(sce) <- cbind("highly_variable" = row.names(sce) %in% top.hvg, m.hvg)
sce <- sce[top.hvg,]

# Multi-level integration
set.seed(1024)
sce <- RunParallelDivisiveICP(object = sce, batch.label = batch.label, threads = 4) # 45.56014 mins

# PCA
set.seed(123)
sce <- RunPCA(object = sce)

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "PCA", umap.method = "uwot",
               n_neighbors = 100, min_dist = 0.5)

# Plot UMAP with cell type annotations
use.colors <- list("stage" = RColorBrewer::brewer.pal(5, "Reds"))
names(use.colors$stage) <- levels(sce$stage)
stage.celltype.weiler.umap <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, use.color = use.colors[[x]], dimred = "UMAP", 
             point.size = 0.75, point.stroke = 0, seed.color = 1024) + 
    theme(legend.position = "bottom",  
          legend.key.size = unit(0.1, 'mm'), 
          legend.justification = "center", 
          plot.margin =  unit(c(0, 0.2, 0, 0), "cm"), 
          legend.text = element_text(size = 10), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    guides(color = guide_legend(title="", nrow = 4, bycol = TRUE, 
                                override.aes = list(size=4))) + 
    ggtitle(ifelse(x == cell.label, "Original annotations<br>(as in Weiler *et al*., 2024)", "Embryo stage")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "UMAP2", "")
    )
})
stage.celltype.weiler.umap.all <- cowplot::plot_grid(plotlist = stage.celltype.weiler.umap, ncol = 2, align = "h", labels = c("D", "E"))
pdf(file = file.path(file.path(res.dir[1], "figure_6_D_E.pdf")), width = 8.2, height = 5)
print(stage.celltype.weiler.umap.all)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 F

## Cell cluster probability score - mean
sce <- SummariseCellClusterProbability(object = sce, icp.round = 4)

# Plot probability scores
prob.weiler.umap <- PlotExpression(object = sce, color.by = "mean_probs", dimred = "UMAP", 
                                   point.size = 0.75, point.stroke = 0, color.scale = "viridis") + 
  scale_color_viridis_c(name = "Mean Cell Cluster Probability", breaks = c(0.3, 0.6, 0.9)) +
  labs(x = "", y = "") + 
  theme(legend.position = "bottom", 
        legend.justification = "center", 
        plot.margin =  unit(c(0, 0.2, 0, 0), "cm"), 
        legend.text = element_text(size = 10), 
        legend.title.align=0.5, 
        legend.title.position = "top", 
        plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
  guides(color = guide_colourbar(theme = theme(
    legend.key.width  = unit(10, "lines"), 
    legend.key.height  = unit(1, "lines")
  ))) + 
  ggtitle("Coralysis cell cluster probability<br>(mean across *L*=50 runs)") 
pdf(file = file.path(file.path(res.dir[1], "figure_6_F.pdf")), width = 3.25, height = 4.25)
print(cowplot::plot_grid(prob.weiler.umap, labels = "F"))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 G

## Cell Probability distribution per cell type per stage
color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
ngroups <- nlevels(sce[[cell.label]])
set.seed(1024)
use.color <- sample(color.palette, ngroups)
cell.probs.dist.celltype.plt <- ggplot(data = as.data.frame(colData(sce)), 
                                       mapping = aes(x = scaled_mean_probs, color = stage)) + 
  stat_density(geom="line",position="identity", linewidth = 0.75) +
  facet_grid(~cell_type) +
  theme_classic() + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  scale_color_manual(name = "Embryo stage", values = use.colors$stage) + 
  labs(x = "Mean cell cluster probability\n(scaled 0-1)", y = "Density") + 
  theme(strip.background = element_rect(color = NA), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 7))
# solution from: https://github.com/tidyverse/ggplot2/issues/2096
g <- ggplot_gtable(ggplot_build(cell.probs.dist.celltype.plt))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- use.color
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 D-G

# Plot altogether
panels_D_F <- cowplot::plot_grid(plotlist = c(stage.celltype.weiler.umap, list(prob.weiler.umap)), 
                                 ncol = 3, align = "h", labels = LETTERS[4:6])
pdf(file = file.path(file.path(res.dir[1], "figure_6_D_F.pdf")), width = 12, height = 5)
print(panels_D_F)
dev.off()
pdf(file = file.path(file.path(res.dir[1], "figure_6_G.pdf")), width = 12, height = 2)
grid::grid.draw(g)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 6 H-Q & Supplementary Figure 14

# Cell states
meta.cells <- BinCellClusterProbability(object = sce, label = cell.label, 
                                        icp.round = 4, funs = "mean", bins = 20)
write.table(x = as.data.frame(colData(meta.cells)),
            file = file.path(res.dir[2], "bins_cell_cluster_probabilities_by_cell_type_embryoid.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Differential expression programs per cell state 
corr.gexp.probs <- CellBinsFeatureCorrelation(object = meta.cells)
write.table(x = cbind("Gene" = row.names(corr.gexp.probs), corr.gexp.probs),
            file = file.path(res.dir[2], "correlated_genes_by_prob_bin_by_cell_type_embryoid.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot top 5 up- & downregulated
heat.plts <- list()
for (cell in levels(sce$cell_type)) {
  tmp <- corr.gexp.probs[,cell, drop = FALSE] %>% filter(!is.na(.data[[cell]]))
  tmp <- tmp %>% arrange(desc(.data[[cell]]))
  genes <- row.names(tmp)
  up.top5.corr <- tmp[1:5,] %>% `names<-`(genes[1:5])
  up.top5.corr <- rev(up.top5.corr)
  down.top5.corr <- tmp[(nrow(tmp)-4):nrow(tmp),] %>% `names<-`(genes[(nrow(tmp)-4):nrow(tmp)])
  down.top5.corr <- rev(down.top5.corr)
  pick.bins <- grepl(cell, colnames(meta.cells))
  data.plt <- t(scale(t(logcounts(meta.cells)[c(names(up.top5.corr), names(down.top5.corr)), pick.bins])))
  row.annot <- rowAnnotation("Pearson" = anno_text(round(c(up.top5.corr, down.top5.corr), 3), gp = gpar(fontsize = 7)))
  col_fun2 <- circlize::colorRamp2(c(0.4, 0.6, 0.9), scales::viridis_pal()(3))
  col.data <-  colData(meta.cells[,pick.bins])[order(meta.cells[,pick.bins]$aggregated_probability_bins),]
  top.annot <- HeatmapAnnotation("Probability" = col.data$aggregated_probability_bins, simple_anno_size = unit(0.25, "cm"), 
                                 col = list("Probability" = col_fun2), show_annotation_name = FALSE, 
                                 show_legend = FALSE)
  col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  heat.plts[[cell]] <- Heatmap(matrix = data.plt[,row.names(col.data)], name = "Row Z-score", cluster_rows = F, cluster_columns = F,
                               row_names_side = "left", row_split = rep(c("up", "down"), each =5), 
                               right_annotation = row.annot, top_annotation = top.annot, col = col_fun,
                               column_title = cell, row_names_gp = gpar(fontsize = 7.5), 
                               column_names_gp = gpar(fontsize = 2), show_heatmap_legend = FALSE, row_title = NULL)
  scatter.plts <- lapply(X = c("mean_probs", names(up.top5.corr), names(down.top5.corr)), FUN = function(x) {
    PlotExpression(sce[,sce$cell_type==cell], color.by = x, 
                   dimred = "UMAP", scale.values = T, 
                   color.scale = ifelse(x == "mean_probs", "viridis", "inferno"), 
                   legend.title = ifelse(x == "mean_probs", "Prob.", x), 
                   point.size = 0.75, point.stroke = 0)
  })
  pdf(file = file.path(res.dir[1], paste0(gsub("\\/", "_", cell), "_scatter_plot_top5_up_down_Pearson_correlated_gene_exp.pdf")), 
      width = 18, height = 4)
  print(cowplot::plot_grid(plotlist = c(scatter.plts[1:6], list(NULL), scatter.plts[7:11]), ncol = 6, align="hv", axis="tblr"))
  dev.off()
}
top.annot.example <- HeatmapAnnotation("Probability" = col.data$aggregated_probability_bins, simple_anno_size = unit(0.25, "cm"), 
                                       col = list("Probability" = col_fun2), show_annotation_name = TRUE, 
                                       show_legend = TRUE, annotation_legend_param = list(Probability = list(direction = "horizontal", title_position = "topcenter")))
heat.example <- Heatmap(matrix = data.plt[,row.names(col.data)], name = "Row Z-score", cluster_rows = F, cluster_columns = F,
                        row_names_side = "left", row_split = rep(c("up", "down"), each =5), 
                        right_annotation = row.annot, top_annotation = top.annot.example, col = col_fun,
                        column_title = cell, row_names_gp = gpar(fontsize = 7), 
                        column_names_gp = gpar(fontsize = 4), show_heatmap_legend = TRUE, row_title = NULL, 
                        heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topcenter"))
pdf(file = file.path(res.dir[1], "example_legend_heatmap_top5_up_down_Pearson_correlated_genes_with_binned_probability.pdf"), 
    width = 8, height = 4)
draw(heat.example)
dev.off()
heat.plts.parsed <- lapply(heat.plts, function(x) draw(x) %>% grid.grabExpr()) 
pdf(file = file.path(res.dir[1], "supp_figure_14.pdf"), 
    width = 12, height = 4)
patchwork::wrap_plots(heat.plts.parsed, ncol = length(heat.plts.parsed)/2)
dev.off()

# Cell bin distribution by stage 
cellbins.tables <- TabulateCellBinsByGroup(object = meta.cells, group = "stage")

# Mesoderm
cell <- "Mesoderm"
cell.bins <- grepl(cell, colnames(meta.cells))
col.data <-  colData(meta.cells[,cell.bins])[order(meta.cells[,cell.bins]$aggregated_probability_bins),]
col.fun.viridis <- circlize::colorRamp2(c(-2, 0, 2), viridis::inferno(3))
colnames(cellbins.tables$Mesoderm) <- paste0("bin_", colnames(cellbins.tables$Mesoderm))
top.annot.meso <- HeatmapAnnotation("Probability" = col.data$aggregated_probability_bins, simple_anno_size = unit(0.25, "cm"), 
                                    col = list("Probability" = col_fun2), show_annotation_name = FALSE, 
                                    show_legend = TRUE, annotation_legend_param = list(Probability = list(direction = "horizontal", title_position = "topcenter")))
stage.cells.meso.heat <- Heatmap(matrix = scale(cellbins.tables$Mesoderm), name = "Column Z-score", cluster_rows = F, cluster_columns = F, 
                                 col = col.fun.viridis, top_annotation = top.annot.meso, row_names_side = "left", 
                                 row_names_gp = gpar(fontsize = 8.5), column_names_gp = gpar(fontsize = 5.5), 
                                 heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topcenter"))
pdf(file = file.path(res.dir[1], "figure_6_H.pdf"), width = 3, height = 2.5)
draw(stage.cells.meso.heat, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

### ESC
cell <- "EN-1"
cell.bins <- grepl(cell, colnames(meta.cells))
col.data <-  colData(meta.cells[,cell.bins])[order(meta.cells[,cell.bins]$aggregated_probability_bins),]
col.fun.viridis <- circlize::colorRamp2(c(-2, 0, 2), viridis::inferno(3))
colnames(cellbins.tables$`EN-1`) <- paste0("bin_", colnames(cellbins.tables$`EN-1`))
top.annot.en1 <- HeatmapAnnotation("Probability" = col.data$aggregated_probability_bins, simple_anno_size = unit(0.25, "cm"), 
                                   col = list("Probability" = col_fun2), show_annotation_name = FALSE, 
                                   show_legend = TRUE, annotation_legend_param = list(Probability = list(direction = "horizontal", title_position = "topcenter")))
stage.cells.en1.heat <- Heatmap(matrix = scale(cellbins.tables$`EN-1`), name = "Column Z-score", cluster_rows = F, cluster_columns = F, 
                                col = col.fun.viridis, top_annotation = top.annot.en1, row_names_side = "left", 
                                row_names_gp = gpar(fontsize = 8.5), column_names_gp = gpar(fontsize = 5.5), 
                                heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topcenter"))

pdf(file = file.path(res.dir[1], "figure_6_M.pdf"), width = 3, height = 2.5)
draw(stage.cells.en1.heat, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
#
#------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" 
saveRDS(object = sce, file = file.path(res.dir[3], "sce_weiler.rds"))
saveRDS(object = meta.cells, file = file.path(res.dir[3], "metacells.rds"))

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
# [6] IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_1.1.0          
# [11] ComplexHeatmap_2.14.0       Coralysis_1.0.0             ggplot2_3.5.1               dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.7.1          colorspace_2.1-0          rjson_0.2.21              class_7.3-20              modeltools_0.2-23         circlize_0.4.15          
# [7] scuttle_1.8.4             markdown_1.5              bluster_1.8.0             XVector_0.38.0            GlobalOptions_0.1.2       BiocNeighbors_1.16.0     
# [13] gridtext_0.1.5            ggtext_0.1.2              clue_0.3-64               rstudioapi_0.14           farver_2.1.1              RSpectra_0.16-1          
# [19] fansi_1.0.4               xml2_1.3.3                codetools_0.2-18          sparseMatrixStats_1.10.0  doParallel_1.0.17         jsonlite_1.8.4           
# [25] Cairo_1.6-0               cluster_2.1.3             png_0.1-8                 uwot_0.1.14               pheatmap_1.0.12           compiler_4.2.1           
# [31] dqrng_0.3.0               basilisk_1.10.2           Matrix_1.6-5              limma_3.54.2              aricode_1.0.2             cli_3.6.1                
# [37] BiocSingular_1.14.0       tools_4.2.1               rsvd_1.0.5                igraph_1.4.1              gtable_0.3.3              glue_1.6.2               
# [43] GenomeInfoDbData_1.2.9    RANN_2.6.1                reshape2_1.4.4            LiblineaR_2.10-22         doRNG_1.8.6               Rcpp_1.0.10              
# [49] vctrs_0.6.1               zellkonverter_1.8.0       iterators_1.0.14          DelayedMatrixStats_1.20.0 xfun_0.37                 stringr_1.5.0            
# [55] beachmat_2.14.0           lifecycle_1.0.3           irlba_2.3.5.1             rngtools_1.5.2            statmod_1.5.0             edgeR_3.40.2             
# [61] zlibbioc_1.44.0           scales_1.3.0              basilisk.utils_1.10.0     doSNOW_1.0.20             parallel_4.2.1            SparseM_1.81             
# [67] RColorBrewer_1.1-3        reticulate_1.34.0         gridExtra_2.3             ggrastr_1.0.1             stringi_1.7.12            foreach_1.5.2            
# [73] ScaledMatrix_1.6.0        scran_1.26.2              filelock_1.0.2            BiocParallel_1.32.6       shape_1.4.6               commonmark_1.9.0         
# [79] rlang_1.1.0               pkgconfig_2.0.3           bitops_1.0-7              lattice_0.20-45           labeling_0.4.2            patchwork_1.2.0          
# [85] cowplot_1.1.1             tidyselect_1.2.0          RcppAnnoy_0.0.20          plyr_1.8.8                magrittr_2.0.3            R6_2.5.1                 
# [91] snow_0.4-4                generics_0.1.3            metapod_1.6.0             DelayedArray_0.24.0       pillar_1.9.0              withr_2.5.0              
# [97] RCurl_1.98-1.10           tibble_3.2.1              dir.expiry_1.6.0          crayon_1.5.2              utf8_1.2.3                GetoptLong_1.0.5         
# [103] viridis_0.6.2             locfit_1.5-9.7            flexclust_1.4-1           digest_0.6.31             munsell_0.5.0             beeswarm_0.4.0           
# [109] viridisLite_0.4.1         vipor_0.4.5          
#
#------------------------------------------------------------------------------#