#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make Figure 1 and Extended Data Figure 2. 
# Date: 20/01/2025
# Last update: 20/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages & Set import global variables

# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("ComplexHeatmap") # v.2.14.0
library("SingleCellExperiment") # v.1.20.1
source("scripts/helper_functions.R") 

# Folders to save results 
analysis <- "01_figure_1"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)

# Colors to use in figures
colors <- list(
  "Sample" = c("V1" = "#beaed4ff", "V2" = "#7fc97fff"),
  "cell_type" = c("aDC" = "#f1ce63ff", "B mem" = "#b6992dff", "B naive" = "#4e79a7ff", 
                  "CD4 mem" = "#f28e2bff", "CD4 naive" = "#ffbe7dff", "CD8 eff" = "#59a14fff", 
                  "CD8 T" = "#8cd17dff", "HSC" = "#79706eff", "Megakaryocyte" = "#bab0acff", 
                  "Monocyte" = "#b07aa1ff", "CD16+ monocyte" = "#499894ff", 
                  "NK" = "#86bcb6ff", "pDC" = "#d37295ff", "Treg" = "#e15759ff")
) 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import & Process data

## Import data
data.dir <- file.path("data", analysis)
# import "Source Data Fig.3" excel file with cell annotations from Harmony paper: https://doi.org/10.1038/s41592-019-0619-0
meta.data <- readxl::read_excel(path = file.path(data.dir, "41592_2019_619_MOESM9_ESM.xlsx")) 
datasets <- c("pbmc6k_3pV1/hg19", "pbmc8k_3pV2/GRCh38")
file.dirs <- file.path(data.dir, datasets)
names(file.dirs) <- names(datasets) <- c("V1", "V2")
sce.list <- lapply(X = setNames(names(datasets), names(datasets)), FUN = function(x) {
  DropletUtils::read10xCounts(samples = file.dirs[x], sample.names = x, col.names = TRUE)
})

# Rename genes (Ensembl 'ID' --> 'Symbol')
rename.duplicates <- function(x, x.append) {
  dup.uniq <- which(duplicated(x))
  dup.names.uniq <- x[dup.uniq] 
  dup.all <- which(x %in% dup.names.uniq) 
  dup.names.all <- x[dup.all]
  dup.renames <- paste(dup.names.all, x.append[dup.all], sep = "_")
  x[dup.all] <- dup.renames
  return(x)
}
convert.ids <- c("ID", "Symbol")
rowdata.cols <- colnames(rowData(sce.list[[1]]))
for (i in names(sce.list)) {
  if (all(convert.ids %in% rowdata.cols)) {
    if (all(row.names(sce.list[[i]])==rowData(sce.list[[i]])$ID)) {
      is.dup <- any(duplicated(rowData(sce.list[[i]])$Symbol))
      new.gene.ids <- rowData(sce.list[[i]])$Symbol
      if (is.dup) {
        old.gene.ids <- rowData(sce.list[[i]])$ID
        new.gene.ids <- rename.duplicates(x = new.gene.ids, x.append = old.gene.ids)
        rowData(sce.list[[i]])$SymbolDeduplicated <- new.gene.ids
      } 
      row.names(sce.list[[i]]) <- new.gene.ids
    }
  } 
}

# Merge SCE objects
counts <- lapply(sce.list, function(x) as.data.frame(as.matrix(counts(x))))
merged.counts <- merge(counts$V1, counts$V2, by = "row.names", all = TRUE)
row.names(merged.counts) <- merged.counts$Row.names
merged.counts <- merged.counts[,colnames(merged.counts)!="Row.names"]
merged.counts <- as(as.matrix(merged.counts), "sparseMatrix")
merged.meta.data <- rbind(colData(sce.list$V1), colData(sce.list$V2)) 
sce <- SingleCellExperiment(assays = list(counts = merged.counts), 
                            colData = merged.meta.data)

# Filter cells based on cell annotated in Harmony paper & add annotations to SCE object
meta.data <- filter(meta.data, dataset %in% c("three_prime_v1", "three_prime_v2"))
stopifnot(all(!duplicated(gsub("-1", "", merged.meta.data$Barcode)))) 
colnames(sce) <- gsub("-1", "", colnames(sce)) # remove "-1" barcode suffix to match 'cell_id' from Harmony excel data annotations
sce <- sce[,meta.data$cell_id] # 13,800 --> 13,189 (removed 611)
stopifnot(all(colnames(sce) == meta.data$cell_id))
colData(sce)[,"cell_subtype"] <- meta.data$cell_subtype
cell.label <- "cell_subtype"

# Remove non-annotated cells, i.e, 'NA' in 'cell_subtype'
pick.cells <- (!is.na(sce[[cell.label]])) # 13,189 --> 13,046 (143 cells removed) 
sce <- sce[,pick.cells]

# Filter non-expressed genes
sce <- sce[!rowAnyNAs(counts(sce)),] # remove genes represented by NAs 36,227 --> 30,205 (removed 6,022) 
sce <- sce[(rowSums(counts(sce))>0),] # remove non-expressed genes 30,205 --> 20,016 (removed 10,189)

# Rename and merge cell types
sce[["cell_type"]] <- as.factor(sce[[cell.label]])
cell.label <- "cell_type"
levels(sce[["cell_type"]]) <- c( # 'cd8mem' & 'cd8naive' cell types were merged into 'CD8 T'
  "adc" = "aDC", "bmem" = "B mem", "bnaive" = "B naive", "cd4mem" = "CD4 mem",  
  "cd4naive" = "CD4 naive", "cd8eff" = "CD8 eff", "cd8mem" = "CD8 T", "cd8naive" = "CD8 T",
  "hsc" = "HSC", "mk" = "Megakaryocyte", "mono14" = "Monocyte", "mono16" = "CD16+ monocyte",  
  "nk" = "NK", "pdc" = "pDC", "treg" = "Treg" 
) 

# Normalize data - function from the 'helper_functions.R' script
sce <- NormalizeData(sce)

# HVG selection
nhvg <- 2000
batch.label <- "Sample"
sce[[batch.label]] <- factor(sce[[batch.label]])
m.hvg <- scran::modelGeneVar(sce, block=sce[[batch.label]])
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg <- row.names(sce)[hvg.ordered[1:nhvg]]
rowData(sce) <- cbind("highly_variable" = row.names(sce) %in% top.hvg, m.hvg)
sce <- sce[top.hvg,]

# Export object
saveRDS(object = sce, file = file.path(res.dir[3], "sce_proc.rds")) # SCE object processed
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 1 A

## Dimensional reduction: before integration

# PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "logcounts", dimred.name = "unintPCA")

# Elbow plot
pdf(file = file.path(res.dir[1], "unintegrated_Elbow_plot.pdf"), width = 4, height = 4)
PCAElbowPlot(sce, dimred.name = "unintPCA")
dev.off()

# t-SNE
pcs <- 1:12
set.seed(123)
sce <- RunTSNE(object = sce, dims = pcs, dimred.type = "unintPCA", dimred.name = "unintTSNE")

# Plot t-SNE
unint.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "unintTSNE", use.color = colors[[x]], 
             point.size = 0.75, point.stroke = 0, plot.theme = theme_void()) + 
    guides(color = guide_legend(title = "", nrow = 6, bycol = TRUE, override.aes = list(size=4))) + 
    theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(0.1, 'mm'), legend.justification = "left",
          legend.text = element_text(size = 10), plot.margin =  unit(c(0, 0, 0, 0), "cm")) 
})
pdf(file = file.path(res.dir[1], "unintegrated_TSNE_plot.pdf"), width = 6.5, height = 4.5)
cowplot::plot_grid(plotlist = unint.plts, ncol = 2, align = "vh")
dev.off()
unint.plts.together <- cowplot::plot_grid((unint.plts$Sample + ggtitle("Unintegrated") +
                                             theme(plot.title = element_text(hjust=0.5))),
                                          unint.plts$cell_type, ncol = 1, align = "vh")

#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Multi-level integration
set.seed(1024)
sce <- RunParallelDivisiveICP(object = sce, batch.label = batch.label, 
                              C = 1, train.k.nn.prop = 0.45, 
                              ari.cutoff = 0.1, allow.free.k = FALSE,
                              threads = 4, verbose = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Dimensional reduction

# Run PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "joint.probability", dimred.name = "intPCA")

# Run t-SNE
set.seed(123)
sce <- RunTSNE(object = sce, dimred.type = "intPCA", dimred.name = "intTSNE")

# Export intermediate SCE object 
saveRDS(object = sce, file = file.path(res.dir[3], "sce_int.rds"))

# Plots
int.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "intTSNE", use.color = colors[[x]],
             point.size = 0.75, legend.nrow = 4, point.stroke = 0, plot.theme = theme_void()) + 
    guides(color = guide_legend(title = "", nrow = 6, bycol = TRUE, override.aes = list(size=4))) + 
    theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(0.1, 'mm'), legend.justification = "left",
          legend.text = element_text(size = 10), plot.margin =  unit(c(0, 0, 0, 0), "cm")) 
})
pdf(file = file.path(res.dir[1], "integrated_TSNE_plot.pdf"), width = 6.5, height = 4.5)
cowplot::plot_grid(plotlist = int.plts, ncol = 2, align = "vh")
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 1 B-D 

## Pick the probability table with 16 clusters with the highest SD

# Find ICP cluster table with the highest standard deviation
by <- log2(metadata(sce)$coralysis$k)
L <- metadata(sce)$coralysis$L
K <- 1:by
pick.k16 <- seq(by, L*by, by)
probs <- GetCellClusterProbability(object = sce, icp.round = 4, concatenate = FALSE)
max.sd <- which.max(unlist(lapply(probs, sd))) # L = 4 (ICP model)

# Plot cluster trees based on previous function
clt.tree.plts <- list()
clt.tree.plts[["median_cluster_probability"]] <- PlotClusterTree(object = sce, icp.run = max.sd)
clt.tree.plts[["batch"]] <- PlotClusterTree(object = sce, color.by = batch.label, icp.run = max.sd, 
                                            legend.title = "Batch", use.color = colors[[batch.label]])
clt.tree.plts[["celltype"]] <- PlotClusterTree(object = sce, color.by = cell.label, icp.run = max.sd, 
                                               legend.title = "Cell Type", use.color = colors[[cell.label]])
pdf(file = file.path(res.dir[1], paste0("median_cluster_probability_tree_L", max.sd, ".pdf")), width = 6, height = 4)
clt.tree.plts[["median_cluster_probability"]]
dev.off()
pdf(file = file.path(res.dir[1], "batch_cluster_tree_L4.pdf"), width = 8, height = 5)
clt.tree.plts[["batch"]]
dev.off()
pdf(file = file.path(res.dir[1], "celltype_cluster_tree_L4.pdf"), width = 8, height = 5)
clt.tree.plts[["celltype"]]
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 1E

# Project every ICP round
plts <- list()
for (k in K) { # ICP rounds
  # Params
  k.round <- as.character(2**k)
  plts[[k.round]] <- list()
  pick.k <- seq(k, L*by, by)
  # Dimreds
  set.seed(123)
  sce <- RunPCA(object = sce, assay.name = "joint.probability", select.icp.tables = pick.k, dimred.name = "pca")
  set.seed(123)
  sce <- RunTSNE(object = sce, dimred.type = "pca", dimred.name = "tsne")
  # Plots
  plts[[k.round]][["batch"]] <- PlotDimRed(object = sce, color.by = batch.label, 
                                           use.color = colors[[batch.label]], 
                                           dimred = "tsne", point.size = 0.75,
                                           point.stroke = 0, plot.theme = theme_void(),
                                           legend.nrow = 4) + 
    ggtitle(paste("Cluster target:", 2**k)) +
    theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=10), 
          plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
    annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", x = Inf, xend = Inf, linewidth = 0.5)
  plts[[k.round]][["celltype"]] <- PlotDimRed(object = sce, color.by = cell.label, 
                                              use.color = colors[[cell.label]], 
                                              dimred = "tsne", point.size = 0.75, 
                                              point.stroke = 0, plot.theme = theme_void(),
                                              legend.nrow = 4) + 
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
    annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", x = Inf, xend = Inf, linewidth = 0.5)
}
pdf(file = file.path(res.dir[1], "ICP_rounds_projected_onto_TSNE_plots.pdf"), width = 8.5, height = 4)
cowplot::plot_grid(plts$`2`$batch, plts$`4`$batch, plts$`8`$batch, plts$`16`$batch, 
                   plts$`2`$celltype, plts$`4`$celltype, plts$`8`$celltype, plts$`16`$celltype, 
                   ncol = 4, align = "vh")
dev.off()

# Plot clustering & probability for the ICP run w/ highest SD
sce <- SummariseCellClusterProbability(object = sce, icp.run = 4, icp.round = 4, scale.funs = FALSE, save.in.sce = TRUE)
sce$icp_run_round_4_4_clusters <- factor(sce$icp_run_round_4_4_clusters, levels = as.character(1:16))
# Cluster
plts[[k.round]][["cluster"]] <- PlotDimRed(object = sce, color.by = "icp_run_round_4_4_clusters", dimred = "tsne", 
                                           use.color = NULL, point.size = 0.5, point.stroke = 0, 
                                           legend.nrow = 4, seed = 1024, plot.theme = theme_void(), 
                                           label = TRUE) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "none") 
plts[[k.round]][["cluster"]] <- PlotDimRed(object = sce, color.by = "icp_run_round_4_4_clusters", dimred = "tsne", 
                                           use.color = NULL, point.size = 0.5, point.stroke = 0, 
                                           legend.nrow = 4, seed = 1024, plot.theme = theme_void(), 
                                           label = FALSE) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "none") +
  geom_label(data = plts[[k.round]][["cluster"]]$data, mapping = aes(label = label), size = 2.5)
# Probability
plts[[k.round]][["prob"]] <- PlotExpression(object = sce, color.by = "icp_run_round_4_4_probs", dimred = "tsne", 
                                            color.scale = "viridis", point.stroke = 0, point.size = 0.5, 
                                            plot.theme = theme_void()) + 
  scale_color_viridis_c(limits = c(0, 1)) + 
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title = element_text(hjust=0.5, size=10)) +
  ggtitle(paste0("Cluster target: ", k, " (L=", max.sd, ")"))

# Assemble figure
pdf(file = file.path(res.dir[1], paste0("ICP_rounds_projected_onto_TSNE_plots_highlighting_L", max.sd, ".pdf")), width = 7.5, height = 3.25)
cowplot::plot_grid(plts$`2`$batch, plts$`4`$batch, plts$`8`$batch, plts$`16`$batch, plts$`16`$prob,
                   plts$`2`$celltype, plts$`4`$celltype, plts$`8`$celltype, plts$`16`$celltype, plts$`16`$cluster,
                   ncol = 5, align = "vh")
dev.off()
reducedDims(sce)[reducedDimNames(sce) %in% c("pca", "tsne")] <- NULL

## Assemble figure 1 A-E
tree.plts <- cowplot::plot_grid(cowplot::plot_grid(clt.tree.plts[["median_cluster_probability"]], NULL, ncol=1, 
                                                   rel_heights = c(0.875,0.125)), 
                                cowplot::plot_grid(clt.tree.plts[["batch"]] + theme(legend.position = "none"), 
                                                   clt.tree.plts[["celltype"]] + theme(legend.position = "none"),
                                                   NULL, ncol=1, labels = c("C", "D", "")), 
                                ncol=2, rel_widths = c(0.4, 0.6), labels = c("B", ""))
icp.rounds.tsne.plots <- cowplot::plot_grid(plts$`2`$batch, plts$`4`$batch, plts$`8`$batch, plts$`16`$batch, plts$`16`$prob,
                                            plts$`2`$celltype, plts$`4`$celltype, plts$`8`$celltype, plts$`16`$celltype, plts$`16`$cluster,
                                            ncol = 5, align = "vh")
fig1AE <- cowplot::plot_grid(unint.plts.together, cowplot::plot_grid(tree.plts, icp.rounds.tsne.plots, 
                                                                     ncol=1, rel_heights = c(0.6, 0.4), 
                                                                     labels = c("", "E")), 
                             ncol=2, rel_widths = c(0.25, 0.75), labels = c("A", ""))
pdf(file = file.path(res.dir[1], "figure1_A_E.pdf"), width = 12, height = 8.75)
print(fig1AE)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 1 F

## Looking into the coefficients for L ICP run with the maximum SD
gene.coeff <- GetFeatureCoefficients(object = sce, icp.run = 4)

# Parse top 10 positive gene coefficients per ICP round per cluster
n.row <- 10 * sum(c(1, 4, 8, 16))
data.coeff <- data.frame("K" = vector(mode = "integer", length = n.row), 
                         "cluster" = vector(mode = "integer", length = n.row), 
                         "coeff" = vector(mode = "numeric", length = n.row), 
                         "gene" = vector(mode = "character", length = n.row))
j <- 1
for (k in seq_along(gene.coeff)) {
  for (i in seq_len(ncol(gene.coeff[[k]][,-1,drop=F]))) {
    tmp <- gene.coeff[[k]][,c(1, i+1)]
    top10 <- head(tmp[order(tmp[,2], decreasing = TRUE),], n = 10)
    h <- (j+9)
    data.coeff[j:h,"K"] <- k
    data.coeff[j:h,"cluster"] <- i
    data.coeff[j:h, "coeff"] <- top10[, 2]
    data.coeff[j:h,"gene"] <- top10$feature
    j <- (h + 1)
  }
}

# Exclude neutral weightsand
data.coeff <- data.coeff[data.coeff$coeff>0,] # 11 excluded - remained 279
write.table(x = data.coeff, 
            file = file.path(res.dir[2], paste0("top10_positive_coefficients_by_cluster_divisive_ICP_L", max.sd, "_table.tsv")), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Gene names that appear multiple times
table(duplicated(data.coeff$gene)) # 110 genes (169 unique)

# Average gene expression by the 16 clusters found
cluster.colors <- c("#F781BF", "#66A61E", "#FC8D62", "#CCCCCC", "#666666", 
                    "#BF5B17", "#666666", "#FDDAEC", "#7570B3", "#A6CEE3", 
                    "#FF7F00", "#FFED6F", "#E5C494", "#A6761D", "#FDCDAC", 
                    "#FDB462")
names(cluster.colors) <- as.character(1:16)
Kannot <- as.factor(data.coeff$K)
levels(Kannot) <- paste0("K", levels(Kannot))
avg.gexp <- Coralysis:::AggregateClusterExpression(mtx = logcounts(sce[data.coeff$gene]), 
                                                   cluster = sce$icp_run_round_4_4_clusters)
colnames(avg.gexp) <- gsub("cluster", "", colnames(avg.gexp))
pick.genes <- data.coeff %>% 
  mutate("idx" = 1:nrow(.)) %>% 
  group_by(., K, cluster) %>% 
  slice_max(., n=3, order_by = coeff)
heat.plt <- Heatmap(matrix = scale(t(avg.gexp)), name = "Row Z-score", 
                    cluster_rows = F, cluster_columns = F, show_row_names = T,  
                    show_column_names = F,
                    col = circlize::colorRamp2(c(-4, -2, 0, 2, 4), c("blue4", "blue1", "white", "red1", "red4")),
                    column_split = Kannot, 
                    left_annotation = HeatmapAnnotation(df = data.frame("Cluster" = as.character(1:16)), 
                                                        col = list("Cluster" = cluster.colors), which = "row", 
                                                        show_annotation_name = TRUE, annotation_name_side = "top", 
                                                        annotation_name_rot = 75,
                                                        show_legend = FALSE), 
                    heatmap_legend_param = list(legend_direction = "horizontal"), 
                    row_names_side = "left", row_names_gp = gpar(fontsize = 10)) 
heat.plt2 <- Heatmap(matrix = scale(t(avg.gexp)), name = "Row Z-score", 
                     cluster_rows = F, cluster_columns = F, show_row_names = T,  
                     show_column_names = F,
                     col = circlize::colorRamp2(c(-4, -2, 0, 2, 4), c("blue4", "blue1", "white", "red1", "red4")),
                     column_split = Kannot, 
                     left_annotation = HeatmapAnnotation(df = data.frame("Cluster" = as.character(1:16)), 
                                                         col = list("Cluster" = cluster.colors), which = "row", 
                                                         show_annotation_name = TRUE, annotation_name_side = "top", 
                                                         annotation_name_rot = 75,
                                                         show_legend = FALSE), 
                     heatmap_legend_param = list(legend_direction = "horizontal"), 
                     row_names_side = "left", row_names_gp = gpar(fontsize = 10), 
                     column_names_side = "bottom", 
                     bottom_annotation = HeatmapAnnotation("Genes" = anno_mark(at = pick.genes$idx, labels = pick.genes$gene, 
                                                                               which = "column", labels_gp = gpar(fontsize = 6.25), 
                                                                               side = "bottom"))) 
heat.plt3 <- grid::grid.grabExpr(draw(heat.plt))
pdf(file = file.path(res.dir[1], paste0("heatmap_avg_gexp_top10_positive_coeff_across_divisive_ICP_L", max.sd, ".pdf")), width = 12, height = 3.75)
print(heat.plt2)
dev.off()

# Add plot to final figure
figure1 <- cowplot::plot_grid(fig1AE, heat.plt3, ncol=1, rel_heights = c(0.75, 0.25), labels = c("", "F"))
pdf(file = file.path(res.dir[1], "figure_1.pdf"), width = 12, height = 11.75)
print(figure1)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export final object
saveRDS(object = sce, file = file.path(res.dir[3], "sce_int.rds"))
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
#   [1] Rtsne_0.16                ggbeeswarm_0.7.1          colorspace_2.1-0          rjson_0.2.21              class_7.3-20             
# [6] modeltools_0.2-23         circlize_0.4.15           scuttle_1.8.4             bluster_1.8.0             XVector_0.38.0           
# [11] GlobalOptions_0.1.2       BiocNeighbors_1.16.0      clue_0.3-64               rstudioapi_0.14           farver_2.1.1             
# [16] ggrepel_0.9.3             scatterpie_0.2.3          RSpectra_0.16-1           fansi_1.0.4               codetools_0.2-18         
# [21] R.methodsS3_1.8.2         sparseMatrixStats_1.10.0  doParallel_1.0.17         polyclip_1.10-4           Cairo_1.6-0              
# [26] cluster_2.1.3             png_0.1-8                 R.oo_1.25.0               pheatmap_1.0.12           ggforce_0.4.1            
# [31] HDF5Array_1.26.0          compiler_4.2.1            dqrng_0.3.0               Matrix_1.6-5              limma_3.54.2             
# [36] aricode_1.0.2             cli_3.6.1                 tweenr_2.0.3              BiocSingular_1.14.0       tools_4.2.1              
# [41] rsvd_1.0.5                igraph_1.4.1              gtable_0.3.3              glue_1.6.2                GenomeInfoDbData_1.2.9   
# [46] RANN_2.6.1                reshape2_1.4.4            LiblineaR_2.10-22         doRNG_1.8.6               Rcpp_1.0.10              
# [51] cellranger_1.1.0          vctrs_0.6.1               rhdf5filters_1.10.0       iterators_1.0.14          DelayedMatrixStats_1.20.0
# [56] stringr_1.5.0             beachmat_2.14.0           lifecycle_1.0.3           irlba_2.3.5.1             rngtools_1.5.2           
# [61] statmod_1.5.0             edgeR_3.40.2              MASS_7.3-57               zlibbioc_1.44.0           scales_1.3.0             
# [66] doSNOW_1.0.20             parallel_4.2.1            rhdf5_2.42.0              SparseM_1.81              RColorBrewer_1.1-3       
# [71] ggfun_0.1.3               ggrastr_1.0.1             stringi_1.7.12            foreach_1.5.2             ScaledMatrix_1.6.0       
# [76] scran_1.26.2              BiocParallel_1.32.6       shape_1.4.6               rlang_1.1.0               pkgconfig_2.0.3          
# [81] bitops_1.0-7              lattice_0.20-45           purrr_1.0.1               Rhdf5lib_1.20.0           labeling_0.4.2           
# [86] cowplot_1.1.1             tidyselect_1.2.0          plyr_1.8.8                magrittr_2.0.3            R6_2.5.1                 
# [91] snow_0.4-4                generics_0.1.3            metapod_1.6.0             DelayedArray_0.24.0       pillar_1.9.0             
# [96] withr_2.5.0               RCurl_1.98-1.10           tibble_3.2.1              crayon_1.5.2              DropletUtils_1.18.1      
# [101] utf8_1.2.3                GetoptLong_1.0.5          locfit_1.5-9.7            readxl_1.4.2              flexclust_1.4-1          
# [106] digest_0.6.31             tidyr_1.3.0               R.utils_2.12.2            munsell_0.5.0             beeswarm_0.4.0           
# [111] viridisLite_0.4.1         vipor_0.4.5              
#
#------------------------------------------------------------------------------#
