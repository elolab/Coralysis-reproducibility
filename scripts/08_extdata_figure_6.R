#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make extended data figure 6 and supplementary figure 12. 
# Date: 24/01/2025
# Last update: 24/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("SingleCellExperiment") # v.1.20.1
source("scripts/helper_functions.R")
#require("scran") # v. 1.26.2
#require("scater") # v. 1.26.1
#require("cowplot") # v.1.1.1

# Folders to save results 
analysis <- "08_extdata_figure_6"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Process data 

# Import data
data.dir <- file.path("data", analysis)
input <- list(
  ref = file.path(data.dir, "tran_exp5_pbmc_batch1_balanced.h5ad"), 
  query = file.path(data.dir, "tran_exp5_pbmc_batch2_balanced.h5ad")
)
sce <- lapply(input, function(x) zellkonverter::readH5AD(file = x, X_name = "counts"))

# Normalization
sce <- lapply(sce, NormalizeData)

# Params
batch.label <- "batch"
cell.label <- "celltype"
nhvg <- 2000
color.palette <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                   "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                   "#8C564B", "#C49C94", "#E377C2") 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 6 A

## Check batch effect 
unint <- cbind(sce$ref, sce$query)

# Dimensional reduction 
set.seed(123)
unint <- scran::fixedPCA(unint, rank=30) %>% 
  scater::runUMAP(., dimred="PCA")

# Plot 
use.cols <- color.palette[seq(1, 13, 2)[1:length(unique(unint$celltype))]]
names(use.cols) <- levels(unint$celltype)
use.cols <- list("celltype" = use.cols, "batch" = c("batch_1" = "#1F77B4", "batch_2" = "#FF7F0E"))
vars.plot <- c(batch.label, cell.label)
names(vars.plot) <- vars.plot
unint.umap <- lapply(X = vars.plot, FUN = function(x) {
  PlotDimRed(object = unint, color.by = x, dimred = "UMAP", use.color = use.cols[[x]], 
             point.size = 0.5, point.stroke = 0, legend.nrow = 1) + 
    theme(legend.title = element_blank(), 
          legend.justification = "left", 
          legend.margin = margin(0, 0, 0, 0, "cm"),
          plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm")) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    guides(color = guide_legend(override.aes = list(size = 3.5), nrow = 1))
})
unint.umap.plt <- cowplot::plot_grid(plotlist = unint.umap, ncol = 1, align = "vh", labels = c("A", "")) 
pdf(file = file.path(res.dir[1], "extdata_figure_6_A.pdf"), width = 2.5, height = 5)
print(unint.umap.plt)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 6 B

# Centroids
centroids <- reducedDim(unint, "PCA") %>% 
  as.data.frame(.) %>% cbind(., colData(unint)[,c(batch.label, cell.label)]) %>% 
  mutate("cell_id" = row.names(.)) %>%
  group_by(batch, celltype) %>% 
  summarise_at(vars(PC1:PC30), mean) %>% 
  as.data.frame(.) %>%
  `rownames<-`(paste(.$batch, .$celltype, sep="-"))
distance <- dist(centroids[,-(1:2)])
# tidy
dist.tidy <- distance %>% 
  as.matrix(.) %>% as.data.frame(.) %>% 
  mutate("rows" = row.names(.)) %>% 
  tidyr::pivot_longer(., cols=`batch_1-B cell`:`batch_2-NK cell`, 
                      names_to = "cols", values_to = "dist") %>% 
  mutate("dups" = !duplicated(unlist(lapply(strsplit(paste(rows, cols), " "), function(x) paste(sort(x), collapse="_"))))) %>% 
  filter(dups & (rows!=cols)) %>% 
  tidyr::separate(rows, c("rows_batch", "rows_celltype"), sep="-") %>% 
  tidyr::separate(cols, c("cols_batch", "cols_celltype"), sep="-") %>% 
  select(-dups) %>% 
  mutate("pair" = ifelse(rows_batch==cols_batch, "within", "across"))

# Lollipop plot data
dist.tidy.plot <- dist.tidy %>% 
  filter(., pair=="across") %>%
  group_by(rows_batch, rows_celltype) %>% 
  mutate("mindist" = dist - min(dist), 
         "alpha" = (1 / rank(mindist)+0.25), 
         "meandist" = mean(mindist)) %>% 
  arrange(meandist) %>%
  ungroup() %>% 
  mutate("rows_celltype" = factor(rows_celltype, levels = unique(as.character(rows_celltype))))
# Plot
centroid.cross.batch.dist.plt <- ggplot(data = dist.tidy.plot, 
                                        mapping = aes(x = mindist, y = rows_celltype, color = cols_celltype)) + 
  geom_segment(mapping = aes(x = 0, xend = mindist, y = rows_celltype, yend = rows_celltype), color = "gray90", lwd = 3) + 
  geom_point(alpha = dist.tidy.plot$alpha, size = 5) + 
  scale_color_manual(name = "Cell types (query)", values = use.cols$celltype) +
  geom_segment(mapping = aes(x = meandist, xend = meandist, 
                             y = as.numeric(rows_celltype)-0.25, yend = as.numeric(rows_celltype)+0.25, 
                             linetype = "mean"), 
               color = "#41AB5D", linewidth = 1) +
  scale_linetype_manual(name = "Distance", values=c("mean"=1)) + 
  theme_minimal() + 
  ylab("Cell types (reference)") +
  xlab("Centroid cross-batch Euclidean distance") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) + 
  guides(colour = guide_legend(override.aes = list(size = 3.5)))
centroid.cross.batch.dist.plt2 <- cowplot::plot_grid(centroid.cross.batch.dist.plt, labels = "B")
pdf(file = file.path(res.dir[1], "extdata_figure_6_B.pdf"), width = 7.75, height = 3)
print(centroid.cross.batch.dist.plt2)
dev.off()

# Export tables
write.table(x = dist.tidy, file = file.path(res.dir[2], "centroids_distance_table.tsv"),
            sep ="\t", row.names = FALSE, quote = FALSE)
write.table(x = dist.tidy.plot, file = file.path(res.dir[2], "centroids_distance_table_used_to_plot.tsv"),
            sep ="\t", row.names = FALSE, quote = FALSE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Perform query-to-reference annotation over imbalance data sets
# Query-to-reference mapping - down sample every cell type on the ref to 10% 
cell.types <- levels(unint[[cell.label]])
down <- 0.1 # down sample every cell to 10%, always on the ref data set 
ncells <- ncol(sce$ref)
res <- list()
for (cell in cell.types) {
  name <- paste(cell, down, sep = "_")
  res[[name]] <- list()
  down.label <- which(sce$ref[[cell.label]] == cell)
  other.labels <- which(sce$ref[[cell.label]] != cell)
  set.seed(1024)
  pick.cells <- sample(x = down.label, size = length(down.label)*down)
  pick.cells <- sort(c(pick.cells, other.labels))
  sce$ref[[name]] <- ((1:ncells) %in% pick.cells)
  ref <- sce$ref[,pick.cells]
  set.seed(2024)
  m.hvg <- scran::modelGeneVar(ref)
  hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
  top.hvg <- row.names(ref)[hvg.ordered[1:nhvg]]
  ref <- ref[top.hvg,]
  set.seed(123)
  ref <- RunParallelDivisiveICP(object = ref, verbose = TRUE) # train
  ref <- RunPCA(object = ref, return.model = TRUE) # PCA
  ref <- RunUMAP(ref, return.model = TRUE) # UMAP
  map <- ReferenceMapping(ref = ref, query = sce$query, ref.label = cell.label, project.umap = TRUE) # query-to-ref mapping
  res[[name]][["refdims"]] <- reducedDims(ref)
  res[[name]][["mapdims"]] <- reducedDims(map)
  res[[name]][["annot"]] <- colData(map)
  rm(list = c("ref", "map")); gc(); 
}
saveRDS(object = res, file = file.path(res.dir[3], "res_downsample.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Supplementary Figure 12

# Plot results
acc_fun <- function(truth, pred) { # function to estimate accuracy
  pred <- as.character(pred)
  truth <- as.character(truth)
  conf.mtx <- table(pred, truth)
  acc <- (sum(diag(conf.mtx)) / sum(conf.mtx))
  return(acc)
}
all.down.umap.plts <- down.umap.plts <- list()
pdf(file = file.path(res.dir[1], "supp_figure_12_single_plots.pdf"), width = 8, height = 2.75)
for (cell in names(res)) {
  cell.name <- gsub("_0.1", "", cell)
  acc <- acc_fun(truth = res[[cell]]$annot$celltype, pred = res[[cell]]$annot$coral_labels)
  res[[cell]][["conf_mtx"]] <- table(res[[cell]]$annot$coral_labels, res[[cell]]$annot$celltype)  
  res[[cell]][["acc"]] <- acc 
  res[[cell]][["acc_cellspec"]] <- diag(res[[cell]][["conf_mtx"]] / colSums(res[[cell]][["conf_mtx"]]))
  down.umap.plts[[cell]] <- list()
  down.umap.plts[[cell]][["ref"]] <- res[[cell]]$refdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    merge(x = ., y = as.data.frame(colData(sce$ref[,sce$ref[[cell]]])), by = "row.names", all = TRUE) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = celltype)) +
    geom_point(size = 0.05, alpha = 0.5) +
    ggforce::geom_mark_ellipse(mapping = aes(color = celltype,
                                             filter = celltype == cell.name), 
                               expand = unit(1, "mm")) +
    scale_color_manual(values = use.cols$celltype) +
    theme_minimal() +
    theme(plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
          legend.position = "none", 
          plot.title = element_text(size=9.5, hjust=0.5), 
          axis.title = element_text(size = 10)) +
    ggtitle("reference - ground-truth")
  down.umap.plts[[cell]][["map_truth"]] <- res[[cell]]$mapdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    cbind(., res[[cell]]$annot) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = celltype)) +
    geom_point(size = 0.05, alpha = 0.5) +
    ggforce::geom_mark_ellipse(mapping = aes(color = celltype,
                                             filter = celltype == cell.name), 
                               expand = unit(1, "mm")) +
    scale_color_manual(values = use.cols$celltype, breaks = names(use.cols$celltype), 
                       labels = paste(names(use.cols$celltype), 
                                      paste0("(", round(res[[cell]][["acc_cellspec"]][names(use.cols$celltype)], 3), ")"), 
                                      sep ="\n")) +
    theme_minimal() +
    theme(legend.key.size = unit(0.4, "cm"), 
          plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
          legend.position = "bottom", legend.title = element_blank(), 
          axis.title = element_blank(), 
          legend.text = element_text(size=10), 
          plot.title = element_text(size=9.5, hjust=0.5), 
          legend.margin = margin(0, 0, 0, 0, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3.5, linetype = 0), nrow = 1)) + 
    ggtitle("query - ground-truth")
  down.umap.plts[[cell]][["map_preds"]] <- res[[cell]]$mapdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    cbind(., res[[cell]]$annot) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = coral_labels)) +
    geom_point(size = 0.05, alpha = 0.5) +
    ggforce::geom_mark_ellipse(mapping = aes(color = celltype,
                                             filter = celltype == cell.name), 
                               expand = unit(1, "mm")) +
    scale_color_manual(values = use.cols$celltype) +
    theme_minimal() +
    theme(plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
          legend.position = "none", axis.title = element_blank(), 
          plot.title = element_text(size=9.5, hjust=0.5)) +
    ggtitle("query - predictions")
  down.umap.plts[[cell]][["map_conf"]] <- res[[cell]]$mapdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    cbind(., res[[cell]]$annot) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = coral_probability)) +
    geom_point(size = 0.05) +
    theme_minimal() +
    theme(legend.key.size = unit(0.5, "cm")) +
    scale_colour_continuous(limits = c(0, 1), type = "viridis") +
    theme(legend.position = "bottom", legend.title = element_blank(), 
          legend.text = element_text(size=8), axis.title = element_blank(), 
          plot.title = element_text(size=9.5, hjust=0.5), 
          plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
          legend.margin = margin(0, 0, 0, 0, "cm")) + 
    ggtitle("query - confidence")
  title.name <- paste("Coralysis prediction accuracy (10%", paste0(cell.name, "):"), round(acc, 3))
  all.down.umap.plts[[cell]] <- cowplot::plot_grid(plotlist = down.umap.plts[[cell]], ncol = 4, align = "vh", #axis = "tbrl", 
                                                   labels = title.name, label_x = 0.15, hjust = 0)
  print(
    all.down.umap.plts[[cell]]
  )
}
dev.off()
# Plot altogether
pdf(file = file.path(res.dir[1], "supp_figure_12.pdf"), width = 8, height = 16)
cowplot::plot_grid(plotlist = all.down.umap.plts, labels = c("A", "B", "C", "D", "E", "F"), nrow = 6)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Plot accuracy cell-specific scores
acc.cell.spec.scores <- lapply(res, function(x) x$acc_cellspec)
acc.cell.spec.scores <- data.frame(acc.cell.spec.scores)
acc.cell.spec.scores.tidy <- acc.cell.spec.scores %>% 
  mutate("Cell_type" = row.names(.)) %>% 
  tidyr::pivot_longer(data = ., cols = `B.cell_0.1`:`NK.cell_0.1`, names_to = "Downsample", values_to = "Accuracy") %>%
  group_by(Cell_type) %>% 
  mutate("Downsample" = gsub("0 1", "0.1", gsub("\\.", " ", Downsample)), 
         "alpha" = ifelse(Cell_type==gsub("_0.1", "", Downsample), 1, 0.25)) %>% ungroup(.) %>% 
  left_join(x = ., y = data.frame("Downsample" = names(round(unlist(lapply(res, function(x) x$acc)), 3)), 
                                  "Overall" = round(unlist(lapply(res, function(x) x$acc)), 3)), by = "Downsample") 
cell.spec.acc.scores.plt <- ggplot(data = acc.cell.spec.scores.tidy, mapping = aes(x = Cell_type, y = Accuracy, group = Downsample)) + 
  geom_line(mapping = aes(color = Downsample)) + 
  geom_point(mapping = aes(fill = Cell_type), shape = 21, size = 5, stroke = 1, alpha = acc.cell.spec.scores.tidy$alpha) +
  scale_color_manual(name = "Downsample 10%", values = use.cols$celltype %>% `names<-`(paste0(names(.), "_0.1")), 
                     labels = names(use.cols$celltype), breaks = paste0(names(use.cols$celltype), "_0.1")) +
  scale_fill_manual(values = use.cols$celltype) +
  scale_y_continuous(limits = c(0,1)) + 
  geom_hline(aes(yintercept = Overall, linetype = "dashed", color = Downsample), 
             alpha = 0.75) + 
  scale_linetype_manual(name = "", values = "dashed", labels = "overall accuracy") +
  annotate(geom = "text", x = 6.4, y = range(round(unlist(lapply(res, function(x) x$acc)), 3))+0.03, 
           label = range(round(unlist(lapply(res, function(x) x$acc)), 3)), 
           color = use.cols$celltype[c("CD8 T cell", "B cell")]) + 
  theme_classic() + 
  guides(fill = "none", color = guide_legend(override.aes = list(linewidth = 1.5))) + 
  theme(axis.text = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=15, hjust = 1), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 12)) + 
  ylab("Cell-specific accuracy scores")
cell.spec.acc.scores.plt2 <- cowplot::plot_grid(cell.spec.acc.scores.plt, labels = "C")
pdf(file = file.path(res.dir[1], "extdata_figure_6_C.pdf"), width = 8.5, height = 4.65)
print(cell.spec.acc.scores.plt2)
dev.off()

# Export table
write.table(x = acc.cell.spec.scores.tidy, file = file.path(res.dir[2], "cell_specific_accuracy_scores_table.tsv"), 
            sep ="\t", row.names = FALSE, quote = FALSE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 6 

# Plot altogether 
p1 <- cowplot::plot_grid(unint.umap[[1]] + 
                           theme(legend.key.size = unit(0.4, "cm"), 
                                 plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
                                 legend.position = "bottom", legend.title = element_blank(), 
                                 legend.text = element_text(size=10), 
                                 plot.title = element_text(size=9.5, hjust=0.5), 
                                 legend.margin = margin(0.5, 0, 0, 0, "cm")) +
                           scale_color_manual(values = color.palette[c(1,3)], 
                                              labels = c("reference", "query"), breaks = c("batch_1", "batch_2")) + 
                           guides(color = guide_legend(override.aes = list(size = 3.5), title = NULL, nrow=1)), 
                         unint.umap[[2]] + 
                           theme(legend.key.size = unit(0.4, "cm"), 
                                 plot.margin = margin(1, 0.1, 0.1, 0.1, "cm"), 
                                 legend.position = "bottom", legend.title = element_blank(), 
                                 legend.text = element_text(size=10), 
                                 plot.title = element_text(size=9.5, hjust=0.5), 
                                 legend.margin = margin(0.5, 0, 0, 0, "cm")) +
                           guides(color = guide_legend(override.aes = list(size = 3.5), title = NULL, nrow=1)), 
                         ncol=1, labels = "A")
p2 <- cowplot::plot_grid(NULL, 
                         centroid.cross.batch.dist.plt + 
                           theme(axis.text = element_text(size=10), 
                                 axis.title = element_text(size=12),
                                 legend.title = element_text(size=11),
                                 legend.text = element_text(size=10), 
                                 legend.spacing.y = unit(0, "cm"), 
                                 legend.margin=margin(0, 0, 0.2, 0, unit='cm')),# + 
                         #guides(colour = guide_legend(override.aes = list(size = 3), byrow = TRUE)), 
                         NULL, NULL, NULL, 
                         cell.spec.acc.scores.plt + 
                           theme(axis.text = element_text(size = 10), 
                                 axis.title.x = element_blank(), 
                                 axis.text.x = element_text(angle=11, hjust = 1), 
                                 axis.title = element_text(size = 12), 
                                 legend.text = element_text(size=10), 
                                 legend.title = element_text(size = 11), 
                                 legend.spacing.y = unit(0, "cm"), 
                                 legend.margin=margin(0, 0, 0.2, 0, unit='cm'), 
                                 legend.position = c(0.88, 0.4)), 
                         NULL, NULL,
                         rel_heights = c(0.25, 0.001, 0.3, 0.02), rel_widths = c(0.02, 0.8),
                         ncol = 2, labels = c("", "B", "", "", "", "C", "", ""), 
                         scale = 1)
p3 <- cowplot::plot_grid(all.down.umap.plts$`CD8 T cell_0.1`, labels = "D")
p4 <- cowplot::plot_grid(all.down.umap.plts$`B cell_0.1`, labels = "E")
extdata.fig.6 <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(0.25, 0.75)), 
                                    p3, p4, rel_heights = c(0.5, 0.25, 0.25), nrow = 3) 
pdf(file.path(res.dir[1], "extdata_figure_6.pdf"), width = 10, height = 13)
print(extdata.fig.6)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" 
saveRDS(object = sce, file = file.path(res.dir[3], "sce.rds"))

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.1               SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [5] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
# [9] BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_1.1.0           Coralysis_1.0.0            
# [13] ggplot2_3.5.1               dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.3                spatstat.explore_3.1-0    reticulate_1.34.0         R.utils_2.12.2            tidyselect_1.2.0         
# [6] htmlwidgets_1.6.2         grid_4.2.1                BiocParallel_1.32.6       Rtsne_0.16                DropletUtils_1.18.1      
# [11] zellkonverter_1.8.0       munsell_0.5.0             ScaledMatrix_1.6.0        ica_1.0-3                 codetools_0.2-18         
# [16] statmod_1.5.0             scran_1.26.2              umap_0.2.10.0             miniUI_0.1.1.1            future_1.32.0            
# [21] withr_2.5.0               spatstat.random_3.1-4     colorspace_2.1-0          progressr_0.13.0          filelock_1.0.2           
# [26] rstudioapi_0.14           Seurat_4.3.0              ROCR_1.0-11               tensor_1.5                listenv_0.9.0            
# [31] labeling_0.4.2            GenomeInfoDbData_1.2.9    polyclip_1.10-4           farver_2.1.1              pheatmap_1.0.12          
# [36] rhdf5_2.42.0              basilisk_1.10.2           parallelly_1.35.0         vctrs_0.6.1               generics_0.1.3           
# [41] R6_2.5.1                  ggbeeswarm_0.7.1          rsvd_1.0.5                locfit_1.5-9.7            spatstat.utils_3.1-1     
# [46] bitops_1.0-7              rhdf5filters_1.10.0       DelayedArray_0.24.0       promises_1.2.0.1          scales_1.3.0             
# [51] beeswarm_0.4.0            gtable_0.3.3              beachmat_2.14.0           Cairo_1.6-0               globals_0.16.2           
# [56] goftest_1.2-3             spam_2.10-0               rlang_1.1.0               GlobalOptions_0.1.2       splines_4.2.1            
# [61] lazyeval_0.2.2            spatstat.geom_3.1-0       abind_1.4-5               reshape2_1.4.4            httpuv_1.6.9             
# [66] gridtext_0.1.5            tools_4.2.1               ellipsis_0.3.2            RColorBrewer_1.1-3        ggridges_0.5.4           
# [71] Rcpp_1.0.10               plyr_1.8.8                sparseMatrixStats_1.10.0  zlibbioc_1.44.0           purrr_1.0.1              
# [76] RCurl_1.98-1.10           basilisk.utils_1.10.0     deldir_1.0-6              openssl_2.0.6             pbapply_1.7-0            
# [81] viridis_0.6.2             zoo_1.8-11                SeuratObject_5.0.2        ggrepel_0.9.3             cluster_2.1.3            
# [86] magrittr_2.0.3            scattermore_0.8           data.table_1.14.8         RSpectra_0.16-1           SparseM_1.81             
# [91] circlize_0.4.15           lmtest_0.9-40             RANN_2.6.1                fitdistrplus_1.1-8        flexclust_1.4-1          
# [96] patchwork_1.2.0           xtable_1.8-4              mime_0.12                 readxl_1.4.2              gridExtra_2.3            
# [101] shape_1.4.6               compiler_4.2.1            scater_1.26.1             tibble_3.2.1              KernSmooth_2.23-20       
# [106] crayon_1.5.2              htmltools_0.5.4           R.oo_1.25.0               later_1.3.0               snow_0.4-4               
# [111] tidyr_1.3.0               ggtext_0.1.2              DBI_1.1.3                 tweenr_2.0.3              MASS_7.3-57              
# [116] LiblineaR_2.10-22         Matrix_1.6-5              cli_3.6.1                 R.methodsS3_1.8.2         parallel_4.2.1           
# [121] metapod_1.6.0             dotCall64_1.1-1           igraph_1.4.1              pkgconfig_2.0.3           dir.expiry_1.6.0         
# [126] sp_1.6-0                  spatstat.sparse_3.0-1     plotly_4.10.1             scuttle_1.8.4             xml2_1.3.3               
# [131] foreach_1.5.2             vipor_0.4.5               dqrng_0.3.0               rngtools_1.5.2            XVector_0.38.0           
# [136] doRNG_1.8.6               stringr_1.5.0             digest_0.6.31             sctransform_0.3.5         RcppAnnoy_0.0.20         
# [141] spatstat.data_3.0-1       leiden_0.4.3              cellranger_1.1.0          uwot_0.1.14               edgeR_3.40.2             
# [146] DelayedMatrixStats_1.20.0 shiny_1.7.4               modeltools_0.2-23         nlme_3.1-157              lifecycle_1.0.3          
# [151] aricode_1.0.2             jsonlite_1.8.4            Rhdf5lib_1.20.0           BiocNeighbors_1.16.0      viridisLite_0.4.1        
# [156] askpass_1.1               limma_3.54.2              fansi_1.0.4               pillar_1.9.0              lattice_0.20-45          
# [161] fastmap_1.1.1             httr_1.4.5                ggrastr_1.0.1             survival_3.3-1            glue_1.6.2               
# [166] FNN_1.1.3.2               png_0.1-8                 iterators_1.0.14          bluster_1.8.0             ggforce_0.4.1            
# [171] class_7.3-20              stringi_1.7.12            HDF5Array_1.26.0          BiocSingular_1.14.0       doSNOW_1.0.20            
# [176] irlba_2.3.5.1             future.apply_1.10.0   
#
#------------------------------------------------------------------------------#