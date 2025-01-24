#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make extended data figure 7 and supplementary figure 13. 
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
analysis <- "09_extdata_figure_7"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Process data 

# Import data
input <- "results/08_extdata_figure_6/objects/sce.rds"
sce <- readRDS(file = input)

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
## Extended Data Figure 7 A

## Check batch effect 
colData(sce$ref) <- colData(sce$ref)[,1:10]
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
pdf(file = file.path(res.dir[1], "extdata_figure_7_A.pdf"), width = 2.5, height = 5)
print(unint.umap.plt)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 7 B

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
# lollipop plot
dist.tidy.plot <- dist.tidy %>% 
  filter(., pair=="across") %>%
  group_by(rows_batch, rows_celltype) %>% 
  mutate("mindist" = dist - min(dist), 
         "alpha" = (1 / rank(mindist)+0.25), 
         "meandist" = mean(mindist)) %>% 
  arrange(meandist) %>%
  ungroup() %>% 
  mutate("rows_celltype" = factor(rows_celltype, levels = unique(as.character(rows_celltype))))
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
pdf(file = file.path(res.dir[1], "extdata_figure_7_B.pdf"), width = 7.75, height = 3)
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
## Perform query-to-reference annotation over data sets with unshared/ablated cell types
# Query-to-reference mapping - down sample every cell type on the ref to 10% 
cell.types <- levels(unint[[cell.label]])
ncells <- ncol(sce$ref)
res <- list()
for (cell in cell.types) {
  name <- paste(cell, "ablated", sep = "_")
  res[[name]] <- list()
  pick.cells <- which(sce$ref[[cell.label]] != cell)
  sce$ref[[name]] <- ((1:ncells) %in% pick.cells)
  ref <- sce$ref[,pick.cells]
  set.seed(2024)
  m.hvg <- scran::modelGeneVar(ref)
  hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
  top.hvg <- row.names(ref)[hvg.ordered[1:nhvg]]
  ref <- ref[top.hvg,]
  set.seed(123)
  ref <- RunParallelDivisiveICP(object = ref, threads = 4) # train
  ref <- RunPCA(object = ref, return.model = TRUE) # PCA
  ref <- RunUMAP(ref, return.model=TRUE) # UMAP
  map <- ReferenceMapping(ref = ref, query = sce$query, ref.label = cell.label, project.umap = TRUE) # query-to-ref mapping
  res[[name]][["refdims"]] <- reducedDims(ref)
  res[[name]][["mapdims"]] <- reducedDims(map)
  res[[name]][["annot"]] <- colData(map)
  rm(list = c("ref", "map")); gc(); 
}
saveRDS(object = res, file = file.path(res.dir[3], "res_ablated.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Supplementary Figure 12

# Plot results
acc_fun <- function(truth, pred) { # function to estimate accuracy
  pred <- as.character(pred)
  truth <- as.character(truth)
  uniq.truth <- unique(truth)
  uniq.pred <- unique(pred)
  labels <- union(uniq.truth, uniq.pred)
  truth <- factor(truth, levels = labels)
  pred <- factor(pred, levels = labels)
  conf.mtx <- table(pred, truth)
  acc <- (sum(diag(conf.mtx)) / sum(conf.mtx))
  return(acc)
}
all.ablated.umap.plts <- ablated.umap.plts <- list()
pdf(file = file.path(res.dir[1], "supp_figure_13_single_plots.pdf"), width = 8, height = 2.75)
for (cell in names(res)) {
  cell.name <- gsub("_ablated", "", cell)
  acc <- acc_fun(truth = res[[cell]]$annot$celltype, pred = res[[cell]]$annot$coral_labels)
  res[[cell]][["conf_mtx"]] <- table(res[[cell]]$annot$coral_labels, res[[cell]]$annot$celltype)  
  res[[cell]][["acc"]] <- acc 
  res[[cell]][["acc_cellspec"]] <- diag(res[[cell]][["conf_mtx"]] / colSums(res[[cell]][["conf_mtx"]]))
  ablated.umap.plts[[cell]] <- list()
  ablated.umap.plts[[cell]][["ref"]] <- res[[cell]]$refdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    merge(x = ., y = as.data.frame(colData(sce$ref[,sce$ref[[cell]]])), by = "row.names", all = TRUE) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = celltype)) +
    geom_point(size = 0.05, alpha = 0.5) +
    scale_color_manual(values = use.cols$celltype) +
    theme_minimal() +
    theme(plot.margin = margin(0.2, 0, 0, 0.25, "cm"),  
          legend.position = "none", 
          plot.title = element_text(size=9.5, hjust=0.5), 
          axis.title = element_text(size = 10)) +
    ggtitle("reference - ground-truth")
  ablated.umap.plts[[cell]][["map_truth"]] <- res[[cell]]$mapdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    cbind(., res[[cell]]$annot) %>% 
    mutate("filter" = celltype == cell.name) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = celltype)) +
    geom_point(size = 0.05, alpha = 0.5) +
    ggforce::geom_mark_ellipse(mapping = aes(color = celltype,
                                             filter = filter), 
                               expand = unit(1, "mm")) +
    scale_color_manual(values = use.cols$celltype, breaks = names(use.cols$celltype), 
                       labels = paste(names(use.cols$celltype), 
                                      paste0("(", round(res[[cell]][["acc_cellspec"]][names(use.cols$celltype)], 3), ")"), 
                                      sep ="\n")) +
    theme_minimal() +
    theme(legend.key.size = unit(0.4, "cm"), 
          plot.margin = margin(0.2, 0, 0, 0.25, "cm"), 
          legend.position = "bottom", legend.title = element_blank(), 
          axis.title = element_blank(), 
          legend.text = element_text(size=10), 
          plot.title = element_text(size=9.5, hjust=0.5), 
          legend.margin = margin(0, 0, 0, 0, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3.5, linetype = 0), nrow = 1)) + 
    ggtitle("query - ground-truth")
  bar.plt <- res[[cell]]$conf_mtx[,cell.name] %>% 
    data.frame("cell_type" = names(.), "freq" = ., "type" = cell.name) %>% 
    filter(., freq != 0) %>% 
    ggplot(data = ., mapping = aes(x = type, y = (freq / sum(freq) * 100), fill = cell_type)) + 
    geom_col(alpha=0.75) + 
    scale_fill_manual(values = use.cols$celltype) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_classic() + 
    coord_flip() + ylab("") + 
    theme(axis.line = element_blank(), axis.title = element_blank(), 
          axis.text.y = element_blank(), legend.position = "none", 
          axis.ticks.y = element_blank(), axis.text.x = element_text(size=6), 
          plot.background = element_rect(fill = "transparent", colour = NA), 
          plot.margin = margin(0,0,0,0, unit="cm")) 
  ablated.umap.plts[[cell]][["map_preds"]] <- res[[cell]]$mapdims$UMAP %>% 
    as.data.frame(.) %>% `colnames<-`(c(paste0("UMAP", 1:2))) %>% 
    cbind(., res[[cell]]$annot) %>% 
    mutate("filter" = celltype == cell.name) %>% 
    ggplot(data = .,
           mapping = aes(x = UMAP1, y = UMAP2, color = coral_labels)) +
    geom_point(size = 0.05, alpha = 0.5) +
    ggforce::geom_mark_ellipse(mapping = aes(color = celltype,
                                             filter = filter), 
                               expand = unit(1, "mm")) +
    scale_color_manual(values = use.cols$celltype) +
    theme_minimal() +
    theme(plot.margin = margin(0.2, 0, 0, 0.25, "cm"), 
          legend.position = "none", axis.title = element_blank(), 
          plot.title = element_text(size=9.5, hjust=0.5), 
          plot.background = element_rect(fill = "transparent", colour = NA)) +
    ggtitle("query - predictions")
  ablated.umap.plts[[cell]][["map_preds"]] <- ablated.umap.plts[[cell]][["map_preds"]] + 
    patchwork::inset_element(bar.plt, left = 0.6, bottom = 0.8, right = 1, top = 1) 
  ablated.umap.plts[[cell]][["map_conf"]] <- res[[cell]]$mapdims$UMAP %>% 
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
          plot.margin = margin(0.2, 0, 0, 0.25, "cm"), 
          legend.margin = margin(0, 0, 0, 0, "cm")) + 
    ggtitle("query - confidence")
  title.name <- paste("Coralysis prediction accuracy", paste0("(", cell.name, " ablated):"), round(acc, 3))
  all.ablated.umap.plts[[cell]] <- (ablated.umap.plts[[cell]]$ref + ablated.umap.plts[[cell]]$map_truth + 
                                      ablated.umap.plts[[cell]]$map_preds + ablated.umap.plts[[cell]]$map_conf + 
                                      patchwork::plot_annotation(title = title.name, theme = theme(plot.title = element_text(hjust=0.1))) + 
                                      patchwork::plot_layout(ncol = 4, widths = 1))
  print(
    all.ablated.umap.plts[[cell]]
  )
}
dev.off()
pdf(file = file.path(res.dir[1], "supp_figure_13.pdf"), width = 8, height = 16)
print(cowplot::plot_grid(plotlist = all.ablated.umap.plts, 
                         labels = c("A", "B", "C", "D", "E", "F"), nrow = 6))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 7 C

## Plot accuracy cell-specific scores
acc.cell.spec.scores <- lapply(res, function(x) x$acc_cellspec)
acc.cell.spec.scores <- data.frame(acc.cell.spec.scores)
acc.cell.spec.scores.tidy <- acc.cell.spec.scores %>% 
  mutate("Cell_type" = row.names(.)) %>% 
  tidyr::pivot_longer(data = ., cols = `B.cell_ablated`:`NK.cell_ablated`, names_to = "Ablated", values_to = "Accuracy") %>%
  group_by(Cell_type) %>% 
  mutate("Ablated" = gsub("\\.", " ", Ablated), 
         "alpha" = ifelse(Cell_type==gsub("_ablated", "", Ablated), 1, 0.25)) %>% ungroup(.) %>% 
  left_join(x = ., y = data.frame("Ablated" = names(round(unlist(lapply(res, function(x) x$acc)), 3)), 
                                  "Overall" = round(unlist(lapply(res, function(x) x$acc)), 3)), by = "Ablated") 
cell.spec.acc.scores.plt <- ggplot(data = acc.cell.spec.scores.tidy, mapping = aes(x = Cell_type, y = Accuracy, group = Ablated)) + 
  geom_line(mapping = aes(color = Ablated)) + 
  geom_point(mapping = aes(fill = Cell_type), shape = 21, size = 5, stroke = 1, alpha = acc.cell.spec.scores.tidy$alpha) +
  scale_color_manual(name = "Ablated", values = use.cols$celltype %>% `names<-`(paste0(names(.), "_ablated")), 
                     labels = names(use.cols$celltype), breaks = paste0(names(use.cols$celltype), "_ablated")) +
  scale_fill_manual(values = use.cols$celltype) +
  scale_y_continuous(limits = c(0,1)) + 
  geom_hline(aes(yintercept = Overall, linetype = "dashed", color = Ablated), 
             alpha = 0.75) + 
  scale_linetype_manual(name = "", values = "dashed", labels = "overall accuracy") +
  annotate(geom = "text", x = 6.4, y = range(round(unlist(lapply(res, function(x) x$acc)), 3))+c(-0.03,0.03),
           label = range(round(unlist(lapply(res, function(x) x$acc)), 3)),
           color = use.cols$celltype[c("B cell", "Monocyte_CD14")]) +
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
pdf(file = file.path(res.dir[1], "extdata_figure_7_C.pdf"), width = 8.5, height = 4.65)
print(cell.spec.acc.scores.plt2)
dev.off()

# Export table
write.table(x = acc.cell.spec.scores.tidy, file = file.path(res.dir[2], "cell_specific_accuracy_scores_table.tsv"), 
            sep ="\t", row.names = FALSE, quote = FALSE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 7

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
                                 legend.position = c(0.88, 0.4), 
                                 legend.key = element_blank(), 
                                 legend.background = element_blank()), 
                         NULL, NULL,
                         rel_heights = c(0.25, 0.001, 0.3, 0.02), rel_widths = c(0.02, 0.8),
                         ncol = 2, labels = c("", "B", "", "", "", "C", "", ""), 
                         scale = 1)
p3 <- cowplot::plot_grid(all.ablated.umap.plts$`B cell_ablated`, labels = "D")
p4 <- cowplot::plot_grid(all.ablated.umap.plts$Monocyte_CD14_ablated, labels = "E")
pdf(file.path(res.dir[1], "extdata_figure_7.pdf"), width = 10, height = 13)
cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(0.25, 0.75)), 
                   p3, p4, rel_heights = c(0.5, 0.25, 0.25), nrow = 3) 
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
#   [1] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
# [5] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
# [9] MatrixGenerics_1.10.0       matrixStats_1.1.0           Coralysis_1.0.0             ggplot2_3.5.1              
# [13] dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.7.1          colorspace_2.1-0          class_7.3-20              modeltools_0.2-23         scuttle_1.8.4            
# [6] bluster_1.8.0             XVector_0.38.0            BiocNeighbors_1.16.0      rstudioapi_0.14           farver_2.1.1             
# [11] ggrepel_0.9.3             RSpectra_0.16-1           fansi_1.0.4               codetools_0.2-18          sparseMatrixStats_1.10.0 
# [16] scater_1.26.1             polyclip_1.10-4           jsonlite_1.8.4            Cairo_1.6-0               umap_0.2.10.0            
# [21] cluster_2.1.3             png_0.1-8                 pheatmap_1.0.12           uwot_0.1.14               ggforce_0.4.1            
# [26] compiler_4.2.1            dqrng_0.3.0               Matrix_1.6-5              limma_3.54.2              aricode_1.0.2            
# [31] cli_3.6.1                 tweenr_2.0.3              BiocSingular_1.14.0       tools_4.2.1               rsvd_1.0.5               
# [36] igraph_1.4.1              gtable_0.3.3              glue_1.6.2                GenomeInfoDbData_1.2.9    RANN_2.6.1               
# [41] reshape2_1.4.4            LiblineaR_2.10-22         doRNG_1.8.6               Rcpp_1.0.10               vctrs_0.6.1              
# [46] iterators_1.0.14          DelayedMatrixStats_1.20.0 stringr_1.5.0             beachmat_2.14.0           lifecycle_1.0.3          
# [51] irlba_2.3.5.1             rngtools_1.5.2            statmod_1.5.0             edgeR_3.40.2              zlibbioc_1.44.0          
# [56] MASS_7.3-57               scales_1.3.0              doSNOW_1.0.20             parallel_4.2.1            SparseM_1.81             
# [61] RColorBrewer_1.1-3        reticulate_1.34.0         gridExtra_2.3             ggrastr_1.0.1             stringi_1.7.12           
# [66] foreach_1.5.2             ScaledMatrix_1.6.0        scran_1.26.2              BiocParallel_1.32.6       rlang_1.1.0              
# [71] pkgconfig_2.0.3           bitops_1.0-7              lattice_0.20-45           purrr_1.0.1               patchwork_1.2.0          
# [76] labeling_0.4.2            cowplot_1.1.1             tidyselect_1.2.0          plyr_1.8.8                magrittr_2.0.3           
# [81] R6_2.5.1                  snow_0.4-4                generics_0.1.3            metapod_1.6.0             DelayedArray_0.24.0      
# [86] pillar_1.9.0              withr_2.5.0               RCurl_1.98-1.10           tibble_3.2.1              crayon_1.5.2             
# [91] utf8_1.2.3                viridis_0.6.2             locfit_1.5-9.7            grid_4.2.1                flexclust_1.4-1          
# [96] FNN_1.1.3.2               digest_0.6.31             tidyr_1.3.0               openssl_2.0.6             munsell_0.5.0            
# [101] beeswarm_0.4.0            viridisLite_0.4.1         vipor_0.4.5               askpass_1.1    
#
#------------------------------------------------------------------------------#