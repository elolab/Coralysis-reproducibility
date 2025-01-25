#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make figure 5. 
# Date: 25/01/2025
# Last update: 25/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("SingleCellExperiment") # v.1.20.1
#require("scran") # v. 1.26.2
#require("ILoReg") # v.0.99.6
#require("scater") # v.1.26.1
#require("cowplot") # v.1.1.1

# Folders to save results 
analysis <- "11_figure_5"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)

# Import data 
input <- file.path("data", analysis, "sce.rds")
sce <- readRDS(file = input)

# Params
batch.label <- "batch"
cell.label <- "cell_type"
nhvg <- 2000
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Process data 

# Split object into 'ref' ('10x Chromium (v2) A)' versus 'query' (remaining batches)
sce[[batch.label]] <- factor(sce[[batch.label]])
batches <- levels(sce[[batch.label]])
ref.cells <- (sce[[batch.label]] %in% batches[2]) 
query.cells <- (!ref.cells)
ref <- sce[,ref.cells]
query <- sce[,query.cells]
ref[[batch.label]] <- factor(x = ref[[batch.label]], 
                             levels = batches[2])
query <- list(
  "chromium10xV2" =  query[,query$batch==batches[1]],
  "chromium10xV2B" =  query[,query$batch==batches[3]],
  "chromium10xV3" =  query[,query$batch==batches[4]],
  "celseq2" =  query[,query$batch==batches[5]], 
  "dropseq" =  query[,query$batch==batches[6]],
  "inDrops" =  query[,query$batch==batches[7]],
  "seqwell" =  query[,query$batch==batches[8]], 
  "smartseq2" =  query[,query$batch==batches[9]]
)
query[[1]][[batch.label]] <- factor(x = query[[1]][[batch.label]], levels = batches[1])
query[[2]][[batch.label]] <- factor(x = query[[2]][[batch.label]], levels = batches[3])
query[[3]][[batch.label]] <- factor(x = query[[3]][[batch.label]], levels = batches[4])
query[[4]][[batch.label]] <- factor(x = query[[4]][[batch.label]], levels = batches[5])
query[[5]][[batch.label]] <- factor(x = query[[5]][[batch.label]], levels = batches[6])
query[[6]][[batch.label]] <- factor(x = query[[6]][[batch.label]], levels = batches[7])
query[[7]][[batch.label]] <- factor(x = query[[7]][[batch.label]], levels = batches[8])
query[[8]][[batch.label]] <- factor(x = query[[8]][[batch.label]], levels = batches[9])

# Colors
tableau10medium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                     "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                     "#CDCC5D", "#6DCCDA")
use.cols <- c(tableau10medium[seq_along(batches)], tableau10medium[1:length(unique(sce$cell_type))])
names(use.cols) <- c(batches, unique(sce$cell_type))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Train reference

# HVG selection - reference object
m.hvg <- scran::modelGeneVar(ref)
hvg.ordered <- order(m.hvg[["bio"]], decreasing=TRUE)
top.hvg <- row.names(ref)[hvg.ordered[1:nhvg]]
rowData(ref) <- cbind("highly_variable" = row.names(ref) %in% top.hvg, m.hvg)
ref <- ref[top.hvg,]

# Train (integrate) reference
set.seed(1024)
ref <- RunParallelDivisiveICP(object = ref, threads = 4)

# PCA
ref <- RunPCA(object = ref, return.model = TRUE)

set.seed(123)
ref <- RunUMAP(ref, return.model=TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Reference mapping 
query.map <- lapply(X = query, FUN = function(x) {
  ReferenceMapping(ref = ref, query = x, ref.label = cell.label, project.umap = TRUE)
})
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Figure 5 A-C

# Concatenate SCE objects
ref.query <- do.call(cbind, query.map)
colData(ref)[,c("coral_labels", "coral_probability")] <- NA
ref$coral_labels <- ref$cell_type
ref.query <- cbind(ref, ref.query[row.names(ref),])
ref.query[["type"]] <- as.character(ref.query$batch)
ref.query[["type"]][ref.query[["type"]]==batches[2]] <- "reference"
ref.query[["type"]] <- factor(ref.query[["type"]], levels = c("reference", levels(ref.query$batch)[levels(ref.query$batch) != batches[2]]))
ref.query$cell_type <- as.factor(ref.query$cell_type)
stopifnot(all(sort(levels(ref.query$cell_type)) == sort(levels(ref.query$coral_labels))))
ref.query$coral_labels <- factor(ref.query$coral_labels, levels = levels(ref.query$cell_type))
ref.query$Type <- ifelse(ref.query$type == "reference", "Reference", "Queries")
ref.query$Type <- factor(ref.query$Type, levels = c("Reference", "Queries"))

# Plot
vars.plot <- c("type", cell.label, "coral_labels")
names(vars.plot) <- vars.plot
title <- c("Reference-Queries", "Ground-truth", "Transfer labels")
names(title) <- vars.plot
ref.query.umaps <- lapply(X = vars.plot, FUN = function(x) {
  p <- PlotDimRed(object = ref.query, color.by = x, dimred = "UMAP", point.size = 0.5, point.stroke = 0) + 
    facet_grid(~Type) + 
    theme(strip.background = element_rect(color = "white"), 
          strip.text.x = element_text(size = 12))
  if (x == "type") {
    p <- p + scale_color_manual(values = RColorBrewer::brewer.pal(9, "Paired")[1:9]) 
  } else {
    p <- p + scale_color_manual(values = use.cols[levels(ref.query$cell_type)])
  }
  p + guides(color = guide_legend(title = title[x], override.aes = list(size = 3.5), nrow = 5))
})
panels_A_C <- cowplot::plot_grid(plotlist = ref.query.umaps, align = "vh", labels = LETTERS[1:3], nrow = 1)
pdf(file = file.path(res.dir[1], "final_5_A_C.pdf"), width = 11.5, height = 3.75)
print(panels_A_C)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 5 D

## Explore misclassifications
query.coldata <- as.data.frame(colData(ref.query)) %>%
  filter(batch %in% batches[-2]) %>% # remove reference
  mutate("Correct_classification" = (cell_type == as.character(coral_labels))) %>%
  mutate("Correct_classification" = factor(ifelse(Correct_classification, "Correct", "Incorrect")))

# Accuracy and precision
query.coldata$cell_type <- factor(query.coldata$cell_type, levels = levels(query.coldata$coral_labels))
acc.prec.data <- split(x = query.coldata[,c("cell_type", "coral_labels")],
                       f = as.character(query.coldata$batch))
conf.mtx <- lapply(X = acc.prec.data, FUN = function(x) {
  table(x$coral_labels, x$cell_type)
})
acc <- lapply(X = conf.mtx, FUN = function(x) caret::confusionMatrix(x)$overall[1])

# Export table
write.table(x = cbind("cell_id" = row.names(query.coldata), query.coldata), 
            file = file.path(res.dir[2], "query_coldata.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE)

## Barplot of accuracy per query dataset
acc.query.data <- data.frame("Query" = names(acc), "Accuracy" = unlist(acc)) %>% 
  mutate("Query" = factor(Query, levels = gsub(".Accuracy", "", names(sort(unlist(acc), decreasing = TRUE)))))
acc.query.barplot <-  ggplot(data = acc.query.data, mapping = aes(x = Query, y = Accuracy)) + 
  geom_col(fill = RColorBrewer::brewer.pal(9, "Paired")[2:9], alpha = 0.75) + 
  geom_hline(mapping = aes(yintercept = mean(Accuracy), linetype = "Mean"), color = "black", linewidth = 0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, angle = 25, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 14)) + 
  scale_linetype_manual(name = "", values = 2) +
  ylim(0, 1) + 
  ylab("Overall accuracy\nper query dataset") + 
  annotate(geom = "text", x = 8.25, y = round(mean(acc.query.data$Accuracy), 2)+0.07, 
           label = round(mean(acc.query.data$Accuracy), 2))
pdf(file = file.path(res.dir[1], "figure_5_D.pdf"), width = 10, height = 5)
cowplot::plot_grid(acc.query.barplot, labels = "D")
dev.off()
write.table(x = acc.query.data, file = file.path(res.dir[2], "accuracy_scores_per_query_dataset_table.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 5 D

## Density plots showing the confidence score per query dataset per classification
class.query.density.data <- query.coldata %>% 
  mutate("batch" = factor(batch, levels = levels(acc.query.barplot$data$Query)))
class.query.density.plot <- ggplot(data = class.query.density.data, 
                                   mapping = aes(x = coral_probability, fill = Correct_classification)) + 
  geom_density(alpha = 0.75) +
  facet_grid(~batch) + 
  scale_fill_brewer(name = "Classification", palette = "Accent") + 
  theme_minimal() + 
  ylab("Density") + 
  xlab("Confidence scores<br />
       (proportion of *K* neighbors from the winning class)") +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        axis.title.x = ggtext::element_markdown(), 
        panel.spacing = unit(1, "lines")) + 
  scale_x_continuous(breaks = c(0.0, 0.5, 1), limits = c(0, 1))
pdf(file = file.path(res.dir[1], "figure_5_E.pdf"), width = 13, height = 2.25)
cowplot::plot_grid(class.query.density.plot, labels = "E")
dev.off()
write.table(x = cbind("Cell_ID" = row.names(class.query.density.data), class.query.density.data), 
            file = file.path(res.dir[2], "metadata_table.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 5 F

## Dotplot of PBMC marker genes for predicted cell types - only for seqwell dataset
pbmc.markers <- list(
  "B cell" = c("CD19", "MS4A1"), 
  "CD14+ monocytes" = c("CD14"),
  "CD16+ monocytes" = c("FCGR3A"), 
  "CD4+ T cell" = c("CD3D", "CD4"), 
  "Cytotoxic T cell" = c("CD8A", "TRAC", "TRBC1", "TRBC2"), 
  "Dendritic cell" = c("ITGAX", "HLA-DRA", "HLA-DRB1", "HLA-DRB5"),
  "Megakaryocyte" =  c("ITGA2B", "ITGB3"),
  "Natural killer" = c("KIR2DL1", "KIR2DL3", "KLRK1", "FCGR3A"), 
  "Plasmacytoid dendritic cell" = c("IL3RA", "CLEC4C")
)
query.map$seqwell$coral_labels <- factor(query.map$seqwell$coral_labels, levels = rev(levels(query.map$seqwell$coral_labels)))
pbmc.pred.cell_type.markers.dotplot <- scater::plotDots(object = query.map$seqwell, 
                                                        features = unique(unlist(pbmc.markers)), 
                                                        group = "coral_labels", scale = TRUE, center = TRUE) + 
  coord_flip() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1), 
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 14)) + 
  xlab("Predicted cell types") + 
  ylab("PBMC marker genes") + 
  ggtitle("Query seqwell") +
  guides(size = guide_legend(title = "Proportion\nExpressed"))
pdf(file = file.path(res.dir[1], "figure_5_F.pdf"),  width = 11, height = 3.75)
cowplot::plot_grid(pbmc.pred.cell_type.markers.dotplot, labels = "F")
dev.off()
query.map$seqwell$coral_labels <- factor(query.map$seqwell$coral_labels, levels = rev(levels(query.map$seqwell$coral_labels)))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 5 I

## Plot marker genes distinguishing CD14 vs CD16 monocytes and cytotoxic T cells vs NK cells (seqwell data)
marker.2.plt <- c("CD14", "FCGR3A", "CD3D", "KIR2DL3") 
markers.plt <- ILoReg::GeneScatterPlot(object = query.map$seqwell, genes = marker.2.plt, 
                                       return.plot = TRUE, dim.reduction.type = "umap", 
                                       point.size = 0.01, plot.expressing.cells.last = TRUE, ncol = 4)
markers.plt2 <- lapply(X = marker.2.plt, FUN = function(x) {
  PlotExpression(object = query.map$seqwell, color.by = x, point.size = 0.25, point.stroke = 0.25)
}) 
pdf(file = file.path(res.dir[1], "figure_5_I_alternative.pdf"), width = 12, height = 2.1)
cowplot::plot_grid(markers.plt, labels = "I")
dev.off()
pdf(file = file.path(res.dir[1], "figure_5_I.pdf"), width = 12, height = 2.1)
cowplot::plot_grid(plotlist = markers.plt2, labels = c("I", "", "", ""), ncol = 4)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 5 G-H

## Project seqwell data alone
query.map$seqwell$cell_type <- factor(query.map$seqwell$cell_type, levels = levels(ref.query$cell_type)[levels(ref.query$cell_type) %in% unique(query.map$seqwell$cell_type)])
seqwell.ground.truth.umap <- PlotDimRed(object = query.map$seqwell, color.by = cell.label, 
                                        dimred = "UMAP", point.size = 0.5, point.stroke = 0.25, 
                                        use.color = use.cols[levels(query.map$seqwell$cell_type)]) + 
  theme(legend.position = "right", 
        legend.key.size = unit(0.4, "cm"),
        legend.text=element_text(size=10.5),
        legend.title = element_text(size = 12.5, hjust = 0.5)) +
  ggtitle("Query seqwell: original cell labels") +
  guides(color = guide_legend(title = "Ground-truth", nrow = 5, override.aes = list(size = 3.5)))
seqwell.preds.umap <- PlotDimRed(object = query.map$seqwell, color.by = "coral_labels", 
                                 dimred = "UMAP", point.size = 0.5, point.stroke = 0.25, 
                                 use.color = use.cols[levels(query.map$seqwell$coral_labels)]) + 
  theme(legend.position = "right", 
        legend.key.size = unit(0.4, "cm"),
        legend.text=element_text(size=10.5),
        legend.title = element_text(size = 12.5, hjust = 0.5)) +
  ggtitle("Query seqwell: Coralysis's cell label predictions") +
  guides(color = guide_legend(title = "Transfer labels", nrow = 5, override.aes = list(size = 3.5)))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Assemble final figure 
panels_D_F <- cowplot::plot_grid(acc.query.barplot, class.query.density.plot, pbmc.pred.cell_type.markers.dotplot, 
                                 ncol = 1, labels = LETTERS[4:6], rel_heights = c(0.35, 0.225, 0.425))
panels_G_I <- cowplot::plot_grid(cowplot::plot_grid(seqwell.ground.truth.umap,  
                                                    seqwell.preds.umap,  
                                                    ncol = 2,  align = "hv", axis = "tblr", 
                                                    rel_widths = c(0.475, 0.475),
                                                    labels = LETTERS[7:8]), 
                                 markers.plt, rel_heights = c(0.525, 0.475),
                                 nrow = 2, align = "hv", axis = "tblr", labels = c("", LETTERS[9]))
pdf(file = file.path(res.dir[1], "figure_5.pdf"), width = 11.5, height = 17)
cowplot::plot_grid(
  panels_A_C, panels_D_F, panels_G_I, 
  ncol = 1, align = "hv", axis = "tblr", 
  rel_heights = c(0.235, 0.54, 0.225)
)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" 
saveRDS(object = ref, file = file.path(res.dir[3], "ref.rds"))
saveRDS(object = query, file = file.path(res.dir[3], "query.rds"))
saveRDS(object = query.map, file = file.path(res.dir[3], "query_map.rds"))
saveRDS(object = ref.query, file = file.path(res.dir[3], "ref_query.rds"))

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
# [5] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
# [9] MatrixGenerics_1.10.0       matrixStats_1.1.0           Coralysis_1.0.0             ggplot2_3.5.1              
# [13] dplyr_1.1.1                
# 
# loaded via a namespace (and not attached):
#   [1] snow_0.4-4                readxl_1.4.2              parallelDist_0.2.6        plyr_1.8.8                igraph_1.4.1             
# [6] splines_4.2.1             BiocParallel_1.32.6       listenv_0.9.0             scater_1.26.1             digest_0.6.31            
# [11] foreach_1.5.2             viridis_0.6.2             fansi_1.0.4               magrittr_2.0.3            ScaledMatrix_1.6.0       
# [16] cluster_2.1.3             limma_3.54.2              fastcluster_1.2.3         recipes_1.1.0             globals_0.16.2           
# [21] gower_1.0.1               RcppParallel_5.1.7        askpass_1.1               hardhat_1.4.0             colorspace_2.1-0         
# [26] ggrepel_0.9.3             xfun_0.37                 jsonlite_1.8.4            RCurl_1.98-1.10           Exact_3.2                
# [31] survival_3.3-1            iterators_1.0.14          glue_1.6.2                gtable_0.3.3              ipred_0.9-15             
# [36] zlibbioc_1.44.0           XVector_0.38.0            DelayedArray_0.24.0       BiocSingular_1.14.0       future.apply_1.10.0      
# [41] SparseM_1.81              scales_1.3.0              pheatmap_1.0.12           mvtnorm_1.1-3             rngtools_1.5.2           
# [46] edgeR_3.40.2              Rcpp_1.0.10               viridisLite_0.4.1         gridtext_0.1.5            aricode_1.0.2            
# [51] reticulate_1.34.0         dqrng_0.3.0               rsvd_1.0.5                proxy_0.4-27              umap_0.2.10.0            
# [56] lava_1.8.0                prodlim_2024.06.25        httr_1.4.5                metapod_1.6.0             RColorBrewer_1.1-3       
# [61] modeltools_0.2-23         farver_2.1.1              pkgconfig_2.0.3           scuttle_1.8.4             nnet_7.3-17              
# [66] locfit_1.5-9.7            utf8_1.2.3                caret_6.0-94              labeling_0.4.2            tidyselect_1.2.0         
# [71] rlang_1.1.0               reshape2_1.4.4            munsell_0.5.0             cellranger_1.1.0          tools_4.2.1              
# [76] cli_3.6.1                 generics_0.1.3            stringr_1.5.0             ModelMetrics_1.2.2.2      RANN_2.6.1               
# [81] purrr_1.0.1               dendextend_1.16.0         rootSolve_1.8.2.3         doRNG_1.8.6               future_1.32.0            
# [86] nlme_3.1-157              sparseMatrixStats_1.10.0  ggrastr_1.0.1             scran_1.26.2              xml2_1.3.3               
# [91] flexclust_1.4-1           compiler_4.2.1            rstudioapi_0.14           png_0.1-8                 beeswarm_0.4.0           
# [96] e1071_1.7-13              tibble_3.2.1              statmod_1.5.0             DescTools_0.99.48         stringi_1.7.12           
# [101] RSpectra_0.16-1           lattice_0.20-45           bluster_1.8.0             Matrix_1.6-5              commonmark_1.9.0         
# [106] markdown_1.5              vctrs_0.6.1               pillar_1.9.0              lifecycle_1.0.3           LiblineaR_2.10-22        
# [111] BiocNeighbors_1.16.0      data.table_1.14.8         cowplot_1.1.1             bitops_1.0-7              irlba_2.3.5.1            
# [116] lmom_2.9                  R6_2.5.1                  gridExtra_2.3             vipor_0.4.5               parallelly_1.35.0        
# [121] gld_2.6.6                 codetools_0.2-18          boot_1.3-28               MASS_7.3-57               openssl_2.0.6            
# [126] withr_2.5.0               GenomeInfoDbData_1.2.9    doSNOW_1.0.20             expm_0.999-7              parallel_4.2.1           
# [131] ggtext_0.1.2              ILoReg_0.99.6             grid_4.2.1                rpart_4.1.16              beachmat_2.14.0          
# [136] timeDate_4032.109         class_7.3-20              DelayedMatrixStats_1.20.0 Cairo_1.6-0               Rtsne_0.16               
# [141] pROC_1.18.4               lubridate_1.8.0           ggbeeswarm_0.7.1
#
#------------------------------------------------------------------------------#