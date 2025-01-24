#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to to make Figure 3. 
#Integration of single-cell proteomics data sets: 
#(1) antibody-based sequenced data set from Hao et al., 2021, Cell (PMID: 34062119) 
#(2) mass cytometry (CyTOF) based data set from Rahil et al., 2020, J. Clin. Invest. 
#(PMID: 33044226) & Bjornson-Hooper et al., 2022, Frontiers Immunol. (PMID: 35359965).
#Both data sets were imported from the repository https://github.com/single-cell-proteomic/SCPRO-HI ,
#Koca and Sevilgen, 2024, Proteomics (PMID: 38135888).
# Date: 26/09/2024
# Last update: 22/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("SingleCellExperiment") # v.1.20.1
#require("zellkonverter") # v.1.8.0
#require("Seurat") # v.4.3.0
#require("scran") # v.1.26.2
#require("bluster") # v.1.8.0
#require("cowplot") # v.1.1.1
 
# Folders to save results 
analysis <- "05_figure_3"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Process data 

# Download manually data from the GitHub repository from Koca and Sevilgen, 2024, 
#Proteomics (PMID: 38135888): https://github.com/single-cell-proteomic/SCPRO-HI

# Import the proteomic data set: 
#(1) antibody-based sequenced data set from Hao et al., 2021, Cell (PMID: 34062119)
sample.files <- list.files(path = res.dir[3], pattern = "GSE164378", full.names = TRUE)
names(sample.files) <- paste0("batch", 1:2)
sce.list <- lapply(sample.files, zellkonverter::readH5AD)

# Merge data sets
colnames(sce.list$batch1) <- paste("cell", "batch1", colnames(sce.list$batch1), sep = "_") 
colnames(sce.list$batch2) <- paste("cell", "batch2", colnames(sce.list$batch2), sep = "_") 
stopifnot(all(row.names(sce.list$batch1) == row.names(sce.list$batch2)))
sce <- do.call(cbind, sce.list)

# Seurat normalization/transformation
assayNames(sce) <- "counts"
logcounts(sce) <- counts(sce)
seu <- Seurat::as.Seurat(sce)
seu <- Seurat::NormalizeData(seu, normalization.method = 'CLR', margin = 2)
sce <- Seurat::as.SingleCellExperiment(seu)

# Export SCE object processed 
saveRDS(object = sce, file = file.path(res.dir[3], "sce_adt_proc.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 3 A-B

## Dimensional reduction: before integration

# Main variables
batch.label <- "Batch"
cell.label <- "cluster_s"

# PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "logcounts", dimred.name = "unintPCA")

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "unintPCA", dimred.name = "unintUMAP", 
               umap.method = "uwot")

# Plot 
seeds <- c("Batch" = 4297, "cluster_s" = 705)
unint.adt.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "unintUMAP", seed.color = seeds[x], 
               point.size = 0.1, point.stroke = 0, legend.nrow = 3, rasterise = TRUE, 
             rasterise.dpi = 300) + 
    theme(legend.position = "bottom",
          legend.key.size = unit(0.1, 'mm'),
          legend.justification = "left",
          plot.margin =  unit(c(0, 0, 0, 0), "cm"),
          legend.text = element_text(size = 10),
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) +
    guides(color = guide_legend(title="", nrow = 3, bycol = TRUE,
                                override.aes = list(size=4))) +
    ggtitle(ifelse(x == cell.label, "cell type", "batch")) +
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "Unintegrated\n\nUMAP2", "")
    )
})
unint.adt.plts.all <- cowplot::plot_grid(plotlist = c(unint.adt.plts, list(NULL)), 
                                         ncol = 3, align = "vh", rel_widths = c(1, 1, 0.05), 
                                         labels = c("A", "B", ""))
pdf(file = file.path(file.path(res.dir[1], "figure_3_A_B.pdf")), 
    width = 8, height = 4.5)
print(unint.adt.plts.all)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Multi-level integration
set.seed(1024)
sce <- RunParallelDivisiveICP(object = sce, batch.label = batch.label, 
                              build.train.params = list(nhvg = nrow(sce)), 
                              threads = 4) # 1.730412 hours
saveRDS(sce, file.path(res.dir[3], "sce_adt_int.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 3 C-D

# Dimensional reduction: after integration

# PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "joint.probability", dimred.name = "intPCA")

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "intPCA", dimred.name = "intUMAP", 
               umap.method = "uwot")

# Plots
int.adt.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "intUMAP", seed.color = seeds[x], 
             point.size = 0.1, point.stroke = 0, legend.nrow = 3, rasterise = TRUE, 
             rasterise.dpi = 300) + 
    theme(legend.position = "bottom", 
          legend.key.size = unit(0.1, 'mm'), 
          legend.justification = "left", 
          plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
          legend.text = element_text(size = 10), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    guides(color = guide_legend(title="", nrow = 3, byrow = TRUE, 
                                override.aes = list(size=4))) + 
    ggtitle(ifelse(x == cell.label, "cell type", "batch")) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "Integrated with Coralysis\n\nUMAP2", "")
    )
})
int.adt.plts.all <- cowplot::plot_grid(plotlist = c(int.adt.plts, list(NULL)), 
                                       ncol = 3, align = "vh", rel_widths = c(1, 1, 0.05), 
                                       labels = c("C", "D", ""))
pdf(file = file.path(file.path(res.dir[1], "figure_3_C_D.pdf")), 
    width = 8, height = 4.5)
print(int.adt.plts.all)
dev.off()
plot.title <- title <- cowplot::ggdraw() + 
  cowplot::draw_label("PBMC ADT data from Hao et al., 2021", 
                      x = 0, hjust = 0, size = 11.5) + 
  theme(plot.margin = margin(0, 0, 0, 7))
int.adt.plts.all <- cowplot::plot_grid(plot.title, 
                                       cowplot::plot_grid(unint.adt.plts.all, NULL, int.adt.plts.all, 
                                                          nrow = 1, rel_widths = c(1, 0.05, 1)), 
                                       ncol = 1, rel_heights = c(0.1, 1))
pdf(file = file.path(file.path(res.dir[1], "figure_3_A_D.pdf")), 
    width = 12.75, height = 4)
print(int.adt.plts.all)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 3 E-G

## Cell states integrated
# Clustering
set.seed(1024)
sce$clusters <- scran::clusterCells(sce, use.dimred="intPCA",  
                                    BLUSPARAM= bluster::NNGraphParam(cluster.fun="louvain"))

# Plot
int.cell.states.adt.plts <- lapply(setNames(c(cell.label, "celltype.l2", "clusters"), c(cell.label, "celltype.l2", "clusters")), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "intUMAP", 
             seed.color = c(seeds, "celltype.l2" = 121765, "clusters" = 123)[x], 
             point.size = 0.1, point.stroke = 0, legend.nrow = 8, rasterise = TRUE, 
             rasterise.dpi = 300, label = TRUE) + 
    theme(legend.position = "none", 
          legend.key.size = unit(0.1, 'mm'), 
          legend.justification = "left", 
          plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
          legend.text = element_text(size = 10), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    labs(
      title = ifelse(x == cell.label, "cell type<br>(as in Koca & Sevilgen, 2024)", 
                     ifelse(x != "clusters", "cell type annotated at level 2 of granularity<br>(as in Hao _et al._, 2021)", 
                            "Louvain clusters on Coralysis integrated embedding")),
      x = ifelse(x == cell.label, "UMAP1", ""),
      y = ifelse(x == cell.label, "Integrated with Coralysis\n\nUMAP2", "")
    )
})
plot.title <- title <- cowplot::ggdraw() + 
  cowplot::draw_label("Coralysis embedding recapitulates main cell populations and states originally annotated with CITE-seq", 
                      x = 0, hjust = 0, size = 11.5) + 
  theme(plot.margin = margin(0, 0, 0, 7))
int.cell.states.adt.plts.all <- cowplot::plot_grid(plot.title, 
                                                   cowplot::plot_grid(plotlist = c(int.cell.states.adt.plts, list(NULL)), 
                                                                      ncol = 4, align = "vh", rel_widths = c(1, 1, 1, 0.05), 
                                                                      labels = c("E", "F", "G", "")), 
                                                   ncol = 1, rel_heights = c(0.1, 1))
pdf(file = file.path(file.path(res.dir[1], "figure_3_E_G.pdf")), 
    width = 12.75, height = 4.25)
print(int.cell.states.adt.plts.all)
dev.off()
adt.final.figure <- cowplot::plot_grid(int.adt.plts.all, int.cell.states.adt.plts.all, 
                                       ncol = 1, rel_heights = c(0.5, 0.5))
pdf(file = file.path(file.path(res.dir[1], "figure_3_A_G.pdf")), 
    width = 12.75, height = 8)
print(adt.final.figure)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" and ".h5ad" (for "scib-metrics" package)
saveRDS(sce, file.path(res.dir[3], "sce_adt_int.rds"))
h5ad.obj <- SingleCellExperiment(assays = list(logcounts = as(logcounts(sce), "sparseMatrix")), 
                                 colData = colData(sce),
                                 reducedDims = reducedDims(sce)) 
zellkonverter::writeH5AD(sce = h5ad.obj, file = file.path(res.dir[3], "sce_adt_int.h5ad"), 
                         X_name = "logcounts")

# Clean environment variables
rm(list=ls())
gc()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Import the proteomic data set: 
#(2) mass cytometry (CyTOF) based data set from Rahil et al., 2020, J. Clin. Invest. (PMID: 33044226)
analysis <- "05_figure_3"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
file.paths <- list.files(path = res.dir[3], pattern = "^Human", full.names = TRUE)
names(file.paths) <- c("ifgn", "h1n1")
sce.list <- lapply(X = file.paths, function(x) zellkonverter::readH5AD(file = x))

# Filter the datasets for the no. of proteins shared and merged them into one SCE
table(row.names(sce.list$h1n1) %in% row.names(sce.list$ifgn)) # 39 proteins shared between datasets
shared.proteins <- row.names(sce.list$h1n1)[row.names(sce.list$h1n1) %in% row.names(sce.list$ifgn)] 
sce.list$ifgn[["batch"]] <- "ifgn"
sce.list$h1n1[["batch"]] <- "h1n1"
sce.list <- lapply(X = sce.list, function(x) x[shared.proteins,])
colnames(sce.list$ifgn) <- paste("cell", "ifgn", colnames(sce.list$ifgn), sep = "_")
colnames(sce.list$h1n1) <- paste("cell", "h1n1", colnames(sce.list$h1n1), sep = "_")
sce <- do.call(cbind, sce.list)

# Create logcounts 
logcounts(sce) <- assay(sce, "X")
assay(sce, "X") <- NULL

# Export SCE object processed
saveRDS(object = sce, file = file.path(res.dir[3], "sce_cytof_proc.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 3 I-J

# Dimensional reduction: before integration
batch.label <- "batch"
cell.label <- "cluster_s"

# PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "logcounts", dimred.name = "unintPCA", p = 10)

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "unintPCA", dimred.name = "unintUMAP", 
               umap.method = "uwot")

# Plot unintegrated UMAP
seeds <- c("batch" = 1994, "cluster_s" = 56)
unint.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "unintUMAP", seed.color = seeds[x], 
             point.size = 0.1, point.stroke = 0, legend.nrow = 3, rasterise = TRUE, 
             rasterise.dpi = 300) + 
    theme(legend.position = "bottom",  
          legend.key.size = unit(0.1, 'mm'), 
          legend.justification = "left", 
          plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
          legend.text = element_text(size = 10), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    guides(color = guide_legend(title="", nrow = 3, bycol = TRUE, 
                                override.aes = list(size=4))) + 
    ggtitle(ifelse(x == cell.label, "cell type", x)) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "Unintegrated\n\nUMAP2", "")
    )
})
unint.plts.all <- cowplot::plot_grid(plotlist = c(unint.plts, list(NULL)), ncol = 3, align = "vh", rel_widths = c(1, 1, 0.05), 
                                     labels = c("I", "J", ""))
pdf(file = file.path(file.path(res.dir[1], "figure_3_I_J.pdf")), 
    width = 8, height = 4.5)
print(unint.plts.all)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Multi-level integration
set.seed(1024)
sce <- RunParallelDivisiveICP(object = sce, batch.label = batch.label, 
                              scale.by = "feature", 
                              build.train.params = list(nhvg = nrow(sce), p = 7L),
                              threads = 4) # 43.05687 mins
saveRDS(object = sce, file = file.path(res.dir[3], "sce_cytof_int.rds"))
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Figure 3 K-L 

# Dimensional reduction

# PCA
set.seed(123)
sce <- RunPCA(object = sce, assay.name = "joint.probability", dimred.name = "intPCA")

# UMAP
set.seed(123)
sce <- RunUMAP(object = sce, dimred.type = "intPCA", dimred.name = "intUMAP", 
               umap.method = "uwot")

# Plots
int.plts <- lapply(setNames(c(batch.label, cell.label), c(batch.label, cell.label)), function(x) {
  PlotDimRed(object = sce, color.by = x, dimred = "intUMAP", seed.color = seeds[x], 
             point.size = 0.1, point.stroke = 0, legend.nrow = 3, rasterise = TRUE, 
             rasterise.dpi = 300) + 
    theme(legend.position = "bottom", 
          legend.key.size = unit(0.1, 'mm'), 
          legend.justification = "left", 
          plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
          legend.text = element_text(size = 10), 
          plot.title = ggtext::element_markdown(size = 9.5, hjust = 0.5)) + 
    guides(color = guide_legend(title="", nrow = 3, byrow = TRUE, 
                                override.aes = list(size=4))) + 
    ggtitle(ifelse(x == cell.label, "cell type", x)) + 
    labs(
      x = ifelse(x == batch.label, "UMAP1", ""),
      y = ifelse(x == batch.label, "Integrated with Coralysis\n\nUMAP2", "")
    )
})
int.plts.all <- cowplot::plot_grid(plotlist = c(int.plts, list(NULL)), ncol = 3, align = "vh", rel_widths = c(1, 1, 0.05), 
                                   labels = c("K", "L", ""))
pdf(file = file.path(file.path(res.dir[1], "figure_3_K_L.pdf")), 
    width = 8, height = 4.5)
print(int.plts.all)
dev.off()
plot.title <- title <- cowplot::ggdraw() + 
  cowplot::draw_label("Whole blood CyTOF data from Rahil et al., 2020 (h1n1) & Bjornson-Hooper et al., 2022 (ifgn)", 
                      x = 0, hjust = 0, size = 11.5) + 
  theme(plot.margin = margin(0, 0, 0, 7))
pdf(file = file.path(file.path(res.dir[1], "figure_3_I_L.pdf")), 
    width = 12.75, height = 4)
cowplot::plot_grid(plot.title, cowplot::plot_grid(unint.plts.all, NULL, int.plts.all, nrow = 1, rel_widths = c(1, 0.05, 1)), 
                   ncol = 1, rel_heights = c(0.1, 1))
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE object as ".rds" and ".h5ad" (for "scib-metrics" package)
saveRDS(sce, file.path(res.dir[3], "sce_cytof_int.rds"))
h5ad.obj <- SingleCellExperiment(assays = list(logcounts = as(logcounts(sce), "sparseMatrix")), 
                                 colData = colData(sce),
                                 reducedDims = reducedDims(sce)) 
zellkonverter::writeH5AD(sce = h5ad.obj, file = file.path(res.dir[3], "sce_cytof_int.h5ad"), 
                         X_name = "logcounts")

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
#   [1] utf8_1.2.3                spatstat.explore_3.1-0    reticulate_1.34.0         R.utils_2.12.2            tidyselect_1.2.0         
# [6] htmlwidgets_1.6.2         grid_4.2.1                BiocParallel_1.32.6       Rtsne_0.16                DropletUtils_1.18.1      
# [11] zellkonverter_1.8.0       munsell_0.5.0             ScaledMatrix_1.6.0        codetools_0.2-18          ica_1.0-3                
# [16] statmod_1.5.0             scran_1.26.2              future_1.32.0             miniUI_0.1.1.1            withr_2.5.0              
# [21] spatstat.random_3.1-4     colorspace_2.1-0          progressr_0.13.0          filelock_1.0.2            rstudioapi_0.14          
# [26] Seurat_4.3.0              ROCR_1.0-11               tensor_1.5                listenv_0.9.0             labeling_0.4.2           
# [31] GenomeInfoDbData_1.2.9    polyclip_1.10-4           farver_2.1.1              pheatmap_1.0.12           rhdf5_2.42.0             
# [36] basilisk_1.10.2           parallelly_1.35.0         vctrs_0.6.1               generics_0.1.3            xfun_0.37                
# [41] markdown_1.5              R6_2.5.1                  ggbeeswarm_0.7.1          rsvd_1.0.5                locfit_1.5-9.7           
# [46] bitops_1.0-7              rhdf5filters_1.10.0       spatstat.utils_3.1-1      DelayedArray_0.24.0       promises_1.2.0.1         
# [51] scales_1.3.0              beeswarm_0.4.0            gtable_0.3.3              beachmat_2.14.0           Cairo_1.6-0              
# [56] globals_0.16.2            goftest_1.2-3             spam_2.10-0               rlang_1.1.0               GlobalOptions_0.1.2      
# [61] splines_4.2.1             lazyeval_0.2.2            spatstat.geom_3.1-0       abind_1.4-5               reshape2_1.4.4           
# [66] httpuv_1.6.9              gridtext_0.1.5            tools_4.2.1               ellipsis_0.3.2            RColorBrewer_1.1-3       
# [71] ggridges_0.5.4            Rcpp_1.0.10               plyr_1.8.8                sparseMatrixStats_1.10.0  zlibbioc_1.44.0          
# [76] purrr_1.0.1               RCurl_1.98-1.10           basilisk.utils_1.10.0     deldir_1.0-6              pbapply_1.7-0            
# [81] cowplot_1.1.1             zoo_1.8-11                SeuratObject_5.0.2        ggrepel_0.9.3             cluster_2.1.3            
# [86] magrittr_2.0.3            data.table_1.14.8         RSpectra_0.16-1           scattermore_0.8           SparseM_1.81             
# [91] circlize_0.4.15           lmtest_0.9-40             RANN_2.6.1                fitdistrplus_1.1-8        flexclust_1.4-1          
# [96] patchwork_1.2.0           mime_0.12                 xtable_1.8-4              readxl_1.4.2              gridExtra_2.3            
# [101] shape_1.4.6               compiler_4.2.1            tibble_3.2.1              KernSmooth_2.23-20        R.oo_1.25.0              
# [106] htmltools_0.5.4           later_1.3.0               snow_0.4-4                tidyr_1.3.0               ggtext_0.1.2             
# [111] DBI_1.1.3                 MASS_7.3-57               LiblineaR_2.10-22         Matrix_1.6-5              cli_3.6.1                
# [116] R.methodsS3_1.8.2         parallel_4.2.1            metapod_1.6.0             dotCall64_1.1-1           igraph_1.4.1             
# [121] pkgconfig_2.0.3           dir.expiry_1.6.0          sp_1.6-0                  spatstat.sparse_3.0-1     plotly_4.10.1            
# [126] scuttle_1.8.4             xml2_1.3.3                foreach_1.5.2             vipor_0.4.5               dqrng_0.3.0              
# [131] rngtools_1.5.2            XVector_0.38.0            doRNG_1.8.6               stringr_1.5.0             digest_0.6.31            
# [136] sctransform_0.3.5         RcppAnnoy_0.0.20          spatstat.data_3.0-1       cellranger_1.1.0          leiden_0.4.3             
# [141] uwot_0.1.14               edgeR_3.40.2              DelayedMatrixStats_1.20.0 commonmark_1.9.0          shiny_1.7.4              
# [146] modeltools_0.2-23         nlme_3.1-157              lifecycle_1.0.3           aricode_1.0.2             jsonlite_1.8.4           
# [151] Rhdf5lib_1.20.0           BiocNeighbors_1.16.0      viridisLite_0.4.1         limma_3.54.2              fansi_1.0.4              
# [156] pillar_1.9.0              lattice_0.20-45           ggrastr_1.0.1             fastmap_1.1.1             httr_1.4.5               
# [161] survival_3.3-1            glue_1.6.2                png_0.1-8                 iterators_1.0.14          bluster_1.8.0            
# [166] class_7.3-20              stringi_1.7.12            HDF5Array_1.26.0          BiocSingular_1.14.0       doSNOW_1.0.20            
# [171] irlba_2.3.5.1             future.apply_1.10.0            
#
#------------------------------------------------------------------------------#