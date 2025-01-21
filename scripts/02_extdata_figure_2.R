#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make Extended Data Figure 2. 
# Date: 21/01/2025
# Last update: 21/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages & Set import global variables

# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("Coralysis") # v.1.0.0
library("SingleCellExperiment") # v.1.20.1
source("scripts/helper_functions.R") 

# Folders to save results 
analysis <- "02_exdata_figure_2"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import data 

# SCE object
sce <- readRDS(file = file.path("results/01_figure_1/objects", "sce_int.rds"))

# Top 3 coefficients
top.coeff <- read.table(file = "results/01_figure_1/tables/top10_positive_coefficients_by_cluster_divisive_ICP_L4_table.tsv", 
                        header = TRUE, sep = "\t")
pick.genes <- top.coeff %>% 
  mutate("idx" = 1:nrow(.)) %>% 
  group_by(., K, cluster) %>% 
  slice_max(., n=3, order_by = coeff)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 2

## Plotting the top 3 coefficients per cluster for every ICP clustering round in TSNE
by <- log2(metadata(sce)$coralysis$k)
L <- metadata(sce)$coralysis$L
K <- 1:by
coef.exp.tsne <- list()
for (k in K) {
  coef.exp.tsne[[k]] <- list()
  k.round <- as.character(2**k)
  pick.k <- seq(k, L*by, by)
  set.seed(123)
  sce <- RunPCA(object = sce, assay.name = "joint.probability", select.icp.tables = pick.k, dimred.name = "pca")
  set.seed(123)
  sce <- RunTSNE(object = sce, dimred.name = "TSNE", dimred.type = "pca")
  genes.k <- filter(pick.genes, K == k)
  genes.clt.k <- split(x = genes.k, f = genes.k$cluster)
  for (clt in names(genes.clt.k)) {
    coef.exp.tsne[[k]][[clt]] <- CustomGeneScatterPlot(object = sce, genes = genes.clt.k[[clt]]$gene, 
                                                       return.plot = TRUE, point.size = 0.25,  
                                                       dim.reduction.type = "tsne", title = clt, 
                                                       plot.expressing.cells.last = TRUE, nrow = 1)
  }
}
pdf(file.path(res.dir[1], "extended_data_figure_2_A_C.pdf"), width = 9, height = 8.75)
cowplot::plot_grid(cowplot::plot_grid(coef.exp.tsne[[1]]$`1`, NULL, ncol = 2, align = "vh", labels = c("A", NULL)),
                   cowplot::plot_grid(coef.exp.tsne[[2]]$`1`, coef.exp.tsne[[2]]$`2`, ncol = 2, align = "vh", labels = c("B", NULL)),
                   cowplot::plot_grid(coef.exp.tsne[[2]]$`3`, coef.exp.tsne[[2]]$`4`, ncol = 2, align = "vh"), 
                   cowplot::plot_grid(coef.exp.tsne[[3]]$`1`, coef.exp.tsne[[3]]$`2`, ncol = 2, align = "vh", labels = c("C", NULL)), 
                   cowplot::plot_grid(coef.exp.tsne[[3]]$`3`, coef.exp.tsne[[3]]$`4`, ncol = 2, align = "vh"), 
                   cowplot::plot_grid(coef.exp.tsne[[3]]$`5`, coef.exp.tsne[[3]]$`6`, ncol = 2, align = "vh"), 
                   cowplot::plot_grid(coef.exp.tsne[[3]]$`7`, coef.exp.tsne[[3]]$`8`, ncol = 2, align = "vh"), 
                   ncol = 1, align = "vh")
dev.off()
pdf(file.path(res.dir[1], "extended_data_figure_2_D.pdf"), width = 9, height = 8.75)
cowplot::plot_grid(plotlist = coef.exp.tsne[[4]], ncol = 2, align = "vh", labels = "D")
dev.off()
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
#   [1] bitops_1.0-7              RColorBrewer_1.1-3        tools_4.2.1               doRNG_1.8.6               utf8_1.2.3               
# [6] R6_2.5.1                  irlba_2.3.5.1             HDF5Array_1.26.0          colorspace_2.1-0          rhdf5filters_1.10.0      
# [11] withr_2.5.0               tidyselect_1.2.0          compiler_4.2.1            flexclust_1.4-1           cli_3.6.1                
# [16] BiocNeighbors_1.16.0      SparseM_1.81              DelayedArray_0.24.0       labeling_0.4.2            scales_1.3.0             
# [21] stringr_1.5.0             digest_0.6.31             R.utils_2.12.2            aricode_1.0.2             XVector_0.38.0           
# [26] pkgconfig_2.0.3           sparseMatrixStats_1.10.0  limma_3.54.2              GlobalOptions_0.1.2       rlang_1.1.0              
# [31] readxl_1.4.2              rstudioapi_0.14           DelayedMatrixStats_1.20.0 shape_1.4.6               generics_0.1.3           
# [36] farver_2.1.1              BiocParallel_1.32.6       R.oo_1.25.0               RCurl_1.98-1.10           magrittr_2.0.3           
# [41] BiocSingular_1.14.0       modeltools_0.2-23         GenomeInfoDbData_1.2.9    scuttle_1.8.4             Matrix_1.6-5             
# [46] Rhdf5lib_1.20.0           Rcpp_1.0.10               munsell_0.5.0             fansi_1.0.4               R.methodsS3_1.8.2        
# [51] lifecycle_1.0.3           stringi_1.7.12            edgeR_3.40.2              zlibbioc_1.44.0           rhdf5_2.42.0             
# [56] Rtsne_0.16                plyr_1.8.8                grid_4.2.1                LiblineaR_2.10-22         parallel_4.2.1           
# [61] dqrng_0.3.0               crayon_1.5.2              doSNOW_1.0.20             lattice_0.20-45           cowplot_1.1.1            
# [66] beachmat_2.14.0           circlize_0.4.15           locfit_1.5-9.7            metapod_1.6.0             pillar_1.9.0             
# [71] igraph_1.4.1              rngtools_1.5.2            reshape2_1.4.4            codetools_0.2-18          ScaledMatrix_1.6.0       
# [76] glue_1.6.2                scran_1.26.2              vctrs_0.6.1               foreach_1.5.2             cellranger_1.1.0         
# [81] gtable_0.3.3              RANN_2.6.1                DropletUtils_1.18.1       rsvd_1.0.5                RSpectra_0.16-1          
# [86] class_7.3-20              tibble_3.2.1              snow_0.4-4                pheatmap_1.0.12           iterators_1.0.14         
# [91] cluster_2.1.3             bluster_1.8.0             statmod_1.5.0               
#
#------------------------------------------------------------------------------#