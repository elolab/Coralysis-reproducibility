#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make Figure 2. 
# Date: 28/01/2025
# Last update: 29/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages & Set import global variables

# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
source("scripts/helper_functions.R")
#require("cowplot") # v.1.1.1

# Folders to save results 
analysis <- "03_figure_2"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)

## Import data
data.dir <- file.path("data", analysis)
pick.methods <- c("unscaled/full_feature/coralysis_embed.csv", # best for each method selected based on scib-pipeline benchmark ranking
                  "scaled/full_feature/scanorama_embed.csv",
                  "unscaled/full_feature/seuratrpca_full.csv",
                  "unscaled/full_feature/fastmnn_embed.csv",
                  "unscaled/full_feature/seurat_full.csv",
                  "scaled/full_feature/harmony_embed.csv",
                  "unscaled/full_feature/scvi_embed.csv")
data.type <- c("scaled", "unscaled")
feature.type <- c("full_feature", "hvg")

# Integrated embeddings 
embs.dir <- file.path(data.dir, pick.methods)
names(embs.dir) <- gsub("_embed.csv|_full.csv", "", basename(pick.methods))
embs <- lapply(embs.dir, function(x) { 
  read.table(file = x, header = TRUE, sep = ",", colClasses = c("character", "factor", "factor", "numeric", "numeric"))
})

# Unintegrated embeddings
unint.embs <- read.table(file = "data/03_figure_2/unscaled/full_feature/unintegrated_full.csv",
                         header = TRUE, sep = ",", colClasses = c("character", "factor", "factor", "numeric", "numeric"))

# Fix batch factor labels
levels(embs$scanorama$stim) <- levels(embs$harmony$stim) <- c("CTRL", "STIM") 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Plot ILoReg2 integrated embeddings
# Define colors for cell types and batches
color.palette <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                   "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                   "#8C564B", "#C49C94", "#E377C2") 
colors <- list(
  "batch" = c("#8C564B", "#FFBB78"), #color.palette[c(1, 3)], 
  "cell_type" = color.palette[1:nlevels(unint.embs$seurat_annotations)]
)
names(colors$batch) <- levels(unint.embs$stim)
names(colors$cell_type) <- levels(unint.embs$seurat_annotations)
labels <- c("Coralysis_full(-)", "Scanorama_full(+)", "Seurat v4 RPCA_full(-)", "fastMNN_full(-)",
            "Seurat v4 CCA_full(-)", "Harmony_full(+)", "scVI_full(-)")
names(labels) <- names(embs)

# Plot all methods
# unintegrated 
plts <- list()
plts[["unint"]] <- cowplot::plot_grid(plotlist = lapply(X = c("batch", "cell_type"), FUN = function(x) {
  if (x=="batch") title <- c("Unintegrated", "full(-)") else title <- NULL
  plot_umap(x = unint.embs, by = x, size = 0.35, title = title,
            use.color = colors[[x]], draw.line = TRUE) 
}), ncol = 1,  rel_heights = c(0.55, 0.45))
for (m in names(embs)) {
  plts[[m]] <- cowplot::plot_grid(plotlist = lapply(X = c("batch", "cell_type"), FUN = function(x) {
    if (x=="batch") title <- unlist(strsplit(labels[m], "_")) else title <- NULL
    draw.line <- ifelse(m %in% c("seuratrpca", "scvi"), FALSE, TRUE) 
    plot_umap(x = embs[[m]], by = x, size = 0.35, title = title, 
              use.color = colors[[x]], draw.line = draw.line)
  }), ncol = 1, rel_heights = c(0.55, 0.45))
}

# Print unintegrated
pdf(file = file.path(res.dir[1], "ifnb_unintegrated_UMAP_plots.pdf"), width = 4, height = 8)
print(
  cowplot::plot_grid(plotlist = lapply(X = c("batch", "cell_type"), FUN = function(x) {
    if (x=="batch") title <- c("Unintegrated", "full(-)") else title <- NULL
    plot_umap(x = unint.embs, by = x, size = 0.5, title = title,
              use.color = colors[[x]], draw.line = TRUE, no_legend = FALSE) +
      theme(legend.text = element_text(size = 12)) + 
      guides(color = guide_legend(title="", nrow = 4, bycol = TRUE, 
                                  override.aes = list(size=5)))
  }), ncol = 1,  rel_heights = c(0.55, 0.45))
)
dev.off()

# Plot altogether
pdf(file = file.path(res.dir[1], "figure_2_supp_figure_9.pdf"), width = 8, height = 7)
print(
  cowplot::plot_grid(plotlist = plts, ncol = 4, align = "vh")
)
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.5.1 dplyr_1.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.4      withr_2.5.0      utf8_1.2.3       grid_4.2.1       R6_2.5.1         lifecycle_1.0.3  gtable_0.3.3     magrittr_2.0.3  
# [9] scales_1.3.0     pillar_1.9.0     rlang_1.1.0      cli_3.6.1        farver_2.1.1     rstudioapi_0.14  generics_0.1.3   vctrs_0.6.1     
# [17] cowplot_1.1.1    labeling_0.4.2   tools_4.2.1      glue_1.6.2       munsell_0.5.0    compiler_4.2.1   pkgconfig_2.0.3  colorspace_2.1-0
# [25] tidyselect_1.2.0 tibble_3.2.1    
#
#------------------------------------------------------------------------------#