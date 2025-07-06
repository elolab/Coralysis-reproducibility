#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make Extended Data Figure 3-5 and Supplementary Figures 7-8.
# Date: 30/01/2025
# Last update: 30/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages & Set import global variables

# Packages
library("dplyr") # v.1.1.1
library("ggplot2") # v.3.5.1 
library("cowplot") # v.1.1.1
library("ComplexHeatmap") # v.2.14.0
source("scripts/helper_functions.R") 

# Folders to save results 
analysis <- "04_extdata_figure_3_4_5"
res.dir <- file.path("results", analysis)
res.dir <- file.path(res.dir, c("plots", "tables", "objects"))
for (f in res.dir) if (!dir.exists(f)) dir.create(path = f, recursive = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import data
results.folder <- file.path("data", analysis)
bench.datasets <- c("immune_cell_hum", "immune_cell_hum_mou", "lung_atlas", "pancreas", "simulations_1_1", "simulations_2")
data.sets <- list.files(path = results.folder, 
                        full.names = FALSE, recursive = FALSE, 
                        pattern = paste(bench.datasets, collapse="|"))
data.type <- c("scaled", "unscaled")
feature.type <- c("full_feature", "hvg")

# Integrated embeddings (Coralysis)
embs.dir <- list.files(path = file.path(results.folder, data.sets), 
                       full.names = TRUE, recursive = TRUE, 
                       pattern = "coralysis_embed.csv") # 24 (4 int. tasks for every of 6 data sets)
embs <- lapply(setNames(data.sets, data.sets), 
               function(x) lapply(setNames(data.type, data.type), 
                                  function(y) lapply(setNames(feature.type, feature.type), 
                                                     function(z) {
                                                       read.table(file = grep(pattern = file.path(x, y, z), x = embs.dir, value = TRUE), 
                                                                  header = TRUE, sep = ",", colClasses = c("character", "factor", "factor", "numeric", "numeric"))
                                                     })))

# Remove dummy embedding for Coralysis task that failed: simulations_1_1 w/ input unscaled & hvg
embs$simulations_1_1$unscaled$hvg <- NULL

# Unintegrated embeddings
unint.embs.dir <- list.files(path = file.path(results.folder, data.sets), 
                             full.names = TRUE, recursive = TRUE, 
                             pattern = "unintegrated_full.csv") # 24 (4 int. tasks for every of 6 data sets)
unint.embs <- lapply(setNames(data.sets, data.sets), 
                     function(x) lapply(setNames(data.type[2], data.type[2]), 
                                        function(y) lapply(setNames(feature.type[1], feature.type[1]), 
                                                           function(z) {
                                                             read.table(file = grep(pattern = file.path(x, y, z), x = unint.embs.dir, value = TRUE),
                                                                        header = TRUE, sep = ",", colClasses = c("character", "factor", "factor", "numeric", "numeric"))
                                                           })))

# Import benchmark scib metrics
scib.file <- file.path(results.folder, "metrics.csv")
scib <- read.table(file = scib.file, header = TRUE, sep = ",", row.names = 1)
pick.bench.data <- grepl(pattern = paste(bench.datasets, collapse="|"), x = row.names(scib))
scib <- scib[pick.bench.data,]
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 4-5

# Plot Coralysis integrated embeddings
plot.all <- list()
data.sets.full <- c("Immune (human)", "Immune (human/mouse)", "Lung", "Pancreas", 
                    "Simulation 1", "Simulation 2")
names(data.sets.full) <- data.sets
for (d in data.sets) {
  if (d == "immune_cell_hum_mou") { 
    size <- 0.01
  } else if (d %in% c("Simulation 1", "Simulation 2")) {
    size <- 0.5
  } else {
    size <- 0.2
  }
  if (d != "simulations_1_1") {
    plot.all[[d]] <- plot_umap_dataset(emb.list = list(embs[[d]]$unscaled$full_feature, 
                                                       embs[[d]]$scaled$full_feature, 
                                                       embs[[d]]$unscaled$hvg, 
                                                       embs[[d]]$scaled$hvg), 
                                       dataset = data.sets.full[d], size = size)
  } else {
    plot.all[[d]] <- plot_umap_dataset(emb.list = list(embs[[d]]$unscaled$full_feature, 
                                                       embs[[d]]$scaled$full_feature, 
                                                       embs[[d]]$scaled$hvg, # dummy plot to print just the batch legend - remove it manually
                                                       embs[[d]]$scaled$hvg), 
                                       dataset = data.sets.full[d], size = size)
  }

  pdf(file = file.path(res.dir[1], paste0(d, "_coralysis_integrated_UMAPs_diff_inputs.pdf")), width = 5, height = 4)
  print(plot.all[[d]])
  dev.off()
}

## Plot altogether 
# real data sets
pdf(file = file.path(res.dir[1], "extdata_figure_4.pdf"), width = 10, height = 8)
print(
  plot_grid(plot.all$pancreas, plot.all$lung_atlas, plot.all$immune_cell_hum, plot.all$immune_cell_hum_mou, 
            ncol = 2, align = "vh")
) 
dev.off()
# simulations
pdf(file = file.path(res.dir[1], "extdata_figure_5.pdf"), width = 10, height = 4)
print(
  plot_grid(plot.all$simulations_1_1, plot.all$simulations_2, 
            ncol = 2, align = "vh")
) 
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Supplementary Figure 7-8

## Plot unintegrated embeddings
plot.unint.all <- list()
for (d in data.sets) {
  plot.unint.all[[d]] <- plot_grid(plotlist = lapply(c("batch", "cell"), function(x) {
    if (x=="batch") {
      plot_umap(x = unint.embs[[d]]$unscaled$full_feature, by = x, size = 0.2, legend.nrow = 10) + 
        theme(plot.margin = unit(c(0.5, 0, 3, 0), "cm"), 
              legend.key.spacing.x = unit(0.5, "mm"), 
              legend.key.spacing.y = unit(0.5, "mm"), 
              legend.text = element_text(margin = margin(l = unit(0.5, "mm"))),
              legend.justification = "left", 
              legend.position = c(0, -0.175)) + 
        labs(subtitle = "full(-)")
    } else {
      plot_umap(x = unint.embs[[d]]$unscaled$full_feature, by = x, size = 0.2, legend.nrow = 10) + 
        theme(plot.margin = unit(c(0.5, 0, 3, 0), "cm"), 
              legend.key.spacing.x = unit(0.5, "mm"), 
              legend.key.spacing.y = unit(0.5, "mm"), 
              legend.text = element_text(margin = margin(l = unit(0.5, "mm"))),
              legend.justification = "left", 
              legend.position = c(0, -0.175)) + 
        annotate(geom = "segment", y = -Inf, yend = Inf, color = "gray", 
                 x = Inf, xend = Inf, linewidth = 0.5) 
    }
  }), 
  ncol=2, align = "vh", label_x = 0, label_y = 0.96,
  hjust = 0, vjust = 0, labels = c(data.sets.full[d], "")
  )
  pdf(file = file.path(res.dir[1], paste0(d, "_coralysis_unintegrated_UMAPs.pdf")), width = 5.5, height = 4.25)
  print(plot.unint.all[[d]])
  dev.off()
}

## Plot altogether
# real data sets
pdf(file = file.path(res.dir[1], "supp_figure_7.pdf"), width = 10, height = 8)
print(
  plot_grid(plot.unint.all$pancreas, plot.unint.all$lung_atlas, 
            plot.unint.all$immune_cell_hum, plot.unint.all$immune_cell_hum_mou, 
            ncol = 2, align = "vh")
) 
dev.off()

# simulation data sets
pdf(file = file.path(res.dir[1], "supp_figure_8.pdf"), width = 10, height = 4)
print(
  plot_grid(plot.unint.all$simulations_1_1, plot.unint.all$simulations_2, 
            ncol = 2, align = "vh")
) 
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 3 A, C, D

## Plot scatter plot of batch correction vs bio-conservation
# Parse data
ids <- row.names(scib)
scib.data <- scib %>% 
  mutate("Dataset" = unlist(lapply(ids, function(x) unlist(strsplit(x, split = "/"))[2])), 
         "Method" = unlist(lapply(ids, function(x) unlist(strsplit(x, split = "/"))[6])),
         "Scaling" = unlist(lapply(ids, function(x) unlist(strsplit(x, split = "/"))[4])), 
         "Features" = unlist(lapply(ids, function(x) unlist(strsplit(x, split = "/"))[5]))) %>% 
  tidyr::separate(., Method, c("Method", "Output")) %>% 
  select(Dataset:Features, `NMI_cluster.label`:trajectory) 
# Set as NAs methods that failed to run
# failed.tasks <- c(
#   "/lung_atlas/metrics/scaled/full_feature/seurat_full", 
#   "/immune_cell_hum_mou/metrics/scaled/full_feature/seuratrpca_full",
#   "/immune_cell_hum_mou/metrics/scaled/hvg/seuratrpca_full",
#   "/immune_cell_hum_mou/metrics/unscaled/full_feature/seuratrpca_full",
#   "/immune_cell_hum_mou/metrics/unscaled/hvg/seuratrpca_full",
#   "/immune_cell_hum_mou/metrics/scaled/full_feature/seurat_full"
# )
# Remove failed tasks 
failed.tasks <- c("/immune_cell_hum_mou/metrics/scaled/full_feature/seurat_full",
                  "/immune_cell_hum_mou/metrics/scaled/full_feature/seuratrpca_full",
                  "/immune_cell_hum_mou/metrics/scaled/hvg/seuratrpca_full",
                  "/lung_atlas/metrics/scaled/full_feature/seurat_full",
                  "/simulations_2/metrics/unscaled/hvg/seuratrpca_full",
                  "/simulations_2/metrics/unscaled/full_feature/seuratrpca_full",
                  "/immune_cell_hum_mou/metrics/unscaled/full_feature/seuratrpca_full",
                  "/immune_cell_hum_mou/metrics/unscaled/hvg/seuratrpca_full",
                  "/simulations_2/metrics/scaled/hvg/seuratrpca_full",
                  "/simulations_2/metrics/scaled/full_feature/seuratrpca_full",
                  "/simulations_1_1/metrics/unscaled/hvg/coralysis_embed"
)
stopifnot(all(failed.tasks %in% row.names(scib)))
scib.data[failed.tasks, !colnames(scib.data) %in% c("Dataset", "Method", "Scaling", "Features")] <- NA

# Export table without tasks that failed
scib.wo.failed.tasks <- filter(scib, !row.names(scib) %in% failed.tasks) 
write.table(x = cbind(" " = row.names(scib.wo.failed.tasks), scib.wo.failed.tasks), 
            file = file.path(res.dir[2], "metrics_without_failed_tasks.csv"), 
            quote = FALSE, sep = ",", row.names = FALSE)

## Plot scatter
metrics <- list(
  "Batch correction" = c("PCR_batch" = "PCR batch", "ASW_label.batch" = "Batch ASW", 
                         "iLISI" = "Graph iLISI", "graph_conn" = "Graph connectivity", 
                         "kBET" = "kBET"), 
  "Bio conservation" = c("NMI_cluster.label" = "NMI cluster/label", "ARI_cluster.label" = "ARI cluster/label", 
                         "ASW_label" = "Cell type ASW", "isolated_label_F1" = "Isolated label F1", 
                         "isolated_label_silhouette" = "Isolated label silhouette",
                         "cLISI" = "Graph cLISI", "cell_cycle_conservation" = "Cell cycle conservation", 
                         "hvg_overlap" = "HVG conservation", "trajectory" = "Trajectory conservation")
)
color.palettes <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
color.palette <- unlist(mapply(RColorBrewer::brewer.pal, color.palettes$maxcolors, rownames(color.palettes)))
set.seed(123)
use.color <- sample(color.palette, 35)
scib.data.scatter <- scib.data %>% 
  rowwise() %>%
  mutate("Batch_correction" = mean(c(PCR_batch, `ASW_label.batch`, iLISI, graph_conn, kBET), na.rm = TRUE), 
         "Bio_conservation" = mean(c(`NMI_cluster.label`, `ARI_cluster.label`, ASW_label, isolated_label_F1, 
                                     isolated_label_silhouette, cLISI, cell_cycle_conservation, hvg_overlap, 
                                     trajectory), na.rm = TRUE), 
         "Overall" = weighted.mean(c(Batch_correction, Bio_conservation), c(0.4, 0.5))) %>% 
  filter(!is.na(Overall))
benchmark.dataset.plt.data <- scib.data.scatter %>% 
  mutate("Task" = paste(Method, Output, Features, Scaling, sep="_"), 
         "Type" = ifelse(Dataset %in% paste0("simulations_", c("1_1", "2")), "Simulations", "Real"), 
         "Dataset" = factor(Dataset, levels = c("pancreas", "lung_atlas", "immune_cell_hum", 
                                                "immune_cell_hum_mou", "simulations_1_1", "simulations_2")), 
         "Output-Input" = paste0(ifelse(Output=="embed", "Embedding", "Genes"), " - ", paste(gsub("_feature", "", Features), ifelse(Scaling=="scaled", "(+)", "(-)"), sep="")), 
         "Method" = ifelse(Method == "coralysis", "Coralysis", ifelse(Method == "seurat", "Seurat v4 CCA", 
                                                                  ifelse(Method == "seuratrpca", "Seurat v4 RPCA", 
                                                                         ifelse(Method == "scvi", "scVI", ifelse(Method == "fastmnn", "fastMNN", 
                                                                                                                 stringr::str_to_title(Method))))))) 
benchmark.dataset.plt.data$Method <- factor(benchmark.dataset.plt.data$Method, 
                                            levels = c("Coralysis", "scVI", "Scanorama", "fastMNN", "Harmony", 
                                                       "Seurat v4 CCA", "Seurat v4 RPCA", "Unintegrated"))
benchmark.dataset.plt <- ggplot(data = benchmark.dataset.plt.data, 
                                mapping = aes(x = Batch_correction, y = Bio_conservation, 
                                              color = `Output-Input`, shape = Method, 
                                              fill = `Output-Input`)) + 
  geom_point(#fill = ifelse(benchmark.dataset.plt.data$Method=="Coralysis", "gray", "white"), 
    stroke = ifelse(benchmark.dataset.plt.data$Method=="Coralysis", 2, 1), 
    size = 2) +  
  scale_shape_manual(values = c(24, 0:1, 3:8)) +
  facet_wrap(~Dataset, ncol = 4, labeller = labeller(Dataset = data.sets.full)) + 
  theme_minimal() +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Reds")[c(3, 5, 7, 9)], 
                                RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)])) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Reds")[c(3, 5, 7, 9)], 
                               RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)])) +
  scale_x_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  theme(#legend.position = "bottom", 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14), 
    legend.text = element_text(size = 12, margin = margin(l = unit(2, "mm"), r = unit(5, "mm"))),
    legend.key.spacing.x = unit(0.5, "mm"), 
    legend.key.spacing.y = unit(2, "mm"), 
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.key.size = unit(1, 'mm'), 
    legend.justification = "left", 
    legend.title.position = "top", 
    legend.title = element_text(size = 14),
    #legend.title.align = 0.5, 
    strip.text.x = element_text(size = 14), 
    legend.position = c(0.55, 0.2), 
    legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 10, bycol = TRUE, 
                              override.aes = list(size=3.5)), 
         shape = guide_legend(nrow = 10, bycol = TRUE, 
                              override.aes = list(size=3.5))) +
  xlab("Batch correction") + 
  ylab("Bio-conservation")
pdf(file = file.path(res.dir[1], "extdata_figure_3_A.pdf"), width = 12, height = 6.5)
print(benchmark.dataset.plt)
dev.off()

# Average metrics across data sets
benchmark.plt.data <- benchmark.dataset.plt.data %>% 
  group_by(Method, `Output-Input`) %>% 
  summarise("Bio_conservation_mean" = mean(Bio_conservation), 
            "Bio_conservation_sd" = sd(Bio_conservation), 
            "Batch_correction_mean" = mean(Batch_correction), 
            "Batch_correction_sd" = sd(Batch_correction), 
            .groups = "drop") 

benchmark.plt <- ggplot(data = benchmark.plt.data, 
                        mapping = aes(x = Batch_correction_mean, y = Bio_conservation_mean, 
                                      shape = Method, color = `Output-Input`, fill = `Output-Input`)) +
  geom_point(#fill = ifelse(benchmark.plt.data$Method=="Coralysis", "gray", "white"), 
    size = 5, stroke = ifelse(benchmark.plt.data$Method=="Coralysis", 2, 1)) +  
  # geom_linerange(aes(ymin=Bio_conservation_mean-Bio_conservation_sd, ymax=Bio_conservation_mean+Bio_conservation_sd, 
  #                    xmin=Batch_correction_mean-Batch_correction_sd, xmax=Batch_correction_mean+Batch_correction_sd)) + 
  scale_shape_manual(values = c(24, 0:1, 3:8)) +
  theme_minimal() +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Reds")[c(3, 5, 7, 9)], 
                                RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)])) + 
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Reds")[c(3, 5, 7, 9)], 
                               RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)])) +
  scale_x_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  theme(#legend.position = "bottom", 
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16), 
    legend.text = element_text(size = 12, margin = margin(l = unit(2, "mm"), r = unit(5, "mm"))),
    legend.key.spacing.x = unit(0.5, "mm"), 
    legend.key.spacing.y = unit(2, "mm"), 
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.key.size = unit(1, 'mm'), 
    legend.justification = "left", 
    legend.title.position = "top", 
    legend.title = element_text(size = 14),
    #legend.title.align = 0.5, 
    strip.text.x = element_text(size = 14), 
    legend.position = c(0.05, 0.25),
    legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 10, bycol = TRUE, 
                              override.aes = list(size=3.5)), 
         shape = guide_legend(nrow = 10, bycol = TRUE, 
                              override.aes = list(size=3.5))) +
  xlab("Batch correction") + 
  ylab("Bio-conservation")
pdf(file = file.path(res.dir[1], "extdata_figure_3_C.pdf"), width = 7, height = 6)
print(benchmark.plt)
dev.off()

## Boxplot of variance given input-outputs per dataset per method
# Estimating the Median Performance per Dataset using the Overall 
#performance score per method per dataset
median.performance.dataset <- benchmark.dataset.plt.data %>% 
  filter(Method != "Unintegrated") %>% 
  group_by(Method, Dataset) %>% 
  mutate("Median_Performance" = median(Overall)) %>% 
  ungroup(.) %>% 
  add_row(Dataset = "simulations_2", Method = "Seurat v4 RPCA")

# Add Median Performance for the 'Method' "Seurat v4 RPCA"
#in the 'Dataset' "simulations_2" as this method failed 
#for all the tasks related with this dataset. Performing 
#the mean of the median performance for the Seurat v4 RPCA
#method without considering the Simulations 2 dataset would
#overestimate its performance, and its inclusion (as 0) would 
#underestimate its performance. Therefore, it was estimated that 
#its performance would be the mean of the median performance 
#from all the other methods for the dataset Simulations 2, 
#yielding a Median Performance of 0.5290003.
median.performance.dataset[nrow(median.performance.dataset), "Median_Performance"] <- median.performance.dataset %>% 
  filter(Dataset == "simulations_2" & Method != "Seurat v4 RPCA") %>% 
  pull(Median_Performance) %>% 
  mean(.) 

# Adding the Mean of the Median Performance per Dataset
mean.median.performance <- median.performance.dataset[,c(1:2,26)] %>% 
  distinct(.) %>%
  group_by(Method) %>% 
  mutate("Mean_Median_Performance" = mean(Median_Performance)) %>%
  ungroup(.) %>%
  select(Method, Mean_Median_Performance) %>% 
  distinct(.) %>%
  arrange(tolower(Method)) %>% 
  mutate("Rank" = dense_rank(desc(Mean_Median_Performance)), 
         "x" = row_number()-0.4, 
         "xend" = row_number()+0.4, 
         "color" = alpha("#117733", 1.1 - (Rank / 10)))

# Add metrics to overall data frame
benchmark.dataset.plt.data.metrics <- left_join(median.performance.dataset, mean.median.performance) 
write.table(x = as.data.frame(benchmark.dataset.plt.data.metrics[,!colnames(benchmark.dataset.plt.data.metrics) %in% c("x", "xend", "color")]), 
            file = file.path(res.dir[2], "metrics_without_failed_tasks_with_additional_metrics.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

# Plot
benchmark.plt.boxplot <- benchmark.dataset.plt.data.metrics %>% 
  ggplot(data = ., mapping = aes(x = Method, y = Overall)) + 
  geom_hline(mapping = aes(yintercept = Mean_Median_Performance), 
             color = "gray") + 
  geom_boxplot(mapping = aes(fill=Dataset)) + 
  geom_segment(data = mean.median.performance, mapping = aes(x = x, xend = xend, 
                             y = Mean_Median_Performance),
               color = mean.median.performance$color,
               linewidth = 1.5) +
  scale_fill_brewer(palette="Spectral", labels = data.sets.full) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0.4,1)) + 
  scale_alpha(guide = "none") + 
  theme(#legend.position = "bottom", 
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
    legend.text = element_text(size = 12, margin = margin(l = unit(2, "mm"), r = unit(5, "mm"))),
    legend.key.spacing.x = unit(0.5, "mm"), 
    legend.key.spacing.y = unit(2, "mm"), 
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.key.size = unit(1, 'mm'), 
    legend.justification = "left", 
    legend.title.position = "top", 
    legend.position = c(0.75, 0.8),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14), 
    title = element_text(size = 16), 
    legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(nrow = 10, bycol = TRUE, override.aes = list(size=3.5))) +
  ylab("Performance of integration:\n((0.6*Bio-conservation) + (0.4*Batch correction))") + 
  ggtitle("Variance in performance for different input/output(s)")
pdf(file = file.path(res.dir[1], "extdata_figure_3_D.pdf"), width = 10, height = 6.5)
print(benchmark.plt.boxplot)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Extended Data Figure 3 B

# Plot Heatmap 
heat.data <- benchmark.dataset.plt.data %>% 
  arrange(Dataset, Method, Output)  
heat.data2 <- select(heat.data, `NMI_cluster.label`:trajectory) 
heat.data2 <- as.data.frame(heat.data2)
row.names(heat.data2) <- paste(heat.data$Dataset, 
                               heat.data$Method, 
                               heat.data$`Output-Input`, 
                               sep = "_")
heat.data.set <- split(x = heat.data2, f = heat.data$Dataset)
scl.row.heat.data.set <- scl.col.heat.data.set <- list()
for (i in names(heat.data.set)) {
  #heat.data.set[[i]][is.na(heat.data.set[[i]])] <- 0 # remove NAs
  scl.col.heat.data.set[[i]] <- scale(heat.data.set[[i]]) # scale by col
  scl.col.heat.data.set[[i]][is.na(scl.col.heat.data.set[[i]])] <- 0 # remove NAs
  scl.row.heat.data.set[[i]] <- t(scale(t(heat.data.set[[i]]))) # scale by row
  scl.row.heat.data.set[[i]][is.na(scl.row.heat.data.set[[i]])] <- 0 # remove NAs
}

# Plot heatmaps
heat.list[["row"]] <- heat.list[["col"]] <- heat.list <- list()
for (j in c("row", "col")) {
  if (j == "row") {
    tmp.data <- scl.row.heat.data.set
    legend.title <- "Row Z-score"
  } else {
    tmp.data <- scl.col.heat.data.set
    legend.title <- "Column Z-score"
  }
  for (i in names(heat.data.set)) {
    if (i == names(heat.data.set)[length(heat.data.set)]) {
      heat.plt <- plot_scib_metrics_heatmap(scaled_data = tmp.data[[i]], 
                                            meta_data = benchmark.dataset.plt.data, 
                                            metrics_list = metrics, 
                                            column_title = data.sets.full[i], 
                                            legend_title = legend.title)
    } else {
      heat.plt <- plot_scib_metrics_heatmap(scaled_data = tmp.data[[i]], 
                                            meta_data = benchmark.dataset.plt.data, 
                                            metrics_list = metrics, 
                                            show_annot_legend = FALSE,
                                            show_heatmap_legend = FALSE, 
                                            column_title = data.sets.full[i], 
                                            legend_title = legend.title)
    }
    heat.list[[j]][[i]] <- heat.plt
  }
}
heat.row.scib.metrics.plts <- cowplot::plot_grid(plotlist = heat.list$row, ncol = 6, rel_widths = c(rep(0.145, 5), 0.27))
heat.col.scib.metrics.plts <- cowplot::plot_grid(plotlist = heat.list$col, ncol = 6, rel_widths = c(rep(0.145, 5), 0.27))
pdf(file = file.path(res.dir[1], "scib_metrics_scaled_within_method_heatmap_plot.pdf"), w=12.5, h=4.75) 
heat.row.scib.metrics.plts
dev.off()
pdf(file = file.path(res.dir[1], "extdata_figure_3_B.pdf"), w=12.5, h=4.75) 
heat.col.scib.metrics.plts
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Assemble figure altogether
pdf(file = file.path(res.dir[1], "extdata_figure_3.pdf"), width = 12.5, height = 16.5) 
cowplot::plot_grid(benchmark.dataset.plt, heat.col.scib.metrics.plts, 
                   cowplot::plot_grid(benchmark.plt, benchmark.plt.boxplot, ncol=2, align = "vh", 
                                      labels = c("C", "D"), rel_widths = c(0.45, 0.55)), 
                   ncol = 1, align = "vh", labels = c("A", "B", ""), rel_heights = c(0.35, 0.3, 0.33))
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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.14.0 cowplot_1.1.1         ggplot2_3.5.1         dplyr_1.1.1          
# 
# loaded via a namespace (and not attached):
#   [1] vipor_0.4.5         ggrastr_1.0.1       pillar_1.9.0        compiler_4.2.1      RColorBrewer_1.1-3  iterators_1.0.14    tools_4.2.1        
# [8] digest_0.6.31       clue_0.3-64         lifecycle_1.0.3     tibble_3.2.1        gtable_0.3.3        pkgconfig_2.0.3     png_0.1-8          
# [15] rlang_1.1.0         foreach_1.5.2       cli_3.6.1           rstudioapi_0.14     parallel_4.2.1      beeswarm_0.4.0      cluster_2.1.3      
# [22] withr_2.5.0         stringr_1.5.0       S4Vectors_0.36.2    IRanges_2.32.0      GlobalOptions_0.1.2 generics_0.1.3      vctrs_0.6.1        
# [29] stats4_4.2.1        tidyselect_1.2.0    glue_1.6.2          R6_2.5.1            GetoptLong_1.0.5    fansi_1.0.4         ggbeeswarm_0.7.1   
# [36] farver_2.1.1        tidyr_1.3.0         purrr_1.0.1         magrittr_2.0.3      BiocGenerics_0.44.0 scales_1.3.0        codetools_0.2-18   
# [43] matrixStats_1.1.0   shape_1.4.6         circlize_0.4.15     colorspace_2.1-0    labeling_0.4.2      utf8_1.2.3          stringi_1.7.12     
# [50] munsell_0.5.0       doParallel_1.0.17   rjson_0.2.21        crayon_1.5.2        Cairo_1.6-0 
#
#------------------------------------------------------------------------------#

