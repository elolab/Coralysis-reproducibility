#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to make figure 3 H,M and supplementary figure 11. 
# Date: 22/01/2025
# Last update: 22/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages, define variables & import data

# Import packages
import os 
import numpy as np # v.1.26.4
import scanpy as sc # v.1.10.3
from scib_metrics.benchmark import Benchmarker, BioConservation # v.0.5.1

# Import data
save_dir = "results/05_figure_3"
outfile = "supp_figure_11" # "figure_3_M" # or "figure_3_H" or "supp_figure_11"
infile = "sce_adt_int.h5ad" # "sce_cytof_int.h5ad" # or "sce_adt_int.h5ad"
adata = sc.read_h5ad(os.path.join(save_dir, "objects", infile))
batch_label = "Batch" # "batch" # or "Batch"
cell_label = "celltype.l2" # "cluster_s" # or "celltype.l2"
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Parse data & benchmark

# Parse data 
import time
adata.obsm["unintPCA"] = adata.obsm["unintPCA"].to_numpy(dtype = "float32")
adata.obsm["Unintegrated"] = adata.obsm["unintPCA"]
adata.obsm["Coralysis"] = adata.obsm["intPCA"].to_numpy(dtype = "float32")
biocons = BioConservation(isolated_labels=False)

# Benchmark
start = time.time()
bm = Benchmarker(
    adata,
    batch_key=batch_label,
    label_key=cell_label,
    embedding_obsm_keys=["Unintegrated", "Coralysis"],
    pre_integrated_embedding_obsm_key="unintPCA",
    bio_conservation_metrics=biocons,
    n_jobs=-1,
)
bm.prepare()
bm.benchmark()
end = time.time()
print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec") # ~52 min CyTOF / ~43 min ADT 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Save output
bm.get_results(min_max_scale=False).to_csv(os.path.join(save_dir, "tables", outfile + ".tsv"), sep = "\t")
bm.plot_results_table(min_max_scale=False, save_dir=os.path.join(save_dir, "plots"))
os.rename(os.path.join(save_dir, "plots", "scib_results.svg"), os.path.join(save_dir, "plots", outfile + ".svg"))
#
#------------------------------------------------------------------------------#
