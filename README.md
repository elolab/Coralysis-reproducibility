
<!-- README.md is generated from README.Rmd. Please edit that file -->

# :recycle: Coralysis-reproducibility

Code to reproduce the analyses from the
[Coralysis](https://github.com/elolab/Coralysis) manuscript:

> **António GG Sousa, Johannes Smolander, Sini Junttila, Laura L Elo**
> (2025).  
> *Coralysis enables sensitive identification of imbalanced cell types
> and states in single-cell data via multi-level integration.*
> *bioRxiv*. <https://doi.org/10.1101/2025.02.07.637023>

<br>

## :card_index_dividers: Scripts & Datasets

All the scripts required to reproduce the **Coralysis** manuscript
figures described in the table below are available in this repository
under the folder [`scripts`](./scripts/), with the exception of the code
from the last row, which has its own repository and respective
documentation at
[`elolab/scib-pipeline`](https://github.com/elolab/scib-pipeline).

All the datasets used in the manuscript are publicly available online,
and the links are provided in the table below. A more detailed
description of the datasets and their respective references can be found
at [`data/code_datasets_description.tsv`](./data).

##### 🔗 **Summary of scripts used to reproduce the Coralysis manuscript figures.**

| Scripts                                                           | Description                                                  | Data                                                                                                                                                                                                                  |
|-------------------------------------------------------------------|--------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `01_figure_1.R`                                                   | R script to make Figure 1                                    | [Assay V1](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc6k); [Assay V2](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k)                          |
| `02_extdata_figure_2.R`                                           | R script to make Supplementary Figure S2                     | [Assay V1](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc6k); [Assay V2](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k)                          |
| `03_figure_2.R`                                                   | R script to make Figure 2B & Supplementary Figure S14        | `ifnb` from [SeuratData](https://github.com/satijalab/seurat-data) ([parsed](https://github.com/elolab/scib-pipeline/blob/main/data/create_ifnb_dataset.R))                                                           |
| `04_extdata_figure_3_4_5.R`                                       | R script to make Supplementary Figure S3,S10-S13             | [Figshare](https://doi.org/10.6084/m9.figshare.12420968)                                                                                                                                                              |
| `05_figure_3.R`                                                   | R script to to make Figure 3A-G,I                            | [GitHub](https://github.com/single-cell-proteomic/SCPRO-HI/tree/main/Data)                                                                                                                                            |
| `06_figure_3.py`                                                  | Python script to make Figure 3H,M & Supplementary Figure S16 | [GitHub](https://github.com/single-cell-proteomic/SCPRO-HI/tree/main/Data)                                                                                                                                            |
| `07_figure_4.R`                                                   | R script to make Figure 4                                    | `panc8` from [SeuratData](https://github.com/satijalab/seurat-data)                                                                                                                                                   |
| `08_extdata_figure_6.R`                                           | R script to make Supplementary Figure S17,S18                | [Figshare](https://doi.org/10.6084/m9.figshare.24625302.v1)                                                                                                                                                           |
| `09_extdata_figure_7.R`                                           | R script to make Supplementary Figure S19,S20                | [Figshare](https://doi.org/10.6084/m9.figshare.24625302.v1)                                                                                                                                                           |
| `10_extdata_figure_8.R`                                           | R script to make Supplementary Figure S21                    | `ifnb` from [SeuratData](https://github.com/satijalab/seurat-data)                                                                                                                                                    |
| `11_figure_5.R`                                                   | R script to make Figure 5                                    | `pbmcsca` from [SeuratData](https://github.com/satijalab/seurat-data)                                                                                                                                                 |
| `12_figure_6.R`                                                   | R script to make Figure 6 & Supplementary Figure S22         | [Zenodo](https://zenodo.org/records/6383269/files/cd34_multiome_rna.h5ad); [Figshare](https://figshare.com/ndownloader/files/41674599)                                                                                |
| `helper_functions.R`                                              | R script with custom R functions                             |                                                                                                                                                                                                                       |
| [`elolab/scib-pipeline`](https://github.com/elolab/scib-pipeline) | Snakemake workflow to benchmark Coralysis                    | [Figshare](https://doi.org/10.6084/m9.figshare.12420968); `ifnb` from [SeuratData](https://github.com/satijalab/seurat-data) ([parsed](https://github.com/elolab/scib-pipeline/blob/main/data/create_ifnb_dataset.R)) |

<br>

## :hammer_and_wrench: Reproducibility

All `R` scripts were run with `R` version 4.2.1 under the RStudio Server
environment (v.2022.07.2 Build 576), deployed through the Docker image
[`elolab/sctoolkit`](https://hub.docker.com/r/elolab/sctoolkit).

The analyses can be reproduced using the
[Coralysis](https://github.com/elolab/Coralysis) version corresponding
to the commit `47f1b3415663ee895df188f264ac4d8ad8d24c11`, which can be
installed as follows:

``` r
devtools::install_github("elolab/Coralysis", ref = "47f1b3415663ee895df188f264ac4d8ad8d24c11")
```

The remaining `R` package dependencies and their respective versions can
be found at the end of every `R` script (obtained with the `R` command
`sessionInfo()`). The `R` scripts were run in the order specified by the
prefix in each file name — i.e., `01_figure_1.R`, followed by
`02_extdata_figure_2.R`, and so on.

The `Python` script (`06_figure_3.py`) was run with `Python` version
3.9.16, along with the packages `numpy` (v.1.26.4), `scanpy` (v.1.10.3),
and `scib_metrics` (v.0.5.1).

The [`elolab/scib-pipeline`](https://github.com/elolab/scib-pipeline)
benchmark was run with `Snakemake` (v.7.25.2) in a cluster environment
using `Slurm` (v.23.02.6) with 8 threads and 354 GB of RAM. The
respective `conda` environments can be found under `envs`
([`scib-R4.1.yml`](https://github.com/elolab/scib-pipeline/blob/main/envs/scib-R4.1.yml)
and
[`scib-pipeline-R4.1.yml`](https://github.com/elolab/scib-pipeline/blob/main/envs/scib-pipeline-R4.1.yml)).
The configuration file used for the benchmark is available here:
[`configs/benchmark_coralysis.yaml`](https://github.com/elolab/scib-pipeline/blob/main/configs/benchmark_coralysis.yaml).

After activating the `Snakemake` `conda` environment, the benchmark was
initiated with the following command:

``` bash
snakemake --configfile configs/benchmark_coralysis.yaml --cores 8 
```

<br>

## :classical_building: Funding

> This project has received funding from the European Union’s Horizon
> 2020 research and innovation programme under the Marie
> Skłodowska-Curie grant agreement no.: 955321.

<br>

<br>
