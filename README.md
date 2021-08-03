
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scBSC (single cell Bivariate Spatial Correlation)

<!-- badges: start -->
<!-- badges: end -->

scBSC helps users calculate bivariate spatial correlation statistics as
from [Lee S., J Geography Syst
(2001)](https://link.springer.com/article/10.1007/s101090100064).

## Installation

You can install the released version of scBSC from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("scBSC")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("soymintc/scBSC")
```

## Example

The following calculates spatial correlation between expression of two
genes, Penk and Cck, from a Seurat-format mouse brain spatial scRNA-seq
data from Visium 10x.

#### Using pre-made SeuratData

``` r
# Import libraries
library(Seurat)
library(SeuratData)
library(scBSC)

# Import Seurat-format mouse brain data
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# Calculate bivariate spatial correlation statistics
conn_mat <- make_conn_mat(brain)
bsc <- calc_bsc("Penk", "Cck", brain, conn_mat)
print_ <- sprintf("%.3f %.3f %.3f %.3f %.3f",
  bsc$L_XX, # How well is Penk expression spatially clustered
  bsc$L_YY, # How well is Cck expression spatially clustered
  bsc$r_sm, # Correlation between spatially smoothened expression of two genes
  bsc$r, # Correlation between expression of two genes
  bsc$L_XY) # Bivariate spatial correlation (clustering effect & expression correlation)
print(print_) # result: 0.869 0.865 -0.809 -0.729 -0.702
```

#### Using Visium 10x data, imported by Seurat Load10X\_Spatial

Of course scBSC also supports your own spatial scRNA-seq, as long as the
data is in Seurat spatial object format. Link to raw data: [Visium 10x
Mouse Brain Serial
Section](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Anterior)

``` r
# Import libraries
library(Seurat)
library(SeuratData)
library(scBSC)

# Import Visium 10x data using Seurat
setwd('~/path/to/mouse_brain_data')
# Use Load10X_Spatial to load data from the following file structure.
# Data used was Visium 10x Mouse Brain Serial Section 1 (Anterior-Sagittal)
# data directory (e.g. "anterior1" as written below)
# ├── filtered_feature_bc_matrix
# │   ├── barcodes.tsv.gz
# │   ├── features.tsv.gz
# │   └── matrix.mtx.gz
# ├── filtered_feature_bc_matrix.h5
# └── spatial
#     ├── aligned_fiducials.jpg
#     ├── detected_tissue_image.jpg
#     ├── scalefactors_json.json # [req] high -> low-res image scaling info
#     ├── tissue_hires_image.png
#     ├── tissue_lowres_image.png # [req]
#     └── tissue_positions_list.csv # [req] the coordinates for each spot
brain <- Load10X_Spatial(data.dir="./anterior1",
                         filename="filtered_feature_bc_matrix.h5",
                         assay='Spatial',
                         slice="anterior1",
                         filter.matrix = TRUE)
brain@meta.data$orig.ident <- "anterior1" # give cells identity; default: SeuratProject
Idents(object = brain) <- "anterior1" # label experiment
Project(object = brain) <- "brain" # label project

# Remove some unwanted genes for your downstream experiments
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")
brain = brain[, brain$nFeature_Spatial > 500 &
                brain$percent_mito < 25 &
                brain$percent_hb < 20]
brain <- brain[!grepl("Bc1", rownames(brain)), ] # filter Bl1
brain <- brain[!grepl("^mt-", rownames(brain)), ] # filter Mitocondrial
brain <- brain[!grepl("^Hb.*-", rownames(brain)), ] # filter Hemoglobin gene (optional if Hb genes needed)

# Normalize, Find HVGs, Scale data
brain <- SCTransform(brain, assay = "Spatial", verbose = TRUE, method = "poisson")

# Calculate bivariate spatial correlation statistics
conn_mat <- make_conn_mat(brain)
bsc <- calc_bsc("Penk", "Cck", brain, conn_mat)
print_ <- sprintf("%.3f %.3f %.3f %.3f %.3f",
  bsc$L_XX, # How well is Penk expression spatially clustered
  bsc$L_YY, # How well is Cck expression spatially clustered
  bsc$r_sm, # Correlation between spatially smoothened expression of two genes
  bsc$r, # Correlation between expression of two genes
  bsc$L_XY) # Bivariate spatial correlation (clustering effect & expression correlation)
print(print_) # result: 0.869 0.865 -0.809 -0.729 -0.702
```

Plotting the relative gene expression between the two genes will show
the following result. Below includes a highly negative, highly positive,
and close-to-zero Lxy cases.

``` r
# Highly positive spatial correlation example with Lxy 0.869, Lyy 0.865, rsm -0.809, r -0.729, Lxy -0.702
SpatialFeaturePlot(brain, features = c("Penk", "Cck"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
```

![Expression between two
genes](figures/Lxy_negative.png?raw=true "Highly negative spatial correlation")

``` r
# Highly negative spatial correlation with Lxy 0.869, Lyy 0.864, rsm 0.936, r 0.855, Lxy 0.813
SpatialFeaturePlot(brain, features = c("Gpr88", "Ppp1r1b"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
```

![Expression between two
genes](figures/Lxy_positive.png?raw=true "Highly positive spatial correlation")

``` r
# Close-to-zero spatial correlation with Lxy 0.461, Lyy 0.569, rsm 0.075, r 0.059, Lxy 0.038
SpatialFeaturePlot(brain, features = c("Rpl34", "Ptn"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
```

![Expression between two
genes](figures/Lxy_neutral.png?raw=true "Close-to-zero spatial correlation")
