
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scBSC (single cell Bivariate Spatial Correlation)

<!-- badges: start -->
<!-- badges: end -->

scBSC help users calculate bivariate spatial correlation statistics as
from [Lee S. J Geography Syst
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
genes, Gpr88 and Penk, from a Seurat-format mouse brain spatial
scRNA-seq data from Visium 10x.

``` r
# Import libraries
library(Seurat)
library(SeuratData)
library(dplyr)
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
  bsc$L_XY) # Bivariate spatial correlation (clustering effect + expression correlation)
print(print_) # result: 0.869 0.865 -0.809 -0.729 -0.702
```

Plotting the relative gene expression between the two genes will show
the following result. Below includes a highly negative, highly positive,
and close-to-zero cases.

``` r
SpatialFeaturePlot(brain, features = c("Penk", "Cck"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
SpatialFeaturePlot(brain, features = c("Gpr88", "Ppp1r1b"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
SpatialFeaturePlot(brain, features = c("Rpl34", "Ptn"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
```

![Expression between two
genes](figures/Lxy_negative.png?raw=true "Highly negative Lxy case")
![Expression between two
genes](figures/Lxy_positive.png?raw=true "Highly positive Lxy case")
![Expression between two
genes](figures/Lxy_neutral.png?raw=true "Close-to-zero Lxy case")
