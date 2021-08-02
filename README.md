
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scBSC

<!-- badges: start -->
<!-- badges: end -->

scBSC help users calculate bivariate spatial correlation statistics as
from Lee S. J Geography Syst (2001).

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
bsc <- calc_bsc("Vxn", "Dkk3", brain, conn_mat)
print_line = sprintf("%.3f %.3f %.3f %.3f %.3f",
  bsc$L_XX, # How well is Vxn expression spatially clustered
  bsc$L_YY, # How well is Dkk3 expression spatially clustered
  bsc$r_sm, # Correlation between spatially smoothened expression of two genes
  bsc$r, # Correlation between expression of two genes
  bsc$L_XY) # Bivariate spatial correlation (clustering effect + expression correlation)
print(print_line) # result: 0.802 0.696 0.857 0.668 0.64
```

Plotting the relative gene expression between the two genes will show
the following result.

``` r
SpatialFeaturePlot(brain, features = c("Vxn", "Dkk3"), ncol = 2, alpha = c(0.1, 1), max.cutoff = 5)
![Expression between two genes](figures/Vxn_Dkk3.png?raw=true "Expression of Vxn and Dkk3")
```
