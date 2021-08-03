#' Returns a list of bivariate spatial correlation statistics
#' according to Lee S. J Geograph Syst (2001) 3:369-385.
#'
#' @title Calculate Bivariate Spatial Correlation
#' @name calc_bsc
#' @description This function converts input temperatures in Fahrenheit to Celsius.
#' @param feature1 a gene symbol as in Seurat object
#' @param feature2 another gene symbol as in Seurat object
#' @param sdata Seurat formatted spatial 10X scRNA expression data
#' @param conn_mat connectivity matrix (conn_mat) from make_conn_mat
#' @param assay default to 'SCT' but can be 'Spatial' (raw) or 'SCT'
#' @return A list containing $r (Pearson's r),
#' $r_sm (spatially smoothened r),
#' $L_XX (clustering factor of feature1),
#' $L_YY (clustering factor of feature2),
#' $L_XY (L_X,Y bivariate spatial correlation from)
#' @import Matrix
#' @import dplyr
#' @export
#' @examples # Refer to https://github.com/soymintc/scBSC
calc_bsc = function(feature1, feature2, sdata, conn_mat, assay='SCT') {
  # Set X, Y variables from two features
  # Seurat - sdata@assays$SCT@data: row=feature, col=cell matrix
  # feature1 = 'Gpr88' # debug
  # feature2 = 'Penk' # debug
  mrna = sdata@assays[[assay]]@data # default to assay='SCT'
  X_values = mrna[feature1, conn_mat$barcodes_in_tissue] # 1:X
  Y_values = mrna[feature2, conn_mat$barcodes_in_tissue] # 2:Y
  X = data.frame(value=X_values)
  Xmean = mean(X$value)
  Y = data.frame(value=Y_values)
  Ymean = mean(Y$value)

  # Smoothened values
  X['smooth'] = conn_mat$W %*% X[,'value'] # Xsm = W * X
  Y['smooth'] = conn_mat$W %*% Y[,'value'] # Ysm = W * Y
  Xmean_sm = mean(X$smooth) # muX
  Ymean_sm = mean(Y$smooth) # muY

  # Calculate Peason's r(X,Y), r(smooth), L_XX, L_YY, L_XY as in Lee S (2001)
  r = (sum((X$value - Xmean) * (Y$value - Ymean))
       / (sqrt(sum((X$value - Xmean)^2)) * sqrt(sum((Y$value - Ymean)^2))))
  r_sm = (sum((X$smooth - Xmean_sm) * (Y$smooth - Ymean_sm))
          / (sqrt(sum((X$smooth - Xmean_sm)^2)) * sqrt(sum((Y$smooth - Ymean_sm)^2))))
  L_XX = (sum((X$smooth - Xmean)^2) / sum((X$value - Xmean)^2))
  L_YY = (sum((Y$smooth - Ymean)^2) / sum((Y$value - Ymean)^2))
  L_XY = sqrt(L_XX) * sqrt(L_YY) * r_sm

  # Group into Bivariate Spatial Correlation values list
  bsc = list('r' = r, 'r_sm' = r_sm,
             'L_XX' = L_XX, 'L_YY' = L_YY, 'L_XY' = L_XY)

  return(bsc)
}
