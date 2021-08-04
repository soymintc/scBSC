#' Returns a list of spatial connectivity matrix variables (C, W)
#' according to Lee S. J Geograph Syst (2001) 3:369-385.
#'
#' @title Make connectivity matrix
#' @name make_conn_mat
#' @description This function calculates the connectivity of a spatial cell to other cells.
#' @param sdata Seurat formatted spatial 10X scRNA expression data
#' @param orig.ident A unique orig.ident for creating connectivity matrix.
#' @param verbose Print intermediate value checks as stderr()
#' @return A list containing $barcodes_in_tissue (barcodes of spatial cells with 'tissue=TRUE'),
#' $nbarcodes_in_tissue (length of barcodes_in_tissue),
#' $C (raw connectivity/adjacency matrix),
#' $W (row-sum divided connectivity matrix),
#' $L_estimate_divR (L estimate, divided by Pearson's R (used later))
#' @import Matrix
#' @import dplyr
#' @import methods
#' @export
#' @examples #' # Refer to https://github.com/soymintc/scBSC

make_conn_mat = function(sdata, orig.ident=NULL, verbose=F) { # input Seurat spatial data
  # Seurat - sdata@images$anterior1@coordinates format:
  #                    tissue row col imagerow imagecol
  # AAACAAGTATCTCCCA-1      1  50 102     7474     8500
  # AAACACCAATAACTGC-1      1  59  19     8552     2788
  # AAACAGAGCGACTCCT-1      1  14  94     3163     7950
  # AAACAGCTTTCAGAAG-1      1  43   9     6636     2100
  if (!.hasSlot(sdata, "meta.data")) stop("[ERROR] data class does not contain @meta.data")
  if (is.null(sdata@meta.data$orig.ident)) stop(paste("[ERROR] @meta.data do not contain $orig.ident"))
  if (length(unique(sdata@meta.data$orig.ident)) > 1) stop("[ERROR] $orig.ident is not unique.")
  if (is.null(orig.ident)) orig.ident = unique(sdata@meta.data$orig.ident)[1]
  if (!.hasSlot(sdata@images[[orig.ident]], "coordinates")) stop("[ERROR] data class does not contain @coordinates")
  if (dim(sdata@images[[orig.ident]]@coordinates)[1] == 0) stop("[ERROR] @coordinates count is zero.")
  coordinates = sdata@images[[orig.ident]]@coordinates
  positions_in_tissue = coordinates %>% filter(tissue==1)
  if (verbose) print(sprintf('dim(positions_in_tissue): %s', dim(positions_in_tissue)))
  barcodes_in_tissue = rownames(positions_in_tissue)
  if (verbose) print(sprintf('length(barcodes_in_tissue): %d', length(barcodes_in_tissue)))
  nbarcodes_in_tissue = length(barcodes_in_tissue)
  if (verbose) print(sprintf('nbarcodes_in_tissue: %d', nbarcodes_in_tissue))

  C = Matrix(nrow=nbarcodes_in_tissue, ncol=nbarcodes_in_tissue, data=0, sparse=T) # C: (sparse) connectivity matrix
  if (!exists("C")) stop("[ERROR] inner matrix C does not exist")
  if (dim(C)[1] == 0) stop("[ERROR] row length of C == 0")
  if (verbose) print(sprintf('dim(C): %s', dim(C)))
  rownames(C) = barcodes_in_tissue
  colnames(C) = barcodes_in_tissue

  for (barcode in barcodes_in_tissue) {
    conn_mat = list() # data holder for connectivity matrices
    row_i = positions_in_tissue[barcode, 'row']
    col_i = positions_in_tissue[barcode, 'col'] # sdata$counts[['A']][[barcode]] later
    # Set nearby nodes and check connectivity (only in_tissue)
    neighbors = subset(coordinates,
                       tissue==1 &
                         (
                           ((row==row_i-1) & (col==col_i-1)) |
                             ((row==row_i-1) & (col==col_i+1)) |
                             ((row==row_i) & (col==col_i-2)) |
                             ((row==row_i) & (col==col_i+2)) |
                             ((row==row_i+1) & (col==col_i-1)) |
                             ((row==row_i+1) & (col==col_i+1))
                         ))
    neighbor_barcodes = rownames(neighbors)
    if (length(neighbor_barcodes) > 0) C[barcode, neighbor_barcodes] = 1
  }
  if (verbose) print(sprintf('dim(C): %s', dim(C)))
  if (verbose) print(sprintf('head(summary(C)): %s', head(summary(C))))

  W = C / rowSums(C) # W: weighted connectivity matrix
  if (!exists("W")) stop("[ERROR] inner matrix W does not exist")
  if (dim(W)[1] == 0) stop("[ERROR] row length of W == 0")
  if (verbose) print(sprintf('head(summary(W)): %s', head(summary(W))))
  W[is.na(W)] = 0

  conn_mat[['L_estimate_divR']] = (sum(diag(t(W %*% W)))
                                   / (nbarcodes_in_tissue - 1))

  conn_mat[['barcodes_in_tissue']] = barcodes_in_tissue
  conn_mat[['nbarcodes_in_tissue']] = nbarcodes_in_tissue
  conn_mat[['W']] = W
  conn_mat[['C']] = C

  return(conn_mat)
}
