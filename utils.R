#' Create a panel of plots, each showing deconvolved fractions of a particular population.
#' 
#' @param rctd An RCTD object, from the spacexr package.
#' @return A data frame whose rows are spots, whose columns are deconvolved populations, and whose 
#'                    entries are the (predicted) fraction of a population in a given spot.
#'                    Also adds columns x and y, holding spatial coordinates of spots.
format.rctd.output_ <- function(rctd) {
  barcodes <- colnames(rctd@spatialRNA@counts)
  weights <- rctd@results$weights
  norm_weights <- normalize_weights(weights)
  df <- as.data.frame(norm_weights)
  df$x <- rctd@spatialRNA@coords$x
  df$y <- rctd@spatialRNA@coords$y
  df
}