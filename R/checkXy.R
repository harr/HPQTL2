#' Check input values.
#'
#' @param X scaled genotype matrix
#' @param y phenotype vector
#' @param map map of markers
#' @param covar covariance matrix
#' @return Nothing if everything is ok. Error or warning otherwise.
#' @export


checkXy <- function(X,y=NULL,map=NULL,covar=NULL) {
  if (is.null(y)) y <- cbind(rep(1,nrow(X)))
  stopifnot(is.matrix(X))
  stopifnot(!any(is.na(X)) & !any(is.na(y))) 
  stopifnot(all(is.finite(X)) & all(is.finite(y)))
  stopifnot(nrow(X)==NROW(y))
  stopifnot(NCOL(y)==1)
  if (!is.null(map)) {
    stopifnot("chr" %in% names(map))
    stopifnot("pos" %in% names(map))
    stopifnot(!any(is.na(map)))
    stopifnot(ncol(X)==nrow(map))
    stopifnot(is.numeric(map$pos))
  }
  if (!is.null(covar)) {
    stopifnot(nrow(X)==nrow(covar))
    stopifnot(is.matrix(covar))
    stopifnot(!any(is.na(covar)) & all(is.finite(covar))) 
  }
}