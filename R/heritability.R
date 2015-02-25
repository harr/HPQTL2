#' Estimate heritability.
#'
#' @param y phenotype vector
#' @param covar matrix with covariates
#' @param G genetic similarity matrix
#' @return numeric 
#' @export

heritability <- function(y, G, covar = NULL) {
  
  if (is.null(covar))
    covar <- rep(1, nrow(G))
  
  # to get reasonable interpretation, normalize matrix
  G <- normalize.matrix(G)
  
  # fit the mixed model
  rg.fit <- regress(y~covar, ~G, pos=c(TRUE,TRUE))
  
  # estimate heritability
  h2 <- as.numeric(rg.fit$sigma[1] / sum(rg.fit$sigma))
  
  h2
}