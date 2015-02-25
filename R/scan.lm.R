#' One Locus Scan for Linear Model.
#'
#' @param y phenotype vector
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param covar matrix with covariates
#' @return R/rql scanone object with LOD scores 
#' @export

scan.lm <- function(y, X, map, covar = NULL) {
  
  checkXy(X, y, map, covar)
  output <- scanoneTemplate(map)
  n <- nrow(X) # number of animals
  Z <- cbind(rep(1, n), covar) # covariates and intercept
  
  for (i in 1:ncol(X)) {
    rss0 <- sum(lsfit(y=y, x=Z, intercept=FALSE)$residuals^2)
    rss1 <- sum(lsfit(y=y, x=cbind(Z, X[,i]), intercept=FALSE)$residuals^2)
    output$lod[i] <- n/2 * (log10(rss0) - log10(rss1))
  }  
  
  output
}