#' One Locus Scan for Linear Mixed Model.
#'
#' @param y phenotype vector
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param covar matrix with covariates
#' @param G genetic similarity matrix
#' @return R/rql scanone object with LOD scores 
#' @export

scan.lmm <- function(y, X, map, covar = NULL, G = NULL) {
  
  checkXy(X, y, map, covar)
  output <- scanoneTemplate(map)
  n <- nrow(X)
  Z <- cbind(rep(1,n), covar)
  
  
  if (is.null(G)) G <- gensim(X, map, method="LMM")
  G <- normalize.matrix(G)
  fit <- regress(y~Z, ~G, pos=c(TRUE,TRUE))
  V <- fit$sigma["G"]*G + fit$sigma["In"]*diag(n)
  A <- half.inv(V)
  
  y.rot <- A %*% cbind(y)
  Z.rot <- A %*% Z
  X.rot <- A %*% X
  
  for (i in 1:ncol(X)) {
    rss0 <- sum(lsfit(y=y.rot, x=Z.rot, intercept=FALSE)$residuals^2)
    rss1 <- sum(lsfit(y=y.rot, x=cbind(Z.rot, X.rot[,i]), intercept=FALSE)$residuals^2)
    output$lod[i] <- n/2 * (log10(rss0) - log10(rss1))
  }  
  
  output
}