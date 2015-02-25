#' One Locus Scan for Project One Marker Out Method.
#'
#' @param y phenotype vector
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param covar matrix with covariates
#' @param G genetic similarity matrix
#' @param faster overall estimate of heritability?
#' @return R/rql scanone object with LOD scores 
#' @export

scan.pomo <- function(y, X, map, covar = NULL, G = NULL, faster=TRUE) {
  
  checkXy(X, y, map, covar)
  output <- scanoneTemplate(map)
  n <- nrow(X)
  Z <- cbind(rep(1,n), covar)
  
  
  if (is.null(G)) G <- gensim(X, map, method="LMM")
  G <- normalize.matrix(G)
  if (faster) fit <- regress(y~Z, ~G, pos=c(TRUE,TRUE))
  
  for (i in 1:ncol(X)) {
    P <- diag(n) - X[,i] %*% t(X[,i]) / sum(X[,i]^2)
    Gp <- P %*% G %*% P
    if (!faster) fit <- regress(y~Z, ~Gp, pos=c(TRUE,TRUE))
    
    V <- fit$sigma[1]*Gp +  fit$sigma["In"]*diag(n)
    A <- half.inv(V)
    y.rot <- A %*% cbind(y)
    Z.rot <- A %*% Z
    
    rss0 <- sum(lsfit(y=y.rot, x=Z.rot, intercept=FALSE)$residuals^2)
    rss1 <- sum(lsfit(y=y.rot, x=cbind(Z.rot, A%*%X[,i,drop=FALSE]), intercept=FALSE)$residuals^2)
    output$lod[i] <- n/2 * (log10(rss0) - log10(rss1))
  }  
  
  output
}