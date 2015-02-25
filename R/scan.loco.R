#' One Locus Scan for for Leave One Chromosome Out Method.
#'
#' @param y phenotype vector
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param covar matrix with covariates
#' @param Glist list of genetic similarity matrices
#' @return R/rql scanone object with LOD scores 
#' @export

scan.loco <- function(y, X, map, covar = NULL, Glist = NULL) {
  
  checkXy(X, y, map, covar)
  output <- scanoneTemplate(map)
  n <- nrow(X)
  Z <- cbind(rep(1,n), covar)
  
  if (is.null(Glist)) Glist <- gensim(X, map, method="LOCO")
  last.chr <- ""
  
  for (i in 1:ncol(X)) {
    
    if (map$chr[i] != last.chr) {
      last.chr <- map$chr[i]
      G <- normalize.matrix(Glist[[last.chr]])
      fit <- regress(y~Z, ~G, pos=c(TRUE,TRUE))
      V <- fit$sigma[1]*G +  fit$sigma[2]*diag(n)
      A <- half.inv(V)
      
      y.rot <- A %*% cbind(y)
      Z.rot <- A %*% Z
      
      chr.ids <- map$chr == last.chr
      X.rot <- A %*% X[,chr.ids]
      idx <- rep(0, ncol(X))
      idx[chr.ids] <- 1:sum(chr.ids)
    }  
    
    rss0 <- sum(lsfit(y=y.rot, x=Z.rot, intercept=FALSE)$residuals^2)
    rss1 <- sum(lsfit(y=y.rot, x=cbind(Z.rot, X.rot[,idx[i]]), intercept=FALSE)$residuals^2)
    output$lod[i] <- n/2 * (log10(rss0) - log10(rss1))
  }  
  
  output
}