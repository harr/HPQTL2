#' One Locus Scan for Project One Marker Out Method.
#'
#' @param y phenotype vector
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param covar matrix with covariates
#' @param Glist list of genetic similarity matrices
#' @param faster overall estimate of heritability?
#' @return R/rql scanone object with LOD scores 
#' @export

scan.pomoloco <- function(y, X, map, covar = NULL, Glist = NULL, faster=TRUE) {
  
  checkXy(X, y, map, covar)
  output <- scanoneTemplate(map)
  n <- nrow(X)
  Z <- cbind(rep(1,n), covar)
  
  if (is.null(Glist)) Glist <- gensim(X, map, method="LOCO")
  nchr <- length(Glist)
  last.chr <- ""
  
  for (i in 1:ncol(X)) {
    
    if (map$chr[i] != last.chr) {
      last.chr <- map$chr[i]
      Gchr <- Glist[[last.chr]]
      Grest <- Reduce("+", Glist[names(Glist)!=last.chr])
      Gp <- normalize.matrix(Grest) 
      if (faster) fit <- regress(y~Z, ~Gp, pos=c(TRUE,TRUE))  
    }    
    
    P  <- diag(n) - X[,i] %*% t(X[,i]) / sum(X[,i]^2)
    Gp <- normalize.matrix(P %*% Gchr %*% P + Grest) 
    if (!faster) fit <- regress(y~Z, ~Gp, pos=c(TRUE,TRUE)) 
    V <- fit$sigma[1]*Gp +  fit$sigma[2]*diag(n)
    A <- half.inv(V)
    
    y.rot <- A %*% cbind(y)
    Z.rot <- A %*% Z
    Xi.rot <- A %*% X[,i]
    
    rss0 <- sum(lsfit(y=y.rot, x=Z.rot, intercept=FALSE)$residuals^2)
    rss1 <- sum(lsfit(y=y.rot, x=cbind(Z.rot, Xi.rot), intercept=FALSE)$residuals^2)
    output$lod[i] <- n/2 * (log10(rss0) - log10(rss1))
  }  
  
  output
}