#' Simulate phenotype with the given background heritability, add QTLs.
#'
#' @param y phenotype vector
#' @param X matrix with covariates
#' @param h2 background heritability
#' @param markers which markers 
#' @return numeric vector
#' @export

sim.pheno <- function(X, G, h2, markers = NULL, effects = NULL) {
  
  # to get reasonable interpretation, normalize matrix
  G <- normalize.matrix(G)
  n <- ncol(G)
  
  # non-normalized noise 
  noise <- rnorm(n)
  
  # non-normalized background
  background <-  mvrnorm(1, rep(0,n), G)
  
  # combine noise and background
  output <- preoutput <- noise + sqrt(h2/(1-h2))*background
  
  # add qtls
  if (!is.null(markers))
    
    if (!is.list(markers)) {
      for (i in 1:length(markers))
        output <- output + sqrt(effects[i]) * scale(X[,markers[i]])
    } else {
      output <- NULL
      for (j in 1:length(markers)) {
        tmp <- preoutput
        for (i in 1:length(markers[[j]]))
          tmp <- tmp + sqrt(effects[[j]][i]) * scale(X[,markers[[j]][i]])
        output <- cbind(output, tmp)
      }
    }
  output
}