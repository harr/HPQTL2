#' Genetic similarity matrix for scaled genotype X.
#'
#' @param X scaled genotype matrix
#' @param map map of markers
#' @param method for which method is gensim. matrix calculated
#' @return a square matrix 
#' @export


gensim <- function(X, map, method=c("LMM", "LOCO", "POMOLOCO")) {
  method = match.arg(method)
  checkXy(X, map=map)
  
  if (method=="LMM") {
    return(tcrossprod(X))
  } 
  
  chrs <- unique(map$chr)
  Glist <- list()
  for (i in 1:length(chrs)) {
    chrs.ids <- map$chr == chrs[i]
    Glist[[i]] <- tcrossprod(X[,chrs.ids, drop=FALSE]) 
  }
  names(Glist) <- chrs
  
  if (method=="POMOLOCO") return(Glist)
  
  
  if (method == "LOCO") {
    Gloco <- list()
    for (i in 1:length(chrs)) {
      Gloco[[i]] <- Reduce("+", Glist[chrs!=chrs[i]]) 
    }
    return(Gloco)
  }
}