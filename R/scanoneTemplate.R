#' Create R/qtl scanone object
#'
#' @param map map of markers
#' @export

scanoneTemplate <- function(map) {
  output <- data.frame(chr=as.character(map$chr), 
                       pos=map$pos, 
                       lod=rep(0, nrow(map)), 
                       row.names=map$marker,
                       stringsAsFactors = FALSE)
  class(output) <- c("scanone", "data.frame")
  output
}