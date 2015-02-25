#' Simulate F2 cross genotypes
#'
#' @param n.ind number of animals
#' @param chrlen length of chromosomes
#' @param nchr number of chromosomes
#' @param nmar number of markers per chromosome 
#' @return list with genotype probabilities and map of markers 
#' @export

sim.geno.f2 <- function(n.ind, chrlen=100, nchr=20, nmar=30) {
  
  map <- sim.map(len=rep(chrlen, nchr), n.mar=nmar)
  cross <- sim.cross(map, type="f2", n.ind=n.ind, model = NULL)
  cross <- calc.genoprob(cross)
  
  chrtype <- sapply(cross$geno, class)
  if(any(chrtype =="X")) {
    for(i in which(chrtype=="X"))
      cross$geno[[i]]$prob <- qtl:::reviseXdata(class(cross)[1], "simple", cross$pheno, prob=cross$geno[[i]]$prob)
  }
  
  subjects = as.character(1:nrow(cross$pheno))
  markers <- data.frame(marker = do.call("c", lapply(cross$geno, function(x) names(x$map))),
                        chr = do.call("c", mapply(rep, x = names(nmar(cross)), each = nmar(cross), SIMPLIFY = FALSE)),
                        pos = do.call("c", lapply(cross$geno, function(x) as.vector(x$map))))
  
  list.of.probs <- lapply(cross$geno, function(x) aperm(x$prob, c(1,3,2)))
  if (length(unique(sapply(list.of.probs, ncol)))!=1) {
    max.calls <- max(sapply(list.of.probs, ncol))
    for (i in which(sapply(list.of.probs, ncol) < max.calls)) {
      miss.calls <- max.calls - dim(list.of.probs[[i]])[2] 
      zeros <- array(rep(0,dim(list.of.probs[[i]])[1] * miss.calls * dim(list.of.probs[[i]])[3]))
      dim(zeros) <- c(dim(list.of.probs[[i]])[1], miss.calls,  dim(list.of.probs[[i]])[3])
      list.of.probs[[i]] <- abind(list.of.probs[[i]], zeros, along = 2)
    }
  }
  
  geno <- list(probs = abind(list.of.probs, along=3), 
               subjects = subjects,    
               calls = dimnames(cross$geno[[1]]$prob)[[3]],
               markers = markers,
               chromosomes = data.frame(chr = names(cross$geno), type = as.vector(chrtype)))
  
  geno$X <- scale(apply(geno$probs, c(1,3), function(x) x[1] - x[2]))
  
  return(geno)
}