context('grav2 dataset')

grav2 <- read_cross2(paste0(.libPaths(),"/qtl2/extdata/grav2.zip"))
grav2.probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
grav2.aprobs <- genoprob_to_alleleprob(grav2.probs)

grav2.X <- NULL
grav2.map <- NULL
for (i in 1:length(grav2.aprobs)) {
  tmp <- apply(grav2.aprobs[[i]], c(1,2), function(x) x[1] - x[2])
  grav2.X <- cbind(grav2.X, scale(tmp))
  grav2.map <- rbind(grav2.map, 
                     data.frame(chr=i, pos=1:ncol(tmp),
                                name=colnames(tmp)))
}  
grav2.y <- grav2$pheno[,241]

G <- gensim(grav2.X, grav2.map, method="LMM")
Glist <- gensim(grav2.X, grav2.map, method="LOCO")
Glist2 <- gensim(grav2.X, grav2.map, method="POMOLOCO")

test_that("scan.lm",{
  
  fit.lm <- scan.lm(grav2.y, grav2.X, grav2.map)
  
  # compare to saved version
  correct_lod <- readRDS("grav2_lm.rds") # saveRDS(fit.lm$lod, file="grav2_lm.rds")
  expect_equal(fit.lm$lod, correct_lod, tolerance=0.01)
})


test_that("scan.lmm",{
  
  fit.lmm <- scan.lmm(grav2.y, grav2.X, grav2.map)
  
  # compare to saved version
  correct_lod <- readRDS("grav2_lmm.rds") # saveRDS(fit.lmm$lod, file="grav2_lmm.rds")
  expect_equal(fit.lmm$lod, correct_lod, tolerance=0.01)
})

test_that("scan.loco",{
  
  fit.loco <- scan.loco(grav2.y, grav2.X, grav2.map)
  
  # compare to saved version
  correct_lod <- readRDS("grav2_loco.rds") # saveRDS(fit.loco$lod, file="grav2_loco.rds")
  expect_equal(fit.loco$lod, correct_lod, tolerance=0.01)
})





context('iron dataset')

iron <- read_cross2(paste0(.libPaths(),"/qtl2/extdata/iron.zip"))
iron.probs   <-  calc_genoprob(iron, step=1, error_prob=0.002)
iron.aprobs <- genoprob_to_alleleprob(iron.probs)

iron.X <- NULL
iron.map <- NULL
for (i in 1:length(iron.aprobs)) {
  tmp <- apply(iron.aprobs[[i]], c(1,2), function(x) x[1] - x[2])
  iron.X <- cbind(iron.X, scale(tmp))
  iron.map <- rbind(iron.map, 
                     data.frame(chr=i, pos=1:ncol(tmp),
                                name=colnames(tmp)))
}  
covar <- as.matrix(model.matrix(~sex+cross_direction, iron$covar)[,-1])
iron.y <- iron$pheno[,"spleen"]

G <- gensim(iron.X, iron.map, method="LMM")
Glist <- gensim(iron.X, iron.map, method="LOCO")
Glist2 <- gensim(iron.X, iron.map, method="POMOLOCO")

test_that("scan.lm",{
  
  fit.lm <- scan.lm(iron.y, iron.X, iron.map, covar)
  
  # compare to saved version
  correct_lod <- readRDS("iron_lm.rds") # saveRDS(fit.lm$lod, file="iron_lm.rds")
  expect_equal(fit.lm$lod, correct_lod, tolerance=0.01)
})


test_that("scan.lmm",{
  
  fit.lmm <- scan.lmm(iron.y, iron.X, iron.map, covar)
  
  # compare to saved version
  correct_lod <- readRDS("iron_lmm.rds") # saveRDS(fit.lmm$lod, file="iron_lmm.rds")
  expect_equal(fit.lmm$lod, correct_lod, tolerance=0.01)
})

test_that("scan.loco",{
  
  fit.loco <- scan.loco(iron.y, iron.X, iron.map, covar)
  
  # compare to saved version
  correct_lod <- readRDS("iron_loco.rds") # saveRDS(fit.loco$lod, file="iron_loco.rds")
  expect_equal(fit.loco$lod, correct_lod, tolerance=0.01)
})