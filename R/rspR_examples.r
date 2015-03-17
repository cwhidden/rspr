## Testing rspr() and rspr.matrix()
source("rspR.r")
library(phangorn)
set.seed(4521435)

## Two-trees comparison
t1 <-  rtree(50, rooted = FALSE)
t2 <-  phangorn::rSPR(t1, 3) ## lets perform an SPR operation on t1
dist.topo(t1, t2)
rspr(t1, t2)
rspr(t1, t2, exact = FALSE)

## Matrix computation
trees <- rmtree(N = 10, n = 20)
trees[[2]] <- unroot(trees[[2]]) # introducing an unrooted black sheep
( radius <- rspr.matrix(trees, type = "restricted", maxdist = 1) )
isSymmetric(radius)
( radius5 <- rspr.matrix(trees, type = "restricted", maxdist = 5) )
( first.entry <- rspr.matrix(trees, type = "onerow") )
( full <- rspr.matrix(trees, type = "full") )

### Big matrix test
## Warning: this takes a long time. Don't run unless you really want to see if the code can handle big matrices ;0)
#big.trees <- rmtree(N = 4096, n = 20)
#( big.test <- rspr.matrix(big.trees, type = "restricted") )
#isSymmetric(big.test)
